#!/usr/bin/env python3

"""
Get data from BioCyc PGDB

Before running this script, launch pathway-tools python API with
   pathway-tools -lisp -python-local-only-non-strict

"""

import logging
import argparse
from typing import Iterable

from tqdm import tqdm


logger = logging.getLogger("pan2met:io:metacyc")
logger.setLevel(logging.DEBUG)

PGDB_ID: str = "META"

def is_spontaneous(pgdb, reaction: str) -> bool:
    """
    Return True if and only if reaction spontaneous_p predicate is True, otherwise returns False.
    """
    pgdb[reaction].spontaneous_p is not None and pgdb[reaction].spontaneous_p


def is_orphan(pgdb, reaction: str) -> bool:
    """
    Return true of value 'orphan_p' is either "YES-CONFIRMED" or "YES-PUTATIVE", False if orphan_p is "NO".

    Raises a value error if orphan_p is neither of these three cases. If orphan_p is None, return None.
    """
    if pgdb[reaction].orphan_p is None:
        logger.warning(f"{reaction} reaction's 'orphan_p' is None")
        return None
    orphan_p = pgdb[reaction].orphan_p[0]
    if orphan_p in ["|YES-CONFIRMED|", "|YES-PUTATIVE|"]:
        return True
    elif orphan_p == "|NO|":
        return False
    else:
        raise ValueError(
            "'orphan_p' predicate is neither YES-CONFIRMED, YES-PUTATIVE nor NO: "
            + orphan_p
        )


def remove_pipes(identifier: str) -> str:
    return identifier.replace("|", "")


def get_pgdb_reactions(pgdb) -> list[str]:
    return pgdb.all_rxns(type_of_reactions=":all")


def get_pgdb_pathways(pgdb) -> list[str]:
    return pgdb.all_pathways()


def enzymes_of_reaction(pgdb, reaction: str) -> list[str]:
    return pgdb.enzymes_of_reactions(reaction)


def get_monomers_of_enzyme(pgdb, enzyme: str) -> list[str]:
    return pgdb.monomers_of_protein(enzyme, unmodify=True)


def is_pathway(pgdb, pathway: str) -> bool:
    return pgdb[pathway]["instance_name_template"] == "PWY-*"


def get_reactions_of_pathway(pgdb, pathway: str) -> list[str]:
    frame_object = pgdb.get_frame_objects([pathway])[0]
    reaction_list = frame_object["reaction_list"]
    # reaction_list can contain subpathways, so we return also the reactions of the subpathays:
    # FIXME: this might be not the best way to deal with such situations.
    mask = [is_pathway(pgdb, reaction) for reaction in reaction_list]
    subpathways = [
        pathway for index, pathway in enumerate(reaction_list) if mask[index]
    ]
    reactions = [
        pathway for index, pathway in enumerate(reaction_list) if not mask[index]
    ]
    for subpathway in subpathways:
        subreactions = get_reactions_of_pathway(pgdb, subpathway)
        reactions = reactions + subreactions
    return reactions


def insert_pathway(
    conn, queries, metacyc_pathway_id: str, metacyc_pathway_name: str
) -> str:
    """
    Insert a pathway in SQL database, and return the created record id.
    """
    db_pathway_id = queries.insert_pathway(conn, name=metacyc_pathway_name)
    queries.insert_crossref(
        conn,
        local_table="pathway",
        local_id=db_pathway_id,
        reference_database="MetaCyc",
        reference_id=metacyc_pathway_id,
    )
    return db_pathway_id


def insert_reaction(
    conn, queries, metacyc_reaction_id: str, metacyc_reaction_name: str
) -> str:
    """
    Insert a reaction in SQL database, and return the created record id.
    """
    db_reaction_id = queries.insert_reaction(conn, name=metacyc_reaction_name)
    queries.insert_crossref(
        conn,
        local_table="reaction",
        local_id=db_reaction_id,
        reference_database="MetaCyc",
        reference_id=metacyc_reaction_id,
    )
    return db_reaction_id


def add_reaction_to_pathway(conn, queries, db_pathway_id: int, db_reaction_id: int):
    """
    Insert a relation pathway - reaction in SQL database
    """
    queries.add_reaction_to_pathway(
        conn, pathway_id=db_pathway_id, reaction_id=db_reaction_id
    )


def name_of_frame(pgdb, frame_id: str) -> str:
    frame_object = pgdb.get_frame_objects([frame_id])[0]
    if frame_object["names"] is not None and len(frame_object["names"]) >= 1:
        return frame_object["names"][0]
    else:
        return f"unnamed: {frame_id}"


def load_pathway_data_into_sqlite(conn, queries, pgdb):
    logger.info("Loading pathways data into sqlite")
    # Load pathways
    pathways = get_pgdb_pathways(pgdb)
    for pathway_id in tqdm(pathways):
        pathway_id_norm = remove_pipes(pathway_id)
        pathway_name = name_of_frame(pgdb, pathway_id)
        db_pathway_id = insert_pathway(conn, queries, pathway_id_norm, pathway_name)
        for reaction_id in get_reactions_of_pathway(pgdb, pathway_id):
            reaction_id_norm = remove_pipes(reaction_id)
            reaction_name = name_of_frame(pgdb, reaction_id)
            db_reaction_id = insert_reaction(
                conn, queries, reaction_id_norm, reaction_name
            )
            add_reaction_to_pathway(conn, queries, db_pathway_id, db_reaction_id)
    logger.info("Loading done.")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Load MetaCyc data into a SQLite database."
    )
    parser.add_argument("-o", "--outdb", help="Output SQLite database")
    parser.add_argument(
        "-f",
        "--force",
        action=argparse.BooleanOptionalAction,
        help="Overwrite existing database SQLite file, if any.",
    )
    return parser.parse_args(), parser


def list_all_monomers(pgdb):
    return list(
        filter(
            lambda monomer: "MONOMER" in monomer,
            pgdb.get_class_all_instances("|Polypeptides|"),
        )
    )


def is_proteic_complex(pgdb, enzyme: str) -> bool:
    """
    Assert True if `enzyme` is a protein complex.
    """
    return pgdb.complex_p(enzyme)


def proteic_complex_subunits(pgdb, enzyme: str) -> Iterable[str]:
    """
    Iter over enzyme protein complex subunits
    """
    enzyme_frame_objects = pgdb.get_frame_objects([enzyme])[0]
    components = enzyme_frame_objects["components"]
    return components


def is_homomeric(pgdb, polymer: str) -> bool:
    components = pgdb.get_frame_objects([polymer])[0]["components"]
    return (
        is_proteic_complex(pgdb, polymer)
        and components is not None
        and len(components) == 1
    )


class SingletonMetaCycPGDB:
    def __new__(cls):
        if not hasattr(cls, "instance"):
            cls.instance = super(SingletonMetaCycPGDB, cls).__new__(cls)
        return cls.instance
