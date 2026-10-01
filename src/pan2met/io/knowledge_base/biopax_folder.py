"""
A BioPAX backed metabolic pathway knowledge base.

Warning: As it is currently implemented, the list of identifiers are the IRI of BioPAX entities.
I will have to make sure this is concordant with the list of provided reaction identifiers.
For instance, when the identifier type is MetaCyc, it will not work in conjonction with BioPAX files,
which not necessarily includes all the required cross references.
In pathbank, for instance reaction would most probably contain only links to Rhea identifiers.
"""

import glob
import importlib.resources
import os
from pathlib import Path

import duckdb
import ouisparql

import pan2met
from pan2met.io.knowledge_base import KnowledgeBase, unsupported
from pan2met.utils import logger


def load_duckdb_sparql_extension(conn):
    conn.sql("INSTALL rdf FROM community; LOAD rdf;")


def duckdb_owl_to_parquet(conn, owl_filepath: Path, parquet_filepath: Path):
    """
    Dump a OWL file into a parquet file with duckdb.
    """
    load_duckdb_sparql_extension(conn)

    conn.sql(
        f"""
        COPY (
            SELECT * FROM read_rdf('{owl_filepath.as_posix()}', file_type='xml')
        ) TO '{parquet_filepath.as_posix()} (FORMAT PARQUET);
        """
    )


def duckdb_set_view(conn, table_name: str, new_view_name: str):
    conn.sql(f"""
            DROP VIEW IF EXISTS {new_view_name};
            CREATE VIEW {new_view_name} AS (
                SELECT * FROM {table_name}
            );
            """)


def load_biopax_as_view(conn, biopax_filename: str, view_name: str):
    conn.sql(f"""
            DROP VIEW IF EXISTS {view_name};
            CREATE VIEW {view_name} AS
            SELECT * FROM read_rdf("{biopax_filename}", file_type='xml');
            """)


def basename_from_filename(pathway_filename: str) -> str:
    return os.path.splitext(os.path.basename(pathway_filename))[0]


# we would type it as tuple_extract_value[T](tuples: Iterable[Tuple[T, ...]], index: int) -> Iterable[T] starting from Python3.12
def tuple_extract_value(tuples, index: int):
    for tuple in tuples:
        yield tuple[index]


class BioPAXFolderKnowledgeBase(KnowledgeBase):
    """
    This KnowledgeBase leverages DuckDB rdf extension to manage a metabolic knowledgebase backed by a set of BioPAX files.
    One example of application of such a knowledgebase would be to use it on the BioPAX dump of the PathBank database (https://pathbank.org/).
    """

    def __init__(self, config):
        self.biopax_folder = Path(config["reference"]["biopax_folder"]["path"])
        # Create a DuckDB connection in memory and load DuckDB's rdf extension
        self.conn = duckdb.connect(":memory:")
        load_duckdb_sparql_extension(self.conn)

        # Load the SPARQL queries
        with importlib.resources.path(
            pan2met, "sparql/biopax/queries.rq"
        ) as sparql_queries_filepath:
            self.queries = ouisparql.from_path(
                sparql_queries_filepath, driver_adapter="duckdb_rdf"
            )

        # Load all the triplets from the OWL files available in the directory
        biopax_glob_path = self.biopax_folder / "*.owl"
        logger.debug(
            f"Loading BioPAX OWL file from glob path {self.biopax_folder.as_posix()}"
        )
        self.conn.sql(
            f"""
            create table triplets AS select * from read_rdf('{biopax_glob_path.as_posix()}', file_type='xml')
            """
        )

    # Helper functions

    def duckdb_checkout(self, table_name: str):
        duckdb_set_view(self.conn, table_name, "triplets")

    def list_pathway_folder(self) -> list[str]:
        return glob.glob(os.path.join(self.biopax_folder, "*.owl"))

    # KnowledgeBase API:

    def monomers(self) -> list[str]:
        """
        List monomer polypeptides
        """
        results = self.queries.get_protein(self.conn)
        return tuple_extract_value(results, index=0)

    def pathways(self) -> list[str]:
        """
        List the pathways referenced by the knowledge base
        """
        results = self.queries.get_all_pathways(self.conn)
        return tuple_extract_value(results, 0)

    def reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List reactions of a pathway
        """
        # As we are not certain of the proper location of the pathway definition
        pathway_uri = f"<{pathway_id}>"
        results = self.queries.get_reactions_of_pathway(
            self.conn, pathway_uri=pathway_uri
        )
        logger.warning(
            "We listed the URI of reactions, maybe not the reaction identifiers used in the protein-reaction mapping?"
        )
        return tuple_extract_value(results, index=0)

    def name_of_pathway(self, pathway_id: str) -> str:
        """
        Get the name of a pathway
        """
        pathway_uri = f"<{pathway_id}>"
        results = self.queries.get_name_of_pathway(self.conn, pathway_uri=pathway_uri)
        return tuple_extract_value(results, index=0)[0]

    def spontaneous_reactions(self) -> list[str]:
        """
        List spontaneous reactions
        """
        results = self.queries.get_spontaneous_reactions(self.conn)
        return tuple_extract_value(results, index=0)

    def orphan_reactions(self) -> list[str]:
        """
        List orphan reactions

        TODO: think about how to handle this.
        The list of orphan reactions, considering only the provided
        list of BioPAX might be different from the list of orphan reactions all available knowledge considered
        (MetaCyc and Rhea/Uniprot for instance)
        """
        results = self.queries.get_orphan_reactions(self.conn)
        return tuple_extract_value(results, index=0)

    def non_spontaneous_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        Non-spontaneous reactions of a pathway
        """
        pathway_uri = f"<{pathway_id}>"
        results = self.queries.get_spontaneous_reactions_of_pathway(
            self.conn, pathway_uri
        )
        return tuple_extract_value(results, index=0)

    def non_orphan_non_spontaneous_reactions_of_pathway(
        self, pathway_id: str
    ) -> list[str]:
        """
        Non-orphan and non-spontaneous reactions of a pathway
        """
        pathway_uri = f"<{pathway_id}>"
        results = self.queries.get_non_orphan_non_spontaneous_reactions_of_pathway(
            self.conn, pathway_uri
        )
        return tuple_extract_value(results, index=0)

    def enzymes_of_reaction(self, reaction_id: str) -> list[str]:
        """
        List enzymes catalyzing a reaction
        """
        reaction_uri = f"<{reaction_id}>"
        results = self.queries.get_enzymes_of_reaction(
            self.conn, reaction_uri=reaction_uri
        )
        return tuple_extract_value(results, index=0)

    def reactions_by_ec_number(self, ec_number: str) -> list[str]:
        """
        List reactions annotated with given EC-number

        The bare EC-number should be given without the "EC-" prefix we can find in BioCyc databases.
        """
        # Wrap the ec_number string in double quotes
        ec_number = f'"{ec_number}"'
        results = self.queries.get_reactions_by_ec_number(
            self.conn, ec_number=ec_number
        )
        return tuple_extract_value(results, 0)

    @unsupported
    def pathway_taxonomic_range(self, pathway_id: str) -> int | None:
        """
        Return a NCBI-Taxonomy Taxonomy Identifier number
        """

    @unsupported
    def reaction_is_key(self, pathway_id: str, reaction_id: str) -> bool:
        """
        Return True if the reaction is a key reaction of the pathway
        """

    @unsupported
    def key_reactions_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List all key reactions of a pathway.
        """

    @unsupported
    def pathway_reaction_order(self, pathway_id: str) -> list[tuple[str, str]]:
        """
        Return a list of (predecessor, successor) reactions of a pathway.
        """

    def pathways_with_reaction(self, reaction_id: str) -> list[str]:
        """
        List all pathways with the given reaction identifier.
        """
        reaction_uri = f"<{reaction_id}>"
        results = self.queries.pathways_with_reaction(
            self.conn, reaction_uri=reaction_uri
        )
        return tuple_extract_value(results, index=0)

    @unsupported
    def variants_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the variants of a pathway.
        """

    @unsupported
    def species_evidence_of_pathway(self, pathway_id: str) -> list[int]:
        """
        List species NCBI Taxonomy identifiers where pathway presence evidence were found.
        """

    @unsupported
    def ontology_parent_class_of_pathway(self, pathway_id: str) -> list[str]:
        """
        List the pathway ontology parents of a pathway
        """

    def complex(self) -> list[str]:
        """
        List the protein complex
        """
        results = self.queries.get_complex(self.conn)
        return tuple_extract_value(results, index=0)

    def components_of_complex(self, complex_id: str) -> list[str]:
        """
        List the components of a complex
        """
        complex_uri = f"<{complex_id}>"
        results = self.queries.get_complex_components(self.conn, complex_uri)
        return tuple_extract_value(results, index=1)
