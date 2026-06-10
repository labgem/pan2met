"""
Some helper functions to deal with the taxonomic range.

The taxonomic "ground truth" is taken from NCBI Taxonomy.
"""

from typing import Dict, Set
from pathlib import Path
import csv


class NCBITaxonomyTree:
    """
    A parser and utilitary methods on NCBI taxonomy tree.
    """

    def __init__(self, dump_path: str):
        """
        :param dump_path: path to the folder containing the dump of the NCBI Taxonomy database.
        """
        self.parse_parent_dict(dump_path)
        self.root_tax_id = 1

    def parse_parent_dict(self, dump_path: str):
        """
        :param dump_path: path to the directory with the dump of the NCBI Taxonomy database.
        load a dictionnary linking a taxid to its direct parent in the NCBI Taxonomy tree.
        """
        self.parent_dict: Dict[int, int] = {}
        self.tax_rank: Dict[int, str] = {}
        self.species_taxa: Set[int] = set()

        with open(Path(dump_path) / "nodes.dmp", "r") as nodes_file:
            reader = csv.DictReader(
                nodes_file,
                fieldnames=[
                    "tax_id",
                    "parent tax_id",
                    "rank",
                    "embl code",
                    "division id",
                    "inherited div flag",
                    "genetic code id",
                    "inherited GC flag",
                    "mitochondrial genetic code id",
                    "inherited MGC flag",
                    "GenBank hidden flag",
                    "hidden subtree root flag",
                    "comments",
                ],
                delimiter="|",
            )
            for row in reader:
                tax_id = int(row["tax_id"].replace("\t", ""))
                parent_tax_id = int(row["parent tax_id"].replace("\t", ""))
                if row["rank"] != "no rank":
                    if row["rank"] == "species":
                        self.species_taxa.add(tax_id)

                self.parent_dict[tax_id] = parent_tax_id
                self.tax_rank[tax_id] = row["rank"].replace("\t", "")

    def is_child_of_parent_tax_id(self, parent_tax_id: int, child_tax_id: int) -> bool:
        """

        :param parent_tax_id: NCBI Taxonomy identifier
        :param child_tax_id:  NCBI Taxonomy identifier
        :return: True if parent_tax_id is a parent of child_tax_id in the NCBI Taxonomy tree
        """
        if parent_tax_id not in self.parent_dict:
            raise ValueError(
                f"NCBI Taxonomy identifier {parent_tax_id} not found in parsed taxonomy tree."
            )
        if child_tax_id not in self.parent_dict:
            raise ValueError(
                f"NCBI Taxonomy identifier {child_tax_id} not found in parsed taxonomy tree."
            )
        current_tax_id = child_tax_id
        if current_tax_id == parent_tax_id:
            return True
        while current_tax_id != self.root_tax_id:
            if current_tax_id == parent_tax_id:
                return True
            current_tax_id = self.parent_dict[current_tax_id]
        return False

    def path_to_root(self, tax_id: int) -> list[int]:
        """
        Trace back the path to the root of the tree

        :param tax_id:
        :return: The list of NCBI Taxonomy identifier from tax_id to the root of the NCBI Taxonomy
        """
        if tax_id not in self.parent_dict:
            raise ValueError(
                f"NCBI Taxonomy identifier {tax_id} not found in parsed taxonomy tree."
            )
        path = []
        current = tax_id
        while current != self.root_tax_id:
            path.append(current)
            current = self.parent_dict[current]
        return path

    def last_common_ancestor(self, taxid_1: int, taxid_2: int) -> int:
        """
        Return the NCBI-Taxonomy tax id that is the shallowest
        in the NCBI-Taxonomy tree that is a parent of both tax id

        Naive implementation: traverse the tree up to the root from both nodes
        and then go down the tree on the two paths from the root,
        while the node is a shared ancestor.
        It is expected to take a time in the order of the depth of the tree.

        :param taxid_1: an NCBI Taxonomy identifier
        :param taxid_2: an other NCBI Taxonomy identifier
        :return: the deepest taxonomy identifier in the NCBI Taxonomy tree that is a common ancestor of both tax id.
        """
        if taxid_1 == taxid_2:
            return taxid_1
        path1 = self.path_to_root(taxid_1)
        path2 = self.path_to_root(taxid_2)

        path1 = path1[::-1]
        path2 = path2[::-1]

        common_ancestor = self.root_tax_id
        i = 0
        while i < len(path1) and i < len(path2) and path1[i] == path2[i]:
            common_ancestor = path1[i]
            i += 1
        return common_ancestor

    def is_under_same_species(self, taxid_1: int, taxid_2: int) -> bool:
        """
        Check whether two NCBI Taxonomy identifiers are under the same species subtree.

        We arbitrarily assume that a species subtree is at most 5 items deep, and do not explore it further.

        :param taxid_1: a NCBI Taxonomy identifier
        :param taxid_2: a NCBI Taxonomy identifier
        :return: True when a parent of the NCBI Taxonomy identifiers is of rank species
        """
        if taxid_1 == taxid_2:
            return True
        MAX_RECURSE: int = 5
        common_ancestor: int = self.last_common_ancestor(taxid_1, taxid_2)
        current = common_ancestor
        depth = 0  # could be even already greater than that
        while depth < MAX_RECURSE:
            if current in self.species_taxa:
                return True
            current = self.parent_dict[current]
            depth += 1
        return False
