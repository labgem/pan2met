"""
In PathwayTools PathoLogic,
expected Taxonomic range is a particularly important
predictor used in reaction / pathway inference prediction.

This module introduces some helper functions to deal with the taxonomic range.

The taxonomic "ground truth" is taken from NCBI Taxonomy, as in PathwayTools.
"""

from pathlib import Path
import csv


class NCBITaxonomyTree:
    def __init__(self, dump_path: str):
        self.parse_parent_dict(dump_path)
        self.root_tax_id = 1

    def parse_parent_dict(self, dump_path: str):
        """
        Return a dictionnary linking a taxid to its direct parent in the NCBI Taxonomy tree.
        """
        self.parent_dict: dict[int, int] = {}
        self.tax_rank: dict[int, str] = {}

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
                self.parent_dict[tax_id] = parent_tax_id
                self.tax_rank[tax_id] = row["rank"].replace("\t", "")

    def is_child_of_parent_tax_id(self, parent_tax_id: int, child_tax_id: int) -> bool:
        current_tax_id = child_tax_id
        if current_tax_id == parent_tax_id:
            return True
        while current_tax_id != self.root_tax_id:
            if current_tax_id == parent_tax_id:
                return True
            current_tax_id = self.parent_dict[current_tax_id]
        return False

    def path_to_root(self, taxid: int) -> list[int]:
        """
        Trace back the path to the root of the tree
        """
        path = []
        current = taxid
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
