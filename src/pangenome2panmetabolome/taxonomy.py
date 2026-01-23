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
        self.parent_dict = self.parse_parent_dict_even_faster(dump_path)
        self.root_tax_id = 1

    def parse_parent_dict(self, dump_path: str) -> dict[int, int]:
        """
        Return a dictionnary linking a taxid to its direct parent in the NCBI Taxonomy tree.
        """
        parent_dict: dict[int, int] = {}

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
                    "inherited GC  flag",
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
                parent_dict[tax_id] = parent_tax_id
        return parent_dict

    def parse_parent_dict_faster(self, dump_path: str) -> dict[int, int]:
        parent_dict: dict[int, int] = {}

        with open(Path(dump_path) / "nodes.dmp", "r") as nodes_file:
            reader = csv.reader(
                nodes_file,
                delimiter="|",
            )
            for row in reader:
                tax_id = int(row[0].replace("\t", ""))
                parent_tax_id = int(row[1].replace("\t", ""))
                parent_dict[tax_id] = parent_tax_id
        return parent_dict

    def parse_parent_dict_even_faster(self, dump_path: str) -> dict[int, int]:
        parent_dict: dict[int, int] = {}

        with open(Path(dump_path) / "nodes.dmp", "r") as nodes_file:
            for row in nodes_file:
                row = row[:20].split("\t|\t")
                tax_id = int(row[0])
                parent_tax_id = int(row[1])
                parent_dict[tax_id] = parent_tax_id
        return parent_dict

    def is_child_of_parent_tax_id(self, parent_tax_id: int, child_tax_id: int) -> bool:
        current_tax_id = child_tax_id
        if current_tax_id == parent_tax_id:
            return True
        while current_tax_id != self.root_tax_id:
            if current_tax_id == parent_tax_id:
                return True
            current_tax_id = self.parent_dict[current_tax_id]
        return False
