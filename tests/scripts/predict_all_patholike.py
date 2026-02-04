"""
Predict for all generated test cases and make a report
"""

import os

from pangenome2panmetabolome.utils import read_list, write_output
from pangenome2panmetabolome.inference.metabolome import infer_metabolome


TEST_CASES_FOLDER = "tests/cases/generated/"
ECOLI_TAXID = 562


def main():
    test_cases = os.listdir(TEST_CASES_FOLDER)
    for test_case in test_cases:
        if os.path.isdir(os.path.join(TEST_CASES_FOLDER, test_case)):
            print(f"Running prediction for {test_case}")
            reactome_filename = os.path.join(TEST_CASES_FOLDER, test_case, "reactome")
            reactome = set(read_list(reactome_filename))
            output_filename = os.path.join(
                TEST_CASES_FOLDER, test_case, "pangenome2panmetabolome.pathways.list"
            )
            metabolome = infer_metabolome(reactome, ECOLI_TAXID, None)
            write_output(output_filename, metabolome)
    print("done.")


if __name__ == "__main__":
    main()
