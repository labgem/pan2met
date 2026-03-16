"""
Predict for all generated test cases and make a report
"""

import os

from pan2met.utils import read_list, write_output
from pan2met.inference.metabolome import infer_metabolome


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
                TEST_CASES_FOLDER, test_case, "pan2met.pathways.list"
            )
            reason_filename = os.path.join(
                TEST_CASES_FOLDER, test_case, "reason.pan2met.log"
            )
            metabolome = infer_metabolome(reactome, ECOLI_TAXID, reason_filename)
            write_output(output_filename, metabolome)
    print("done.")


if __name__ == "__main__":
    main()
