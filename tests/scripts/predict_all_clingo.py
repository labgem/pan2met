"""
Predict for all generated test cases and make a report
"""

import os

from pan2met.utils import read_list, write_output
from pan2met.inference.constraint.metabolome import ASPPathwayInference
from pan2met.io.knowledge_base import select_kb

TEST_CASES_FOLDER = "tests/cases/generated/"
KB_ASP = "tmp/metabiantes_pathway_asp.lp"


def infer_metabolome(reactome: set[str]) -> list[str]:
    kb = select_kb("metabiantes")
    inference = ASPPathwayInference(kb, reactome)
    inferred_pathways: list[str] = inference.inferred_pathways(KB_ASP)
    return inferred_pathways


def main():
    test_cases = os.listdir(TEST_CASES_FOLDER)
    for test_case in test_cases:
        if os.path.isdir(os.path.join(TEST_CASES_FOLDER, test_case)):
            print(f"Running prediction for {test_case}")
            reactome_filename = os.path.join(TEST_CASES_FOLDER, test_case, "reactome")
            reactome = set(read_list(reactome_filename))
            output_filename = os.path.join(
                TEST_CASES_FOLDER,
                test_case,
                "pan2met-asp.pathways.list",
            )
            metabolome = infer_metabolome(reactome)
            write_output(output_filename, metabolome)
    print("done.")


if __name__ == "__main__":
    main()
