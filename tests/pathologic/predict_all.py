"""
Predict for all generated test cases and make a report
"""

import os

from pangenome2panmetabolome.utils import read_list, write_output
from pangenome2panmetabolome.io import knowledge_base
from pangenome2panmetabolome.config import config
from pangenome2panmetabolome.inference.metabolome import PathwayInference


TEST_CASES_FOLDER = "tests/cases/generated2/"


def run_prediction(kb, input_reactome, output_filename):
    inference = PathwayInference(kb, input_reactome, taxon_id=1)
    metabolome = inference.inferred_pathways()
    write_output(output_filename, metabolome)


def main():
    kb = knowledge_base.select_kb(config["reference"]["source"])
    test_cases = os.listdir(TEST_CASES_FOLDER)
    for test_case in test_cases:
        if os.path.isdir(os.path.join(TEST_CASES_FOLDER, test_case)):
            print(f"Running prediction for {test_case}")
            reactome_filename = os.path.join(TEST_CASES_FOLDER, test_case, "reactome")
            reactome = set(read_list(reactome_filename))
            output_filename = os.path.join(
                TEST_CASES_FOLDER, test_case, "pangenome2panmetabolome.pathways.list"
            )
            run_prediction(kb, reactome, output_filename)
    print("done.")


if __name__ == "__main__":
    main()
