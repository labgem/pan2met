"""
Extract a simplified table with columns pathway, decision, reason from the Lisp dump of PathoLogic decisions.

This Lisp dump is written to file: ptools-local/pgdbs/user/testcyc/1.0/reports/pwy-inference-description.data,
when testcyc is the name of the generated PGDB.
"""

import sys
from typing import Iterable
from typing import Optional


def parse_lisp_token(lisp: str) -> Iterable[str]:
    acc = ""
    for letter in lisp:
        if letter in "()":
            if acc != "":
                yield acc
                acc = ""
            yield letter
        elif letter == " ":
            if acc != "":
                yield acc
                acc = ""
        else:
            acc += letter


def first_order_lisp_objects(lisp: str) -> Iterable[list[str]]:
    depth = 0
    objects = []
    current_object = []
    for token in parse_lisp_token(lisp):
        if token == "(":
            depth += 1
        elif token == ")":
            depth -= 1
        current_object.append(token)
        if depth == 0 and current_object != []:
            objects.append(current_object)
            current_object = []
    return objects


def find_row_by_start(rows: list[str], start: str) -> int:
    for index, item in enumerate(rows):
        if item.startswith(start):
            return index


def report(filename: str, output_filename):
    with open(filename, "r") as f:
        lines = f.read().splitlines()
        begin_index = find_row_by_start(lines, "::: Pathway Inference Report written")
        end_index = find_row_by_start(lines, "List of pathways pruned")
        lisp_rows = lines[begin_index + 2 : end_index]
        lisp_code = " ".join(lisp_rows)
        print(lisp_code)


def direct_report_inference_description(
    inference_description_file: str, output_filename: str
):
    with open(inference_description_file, "r") as input_file:
        with open(output_filename, "w") as output_file:
            lisp_code = input_file.read().replace("\n", " ")
            for lisp_object in first_order_lisp_objects(lisp_code):
                print(lisp_object)
                pathway: Optional[str] = None
                explanation: Optional[str] = None
                keep: Optional[bool] = None
                i = 0
                while i < len(lisp_object):
                    if lisp_object[i] == ":PATHWAY":
                        pathway = lisp_object[i + 1]
                    elif lisp_object[i] == ":EXPLANATION-CODE":
                        explanation = lisp_object[i + 2]
                    elif lisp_object[i] == ":KEEP?":
                        keep = True if lisp_object[i + 1] == "T" else False
                    i += 1
                if pathway is not None and explanation is not None and keep is not None:
                    output_file.write(
                        "\t".join([pathway, str(keep), explanation]) + "\n"
                    )


def main():
    inference_description_file = sys.argv[
        1
    ]  # "./tmp/pathologic_output/pathways1000_frac100_seed1/tmp.sortion.ptools_local.SWAImfavJd/ptools-local/pgdbs/user/testcyc/1.0/reports/pwy-inference-description.data"
    direct_report_inference_description(
        inference_description_file, sys.argv[2]
    )  # "./tmp/pathways1000_frac100_seed1.pathologic_report.tsv")


if __name__ == "__main__":
    main()
