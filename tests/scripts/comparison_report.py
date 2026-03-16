"""
Basic comparison report between PathoLogic metabolome output and the output of pan2met
"""

import datetime
import subprocess
import os

from pan2met.utils import read_list


def get_git_revision():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )


git_revision = get_git_revision()

TEST_CASES_FOLDER = "tests/cases/generated/"
PATHOLOGIC_FILENAME = "pathologic.pathways.list"
PANGENOME2PANMETABOLOME_FILENAME = "pan2met.pathways.list"

test_cases = os.listdir(TEST_CASES_FOLDER)


def report(pathologic_metabolome: set[str], pan2met_metabolome: set[str]) -> list[int]:
    union = pathologic_metabolome | pan2met_metabolome
    unique_pathologic = pathologic_metabolome - pan2met_metabolome
    unique_pan2met = pan2met_metabolome - pathologic_metabolome
    intersection = pathologic_metabolome & pan2met_metabolome
    return [
        len(pathologic_metabolome),
        len(pan2met_metabolome),
        len(union),
        len(intersection),
        len(unique_pathologic),
        len(unique_pan2met),
    ]


def main():
    with open(
        f"./tmp/report_{datetime.datetime.today().strftime('%Y-%m-%d_%H-%M')}_gitrev:{git_revision}.tsv",
        "w",
    ) as report_file:
        report_file.write(
            "\t".join(
                [
                    "test_case",
                    "pathologic",
                    "pan2met",
                    "union",
                    "intersection",
                    "unique_pathologic",
                    "unique_pan2met",
                ]
            )
        )
        for test_case in test_cases:
            if os.path.isdir(os.path.join(TEST_CASES_FOLDER, test_case)):
                pathologic_metabolome = set(
                    read_list(
                        os.path.join(TEST_CASES_FOLDER, test_case, PATHOLOGIC_FILENAME)
                    )
                )
                pan2met_metabolome = set(
                    read_list(
                        os.path.join(
                            TEST_CASES_FOLDER,
                            test_case,
                            PANGENOME2PANMETABOLOME_FILENAME,
                        )
                    )
                )
                report_file.write(
                    "\t".join(
                        [test_case]
                        + list(
                            map(
                                str,
                                report(
                                    pathologic_metabolome,
                                    pan2met_metabolome,
                                ),
                            )
                        )
                    )
                    + "\n"
                )


if __name__ == "__main__":
    main()
