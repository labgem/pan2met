"""
Basic comparison report between PathoLogic metabolome output and the output of pangenome2panmetabolome

"""

import datetime
import subprocess
import os

from pangenome2panmetabolome.utils import read_list


def get_git_revision():
    return (
        subprocess.check_output(["git", "rev-parse", "--short", "HEAD"])
        .decode("ascii")
        .strip()
    )


git_revision = get_git_revision()

TEST_CASES_FOLDER = "tests/cases/generated/"
PATHOLOGIC_FILENAME = "pathologic.pathways.list"
PANGENOME2PANMETABOLOME_FILENAME = "pangenome2panmetabolome.pathways.list"

test_cases = os.listdir(TEST_CASES_FOLDER)


def report(
    pathologic_metabolome: set[str], pangenome2panmetabolome_metabolome: set[str]
) -> set:
    union = pathologic_metabolome | pangenome2panmetabolome_metabolome
    unique_pathologic = pathologic_metabolome - pangenome2panmetabolome_metabolome
    unique_pangenome2panmetabolome = (
        pangenome2panmetabolome_metabolome - pathologic_metabolome
    )
    intersection = pathologic_metabolome & pangenome2panmetabolome_metabolome
    return [
        len(pathologic_metabolome),
        len(pangenome2panmetabolome_metabolome),
        len(union),
        len(intersection),
        len(unique_pathologic),
        len(unique_pangenome2panmetabolome),
    ]


def main():
    with open(
        f"report_{datetime.datetime.today().strftime('%Y-%m-%d')}_gitrev{git_revision}.tsv",
        "w",
    ) as report_file:
        report_file.write(
            "\t".join(
                [
                    "test_case",
                    "pathologic",
                    "pangenome2panmetabolome",
                    "union",
                    "intersection",
                    "unique_pathologic",
                    "unique_pangenome2panmetabolome",
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
                pangenome2panmetabolome_metabolome = set(
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
                                    pangenome2panmetabolome_metabolome,
                                ),
                            )
                        )
                    )
                    + "\n"
                )


if __name__ == "__main__":
    main()
