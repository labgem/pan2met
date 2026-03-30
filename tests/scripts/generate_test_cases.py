import random
import os
import shutil

import pan2met.io.metabiantes_kb.kb

TEST_FOLDER = "./tests/cases/generated2/"

kb = pan2met.io.metabiantes_kb.kb.MetabiantesKnowledgeBase()

pathways = kb.pathways()

ecoli_tax_id = 562
thermococcales_tax_id = 2258
human_tax_id = 9606

tax_ids = [ecoli_tax_id, thermococcales_tax_id, human_tax_id]
test_case_sizes = [1, 5, 10, 20, 50, 1000, 2000, 3000]
pathway_fractions = [0.1, 0.35, 0.5, 0.75, 1.0]
seeds = list(range(5))


def random_protein_id() -> str:
    return "RAND_" + "".join(
        map(str, list(map(lambda _: random.randint(0, 10), range(10))))
    )


def write_pathologic_file(file, reactions: list[str]):
    for reaction in reactions:
        protein_id = random_protein_id()
        file.writelines(
            [
                f"ID\t{protein_id}\n",
                f"NAME\t{protein_id}\n",
                "STARTBASE\t1\n",
                "ENDBASE\t9\n",
                "PRODUCT-TYPE\tP\n",
                f"METACYC\t{reaction}\n",
                "//\n",
            ]
        )


def write_organism_params_file(file, tax_id):
    file.writelines(
        [
            "ID      test\n",
            "STORAGE FILE\n",
            "NAME    test\n",
            "ABBREV-NAME test\n",
            "STRAIN  test\n",
            "CREATE? T\n",
            f"DOMAIN  TAX-{tax_id}\n",
            "RANK    species\n",
            "AUTHOR  pan2met developer\n",
        ]
    )


os.makedirs(TEST_FOLDER, exist_ok=True)
with open(os.path.join(TEST_FOLDER, "taxon_id.tsv"), "w") as taxon_id_file:
    taxon_id_file.write("species\ttaxon_id\n")
    for seed in seeds:
        for tax_id in tax_ids:
            for test_case_size in test_case_sizes:
                for pathway_fraction in pathway_fractions:
                    organism_id = f"taxid{tax_id}_pathways{test_case_size}_frac{int(100 * pathway_fraction)}_seed{seed}"
                    test_folder = os.path.join(TEST_FOLDER, organism_id)
                    taxon_id_file.write(f"{organism_id}\t1\n")
                    os.makedirs(test_folder, exist_ok=True)
                    shutil.copy2(
                        "tests/pathologic/pathologic_template/genetic-elements.dat",
                        test_folder,
                    )
                    # shutil.copy2(
                    #    "tests/pathologic/pathologic_template/organism-params.dat", test_folder
                    # )
                    with open(
                        os.path.join(test_folder, "organism-params.dat"), "w"
                    ) as organism_params_file:
                        write_organism_params_file(organism_params_file, tax_id)
                    shutil.copy2(
                        "tests/pathologic/pathologic_template/test.fasta", test_folder
                    )
                    test_pathways = random.choices(pathways, k=test_case_size)
                    with open(
                        os.path.join(test_folder, "reactome"), "w"
                    ) as reactome_file:
                        with open(os.path.join(test_folder, "test.pf"), "w") as file:
                            for pathway in test_pathways:
                                reactions = kb.reactions_of_pathway(pathway)
                                n_reactions: int = int(
                                    len(reactions) * pathway_fraction
                                )
                                test_reactions = random.choices(
                                    reactions, k=n_reactions
                                )
                                write_pathologic_file(file, test_reactions)
                                reactome_file.writelines(
                                    reaction + "\n" for reaction in test_reactions
                                )
