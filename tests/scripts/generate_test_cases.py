import random
import os
import shutil

import pangenome2panmetabolome.io.metabiantes.kb

SEED = 1
TAX_ID = 562

kb = pangenome2panmetabolome.io.metabiantes.kb.MetabiantesKnowledgeBase()

pathways = kb.pathways()

test_case_sizes = [1, 5, 10, 20, 100, 1000]
pathway_fractions = [0.35, 0.5, 0.75, 1.0]


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


os.makedirs("./tests/cases/generated/", exist_ok=True)
with open("./tests/cases/generated/taxon_id.tsv", "w") as taxon_id_file:
    taxon_id_file.write("species\ttaxon_id\n")
    for test_case_size in test_case_sizes:
        for pathway_fraction in pathway_fractions:
            random.seed(SEED)
            organism_id = (
                f"pathways{test_case_size}_frac{int(100 * pathway_fraction)}_seed{SEED}"
            )
            test_folder = f"./tests/cases/generated/{organism_id}"
            taxon_id_file.write(f"{organism_id}\t1\n")
            os.makedirs(test_folder, exist_ok=True)
            shutil.copy2(
                "tests/pathologic/pathologic_template/genetic-elements.dat", test_folder
            )
            shutil.copy2(
                "tests/pathologic/pathologic_template/organism-params.dat", test_folder
            )
            shutil.copy2("tests/pathologic/pathologic_template/test.fasta", test_folder)
            test_pathways = random.choices(pathways, k=test_case_size)
            with open(os.path.join(test_folder, "reactome"), "w") as reactome_file:
                with open(os.path.join(test_folder, "test.pf"), "w") as file:
                    for pathway in test_pathways:
                        reactions = kb.reactions_of_pathway(pathway)
                        n_reactions: int = int(len(reactions) * pathway_fraction)
                        test_reactions = random.choices(reactions, k=n_reactions)
                        write_pathologic_file(file, test_reactions)
                        reactome_file.writelines(
                            reaction + "\n" for reaction in test_reactions
                        )
