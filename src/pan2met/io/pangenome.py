from typing import Dict, Set, Literal
from pathlib import Path
from collections import defaultdict


def read_pangenome_rtab(rtab_filename: Path) -> Dict[str, Set[str]]:
    """
    Read a .Rtab file
    :param rtab_filename: path to the .Rtab file
    :return: a dictionary mapping each strain identifier to the set of gene family identifiers present in the strain
    """
    with open(rtab_filename, "r") as rtab_file:
        header = rtab_file.readline().rstrip()
        strains = header.split()[1:]
        strain_to_families = defaultdict(set)
        for strain in strains:
            strain_to_families[strain] = set()
            for line in rtab_file:
                i = 0
                for field in line.rstrip().split("\t"):
                    if i == 0:
                        family_id = field
                    elif field == "1":
                        strain_to_families[strains[i - 1]].add(family_id)
                    i += 1
    return dict(strain_to_families)


def read_partition_files(
    persistent_filename: Path, shell_filename: Path, cloud_filename: Path
) -> Dict[Literal["persistent", "shell", "cloud"], Set[str]]:
    """
    Read persistent.txt, shell.txt and cloud.txt files, with lists of gene families identifiers.
    :param persistent_filename: Path to a text file listing gene families of the persistent partition
    :param shell_filename: Path to a text file listing gene families of the shell partition
    :param cloud_filename: Path to a text file listing gene families of the cloud partition
    :return: A dictionnary mapping the name of the partition ("persistent", "shell", "cloud") to the set of gene families identifiers.
    """
    partitions_with_families = dict()
    with open(persistent_filename) as f:
        partitions_with_families["persistent"] = {line.rstrip().upper() for line in f}
    with open(shell_filename) as f:
        partitions_with_families["shell"] = {line.rstrip().upper() for line in f}
    with open(cloud_filename) as f:
        partitions_with_families["cloud"] = {line.rstrip().upper() for line in f}
    return partitions_with_families


def read_functional_module_file(module_filename: Path) -> Dict[str, Set[str]]:
    """
    Read a pangenome module file.
    :param module_filename: Path to the module TSV file with two columns (module identifier, gene families identifier)
    :return: a dictionnary mapping the name of the module to the set of gene families identifiers.
    """
    modules_with_families = dict()
    with open(module_filename, "r") as module_file:
        _header = module_file.readline().rstrip()
        for line in module_file:
            (module_id, fam_id) = line.rstrip().split("\t")
            if module_id not in modules_with_families:
                modules_with_families[module_id] = set()
            modules_with_families[module_id].add(fam_id.upper())
    return modules_with_families
