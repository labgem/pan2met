#!/usr/bin/env python3

"""
Export a TSV from BioPAX OWL reference
for reactions, containing the following columns:
- Reaction ID
- Reaction Name
- Spontaneous ?
- Orphan ?

"""

import csv

from typing import Iterable
import argparse

import pybiopax
from pybiopax.biopax import BioPaxModel
from pybiopax.biopax.interaction import (Catalysis, Conversion, BiochemicalReaction)
from pybiopax.biopax.physical_entity import (Protein, Complex)




def gibbs_free_energy(reaction: Catalysis) -> float:
    r"""
    Gibbs free energy of a biochemical reaction $\Delta G^{\degree}$
    """
    assert len(reaction.delta_g) == 1, "There should be one value of Gibbs free energy ΔG°"
    return reaction.delta_g[0].delta_g_prime0

def iter_catalysis(model: BioPaxModel, reaction: BiochemicalReaction) -> Iterable[Catalysis]:
    for obj in model.get_objects_by_type(Catalysis):
        if (isinstance(obj.controlled, BiochemicalReaction)
            and obj.controlled.uid == reaction.uid):
            yield obj
            
def iter_enzymes(model: BioPaxModel, reaction: BiochemicalReaction) -> Iterable[Protein | Complex]:
    for catalysis in iter_catalysis(model, reaction):
        yield catalysis.controller
        

def is_orphan(model: BioPaxModel, reaction: BiochemicalReaction) -> bool:
    """
    Assume an enzymatic reaction is orphan, if no evidence is found for gene coding an enzyme catalyzing this reaction
    """
    enzymes = iter_enzymes(model, reaction)
    return next(enzymes, None) is None


def reaction_table(model: BioPaxModel):
    reactions = filter(lambda obj: issubclass(obj.__class__, Conversion), model.objects.values())
    for reaction in reactions:
        uid = reaction.uid
        name = reaction.name[0]
        spontaneous = reaction.spontaneous
        orphan = None
        ec_number = None
        if issubclass(reaction.__class__, BiochemicalReaction):
            orphan = is_orphan(model, reaction)
            if isinstance(reaction.e_c_number, list):
                ec_number = reaction.e_c_number[0] 
        yield [uid, name, spontaneous, orphan, ec_number]


def parse_arguments():
    parser = argparse.ArgumentParser(description="BioPAX to reaction list")
    parser.add_argument('-i', '--input',
                        help="Input BioPAX file in owl format or owl.gz format")
    parser.add_argument('-o', '--output',
                        help="Output TSV with columns reaction id, reaction name, spontaneous?, orphan? and ec-number")
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    model = pybiopax.model_from_owl_file(args.input)    
   
    fieldnames = ["uid", "name", "spontaneous?", "orphan?", "ec-number"]
   
    with open(args.output, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(fieldnames)
        for row in reaction_table(model):
            writer.writerow(row)
 

if __name__ == "__main__":
    main()
