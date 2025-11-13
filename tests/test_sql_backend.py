from pangenome2panmetabolome.backend import PathwayKnowledge, SQLiteBioPAXPathwayKnowledge 

def test_listing_all_pathways():

    kb: PathwayKnowledge = SQLiteBioPAXPathwayKnowledge("metacyc-biopax.db")

    print(kb.list_pathways())

    print(kb.reactions_of_pathway("bp:Pathway683791"))
