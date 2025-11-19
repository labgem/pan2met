import pytest
import ouisparql
from pangenome2panmetabolome.io.sparql import SPARQL

@pytest.fixture
def sparql_wrapper():
    sparql = SPARQL().wrapper
    return sparql

def test_simple_sparql_query():
    sparql = SPARQL().wrapper
    sparql.setQuery(
        """
        PREFIX bp: <http://www.biopax.org/release/biopax-level3.owl#>
        SELECT DISTINCT ?predicate
        WHERE {
             ?subject ?predicate ?object
        }
        ORDER BY ?predicate
        LIMIT 3
        """
    )
    ret = sparql.queryAndConvert()
    for r in ret["results"]["bindings"]:
        print(r)

def test_ouisparql_sparql_query(sparql_wrapper):
    queries = ouisparql.from_path("src/sparql/queries.rq", "sparql_wrapper")
    queries.get_all_predicate(sparql_wrapper)
    
