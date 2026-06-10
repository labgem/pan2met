# import itertools

# import pytest
# import ouisparql

# from pan2met.io.sparql import SPARQL


# @pytest.fixture
# def sparql_wrapper():
#     sparql = SPARQL().wrapper
#     return sparql


# def test_simple_sparql_query():
#     sparql = SPARQL().wrapper
#     sparql.setQuery(
#         """
#         SELECT DISTINCT ?predicate
#         WHERE {
#              ?subject ?predicate ?object
#         }
#         ORDER BY ?predicate
#         LIMIT 3
#         """
#     )
#     ret = sparql.queryAndConvert()
#     for r in ret["results"]["bindings"]:
#         print(r)


# def test_ouisparql_sparql_query(sparql_wrapper):
#     queries = ouisparql.from_path("src/sparql/queries.rq", "sparql_wrapper")
#     queries.get_all_predicate(sparql_wrapper)


# def test_get_pathway(sparql_wrapper):
#     queries = ouisparql.from_path("src/sparql/queries.rq", "sparql_wrapper")
#     res = queries.get_metacyc_pathways(sparql_wrapper)
#     for r in itertools.islice(res, 5):
#         print(r)
