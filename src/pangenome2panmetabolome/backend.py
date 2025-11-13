#!/usr/bin/env python3

import os
from abc import ABC, abstractmethod

import anosql
import sqlite3

SQL_QUERY_FILE_PATH = os.path.join(os.path.dirname(__file__), '../sql/queries.sql')

class PathwayKnowledge(ABC):
    @abstractmethod
    def list_pathways(self) -> list[str]:
        pass

    def reactions_of_pathways(self) -> list[str]:
        pass


class SQLiteBioPAXPathwayKnowledge(PathwayKnowledge):

    def __init__(self, database: str):
        self.connection = sqlite3.connect(database)
        self.queries = anosql.from_path(SQL_QUERY_FILE_PATH, "sqlite3")
         

    def list_pathways(self) -> list[str]:
        return self.queries.get_all_pathways(self.connection)
        

    def reactions_of_pathway(self, pathway: str) -> list[str]:
        return self.queries.get_biochemical_reactions_of_pathway(self.connection, pathway=pathway)

