#!/usr/bin/env bash

sqlite3 mock_metabiantes.db < ../../submodules/metabiantes/sql/create_schema.sql
sqlite3 mock_metabiantes.db < mock_metabiantes_dump.sql
