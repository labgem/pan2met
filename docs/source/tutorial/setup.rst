Setup pan2met on your machine
=============================

**pan2met** is a Python library with some python dependencies. To install **pan2met** on your machine, you will need a working Python installation.

You have multiple options to choose regarding the installation of **pan2met**.
First, you will have to decide on which tool **pan2met** will rely on to get reference data on gene, reaction and metabolic pathways.

One option is `metabiantes <https://gitlab.com/sortion/metabiantes>`_, the other is `padmet <https://github.com/AuReMe/padmet>`_.

metabiantes is a relational database schema mimicking MetaCyc object data model.
It is provided with a Lisp script that enables an import of MetaCyc pathway knowledge base into a SQL database (currently either postgreSQL or SQLite).
Confer to the `metabiantes repository <https://gitlab.com/sortion/metabiantes>`_, for instruction on how to build a database with metabiantes schema and import data from MetaCyc with the included Lisp script.

padmet is a Python library and custom plain text metabolic database format developped by IRISA/Dyliss team in Rennes, France, also mimicking the MetaCyc object data model.
Refer to `padmet documentation <https://padmet.readthedocs.io/`_ for instructions on how to construct a .padmet file with reference metabolic pathway knowledge.

PADMet with MetaCyc
~~~~~~~~~~~~~~~~~~~

In summary, to create a file `metacyc.padmet` with data from a local pathway-tools installation of MetaCyc, you will proceed in two steps:
First, you will need to dump pathway-tools data into flat files.
To do so, launch `pathway-tools` in Lisp API mode:

.. code:: bash

    pathway-tools -lisp

Then, launch the export to flat file command:

.. code:: lisp

    EC(0): (select-organism :org-id 'meta)
    EC(1): (create-flat-files-for-current-kb)

Secondly, use `padmet pgdb_to_padmet` command to create the padmet file:

.. code:: bash

    padmet pgdb_to_padmet \ --pgdb=~/.local/share/pathway-tools/aic-export/pgdbs/biocyc/metacyc/29.5/data \
    --extract-gene \
    --source="metacyc-29.5" \
    --output=metacyc.padmet

You will probably have to change the path to the pgdb `data` folder to fit to your local installation configuration, and also change the version in the example path '29.5' to your version of MetaCyc.

metabiantes with MetaCyc
~~~~~~~~~~~~~~~~~~~~~~~~

We refer the reader to the `metabiantes repository <https://gitlab.com/sortion/metabiantes>`_ for instruction on how to build a database with metabiantes schema and import data from MetaCyc with the included Lisp script.
In a nutshell, you will have to:

1. Make sure you have a working installation of PostgreSQL or SQLite on your machine.

2. Make sure you have a local installation of PathwayTools and MetaCyc database on your machine, as **metabiantes** data loader relies on the local Lisp PathwayTools API to extract data from PathwayTools.

3. Create a database with the metabiantes schema, and import MetaCyc data into it with the provided Lisp script.

    .. code:: bash

        sudo -u postgres psql

    .. code:: sql

        CREATE DATABASE metabiantes OWNER <user>;

4. Clone the **metabiantes** git repository

    .. code:: bash

        git clone https://gitlab.com/sortion/metabiantes.git
        cd metabiantes/metabiantes/loader


5. Dump the MetaCyc data with the provided Lisp script, or the wrapper shell script:

    .. code:: bash

        bash metabiantes.sh "dump.sql" "meta"


6. Create the metabiantes database schema

   .. code:: bash

        psql -d metabiantes < ../sql/create_schema_pg.sql

7. Load the metabiantes `dump.sql` generated file into the Postgres database with:

   .. code:: bash

        psql -d metabiantes < ./dump.sql


Importing **metabiantes**' '`dump.sql` into a SQL database can take a while, as **metabiantes** uses subqueries in `INSERT INTO` SQL instructions to favor foreign references when possible, in a way that may possibly be not that effective...

Hopefully you will only need to do that once and be done for a while. We advise you to create a 'real' dump once the database is fully created, to allow for a faster load of the data the next time you need this kind of database (e.g, on another of your machine), instead of doing again the slow import of the data from the hand crafted dump generated with **metabiantes**.
