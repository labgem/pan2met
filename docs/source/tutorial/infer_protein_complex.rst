Infer the set of constructible proteic complex from the set of protein monomers
===============================================================================

`pan2met` includes a subcommand to infer the list of protein complex that could be built based on a set of protein monomer components.
To do so, it relies on the on the composition of the protein complex in the reference knowledge base (both metabiantes and padmet uses MetaCyc data as a reference).

For now, complex information is not available in PADMet files, so this functionality will work only with a metabiantes knowledge base backend.

Make sure you have a metabiantes SQLite database available and that you set metabiantes as the source of knowledge in the pan2met configuration file you use (see `<setup>`__)
