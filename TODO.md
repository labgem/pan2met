# Tasks

## Knowledge base input methods

- [ ] QLever repository of BioPAX
- [ ] Plain folder of BioPAX
  - QLever backend or pax2graphml-like backend?
- [ ] Plain folder of SBML

## Input format

- [x] Reaction / Monomer table
- [ ] Reaction list: in this case, we will not be able to use the graph topology.


## Export format

- [ ] PADMet
- [ ] SBML
- [ ] BioPAX
- [x] Plain list of pathway identifiers
- [ ] Plain list of reaction identifiers (the input reactions + the inferred ones?)

## Inference rules

- [ ] Improve the way inference rules are handled
    - [ ] Use a priority queue as in pan2met-rs instead of dumb if then else
