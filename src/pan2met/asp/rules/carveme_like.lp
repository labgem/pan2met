% CarveMe (Machado et al 2018) implements universal metabolic model carving into organism specific models
% with mixed integer linear programming
% Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. DOI: 10.1093/nar/gky537
%
% This ASP program does not attempt to closely reproduce the integer linear programming carving of CarveMe,
% but takes inspiration of this approach.

#maximize { N: gapless_pathway_count(N) }.
gapless_pathway_count(N) :- N=#count { gapless_pathway(Pathway) }.
