This namespace contains movers designed to work with cyclic peptides.  These 
peptides are somewhat difficult to work with because they defy some of the 
assumptions made by rosetta.  For example, fold trees inherently assume linear 
peptides and the PDB loader automatically inserts N- and C-termini.

The classes defined in this namespace try to get around these conventions.  
Note that I am not an experienced rosetta developer, and consequently some of 
these classes are pretty hacky.  In particular, I think someone with a better 
understanding of how AtomTrees work could make this code much more robust.
