(dp0
S'ddg_monomer'
p1
(dp2
S'LIMITATIONS'
p3
S'The high-resolution protocol is designed to limit the amount of conformational flexibility allowed to protein backbones as the amount of noise that seeps into the protocol from increased flexibility tends to drown out the signal that might be gained by searching a larger region of conformation space. The result is that this protocol is probably not well suited to model several mutations simultaneously where backbone motion would be expected.\n'
p4
sS'INPUTS'
p5
S'All PDBs should be renumbered so that their first residue is residue 1 and number consecutively so that, if there are missing residues in the structure (due maybe to missing density) that these residues are simply skipped in the residue numbering. The numbering of all residues in both the distance-restraint file and the mutation-list file should follow this numbering.\\nThere are two main input files: a) Restraint / constraint files for calpha atom pairs, and b) Mutation-list files, describing which set of point mutations to entertain for an input backbone.\\na) Restraint / constraint files for calpha atom pairs  The constraint file should list all calpha atom pairs within 9 Angstroms, giving the distance measured from the original (non-pre-minimized) experimental structure. See RosettaCommons (Maybe a part of the toolkit).\n'
p6
sS'DESCRIPTION'
p7
S'The purpose of this application is to predict the change in stability (the ddG) of a monomeric protein induced by a point mutation. The application takes as input the crystal structure of the wild-type (which must be first pre-minimized), and generates a structural model of the point-mutant. The ddG is given by the difference in rosetta energy between the wild-type structure and the point mutant structure. More precisely, 50 models each of the wild-type and mutant structures should be generated, and the most accurate ddG is taken as the difference between the mean of the top-3-scoring wild type structures and the top-3-scoring point-mutant structures.\\n\\n\n'
p8
sS'ALGORITHM'
p9
S'There are two main ways that this application should be used: a high-resolution and a low-resolution way. They are nearly as accurate as each other with corellation coefficients of 0.69 and 0.68 on a set of 1210 mutations. These are described by the protocols on rows 16 and 3 of [Kellogg,2011], respectively.\\nA) High Resolution Protocol:\\nThis protocol allows a small degree of backbone conformational freedom. The protocol optimizes both the initial input structure for the wild-type and the generated structure for the point mutant in the same way, for the same number of iterations (recommended 50). It begins by optimizing the rotamers at all residues in the protein using Rosetta\'s standard side-chain optimization module (the packer). It follows this initial side-chain optimization with three rounds of gradient-based minimiztion, where the repulsive component of the Lennard-Jones (van der Waals) term is downweighted in the first iteration (10% of its regular strength), weighted at an intermediate value in the second iteration (33% of its regular strength), and weighted at its standard value in the third iteration. This repacking followed by minimiztion is run several times, always starting from the same structure. Scores and optionally PDBs are written out.\\nDistance Restraints: The high-resolution protocol relies on the use of Calpha-Calpha distance restraints as part of the optimization to prevent the backbone from moving too far from the starting conformation. These distance restraints may be generated externally before the protocol may be run. If not specified, constraints will be automatically generated based on the input structure, but the results obtained in the published paper utilized constraints based on the high-resolution crystal structure. The constraints used in the generation of data for row 16 in [Kellogg2011] were given as distance restraints between all Calpha pairs within 9 Angstroms of each other in the wild type structure; for each harmonic restraint, the ideal value for the restraint was taken as the distance in the original crystal structure (not the pre-minimized structure which should be given as input) and the standard-deviation on the harmonic constraint was set to 0.5 Angstroms. For example, the distance restraint between the c-alpha of residue 1 and the c-alpha of residue 2 of the PDB 1hz6 is described in the input constraint file by the line "AtomPair CA 2 CA 1 HARMONIC 3.79007 0.5".\\n\\ndistance restraints can be generated through the use of this shell script:\\n/convert_to_cst_file.sh mincst.log > input.cst\\nthis shell script simply takes the output of the minimization log (from pre-minimization of the input structure) and converts it to the appropriate constraint file format. ( see below to obtain minimization log \\n\\nB) Low Resolution Protocol:\\nThis protocol only allows sidechain conformational flexibility. It optimizes the rotamers for the residues in the neighborhood of the mutation; those with CBeta atoms within 8 Angstroms of the CBeta atom of the mutated residue (or Calpha for glycine). The same set of residues is optimized for both the wild-type and mutant structures. The optimization is performed a recommended 50 times for both the wildtype and point-mutant, and the scores (and optional PDBs) are written out.\n'
p10
sS'AUTHOR'
p11
S'Andrew Leaver-Fay and Elizabeth Kellogg\n'
p12
sS'OUTPUTS'
p13
S"The output of the ddg protocol is a 'ddg_predictions.out' which contains, for each mutation, the total predicted ddg and a breakdown of all the score components which contribute to that total. Furthermore, output structures are dumped either in silent-file or pdb format. If silent-files are output, the following naming convention is used. All wild-type structures are dumped into wt_<WT_AA><RESIDUE_NUM><MUTANT_AA>.out The reason wild-type structures are always dumped is because if local optimization around the site of mutation is being done, the wild-type structures can potentially be different from one another due to different constraint definitions or different packing definitions. Mutant structures follow a similar convention: mut_<WT_AA><RESIDUE_NUM><MUTANT_AA>.out For example, if you made a A to Q mutation at residue 123, you would see two silent-files as output: wt_A123Q.out and mut_A123Q.out this is done regardless of the protocol used.\n"
p14
sS'ANALYSIS'
p15
S'No Comments\n'
p16
sS'REFERENCES'
p17
S'The new algorithm for performing limited relaxation of the backbone was published in\\nE. Kellogg, A. Leaver-Fay, and D. Baker, (2011) "Role of conformational sampling in computing mutation-induced changes in protein structure and stability", Proteins: Structure, Function, and Bioinformatics. V 79, pp 830--838.\\n\\nThe older, fixed-backbone, soft-repulsive scorefunction algorithm (analogous to that described in row 4 of [Kellogg 2011] but with weights trained towards recapitulating alanine-scanning mutation experiments. weights are in minirosetta_database/scoring/weights/ddg_monomer.wts) was published in:\\nKortemme et al. (2002) "A simple physical model for binding energy hot spots in protein-protein complexes", PNAS 22, 14116-21\n'
p18
sS'EXAMPLES'
p19
S'The ddg_monomer application lives in src/apps/public/ddg/ddg_monomer.cc (This file had previously been named "fix_bb_monomer_ddg.cc", but has been renamed since it now moves the backbone). This file houses the main() function. The central subroutines invoked by this file live in the ddGMover class defined in src/protocols/moves/ddGMover.hh and src/protocols/moves/ddGMover.cc.\\nA helper application, minimize_with_cst, lives in src/apps/public/ddg/minimize_with_cst.cc.\\nA helper script for generating one of the input files needed by the ddg_monomer application lives in src/apps/public/ddg/convert_to_cst_file.sh.\\nAn integration test for this application lives in test/integration/tests/ddG_of_mutation/. The test in this directory runs a shortened trajectory for predicting the wild-type and mutant energies. To turn this into a production-run example, set the value for the "-ddg:iterations" flag given in the file "test/integration/tests/ddG_of_mutations/flags" to 50.\n'
p20
sS'TIPS'
p21
S'Make sure to use exactly as described in protocol 13 of the main paper referenced above.\n'
p22
sS'METADATA'
p23
S'The documentation was last updated on 4/7/2011, by Andrew Leaver-Fay. Questions about this documentation should be directed to David Baker: dabaker@u.washington.edu.\n'
p24
ssS'relax'
p25
(dp26
S'LIMITATIONS'
p27
S'Does a great job at minimizing the structure.  button_ be careful.\n'
p28
sS'INPUTS'
p29
S'Relax takes as input one or more structures in silent or PDB format. Or Use the PyRosetta Toolkit.  All JD2 options apply (see JD2 Documentation ofr more details)\n'
p30
sS'DESCRIPTION'
p31
S'The "relax" application in Rosetta carries out the task of simple structural refinement of fullatom Rosetta models. It can also read Centroid models, in which case it will convert the model into a fullatom model and pack the sidechains. Relax does not carry out any extensive refinement and only searches the local conformational space neighbourhood.\\n\\n For virtually all situations it should be sufficient to use either -relax:quick or -relax:thorough and not worry about all the options.\n'
p32
sS'ALGORITHM'
p33
S'The fast relax modes work by running many sidechain repack and minimisation cycles ramping up and down the repulsive weight of the forcefield. This sort of "pulses" the structure in successive collapse and expansion phases which "wiggles" the atoms into a low energy state. No explicit moves are done (this was found not to be useful as most moves are rejected and dont help lowering the nergy). Not that despite that fact, the structure can change up to 2-3 A from the starting conformation during the minimisation cycles.\\n\\nFastRelax is a more modern version of the initial fast relax algotihm which is more flexible and allows definition of a specific script of relax events (how exactly the repack and minimisation cycles are interlaced and what paramters they should take).\n'
p34
sS'AUTHOR'
p35
S'Mike Tyka\n'
p36
sS'OUTPUTS'
p37
S'Relax outputs one or more structures in silent or PDB format. All JD2 output options apply (see JD2 Documentation ofr more details) As always -nstruct regulates the number of outputs per input structure.\n'
p38
sS'ANALYSIS'
p39
S'Typically, one would either cluster similar conformations or simply take the lowest energy score.  You can do this within the toolkit.\n'
p40
sS'EXAMPLES'
p41
S'The example with input files can be found in demo/relax/ Just execute the run.sh script.\n'
p42
sS'TIPS'
p43
S'Relax will increase the rmsd of your protein to a high degree.  Use only if you really want to do this.  A better and quicker way to run this application is to use the built in relax applications for loops, domains, subunits, and proteins within the full PyRosetta Toolkit.  See user guide for more details.\n'
p44
sS'METADATA'
p45
S'This document was edited May 25th 2010 by Mike Tyka. This application in mini was created and documented by Mike Tyka, et al.\n'
p46
ssS'fragment_picker'
p47
(dp48
S'INPUTS'
p49
S'vall: protein structures database, your fragments come from there. http://www.bioshell.pl/rosetta-related/vall.apr24.2008.extended.gz - Mandatory\\n.wghts: defines scoring system for fragment selection - Mandatory\\n.fasta Amino Acid sequence - Mandatory unless .chk file given.\\n.chk: sequence profile created with PSI-Blast with further modifications (pseudocounts added) use make_fragments.pl - any sequence profile - based score, e.g. ProfileScoreL1; mandatory file unless .fasta is given\\n.ss2:secondary structure prediction in PsiPred format - The easiest way is to run make_fragments.pl script. You may also try to run a secondary prediction software on your own and then convert the resulst to the proper format. A script convert_ss_predictions.py can turn TALOS, Juffo, Porter and SAM into ss2.\\n.cst:distance (or dihedral) constraints - Convert your data (distances or torsion angle values) into the proper format.\\n.tab:NMR experiment; examples can be downloaded from BMRB database\\n.pdb:Your structure\n'
p50
sS'REFERENCES'
p51
S'Gront D, Kulp DW, Vernon RM, Strauss CEM and Baker D, "Generalized fragment picking in Rosetta: design, protocols and applications", submitted to PLoS ONE\n'
p52
sS'DESCRIPTION'
p53
S'In brief, the program reads a database file (nicknamed vall), input query sequence or sequence profile and other files and produces fragment files for modeling with Rosetta. \n'
p54
sS'ALGORITHM'
p55
S'Detail of the algorithm are described in Gront D. et al paper. \n'
p56
ssS'backrub'
p57
(dp58
S'LIMITATIONS'
p59
S'The backrub application does not sample either the backbone or side chain of proline residues. As long as a proline residue is specified as flexible, it can be part of backrub segments, but cannot be either the start or end pivot of a segment.\\nOne of the primary differences between this implementation and the previous implementation is that backrub is now atom-centric rather than residue-centric. check_button_ck the references or RosettaCommons for more information on the change.\n'
p60
sS'INPUTS'
p61
S'The starting structures must be in PDB format and can be specified using the -s or -l options. A custom fold tree can be specified on a single line of the PDB file using the silent file format. An overview of that format can be found in the fold tree documentation.\\nSide chain sampling is controlled using the -resfile command line option. If no resfile is specified, then all side chains are made flexibile by default. Please see the resfile documentation for more information about how to create one. There are several things to note when using resfiles: First, because of current limitations with side chain sampling, proline resfiles are not sampled, even if specified in the resfile. Second, while it is possible to sample different amino acids using the backrub application, the fixed temperature Monte Carlo algorithm will bias the selected amino acids towards smaller residues such as alanine. Lastly, the residue numbering in resfiles is based on the residue number and chain letter from the PDB file, which is different from the -pivot_residues option.\\nSimple phi/psi backbone moves can also be enabled by specifying a MoveMap file using the -movemap option and giving a greater than zero value to the -sm_prob option. See the MoveMap documentation for more information.\n'
p62
sS'DESCRIPTION'
p63
S'This application is useful for creating ensembles of protein backbones, modeling protein flexibility, modeling mutations, and detailed refinement of backbone/side chain conformations.\n'
p64
sS'ALGORITHM'
p65
S'The backrub algorithm rotates local segments of the protein backbone as a rigid body about an axis defined by the starting and ending atoms of the segment. It was inspired by observations made by Davis et al (Structure 2006) of alternate side chain/backbone conformations in high resolution crystal structures. Atoms branching of the main chain at the pivot points (side chains, hydrogens, carbonyl oxygens), are updated to minimize the bond angle strain incurred. These moves are accepted or rejected using the Metropolis criterion.\\nIn addition to backrub moves, side chain conformations are sampled directly from the probability distributions described by the Dunbrack rotamer library, and not from a discrete set of chi angles, as is typically done by many side chain sampling algorithms. Side chain moves are also accepted or rejected using the Metropolis criterion.\n'
p66
sS'AUTHOR'
p67
S'Colin A. Smith\n'
p68
sS'OUTPUTS'
p69
S'For each starting structure, an output tag will be generated from the input file name, suffix, prefix, and user tags, if applicable. The backrub application generates two files, output_tag_0001_low.pdb and output_tag_0001_last.pdb. The four digit index is incremented up to the number of structures specified by -nstruct. The "low" file contains the lowest energy structure found during the Monte Carlo simulation. The "last" file contains the last accepted structure during the Monte Carlo simulation. If a custom fold tree was given in the input file, it will be appended to each of the output files.\n'
p70
sS'ANALYSIS'
p71
S'A useful form of post processing is to calculate the RMSD of the output structures to the input structure.\n'
p72
sS'REFERENCES'
p73
S'Smith CA, Kortemme T. Structure-Based Prediction of the Peptide Sequence Space Recognized by Natural and Synthetic PDZ Domains. J Mol Biol. http://dx.doi.org/10.1016/j.jmb.2010.07.032\\nSmith CA, Kortemme T. Backrub-like backbone simulation recapitulates natural protein conformational variability and improves mutant side-chain prediction. J Mol Biol. 2008 Jul 18;380(4):742-56. http://dx.doi.org/10.1016/j.jmb.2008.05.023\\nFriedland GD, Linares AJ, Smith CA, Kortemme T. A simple model of backbone flexibility improves modeling of side-chain conformational variability. J Mol Biol. 2008 Jul 18;380(4):757-74. http://dx.doi.org/10.1016/j.jmb.2008.05.006\\nFriedland GD, Lakomek NA, Griesinger C, Meiler J, Kortemme T. A correspondence between solution-state dynamics of an individual protein and the sequence and conformational diversity of its family. PLoS Comput Biol. 2009 May;5(5):e1000393. http://dx.doi.org/10.1371/journal.pcbi.1000393\\nDavis IW, Arendall WB 3rd, Richardson DC, Richardson JS. The backrub motion: how protein backbone shrugs when a sidechain dances. Structure. 2006 Feb;14(2):265-74. http://dx.doi.org/10.1016/j.str.2005.10.007\\nBetancourt MR. Efficient Monte Carlo trial moves for polypeptide simulations. J check_button_m Phys. 2005 Nov 1;123(17):174905. http://dx.doi.org/10.1063/1.2102896\n'
p74
sS'EXAMPLES'
p75
S'he code for the backrub application is in src/apps/public/backrub.cc. An integration test and demo is located in test/integration/tests/backrub. Backrub moves are made with the BackrubMover. Side chain moves are made with the SidechainMover. Backbone phi/psi moves are made with the SmallMover.\n'
p76
sS'TIPS'
p77
S'To date, typical backrub ensemble generation has used 10,000 Monte Carlo steps at a temperature of 0.6. At this temperature, many structures will unfold if the number of Monte Carlo steps is increased significantly. Many structures remain stable in extended simulations at a temperature of 0.3-0.4.\\nThe 10,000 step backrub simulations for a recent PDZ specificity prediction paper (Smith & Kortemme 2010) took an average of 110 seconds per simulation to generate a single structure. The simulations were each run on a single core of a heterogeneous cluster of 8 core Xeon workstations with E5345, E5430, and E5520 processors.\n'
p78
sS'METADATA'
p79
S'This document was last updated August 10, 2010, by Colin A. Smith. The corresponding principal investigator is Tanja Kortemme <kortemme@cgl.ucsf.edu>.\n'
p80
ssS'loopmodel'
p81
(dp82
S'INPUTS'
p83
S'Start pdbs: The template pdb file and must have real coordinates for all template residues plus the first and last residue of each loop region.  Loop file (Create using the GUI - check export):\n'
p84
sS'DESCRIPTION'
p85
S'Loop modeling is performed by two different algorithms CCD (Cyclic coordinate descent) and KIC (Kinematic closure). Here only CCD is described and the explanation of the latter can be found the corresponding KIC documentation. The goal of both algorithms is to explore the conformational space of the loop using a centroid representation of protein side-chains and explicit backbone representation, followed by a higher-resolution search using explicit representations of all atoms and hydrogen.\\n\\nThe centroid stage of loop-modeling generates loops by performing fragment insertions using Monte Carlo sampling, a score to reward closed chains, and CCD is used to close the loop at the end of the simulation. As the fragments are necessary for the sampling these have to be generated by fragment picker ( c.f. fragment picker documentation) or downloaded from the Robetta web-server.\n'
p86
sS'AUTHOR'
p87
S'Sinisa Bjelic and TJ Brunette\n'
p88
sS'ANALYSIS'
p89
S"For benchmarking purposes, creating a score vs rmsd plot across decoys and looking for near native 'energy funnels' is good way to test the performance of the protocols on a system, and can help to determine whether errors are due to scoring or sampling. For blind prediction and refinement, such plots can still be useful to look for convergence or multiple minima in the energy landscape. Decoys may also be pairwise-clustered to search for well-populated regions of conformational space that may represent alternative low-energy conformations. (from KIC loopclosure)\n"
p90
sS'REFERENCES'
p91
S'Qian, B., Raman, S., Das, R., Bradley, P., McCoy, A.J., Read, R.J. and Baker D. (2007). High resolution protein structure prediction and the crystallographic phase problem. Nature. manuscript accepted.\\nWang, C., Bradley, P. and Baker, D. (2007) Protein-protein docking with backbone flexibility. Journal of Molecular Biology, in press, DOI,http://dx.doi.org/10.1016/j.jmb.2007.07.050\n'
p92
sS'EXAMPLES'
p93
S'See mini/test/integration/tests/loop_modeling for an example loop relax run and input files.\n'
p94
sS'TIPS'
p95
S'For production runs, it is recommended to use the following flags. -loops::remodel quick_ccd -loops::refine refine_kic -loops::relax fastrelax -relax::fastrelax_repeats 8 -loops::extended and to generate at least 1000 models using -nstruct 1000.\\nquick_ccd can also remodel termini. To do this set the cutpoint in the loops file to be equal to the last residue in the chain. For example for a 80 residue protein, if you want to remodel the first 10 residues the loop file would have 1 10 10 0 0\\nquick_ccd does not require constraints, but using constraints from homologs or experimental data can produce more accurate results. Output consists of a pdb and a scorefile. The job concludes with the following command:\\nprotocols.looprelax: '
p96
sS'METADATA'
p97
S'This document was last updated August 11, 2010 by TJ Brunette & Sinisa Bjelic. The corresponding PIs for this application are David Baker <dbaker@u.washington.edu>.\n'
p98
ssS'cluster'
p99
(dp100
S'DESCRIPTION'
p101
S'The "cluster" application in Rosetta carries out a simple clustering of structures (either PDB or silent file format). \n'
p102
sS'ALGORITHM'
p103
S'The algorithm is based on one of Phil Bradley\'s old programs (silent_cluster_c). Starting with a subset of structures (the first 400 structures) the algorithm finds the structure with the largest number of neighbors within the cluster radius and creates a first cluster with that structure as the clsuter center and the neighbors part of and claimed by the cluster. The structures are removed from the pool of "unclaimed" structures. THe algorithm is then repeated untill all structures are assigned a cluster. The remaineder of structures are then assigned to clusters (this avoids having to calculate a full rms matrix) one by one. THe rule is that any structure joins the cluster to who\'s cluster center it is most similar to. If the closest cluster is more then "cluster_radius" away the structure will form a new cluster. This rule is applied to all remaining structures. Clusters can be size limited, sorted by energy etc.. (see options)\n'
p104
sS'AUTHOR'
p105
S'Phil Bradley/Mike Tyka\n'
p106
sS'OUTPUTS'
p107
S'\n'
p108
sS'REFERENCES'
p109
g108
sS'EXAMPLES'
p110
S'cluster.linuxgccrelease @flags > cluster.log\n'
p111
sS'INPUT'
p112
g108
sS'TIPS'
p113
g108
sS'METADATA'
p114
S'This document was edited May 25th 2010 by Mike Tyka. This application in mini was created and documented by Mike Tyka,et al.\n'
p115
ssS'AbinitioRelax'
p116
(dp117
S'DESCRIPTION'
p118
S'Abinitio Structure Prediction - The AbinitioRelax application consists of two main steps. The first step is a coarse-grained fragment-based search through conformational space using a knowledge-based "centroid" score function that favors protein-like features (Abinitio). The second optional step is all-atom refinement using the Rosetta full-atom forcefield (Relax). \n'
p119
sS'ALGORITHM'
p120
S'The "Relax" step is considerably more compute-intensive and time-consuming than the first step. The example run described above in the Code and Demo section takes around 8 minutes to generate one model of a 117 residue protein on a modern computer. A single AbinitioRelax run can generate a user defined number of models via a command line option (see Options section below). For increased conformational sampling, this application is easily parallelized by executing numerous jobs each using a unique random number seed (see Options section below). This is typically done by submitting multiple jobs to a computer cluster or distributed grid. Since the full-atom energy function is very sensitive to imperfect atomic interactions and more noise will exist with insufficient sampling, convergence towards the native structure may require a significant amount of sampling. Additionally, to increase your chance of sampling the correct topology, a diverse set of homologous sequences, preferably with sequence changes that may have a greater impact on sampling like deletions and differences in conserved positions, may also be run since a homologue may converge towards the native structure with significantly less sampling (see Bradley et al reference).\n'
p121
sS'AUTHOR'
p122
S'David E Kim\n'
p123
sS'OUTPUTS'
p124
S'Generates pdb files and an energy file, or a silent output file. Example: rosetta_demos/abinitio/input_files/S_00000001.pdb, rosetta_demos/abinitio/input_files/score.fsc, and rosetta_demos/abinitio/output_files/default.out (silent output file).\n'
p125
sS'ANALYSIS'
p126
S'We recommend generating up to 20,000 to 30,000 models of the target sequence and of up to 10 homologs and then using the Cluster application to identify the most frequently sampled conformations. In a general case, at least one of the top 5-10 clusters by size may have models with the lowest rmsd to the native structure.\\nIn an ideal case, your sequence will have many homologs identified by search tools like PSI-BLAST. Sequence alignments can be extremely helpful in model selection. For example, conserved hydrophobic positions most likely represent the core of the protein so models that have sidechains exposed in such positions may be discarded. The same logic applies to conserved polar positions which are most likely on the surface. Additionally, conserved cysteine pairs may represent disulphides. Tools like Jalview to view alignments and PyMOL to view models are extremely helpful for model selection in this respect.\\nScore versus RMSD plots may be helpful for identifying convergence towards the native structure for the target sequence and homologs. For example, the lowest scoring model can be used for the in:file:native input option when rescoring models with the score.linuxgccrelease score application. A score versus RMSD plot from the resulting score file may show convergence (an energy funnel) towards the lowest scoring model. If an energy funnel exists, the lowest scoring model has a greater chance of being near-native.\\nLowest scoring models that are in a cluster and that have a topology represented in the PDB also have a greater chance of being correct. Structure-structure comparison tools like Dali or Mammoth can be used to search against the PDB database.\n'
p127
sS'REFERENCES'
p128
S'Srivatsan Raman, Robert Vernon, James Thompson, Michael Tyka, Ruslan Sadreyev,Jimin Pei, David Kim, Elizabeth Kellogg, Frank DiMaio, Oliver Lange, Lisa Kinch, Will Sheffler, Bong-Hyun Kim, Rhiju Das, Nick V. Grishin, and David Baker. (2009) Structure prediction for CASP8 with all-atom refinement using Rosetta. Proteins 77 Suppl 9:89-99\\nBradley P, Misura KM, Baker D (2005). Toward high-resolution de novo structure prediction for small proteins. Science 309, 1868-71\\nBonneau R, Strauss CE, Rohl CA, Chivian D, Bradley P, Malmstrom L, Robertson T, Baker D. (2002) De novo prediction of three-dimensional structures for major protein families. J Mol Biol 322(1):65-78.\\nBonneau R, Tsai J, Ruczinski I, Chivian D, Rohl C, Strauss CE, Baker D. (2001) Rosetta in CASP4: progress in ab initio protein structure prediction. Proteins Suppl 5:119-26.\\nSimons KT, Ruczinski I, Kooperberg C, Fox B, Bystroff C, Baker D. (1999) Improved recognition of native-like protein structures using a combination of sequence-dependent and sequence-independent features of proteins. Proteins 34(1) 82-95\\nSimons KT, Kooperberg C, Huang E, Baker, D. (1997) Assembly of protein tertiary structures from fragments with similar local sequences using simulate anealing and Bayesian scoring functions. J Mol Biol 268:209-25.\n'
p129
sS'EXAMPLES'
p130
S'The ab initio executable is in rosetta_source/src/apps/public/AbinitioRelax.cc. The source code for the ab initio protocol is in rosetta_source/src/protocols/abinitio/AbrelaxApplication.cc. See the rosetta_demos/abinitio directory for an example ab initio run which includes input files, expected output files, and an example run log. The example command exists in rosetta_demos/abinitio/readme.txt. Input files exist in rosetta_demos/abinitio/input_files. Expected output files exist in rosetta_demos/abinitio/output_files. An example run log exist in rosetta_demos/abinitio/log.\n'
p131
sS'INPUT'
p132
S'Fasta file. Contains the amino acid protein sequence in fasta format. Example: rosetta_demos/abinitio/input_files/1elwA.fasta. Create using the GUI or by hand.  Fragments files. Generate structural fragment libraries using either the publicly available webserver (http://robetta.bakerlab.org/fragmentsubmit.jsp) or the GUI.  Native structure (optional). The native PDB structure may be used for benchmarking. When used, the RMSD to native is calculated for each model and provided as an extra column in the score line.  Psipred secondary structure prediction psipred_ss2 file (optional). The Psipred secondary structure prediction file is necessary when the -use_filters and -kill_hairpins options are used (see below). Note: the fragment webserver runs Psipred and provides the psipred_ss2 output file. Can make this using the GUI as well.\n'
p133
sS'TIPS'
p134
S'The AbinitioRelax application performs best for small monomeric proteins that are less than 100 residues in length. It is possible to get accurate predictions for some proteins up to around 150 residues, however, larger proteins require significantly more computing as the conformational space is vastly increased. It is not uncommon to sample in the range of 20,000 to 200,000 models in order to converge towards the native structure. The following references provide information relevant to the sampling problem:\\nBradley P, Misura KM, Baker D (2005). Toward high-resolution de novo structure prediction for small proteins. Science 309, 1868-71.\\nKim DE, Blum B, Bradley P, Baker D (2009). Sampling bottlenecks in de novo protein structure prediction. J Mol Biol 393, 249-60.\\n\\n  As stated above, it is beneficial to try to identify homologous sequences to run along with the target sequence (see Bradley et al reference). Homologs can be identified using search tools like PSI-BLAST to search the non-redundant sequence database (NCBI nr database) or Pfam. Using a sequence alignment viewer like Jalview is very useful to help select an optimal set of homologs to run and also to aid in model selection. Typically we look for a diverse set of homologs (up to 10) with differences in conserved positions and deletions which may represent a truncated loop or disordered region. Small changes in sequence can have a large impact on the topologies that are sampled, for example, a polar residue at a conserved hydrophobic position can have a big effect, i.e. the native topology may not be sampled because the full-atom Rosetta score will highly disfavor a polar residue buried in a hydrophobic core. It is also important to identify and trim disordered termini using publicly available programs like Disopred or metaPrDOS. Signal sequences should also be identified and trimmed using publicly available programs like SignalP. This protocol is not developed for membrane proteins. If transmembrane helices are predicted using programs like TMHMM, please refer to our Membrane ab initio application.\n'
p135
sS'METADATA'
p136
S'This document was last updated on November, 2010 by David E Kim <dekim@uw.edu>. The PI is David Baker <dabaker@uw.edu>. The AbinitioRelax application was developed by numerous Rosetta Commons members, primarily:Kim Simons, Richard Bonneau, Kira Misura, Phil Bradley, Oliver Lange, Michael Tyka, Robert Vernon\n'
p137
ssS'AnchoredPDBCreator'
p138
(dp139
S'LIMITATIONS'
p140
S"At this time, the anchor PDB must be a minimum of three residues (2 or 1 will cause a crash or bizzarre other errors). If you ultimately want an shorter anchor, that is fine - use a length of 3 here, then tell AnchoredDesign that it's only 2, and it will happily mutate/move the residue you didn't want as part of the anchor. I intend to fix this bug eventually...\n"
p141
sS'INPUTS'
p142
S"AnchoredPDBCreator needs three PDBs (scaffold for design, target for binding, and an anchor that was derived from the target's old binding partner) and a description of where to insert the anchor into the target. You should prepare the scaffold PDB by deleting any residues you want replaced by the anchor.The scaffold loop specification is a one-line file with four whitespace-delimited values: the chain letter of the scaffold chain you want to insert into, the start of the loop you want to insert into, the residue you want the anchor to be inserted immediately after, and the end of the loop you want to insert into. For example:B 18 23 30  This will treat the scaffold loop from 18-30 flexibly, except the anchor. The anchor will be inserted after 23. Let's say the anchor is 4 residues long, and residues 24-27 are deleted from the scaffold. In this case, PDB residues 18-23 and 28-30 will be treated as one flexible loop with a constant portion in the middle (the anchor). If 24-27 existed in the scaffold then their numbering would be bumped up and the eventual flexible loop would number 18-34.\n"
p143
sS'DESCRIPTION'
p144
S'AnchoredPDBCreator exists to create starting files for AnchoredDesign - it will take a scaffold, target, and the anchor and brew a starting structure for AnchoredDesign. \n'
p145
sS'ALGORITHM'
p146
S"It is implemented separately so that the user can look at the input structure for AnchoredDesign.cc to design its resfile, etc. with the right numbering.  This app assumes you have three structures: an anchor and target in separate PDBs, which are bound to each other relative to the same coordinate origin, and a scaffold PDB (origin unimportant). Generally your anchor and target will be from the same crystal structure. The protocol inserts the anchor's sequence AND COORDINATES into a loop of the scaffold, and runs a fast dirty minimization to try to keep the scaffold and target from crashing into each other. It is not intended to produce low-energy structures...just structures good enough to look at their numbering while designing a resfile, etc. It tries to ensure that the modified loop has at most ONE chainbreak or bad omega angle (which can be fixed in later loop modeling). It does NOT try to prevent eclipsing of the loop onto the scaffold or the target - again, AnchoredDesign will deal with that.  The workhorse is protocols::pose_manipulation::insert_pose_into_pose. This code inserts the anchor into the scaffold and attempts loop closure (this is effectively a domain insertion problem). Interfacing with the target occurs inside AnchoredPDBCreator proper.\n"
p147
sS'AUTHOR'
p148
S'Steven Lewis smlewi@gmail.com\n'
p149
sS'OUTPUTS'
p150
S'The protocol will produce nstruct result structures with the anchor inserted into the scaffold, and the anchor aligned against the target. You should expect the scaffold loop newly containing the anchor to be closed and in a reasonable conformation. You should not expect the scaffold/target interface to be reasonable (or even not clashing; they may overlap). Generating the interface is the job of AnchoredDesign instead.The scorefile will report LAM_total, which is the score from Loop Analyzer Mover. This sums the scores for he ramachandran, omega, chainbreak, and peptide bond scores for the residues in the loops.\n'
p151
sS'ANALYSIS'
p152
S'Essentially the only quality of the results that matters is the quality of loop closure; everything else is meant to be handled by AnchoredDesign. To pick a model for AnchoredDesign, sort by LAM_total, examine the lowest scoring handful of models to see which you like best, and move on.\n'
p153
sS'REFERENCES'
p154
S'Lewis SM, Kuhlman BA. Anchored design of protein-protein interfaces. PLoS One. 2011;6(6):e20872. Epub 2011 Jun 17. Gulyani A, Vitriol E, Allen R, Wu J, Gremyachinskiy D, Lewis S, Dewar B, Graves LM, Kay BK, Kuhlman B, Elston T, Hahn KM. A biosensor generated via high-throughput screening quantifies cell edge Src dynamics. Nat check_button_m Biol. 2011 Jun 12;7(7):437-44. doi: 10.1038/nchembio.585.\n'
p155
sS'EXAMPLES'
p156
S"The code is at mini/src/apps/public/interface_design/anchored_design/AnchoredPDBCreator.cc; there's an integration test at mini/test/integration/tests/AnchoredDesign/. There is a more extensive demo with more documentation at demo/AnchoredDesign, or in the demo section of the release.\n"
p157
sS'TIPS'
p158
S'You will need healthy amounts of sampling. I use 50 nstruct and 100,000 cycles.\n'
p159
sS'METADATA'
p160
S'Code and documentation by Steven Lewis smlewi@gmail.com. This document was last updated 6/24/11 by Steven Lewis. The PI was Brian Kuhlman, bkuhlman@email.unc.edu.\n'
p161
ssS'minimize_with_cst'
p162
(dp163
S'LIMITATIONS'
p164
S'Only known use (by me) is for ddg_monomer\n'
p165
sS'INPUTS'
p166
S'Will only take a LIST of pdbs.  If you run this app within the toolkit and have a directory for input specified, the program will create a list of PDBs and run that.\n'
p167
sS'DESCRIPTION'
p168
S'Used to minimize structure before ddg_monomer.  Use default options.\n'
p169
sS'AUTHOR'
p170
S'Andrew Leaver-Fay and Elizabeth Kellogg\n'
p171
sS'OUTPUTS'
p172
S'Need the Log file, and the minimized structures for ddg_monomer.\n'
p173
sS'ANALYSIS'
p174
S'None that I know of.\n'
p175
sS'EXAMPLES'
p176
S'See Tips\n'
p177
sS'TIPS'
p178
S'Run with these options: -in:file:l lst  -in:file:fullatom -ignore_unrecognized_res -fa_max_dis 9.0 -database /path/to/minirosetta_database/ -ddg::harmonic_ca_tether 0.5 -score:weights standard -ddg::constraint_weight 1.0 -ddg::out_pdb_prefix min_cst_0.5 -ddg::sc_min_only false -score:patch minirosetta_database/scoring/weights/score12.wts_patch > mincst.log\n'
p179
sS'METADATA'
p180
S'Documentation created by Andrew Leaver-Fay, with some extra stuff from me.\n'
p181
ss.