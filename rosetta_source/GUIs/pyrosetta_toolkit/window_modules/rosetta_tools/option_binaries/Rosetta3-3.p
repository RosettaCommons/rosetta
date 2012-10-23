(dp0
S'--Docking--'
p1
(dp2
S'-docking:dock_ppk'
p3
S"Docking prepack mode'. Default='false' [Boolean]\n"
p4
sS'-dock_lowres_filter'
p5
S"<INTERCHAIN_CONTACT CUTOFF> <INTERCHAIN_VDW CUTOFF>. Default values for protein docking are 10.0 and 1.0',n='2' [RealVector]\n"
p6
sS'-docking:dock_lowres_filter'
p7
S'Manually sets the lowres docking filter:\n'
p8
sS'-docking:norepack1'
p9
S"Do not repack the side-chains of partner 1. Default='false' [Boolean]\n"
p10
sS'-docking:dock_mcm_trans_magnitude'
p11
S"The magnitude of the translational perturbation during mcm in docking. Default='0.1'. [Real]\n"
p12
sS'-docking:uniform_trans'
p13
S'Uniform random repositioning within a sphere of the given radius. [Real]\n'
p14
sS'-docking:view'
p15
S"Decide whether to use the viewer (graphical) or not. Default='false'. [Boolean]\n"
p16
sS'-docking:norepack2'
p17
S"Do not repack the side-chains of partner 2. Default='false' [Boolean]\n"
p18
sS'-docking:partners'
p19
S"Defines docking partners by ChainID. Example: docking chains L+H with A is -partners LH_A', Default='_'. [String]\n"
p20
sS'-docking:dock_min'
p21
S"Minimize the final \xef\xbf\xbd\xef\xbf\xbd\xce\xa9fullatom structure.Default='false' [Boolean]\n"
p22
sS'-docking:dock_mcm_rot_magnitude'
p23
S"The magnitude of the rotational perturbation during mcm in docking. Default='5.0'. [Real]\n"
p24
sS'-docking:randomize2'
p25
S"Randomize the second docking partner.Default='false'[Boolean]\n"
p26
sS'-docking:randomize1'
p27
S"Randomize the first docking partner.Default='false'[Boolean]\n"
p28
sS'-docking:docking_centroid_outer_cycles'
p29
S"Outer cycles during docking rigid body adaptive moves. Default='10' [Integer]\n"
p30
sS'-docking:spin'
p31
S"Spin a second docking partner around axes from center of mass of partner1 to partner2. Default='false'. [Boolean]\n"
p32
sS'-docking:docking_local_refine'
p33
S"Do a local refinement of the docking position(high resolution) Default='false'. [Boolean]\n"
p34
sS'-docking:docking_centroid_inner_cycles'
p35
S"Inner cycles during docking rigid body adaptive moves. Default='50' [Integer]\n"
p36
sS'-docking:fake_native'
p37
S"Decide whether to use low res docking filters after centroid mode. Default='false' [Boolean]\n"
p38
sS'-docking:dock_mcm'
p39
S"Do a monte-carlo minimize search.Default='true'[Boolean]\n"
p40
sS'-docking:dock_pert'
p41
S"Do a small perturbation with partner two: -dock_pert DEGREES ANGSTROMS. Good values for protein docking are 3 deg and 8 A.', n='2' [RealVector]\n"
p42
ssS'--Packing--'
p43
(dp44
S'-packing:ex2aro:level'
p45
S'Use extra chi2 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p46
sS'-packing:ex1aro:level'
p47
S'Use extra chi1 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p48
sS'-packing:solvate'
p49
S"Add explicit water, but don't try to place water such that it bridges Hbonds, just put it on every available Hbond donor/acceptor where there's no clash (implies explicit_h2o). [Boolean]\n"
p50
sS'-packing:multi_cool_annealer'
p51
S'Alternate annealer for packing.  Runs multiple quence cycles in a first cooling stage, and tracks the N best network states it observes. It then runs low-temperature rotamer substitutions with repeated quenching starting from each of these N best network states. 10 is recommended.\n'
p52
sS'-packing:ex2:level'
p53
S'Use extra chi2 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p54
sS'-packing:no_his_his_pairE'
p55
S'Set pair term for His-His to zero. [Boolean]\n'
p56
sS'-packing:ex4:operate'
p57
S'Apply special operations (see RotamerOperation class) on ex4 rotamers. [Boolean]\n'
p58
sS'-packing:pack_missing_sidechains'
p59
S'Run packer to fix residues with missing sidechain density at PDB load. default="true"\' [Boolean]\n'
p60
sS'-packing:fix_his_tautomer'
p61
S'Seqpos numbers of his residus whose tautomer should be fixed during repacking. [IntegerVecter]\n'
p62
sS'-packing:linmem_ig'
p63
S" =Force the packer to use the linear memory interaction graph; each RPE may be computed more than once,but recently-computed RPEs are reused.  The integer parameter specifies the number of recent rotamers to store RPEs for.  10 is the recommended size. Memory use scales linearly with the number of rotamers at about 200 bytes per rotamer per recent rotamers to store RPEs for(~4 KB per rotamer by default). default='10' [Integer]\n"
p64
sS'-packing:ex2xaro_exposed:level'
p65
S'Use extra chi2 sub-rotamers for all aromatic residues that pass the extrachi_cutoff.\n'
p66
sS'-packing:extrachi_cutoff'
p67
S'Number of neighbors a residue must have before extra rotamers are used. default: 18 [Integerx]\n'
p68
sS'-packing:ex4'
p69
S'Use extra chi4 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p70
sS'-packing:ex1aro_exposed'
p71
S'Use extra chi1 sub-rotamers for all aromatic residues. [Boolean]\n'
p72
sS'-packing:ex1'
p73
S'Use extra chi1 sub-rotamers for all residues that pass the extrachi_cutoff. [Boolean]\n'
p74
sS'-packing:ex3'
p75
S'Use extra chi3 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p76
sS'-packing:ex2'
p77
S'Use extra chi2 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p78
sS'-packing:ex1aro_exposed:level'
p79
S'Use extra chi1 sub-rotamers for all aromatic residues that pass the extrachi_cutoff.\n'
p80
sS'-packing:ex2aro_exposed'
p81
S'Use extra chi2 sub-rotamers for all aromatic residues. [Boolean]\n'
p82
sS'-packing:ex4:level'
p83
S'Use extra chi4 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p84
sS'-packing:use_input_sc'
p85
S'Use rotamers from input structure in packing. By default, input sidechain coords are NOTincluded in rotamer set but are discarded before the initial pack; with this flag, the inputrotamers will NOT be discarded. Note that once the starting rotamers are replaced by any mechanism, they are no longer included in the rotamer set. (rotamers included by coordinates)\n'
p86
sS'-packing:ex3:operate'
p87
S'Apply special operations (see RotamerOperation class) on ex3 rotamers. [Boolean]\n'
p88
sS'-packing:resfile'
p89
S"Resfile filename(s).  Most protocols use only the first and will ignore the rest; it does not track against -s or -l automatically. Default='resfile' [FileVector]\n"
p90
sS'-packing:ex1aro'
p91
S'Use extra chi1 sub-rotamers for aromatic residues that pass the extrachi_cutoff. [Boolean]\n'
p92
sS'-packing:ex1:level'
p93
S'Use extra chi1 sub-rotamers for all residues that pass the extrachi_cutoff. [Integer]\n'
p94
sS'-packing:ex3:level'
p95
S' Use extra chi3 sub-rotamers for all residues that pass the extrachi_cutoff.\n'
p96
sS'-packing:no_optH'
p97
S' =Do not optimize hydrogen placement at the time of a PDB load. default="true" [Boolean]\n'
p98
sS'-packing:ex2aro'
p99
S'Use extra chi2 sub-rotamers for aromatic residues that pass the extrachi_cutoff.\n'
p100
sS'-packing:ex1:operate'
p101
S'Apply special operations (see RotamerOperation class) on ex2 rotamers. [Boolean]\n'
p102
ssS'--Experiments--'
p103
(dp104
sS'--Run--'
p105
(dp106
S'-run:delay'
p107
S'Wait N seconds before doing anything at all. Useful for cluster job staggering. Default=0. [String]\n'
p108
sS'-run:maxruntime'
p109
S'Maximum runtime in seconds.JobDistributor will signal end if time is exceeded no matter how many jobs were finished. Default=-1. [Integer]\n'
p110
sS'-run:constant_seed '
p111
S'Use a constant seed (1111111 unless specified). [Boolean]\n'
p112
sS'-run:jran'
p113
S'Specify seed (requires -constant_seed) Default=1111111. [Integer]\n'
p114
sS'-run:shuffle'
p115
S"Shuffle job order. Default=false. [Boolean]'\n"
p116
sS'-run:seed_offset'
p117
S'This value will be added to the random number seed. Useful when using time as seed and submitting many jobs to clusters.Using the condor job id will force jobs that are started in the same second to still have different initial seeds".default=0. [Integer]\n'
p118
sS'-run:nodelay'
p119
S'Do not delay launch of Rosetta [Boolean]\n'
p120
sS'-run:rng'
p121
S'Random number generation algorithm: Currently only mt19937 is a accepted here. Default=mt19937. legal=[mt19937]. [String]\n'
p122
sS'-run:random_delay'
p123
S'Wait a random amount of 0..N seconds before doing anything at all. Useful for cluster job staggering." Default=0. [Integer]\n'
p124
sS'-run:rng_seed_device'
p125
S' Obtain random number seed from specified device. Default=/dev/urandom. [String]\n'
p126
sS'-run:use_time_as_seed'
p127
S'Use time as random number seed instead of default rng seed device. [Boolean]\n'
p128
ssS'relax'
p129
(dp130
S'-score:patch'
p131
S'Supply a different patch file (default is score12)\n'
p132
sS'-score:weights'
p133
S'Supply a different weights file (default is standard) (Scoretypes can be found in the FineControl Menu)\n'
p134
sS'-in:file:s'
p135
S'Input pdb file(s)\n'
p136
sS'-nstruct'
p137
S'Make how many decoys per input structure ?\n'
p138
sS'-relax:classic'
p139
S'Do an old old deprecated "classic" relax mode (slow and poor performance)\n'
p140
sS'-relax:fast'
p141
S'Do a simple, small cycle number (5)  fast relax (DEFAULT)\n'
p142
sS'-run:shuffle'
p143
S'Use shuffle mode, treat structures in random order\n'
p144
sS'-relax:thorough'
p145
S'Do a preset, large cycle number (15) fast relax\n'
p146
sS'-in:file:silent'
p147
S'Input silent file\n'
p148
sS'-in:file:fullatom'
p149
S'Read as fullatom input structure\n'
p150
sS'-relax:script'
p151
S'Do custom script relax (you can supply a custom script file) (Use the Relax Script builder in the GUI)\n'
p152
ssS'fixbb'
p153
(dp154
S'-score:weights'
p155
S'Set the weights file to be used.\n'
p156
sS'-ex2aro'
p157
S'Increase chi2 rotamer sampling for buried* aromatic** residues +/- 1 standard deviation\n'
p158
sS'-ex3:level'
p159
S'Increase chi1 sampling for buried* residues to the given sampling level\n'
p160
sS'-score:weights_patch'
p161
S'Set the weights patch file to be used.'
p162
sS'-use_input_sc'
p163
S'Include the side chain from the input pdb.  False by default.  Including the input sidechain is "cheating" if your goal is to measure sequence recovery, but a good idea if your goal is to redesign the input sequence for eventual synthesis. Buried residues are those with #CBeta-neighbors >= threshold (default 18) within 10 Angstroms. This threshold can be controlled by the -extrachi_cutoff flag.\n'
p164
sS'-ex2:level'
p165
S'Increase chi1 sampling for buried* residues to the given sampling level\n'
p166
sS'-ex1aro'
p167
S'Increase chi1 rotamer sampling for buried* aromatic** residues +/- 1 standard deviation\n'
p168
sS'-minimize_sidechains'
p169
S'Follow the packing phase with gradient-based minimization of the sidechains for residues that were either repacked or designed in the packing phase.\n'
p170
sS'-l'
p171
S'A file that lists one or more pdbs to run <code>fixbb</code> upon.\n'
p172
sS'-resfile'
p173
S'The resfile that is to be used for this job\n'
p174
sS'-linmem_ig'
p175
S'Activate the linear-memory interaction graph [Leaver-Fay et al. 2008]\n'
p176
sS'-ex1aro:level'
p177
S'Increase chi1 sampling for buried* aromatic residues to the given sampling level\n'
p178
sS'-ex2aro:level'
p179
S'Increase chi1 sampling for buried* aromatic residues to the given sampling level\n'
p180
sS'-s'
p181
S'A list of one or more pdbs to run <code>fixbb</code> upon.\n'
p182
sS'-nstruct'
p183
S'The number of iterations to perform per input structure; e.g. with 10 input structures and an -nstruct of 10, 100 trajectories will be performed. Interaction Graph (Default is to precompute all rotamer pair energies)\n'
p184
sS'-extrachi_cutoff'
p185
S'Set the number of cbeta neighbors (counting its own) at which a residue is considered buried.  A value of "1" will mean that all residues are considered buried for the purpose of rotamer building.\n'
p186
sS'-overwrite'
p187
S'Overwrite the output files, even if they already exist.\n'
p188
sS'-min_type'
p189
S'When combined with the -minimize_sidechains flag, specifies the line-search algorithm to use in the gradient-based minimization . "dfpmin" by default.\n'
p190
sS'-ex4:level'
p191
S'Increase chi1 sampling for buried* residues to the given sampling level\n'
p192
sS'-lazy_ig'
p193
S'Activate the lazy interaction graph Rotamers\n'
p194
sS'-multi_cool_annealer'
p195
S'Use an alternate annealer that spends more time at low temperatures.  This annealer produces consistently lower energies than the standard annealer.\n'
p196
sS'-ex4'
p197
S'Increase chi4 rotamer sampling for buried* residues +/- 1 standard deviation\n'
p198
sS'-ex2'
p199
S'Increase chi2 rotamer sampling for buried* residues +/- 1 standard deviation\n'
p200
sS'-ex3'
p201
S'Increase chi3 rotamer sampling for buried* residues +/- 1 standard deviation\n'
p202
sS'-ex1'
p203
S'Increase chi1 rotamer sampling for buried* residues +/- 1 standard deviation\n'
p204
sS'-ex1:level'
p205
S'Increase chi1 sampling for buried* residues to the given sampling level\n'
p206
sS'-preserve_input_cb'
p207
S'Do not idealize the CA-CB bond vector -- instead, use the CB coordinates of the input pdb.\n'
p208
ssS'--Score--'
p209
(dp210
S'-score:output_etables'
p211
S'Write out etables to files with given prefix. [String]\n'
p212
sS'-score:optH_weights'
p213
S'Name of weights file (without extension .wts) to use during optH. [String]\n'
p214
sS'-score:weights'
p215
S'Name of weights file (without extension .wts) Default="standard". [String]\n'
p216
sS'-score:etable_lr'
p217
S'Lowers energy well at 6.5A. [Real]\n'
p218
sS'-score:optH_patch'
p219
S'Name of weights file (without extension .wts) to use during optH. [String]\n'
p220
sS'-score:no_smooth_etables'
p221
S'Revert to old style etables. [Boolean]\n'
p222
sS'-score:empty'
p223
S'Make an empty score - i.e. NO scoring. [Boolean]\n'
p224
sS'-score:dun10'
p225
S'Use the 2010 Dunbrack library [Boolean?]\n'
p226
sS'-score:input_etables'
p227
S'Read etables from files with given prefix. [String]\n'
p228
sS'-score:ramaneighbors true'
p229
S'Uses neighbor-dependent ramachandran maps Default=false [Boolean]\n'
p230
sS'-score:fa_Hatr'
p231
S'Turn on Lennard Jones attractive term for hydrogen atoms. [Boolean]\n'
p232
sS'-score:fa_max_dis'
p233
S'How far does the FA pair potential go out to. Default = 6.0. [Real]\n'
p234
sS'-score:rms_target'
p235
S'Target of RMS optimization for RMS_Energy EnergyMethod. Default=0.0 [Real]\n'
p236
sS'-score:dun08'
p237
S'Use the 2008 Dunbrack library instead of the 2002 library. [Boolean]\n'
p238
sS'-score:patch score12'
p239
S'Name of patch file (without extension) Default="". [String]\n'
p240
ssS'--Loops--'
p241
(dp242
S'-loops:combine_rate'
p243
S"Combine successive loops at this rate. Default='0.0' [Real]\n"
p244
sS'-loops:idealize_before_relax'
p245
S"Idealize the structure before going into fa relax. Default='false' [Boolean]\n"
p246
sS'-loops:fast'
p247
S"Reduce number of simulation cycles in loop modeling. Default='false' [Boolean]\n"
p248
sS'-loops:random_grow_loops_by'
p249
S"Randomly grow loops by up to this number of residues on either side.Default='0.0' [Real]\n"
p250
sS'-loops:remodel'
p251
S"Legal values: 'no','perturb_ccd','perturb_kic'. This option controls the centroid-stage loop remodeling (where large perturbations are attempted). Default = 'quick_ccd' [String]\n"
p252
sS'-loops:relax'
p253
S'Legal values: \'no\',\'fastrelax\',\'shortrelax\',\'fullrelax\'    [String] Controls whether a full-structure relax occurs after loop modeling.  Defaults to "no".\n'
p254
sS'-loops:vicinity_sampling'
p255
S"Sample non-pivot torsions within a vicinity of their input values. Default='false' [Boolean]\n"
p256
sS'-loops:build_all_loops'
p257
S"Build all loops(no skip). Default='false' [Boolean]\n"
p258
sS'-loops:kic_max_seglen'
p259
S"Maximum number residues in a KIC move segment. Default='12' [Integer]\n"
p260
sS'-loops:build_attempts'
p261
S"Build attempts per growth attempt. Default='3' [Integer]\n"
p262
sS'-loops:remodel_final_temp'
p263
S'Final temperature for simulated annealing for perturb_kic\n'
p264
sS'-loops:grow_attempts'
p265
S"Total loop growth attempts. Default='7' [Integer]\n"
p266
sS'-loops:random_order'
p267
S"Build in random order. Default='false' [Boolean]\n"
p268
sS'-loops:extended'
p269
S"Force extended on loops, independent of loop input file. Default='false' [Boolean]\n"
p270
sS'-loops:frag_sizes'
p271
S'(non-kic modes only) Length of fragments to be used in loop modeling (in other words, number of amino acids/residues for each segment. For each type of fragment). defaults to 9 3 1 [IntegerVector]\n'
p272
sS'-lopps:input_pdb'
p273
S"Template pdb file. default='input_pdb' [File]\n"
p274
sS'-loops:build_initial'
p275
S"Preceed loopmodelling with a initial round of just removing the missing densities and building simple loop. Default = 'false'. [Boolean]\n"
p276
sS'-loops:accept_aborted_loops'
p277
S"Accept aborted loops. Default='false' [Boolean]\n"
p278
sS'-loops:refine_only'
p279
S"Perform full atom refinement only on loops. Default='false' [Boolean]\n"
p280
sS'-loops:loop_file'
p281
S"Loop definition file. default='loop_file' [File]\n"
p282
sS'-loops:neighbor_dist'
p283
S"Distance in angstroms to optimize side-chains near segment moved by KIC. Default='10.0' [Real]\n"
p284
sS'-loops:fix_natsc'
p285
S"Fix sidechains in template region in loop modeling. Default='false' [Boolean]\n"
p286
sS'-loops:remodel_init_temp'
p287
S'Initial temperature for simulated annealing for perturb_kic\n'
p288
sS'-loops:vicinity_degree'
p289
S" Number of degrees to deviate from (current) non-pivot torsions. Default='1.0' [Real]\n"
p290
sS'-loops:mm_loop_file'
p291
S"Loop definition file. default='loop_file' [File]\n"
p292
sS'-loops:build_specific_loops'
p293
S'Numbers of the loops to be built. [IntegerVector]\n'
p294
sS'-loops:select_best_loop_from'
p295
S"Keep building loops until N and choose best (by score). Default='1' [Integer]\n"
p296
sS'-loops:fa_input'
p297
S"Input structures are in full atom format. Default='false' [Boolean]\n"
p298
sS'-loops:refine_init_temp'
p299
S'Initial temperature for simulated annealing for refine_kic\n'
p300
sS'-loops:frag_files'
p301
S'(non-kic modes only) Path/Name of fragment libraries files, defaults to frag9 frag3 frag1. [FileVector]\n'
p302
sS'-loops:refine'
p303
S"Legal values: 'no','refine_ccd','refine_kic'. This option controls the fullatom refinement stage of loop. modeling. Default = 'no' [String]\n"
p304
sS'-loops:refine_final_temp'
p305
S'Final temperature for simulated annealing for refine_kic\n'
p306
sS'-loops:strict_loops'
p307
S"Do not allow growing of loops. Default='false' [Boolean]\n"
p308
sS'-loops:intermedrelax'
p309
S"Currently not used; Eventually may provide for a full-pose relax between centroid and fullatom modeling. default = 'no' [String]\n"
p310
ssS'ddg_monomer'
p311
(dp312
S'-fa_max_dis'
p313
S'9.0 # optional -- if not given, the default value of 9.0 Angstroms is used.\n'
p314
sS'-ddg::ramp_repulsive true'
p315
S'true # perform three rounds of minimization (and not just the default 1 round) where the weight on the repulsive term is increased from 10% to 33% to 100%\n'
p316
sS'-ignore_unrecognized_res true'
p317
S'true # Recommended\n'
p318
sS'-mute'
p319
S'all # optional -- silence all of the log-file / stdout output generated by this protocol\n'
p320
sS'-ddg::min'
p321
S'false Low-Res, true High-Res # report the minimum energy\n'
p322
sS'-resfile'
p323
S'Resfile to be used instead of mutfile\n'
p324
sS'-ddg::mut_file'
p325
S'the list of point mutations to consider in this run\n'
p326
sS'-ddg::min_cst true'
p327
S'true # use distance restraints (aka constraints) during the backbone minimization phase\n'
p328
sS'-ddg::mean'
p329
S'true Low-Res, false High-Res # do not report the mean energy\n'
p330
sS'-ddg::local_opt_only'
p331
S'true Low-Res, false High-Res : local optimization restricts the sidechain optimization to only the 8 A neighborhood of the mutation (equivalent to row 13)\n'
p332
sS'-ddg:weight_file'
p333
S'soft_rep_design # Use soft-repulsive weights for the initial sidechain optimization stage\n'
p334
sS'-ddg::suppress_checkpointing'
p335
S'true # dont checkpoint\n'
p336
sS'-ddg:minimization_scorefunction'
p337
S'optional -- the weights file to use, if not given, then "score12" will be used (score12 = standard.wts + score12.wts_patch)\n'
p338
sS'-in:file:s'
p339
S'<pdbfile of the preminimized wildtype structure> # the PDB file of the structure on which point mutations should be made\n'
p340
sS'-ddg::iterations 50'
p341
S'50 # 50 is the recommended number of iterations\n'
p342
sS'-ddg::dump_pdbs'
p343
S'true # write out PDB files for the structures, one for the wildtype and one for the pointmutant for each iteration\n'
p344
sS'-ddg::output_silent'
p345
S'true # write output to a silent file\n'
p346
sS'-unmute'
p347
S'core.optimization.LineMinimizer # optional -- unsilence a particular tracer\n'
p348
sS'-in::file::fullatom'
p349
S'read the input PDB file as a fullatom structure\n'
p350
sS'-ddg::sc_min_only false'
p351
S'false # do not minimize only the backbone during the backbone minimization phase\n'
p352
sS'-constraints::cst_file'
p353
S'<cbeta-distance-constraint-file> # the set of constraints to use during minimization which should reflect distances in the original (non-pre-relaxed) structure\n'
p354
sS'-ddg::minimization_patch'
p355
S'optional -- the weight-patch file to apply to the weight file; does not have to be given\n'
p356
ssS'backrub'
p357
(dp358
S'-ignore_unrecognized_res'
p359
S'Ignore unrecognized atom times found in PDB\n'
p360
sS'-min_atoms 3'
p361
S'minimum backrub segment size in atoms (default 3)\n'
p362
sS'-mc_kt 0.6'
p363
S'value of kT for Monte Carlo (default 0.6, 0.3-0.4 recommended for > 100,000 step simulations)\n'
p364
sS'-l'
p365
S'File(s) containing list(s) of PDB files to process\n'
p366
sS'-resfile'
p367
S'Use Resfile to specify which residues can be repacked\n'
p368
sS'-cst_fa_file'
p369
S'constraints filename(s) for fullatom. When multiple files are given a *random* one will be picked.\n'
p370
sS'-sc_prob_withinrot 0.0'
p371
S'probability of sampling within the current rotamer (default 0.0)(For Off-rotamer searches)\n'
p372
sS'-pivot_atoms CA'
p373
S'main chain atoms usable as pivots (default CA)\n'
p374
sS'-movemap'
p375
S'MoveMap file specifying flexible torsions (You can create this using the toolkit)\n'
p376
sS'-sc_prob 0.25'
p377
S'probability of making a side chain move (default 0.25)\n'
p378
sS'-s'
p379
S'Name(s) of single PDB file(s) to process\n'
p380
sS'-backrub:trajectory'
p381
S'record a trajectory\n'
p382
sS'-sc_prob_uniform 0.1'
p383
S'probability of uniformly sampling chi angles (defualt 0.1)\n'
p384
sS'-cst_fa_weight 1'
p385
S'weight of the fullatom constraint terms (defualt 1)\n'
p386
sS'-nstruct'
p387
S'# of structures to output\n'
p388
sS'-sm_prob 0'
p389
S'probability of making a small move (default 0)\n'
p390
sS'-max_atoms 34'
p391
S'maximum backrub segment size in atoms (default 34)\n'
p392
sS'-pivot_residues'
p393
S'residues for which contiguous stretches can contain segments, in absolute residue numbers (defaults to all residues)\n'
p394
sS'-backrub:trajectory_gz'
p395
S'gzip the trajectory\n'
p396
sS'-initial_pack '
p397
S'force a repack/minimization at the beginning regardless of whether mutations are set in the resfile\n'
p398
sS'-minimize_movemap'
p399
S'MoveMap specifying degrees of freedom to be minimized after initial packing, occurs in three stages:1. CHI only  2. CHI+BB  3. CHI+BB+Jump\n'
p400
sS'-backrub:ntrials 10000'
p401
S'number of Monte Carlo trials to run (default 1000, 10000 recommended, 10,000 for production runs)\n'
p402
sS'-backrub:trajectory_stride'
p403
S'write out a trajectory frame every N steps\n'
p404
sS'-mm_bend_weight 1'
p405
S'weight of mm_bend bond angle energy term (default 1)\n'
p406
ssS'loopmodel'
p407
(dp408
S'-out::prefix'
p409
S'Set prefix for output\n'
p410
sS'-loops:kic_recover_last'
p411
S"(Expert) Keep the last sampled conformation at the end of each outer cycle instead of the lowest energy conformation sampled so far. default = 'false'. [Boolean]\n"
p412
sS'-run:test_cycles'
p413
S"(Expert) Sets the number of outer cycles and inner cycles to 3. For extremely quick testing and sanity checks, not for production runs. default = 'false'. [Boolean]\n"
p414
sS'-loops:fast'
p415
S"(Expert) Signifcantly reduces the number of inner cycles. For quick testing, not production runs. default = 'false'. [Boolean]\n"
p416
sS'-extrachi_cutoff 0'
p417
S'Set to 0 to include extra rotamers regardless of number of neighbors\n'
p418
sS'-ex1'
p419
S'Include extra chi1 rotamers (or also chi2, chi3, chi4). Using -ex1 -ex2 improves benchmark performance.\n'
p420
sS'-loops:vicinity_sampling'
p421
S"(Expert) Sample non-pivot torsions within a vicinity of their input values. default = 'false'. For a description of pivot and non-pivot torsions, please see Purpose, above. [Boolean]\n"
p422
sS'-out:pdb_gz'
p423
S'Create compressed output PDB files (using gzip), which saves a lot of space. These files can still be visualized normally with software such as Pymol.\n'
p424
sS'-loops::frag_sizes'
p425
S'9 3 or 1 (Create using Fragment Window...)\n'
p426
sS'-loops:kic_max_seglen'
p427
S"(Expert) Maximum number of residues in a KIC move segment. default = '12'. [Integer]\n"
p428
sS'-loops:optimize_only_kic_region_sidechains_after_move'
p429
S"(Expert) Should rotamer trials and minimization be performed after every KIC move but only within the loops:neighbor_dist of the residues in the moved KIC segment. Speeds up execution when using very large loop definitions (such as when whole chains are used for ensemble generation). default = 'false'. [Boolean]\n"
p430
sS'-loops:remodel_final_temp'
p431
S"(Expert) Final temperature for simulated annealing in 'perturb_kic'. default = '1.0'. [Float]\n"
p432
sS'-extra_res_cen'
p433
S'Path to centroid parameters file for non-protein atoms or ligands\n'
p434
sS'-loops::loop_file'
p435
S'(You can create this using the export menu in the PyRosetta GUI)\n'
p436
sS'-loops:outer_cycles'
p437
S"(Expert) Number of outer cycles for Monte Carlo (described above in Protocol). default = '3'. [Integer]\n"
p438
sS'-loops::input_pdb'
p439
S'You can do most of this within the main GUI (Note that this is not -s or -l like normal)\n'
p440
sS'-loops:neighbor_dist'
p441
S"(Expert) Only optimize side-chains with C-beta atoms within this many angstroms of any loop C-beta atom. default = '10.0'. To speed up runs, try '6.0'. [Float]\n"
p442
sS'-loops:fix_natsc'
p443
S"Don't repack, rotamer trial, or minimize loop residue neighbors. default = 'false'. [Bolean]\n"
p444
sS'-loops:remodel_init_temp'
p445
S"(Expert) Initial temperature for simulated annealing in 'perturb_kic'. default = '2.0'. [Float]\n"
p446
sS'-loops:vicinity_degree'
p447
S"(Expert) Number of degrees allowed to deviate from current non-pivot torsions when using vicinty sampling (smaller number makes tighter sampling). default = '1.0'. [Float]\n"
p448
sS'-loops:max_inner_cycles'
p449
S'(Expert) Maximum number of inner cycles for Monte Carlo (default described above in Protocol). [Integer]\n'
p450
sS'-loops::refine'
p451
S'refine_ccd, refine_alc, refine_kic (default=no)\n'
p452
sS'-loops::random_loop'
p453
S'Randomize loop?\n'
p454
sS'-extra_res_fa'
p455
S'Path to all-atom parameters file for non-protein atoms or ligands\n'
p456
sS'-loops:refine_init_temp'
p457
S"(Expert) Initial temperature for simulated annealing in 'refine_kic'. default = '1.5'. [Float]\n"
p458
sS'-loops:max_kic_build_attempts'
p459
S"(Expert) Number of times to attempt initial closure in 'perturb_kic' protocol. Try increasing to 10000 if initial closure is failing. default = 100. [Integer]\n"
p460
sS'-ex2'
p461
S'Include extra chi2 rotamers (or also chi2, chi3, chi4). Using -ex1 -ex2 improves benchmark performance.\n'
p462
sS'-in::file::fullatom'
p463
S'Full Atom (Use this or else loop modeling will discard the input sidechains and repack!!!!)\n'
p464
sS'-loops::remodel'
p465
S'perturb_ccd (Use This!), perturb_alc, quick_ccd, old_loop_relax, perturb_kic\n'
p466
sS'-loops::relax'
p467
S"fullrelax, shortrelax, fastrelax'\n"
p468
sS'-loops::ccd_closure'
p469
S'Do CCD clusure on loops?\n'
p470
sS'-loops:refine_final_temp'
p471
S"(Expert) Final temperature for simulated annealing in 'refine_kic'. default = '0.5'. [Float]\n"
p472
sS'-loops::frag_files'
p473
S' Frag files 9, 3, 1\n'
p474
ssS'cluster'
p475
(dp476
S'-score:patch'
p477
S'Supply a different patch file (default is score12)\n'
p478
sS'-out:file:silent'
p479
S'Output silent structures instead of PDBs \n'
p480
sS'-in:file:s'
p481
S'Input pdb file(s)  \n'
p482
sS'-cluster:radius'
p483
S'<float>Cluster radius\n'
p484
sS'-cluster:exclude_res'
p485
S'<int> [<int> <int> ..]Exclude residue numbers from structural comparisons              \n'
p486
sS'-score:weights'
p487
S'Supply a different weights file (default is score12)\n'
p488
sS'-nstruct'
p489
S'Make how many decoys per input structure ?\n'
p490
sS'-run:shuffle'
p491
S'Use shuffle mode\n'
p492
sS'-cluster:limit_clusters'
p493
S'<int>Maximal number of clusters\n'
p494
sS'-cluster:input_score_filter'
p495
S'<float>Ignore structures above certain energy \n'
p496
sS'-in:file:silent'
p497
S'Input silent file\n'
p498
sS'-in:file:fullatom'
p499
S'Read as fullatom input structure\n'
p500
sS'-cluster:limit_total_structures'
p501
S'<int>Maximal number of structures in total\n'
p502
sS'-cluster:gdtmm'
p503
S'Cluster by gdtmm instead of rms (Better for loops?)\n'
p504
sS'-cluster:sort_groups_by_energy'
p505
S'Sort clusters by energy.\n'
p506
sS'-cluster:limit_cluster_size'
p507
S'<int>Maximal cluster size\n'
p508
ssS'minimize_with_cst'
p509
(dp510
S'-ignore_unrecognized_res'
p511
S'Use this to keep from crashing\n'
p512
sS'-score:patch'
p513
S'Use minirosetta_database/scoring/weights/score12.wts_patch\n'
p514
sS'>ddgLog.log'
p515
S"Use this to obtain a log file.  Use log file with ./convert_to_cst_file.sh to obtain constraints for crystal structure for DDG Monomer!'\n"
p516
sS'-ddg::harmonic_ca_tether 0.5'
p517
S'Use .5\n'
p518
sS'-score:weights standard'
p519
S'Use standard\n'
p520
sS'-ddg::out_pdb_prefix'
p521
S'Prefix for out structures: Use min_cst_0.5\n'
p522
sS'-ddg::constraint_weight 1.0'
p523
S'Use 1.0\n'
p524
sS'-ddg::sc_min_only false'
p525
S'Use false\n'
p526
sS'-in:file:fullatom'
p527
S'Should be fullatom representation - use this\n'
p528
sS'-fa_max_dis 9.0'
p529
S'Use 9.0\n'
p530
sS'-in:file:l'
p531
S'Application only takes a list of PDB structures...\n'
p532
ssS'AbinitioRelax'
p533
(dp534
S'-kill_hairpins'
p535
S'Setup hairpin killing in score (kill hairpin file or psipred file). This option is useful for all-beta or alpha-beta proteins with predicted strands adjacent in sequence since hairpins are often sampled too frequently.\n'
p536
sS'-out:file:silent'
p537
S'Use silent file output, use filename after this flag, default=default.out\n'
p538
sS'-psipred_ss2'
p539
S'psipred_ss2 secondary structure definition file (required for -use_filters)\n'
p540
sS'-constant_seed'
p541
S'Use a constant seed (1111111 unless specified with -jran)\n'
p542
sS'-relax::fast'
p543
S'Do a fast relax instead of classic\n'
p544
sS'-in:file:native'
p545
S'Native structure (optional)\n'
p546
sS'-seed_offset 10'
p547
S'This value will be added to the random number seed. Useful when using time as seed and submitting many jobs to a cluster.  If jobs are started in the same second they will still have different initial seeds when using a unique offset."\n'
p548
sS'-in:file:frag3'
p549
S'3-residue fragments\n'
p550
sS'-abinitio::rsd_wt_helix 0.5'
p551
S'# to Reweight env,pair,cb for helix residues by this factor (Recommended).\n'
p552
sS'-in:file:frag9'
p553
S'9-residue fragments\n'
p554
sS'-abinitio:relax true'
p555
S'Do a relax after abinitio (abrelax protocol), default=false.\n'
p556
sS'-jran'
p557
S'Specify seed. Should be unique among jobs (requires -constant_seed)\n'
p558
sS'-out:pdb true'
p559
S'Output PDB, default=False\n'
p560
sS'-nstruct'
p561
S'# of output structures\n'
p562
sS'-use_filters true'
p563
S'Use radius of gyration (RG), contact-order, and sheet filters. This option conserves computing by not continuing with refinement if a filter fails. A caveat is that for some sequences, a large  percentage of models may fail a filter. The filters are meant to identify models with non-protein like features. The names of models that fail filters start with F_.\n'
p564
sS'-abinitio::rg_reweight 0.5'
p565
S'# to Reweight contribution of radius of gyration to total score by this scale factor (Recommended).\n'
p566
sS'-out:path'
p567
S'Output Path\n'
p568
sS'-abinitio::rsd_wt_loop 0.5'
p569
S'# to Reweight env,pair,cb for loop residues by this factor (Recommended).\n'
p570
sS'-abinitio::increase_cycles 10'
p571
S'# to increase the number of cycles at each stage by this factor (Recommended).\n'
p572
sS'-in:file:fasta'
p573
S'Protein sequence in fasta format (required if native structure is not provided)\n'
p574
ssS'--Output--'
p575
(dp576
S'-out:file:silent'
p577
S'Use silent file output, use filename after this flag. Default="default.out" [String]\n'
p578
sS'-out:show_accessed_options'
p579
S'At the end of the run show options that has been accessed Default="false". [Boolean]\n'
p580
sS'-unmute all'
p581
S"# unmute all channels - default behavior, the same as just './a.out'\n"
p582
sS'-out:file:o'
p583
S' Output file name. [String]\n'
p584
sS'-out:file:scorefile'
p585
S'Write a scorefile to the provided filename. Default = "default.sc [String]\n'
p586
sS'-mute'
p587
S'Mute specified chanels, specify \'all\' to mute all chanels. Defaule="false" [StringVector]\n'
p588
sS'-out:pdb_gz'
p589
S'Compress (gzip) output pdbs", default="false". [Boolean]\n'
p590
sS'-mute ChannelA ChannelB'
p591
S'# mute ChannelA and  ChannelB\n'
p592
sS'-out:file:fullatom'
p593
S'Enable full-atom output of PDB or centroid structures. Default = "fault" [Boolean]\n'
p594
sS'-mute all -unmute ChannelA ChannelB'
p595
S'# mute all channels, unmute ChannelA and ChannelB\n'
p596
sS'-out:prefix'
p597
S'Prefix for output structure names, like old -series code", Default="". [String]\n'
p598
sS'-out:level 10'
p599
S'# Set priority of output messages to be between 0 and 10\n'
p600
sS'-out:pdb'
p601
S'Output PDBs", default="false". [Boolean]\n'
p602
sS'-out:nstruct'
p603
S'Number of times to process each input PDB", default="1" [Integer]\n'
p604
sS'-out:suffix'
p605
S'Suffix for output structure names Default="". [String]\n'
p606
sS'-out:overwrite'
p607
S"Ignore 'CHECKPOINT' file and the overwrite the PDB file(s). [Boolean]\n"
p608
sS'-out:output'
p609
S'Force outputfiles. Default="false". [Boolean]\n'
p610
sS'-out:user_tag'
p611
S'Add this tag to structure tags: e.g., a process id. Default="". [String]\n'
p612
sS'-out:path:score'
p613
S'Score file output path. [Path]\n'
p614
sS'-out:path:pdb'
p615
S'PDB file output path. [Path]\n'
p616
sS'-unmute'
p617
S'UnMute specified chanels. Default="false" [StringVector]\n'
p618
sS'-out:nooutput'
p619
S'Surpress outputfiles. Default="false". [Boolean]\n'
p620
sS'-chname off'
p621
S'# disable output of channels names\n'
p622
ssS'--Input--'
p623
(dp624
S'-in:file:repair_sidechains'
p625
S"Attempt a repack/minmize to repair sidechain problems. Such as proline geometry and his tautomerization' default = 'false' [Boolean]\n"
p626
sS'-in:file:s'
p627
S'Name(s) of single PDB file(s) to process. [FileVector]\n'
p628
sS'-in:path:database'
p629
S'Database file input search paths", default=[\'~/rosetta_database\']. [PathVector]\n'
p630
sS'-in:file:frag3'
p631
S'Fragments file with residue length of 3 [String]\n'
p632
sS'-in:file:native'
p633
S'Native PDB filename. [File]\n'
p634
sS'-in:file:centroid_input'
p635
S"Enable centroid inputs of PDBs.  default = 'false' [Boolean]\n"
p636
sS'-in:file:frag9'
p637
S'Fragments file with residue length of 9 [String]\n'
p638
sS'-in:file:silent_optH'
p639
S'Call optH when reading a silent file. [Boolean]\n'
p640
sS'-in:file:silent_structure_type'
p641
S"Type of SilentStruct object to use in silent-file input'. Default='protein', [String]\n"
p642
sS'-in:file:silent_list'
p643
S'Silent input filename list(s) - like -l is to -s. [FileVector]\n'
p644
sS'-in:file:residue_type_set'
p645
S"ResidueTypeSet for input files', default = 'fa_standard. [String]\n"
p646
sS'-in:file:silent'
p647
S'Silent input filename(s). [FileVector]\n'
p648
sS'-in:file:fullatom'
p649
S'Enable full-atom input of PDB or centroid structures. [Boolean]\n'
p650
sS'-in:file:native_exclude_res'
p651
S'Residue numbers to be excluded from RMS calculation. [IntegerVector]\n'
p652
sS'-in:file:silent_score_prefix'
p653
S"Prefix that is appended to all scores read in from a silent-file', default='' [String]\n"
p654
sS'-in:file:l'
p655
S'File(s) containing list(s) of PDB files to process. [FileVector]\n'
p656
sS'-in:ignore_unrecognized_res'
p657
S"Do not abort if unknown residues are found in PDB file;  instead, ignore them. default='false' [Boolean].\n"
p658
sS'-in:file:fasta'
p659
S'Fasta-formatted sequence file. [FileVector]\n'
p660
ssS'--Refine--'
p661
(dp662
S'-relax:filter_stage2_beginning'
p663
S"FArelax score filter. default='99999999.00' [Real]\n"
p664
sS'-relax:jump_move'
p665
S"Allow jump to move during relax. default='false' [Boolean]\n"
p666
sS'-relax:stage2_repack_period'
p667
S"Full repack after how many cycles in stage 2. default='100' [Integer]\n"
p668
sS'-relax:fastrelax_rampcycles'
p669
S"Default='5' [Integer]\n"
p670
sS'-relax:fastrelax_start_weight_hi'
p671
S"Higherbound on start weight. default='0.02' [Real]\n"
p672
sS'-relax:stage1_ramp_cycles'
p673
S"Ramp cycles in stage 1. default='8' [Integer]\n"
p674
sS'-relax:ramp_constraints'
p675
S'Ramp constraints during phase1 of relax from full to zero. [Boolean]\n'
p676
sS'-relax:filter_stage2_end'
p677
S"FArelax score filter. default='99999999.00' [Real]\n"
p678
sS'-relax:fastrelax_start_weight_low'
p679
S"Lowerbound on start weight. default='0.02' [Real]\n"
p680
sS'-relax:stage3_cycles'
p681
S"How many stage 3 cycles. (by default its -1 means Nresidues ). default='-1' . [Integer]\n"
p682
sS'-relax:chi_move'
p683
S"Allow sidechain to move during relax. default='true' [Boolean]\n"
p684
sS'-relax:fastrelax_startover'
p685
S"Start from scratch when doing multiple packs in fast relax every this many repeats. default='0' [Integer]\n"
p686
sS'-relax:min_tolerance'
p687
S"Minimizer tolerance. default='0.00025' [Real]\n"
p688
sS'-relax:energycut'
p689
S"Rottrial energycut (per residue!). default='0.01' [Real]\n"
p690
sS'-relax:fast'
p691
S"Do fast repacks only. default='false' [Boolean]\n"
p692
sS'-relax:filter_stage2_quarter'
p693
S"FArelax score filter. default='99999999.00' [Real]\n"
p694
sS'-relax:stage1_ramp_inner_cycles'
p695
S"Inner cycles means how many small shear moves + rottrials. default='1' [Integer]\n"
p696
sS'-relax:fastrelax_repeats'
p697
S"Default='3' [Integer]\n"
p698
sS'-relax:filter_stage2_half'
p699
S"FArelax score filter. default='99999999.00' [Real]\n"
p700
sS'-relax:bb_move'
p701
S"Allow backbone to move during relax. default='true' [Boolean]\n"
p702
sS'-relax:fastrelax_postrelax'
p703
S'Do normal move+rotrial+mc moves after initial repacks. [Boolean]\n'
p704
sS'-relax:stage2_cycles'
p705
S"How many stage 2 cycles. (by default its -1 means Nresidues*4 ). default='-1' [Integer]\n"
p706
ssS'idealize_jd2'
p707
(dp708
S'-coordinate_constraint_weight'
p709
S'Can be used to toggle the coordinate constraint weight. Default value is 0.01.\n'
p710
sS'-atom_pair_constraint_weight .01'
p711
S"Sets an optional non-zero weight for atompair constraints. Using atompair constraints results in a significant slow down, but may lead to better preservation of short-range interactions like hydrogen bonds. A reasonable value for this would be ~0.01, but it's worth experimenting with.\n"
p712
sS'-fast'
p713
S'Fast Protocol: ideal geometry is inserted all at once throughout the structure, which perturbs the conformation leading to a potentially large RMSD to the input structure. Then Rosetta\'s minimizer is called to optimize the dihedral angles and reduce the coordinate deviations from the input structure. In the "slow" protocol (i.e., without the -fast option), ideal geometry is inserted at a single residue at a time, and torsion angles nearby the idealized residue are minimized to bring the structure back into agreement with the starting structure.\n'
p714
ss.