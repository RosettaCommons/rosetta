// Crappy doxygen comment {{{1
/// @details Refines the polypeptide segments in the LoopMover_KIC::loops_ 
/// class variable from their current conformations.  Uses Rosetta's all-atom 
/// representation and high-resolution scoring function (with an up-weighted 
/// chain break term). At the beginning of this stage, unless the flag 
/// -loops:fix_natsc has been set, all residues within the neighbor distance of 
/// a loop (defined by -loops:neighbor_dist) are repacked and then subject to 
/// rotamer trials.  The backbones of all loop residues, and the side-chains of 
/// all loop residues and neighbors are then subject to DFPmin. If 
/// -loops:fix_natsc is set, only the loop residues (and not the neighbors) 
/// will be subject to repacking, rotamer trials, and minimization.  
/// Consequently, if this stage has been preceded by 
/// LoopMover_Perturb_KIC:model_loop, and the -loops:fix_natsc flag is omitted, 
/// the side-chains surrounding the loop will be optimized for the perturbed 
/// loop conformation, rather than the starting conformation that preceded the 
/// call to model_loop.
// }}}1

void LoopMover_Refine_KIC::apply(core::pose::Pose & pose) {
	// Define some variables {{{1
	static int cur_struct=0; // for movie output

	// Apply the loops file (presumably specified on the command line) to this 
	// pose.
	resolve_loop_indices( pose );

	using namespace core;
	using namespace optimization;
	using namespace scoring;
	using namespace basic::options;

	// Set the native pose if provided
	core::pose::Pose native_pose;
	if( get_native_pose() ){
		native_pose = *get_native_pose();
	}else{
		native_pose = pose;
	}

	// 'verbose' should be an option, or better yet Tracer's output filtering capabilities should be used instead
	bool const verbose( true );
	bool const local_debug( option[ OptionKeys::loops::debug ].user() );
	bool const local_movie( false );
	std::ofstream loop_outfile; // for movie

	// prepare for movie output if requested
	if (local_movie) {
		std::string outname_base = option[ OptionKeys::loops::output_pdb ]().name();
		std::string outname_prefix = option[ OptionKeys::out::prefix ];
		std::string outname = outname_prefix + outname_base + "_fullatom_movie_" +
			right_string_of(cur_struct,4,'0') + ".pdb";
		loop_outfile.open(outname.c_str(), std::ios::out | std::ios::binary);
	}

	// Set cutpoint variants for correct chainbreak scoring {{{1
	Size const nres( pose.total_residue() );
	utility::vector1< bool > is_loop( nres, false );

	for( Loops::const_iterator it=loops()->begin(), it_end=loops()->end();
		 it != it_end; ++it ) {
		for ( Size i= it->start(); i<= it->stop(); ++i ) {
			is_loop[i] = true;
		}
		Size const loop_cut(it->cut());

		if ( loop_cut != nres ) { //c-terminal loop
			if ( ! pose.residue(loop_cut).has_variant_type(chemical::CUTPOINT_LOWER) )
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, loop_cut );
			if ( ! pose.residue(loop_cut+1).has_variant_type(chemical::CUTPOINT_UPPER) )
				core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, loop_cut+1 );
		}
	}

	// Set the number of iterations to run {{{1
	int const fast = option[OptionKeys::loops::fast];
	int outer_cycles(3);
	if ( option[ OptionKeys::loops::outer_cycles ].user() ) {
		outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		outer_cycles = 3;
	}
	int max_inner_cycles( 200 );
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() ) {
		max_inner_cycles = option[ OptionKeys::loops::max_inner_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		max_inner_cycles = 3;
	}

	int inner_cycles = std::min( Size(max_inner_cycles), fast ? (int)loops()->loop_size() : 10 * (Size)( loops()->loop_size() ) );
	int repack_period = 20; // should be an option
	if ( option[ OptionKeys::loops::repack_period ].user() ) {
		repack_period = option[ OptionKeys::loops::repack_period ]();
	}

	// Setup the score function {{{1
	scoring::ScoreFunctionOP local_scorefxn;
	scoring::ScoreFunctionOP min_scorefxn;
	if ( scorefxn() != 0 ) local_scorefxn = scorefxn()->clone();
	else {
		local_scorefxn = get_fa_scorefxn();
	}

	// For testing JK's new solvation term. The exact form isn't yet differentiable
	min_scorefxn = local_scorefxn->clone();
	if ( min_scorefxn->get_weight( core::scoring::occ_sol_exact ) > 0.0 ) { // non-default
		min_scorefxn->set_weight( core::scoring::fa_sol, 0.65 );
		min_scorefxn->set_weight( core::scoring::occ_sol_exact, 0.0 );
	}

	loop_mover::loops_set_chainbreak_weight( scorefxn(), 1 );
	loop_mover::loops_set_chainbreak_weight( min_scorefxn, 1 );
	loop_mover::loops_set_chainbreak_weight( local_scorefxn, 1 ); // AS April 25, 2013 -- the local_scorefxn didn't have a chainbreak weight, restoring that
	
	// AS: add rama2b weight in case we're sampling with rama2b
	if ( option[ OptionKeys::loops::kic_rama2b ]() ) { // non-default
		min_scorefxn->set_weight( rama2b, min_scorefxn->get_weight( rama ) );
		min_scorefxn->set_weight( rama, 0);
		
		local_scorefxn->set_weight( rama2b, local_scorefxn->get_weight( rama ) );
		local_scorefxn->set_weight( rama, 0);
	}
	
	// Setup the parameters of the Monte Carlo simulation {{{1
	float const init_temp( option[ OptionKeys::loops::refine_init_temp ]() );
	float const	final_temp( option[ OptionKeys::loops::refine_final_temp ]() );
	float const gamma = std::pow( (final_temp/init_temp), 1.0f/(outer_cycles*inner_cycles) );
	float temperature = init_temp;
	protocols::moves::MonteCarlo mc( pose, *local_scorefxn, temperature );

	// Setup the minimizer {{{1
	// AS Feb 6 2013: rewriting the minimizer section to use the MinMover, which allows seamless integration of cartesian minimization
	protocols::simple_moves::MinMoverOP min_mover;
	
	const std::string min_type = "lbfgs_armijo_nonmonotone";
	core::Real dummy_tol( 0.001 ); 
	bool use_nblist( true ), deriv_check( false ), use_cartmin ( option[ OptionKeys::loops::kic_with_cartmin ]() ); // true ); // false );
	if ( use_cartmin ) runtime_assert( min_scorefxn->get_weight( core::scoring::cart_bonded ) > 1e-3 ); 
	if ( core::pose::symmetry::is_symmetric( pose ) )  {
		min_mover = new simple_moves::symmetry::SymMinMover(); 
	} else {
		min_mover = new protocols::simple_moves::MinMover();
	}
	min_mover->score_function( min_scorefxn );
	min_mover->min_type( min_type );
	min_mover->tolerance( dummy_tol );
	min_mover->nb_list( use_nblist );
	min_mover->deriv_check( deriv_check );
	min_mover->cartesian( use_cartmin );

	// Set up the packer tasks (rotamer trials & repacking) // {{{1
	// (with design if resfile supplied)
	using namespace pack::task;
	if ( task_factory == 0 ) {
		task_factory = new TaskFactory;
		// TaskOperations replace the following kind of code:
		// base_packer_task->initialize_from_command_line().or_include_current( true );
		task_factory->push_back( new operation::InitializeFromCommandline );
		task_factory->push_back( new operation::IncludeCurrent );
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			// Note - resfile is obeyed, so use NATAA as default to maintain protocol behavior
			task_factory->push_back( new core::pack::task::operation::ReadResfile );
			redesign_loop_ = true;
		}
	}

	PackerTaskOP repack_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
	repack_packer_task->set_bump_check( true );
	pack::task::PackerTaskOP rottrials_packer_task;
	if ( !option[ OptionKeys::packing::resfile ].user() && !redesign_loop_ ) {
		// Not designing -- just repack
		repack_packer_task->restrict_to_repacking();
		tr() << "Not designing" << std::endl;
	}
	else tr() << "Activating design" << std::endl;

	// setting redes loop, but no resfile specified. all non-loop positions only repack. loop positions can design.
	if( redesign_loop_ && !option[ OptionKeys::packing::resfile ].user() ) {
		tr() << "Auto-setting loop design for residues:";
		for( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if( !is_loop[i] ) repack_packer_task->nonconst_residue_task( i ).restrict_to_repacking();
			else tr() << " " << i;
		}
	tr() << std::endl;
	}

	// Setup the kinematic closure move itself {{{1
	loop_closure::kinematic_closure::KinematicMover myKinematicMover;

	// AS Oct 3, 2012: create appropriate perturber here, depending on input flags
	if (option[ OptionKeys::loops::vicinity_sampling ]()) {
		loop_closure::kinematic_closure::VicinitySamplingKinematicPerturberOP perturber =
		new loop_closure::kinematic_closure::VicinitySamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );
		myKinematicMover.set_perturber( perturber );
		perturber->set_degree_vicinity( option[ OptionKeys::loops::vicinity_degree ]() );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ].user() || basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() ) { // torsion-restricted sampling
		std::string torsion_bins = basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ]();
		
		// derive torsion string from native/input pose, if requested -- warning: this overwrites the externally provided one
		if ( basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() )
			torsion_bins = torsion_features_string( native_pose ); 
		
		loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturberOP perturber = new loop_closure::kinematic_closure::TorsionRestrictedKinematicPerturber ( &myKinematicMover, torsion_bins );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() ); // to make the code cleaner this should really be a generic function, even though some classes may not use it
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() && basic::options::option[ basic::options::OptionKeys::loops::taboo_in_fa ]() ) {  // TabooSampling with rama2b (neighbor-dependent phi/psi lookup) -- note that Taboo Sampling will only be active during the first half of the full-atom stage, after that we follow the energy landscape / Monte Carlo
		loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturberOP 
		perturber =
        new loop_closure::kinematic_closure::NeighborDependentTabooSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::taboo_in_fa ] ) { // TabooSampling (std rama) -- note that Taboo Sampling will only be active during the first half of the full-atom stage, after that we follow the energy landscape / Monte Carlo
		loop_closure::kinematic_closure::TabooSamplingKinematicPerturberOP perturber = new loop_closure::kinematic_closure::TabooSamplingKinematicPerturber ( &myKinematicMover );
		
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() ); 
		myKinematicMover.set_perturber( perturber );
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() ) { // neighbor-dependent phi/psi selection
		loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturberOP 
		perturber =
        new loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
		myKinematicMover.set_perturber( perturber );        
	} else { // standard KIC
		loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP perturber =
		new loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber( &myKinematicMover );
		perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );  // default = vary angles
		myKinematicMover.set_perturber( perturber );
	}

	myKinematicMover.set_vary_bondangles( true ); // why is this hard-coded?
	myKinematicMover.set_sample_nonpivot_torsions( option[ OptionKeys::loops::nonpivot_torsion_sampling ]());
	myKinematicMover.set_rama_check( true );
	Size kic_start, kic_middle, kic_end; // three pivot residues for kinematic loop closure

	// Tweak the packer setup {{{1
	// AS: adding option to change the number of rotamer trials -- just one (as in the current implementation) may be way too little
	core::Size num_rot_trials = Size( option[ OptionKeys::loops::kic_num_rotamer_trials ]() );

	// Tweak the score function setup {{{1
	// AS: setting up weights for ramping rama[2b] and/or fa_rep
	core::Real orig_local_fa_rep_weight = local_scorefxn->get_weight( fa_rep );
	core::Real orig_min_fa_rep_weight = min_scorefxn->get_weight( fa_rep );
	
	core::Real orig_local_rama_weight = local_scorefxn->get_weight( rama );
	core::Real orig_min_rama_weight = min_scorefxn->get_weight( rama );
	
	core::Real orig_local_rama2b_weight = local_scorefxn->get_weight( rama2b );
	core::Real orig_min_rama2b_weight = min_scorefxn->get_weight( rama2b );
	
	if ( basic::options::option [ basic::options::OptionKeys::loops::ramp_fa_rep ]() ) { // default: don't ramp
		local_scorefxn->set_weight( fa_rep, orig_local_fa_rep_weight/(outer_cycles + 1) );
		min_scorefxn->set_weight( fa_rep, orig_min_fa_rep_weight/(outer_cycles + 1) );
	}
	//AS: ramp rama -- currently only in refine mode
	if ( basic::options::option [ basic::options::OptionKeys::loops::ramp_rama ]() ) { // default: don't ramp
		local_scorefxn->set_weight( rama, orig_local_rama_weight/(outer_cycles + 1) );
		min_scorefxn->set_weight( rama, orig_min_rama_weight/(outer_cycles + 1) );
		
		// ramp rama2b instead if we're using that for sampling		
		if (option[ OptionKeys::loops::kic_rama2b ]() ) {
			local_scorefxn->set_weight( rama2b, orig_local_rama2b_weight/(outer_cycles + 1) );
			min_scorefxn->set_weight( rama2b, orig_min_rama2b_weight/(outer_cycles + 1) );
		}
	}

	// Perform the initial repack {{{1
	utility::vector1<bool> allow_sc_move_all_loops( nres, false );
	(*local_scorefxn)(pose); // update 10A nbr graph, silly way to do this
	// here we'll optimize all side-chains within neighbor_dist_ of any loop (unless fix_natsc_ is true)
	select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);
	core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );
	repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
	core::pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task );

	// setup rottrials packer task here to ensure we're using current sequence if design was active
	rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
	rottrials_packer_task->restrict_to_repacking();
	rottrials_packer_task->set_bump_check( true );
	rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );

	// Minimize after the initial repack {{{1
	std::string move_type = "repack";
	pose.update_residue_neighbors(); // to update 10A nbr graph
	kinematics::MoveMap mm_all_loops; // DJM tmp
	loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
	if(flank_residue_min_){add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
	MoveMapOP mm_all_loops_OP = new kinematics::MoveMap( mm_all_loops ); // AS, required for using MinMover
	min_mover->movemap( mm_all_loops_OP );
	min_mover->score_function( min_scorefxn ); // AS: needs to be adapted in case we ramp any weights (or does this happen automatically through pointers?)
	min_mover->apply( pose );

	// Decide whether or not to accept the minimization? {{{1
	mc.boltzmann( pose, move_type );
	mc.show_scores();
	if ( redesign_loop_ ) {
		tr() << "Sequence after design step: "
		<< pose.sequence() << std::endl;
	}

	// Setup move maps for each loop {{{1
	// Get a vector of move_maps and allow_sc_move vectors, one for each loop. This way if we have multiple loops
	// we can optimize only the side-chains around the loop selected in the inner_cycle
	utility::vector1< kinematics::MoveMap > move_maps ( loops()->size() );
	utility::vector1< utility::vector1< bool > > allow_sc_vectors ( loops()->size() );
	update_movemap_vectors( pose, move_maps );
	update_allow_sc_vectors( pose, allow_sc_vectors );

	// Movie output {{{1
	if (local_movie) {
		// this assumes there is only one loops_.
		Loops::const_iterator one_loop( loops()->begin() );
		Size begin_loop=one_loop->start();
		Size end_loop=one_loop->stop();
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(end_loop - begin_loop + 3);
		for (Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++) {
			indices[j]=i;
		}
		pose.dump_pdb(loop_outfile, indices, "initial_repack");
		loop_outfile << "ENDMDL" << std::endl;
	}
	// }}}1

	core::Size num_kinematic_trials = 2;

	for (int i=1; i<=outer_cycles; ++i) {
		// Ramp the score function {{{1
		// increase CHAINBREAK weight and update monte carlo
		loop_mover::loops_set_chainbreak_weight( scorefxn(), i );
		loop_mover::loops_set_chainbreak_weight( min_scorefxn, i );
		loop_mover::loops_set_chainbreak_weight( local_scorefxn, i ); // AS April 25, 2013
				
		// AS: try ramping the fa_rep over the outer cycles
		if ( option [ OptionKeys::loops::ramp_fa_rep ]() ) { 
			local_scorefxn->set_weight( fa_rep, orig_local_fa_rep_weight/(outer_cycles - i + 1) );
			min_scorefxn->set_weight( fa_rep, orig_min_fa_rep_weight/(outer_cycles - i + 1) );
		}
		
		if ( option [ OptionKeys::loops::ramp_rama ]() ) { 
			local_scorefxn->set_weight( rama, orig_local_rama_weight/(outer_cycles - i + 1) );
			min_scorefxn->set_weight( rama, orig_min_rama_weight/(outer_cycles - i + 1) );
			
			if (option[ OptionKeys::loops::kic_rama2b ]() ) {
				local_scorefxn->set_weight( rama2b, orig_local_rama2b_weight/(outer_cycles - i + 1) );
				min_scorefxn->set_weight( rama2b, orig_min_rama2b_weight/(outer_cycles - i + 1) );
				
			}
		}
		
		// Turn off taboo sampling if more than half way done {{{1
		// more than 50% done? -> switch off Taboo Sampling (if applicable)
		if ( i > outer_cycles/2 && option[ OptionKeys::loops::taboo_in_fa ]() ) {
			// if rama2b is active, generate a new rama2b perturber, otherwise use a standard torsion sampling perturber
			if ( basic::options::option[ basic::options::OptionKeys::loops::kic_rama2b ]() ) { // neighbor-dependent phi/psi selection
				loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturberOP 
				perturber =
				new loop_closure::kinematic_closure::NeighborDependentTorsionSamplingKinematicPerturber( &myKinematicMover );
				perturber->set_vary_ca_bond_angles( ! option[ OptionKeys::loops::fix_ca_bond_angles ]() );
				myKinematicMover.set_perturber( perturber );        
			} else { // standard KIC
				loop_closure::kinematic_closure::TorsionSamplingKinematicPerturberOP perturber =
				new loop_closure::kinematic_closure::TorsionSamplingKinematicPerturber( &myKinematicMover );
				perturber->set_vary_ca_bond_angles(  ! option[OptionKeys::loops::fix_ca_bond_angles ]()  );
				myKinematicMover.set_perturber( perturber );
			}
		}

		// Finish ramping the score function {{{1
		mc.score_function( *local_scorefxn );

		// Optionally recover the lowest scoring pose {{{1
		if ( recover_low_ ) {
			mc.recover_low( pose );
		}
		// }}}1

		for ( int j=1; j<=inner_cycles; ++j ) {
			// Ramp the temperature {{{1
			temperature *= gamma;
			mc.set_temperature( temperature );

			// Pick a random loop {{{1
			Size loop_ind = RG.random_range(1, loops()->size());
			Loops one_loop;
			one_loop.add_loop( ( *loops() )[ loop_ind ] );

			// get loop endpoints
			Size begin_loop=one_loop.begin()->start();
			Size end_loop=one_loop.begin()->stop();
			
			// AS -- for restricted torsion bin sampling, the mover needs to know 
			// about the start of the defined loop, not just the segment that is 
			// sampled in a given move
			myKinematicMover.set_loop_begin_and_end( begin_loop, end_loop ); 
			
			// Update the sequence if necessary {{{1
			// for Taboo Sampling we need to update the sequence here every time, in 
			// case it changes (e.g. when modeling multiple loops)
			utility::vector1< core::chemical::AA > loop_sequence;
			loop_sequence.resize(0); // make sure it is empty
			for (core::Size cur_res = begin_loop; cur_res <= end_loop; cur_res++) {
				loop_sequence.push_back(pose.aa(cur_res)); 
			}
			myKinematicMover.update_sequence( loop_sequence );
			
			// Setup the packers for the chosen loop {{{1
			// get the movemap for minimization and allowed side-chains for rottrials for this loop
			kinematics::MoveMap cur_mm = move_maps[ loop_ind ];
			MoveMapOP cur_mm_OP = new kinematics::MoveMap( cur_mm );
			utility::vector1<bool> cur_allow_sc_move = allow_sc_vectors[ loop_ind ];
			rottrials_packer_task->restrict_to_residues( cur_allow_sc_move );
			// }}}1

			for (core::Size kinematic_trial = 1; kinematic_trial <= num_kinematic_trials; kinematic_trial++) {
				// Pick pivot residues {{{1
				// AS: the previous implementation had a "history bias" towards the 
				// N-terminus of the loop, as the start pivot can be anywhere between 
				// begin_loop and end_loop-2, while the choice of the end pivot depends 
				// on the start pivot 
				if ( option[ OptionKeys::loops::legacy_kic ]() || j % 2 == 0 ) {
					kic_start = RG.random_range(begin_loop,end_loop-2);
					// choose a random end residue so the length is >= 3, <= min(loop_end, start+maxlen)
					kic_end = RG.random_range(kic_start+2, std::min((kic_start+max_seglen_ - 1), end_loop));
					Size middle_offset = (kic_end - kic_start) / 2;
					kic_middle = kic_start + middle_offset;
				} else {
					kic_end = RG.random_range(begin_loop+2,end_loop);
					kic_start = RG.random_range(std::max((kic_end - std::min(max_seglen_, kic_end) + 1), begin_loop), kic_end-2);
					Size middle_offset = (kic_end - kic_start) / 2;
					kic_middle = kic_start + middle_offset;
				}

				// Make a kinematic closure move {{{1
				myKinematicMover.set_pivots(kic_start, kic_middle, kic_end);
				myKinematicMover.set_temperature(temperature);
				
				core::pose::Pose last_accepted_pose = pose; // backup in case of undetected chain breaks
				myKinematicMover.apply( pose );
				// }}}1

				if ( myKinematicMover.last_move_succeeded() ) {
					// Deal with symmetry {{{1
					//fpd symmetrize mm_all_loops
					if ( core::pose::symmetry::is_symmetric( pose ) )  {
						core::pose::symmetry::make_residue_mask_symmetric( pose, cur_allow_sc_move );
					}

					// Optionally restrict which residues are repacked {{{1
					if ( optimize_only_kic_region_sidechains_after_move_ ) {
						set_rottrials_from_kic_segment( pose, rottrials_packer_task, kic_start, kic_end ); // for rottrials
						set_movemap_from_kic_segment( pose, cur_mm, kic_start, kic_end ); // for minimization
					}
					// }}}1
				
					if (not option[OptionKeys::loops::legacy_kic]() and (j % repack_period) == 0) { 
						// Repacking (all loops) {{{1
						update_movemap_vectors( pose, move_maps );
						update_allow_sc_vectors( pose, allow_sc_vectors );
						
						// the repack/design and subsequent minimization within this main_repack_trial block apply
						// to all loops (and their neighbors, if requested)
						loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
						if(flank_residue_min_){add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
						select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);
						
						core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );  //fpd symmetrize res mask -- does nothing if pose is not symm
						repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task ); // design here if resfile supplied
						if ( verbose ) tr() << "energy after design: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
						if ( redesign_loop_ ) { // need to make new rottrials packer task with redesigned sequence
							rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
							rottrials_packer_task->restrict_to_repacking();
							rottrials_packer_task->set_bump_check( true );
							rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
							if ( verbose ) tr() << "energy after design repack: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
						}
						
						// Minimization (optional) {{{1
						// AS: note that this may not be necessary here, as we'll minimize 
						// after the rotamer trials step anyway
						if ( min_after_repack_ ) {
							if ( core::pose::symmetry::is_symmetric( pose ) )  {
								//fpd  minimizing with the reduced movemap seems to cause occasional problems
								//     in the symmetric case ... am looking into this
								loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_ );
								if(flank_residue_min_){add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
							}
							min_mover->movemap( mm_all_loops_OP ); // will symmetry changes automatically be transferred through the pointer?
							min_mover->score_function( min_scorefxn ); 
							min_mover->apply( pose );
						}
						
						// Accept or reject the move {{{1
						std::string move_type;
						if ( redesign_loop_ ) move_type = "repack+design";
						else move_type = "repack";
						mc.boltzmann( pose, move_type );
						mc.show_scores();
						if ( redesign_loop_ ) {
							tr() << "Sequence after design step: "
							<< pose.sequence() << std::endl;
						}
						if ( verbose ) tr() << "energy after repack: " << (*local_scorefxn)(pose) << std::endl;
						// }}}1
					} 

					// Rotamer trials {{{1
					// AS: first do repacking (if applicable), then rotamer trials -- we've shown that this increases recovery of native side chain conformations
					for (Size i = 1; i <= num_rot_trials; i++) {
						pack::rotamer_trials( pose, *local_scorefxn, rottrials_packer_task );
						pose.update_residue_neighbors(); // to update 10A nbr graph
					}
					
					// Minimization {{{1
					// AS: for non-legacy-KIC one might consider putting this into an ELSE branch -- however, min_after_repack_ defaults to false, in which case we wouldn't have minimization here, which could strongly affect the likelihood of acceptance
					min_mover->movemap( cur_mm_OP );
					min_mover->score_function( min_scorefxn );
					min_mover->apply( pose );

					// Accept or reject the move {{{1
					std::stringstream k_trial;
					k_trial << kinematic_trial;
					std::string move_type = "kic_refine_r" + k_trial.str();
					bool accepted = mc.boltzmann( pose, move_type );

					if (accepted) {
						tr() << "RMS to native after accepted kinematic round " << kinematic_trial << " move on loop "
						   << loops()->size() + 1 - loop_ind << ": " // reverse the order so it corresponds with the loop file
						   << loop_rmsd(pose, native_pose, *loops()) << std::endl;
						if (local_movie) {
							loop_outfile << "MODEL" << std::endl;
							utility::vector1<Size> indices(end_loop - begin_loop + 3);
							for (Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++) {
								indices[j]=i;
							}
							pose.dump_pdb(loop_outfile, indices, "refine_r" + k_trial.str());
							loop_outfile << "ENDMDL" << std::endl;
						}
						if ( verbose ) tr() << "energy after accepted move: " << (*local_scorefxn)(pose) << std::endl;
					}
					// }}}1
				}
			}
			// Always repack at the end of an outer cycle.
			if ( ( option[ OptionKeys::loops::legacy_kic ]() and (j % repack_period) == 0 ) or j == inner_cycles ) {
				// Repacking {{{1
				{ // Open up a new scope so we can reuse variable names.  This code is 
					// disgusting.

					// DJM: we have found that updating the move maps and rottrials/ 
					// repack sets once per repack period gives better performance than 
					// updating them after every move, so we do so here for subsequent 
					// cycles until the next repack period
					update_movemap_vectors( pose, move_maps );
					update_allow_sc_vectors( pose, allow_sc_vectors );

					// the repack/design and subsequent minimization within this main_repack_trial block apply
					// to all loops (and their neighbors, if requested)
					loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_);
					if(flank_residue_min_){add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
					select_loop_residues( pose, *loops(), !fix_natsc_, allow_sc_move_all_loops, neighbor_dist_);

					core::pose::symmetry::make_residue_mask_symmetric( pose, allow_sc_move_all_loops );  //fpd symmetrize res mask -- does nothing if pose is not symm
					repack_packer_task->restrict_to_residues( allow_sc_move_all_loops );
					rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );

					pack::pack_rotamers( pose, *local_scorefxn, repack_packer_task ); // design here if resfile supplied
					if ( verbose ) tr() << "energy after design: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
					
					if ( redesign_loop_ ) { // need to make new rottrials packer task with redesigned sequence
						rottrials_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
						rottrials_packer_task->restrict_to_repacking();
						rottrials_packer_task->set_bump_check( true );
						rottrials_packer_task->restrict_to_residues( allow_sc_move_all_loops );
						if ( verbose ) tr() << "energy after design repack: " << (*local_scorefxn)(pose) << std::endl; // DJM remove
					}

				// Minimization (optional) {{{1
					if ( min_after_repack_ ) {
						if ( core::pose::symmetry::is_symmetric( pose ) )  {
							//fpd  minimizing with the reduced movemap seems to cause occasional problems
							//     in the symmetric case ... am looking into this
							loops_set_move_map( pose, *loops(), fix_natsc_, mm_all_loops, neighbor_dist_ );
							if(flank_residue_min_){add_loop_flank_residues_bb_to_movemap(*loops(), mm_all_loops); } // added by JQX
						}
						min_mover->movemap( mm_all_loops_OP );
						min_mover->score_function( min_scorefxn );
						min_mover->apply( pose );
					}

				// Accept or reject the move {{{1
					std::string move_type;
					if ( redesign_loop_ ) move_type = "repack+design";
					else move_type = "repack";
					mc.boltzmann( pose, move_type );
					mc.show_scores();
					if ( redesign_loop_ ) {
						tr() << "Sequence after design step: "
							 << pose.sequence() << std::endl;
					}
					if ( verbose ) tr() << "energy after repack: " << (*local_scorefxn)(pose) << std::endl;
				// Close our disgusting scope.
				}
				// }}}1
			}
		}
	}

	// Pick a pose to return {{{1
	mc.show_counters();
	if ( recover_low_ ) {
		pose = mc.lowest_score_pose();
	}
	else {
		pose = mc.last_accepted_pose(); // I'm forgetting to do this, I think. -KBK
	}

	// Movie output {{{1
	if (local_movie) {
		// this assumes there is only one loops_.
		Loops::const_iterator one_loop( loops()->begin() );
		Size begin_loop=one_loop->start();
		Size end_loop=one_loop->stop();
		loop_outfile << "MODEL" << std::endl;
		utility::vector1<Size> indices(end_loop - begin_loop + 3);
		for (Size i=begin_loop-1, j=1; i<=end_loop+1; i++, j++) {
			indices[j]=i;
		}
		pose.dump_pdb(loop_outfile, indices, "final_refine");
		loop_outfile << "ENDMDL" << std::endl;
	}
	// }}}1
}
