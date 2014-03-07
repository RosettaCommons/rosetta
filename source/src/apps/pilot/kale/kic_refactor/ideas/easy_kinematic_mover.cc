void KinematicMover::apply(core::pose::Pose & pose) {

	// Namespaces and Definitions {{{1
	using protocols::loops::loggers::LoggerOP;
	using namespace std;

	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	using numeric::random::RG;
	using core::id::AtomID;
	last_move_succeeded_=false;

	// make backup copy to restore pose, if all trials fail
	core::pose::Pose start_p = pose;
	
	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	// outputs from loop closure
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, save_t_ang, save_b_len, save_b_ang, R0 (3);
	utility::vector1<conformation::ResidueOP> save_residues;

	Size middle_offset = middle_res_ - start_res_;
	Size seg_len = 1 + end_res_ - start_res_;
	atoms.resize((seg_len + 2) * 3);

	// Handle N-terminal loops {{{1
	if (pose.residue(start_res_).is_lower_terminus()) {
		// make a copy of the lower terminal residue
		conformation::ResidueOP nterm_copy( conformation::ResidueFactory::create_residue(
																						  pose.conformation().residue(start_res_).type(),
																						  pose.conformation().residue(start_res_),
																						  pose.conformation(),
																						  false ) );
		// now create another residue of the same type to prepend
		conformation::ResidueOP pre_nterm( conformation::ResidueFactory::create_residue(
																						pose.residue(start_res_).type() ) );
		// create a new conformation containing the n-term copy
		conformation::Conformation extn;
		extn.append_residue_by_bond( *nterm_copy );
		// attached the pre_nterm residue to the nterm copy
		extn.safely_prepend_polymer_residue_before_seqpos( *pre_nterm, 1, true );
		// set ideal omega angle at junction
		extn.set_torsion( id::TorsionID( 1, id::BB, 3 ),  OMEGA_MEAN_);
		// store the xyz coords for KC from the first 3 atoms of the pre-nterm residue into the first 3 atoms indices
		Size ind=1;
		for (Size j=1; j<=3; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(1).xyz(j).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(1).xyz(j).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(1).xyz(j).z());
			ind++;
		}
	}

	// Handle C-terminal loops {{{1
	if (pose.residue(end_res_).is_upper_terminus()) {
		// make a copy of the upper terminal residue
		conformation::ResidueOP cterm_copy( conformation::ResidueFactory::create_residue(
																						  pose.conformation().residue(end_res_).type(),
																						  pose.conformation().residue(end_res_),
																						  pose.conformation(),
																						  false ) );

		// now create another residue of the same type to append
		conformation::ResidueOP post_cterm( conformation::ResidueFactory::create_residue(
																						 pose.residue(end_res_).type() ) );
		// create a new conformation containing the c-term copy
		conformation::Conformation extn;
		extn.append_residue_by_bond( *cterm_copy );
		// attached the post_cterm residue to the cterm copy
		extn.safely_append_polymer_residue_after_seqpos( *post_cterm, 1, true );
		// set ideal omega angle at junction
		extn.set_torsion( id::TorsionID( 2, id::BB, 3 ),  OMEGA_MEAN_);
		// store the xyz coords for KC from the first 3 atoms of the post-cterm residue into the last 3 atoms indices
		Size ind=atoms.size()-2;
		for (Size j=1; j<=3; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(2).xyz(j).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(2).xyz(j).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(2).xyz(j).z());
			ind++;
		}
	}

	// Copy loop coordinates into matrix structure {{{1
	// Get the current coords of the loop closure segment (if a pivot is terminal we took the coords above so skip)
	Size ind = (pose.residue(start_res_).is_lower_terminus() ? 4 : 1);
	for (Size i =  (pose.residue(start_res_).is_lower_terminus() ? start_res_ : start_res_ - 1);
		 i <= (pose.residue(end_res_).is_upper_terminus() ? end_res_ : end_res_ + 1);
		 i++) {
		conformation::Residue res=pose.residue(i);
		for (Size j=1; j<=3; j++) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
			ind++;
		}
	}

	// More definitions {{{1
	// Choose order to solve half-tangents
	order[1]=1;
	order[2]=2;
	order[3]=3;
	// Set the pivot atoms
	Size pvatom1=5; // second C-alpha
	Size pvatom2=5 + (3 * middle_offset); // middle res C-alpha
	Size pvatom3=(3 * (seg_len+1)) - 1; // second-to-last C-alpha
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;


	// Invoke chainTORS {{{1
	// chainTORS is used to get dt to remove solutions identical to original torsions
	// DJM: it would be quicker to store the last torsions rather than recomputing them,
	// but we use db_ang and db_len below so it's worth it for now. However, we should
	// be able to get the torsions, bond angles, and bond lengths from the atom_tree now.
	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	// Save initial coordinates (why?) {{{1
	// save the previous torsions, lengths and angles to restore if no solution is found
	if (! idealize_loop_first_) {
		save_t_ang.resize(dt_ang.size());
		save_b_ang.resize(db_ang.size());
		save_b_len.resize(db_len.size());
		for (Size i=1; i<=dt_ang.size(); i++) {
			save_t_ang[i] = dt_ang[i];
			save_b_ang[i] = db_ang[i];
			save_b_len[i] = db_len[i];
		}
	}

	// Idealize the loop if requested {{{1
	// Not used with default settings.
	if (idealize_loop_first_) {
		
		// save a backbup of the non-ideal residues in case closure fails -- AS_DEBUG: now trying to simply store the entire pose and write it back
		save_residues.resize(seg_len);
		for (Size seqpos=start_res_, i=1; seqpos<=end_res_; seqpos++, i++) {
			conformation::ResidueOP saveres=pose.conformation().residue( seqpos ).clone();
			save_residues[i]=saveres;
		}

		// set all bond lengths, angles, and omegas for closure to ideal values

		// check sizes to prevent memory crashes:
		for (Size i=1; i<=atoms.size(); i+=3) {
			if( db_ang.size() >= (i) )   db_ang[i]=idl_C_N_CA_;
			if( db_ang.size() >= (i+1) ) db_ang[i+1]=idl_N_CA_C_;
			if( db_ang.size() >= (i+2) ) db_ang[i+2]=idl_CA_C_N_;
			if( db_len.size() >= (i) )   db_len[i]=idl_N_CA_;
			if( db_len.size() >= (i+1) ) db_len[i+1]=idl_CA_C_;
			if( db_len.size() >= (i+2) ) db_len[i+2]=idl_C_N_;
			if( dt_ang.size() >= (i+2) ) dt_ang[i+2]=OMEGA_MEAN_;
		}

		// setting pose omegas here, don't set them later!
		for ( Size i= start_res_; i<= end_res_; ++i ) {
			conformation::idealize_position( i, pose.conformation() ); // this is coupled to above values!
			pose.set_omega(i, OMEGA_MEAN_); 
		}

	}
	// }}}1

	for (Size nits=1; nits <= perturber_->max_sample_iterations(); ++nits ) {

		// Perturb the internal coordinate matrices {{{1
		if( perturber_->perturber_exhausted() ) break;

		perturber_->perturb_chain( pose, dt_ang, db_ang, db_len );

		// Solve the closure problem {{{1
		bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);

		BOOST_FOREACH(LoggerOP logger, loggers_) {
			logger->log_task(pose, "BridgeObjects", nsol != 0, false);
		}

		utility::vector1<Size> pos(nsol);
		for (int i=1; i<=nsol; i++) {
			pos[i]=i;
		}
		
		numeric::random::random_permutation(pos.begin(), pos.end(), RG);
		// }}}1
		
		for (Size i=nsol; i>=1; i--) {
			
			// NaN/Inf filter (probably shouldn't exist) {{{1
			if ( ! ( ( -360.0 <= t_ang[pos[i]][pivots[1]-1] ) && (t_ang[pos[i]][pivots[1]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[1]]   ) && (t_ang[pos[i]][pivots[1]]   <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[2]-1] ) && (t_ang[pos[i]][pivots[2]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[2]]   ) && (t_ang[pos[i]][pivots[2]]   <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[3]-1] ) && (t_ang[pos[i]][pivots[3]-1] <= 360.0) ) ||
				! ( ( -360.0 <= t_ang[pos[i]][pivots[3]]   ) && (t_ang[pos[i]][pivots[3]]   <= 360.0) ) ) {
				TR << "solution " << i << " of " << nsol << " failed NaN / Inf check. Skipping..." << std::endl;
				TR << "failed NaN t_ang matrix: " << std::endl;
				printMatrix(t_ang);
				continue;
			}

			// Ramachandran filter {{{1
			bool rama_ok = perform_rama_check(pose, t_ang[pos[i]], pivots, start_res_, middle_res_, end_res_); // checks only pivot residues
			
			BOOST_FOREACH(LoggerOP logger, loggers_) {
				logger->log_task(pose, "RamaCheck", rama_ok, false);
			}

			if (rama_ok == false) continue;
			
			// Apply solution to pose {{{1
			perturber_->set_pose_after_closure( pose, t_ang[pos[i]], b_ang[pos[i]], b_len[pos[i]], true /* closure successful? */);
			
			insert_sampled_torsion_string_into_taboo_map( torsion_features_string( pose ) );
			
			// Bump filter {{{1
			bool bump_ok = do_hardsphere_bump_check_ && !perform_bump_check(pose, start_res_, end_res_);
			BOOST_FOREACH(LoggerOP logger, loggers_) {
				logger->log_task(pose, "BumpCheck", bump_ok, false);
			}

			if(bump_ok) {
				continue;
			}
			
			// Custom filters {{{1
			if( do_sfxn_eval_every_iteration_ ) (*sfxn_)( pose );

			if( filters_.size() != 0 ){
				bool all_filters_passed( true );
				BOOST_FOREACH(protocols::filters::FilterCOP filter, filters_){
					if( ! filter->apply( pose ) ) {
						all_filters_passed = false;
						break;
					}
				}
				BOOST_FOREACH(LoggerOP logger, loggers_) {
					logger->log_task(pose, "CustomFilters", all_filters_passed, false);
				}
				if ( !all_filters_passed ){
					continue;
				}
			}
			// }}}1
			
			last_move_succeeded_ = true;
			return;
		}
	}

	last_move_succeeded_ = false;
	pose = start_p; 
	return;
}

