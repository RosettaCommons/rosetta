void KinematicMover::apply( core::pose::Pose & pose ) {
	// Variable definitions {{{1
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
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	// length of the closure chain
	Size seg_len = 1 + end_res_ - start_res_;

	// Count the number of residues in the backbone and resize the atom list to 
	// match this.  There should be one extra residue on each side to establish 
	// the geometric frame.
	atoms.resize(count_bb_atoms(pose, start_res_, end_res_));

	// Handle N-terminal loops {{{1
	// KC requires 3 backbone atoms 1 residue N-terminal to start pivot and 1 residue C-terminal to end pivot
	// so check for terminal pivots and compute coordinates from prepended / appended residues if necessary
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

		// Store the xyz coords for KC from the first 3 atoms of the pre-nterm residue into the first 3 atoms indices
		// Note: if this is a beta-amino acid, it will store the first 4 atoms.
		for (Size j=1, numatoms=count_bb_atoms_in_residue(pose, start_res_), ind=1; j<=numatoms; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(1).xyz( get_bb_atoms_for_residue( extn.residue(1), j ) ).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(1).xyz( get_bb_atoms_for_residue( extn.residue(1), j ) ).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(1).xyz( get_bb_atoms_for_residue( extn.residue(1), j ) ).z());
			ind++;
		}
	}

	// Handle C-terminal loops {{{1
	// check for upper terminus
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
		// store the xyz coords for KC from the first 3 atoms of the post-cterm residue into the last 3 atom indices
		// Note: if this is a beta-amino acid, it will store 4 atoms (N, CA, CM, C) in the last 4 atom indices.
		for (Size j=1, numatoms=count_bb_atoms_in_residue(pose, end_res_), ind=atoms.size()-count_bb_atoms_in_residue(pose, end_res_)+1; j<=numatoms; j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (extn.residue(2).xyz( get_bb_atoms_for_residue( extn.residue(2), j ) ).x());
			atoms[ind][2] = static_cast<Real> (extn.residue(2).xyz( get_bb_atoms_for_residue( extn.residue(2), j ) ).y());
			atoms[ind][3] = static_cast<Real> (extn.residue(2).xyz( get_bb_atoms_for_residue( extn.residue(2), j ) ).z());
			ind++;
		}
	}

	// Fill in the "atoms" matrix {{{1
	// Get the current coords of the loop closure segment (if a pivot is terminal we took the coords above so skip)
	Size ind = (pose.residue(start_res_).is_lower_terminus() ? count_bb_atoms_in_residue(pose, start_res_) + 1 : 1); //If the start res is a terminus, then an identical amino acid has been prepended, so the start atom is 1 plus the number of backbone atoms in the prepended residue.
	for (Size i = (pose.residue(start_res_).is_lower_terminus() ? start_res_ : start_res_ - 1);
		 i <= (pose.residue(end_res_).is_upper_terminus() ? end_res_ : end_res_ + 1);
		 i++) {
		conformation::Residue res=pose.residue(i);

		// Updated by VKM -- keep N,CA,C by default, though the 
		// get_bb_atoms_for_residue function will return whatever atoms should be 
		// kept (e.g. N, CA, CM, C for beta-amino acids).  That is, we're no longer 
		// making the protein-only assumption.
		for (Size j=1; j<=count_bb_atoms_in_residue(pose, i); j++) {
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz( get_bb_atoms_for_residue( res, j ) ).x());
			atoms[ind][2] = static_cast<Real> (res.xyz( get_bb_atoms_for_residue( res, j ) ).y());
			atoms[ind][3] = static_cast<Real> (res.xyz( get_bb_atoms_for_residue( res, j ) ).z());
			ind++;
		}
	}

	// Choose pivot atoms {{{1
	// Get alternative conformations from from kinematic loop closure
	// Choose order to solve half-tangents
	order[1]=1;
	order[2]=2;
	order[3]=3;

	// Set the pivot atoms
	core::Size start_minus_one_bb_atom_count = pose.residue(start_res_).is_lower_terminus() ? count_bb_atoms_in_residue(pose, start_res_) : count_bb_atoms_in_residue(pose, start_res_ - 1); //Number of backbone atoms for the first residue in the segment (start_res_ - 1 or the prepended start_res_ if start_res_ is a terminus)
	core::Size end_plus_one_bb_atom_count = pose.residue(end_res_).is_upper_terminus() ? count_bb_atoms_in_residue(pose, end_res_) : count_bb_atoms_in_residue(pose, end_res_ + 1); //Number of backbone atoms for the last residue in the segment (end_res_ + 1 or the appended end_res_ if end_res_ is a terminus)
	core::Size pvatom1 = start_minus_one_bb_atom_count + 2; // Second backbone atom of start_res_ (CA if alpha or beta-amino acid).
	core::Size pvatom2 = start_minus_one_bb_atom_count + 2;
	for(core::Size ir = start_res_; ir<middle_res_; ++ir) pvatom2 += count_bb_atoms_in_residue(pose, ir); //This will set pvatom2 to be the second backbone atom of middle_res_ (CA if alpha or beta-amino acid).
	core::Size pvatom3 = atoms.size() - end_plus_one_bb_atom_count - count_bb_atoms_in_residue(pose, end_res_) + 2; // Second backbone atom of end_res_ (CA if alpha or beta-amino acid).
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;

	// Fill in "angle/length" matrices {{{1
	// chainTORS is used to get dt to remove solutions identical to original 
	// torsions

	// DJM: it would be quicker to store the last torsions rather than 
	// recomputing them, but we use db_ang and db_len below so it's worth it for 
	// now. However, we should be able to get the torsions, bond angles, and bond 
	// lengths from the atom_tree now.

	// VKM, 26 Aug 2013: As far as I can tell, chainTORS calculates bond lengths, 
	// angles, and torsions from the stored backbone atom coordinates, and is 
	// completely agnostic to the relationship between this list of atoms and any 
	// polymer.  This is probably the best (or, at least, safest) way to do this 
	// if we want kinematic closure to be alignly applicable (i.e. useful for 
	// non-protein backbones, too).
	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);
	// }}}1

	if (idealize_loop_first_) {
		// Set "angles/lengths" to ideal values {{{1
		// set all bond lengths, angles, and omegas for closure to ideal values
		// check sizes to prevent memory crashes:
		for (Size i=1, res=start_res_-1; i<=atoms.size(); res+=1) {
			//Check for whether this is a beta-amino acid:
			bool is_beta_residue = false;
			if(res == start_res_ - 1 && pose.residue(start_res_).is_lower_terminus()) is_beta_residue = is_beta_aminoacid(pose.residue(start_res_));
			else if(res == end_res_ + 1 && pose.residue(end_res_).is_upper_terminus()) is_beta_residue = is_beta_aminoacid(pose.residue(end_res_));
			else is_beta_residue = is_beta_aminoacid(pose.residue(res));

			// Idealizing a beta-amino acid
			if(is_beta_residue) {
				if( db_ang.size() >= (i) )   db_ang[i]=idl_beta_C_N_CA_;
				if( db_ang.size() >= (i+1) ) db_ang[i+1]=idl_beta_N_CA_CM_;
				if( db_ang.size() >= (i+2) ) db_ang[i+2]=idl_beta_CA_CM_C_;
				if( db_ang.size() >= (i+3) ) db_ang[i+3]=idl_beta_CM_C_N_;
				if( db_len.size() >= (i) )   db_len[i]=idl_beta_N_CA_;
				if( db_len.size() >= (i+1) ) db_len[i+1]=idl_beta_CA_CM_;
				if( db_len.size() >= (i+2) ) db_len[i+2]=idl_beta_CM_C_;
				if( db_len.size() >= (i+3) ) db_len[i+3]=idl_beta_C_N_;
				if( dt_ang.size() >= (i+3) ) dt_ang[i+3]=( core::chemical::is_D_aa(pose.residue(res).aa()) ? -1.0 : 1.0 ) * OMEGA_MEAN_;				
			}
			// Default case: alpha-amino acid residue
			else {
				if( db_ang.size() >= (i) )   db_ang[i]=idl_C_N_CA_;
				if( db_ang.size() >= (i+1) ) db_ang[i+1]=idl_N_CA_C_;
				if( db_ang.size() >= (i+2) ) db_ang[i+2]=idl_CA_C_N_;
				if( db_len.size() >= (i) )   db_len[i]=idl_N_CA_;
				if( db_len.size() >= (i+1) ) db_len[i+1]=idl_CA_C_;
				if( db_len.size() >= (i+2) ) db_len[i+2]=idl_C_N_;
				if( dt_ang.size() >= (i+2) ) dt_ang[i+2]=OMEGA_MEAN_;
			}

			// Incrementing the atom counter:
			if(res == start_res_ - 1) i += start_minus_one_bb_atom_count; // The padding residue at the start
			else if (res == end_res_ + 1) i += end_plus_one_bb_atom_count; // The padding residue at the end
			else i += count_bb_atoms_in_residue(pose, res); // The loop residues
		}

		// setting pose omegas here, don't set them later!
		for ( Size i= start_res_; i<= end_res_; ++i ) {
			conformation::idealize_position( i, pose.conformation() ); // this is coupled to above values!
			if( is_beta_aminoacid(pose.residue(i))) pose.conformation().set_torsion( id::TorsionID( i, id::BB, 4 ), OMEGA_MEAN_ ); //set_omega doesn't work properly for beta-amino acids (it sets the third torsion rather than the fourth).
			else if (core::chemical::is_D_aa(pose.residue(i).aa())) pose.set_omega(i, -1.0*OMEGA_MEAN_);
			else pose.set_omega(i, OMEGA_MEAN_); //Default case for alpha-amino acids -- set the torsion angle with set_omega.
		}
	// }}}1
	}

	for (Size nits=1; nits <= perturber_->max_sample_iterations(); ++nits ) {
		// Make a perturbation {{{1
		// make sure the perturber can still generate solutions
		if( perturber_->perturber_exhausted() ) break;

		perturber_->perturb_chain( pose, dt_ang, db_ang, db_len );

		// Solve the closure problem {{{1
		// I THINK this doesn't need to know about standard vs. nonstandard 
		// backbones -- i.e. it just takes a list of bond angles, torsion angles, 
		// and lengths and does its thing.
		bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);

		// Randomly shuffle the solutions {{{1
		utility::vector1<Size> pos(nsol);
		for (int i=1; i<=nsol; i++) {
			pos[i]=i;
		}
		
		// Put the solutions in random order, since the first one encountered that 
		// passes the filters is accepted.
		numeric::random::random_permutation(pos.begin(), pos.end(), RG);
		// }}}1
		
		for (Size i=nsol; i>=1; i--) {
			// Sanity filter {{{1
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

			// Rama filter {{{1
			if (!perform_rama_check(pose, t_ang[pos[i]], pivots, start_res_, middle_res_, end_res_)) {
				continue;
			}

			// Apply solution to pose {{{1
			// place the solution into the pose and bump check+eventual filters
			// set the torsions
			perturber_->set_pose_after_closure( pose, t_ang[pos[i]], b_ang[pos[i]], b_len[pos[i]], true /* closure successful? */);

			// AS: fetch & set the torsion string for Taboo Sampling
			insert_sampled_torsion_string_into_taboo_map( torsion_features_string( pose ) ); //VKM, 27 Aug 2013: This won't work properly with beta-amino acids, but won't crash, either...
			
			// Bump filter {{{1
			if( do_hardsphere_bump_check_ && !perform_bump_check(pose, start_res_, end_res_) ){
				continue;
			}
			
			// Evaluate the score function {{{1
			if( do_sfxn_eval_every_iteration_ ) (*sfxn_)( pose );

			// Custom filters {{{1
			if( filters_.size() != 0 ){
				bool all_filters_passed( true );
				foreach(protocols::filters::FilterCOP filter, filters_){
					if( ! filter->apply( pose ) ) {
						all_filters_passed = false;
						break;
					}
				}
				if ( !all_filters_passed ){
					continue;
				}
			}
			// }}}1

			last_move_succeeded_ = true;
			return;
		}

	// Revert to the input pose {{{1
	// prepare to exit here if no solutions
	last_move_succeeded_ = false;

	// AS_DEBUG -- replace pose by stored pose in case the loop closure trials didn't work
	pose = start_p;
	// }}}1
}
