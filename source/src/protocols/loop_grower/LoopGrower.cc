
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Brandon Frenz
/// @author Frank DiMaio

#include <protocols/loop_grower/LoopGrower.hh>
#include <protocols/loop_grower/util.hh>

#include <iostream>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataMap.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/util.hh>
#include <core/fragment/Frame.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/CircularSplineFunc.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/USOGFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/sequence/AnnotatedSequence.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>

#include <protocols/comparative_modeling/coord_util.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/util.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>
#include <protocols/loop_grower/SheetSampler.hh>
#include <protocols/loop_grower/DensSkeleton.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>
#include <boost/algorithm/string/predicate.hpp>

namespace protocols {
namespace loop_grower {

using namespace core;

//static numeric::random::RandomGenerator RG(8403179);

// stupid function
// until we get better dist distribs
double gapdist( int n ) {
	switch (n) {
	case 0 : return 0.5;
	case 1 : return 5.0; // return 5.0
	case 2 : return 8.0; // return 8.5
	case 3 : return 11.0;
	case 4 : return 14.0;
	case 5 : return 17.0;
	case 6 : return 20.0;
	case 7 : return 23.0;
	case 8 : return 26.0;
	case 9 : return 29.0;
	case 10 : return 32.0;
	default : return 999.0;
	}
}


protocols::moves::MoverOP LoopGrower::clone() const { return protocols::moves::MoverOP( new LoopGrower( *this ) ); }
protocols::moves::MoverOP LoopGrower::fresh_instance() const { return protocols::moves::MoverOP( new LoopGrower ); }

void
LoopGrower::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::datacache;

	bool is_nterm = (resstart_ == 0 || pose.fold_tree().is_cutpoint(resstart_-1));
	bool is_cterm = (resstop_ == 0 || pose.fold_tree().is_cutpoint(resstop_));
	bool is_cutpoint = pose.fold_tree().is_cutpoint(resstart_);
	bool dont_cut_jump = false;
	lower_fasta_ = loop_.start();
	upper_fasta_ = loop_.stop();

	int eff_resstart=resstart_, eff_resstop=resstop_, torsionrangelo=0, torsionrangehi=0;
	Size fasta_range_low, fasta_range_hi, pose_range_low, pose_range_hi, cutpoint_pose, cutpoint_fasta;
	if ( is_nterm ) eff_resstart = resstop_;
	if ( is_cterm ) eff_resstop  = resstart_;

	// fragments
	Size maxfrag = fragments_[1]->max_frag_length();
	core::fragment::FragSetOP fragset = fragments_[1];
	for ( Size i=1; i<=fragments_.size(); i++ ) {
		fragset = fragments_[i];
		if ( maxfrag > fragset->max_frag_length() ) {
			maxfrag = fragset->max_frag_length();
		}
	}
	runtime_assert( fragmelt_< maxfrag );


	// chop fragments if short loop
	unsigned loopsize = (loop_.stop() - loop_.start());
	if ( (loopsize < maxfrag) && (maxfrag > 2) ) {
		core::fragment::FragSetOP smallfrags = core::fragment::FragSetOP(new core::fragment::ConstantLengthFragSet( 2 ));
		core::fragment::chop_fragments( *fragset, *smallfrags );
		fragments_.erase (fragments_.begin(),fragments_.end());
		fragments_.push_back(smallfrags);
		fragmelt_ = 1;
		fragset = fragments_[1];
		maxfrag = fragset->max_frag_length();
		update_fragment_library_pointers();
		TRACER << "using shortloop logic" << std::endl;
	}
	maxfrag_ = maxfrag;

	if ( coordfile_ != "" ) {
		read_coordfile();
	}

	if ( skeleton_file_ != "" ) {
		core::Real mapreso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
		core::Real gridspace = 1;
		skeleton_.load_density( skeleton_file_, mapreso, gridspace );
	}


	// count the number of residues  left and right before we encounter a cut
	// This segment of code deals with non-standard loops such as the termini and instances in which the initial fragmelt and minmelt conditions cannot
	// be applied. It then attempts to set those parameters to acceptable values.
	Size initial_melt_left = fragmelt_, initial_melt_right = fragmelt_, minmelt_left = minmelt_+fragmelt_, minmelt_right = minmelt_+fragmelt_, totalmelt = minmelt_+fragmelt_ ;
	if ( !is_nterm ) {
		int i=1;
		while ( i<=(int)fragmelt_ && !pose.fold_tree().is_cutpoint(resstart_-i) ) {
			i++;
		}

		initial_melt_left = i;
	}
	if ( !is_nterm ) {
		int i=1;
		while ( i<(int)totalmelt && !pose.fold_tree().is_cutpoint(resstart_-i) ) {
			if ( pose.fold_tree().is_cutpoint(resstart_-i) ) TRACER << " hit cutpoint melt left " << std::endl;
			i++;
		}

		minmelt_left = i;
		if ( minmelt_left != totalmelt || pose.fold_tree().is_cutpoint(resstart_-minmelt_left) ) {
			if ( minmelt_left > 2 ) {
				minmelt_left -= 2;
			}
		}
	}
	if ( !is_cterm ) {
		int i=1;
		while ( i<=(int)fragmelt_ && !pose.fold_tree().is_cutpoint(resstop_+i-1) ) {
			if ( pose.fold_tree().is_cutpoint(resstop_+i-1) ) TRACER << " hit cutpoint melt right " << std::endl;
			i++;
		}
		initial_melt_right = i;
	}
	if ( !is_cterm ) {
		int i=1;
		while ( i<(int)totalmelt && !pose.fold_tree().is_cutpoint(resstop_+i+1) ) {
			i++;
		}
		minmelt_right = i;
		if ( minmelt_right != totalmelt ) {
			if ( totalmelt != 1 ) {
				minmelt_right -= 1;
			}
			dont_cut_jump = true;
		}
		if ( pose.fold_tree().is_cutpoint(resstop_+minmelt_right-1) ) {
			dont_cut_jump = true;
		}
	}
	if ( is_cterm ) {
		initial_melt_right = 0;
		minmelt_right = 0;
	}
	if ( is_nterm ) { // && lower_fasta_ == 1){
		initial_melt_left  = 0;
		minmelt_left = 0;
	}
	if ( is_nterm && lower_fasta_ != 1 ) {
		int i=1;
		while ( i<=(int)fragmelt_ && !pose.fold_tree().is_cutpoint(eff_resstart-1-i) ) {
			i++;
		}
		initial_melt_left = i;
	}
	if ( is_nterm && lower_fasta_ != 1 ) {
		int i=1;
		while ( i<=(int)totalmelt && !pose.fold_tree().is_cutpoint(eff_resstart-1-i) ) {
			i++;
		}

		minmelt_left = i;
		resstart_ = resstop_-1;

		if ( minmelt_left != totalmelt+1 || pose.fold_tree().is_cutpoint(eff_resstart-1-minmelt_left) ) {
			if ( minmelt_left > 1 ) {
				minmelt_left -= 1;
			}
		}
	}

	// shuffle foldtree while maintaining original cuts
	core::kinematics::FoldTree f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	core::kinematics::FoldTree f_new;
	int ncuts = f_in.num_cutpoint();
	if ( !is_cutpoint ) {
		ncuts +=1;
	}

	ObjexxFCL::FArray1D< Size > fcuts ( ncuts );
	for ( int i=1; i<=(int)f_in.num_cutpoint(); ++i ) {
		fcuts(i) = f_in.cutpoint(i);
	}
	if ( !is_cutpoint ) {
		fcuts(ncuts) = resstart_;
	}
	// reshuffle the foldtree
	ObjexxFCL::FArray2D< Size > fjumps( 2, ncuts );
	for ( int i=1; i<=ncuts; ++i ) {
		int cut_i = (int) fcuts(i);

		fjumps(1,i) = cut_i;
		fjumps(2,i) = cut_i+1;

		if ( cut_i==resstart_ ) {
			fjumps(1,i) = cut_i-minmelt_left;
			if ( !is_cterm ) {
				fjumps(2,i) = cut_i+minmelt_right+2;
			} else {
				fjumps(2,i) = cut_i+minmelt_right+1;
			}
		}
		if ( int(fjumps(1,i)) >resstart_-(int)minmelt_left && int(fjumps(1,i))<=resstart_ ) fjumps(1,i)=resstart_-(int)minmelt_left;
		if ( int(fjumps(2,i)) >resstart_-(int)minmelt_left && int(fjumps(2,i))<=resstart_ ) fjumps(2,i)=resstart_-(int)minmelt_left;
		if ( int(fjumps(1,i)) >resstart_ && int(fjumps(1,i))<resstart_+(int)minmelt_right+1 ) fjumps(1,i)=resstart_+(int)minmelt_right;
		if ( int(fjumps(2,i)) >resstart_ && int(fjumps(2,i))<resstart_+(int)minmelt_right+1 ) fjumps(2,i)=resstart_+(int)minmelt_right;

		// fcuts(i) = f_in.cutpoint(i);
		TRACER << "fjumps and fcuts " << i << " " << fjumps(1,i) << "," << fjumps(2,i) << "," << fcuts(i) << std::endl;
	}
	bool valid_tree = f_new.tree_from_jumps_and_cuts( pose.total_residue(), ncuts, fjumps, fcuts );
	runtime_assert( valid_tree );
	f_new.reorder( f_in.root() );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::symmetrize_fold_tree( pose.conformation(), f_new );
	}
	pose.fold_tree( f_new );
	TRACER << pose.fold_tree() << std::endl;

	int tgt_jump = 0;
	if ( !is_nterm && !is_cterm ) {
		for ( int i=1; i<=(int)pose.fold_tree().num_jump() && tgt_jump==0; ++i ) {
			if ( pose.fold_tree().cutpoint_by_jump(i) == Size( resstart_ ) ) tgt_jump=i;
		}
		runtime_assert (tgt_jump != 0);
	}

	//these are needed for sheet sampling
	numjumps_ = pose.num_jump();
	total_residues_ = pose.total_residue();

	// centroid copy of the pose
	protocols::simple_moves::SwitchResidueTypeSetMoverOP tocen = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid")); // tocen("centroid_rot");
	protocols::simple_moves::SwitchResidueTypeSetMoverOP tocenrot = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid_rot")); // tocen("centroid_rot");
	// if you don't want to use cen and cenrot at the same time be sure to change this to only set the cen pose to cenrot

	//Prep the template poses for constraints
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() ) {
		core::pose::make_pose_from_sequence(cen_seqpose_, seq_->ungapped_sequence(), "centroid", true);
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *cen_sf_ );
	}
	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user() ) {
		core::pose::make_pose_from_sequence(fa_seqpose_, seq_->ungapped_sequence(), "fa_standard", true);
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sf_ );
	}

	core::pose::Pose pose_cen = pose;
	tocen->apply(pose_cen);
	if ( cenrot_ ) {
		tocenrot->apply(pose);
	}
	add_user_csts(pose);
	add_user_csts(pose_cen);
	if ( cenrot_ ) {
		(*cenrot_sf_)(pose);
	} else {
		(*sf_)(pose);
	}
	Real hbond_weight = 0;
	core::scoring::ScoreFunctionOP sf_nodens;
	core::scoring::ScoreFunctionOP censf_nodens = cen_sf_->clone();
	if ( cenrot_ ) {
		sf_nodens = cenrot_sf_->clone();
		hbond_weight = cenrot_sf_->get_weight(scoring::hbond_sr_bb);
	} else {
		hbond_weight = cen_sf_->get_weight(scoring::hbond_sr_bb);
		sf_nodens = sf_->clone();
	}
	sf_nodens->set_weight( core::scoring::elec_dens_fast , 0.0 );
	censf_nodens->set_weight( core::scoring::elec_dens_fast , 0.0 );


	if ( !cenrot_ ) {
		startingscore_ = (*sf_nodens)(pose);
	} else {
		startingscore_ = (*sf_nodens)(pose) - nton3_hbond_score(pose)*hbond_weight*0.7;
	}
	censtartingscore_ = (*censf_nodens)(pose_cen) - nton3_hbond_score(pose_cen)*hbond_weight*0.7;
	pose.remove_constraints();
	pose_cen.remove_constraints();


	// res type sets
	core::chemical::ResidueTypeSetCOP restypeset_cen = pose_cen.residue_type_set_for_pose();
	core::chemical::ResidueTypeCOPs restypes_pose_cen = core::pose::residue_types_from_sequence(seq_->sequence(), *restypeset_cen);
	core::chemical::ResidueTypeSetCOP restypeset = pose.residue_type_set_for_pose();
	core::chemical::ResidueTypeCOPs restypes_pose = core::pose::residue_types_from_sequence(seq_->sequence(), *restypeset);

	// lower and upper positions in the pose and fasta
	//    track with insertions
	int lower_pose = resstart_; // the pose, 1 res before cutpoint
	int upper_pose = eff_resstop;  // the start of the upper loop
	int upper_term = eff_resstop;


	// store top N configurations at each expansion
	LoopPartialSolutionStore solutionset( beamwidth_, rmscutoff_, master_beam_width_, master_beam_cutoff_ ), solutionsetlower( beamwidth_, rmscutoff_, master_beam_width_, master_beam_cutoff_ ); // solutionsetlower( beamwidth_, rmscutoff_ );
	LoopPartialSolutionStore filteronlysolutionset( beamwidth_, rmscutoff_, master_beam_width_, master_beam_cutoff_ );
	bool dumpbeam = dumpbeam_;
	if ( rescorebeams_ ) {
		dumpbeam = false;
	}
	solutionset.setfilterparams( fragmelt_, rmswindow_, parallelcount_, beamscorecutoff_, dumperrors_, dumpbeam, writebeams_, clustercheck_, fafilter_, samplesheets_, filterprevious_, checksymm_,
		asymmdump_, dumpfinalbeam_ );
	solutionsetlower.setfilterparams( fragmelt_, rmswindow_, parallelcount_, beamscorecutoff_, dumperrors_, dumpbeam, writebeams_, clustercheck_, fafilter_, samplesheets_, filterprevious_, checksymm_,
		asymmdump_, dumpfinalbeam_ );
	filteronlysolutionset.setfilterparams( fragmelt_, rmswindow_, parallelcount_, beamscorecutoff_, dumperrors_, dumpbeam, writebeams_, clustercheck_, fafilter_, samplesheets_, filterprevious_,
		checksymm_, asymmdump_, dumpfinalbeam_ );

	if ( minmelt_left != 0 ) {
		torsionrangelo = eff_resstart-minmelt_left+1;
		fasta_range_low = lower_fasta_-minmelt_left;
	} else {

		torsionrangelo = eff_resstart;
		if ( torsionrangelo == 1 ) torsionrangelo -= 1;
		fasta_range_low = lower_fasta_;
	}
	pose_range_low = lower_fasta_;
	if ( minmelt_right != 0 ) {
		torsionrangehi = upper_pose+minmelt_right;
		fasta_range_hi = upper_fasta_+minmelt_right;
	} else {
		torsionrangehi = upper_pose;
		fasta_range_hi =  upper_fasta_;
	}
	pose_range_hi = upper_fasta_;
	if ( is_nterm ) {
		torsionrangelo = eff_resstart;
		fasta_range_low = lower_fasta_;
	}
	if ( is_cterm ) {
		torsionrangehi = upper_pose;
		fasta_range_hi = upper_fasta_-1;
		pose_range_hi = upper_fasta_;
	}

	//add empty solution to solution set
	LoopPartialSolution lps( pose, torsionrangelo, torsionrangehi, 0.0);
	lps.set_bonus_score(0);
	solutionset.push( pose, lps, fragmelt_, beamscorecutoff_, 0, true, false, 0.0, 0, 0 );

	//store side chains
	protocols::moves::MoverOP restore_sc;

	//this sets us up to use the second pass filter with cenrot
	if ( fafilter_ && cenrot_ ) {
		cenrotfilter_ = true;
		cenrot_ = 0;
	}
	if ( fafilter_ ) {
		fafilter_pmcycles_ = pack_min_cycles_;
		pack_min_cycles_ = 0;
	}


	bool done = false;
	int cycle=1;
	Size stepcount = 1;
	rmsrangelo_ = loop_.start()-minmelt_left+1;
	rmsrangehi_ = loop_.stop()+minmelt_right-1;
	if ( is_nterm ) rmsrangelo_ = loop_.start();
	if ( is_cterm ) rmsrangehi_ = loop_.stop();

	if ( parametercheck_ ) {
		lowest_ranked_native_ = 0;
		beamscorecutoff_ = 0;
		if ( master_beam_width_ == 1337 && master_beam_cutoff_ == 1337 ) {
			lowest_ranked_native_ = beamwidth_;
		}
	}
	while ( !done ) {


		int n_to_insert_lower, n_to_insert_upper;
		//set grow directions
		bool growlower = true;
		bool growupper = true;
		if ( direction_ == 1 ) {
			growlower = true;
			growupper = false;
		}
		if ( direction_ == 2 ) {
			growlower = false;
			growupper = true;
		}
		if ( direction_ == 3 ) {
			growlower = (numeric::random::uniform() <0.5);
			growupper = !growlower;
		}
		if ( is_cterm ) {
			growlower = true;
			growupper = false;
		}
		if ( is_nterm ) {
			growlower = false;
			growupper = true;
		}

		//this logic deals with cases where the melts exceeds the length between loops
		core::Size insert_lower_melt, insert_upper_melt;
		if ( initial_melt_left >= fragmelt_ ) {
			insert_lower_melt = fragmelt_;
		} else {
			insert_lower_melt = initial_melt_left;
			if ( insert_lower_melt != 0 ) insert_lower_melt = 0;
		}
		if ( initial_melt_right >= fragmelt_ ) {
			insert_upper_melt = fragmelt_;
		} else {
			insert_upper_melt = initial_melt_right;
			if ( insert_upper_melt != 0 ) insert_upper_melt = 0;
		}

		bool acceptlower = true;
		bool update_pose = false;
		if ( filterbeams_ != "" && filterprevious_ ) {
			read_from_disk(filteronlysolutionset, cycle, acceptlower, true);
			utility::vector1<LoopPartialSolution> filtersolutions = filteronlysolutionset.get_solutions();
			solutionset.store_filteronly_solutions(filtersolutions);
			storehi_ = 0;
			storelow_ = 0;
		}
		solutionsetlower = solutionset;

		n_to_insert_lower = (int)maxfrag-insert_lower_melt;
		n_to_insert_upper = (int)maxfrag-insert_upper_melt;
		//turn off dump beam this is done a precaution to avoid accidentally creating many thousands of pdb files
		if ( steps_ != 0 ) {
			if ( dumpbeam_ ) {
				TR << " warning dump beam set to true in a production run. This will produce many thousands of pdbfiles. This behavior is being turned off" << std::endl;
			}
			dumpbeam_ = false;
		}
		//This logic handles setting up the pose and other inputs when doing filtering using the parallelization through the python script
		if ( readbeams_ && steps_ == 0 && (dumpbeam_ || clustercheck_ || rescorebeams_ || nativegrow_ || (coordfile_ != "") || auto_stop_) ) update_pose = true;
		if ( readbeams_ && steps_ == 0 ) {
			if ( dumpbeam_ ) samplesheets_ = false;
			solutionset.clear();
			read_from_disk(solutionset, cycle, acceptlower, false);
			cycle++;
			stepcount++;
			torsionrangehi += storelow_;
			torsionrangehi += storehi_;
			if ( update_pose ) update_to_stored( pose, pose_cen, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, lower_fasta_, upper_fasta_, storelow_, storehi_, is_nterm, is_cterm);
			total_residues_ = pose.total_residue();
			if ( direction_ == 1 && !is_nterm ) {
				insert_pose_ = lower_pose-((int)maxfrag-1);
			} else {
				insert_pose_ = upper_term;
			}
			if ( !update_pose ) {
				upper_fasta_ -= storehi_;
				lower_fasta_ += storelow_;
			}

			//If coordinate file is provided filter the solutions based on the coordinates.
			if ( coordfile_ != "" ) {
				coordinate_filter(solutionset, pose_cen, acceptlower, lower_fasta_, upper_fasta_, torsionrangelo, torsionrangehi );
			}

			//this line updates lower_pose only if it hasn't already been updated in the update_to_stored function.
			if ( !update_pose ) lower_pose += storelow_;

			//bool lowerstore = storelow_ > storehi_;
			Size total_lower = lower_pose - torsionrangelo+1;

			if ( auto_stop_ ) {
				LoopPartialSolution lps = solutionset[1];
				lps.apply(pose_cen, torsionrangelo, torsionrangehi);
				check_auto_stop( pose_cen, torsionrangelo, lower_pose );
				check_auto_stop( pose_cen, lower_pose+1, upper_pose );
			}

			if ( rescorebeams_ ) {
				add_user_csts( pose );
				add_user_csts( pose_cen);
				if ( !is_nterm && !is_cterm && (upper_fasta_-lower_fasta_+1 <= 10) ) {
					core::Real len = gapdist(upper_fasta_-lower_fasta_+1);
					core::scoring::constraints::BoundFuncOP myfunc(new core::scoring::constraints::BoundFunc( 0.0, len, 1.0, "x"));
					core::scoring::constraints::AtomPairConstraintOP close_cst(
						new core::scoring::constraints::AtomPairConstraint(
						core::id::AtomID(pose.residue(lower_pose).atom_index("C"),lower_pose),
						core::id::AtomID(pose.residue(lower_pose+1).atom_index("N"),lower_pose+1),
						myfunc
						));

					pose.add_constraint( close_cst );
					pose_cen.add_constraint( close_cst );
					//Make sure constraints are actually on
					Real constraint_weight;
					if ( pack_min_cycles_ != 0 ) {
						constraint_weight = sf_->get_weight(core::scoring::atom_pair_constraint);
					} else if ( cenrot_ ) {
						constraint_weight = cenrot_sf_->get_weight(core::scoring::atom_pair_constraint);
					} else {
						constraint_weight = cen_sf_->get_weight(core::scoring::atom_pair_constraint);
					}
					if ( constraint_weight  ==  0 ) {
						TRACER << "WARNING::Constraints on closure are turned off when you need them!" << std::endl;
					}
				}
				// loop is complete
				if ( upper_fasta_<lower_fasta_ ) {
					// add chainbreak vars, turn on chainbreak score
					if ( !is_nterm && !is_cterm ) {
						core::conformation::remove_upper_terminus_type_from_conformation_residue( pose.conformation(), lower_pose );
						core::conformation::remove_lower_terminus_type_from_conformation_residue( pose.conformation(), lower_pose+1 );
						core::conformation::remove_upper_terminus_type_from_conformation_residue( pose_cen.conformation(), lower_pose );
						core::conformation::remove_lower_terminus_type_from_conformation_residue( pose_cen.conformation(), lower_pose+1 );

						pose.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");
						pose_cen.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");

						core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, lower_pose   );
						core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
						core::pose::add_variant_type_to_pose_residue( pose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
						core::pose::add_variant_type_to_pose_residue( pose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );

						sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
						cen_sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
						cenrot_sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
					} else if ( is_nterm ) {
						core::conformation::add_lower_terminus_type_to_conformation_residue( pose.conformation(), lower_pose+1 );
						core::conformation::add_lower_terminus_type_to_conformation_residue( pose_cen.conformation(), lower_pose+1 );
					} else if ( is_cterm ) {
						core::conformation::add_upper_terminus_type_to_conformation_residue( pose.conformation(), lower_pose );
						core::conformation::add_upper_terminus_type_to_conformation_residue( pose_cen.conformation(), lower_pose );
					}
				}
				if ( fafilter_ ) {
					bool minimize = minimize_;
					if ( cenrotfilter_ ) {
						cenrot_ = 1;
					} else {
						pack_min_cycles_ = fafilter_pmcycles_;
						if ( minimize_ == false ) {
							pack_min_cycles_ = 999;
						}
					}
					rescoresolutionset( solutionset, pose, pose_cen, torsionrangelo, torsionrangehi );
					cenrot_ = 0;
					pack_min_cycles_ = 0;
					minimize_ = minimize;
				} else {
					rescoresolutionset( solutionset, pose, pose_cen, torsionrangelo, torsionrangehi );
				}
			}

			if ( nativegrow_ && rescorebeams_ ) {
				if ( fafilter_ ) pack_min_cycles_ = fafilter_pmcycles_;
				addnativesolution(solutionset, pose, pose_cen, loop_.start(), lower_fasta_, upper_fasta_, loop_.stop(), resstart_+1, torsionrangelo, torsionrangehi);
				if ( fafilter_ ) pack_min_cycles_ = 0;
			}
			if ( !parametercheck_ ) {
				if ( pack_min_cycles_ != 0 || fafilter_ || cenrot_ ) {
					solutionset.filter(pose, beamwidth_, total_lower, torsionrangelo, torsionrangehi, cycle);
				} else {
					solutionset.filter(pose_cen, beamwidth_, total_lower, torsionrangelo, torsionrangehi, cycle);
				}
			} else if ( pack_min_cycles_ != 0 && !cenrot_ ) {
				solutionset.parametercheck_filter(pose, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
					lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
			} else {
				solutionset.parametercheck_filter(pose_cen, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
					lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
			}
			write_to_disk( solutionset, cycle-1, storelow_, storehi_, true, acceptlower);

		}

		if ( (upper_fasta_ < lower_fasta_) && (writelps_) ) {
			if ( !update_pose ) lower_pose -= storelow_;
			upper_fasta_ += storehi_;
			lower_fasta_ -= storelow_;
			update_and_writelps(solutionset, pose, pose_cen, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, is_nterm, is_cterm, fasta_range_low, fasta_range_hi,
				pose_range_low, pose_range_hi, torsionrangelo, torsionrangehi, tgt_jump, !update_pose);
		}
		if ( readbeams_ && steps_ == 0 ) {
			exit(0);
		}
		if ( readbeams_ ) {
			solutionset.clear();
			read_from_disk(solutionset, cycle, acceptlower, false);
			cycle +=1;
			solutionsetlower = solutionset;
			readbeams_ = false;
		}
		if ( (storelow_ != 0) || (storehi_ !=0) ) {
			update_to_stored( pose, pose_cen, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, lower_fasta_, upper_fasta_, storelow_, storehi_, is_nterm, is_cterm);
			total_residues_ = pose.total_residue();
			initial_melt_right = fragmelt_;
			initial_melt_left = fragmelt_;
			upper_term += storelow_;
			torsionrangehi += storelow_;
			torsionrangehi += storehi_;
			storelow_ = 0;
			storehi_ = 0;
		}
		if ( upper_fasta_ < lower_fasta_ ) {
			TRACER << "last grow" <<std::endl;
			done = true;
		}

		Real scorelower = 99;
		Real scoreupper = 99;
		core::pose::Pose lowgrowpose = pose;
		core::pose::Pose lowgrowpose_cen = pose_cen;
		core::pose::Pose uppergrowpose = pose;
		core::pose::Pose upperpose_cen = pose_cen;
		//extend the pose
		if ( !done ) {
			if ( growlower ) { scorelower = single_grow( lowgrowpose, lowgrowpose_cen,  solutionsetlower, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, upper_term, lower_fasta_, upper_fasta_,
				torsionrangelo, torsionrangehi, true, initial_melt_left, is_cterm, is_nterm, cycle, n_to_insert_lower );
			}
			if ( growupper ) { scoreupper = single_grow( uppergrowpose, upperpose_cen, solutionset, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, upper_term, lower_fasta_, upper_fasta_,
				torsionrangelo, torsionrangehi, false, initial_melt_left, is_cterm, is_nterm,  cycle, n_to_insert_upper );
			}
			if ( solutionsetlower.size() < solutionset.size() ) {
				acceptlower = true;
			} else {
				acceptlower = false;
			}
			if ( solutionsetlower.size() == solutionset.size() ) {
				if ( scorelower <= scoreupper ) {
					acceptlower = true;
				} else {
					acceptlower = false;
				}
			}
			if ( growlower && !growupper ) {
				acceptlower = true;
			}
			if ( growupper && !growlower ) {
				acceptlower = false;
			}
			if ( ( acceptlower ) ) {
				cycle++;
				stepcount++;
				pose = lowgrowpose;
				pose_cen = lowgrowpose_cen;
				solutionset = solutionsetlower;
				initial_melt_right = fragmelt_;
				lower_pose += n_to_insert_lower;
				upper_pose += n_to_insert_lower;
				upper_term += n_to_insert_lower;
				torsionrangehi += n_to_insert_lower;
				lower_fasta_ += n_to_insert_lower;
			} else {
				cycle++;
				stepcount++;
				pose = uppergrowpose;
				pose_cen = upperpose_cen;
				initial_melt_left = fragmelt_;
				solutionsetlower = solutionset;
				upper_pose += n_to_insert_upper;
				torsionrangehi += n_to_insert_upper;
				upper_fasta_ -= n_to_insert_upper;
			}
			if ( writebeams_ ) {
				Size added_lower = lower_fasta_-loop_.start();
				Size added_upper = loop_.stop()-upper_fasta_;
				write_to_disk( solutionset, cycle-1, added_lower, added_upper, false, acceptlower);
			}
			if ( ((int)steps_ !=0) && (stepcount-1 >= steps_) ) {
				TRACER << "all steps complete" << std::endl;
				exit(0);
			}

		}
		if ( upper_fasta_ < lower_fasta_ ) {
			done = true;
		}
		//if the job is finished write a file simply named "finished" this is used by the python script to detect when loop is complete
		if ( writebeams_ && done ) {
			std::ofstream outbeam;
			outbeam.open("finished.txt", std::ofstream::app);
			outbeam << "finished" << std::endl;
			outbeam.close();
		}

	}

	upper_pose += upper_fasta_-lower_fasta_+1;
	torsionrangehi = upper_pose+minmelt_right;
	if ( is_cterm ) torsionrangehi = upper_pose;
	Real score_best = 1e10;
	core::pose::Pose best_pose;
	core::kinematics::FoldTree ft_final = pose.fold_tree();
	int fctr=0;
	restore_sc = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( pose ) );

	//get the cutpoint for transfer
	cutpoint_fasta = 0;
	if ( !is_nterm && !is_cterm ) {
		cutpoint_pose = pose.fold_tree().cutpoint_by_jump(tgt_jump);
		Size delta_cutpoint = cutpoint_pose - torsionrangelo;
		cutpoint_fasta = fasta_range_low + delta_cutpoint-1;
	}
	if ( is_cterm ) {
		cutpoint_fasta = pose_range_hi;
	}
	if ( is_nterm ) {
		cutpoint_fasta = lower_fasta_-1;
	}


	//creates a final solutionset for the comparator
	LoopPartialSolutionStore solutionsetfinal( beamwidth_, rmscutoff_, master_beam_width_, master_beam_cutoff_ );
	//setting these parameters is probably unncessary but just incase I want to use them later
	solutionsetfinal.setfilterparams( fragmelt_, rmswindow_, parallelcount_, beamscorecutoff_, dumperrors_, dumpbeam_, writebeams_, clustercheck_, fafilter_, samplesheets_, filterprevious_,
		checksymm_, asymmdump_, dumpfinalbeam_ );

	Real bestRMS=999;
	Real bestGDT=0;
	for ( Size i=1; i<=solutionset.size(); i++ ) {
		fctr++;
		LoopPartialSolution lps = solutionset[i];
		if ( !is_nterm && !is_cterm ) {
			pose.fold_tree(ft_final);
			pose_cen.fold_tree(ft_final);
		}
		lps.apply (pose_cen, torsionrangelo, torsionrangehi );
		//if(greedy_) restore_sc->apply( pose );
		if ( dumpbeam_ ) {
			pose_cen.conformation().chains_from_termini();
			core::pose::PDBInfo newpdbinfo(pose_cen, false);
			newpdbinfo.attach_to(pose_cen.conformation());
			pose_cen.dump_pdb("topscoring_final_"+utility::to_string(i)+".pdb");
		}


		// final refine cycle
		if ( !is_nterm && !is_cterm ) {
			// re-add chainbreak vars + cut
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			core::pose::add_variant_type_to_pose_residue( pose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::add_variant_type_to_pose_residue( pose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
			cen_sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
		}
		bool finalfullatomrefine;
		core::pose::Pose * finalpose;
		if ( !greedy_ && pack_min_cycles_ ==0 ) {
			finalfullatomrefine=false;
			finalpose = &pose_cen;
		} else {
			finalfullatomrefine = false;
			finalpose = &pose_cen;
		}

		refine_cycle( pose, pose_cen, torsionrangelo, torsionrangehi, finalfullatomrefine, 0, 0, 0, 0);
		core::scoring::ScoreType elec_dens = core::scoring::score_type_from_name( "elec_dens_fast" );
		core::Real residue_dens = (sf_->get_weight(core::scoring::elec_dens_fast) * finalpose->energies().residue_total_energies(torsionrangelo)[elec_dens]);
		core::Real worst_residue = residue_dens;
		core::Real total_resdens = 0;
		core::Size total_residues = 1;
		for ( core::Size ii = (unsigned)torsionrangelo+1; ii <= (unsigned)torsionrangehi; ii++ ) {
			core::Real residue_dens = ( sf_->get_weight(core::scoring::elec_dens_fast ) * (finalpose->energies().residue_total_energies(ii)[elec_dens]));
			total_resdens += residue_dens;
			total_residues++;
			if ( residue_dens > worst_residue ) {
				worst_residue = residue_dens;
			}
		}
		core::Real adjusted_dens = 0;
		for ( core::Size ii = (unsigned)torsionrangelo+1; ii <= (unsigned)torsionrangehi; ii++ ) {
			core::Real residue_dens = ( sf_->get_weight(core::scoring::elec_dens_fast ) * (finalpose->energies().residue_total_energies(ii)[elec_dens]));
			adjusted_dens -= continuous_weight_ * std::abs(residue_dens-worst_residue);
		}
		core::Real minimal_dens = total_residues * worst_residue;
		//core::Real average_res_dens = total_resdens/total_residues;
		Real score;
		if ( !cenrot_ ) {
			score = (*sf_)(pose) - total_resdens + adjusted_dens + minimal_dens;
		}
		Real cen_score;
		if ( cenrot_ ) {
			cen_score = (*cenrot_sf_)(pose);
			score = cen_score;
		} else {
			cen_score = (*cen_sf_)(pose_cen); std::ofstream outbeam;
			if ( pack_min_cycles_ == 0 ) score = cen_score;
		}


		TRACER << pose.fold_tree() << std::endl;
		if ( !is_nterm && !is_cterm ) {
			core::kinematics::FoldTree ft_closed = pose.fold_tree();
			if ( !dont_cut_jump ) ft_closed.delete_jump_and_intervening_cutpoint(tgt_jump);
			pose.fold_tree(ft_closed);

			// remove chainbreak vars + cut
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			core::pose::remove_variant_type_from_pose_residue( pose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( pose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			sf_->set_weight( core::scoring::chainbreak, 0.0 );
			cen_sf_->set_weight( core::scoring::chainbreak, 0.0 );
		}
		TRACER << "beam score = " << score << std::endl;

		//calculate RMS of loop to native. If native and loop have different numbering number will be junk
		if ( native_ ) {
			Real GDT = GDThatonative(pose, rmsrangelo_, rmsrangehi_, 0, 0, torsionrangelo, 0, 0, torsionrangehi);
			Real RMS = RMStonative(pose, rmsrangelo_, rmsrangehi_, 0, 0, torsionrangelo, 0, 0, torsionrangehi);
			if ( RMS <= bestRMS ) bestRMS = RMS;
			if ( GDT >= bestGDT ) bestGDT = GDT;
			TRACER << "The RMS is ====================================== " << RMS << " GDT is " << GDT << std::endl;
		}
		finalpose->conformation().chains_from_termini();
		core::pose::PDBInfo newpdbinfo(*finalpose, false);
		newpdbinfo.attach_to(finalpose->conformation());
		if ( dumpfinalbeam_ ) {
			finalpose->dump_pdb("finalbeam_"+utility::to_string(score)+"_"+utility::to_string( i )+".pdb" );
		}

		if ( score < score_best ) {
			score_best = score;
			best_pose = pose;
		}
		lps.store_coordinates(*finalpose, torsionrangelo, torsionrangehi, fasta_range_low-1);
		solutionsetfinal.store(lps);
	}

	//store solutions for comparator
	solutionset_ = solutionsetfinal;
	solutionset_.set_fastas(fasta_range_low-1, fasta_range_hi+1);
	solutionset_.set_poses(pose_range_low, pose_range_hi);
	solutionset_.set_cutpoint(cutpoint_fasta);

	pose = best_pose;
	TRACER << "the best score is " << score_best << std::endl;

	//calculate RMS of loop to native.
	if ( native_ ) {
		Real RMS = RMStonative(pose, rmsrangelo_, rmsrangehi_, 0, 0, torsionrangelo, 0, 0, torsionrangehi);
		Real GDT = GDThatonative(pose, rmsrangelo_, rmsrangehi_, 0, 0, torsionrangelo, 0, 0, torsionrangehi);
		TRACER << "The best scoring  RMS is ====================================== " << RMS << " GDT == " << GDT << std::endl;
		TRACER << "The best RMS is =============================================== " << bestRMS << " GDT == " << GDT << std::endl;
		if ( dumprms_ ) {
			solutionset.report_rms_and_scores(loop_.start(), loop_.stop());
		}
	}
	if ( cenrot_ ) {
		protocols::simple_moves::SwitchResidueTypeSetMoverOP tofullatom = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
		tofullatom->apply(pose);
	}

}
void
LoopGrower::addnativesolution(LoopPartialSolutionStore& solutionset, core::pose::Pose& fa_pose, core::pose::Pose& cen_pose, int natlowerstart, int natlowerstop, int natupperstart, int natupperstop,
	int lower_pose, Size torsionrangelo, Size torsionrangehi){

	// centroid copy of the pose
	protocols::simple_moves::SwitchResidueTypeSetMoverOP tocen;
	if ( cenrot_ ) {
		tocen = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid_rot")); // tocen("centroid_rot");
	} else {
		tocen = protocols::simple_moves::SwitchResidueTypeSetMoverOP( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid")); // tocen("centroid_rot");
	}
	core::pose::Pose natpose_cen = *native_;
	tocen->apply(natpose_cen);

	utility::vector1< numeric::xyzVector<core::Real> > positions;
	utility::vector1< numeric::xyzVector<core::Real> > fapositions;
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< core::id::AtomID > faids;
	Size k = lower_pose;
	for ( int i=natlowerstart; i<=natupperstop; ++i ) {
		if ( i >= natlowerstop && i <= natupperstart ) continue;
		for ( Size j=1; j<=natpose_cen.residue_type(i).natoms(); j++ ) {
			core::id::AtomID atomid = core::id::AtomID(j,k);
			ids.push_back( atomid );
			positions.push_back( natpose_cen.xyz(core::id::AtomID(j,i)));
		}
		k++;
	}
	k = lower_pose;
	for ( int i=natlowerstart; i<=natupperstop; ++i ) {
		if ( i >= natlowerstop && i <= natupperstart ) continue;
		for ( Size j=1; j<=native_->residue_type(i).natoms(); j++ ) {
			core::id::AtomID atomid = core::id::AtomID(j,k);
			faids.push_back( atomid );
			fapositions.push_back( native_->xyz(core::id::AtomID(j,i)));
		}
		k++;
	}
	cen_pose.batch_set_xyz(ids, positions);
	fa_pose.batch_set_xyz(faids, fapositions);
	if ( debug_ ) {
		cen_pose.conformation().chains_from_termini();
		core::pose::PDBInfo newpdbinfo(cen_pose, false);
		newpdbinfo.attach_to(cen_pose.conformation());
		cen_pose.dump_pdb("aftercopyfromnative.pdb");
	}
	Real prerms = GDThatonative(cen_pose, rmsrangelo_, rmsrangehi_, lower_fasta_-1, upper_fasta_+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
	TRACER << "rms before refine is " << prerms << std::endl;
	refine_cycle( fa_pose, cen_pose, torsionrangelo, torsionrangehi, false, 0, 0, 0, 0);
	Real modscore = modifiedscore(fa_pose, cen_pose, torsionrangelo, torsionrangehi);
	LoopPartialSolution nativelps= LoopPartialSolution( fa_pose, torsionrangelo, torsionrangehi, modscore );
	Real RMS;
	Real GDT;
	if ( pack_min_cycles_ !=0 || cenrot_ ) {
		GDT = GDThatonative(fa_pose, rmsrangelo_, rmsrangehi_, lower_fasta_-1, upper_fasta_+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
		RMS = RMStonative(fa_pose, rmsrangelo_, rmsrangehi_, lower_fasta_-1, upper_fasta_+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
	} else {
		GDT = GDThatonative(cen_pose, rmsrangelo_, rmsrangehi_, lower_fasta_-1, upper_fasta_+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
		RMS = RMStonative(cen_pose, rmsrangelo_, rmsrangehi_, lower_fasta_-1, upper_fasta_+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
	}
	if ( debug_ && pack_min_cycles_ != 0 ) {
		fa_pose.conformation().chains_from_termini();
		core::pose::PDBInfo newpdbinfo(fa_pose, false);
		newpdbinfo.attach_to(fa_pose.conformation());
		fa_pose.dump_pdb("faafternative.pdb");
	}
	nativelps.set_rms(RMS);
	nativelps.set_gdt(GDT);
	TRACER << " native solution rms is " << RMS << " and the score is " << modscore << std::endl;
	solutionset.store(nativelps);
}
void
LoopGrower::update_and_writelps(LoopPartialSolutionStore & solutionset, core::pose::Pose & fa_pose, core::pose::Pose & pose_cen, core::chemical::ResidueTypeCOPs & restypes_pose,
	core::chemical::ResidueTypeCOPs & restypes_pose_cen, int lower_pose, int upper_pose, bool is_nterm, bool is_cterm, Size fasta_range_low, Size fasta_range_hi, Size pose_range_low,
	Size pose_range_hi, int torsionrangelo, int torsionrangehi, int tgt_jump, bool update_pose){

	if ( update_pose ) {
		update_to_stored( fa_pose, pose_cen, restypes_pose, restypes_pose_cen, lower_pose, upper_pose, lower_fasta_, upper_fasta_,
			storelow_, storehi_, is_nterm, is_cterm);
	}
	int cutpoint_fasta = 0;
	if ( !is_nterm && !is_cterm ) {
		int cutpoint_pose = fa_pose.fold_tree().cutpoint_by_jump(tgt_jump);
		int delta_cutpoint = cutpoint_pose - torsionrangelo;
		cutpoint_fasta = fasta_range_low + delta_cutpoint-1;
	}
	if ( is_cterm ) {
		cutpoint_fasta = pose_range_hi;
	}
	if ( is_nterm ) {
		cutpoint_fasta = lower_fasta_-1;
	}
	LoopPartialSolutionStore newsolutionset;
	newsolutionset.set_fastas(fasta_range_low, fasta_range_hi);
	newsolutionset.set_poses(pose_range_low, pose_range_hi);
	newsolutionset.set_cutpoint(cutpoint_fasta);
	for ( Size i=1; i<=solutionset.size(); i++ ) {
		LoopPartialSolution lps = solutionset[i];
		if ( pack_min_cycles_ == 0 && !cenrot_ ) {
			lps.apply( pose_cen, torsionrangelo, torsionrangehi );
			lps.store_coordinates(pose_cen, torsionrangelo, torsionrangehi, fasta_range_low);
		} else {
			lps.apply( fa_pose, torsionrangelo, torsionrangehi );
			lps.store_coordinates(fa_pose, torsionrangelo, torsionrangehi, fasta_range_low);
		}
		newsolutionset.store(lps);
	}
	solutionset = newsolutionset;
	solutionset.writelpsstore(loopnumber_, parallelcount_);
	std::ofstream outbeam;
	outbeam.open("finished.txt", std::ofstream::app);
	outbeam << "finished" << std::endl;
	outbeam.close();

}

void
LoopGrower::full_atom_beam( LoopPartialSolutionStore& solutionset, core::pose::Pose & fa_pose, core::pose::Pose & cen_pose, Size lower_pos, Size upper_pos){
	Size initialpackmin = pack_min_cycles_;
	LoopPartialSolutionStore fasolutionset = solutionset;
	fasolutionset.clear();
	for ( Size i=1; i<=solutionset.size(); i++ ) {
		LoopPartialSolution lps = solutionset[i];
		lps.apply( cen_pose, lower_pos, upper_pos );
		Real censcore = lps.get_score();
		refine_cycle( fa_pose, cen_pose, lower_pos, upper_pos, true, 0, 0, 0, 0);
		Real fa_score = (*sf_)(fa_pose);
		core::scoring::ScoreType elec_dens = core::scoring::score_type_from_name( "elec_dens_fast" );
		core::Real residue_dens = (sf_->get_weight(core::scoring::elec_dens_fast) * fa_pose.energies().residue_total_energies(lower_pos)[elec_dens]);
		core::Real total_resdens = residue_dens;
		for ( core::Size ii = lower_pos+1; ii <= upper_pos; ii++ ) {
			core::Real residue_dens = ( sf_->get_weight(core::scoring::elec_dens_fast ) * (fa_pose.energies().residue_total_energies(ii)[elec_dens]));
			total_resdens += residue_dens;
		}
		LoopPartialSolution fa_solution = LoopPartialSolution( fa_pose, lower_pos, upper_pos, fa_score );
		fa_solution.set_score(censcore);
		fa_solution.set_bonus_score(total_resdens);
		fasolutionset.store(fa_solution);
	}
	pack_min_cycles_ = initialpackmin;
	solutionset = fasolutionset;
}
void
LoopGrower::fafilter( LoopPartialSolutionStore &solutionset, core::pose::Pose &fapose, core::pose::Pose &cenpose, Size total_lower, Size torsionrangelo, Size torsionrangehi,
	Size cycle, Size lower_fasta, Size upper_fasta, Size lower_pose){

	Size maxbeam = beamwidth_ * 2 ;
	bool oldfilterprevious = filterprevious_;
	solutionset.set_filterprevious(false);
	if ( parametercheck_ ) {
		solutionset.parametercheck_filter(fapose, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
			lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
	} else {
		solutionset.filter(cenpose, maxbeam, total_lower, torsionrangelo, torsionrangehi, cycle );
	}
	if ( cenrotfilter_ ) {
		cenrot_ = 1;
	} else {
		pack_min_cycles_ = fafilter_pmcycles_;
	}
	rescoresolutionset( solutionset, fapose, cenpose, torsionrangelo, torsionrangehi );
	pack_min_cycles_ = 0;
	cenrot_ = 0;
	if ( native_ ) {
		LoopPartialSolution rmslps;
		LoopPartialSolutionStore updaterms = solutionset;
		updaterms.clear();
		for ( Size i=1; i<=solutionset.size(); i++ ) {
			rmslps = solutionset[i];
			rmslps.apply(fapose, torsionrangelo, torsionrangehi);
			Real RMS = RMStonative(fapose, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
			Real GDT = GDThatonative(fapose, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
			rmslps.set_rms(RMS);
			rmslps.set_gdt(GDT);
			updaterms.store(rmslps);
		}
		solutionset = updaterms;
	}
	solutionset.set_filterprevious(oldfilterprevious);
	if ( parametercheck_ ) {
		solutionset.parametercheck_filter(fapose, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
			lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
	} else {
		solutionset.filter(fapose, beamwidth_, total_lower, torsionrangelo, torsionrangehi, cycle );
	}
}

void
LoopGrower::refine_cycle( core::pose::Pose & refinepose, core::pose::Pose & refinepose_cen, int loop_start, int loop_end, bool finalrefinement, int cycle, int beam, int fragment,
	Size is_lower) {

	//find gap
	Size cutpoint = 0;
	for ( int i=loop_start; i<=loop_end; i++ ) {
		if ( refinepose_cen.fold_tree().is_cutpoint(i) ) cutpoint = i;
	}
	if ( loop_start != 1 ) {
		runtime_assert( cutpoint != 0 );
	}
	Size lower_min = cutpoint-minmelt_-1;
	Size upper_min = cutpoint+minmelt_;
	if ( (int)lower_min < loop_start ) lower_min = loop_start;
	if ( (int)upper_min > loop_end ) upper_min = loop_end;

	// setup minimizer
	core::optimization::AtomTreeMinimizerOP minimizer;
	if ( core::pose::symmetry::is_symmetric( refinepose ) ) {
		core::optimization::AtomTreeMinimizerOP symm_min( new core::optimization::symmetry::SymAtomTreeMinimizer);
		minimizer = utility::pointer::dynamic_pointer_cast< core::optimization::AtomTreeMinimizer >(symm_min);
	} else {
		minimizer = core::optimization::AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer );
	}

	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-4, true, false, false);
	options.max_iter(200);
	core::kinematics::MoveMap mm;
	for ( Size i=lower_min; i<=upper_min; ++i ) {
		mm.set_bb ( i, true ); mm.set_chi ( i, true );
	}
	//Allow residues from sheet sampler to refine
	for ( Size i=total_residues_+1; i<=refinepose_cen.total_residue(); ++i ) {
		mm.set_bb ( i, true ); mm.set_chi ( i, true );
	}
	if ( core::pose::symmetry::is_symmetric( refinepose ) ) {
		core::pose::symmetry::make_symmetric_movemap( refinepose, mm );
	}

	// centroid min
	if ( minimize_ ) {
		(*cen_sf_)(refinepose_cen);
		minimizer->run( refinepose_cen, mm, *cen_sf_, options );
	}

	// copy torsions cen->fa
	LoopPartialSolution lps( refinepose_cen, loop_start, loop_end, 0.0);
	lps.apply( refinepose, loop_start, loop_end );

	//apply sheets
	if ( samplesheets_ ) {
		store_sheets(refinepose_cen, lps);
		lps.apply_sheets(refinepose);
	}

	// Full Atom Only
	//finalrefinement = false;
	if ( (pack_min_cycles_ !=0 || finalrefinement || cenrot_) ) {
		if ( finalrefinement && pack_min_cycles_ == 0 ) { pack_min_cycles_ = 1; }

		//  score the pose since RestrictToLoopsAndNeighbors expects it
		if ( cenrot_ ) {
			(*cenrot_sf_)(refinepose);
		} else {
			(*sf_)(refinepose);
		}

		// setup task
		core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
		protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP repackaround( new
			protocols::toolbox::task_operations::RestrictToLoopsAndNeighbors);
		protocols::loops::LoopsOP pack_core( new protocols::loops::Loops );

		pack_core->add_loop( protocols::loops::Loop( lower_min, upper_min ) );

		protocols::simple_moves::PackRotamersMoverOP pack_mover;
		if ( core::pose::symmetry::is_symmetric( refinepose ) ) {
			protocols::simple_moves::PackRotamersMoverOP symmpack( new protocols::simple_moves::symmetry::SymPackRotamersMover );
			pack_mover = utility::pointer::dynamic_pointer_cast< protocols::simple_moves::PackRotamersMover >(symmpack);
		} else {
			pack_mover = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover );
		}
		repackaround->set_loops(pack_core);
		main_task_factory->push_back( repackaround );
		main_task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::IncludeCurrent) );
		main_task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		core::pack::task::PackerTaskOP ptask( main_task_factory->create_task_and_apply_taskoperations( refinepose ) );
		if ( core::pose::symmetry::is_symmetric(refinepose) ) {
			core::pack::make_new_symmetric_PackerTask_by_requested_method(refinepose, ptask); // NK 110621
		}
		if ( pack_min_cycles_ == 999 ) {
			// prevent loop from repacking
			for ( int i=loop_start; i<=loop_end; i++ ) {
				core::pack::task::ResidueLevelTask & restask( ptask->nonconst_residue_task( i ) );
				restask.prevent_repacking();
			}
		}
		pack_mover->task_factory( main_task_factory );

		for ( int i=1; i<=(int)refinepose.total_residue(); ++i ) {
			if ( !ptask->pack_residue(i) ) continue;
			mm.set_chi ( i, true );
			if ( cenrot_ ) {
				core::conformation::Residue const &res_i = refinepose.residue(i);
				if ( res_i.aa()!=chemical::aa_gly && res_i.aa()!=chemical::aa_ala && res_i.type().has( "CEN") ) {
					mm.set( core::id::DOF_ID( core::id::AtomID( res_i.atom_index("CEN"), i ), core::id::D ), true );
					mm.set( core::id::DOF_ID( core::id::AtomID( res_i.atom_index("CEN"), i ), core::id::THETA ), true );
				}

			}
		}
		//Turn off bb min of full atom if requested
		if ( !famin_ ) {
			for ( int i=loop_start; i<=loop_end; ++i ) {
				mm.set_bb ( i, false );
			}
		}

		// run repack and min
		if ( cenrot_ ) {
			pack_mover->score_function(cenrot_sf_);
			pack_mover->apply( refinepose );
			if ( minimize_ ) minimizer->run( refinepose, mm, *cenrot_sf_, options );
		} else {
			core::Real fa_rep_base = sf_->get_weight( core::scoring::fa_rep );
			for ( int i = 1; i<=(int)pack_min_cycles_; ++i ) {
				if ( pack_min_cycles_>1 ) {
					sf_->set_weight( core::scoring::fa_rep, fa_rep_base);
				}
				options.max_iter(200);
				core::pack::pack_rotamers( refinepose, *sf_, ptask);
				if ( pack_min_cycles_ == 999 ) break;
				if ( minimize_ ) minimizer->run( refinepose, mm, *sf_, options );
			}
		}
		if ( debug_ ) {
			refinepose.conformation().chains_from_termini();
			core::pose::PDBInfo newpdbinfo(refinepose, false);
			newpdbinfo.attach_to(refinepose.conformation());
			refinepose.dump_pdb("faafterrefine_"+utility::to_string( cycle )+"."+utility::to_string(parallelcount_)+"."+utility::to_string( beam )+"."+utility::to_string( fragment )+"_"+utility::to_string( is_lower )+".pdb" );
		}
	}
	// copy torsions fa->cen
	LoopPartialSolution falps( refinepose, loop_start, loop_end, 0.0);
	falps.apply( refinepose_cen, loop_start, loop_end );
	if ( debug_ ) {
		refinepose_cen.conformation().chains_from_termini();
		core::pose::PDBInfo newpdbinfo(refinepose_cen, false);
		newpdbinfo.attach_to(refinepose_cen.conformation());
		refinepose_cen.dump_pdb("cenafterrefine_"+utility::to_string( cycle )+"."+utility::to_string(parallelcount_)+"."+utility::to_string( beam )+"."+utility::to_string( fragment )+"_"+utility::to_string( is_lower )+".pdb" );
	}
}


Real
LoopGrower::nton3_hbond_score(core::pose::Pose& pose){

	Real hbond_energies = 0;
	//populate hbond set:
	core::scoring::hbonds::HBondSet set1;
	if ( pose.conformation().residue_typeset_mode( true ) == core::chemical::CENTROID_ROT_t ) {
		(*cenrot_sf_)(pose);
	} else {
		(*cen_sf_)(pose);
	}
	pose.update_residue_neighbors();
	set1.setup_for_residue_pair_energies( pose, false, false );


	//query hbond set, get energies:
	for ( Size i = 1; i<= set1.nhbonds(); i++ ) {

		core::scoring::hbonds::HBond bond = set1.hbond( i );
		//need to access donor and acc Residues as well as ints
		int accResNum = bond.acc_res();
		int donResNum = bond.don_res();
		//get acc and donor residues from sequence numbers
		core::conformation::Residue accRes = pose.residue( accResNum );
		core::conformation::Residue donRes = pose.residue( donResNum );
		if ( std::abs( accResNum - donResNum) != 3 ) continue;
		Real energy = bond.energy();
		Real weight = bond.weight();

		hbond_energies +=  weight*energy;

	}
	return hbond_energies;


}

Real
LoopGrower::get_resrange_hbond_energy(core::pose::Pose& pose, Size lower, Size upper){

	Real hbond_energy =0;
	//populate hbond set:
	core::scoring::hbonds::HBondSet set1;
	if ( pose.conformation().residue_typeset_mode(true) ==  core::chemical::CENTROID_ROT_t ) {
		(*cenrot_sf_)(pose);
	} else if ( pose.conformation().is_centroid() ) {
		(*cen_sf_)(pose);
	} else {
		(*sf_)(pose);
	}
	pose.update_residue_neighbors();
	set1.setup_for_residue_pair_energies( pose, false, false );

	//query hbond set, get energies:
	for ( Size i = 1; i<= set1.nhbonds(); i++ ) {
		core::scoring::hbonds::HBond bond = set1.hbond( i );
		//need to access donor and acc Residues as well as ints
		int accResNum = bond.acc_res();
		int donResNum = bond.don_res();
		//get acc and donor residues from sequence numbers
		core::conformation::Residue accRes = pose.residue( accResNum );
		core::conformation::Residue donRes = pose.residue( donResNum );
		if ( (accResNum >= (int)lower && accResNum <= (int)upper) || (donResNum >= (int)lower && donResNum <= (int)upper) ) {
			Real energy = bond.energy();
			Real weight = bond.weight();
			hbond_energy +=  weight*energy;
		} else {
			continue;
		}
	}
	return hbond_energy;

}
void
LoopPartialSolution::add_sheets(core::pose::Pose & pose, Size takeoffres, Size totalres){

	std::string chemicaltype;
	if ( pose.is_fullatom() ) {
		chemicaltype = core::chemical::FA_STANDARD;
	} else if ( pose.is_centroid() ) {
		chemicaltype = core::chemical::CENTROID;
	} else {
		chemicaltype = core::chemical::CENTROID_ROT;
	}
	core::chemical::ResidueType const &ala_type = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemicaltype )->name_map("ALA");
	core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(ala_type) );
	pose.conformation().append_residue_by_jump(*newres, takeoffres );

	Size newlower = totalres/2;
	Size newupper = totalres - newlower;

	core::Size midres = pose.total_residue();
	Size count = 0;
	for ( Size i=1; i<=newlower; i++ ) {
		pose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, pose.total_residue()-count, true);
		midres++;
		count++;
	}
	for ( Size i=1; i<newupper; i++ ) {
		pose.conformation().safely_append_polymer_residue_after_seqpos(*newres, pose.total_residue(), true);
	}

}
void
LoopPartialSolution::write_beam_lps(std::ofstream &lpsfile){
	Size n = ids_.size();
	lpsfile << score_ << std::endl;
	lpsfile << n << std::endl;
	for ( Size i=1; i<=ids_.size(); i++ ) {
		Size atomno = ids_[i].atomno();
		Size rsd = ids_[i].rsd();
		lpsfile << atomno << " " << rsd << " ";
	}
	lpsfile << std::endl;
	for ( Size i=1; i<=positions_.size(); i++ ) {
		lpsfile << positions_[i];
	}
	lpsfile << std::endl;
}
void
LoopPartialSolution::write_beam(std::ofstream &outbeam){
	for ( Size i=1; i<=residues_.size(); i++ ) {
		ResTorsions storeres = residues_[i];
		outbeam << storeres.phi_ << " ";
		outbeam << storeres.psi_ << " ";
		outbeam << storeres.omega_ << " ";
		outbeam << storeres.nchi_ << " ";
		if ( storeres.nchi_ > 0 ) outbeam << storeres.chi1_ << " ";
		if ( storeres.nchi_ > 1 ) outbeam << storeres.chi2_ << " " ;
		if ( storeres.nchi_ > 2 ) outbeam << storeres.chi3_ << " ";
		if ( storeres.nchi_ > 3 ) outbeam << storeres.chi4_ << " ";
	}
	outbeam << std::endl;
	for ( Size j=1; j<=backbones_.size(); j++ ) {
		outbeam << backbones_[j].size() << " ";
		for ( Size k=1; k<=backbones_[j].size(); k++ ) {
			//outbeam << backbones_[j];
			outbeam << backbones_[j][k][0] << " " << backbones_[j][k][1] << " " << backbones_[j][k][2] << " ";
		}
	}
	outbeam << std::endl << sheets_.size() << " " << bonus_score_ << std::endl;
	for ( Size i=1; i<=sheets_.size(); i++ ) {
		outbeam << sheets_[i].baseres_ << " " << sheets_[i].jumpid_ << std::endl;
		for ( Size ii=1; ii<=3; ii++ ) {
			for ( Size j=1; j<=3; j++ ) {
				outbeam << sheets_[i].rotation_(ii,j) << " ";
			}
		}
		outbeam << std::endl;
		for ( Size j=0; j<3; j++ ) {
			outbeam << sheets_[i].translation_[j] << " " ;
		}
		outbeam << std::endl;
		outbeam << sheets_[i].residues_.size() << " ";
		for ( Size j=1; j<=sheets_[i].residues_.size(); j++ ) {
			ResTorsions storeres = sheets_[i].residues_[j];
			outbeam << storeres.phi_ << " ";
			outbeam << storeres.psi_ << " ";
			outbeam << storeres.omega_ << " ";
			outbeam << storeres.nchi_ << " ";
			if ( storeres.nchi_ > 0 ) outbeam << storeres.chi1_ << " ";
			if ( storeres.nchi_ > 1 ) outbeam << storeres.chi2_ << " " ;
			if ( storeres.nchi_ > 2 ) outbeam << storeres.chi3_ << " ";
			if ( storeres.nchi_ > 3 ) outbeam << storeres.chi4_ << " ";
		}
		outbeam << std::endl;
	}
	outbeam << score_ << " " << RMS_ << " " << GDT_ << std::endl;
}

void
LoopPartialSolution::apply_sheets( core::pose::Pose & pose ){
	for ( Size i=1; i<= sheets_.size(); i++ ) {
		Size numjumps = pose.num_jump();
		Size totalres = pose.total_residue();
		Size newres = sheets_[i].residues_.size();
		Size takeoffres = sheets_[i].baseres_;
		add_sheets(pose, takeoffres, newres);
		core::kinematics::Jump newjump = pose.jump(numjumps+1);
		newjump.set_rotation( sheets_[i].rotation_ );
		newjump.set_translation( sheets_[i].translation_ );
		pose.set_jump(numjumps+1, newjump);
		for ( Size j=1; j<=sheets_[i].residues_.size(); j++ ) {
			ResTorsions const &res_j = sheets_[i].residues_[j];
			pose.set_phi( totalres+j, res_j.phi_ );
			pose.set_psi( totalres+j, res_j.psi_ );
			pose.set_omega( totalres+j, res_j.omega_ );
			if ( !pose.is_centroid() && pose.conformation().residue_typeset_mode(true) != core::chemical::CENTROID_ROT_t ) {
				if ( res_j.nchi_ > 0 ) pose.set_chi( 1, totalres+j, res_j.chi1_ );
				if ( res_j.nchi_ > 1 ) pose.set_chi( 2, totalres+j, res_j.chi2_ );
				if ( res_j.nchi_ > 2 ) pose.set_chi( 3, totalres+j, res_j.chi3_ );
				if ( res_j.nchi_ > 3 ) pose.set_chi( 4, totalres+j, res_j.chi4_ );
			}

		}
	}
}

Real LoopPartialSolution::max_calpha_distance( const LoopPartialSolution & newlps, Size fragmelt, Size total_lower, bool takeoffonly, bool full_loop, bool is_lower ) const {
	Real largest_distance = 0.0;
	Real new_distance = 0.0;
	utility::vector1< core::Vector > newcalphas = newlps.get_calphas();
	utility::vector1< utility::vector1< core::Vector > > newbackbone = newlps.get_backbone();
	runtime_assert( calphas_.size() == newcalphas.size() );
	if ( takeoffonly ) {
		Size i;
		if ( is_lower ) {
			i = (total_lower-fragmelt-1) * 3;
		} else {
			i = (total_lower+fragmelt) * 3 + 1;
		}
		for ( Size j=1; j<=5; j++ ) {
			i++;
		}

	} else {
		for ( Size i=1; i<=newcalphas.size(); i++ ) {
			if ( !full_loop ) {
				if ( (i <= total_lower && i > total_lower-fragmelt+1) || (i > total_lower && i < total_lower+fragmelt) ) {
					continue;
				}
			}
			if ( takeoffonly ) {
				if ( is_lower ) {
					if ( i != total_lower-fragmelt ) {
						continue;
					}
				} else {
					if ( i != total_lower+fragmelt+1 ) {
						continue;
					}
				}
			}
			new_distance = (calphas_[i]-newcalphas[i]).length();
			if ( new_distance > largest_distance ) {
				largest_distance = new_distance;
			}
		}
	}
	return largest_distance;
}

Real LoopPartialSolution::partialrms(const LoopPartialSolution & newlps, int fragmelt, int total_lower, bool takeoffonly, int rmswindow ) const {
	Real worstresrms=-1;
	utility::vector1< core::Vector > newcalphas = newlps.get_calphas();
	utility::vector1< utility::vector1< core::Vector > > newbackbones = newlps.get_backbone();
	runtime_assert( backbones_.size() == newbackbones.size() );
	//if rmswindow set to 0 loop over all the residues
	if ( rmswindow == 0 ) rmswindow = residues_.size();
	int rmswindowupper = rmswindow;
	int rmswindowlower = rmswindow;
	int lowercount = total_lower-fragmelt+2;
	int uppercount = residues_.size() - (total_lower+fragmelt);
	Size upperstart = total_lower+fragmelt;
	if ( takeoffonly ) {
		lowercount -= 2;
		uppercount -= 2;
		upperstart += 2;
	}
	if ( rmswindowupper > uppercount ) {
		rmswindowupper = uppercount;
		if ( rmswindowupper < 1 ) rmswindowupper=1;
	}
	if ( rmswindowlower > lowercount ) {
		rmswindowlower = lowercount;
		if ( rmswindowlower < 1 ) rmswindowlower=1;
	}
	//TRACER << "total lower " << total_lower << std::endl;
	for ( int i=1; i<=lowercount-rmswindowlower+1; i++ ) {
		Real err = 0;
		Size atmcount=0;
		for ( int ii=1; ii<=rmswindowlower; ii++ ) {
			int resi = i+ii-1;
			if ( resi > (int)residues_.size() ) break;
			Size j=1;
			while ( j<=backbones_[resi].size() ) {
				Real new_distance = (backbones_[resi][j]-newbackbones[resi][j]).length_squared();
				err += new_distance;
				atmcount++;
				j++;
			}
			//TRACER << "using resi " << resi << " for lower " << std::endl;
		}
		Real windowRMS = sqrt(err/atmcount);
		//TRACER << " window rms " << windowRMS << std::endl;
		if ( windowRMS>worstresrms ) {
			worstresrms=windowRMS;
		}
	}
	for ( int i=upperstart; i<=(int)residues_.size()-rmswindowupper; i++ ) {
		Real err = 0;
		Size atmcount=0;
		for ( int ii=1; ii<=rmswindowupper; ii++ ) {
			int resi = i+ii-1;
			if ( resi > (int)residues_.size() ) break;
			Size j=1;
			while ( j<=backbones_[resi].size() ) {
				Real new_distance = (backbones_[resi][j]-newbackbones[resi][j]).length_squared();
				err += new_distance;
				atmcount++;
				j++;
			}
			//TRACER << "using resi " << resi << " for upper " << std::endl;
		}
		Real windowRMS = sqrt(err/atmcount);
		if ( windowRMS>worstresrms ) {
			worstresrms=windowRMS;
		}
	}
	return worstresrms;
}

void LoopPartialSolution::set_centerofmass(){
	Real totalx = 0;
	Real totaly = 0;
	Real totalz = 0;
	Size totalatoms = 0;
	for ( Size i=1; i<=backbones_.size(); i++ ) {
		for ( Size j=1; j<=backbones_[i].size(); j++ ) {
			totalx += backbones_[i][j][0];
			totaly += backbones_[i][j][1];
			totalz += backbones_[i][j][2];
			totalatoms++;
		}
	}
	Real centerx = totalx/totalatoms;
	Real centery = totaly/totalatoms;
	Real centerz = totalz/totalatoms;
	centerofmass_.assign(centerx,centery,centerz);

}
void LoopPartialSolutionStore::push( core::pose::Pose& pose, LoopPartialSolution const &new_entry, Size fragmelt, Real beamscorecutoff, Size total_lower, bool lower, bool dump_errors,
	Real nativeRMS, Size lower_res, Size upper_res){

	Real pre_best_RMS = nativeRMS;
	Real post_best_RMS = 999;
	if ( dump_errors ) { old_solutions_ = solutions_; }
	Real newscore = new_entry.get_score();
	Real worstscore = -9999;
	Real bestscore = 9999;
	Size worstindex = 1;
	bool getnewscores = false;
	//check the stored beams for best and worst score
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		Real oldscore = solutions_[i].get_score();
		Real oldRMS = solutions_[i].get_rms();
		if ( oldRMS <= pre_best_RMS ) {
			pre_best_RMS = oldRMS;
		}
		if ( oldscore > worstscore || worstscore == -9999 ) {
			worstscore = oldscore;
			worstindex = i;
		}
		if ( oldscore < bestscore || bestscore == 9999 ) {
			bestscore = oldscore;
		}
	}
	if ( newscore < bestscore && bestscore != 9999 ) bestscore = newscore;

	//removes any beams above the cutoff
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		Real oldscore = solutions_[i].get_score();
		if ( oldscore > bestscore+beamscorecutoff && bestscore != 9999 ) {
			solutions_.erase(solutions_.begin() +i-1 );
			getnewscores = true;
		}
	}


	bool considered = false, checked = false;
	if ( newscore > worstscore ) {
		considered = true;
	}
	if ( newscore > bestscore+beamscorecutoff && bestscore != 9999 ) {
		checked = true;
		considered = true;
	} else {
		bool delete_close_beams = true;
		utility::vector1<Size> worsesolutions;
		for ( Size i=1; i<=solutions_.size(); ++i ) {
			Real RMS = new_entry.max_calpha_distance( solutions_[i], fragmelt, total_lower, false, false, lower );
			Real oldscore = solutions_[i].get_score();
			if ( RMS < rmscutoff_ ) {
				if ( newscore > oldscore ) {
					considered = true;
					checked = true;
					delete_close_beams = false;
				} else {
					worsesolutions.push_back(i);
					considered = true;
					checked = true;
				}
			}
		}
		if ( delete_close_beams ) {
			for ( Size i=1; i<=worsesolutions.size(); ++i ) {
				Size j = worsesolutions[i];
				solutions_.erase(solutions_.begin() +j-1 );
				if ( i==worsesolutions.size() ) {
					solutions_.push_back(new_entry);
					Size newindex = solutions_.size();
					solutions_[newindex].set_rms(nativeRMS);
				}
			}
		}
	}
	core::Size cluster_count;
	core::Real worst_cluster_score = -99999;
	core::Size worst_cluster_index;
	bool finished_master_beam = false;
	while ( !checked && !finished_master_beam ) {
		worst_cluster_score = -99999;
		cluster_count = 0;
		for ( Size i=1; i<=solutions_.size(); ++i ) {
			Real max_distance = new_entry.max_calpha_distance( solutions_[i], fragmelt, total_lower, false, false, lower );
			Real oldscore = solutions_[i].get_score();
			if ( max_distance < master_beam_cutoff_ ) {
				cluster_count++;
				if ( oldscore >= worst_cluster_score ) {
					worst_cluster_score = oldscore;
					worst_cluster_index = i;
				}
			}
		}
		if ( cluster_count >= master_beam_width_ ) {
			if ( newscore < worst_cluster_score ) {
				solutions_.erase(solutions_.begin() +worst_cluster_index-1 );
				getnewscores = true;
				if ( cluster_count == master_beam_width_ ) {
					solutions_.push_back(new_entry);
					Size newindex = solutions_.size();
					solutions_[newindex].set_rms(nativeRMS);
					considered = true;
					checked = true;
					finished_master_beam = true;
				}
			} else {
				considered = true;
				checked = true;
				finished_master_beam = true;
			}
		} else {
			finished_master_beam = true;
		}
	}
	//If beams were removed rechecks for best and worst score
	if ( getnewscores ) {
		worstscore = -9999;
		bestscore = 9999;
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			Real oldscore = solutions_[i].get_score();
			if ( oldscore > worstscore || worstscore == -9999 ) {
				worstscore = oldscore;
				worstindex = i;
			}
			if ( oldscore < bestscore || bestscore == 9999 ) {
				bestscore = oldscore;
			}
		}
	}
	if ( solutions_.size() < maxelts_ && !checked ) {
		solutions_.push_back( new_entry );
		Size newindex = solutions_.size();
		solutions_[newindex].set_rms(nativeRMS);
		considered = true;
	} else if ( !considered ) {
		solutions_[worstindex] = new_entry;
		solutions_[worstindex].set_rms(nativeRMS);
	}
	//If we've moved too far from native dump structure (debugging purposes)
	if ( dump_errors ) {
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			Real storedRMS = solutions_[i].get_rms();
			if ( storedRMS <= post_best_RMS ) {
				post_best_RMS = storedRMS;
			}
		}
		TRACER << "pre best is " << pre_best_RMS << " post best is " << post_best_RMS << std::endl;
		if ( pre_best_RMS <= 2.0 && post_best_RMS >= 2.0 ) {
			pose.conformation().chains_from_termini();
			core::pose::PDBInfo newpdbinfo(pose, false);
			newpdbinfo.attach_to(pose.conformation());
			for ( Size j=1; j<=old_solutions_.size(); j++ ) {
				Real oldRMS;
				oldRMS = old_solutions_[j].get_rms();
				old_solutions_[j].apply(pose, (int)lower_res, (int)upper_res );
				pose.dump_pdb("after_error_"+utility::to_string(j)+"_"+utility::to_string(newscore)+"_"+utility::to_string(oldRMS)+".pdb");
			}
			LoopPartialSolution new_entry_dump = new_entry;
			new_entry_dump.apply(pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("after_error_new_"+utility::to_string(newscore)+"_"+utility::to_string(nativeRMS)+".pdb");
			TRACER << " Near natives reject dumping beam " << std::endl;
		}
	}
}
void
LoopPartialSolutionStore::report_rms_and_scores(Size rangelo, Size rangehi){
	Real bestrms = solutions_[1].get_rms();
	Real bestscore = solutions_[1].get_score();
	Real bestgdt = solutions_[1].get_gdt();
	std::ofstream outputfile;
	outputfile.open("rms_and_scores.txt", std::ofstream::app);
	for ( Size i=1; i<= solutions_.size(); i++ ) {
		Real beamrms = solutions_[i].get_rms();
		Real beamgdt = solutions_[i].get_gdt();
		Real beamscore = solutions_[i].get_score();
		outputfile << "The RMS and GDT for loop " << rangelo << "-" << rangehi << " with score " << beamscore << " is " << beamrms << " and " << beamgdt << std::endl;
		if ( beamrms < bestrms ) {
			bestrms = beamrms;
			bestscore = beamscore;
			bestgdt = beamgdt;
		}
	}
	outputfile << "The Best RMS with GDT for loop " << rangelo << "-" << rangehi << " with score " << bestscore << " is " << bestrms << " with " << bestgdt << std::endl;
	outputfile.close();

}

bool
LoopPartialSolutionStore::filterprevious( LoopPartialSolution lps, Size total_lower, Size fragmelt, Real rmscutoff, Size masterbeamwidth){

	Size masterbeamcount = 0;
	Real score = lps.get_score();
	for ( Size i=1; i<=filteronly_solutions_.size(); i++ ) {
		LoopPartialSolution oldlps = filteronly_solutions_[i];
		Real oldscore = oldlps.get_score();
		Real RMS = lps.partialrms(oldlps, fragmelt, total_lower, false, 1);
		Real mbRMS = lps.partialrms(oldlps, fragmelt, total_lower, false, rmswindow_);
		if ( RMS < rmscutoff ) return false;
		if ( mbRMS < master_beam_cutoff_ && score > oldscore ) masterbeamcount++;
		if ( masterbeamcount > masterbeamwidth ) return false;
	}
	return true;

}
void
LoopPartialSolutionStore::diversityfilter(Size maxbeamsize, Size total_lower ){

	LoopPartialSolutionStore filteredsolutions;
	sort();
	Real bestscore = solutions_[1].get_score();
	std::sort(filteronly_solutions_.begin(), filteronly_solutions_.end());
	Real filteredbestscore = 999999;
	Size includeoldbeams = (maxbeamsize/10);
	if ( filterprevious_ ) {
		filteredbestscore = filteronly_solutions_[1].get_score();
		if ( bestscore > filteredbestscore ) bestscore = filteredbestscore;
		maxbeamsize-=includeoldbeams;
	}

	bool done = false;
	Size masterbeamwidth = 1;
	bool reported = false;
	Real bestRMS = 999;
	Size bestrmsrank = 0;

	//remove bad beams and those in previous set before starting normal filter.
	LoopPartialSolutionStore possiblesolutions;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		LoopPartialSolution lps = solutions_[i];
		Real beamscore = lps.get_score();
		Real absolutebestscore = std::abs(bestscore);
		Real standardizedscore = 0;
		Real denominator = absolutebestscore + std::abs(bestscore-beamscore);
		if ( denominator != 0 ) {
			standardizedscore = absolutebestscore/denominator;
		}
		if ( standardizedscore < beamscorecutoff_ ) {
			if ( !reported ) {
				reported = true;
			}
			break;
		}
		if ( filterprevious_ ) {
			bool passedprevious = filterprevious(lps, total_lower, fragmelt_, rmscutoff_, master_beam_width_);
			if ( !passedprevious ) {
				continue;
			}
		}
		possiblesolutions.store(lps);
	}
	solutions_ = possiblesolutions.get_solutions();

	if ( solutions_.size() < 1 ) done = true;
	while ( !done ) {
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			LoopPartialSolution lps = solutions_[i];
			bool addlps = true;
			Size masterbeamcount = 0;
			if ( filteredsolutions.size() >= maxbeamsize ) {
				done = true;
				TRACER << "hit max beam " << std::endl;
				break;
			}
			for ( Size j=1; j<=filteredsolutions.size(); j++ ) {
				LoopPartialSolution filteredlps = filteredsolutions[j];
				Real masterbeamrms = lps.partialrms(filteredlps, fragmelt_, total_lower, false, rmswindow_ );
				Real localbeamrms = lps.partialrms(filteredlps, fragmelt_, total_lower, false, 1 );

				if ( masterbeamrms < master_beam_cutoff_ ) {
					masterbeamcount++;
				}
				if ( masterbeamcount >= masterbeamwidth ) {
					addlps=false;
					break;
				}
				if ( localbeamrms <= rmscutoff_ ) {
					addlps = false;
					break;
				}
			}
			if ( addlps ) {
				filteredsolutions.store(lps);
				if ( dump_errors_ ) {
					Real RMS = lps.get_rms();
					if ( RMS <= 1.5 ) {
						rank_ = filteredsolutions.size();
					}
					if ( RMS < bestRMS ) {
						bestRMS = RMS;
						bestrmsrank = filteredsolutions.size();
					}
				}
			}
			if ( masterbeamwidth == master_beam_width_ ) {
				done = true;
			}
		}
		masterbeamwidth++;
	}

	if ( bestRMS > 1.5 ) {
		rank_ = bestrmsrank;
	}

	maxbeamsize+=includeoldbeams;
	//there is probably an unnecessary amount of code duplication here that should be removed but in the meantime...
	if ( filterprevious_ && filteredsolutions.size() < maxbeamsize ) {
		std::sort(filteronly_solutions_.begin(), filteronly_solutions_.end());
		Size filtermasterbeamwidth = 1;
		done = false;
		while ( !done ) {
			for ( Size i=1; i<=filteronly_solutions_.size(); i++ ) {
				bool addlps = true;
				reported = false;
				LoopPartialSolution lps = filteronly_solutions_[i];
				Size masterbeamcount = 0;
				if ( filteredsolutions.size() >= maxbeamsize ) {
					done = true;
					break;
				}
				Real beamscore = lps.get_score();
				Real absolutebestscore = std::abs(bestscore);
				Real standardizedscore = absolutebestscore/(absolutebestscore + std::abs(bestscore-beamscore));
				//TRACER << "absolutebestscore is " << absolutebestscore << std::endl;
				if ( standardizedscore < beamscorecutoff_ ) {
					if ( !reported ) {
						reported = true;
					}
					break;
				}
				for ( Size j=1; j<=filteredsolutions.size(); j++ ) {
					LoopPartialSolution storedsolution = filteredsolutions[j];
					Real masterbeamrms = lps.partialrms(storedsolution, fragmelt_, total_lower, false, rmswindow_ );
					Real localbeamrms = lps.partialrms(storedsolution, fragmelt_, total_lower, false, 1 );
					if ( masterbeamrms < master_beam_cutoff_ ) {
						masterbeamcount++;
					}
					//TRACER << " looping over filteredsolutions " << std::endl;
					if ( masterbeamcount >= filtermasterbeamwidth ) {
						addlps=false;
						break;
					}
					if ( localbeamrms <= rmscutoff_ ) {
						addlps = false;
						break;
					}
				}
				if ( addlps ) {
					//TRACER << "adding filteronly lps " << std::endl;
					filteredsolutions.store(lps);
				}
			}
			filtermasterbeamwidth++;
			if ( filtermasterbeamwidth == master_beam_width_ ) done = true;
		}
	}
	solutions_ = filteredsolutions.get_solutions();


}
void
LoopPartialSolutionStore::filter( core::pose::Pose& pose, Size maxbeam, Size total_lower, Size lower_res, Size upper_res, Size cycle ){

	sort();

	Real pre_best_RMS = solutions_[1].get_rms();
	Real pre_best_GDT = solutions_[1].get_gdt();
	LoopPartialSolution nearnativelps = solutions_[1];
	if ( dump_errors_ ) {
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			Real oldRMS = solutions_[i].get_rms();
			Real oldGDT = solutions_[i].get_gdt();

			if ( oldGDT >= pre_best_GDT ) {
				pre_best_GDT = oldGDT;
			}
			if ( oldRMS <= pre_best_RMS ) {
				pre_best_RMS = oldRMS;
				nearnativelps = solutions_[i];
			}
		}
	}

	diversityfilter( maxbeam, total_lower );


	// The rest of this function is just debugging tools
	if ( clustercheck_ && dump_errors_ ) {
		cluster_check(nearnativelps, pose, total_lower, lower_res, upper_res);
	}

	Real post_best_RMS = solutions_[1].get_rms();
	Real post_best_GDT = solutions_[1].get_gdt();
	if ( dump_beam_ ) {
		core::pose::Pose scpose = pose;
		for ( Size j=1; j<=solutions_.size(); j++ ) {
			Real oldRMS = solutions_[j].get_rms();
			std::string oldid = solutions_[j].get_id();
			Real newscore = solutions_[j].get_score();
			solutions_[j].apply(pose, (int)lower_res, (int)upper_res );
			if ( samplesheets_ ) solutions_[j].apply_sheets(pose);
			if ( checksymm_ ) {
				core::pose::Pose posecopy = pose;
				transform_to_closest_symmunit(pose,posecopy,total_lower);
			}
			pose.conformation().chains_from_termini();
			core::pose::PDBInfo newpdbinfo(pose, false);
			newpdbinfo.attach_to(pose.conformation());
			std::string name = "after_filter_"+utility::to_string(cycle)+"_"+oldid+"_"+utility::to_string(oldRMS)+"_"+utility::to_string(newscore)+".pdb";
			if ( dumpfinalbeam_ ) name = "final_"+utility::to_string(j)+".pdb";
			if ( asymmdump_ && core::pose::symmetry::is_symmetric( pose ) ) {
				core::pose::Pose asymm_pose;
				core::pose::symmetry::extract_asymmetric_unit(pose, asymm_pose, false);
				asymm_pose.dump_pdb(name);
			} else {
				pose.dump_pdb(name);
			}
			pose = scpose;
		}
		if ( dump_errors_ ) {
			std::string nearnativeid = nearnativelps.get_id();
			nearnativelps.apply(pose, (int)lower_res, (int)upper_res );
			nearnativelps.apply_sheets(pose);
			pose.conformation().chains_from_termini();
			core::pose::PDBInfo newpdbinfo(pose, false);
			newpdbinfo.attach_to(pose.conformation());
			if ( asymmdump_ && core::pose::symmetry::is_symmetric( pose ) ) {
				core::pose::Pose asymm_pose;
				core::pose::symmetry::extract_asymmetric_unit(pose, asymm_pose, false);
				asymm_pose.dump_pdb("after_filter_nn_"+utility::to_string(cycle)+"_"+nearnativeid+"_"+utility::to_string(pre_best_RMS)+".pdb");
			} else {
				pose.dump_pdb("after_filter_nn_"+utility::to_string(cycle)+"_"+nearnativeid+"_"+utility::to_string(pre_best_RMS)+".pdb");
			}
			pose = scpose;
		}
	}

	if ( dump_errors_ ) {
		sort();
		Size rank = 0;
		Size ranknat = 0;
		Real scorenearnat = 0;
		Real natrms = 0;
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			Real oldRMS = solutions_[i].get_rms();
			if ( oldRMS <= 1.5 && ranknat ==0 ) {
				ranknat = i;
				scorenearnat = solutions_[i].get_score();
				natrms = oldRMS;
			}
			if ( oldRMS <= post_best_RMS ) {
				post_best_RMS = oldRMS;
				rank = i;
			}
			Real oldGDT = solutions_[i].get_gdt();
			if ( oldGDT >= post_best_GDT ) post_best_GDT = oldGDT;
		}
		//if there is no near native solution report the rank of the best option
		if ( ranknat == 0 ) {
			ranknat = rank;
		}

		std::string nearnativeid = nearnativelps.get_id();
		TRACER << "best score is " << solutions_[1].get_score();
		TRACER << " pre best RMS, GDT, and score " << nearnativeid << " " << pre_best_RMS << " " << post_best_GDT << " " << nearnativelps.get_score() << " post_best_RMS " << post_best_RMS <<
			" post best_GDT " << post_best_GDT << " post_best ranked " << rank << std::endl;
		TRACER << " 1.5RMS rank is " << ranknat << " 1.5 score is " << scorenearnat << " natrms " << natrms << " total beams " << solutions_.size()
			<< std::endl;
		if ( dump_errors_ ) {
			std::ofstream outbeam;
			outbeam.open("rank.txt", std::ofstream::app);
			if ( parallelcount_ == 0 ) {
				cycle -= 1;
			}
			outbeam << "cycle " << cycle << "." << parallelcount_ << " pre best rms and gdt " << pre_best_RMS << " " << pre_best_GDT << " post best " << post_best_RMS << " " << post_best_GDT << " rank "
				<< ranknat << " accepted at " << rank_ << " total " << solutions_.size() <<std::endl;
			outbeam.close();
		}

		if ( pre_best_RMS < 1.2 && post_best_RMS > 1.2 ) {
			if ( !dump_beam_ && !writebeams_ ) {
				core::pose::Pose scpose = pose;
				pose.conformation().chains_from_termini();
				core::pose::PDBInfo newpdbinfo(pose, false);
				newpdbinfo.attach_to(pose.conformation());
				for ( Size j=1; j<=solutions_.size(); j++ ) {
					Real oldRMS;
					oldRMS = solutions_[j].get_rms();
					std::string oldid = solutions_[j].get_id();
					solutions_[j].apply(pose, (int)lower_res, (int)upper_res );
					solutions_[j].apply_sheets(pose);
					pose.dump_pdb("after_error_"+utility::to_string(j)+"_"+nearnativeid+"_"+oldid+"_"+utility::to_string(oldRMS)+".pdb");
					pose = scpose;
				}
				nearnativelps.apply(pose, (int)lower_res, (int)upper_res );
				pose.dump_pdb("after_error_new_"+nearnativeid+"_"+utility::to_string(pre_best_RMS)+".pdb");
				pose = scpose;
			}
			TRACER << " Near natives rejected dumping beam " << std::endl;
		}
		if ( writebeams_ && ((pre_best_RMS < 1.2 && post_best_RMS > 1.2 ) || ( pre_best_RMS < 1.5 && post_best_RMS > 1.5)) ) {
			std::ofstream writeerrors;
			writeerrors.open("errorpoints.txt", std::ofstream::app);
			writeerrors << "error at cycle " << utility::to_string(cycle) << "." << parallelcount_ << " pre and post best rms " << utility::to_string( pre_best_RMS) << " " << utility::to_string( post_best_RMS)
				<< std::endl;
			writeerrors.close();
		}
	}
}

void
LoopPartialSolutionStore::skeleton_filter( core::pose::Pose & pose, DensSkeleton & skeleton, Size start_res, Size stop_res, Size lower_term, Size res_gap ){

	TR << " running skeleton filter " << solutions_.size() << "solutions going in " << std::endl;
	utility::vector1<LoopPartialSolution> new_solutions;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		LoopPartialSolution lps = solutions_[i];
		lps.apply(pose,start_res,stop_res);
		numeric::xyzVector< core::Real > start_xyz = pose.residue(lower_term).atom("C").xyz();
		numeric::xyzVector< core::Real > stop_xyz = pose.residue(lower_term+1).atom("N").xyz();
		//The maximum distance a number of residues could likely stretch
		core::Real max_dist = 4.5*res_gap;
		bool hit_max = false;
		core::Size maxque = 3e4;
		core::Size max_grid = 5;
		core::Real min_dist = skeleton.shortest_path_bfs( start_xyz, stop_xyz, max_dist, maxque, hit_max, max_grid );
		TR << " min dist is " << utility::to_string(min_dist) << std::endl;
		//min_dist == -1 means the the bfs hit the maximum allowed que size and so was unable to complete
		if ( min_dist >= -1 ) {
			new_solutions.push_back(lps);
		}
	}
	//if the new_solutions are empty don't filter at all.
	if ( new_solutions.size() != 0 ) {
		TR << " no solutions exist after the skeleton filter. Keeping the original" << std::endl;
		solutions_ = new_solutions;
	}
	TR << " finished skeleton filter " << solutions_.size() << "solutions coming out " << std::endl;
}

bool
LoopGrower::check_auto_stop(core::pose::Pose & pose, Size lower_res, Size upper_res){

	core::Real density_weight = cen_sf_->get_weight(scoring::elec_dens_fast);
	utility::vector1< core::Real > dens_scores;
	for ( Size i=lower_res; i<=upper_res; i++ ) {
		core::conformation::Residue currentres = pose.residue(i);
		Real resdens = density_weight*-core::scoring::electron_density::getDensityMap().matchResFast( i, currentres, pose, NULL );
		dens_scores.push_back(resdens);
	}
	Real minE = 1e5;
	for ( Size i=1; i<=dens_scores.size(); i++ ) {
		if ( minE > dens_scores[i] ) minE = dens_scores[i];
	}
	Real sum = 0;
	Real sum2 = 0;
	for ( Size i=1; i<=dens_scores.size(); i++ ) {
		Real score = dens_scores[i];
		sum +=(score - minE);
		sum2 += (score -minE)*(score - minE);
	}
	Real mean = sum/dens_scores.size() +minE;
	Real stddev = std::sqrt((dens_scores.size()*sum2 - sum*sum)/dens_scores.size());
	for ( Size i=1; i<=dens_scores.size(); i++ ) {
		Real score = dens_scores[i];
		Real zscore = (score - mean)/stddev;
		TRACER << "zscore of residue " << lower_res+i-1 << " is " << zscore << " dens score is " << score << std::endl;
	}

	return false;
}

void
LoopPartialSolutionStore::zscoretransform(){
	Real minE = solutions_[1].get_score();
	Size N = 0;
	Real sum = 0;
	Real sum2 = 0;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		N++;
		Real score = solutions_[i].get_score();
		sum += (score - minE);
		sum2 += (score - minE)*(score - minE);
	}
	Real mean = sum/N + minE;
	Real stddev = std::sqrt(N*sum2 - sum*sum)/N;
	TRACER << "mean and stddev" << mean << " " << stddev << " sum 2 " << sum2 << std::endl;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		Real score = solutions_[i].get_score();
		Real zscore = (score - mean)/stddev;
		TRACER << "zscore of  " << i << " is " << zscore << std::endl;
		solutions_[i].set_score(zscore);
	}
}
void
LoopPartialSolutionStore::cluster_check(LoopPartialSolution nearnativelps, core::pose::Pose& pose, Size total_lower, Size lower_res, Size upper_res ){

	pose.conformation().chains_from_termini();
	core::pose::PDBInfo newpdbinfo(pose, false);
	newpdbinfo.attach_to(pose.conformation());
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		LoopPartialSolution beam = solutions_[i];
		Real rmsdist = nearnativelps.partialrms(beam, fragmelt_, total_lower, false, 1);
		Real masterbeamrmsdist = nearnativelps.partialrms(beam, fragmelt_, total_lower, true, rmswindow_);
		Real score = beam.get_score();
		if ( rmsdist < rmscutoff_ ) {
			beam.apply( pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("afterclustercheck_"+utility::to_string(score)+"_localrmsdist_"+utility::to_string(rmsdist)+"_.pdb");
		}
		if ( masterbeamrmsdist <= master_beam_cutoff_ && rmsdist > rmscutoff_ ) {
			beam.apply( pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("afterclustercheck_"+utility::to_string(score)+"_masterbeamwidth_"+utility::to_string(masterbeamrmsdist)+"_.pdb");
		}
	}
	nearnativelps.apply(pose, (int)lower_res, (int)upper_res);
	Real score = nearnativelps.get_score();
	Real nnrms = nearnativelps.get_rms();
	pose.dump_pdb("afterclustercheck_"+utility::to_string(score)+"_native_"+utility::to_string(nnrms)+"_.pdb");
}
void
LoopPartialSolutionStore::parametercheck_filter(core::pose::Pose& pose, Size fragmelt, Size total_lower, Size lower_res, Size upper_res, Size rmswindow, bool dump_beam,
	Size& totalreqbeam, Size& totalreqmaxbeam, Real& totalreqscorecut, bool writebeams, Size cycle, Size parallelcount ){

	LoopPartialSolutionStore filteredsolutions;
	sort();
	Real bestscore = solutions_[1].get_score();
	TRACER << "The best score is " << bestscore << std::endl;

	Real pre_best_RMS = solutions_[1].get_rms();
	LoopPartialSolution nearnativelps = solutions_[1];
	LoopPartialSolution nativeranklps = solutions_[1];
	Size rank =999999;
	//std::map< Size, Real> bestsolutions;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		Real oldRMS = solutions_[i].get_rms();
		if ( oldRMS <= pre_best_RMS ) {
			pre_best_RMS = oldRMS;
			nearnativelps = solutions_[i];
		}
	}
	if ( master_beam_width_ == 1337 && master_beam_cutoff_ == 1337 ) {
		sortrms();
		for ( Size i=1; i<=solutions_.size(); i++ ) {
			LoopPartialSolution incominglps = solutions_[i];
			bool addlps = true;
			for ( Size j=1; j<=filteredsolutions.size(); j++ ) {
				LoopPartialSolution storedlps = filteredsolutions[j];
				Real rmsdist = incominglps.partialrms(storedlps, fragmelt, total_lower, false, 1);
				if ( rmsdist < rmscutoff_ ) {
					addlps = false;
				}
			}
			if ( addlps ) filteredsolutions.store(incominglps);
			if ( filteredsolutions.size() >= totalreqbeam ) break;
		}
		solutions_ = filteredsolutions.get_solutions();
		return;
	}

	for ( Size i=1; i<=solutions_.size(); i++ ) {
		LoopPartialSolution rmscheck = solutions_[i];
		Real rmsdist = nearnativelps.partialrms(rmscheck, fragmelt, total_lower, false, 1);
		if ( rmsdist <= rmscutoff_ ) {
			rank = i;
			nativeranklps = solutions_[i];
			TRACER << "rmsdist " << rmsdist << std::endl;
			if ( dump_beam ) {
				nearnativelps.apply(pose, (int)lower_res, (int)upper_res );
				pose.dump_pdb("after_natonly_" + nearnativelps.get_id()+".pdb");
				nativeranklps.apply(pose, (int)lower_res, (int)upper_res );
				pose.dump_pdb("after_natrank_" + nativeranklps.get_id()+".pdb");
			}
			break;
		}
	}
	if ( master_beam_width_ == 1337 ) {
		solutions_.clear();
		solutions_.push_back(nativeranklps);
		TRACER << "accepting only the nearnative. pre_best_RMS is " << pre_best_RMS << " nr rms is " << nativeranklps.get_rms() << " it's rank is " << rank << std::endl;
		if ( dump_beam ) {
			nearnativelps.apply(pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("after_natonly_" + nearnativelps.get_id()+".pdb");
			nativeranklps.apply(pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("after_natrank_" + nativeranklps.get_id()+".pdb");
		}
		return;
	}
	Real nativerankscore = nativeranklps.get_score();
	Real scorecutoff = -(bestscore-nativerankscore)/(upper_res-lower_res);

	Size reqmaxbeam = 1;
	rank += 30;
	if ( rank > solutions_.size() ) rank = solutions_.size();

	TRACER << "master beam width " << master_beam_width_ << std::endl;
	for ( Size i=1; i<=rank; i++ ) {
		LoopPartialSolution lps = solutions_[i];
		bool addlps = true;
		Size maxbeamcount = 0;
		for ( Size j=1; j<=filteredsolutions.size(); j++ ) {
			LoopPartialSolution filteredlps = filteredsolutions[j];
			Real partialrms = lps.partialrms(filteredlps, fragmelt, total_lower, false, 1 );
			Real masterbeamrms = lps.partialrms(filteredlps, fragmelt, total_lower, false, rmswindow );
			if ( partialrms < rmscutoff_ ) {
				addlps = false;
				TRACER << lps.get_id() << " is rejected " << filteredlps.get_id() << " is too similar with a distance of " << partialrms << std::endl;
				break;
			}
			if ( masterbeamrms <= master_beam_cutoff_ ) {
				maxbeamcount++;
			}
			if ( maxbeamcount >= master_beam_width_ ) {
				addlps = false;
				TRACER << lps.get_id() << " is rejected. Master beam width reached " << std::endl;
				break;
			}
		}
		if ( i==rank ) reqmaxbeam = maxbeamcount;
		if ( addlps ) {
			filteredsolutions.store(lps);
			TRACER << "adding " << lps.get_id() << std::endl;
		}
	}
	solutions_ = filteredsolutions.get_solutions();

	Real post_best_RMS = solutions_[1].get_rms();
	std::string nearnativeid = nearnativelps.get_id();
	if ( dump_beam ) {
		for ( Size j=1; j<=solutions_.size(); j++ ) {
			Real oldRMS = solutions_[j].get_rms();
			std::string oldid = solutions_[j].get_id();
			solutions_[j].apply(pose, (int)lower_res, (int)upper_res );
			pose.dump_pdb("after_filter_"+utility::to_string(cycle)+"_"+oldid+"_"+utility::to_string(oldRMS)+".pdb");
		}
		nearnativelps.apply(pose, (int)lower_res, (int)upper_res );
		pose.dump_pdb("after_filter_nn_"+utility::to_string(cycle)+"_"+nearnativeid+"_"+utility::to_string(pre_best_RMS)+".pdb");
		nativeranklps.apply(pose, (int)lower_res, (int)upper_res );
		std::string nativerankid = nativeranklps.get_id();
		Real nativerankrms = nativeranklps.get_rms();
		pose.dump_pdb("after_filter_nr_"+utility::to_string(cycle)+"_"+nativerankid+"_"+utility::to_string(nativerankrms)+".pdb");
	}

	rank = 0;
	bool rankset = false;
	for ( Size i=1; i<=solutions_.size(); i++ ) {
		Real oldRMS = solutions_[i].get_rms();
		//this will set the rank to the first structure to be lower than 1.5 RMS
		if ( oldRMS <=1.5 && !rankset ) {
			rank = i;
			rankset = true;
		}
		if ( oldRMS <= post_best_RMS ) {
			post_best_RMS = oldRMS;
			if ( !rankset ) {
				rank = i;
			}
		}
	}

	if ( totalreqbeam < filteredsolutions.size() ) totalreqbeam = filteredsolutions.size();
	//if((totalreqdist < reqdist) && (reqdist != 999)) totalreqdist = reqdist;
	if ( totalreqmaxbeam < reqmaxbeam ) totalreqmaxbeam = reqmaxbeam;
	if ( totalreqscorecut < scorecutoff ) totalreqscorecut = scorecutoff;

	Real nativerankrms = nativeranklps.get_rms();
	TRACER << "pre_best_RMS " << nearnativeid << " " << pre_best_RMS << " nativerankrms " << nativerankrms << " post_best_RMS " << post_best_RMS << std::endl;
	TRACER << "The necessary parameters for this step are: " <<  std::endl;
	//TRACER << "rmscutoff = " << reqdist << " overall = " << totalreqdist << std::endl;
	TRACER << "masterbeamwidth = " << reqmaxbeam << " overall = " << totalreqmaxbeam << std::endl;
	TRACER << "score cutoff = " << scorecutoff << " overall = " << totalreqscorecut << std::endl;
	TRACER << "totalbeam = " << filteredsolutions.size() << " overall = " << totalreqbeam << std::endl;
	if ( dump_errors_ ) {
		std::ofstream outbeam;
		outbeam.open("rank.txt", std::ofstream::app);
		outbeam << "cycle " << cycle << "." << parallelcount << " pre and post best rms " << pre_best_RMS << " " << post_best_RMS << " rank " << rank << std::endl;
		outbeam.close();
	}

	if ( (pre_best_RMS < 1.2 && post_best_RMS > 1.3) || (solutions_.size()>maxelts_) ) {
		if ( !dump_beam && !writebeams ) {
			for ( Size j=1; j<=solutions_.size(); j++ ) {
				Real oldRMS;
				oldRMS = solutions_[j].get_rms();
				std::string oldid = solutions_[j].get_id();
				solutions_[j].apply(pose, (int)lower_res, (int)upper_res);
				pose.dump_pdb("after_error_"+utility::to_string(j)+"_"+nearnativeid+"_"+oldid+"_"+utility::to_string(oldRMS)+".pdb");
			}
			nearnativelps.apply(pose, (int)lower_res, (int)upper_res);
			pose.dump_pdb("after_error_new_"+nearnativeid+"_"+utility::to_string(pre_best_RMS)+".pdb");
			nativeranklps.apply(pose, (int)lower_res, (int)upper_res);
			std::string nativerankid = nativeranklps.get_id();
			Real nativerankrms = nativeranklps.get_rms();
			pose.dump_pdb("after_error_nr_"+nativerankid+"_"+utility::to_string(nativerankrms)+".pdb");
		}
		TRACER << " Near natives rejected dumping beam " << std::endl;
	}
	if ( (post_best_RMS >= 3.0 && !writebeams) || solutions_.size()>maxelts_ ) {
		TRACER << "no natives remain exiting protocol or they score too poorly" << std::endl;
		if ( writebeams ) {
			std::ofstream outbeam;
			outbeam.open("finished.txt", std::ofstream::app);
			outbeam << "finished" << std::endl;
			outbeam.close();
		}
		exit(1);
	}
}

void LoopGrower::read_coordfile(){
	std::ifstream coordstream;
	coordstream.precision(8);
	if ( coordfile_ == "" ) {
		TRACER << " No coordfile name provided " << std::endl;
		return;
	}
	coordstream.open(coordfile_.c_str());
	core::Real scoreweight;
	bool firstline = true;
	if ( coordstream.peek() == std::ifstream::traits_type::eof() ) {
		TRACER << "THE COORDINATE FILE YOU SPECIFIED IS EMPTY OR DOES NOT EXIST. EITHER PROVIDE CORRECT COORDINATES OR REMOVE THE FILENAME FROM THE OPTIONS TO RUN WITHOUT THEM" << std::endl;
		exit(0);
	}
	while ( coordstream.peek() != std::ifstream::traits_type::eof() ) {
		if ( firstline ) {
			coordstream >> scoreweight;
			firstline = false;
		}
		if ( coordstream.peek() == std::ifstream::traits_type::eof() ) break;
		Size res = 0;
		Size atom = 0;
		Size totalpositions = 0;
		coordstream >> res >> atom >> totalpositions;
		utility::vector1< numeric::xyzVector< core::Real > > positions;
		for ( Size i=1; i<=totalpositions; i++ ) {
			if ( coordstream.peek() == std::ifstream::traits_type::eof() ) break;
			numeric::xyzVector<core::Real> position(0,0,0);
			coordstream >> position[0] >> position[1] >> position[2];
			positions.push_back(position);
			//TRACER << res << "." << atom << " " << i << " " << position[0] << " " << position[1] << " " << position[2] << std::endl;
			if ( coordstream.peek() == std::ifstream::traits_type::eof() ) break;
		}
		if ( coordstream.peek() == std::ifstream::traits_type::eof() ) break;
		std::pair< Size, Size > atomid = std::make_pair(res,atom);
		scoringcoords_[atomid] = positions;
	}
	/*std::map< std::pair<Size, Size>, utility::vector1< numeric::xyzVector< core::Real > > >::iterator iter;
	Size count = 1;
	for( iter = scoringcoords_.begin(); iter != scoringcoords_.end(); iter++){
	count++;
	for( Size i=1; i<=iter->second.size(); i++){
	TRACER << " setting " << iter->first.first << "." << iter->first.second << " " << iter->second.size() << " " << iter->second[i][0] << " " << iter->second[i][1] << " " << iter->second[i][2] << std::endl;
	}
	}*/
}

bool
LoopGrower::is_beta(Real phi, Real psi){
	if ( phi > -180 && phi <= 0 && psi > 45 && psi < 180 ) {
		return true;
	} else {
		return false;
	}
}

void
LoopGrower::coordinate_filter(LoopPartialSolutionStore& solutionset, core::pose::Pose pose, bool lower, Size lower_fasta, Size upper_fasta, Size rangelo, Size rangehi){

	bool coordmatch = false;
	//find the cutpoint
	Size cutpoint = rangelo;
	while ( !pose.fold_tree().is_cutpoint(cutpoint) ) {
		cutpoint++;
		if ( cutpoint > rangehi ) {
			cutpoint = rangelo;
			break;
		}
	}
	while ( pose.residue(cutpoint).name3() == "XXX" ) cutpoint--;
	//get the range via the cutpoint
	Size lowerrange;
	Size upperrange;
	Size fasta_start;
	if ( lower ) {
		lowerrange = cutpoint-4;
		upperrange = cutpoint;
		fasta_start = lower_fasta-5;
	} else {
		lowerrange = cutpoint+1;
		upperrange = cutpoint+6;
		fasta_start = upper_fasta+1;
	}
	if ( lowerrange < rangelo ) {
		TRACER << " setting lower range to torsionrangelo " << std::endl;
		lowerrange = rangelo;
	}
	if ( upperrange > rangehi ) {
		TRACER << " setting upper range to torsionrangehi " << std::endl;
		upperrange = rangehi;
	}
	LoopPartialSolutionStore matchingsolutions = solutionset;
	matchingsolutions.clear();
	//check if any solutions match coordinates
	for ( Size i=1; i<=solutionset.size(); i++ ) {
		LoopPartialSolution lps = solutionset[i];
		lps.apply(pose,rangelo,rangehi);
		Size radiuscluster = check_coordinates(pose, lowerrange, upperrange, fasta_start, 3.0);
		if ( radiuscluster > 4 ) {
			coordmatch = true;
			matchingsolutions.store(lps);
		}

	}
	if ( coordmatch ) {
		solutionset = matchingsolutions;
	}
	//if one solution is a match filter all those that don't match.

}

void LoopGrower::write_to_disk( LoopPartialSolutionStore solutionset, Size step, Size added_lower, Size added_upper, bool filteronly, bool lower){

	std::ofstream outbeam;
	outbeam.precision(12);
	std::string filename = "beam_"+utility::to_string(step)+"."+utility::to_string(parallelcount_)+".txt";
	if ( filteronly ) filename = "beam_0.txt";
	if ( !filteronly ) {
		outbeam.open (filename.c_str());
	} else {
		outbeam.open (filename.c_str());
	}
	LoopPartialSolution outlps;
	for ( Size i=1; i<= solutionset.size(); i++ ) {
		outlps = solutionset[i];
		outbeam << outlps.size() << " " << added_lower << " " << added_upper << " " << step << " " << lower << std::endl;;
		outlps.write_beam(outbeam);
	}
	outbeam.close();

}

void LoopGrower::read_from_disk(LoopPartialSolutionStore & solutionset, int & cycle, bool & lower, bool filterbeams ){
	std::ifstream beamread;
	std::string beamfile;
	if ( filterbeams ) {
		beamfile = filterbeams_;
	} else {
		beamfile = storedbeams_;
	}
	beamread.open(beamfile.c_str());
	if ( beamread.peek() == std::ifstream::traits_type::eof() ) {
		TRACER << " The beam file " << beamfile << " is empty exiting the protocol " << std::endl;
		runtime_assert( beamread.peek() != std::ifstream::traits_type::eof() );
	}
	TRACER << "reading beams from " << beamfile << std::endl;
	// Size countertest = 0;
	while ( !beamread.eof() ) {
		Size rescount;
		beamread >> rescount >> storelow_ >> storehi_ >> cycle >> lower;
		if ( beamread.eof() ) break;
		LoopPartialSolution storedbeam;
		for ( Size j=1; j<=rescount; j++ ) {
			Real phi, psi, omega, chi1, chi2, chi3, chi4;
			phi = psi = omega = chi1 = chi2 = chi3 = chi4 = 0;
			Size nchi = 0;
			beamread >> phi >> psi >> omega >> nchi;
			if ( nchi >=1 ) beamread >> chi1;
			if ( nchi >=2 ) beamread >> chi2;
			if ( nchi >=3 ) beamread >> chi3;
			if ( nchi >=4 ) beamread >> chi4;
			ResTorsions residue(phi, psi, omega, nchi, chi1, chi2, chi3, chi4);
			storedbeam.push_back_restorsions(residue);
			if ( beamread.eof() ) break;
		}
		for ( Size j=1;  j<=rescount; j++ ) {
			Size atomcount;
			utility::vector1< core::Vector > resbackbone;
			beamread >> atomcount;
			if ( beamread.eof() ) break;
			for ( Size k=1; k<=atomcount; k++ ) {
				numeric::xyzVector<core::Real> atomxyzvector;
				beamread >> atomxyzvector[0] >> atomxyzvector[1] >> atomxyzvector[2];
				resbackbone.push_back(atomxyzvector);
				if ( beamread.eof() ) break;
			}
			storedbeam.push_back_backbone(resbackbone);
		}
		Size sheets;
		Real bonus_score;
		beamread >> sheets >> bonus_score;
		storedbeam.set_bonus_score(bonus_score);
		for ( Size i=1; i<=sheets; i++ ) {
			Size baseres, jumpid;
			beamread >> baseres >> jumpid;
			numeric::xyzMatrix< Real > rotation;
			numeric::xyzVector< Real > translation;
			for ( Size ii=1; ii<=3; ii++ ) {
				for ( Size j=1; j<=3; j++ ) {
					Real rotationdata;
					beamread >> rotationdata;
					rotation(ii,j) = rotationdata;
					if ( beamread.eof() ) break;
				}
				if ( beamread.eof() ) break;
			}
			beamread >> translation[0] >> translation[1] >> translation[2];
			Size sheetrescount;
			beamread >> sheetrescount;
			utility::vector1< ResTorsions > sheetresidues;
			for ( Size ii=1; ii<=sheetrescount; ii++ ) {
				Real phi, psi, omega, chi1, chi2, chi3, chi4;
				phi = psi = omega = chi1 = chi2 = chi3 = chi4 = 0;
				Size nchi = 0;
				beamread >> phi >> psi >> omega >> nchi;
				if ( nchi >=1 ) beamread >> chi1;
				if ( nchi >=2 ) beamread >> chi2;
				if ( nchi >=3 ) beamread >> chi3;
				if ( nchi >=4 ) beamread >> chi4;
				ResTorsions residue(phi, psi, omega, nchi, chi1, chi2, chi3, chi4);
				sheetresidues.push_back(residue);
				if ( beamread.eof() ) break;
			}
			SheetPositions sheet( sheetresidues, rotation, translation, jumpid, baseres);
			storedbeam.push_back_sheet(sheet);
		}
		Real score, rms, gdt;
		beamread >> score >> rms >> gdt;
		storedbeam.set_rms(rms);
		storedbeam.set_gdt(gdt);
		storedbeam.set_score(score);
		solutionset.store(storedbeam);
	}
	TRACER << "finished reading" << beamfile << std::endl;
	beamread.close();

}
void
LoopGrower::update_fragment_library_pointers() {
	core::Size nfragsets = fragments_.size();

	// map positions to fragments
	library_.resize( nfragsets );
	for ( int i=1; i<=(int)nfragsets; ++i ) {
		for ( core::fragment::ConstFrameIterator j = fragments_[i]->begin(); j != fragments_[i]->end(); ++j ) {
			core::Size position = (*j)->start();
			library_[i][position] = **j;
		}
	}
}
Real
LoopGrower::GDThatonative( core::pose::Pose const &pose, int natlow, int nathi, int natstoplow, int natstarthi, int startlow, int stoplow, int starthi, int stophi){

	if ( !native_ ) {
		TRACER << "no native exists cannot determine GDTha " << std::endl;
		runtime_assert( native_);
	}

	Size GDTha = 0;
	int rescount = startlow;
	if ( (resstart_ == 0 || pose.fold_tree().is_cutpoint(startlow-1)) ) {
		rescount-=1;
	}
	Real looplength = 0;
	std::string atomnames[5] = {"CA","O","C","N","CB"};
	for ( int i=natlow; i<=nathi; i++ ) {
		if ( (i > natstoplow) && (i < natstarthi) ) continue;
		rescount++;
		if ( (rescount > stoplow) && (rescount < starthi) ) continue;
		if ( (rescount > stophi) && (stophi !=0) ) break;

		for ( Size ii=0; ii < 5; ii++ ) {
			std::string atomname = atomnames[ii];
			if ( pose.residue(rescount).name3() == "GLY" && atomname == "CB" ) {
				continue;
			}
			core::Vector newatom = pose.residue(rescount).atom(atomname).xyz();
			core::Vector natatom = native_->residue(i).atom(atomname).xyz();
			atomGDTha(natatom, newatom, GDTha);
			looplength++;
		}
	}
	Real GDThatotal = GDTha/(looplength*4);
	return GDThatotal;
}

void
LoopGrower::atomGDTha(core::Vector atom1, core::Vector atom2, Size &GDTha){
	Real distance = (atom1-atom2).length_squared();
	if ( distance < 0.5 ) GDTha++;
	if ( distance < 1.0 ) GDTha++;
	if ( distance < 2.0 ) GDTha++;
	if ( distance < 4.0 ) GDTha++;
}

Real
LoopGrower::RMStonative( core::pose::Pose const &pose, int natlow, int nathi, int natstoplow, int natstarthi, int startlow, int stoplow, int starthi, int stophi){

	if ( !native_ ) {
		TRACER << "no native exists cannot determine RMS " << std::endl;
		runtime_assert( !native_ );
	}


	Real err = 0;
	int rescount = startlow;
	if ( (resstart_ == 0 || pose.fold_tree().is_cutpoint(startlow-1)) ) {
		rescount-=1;
	}
	Real looplength = 0;
	for ( int i=1; i<=nathi; i++ ) {}
	for ( int i=natlow; i<=nathi; i++ ) {
		if ( (i > natstoplow) && (i < natstarthi) ) continue;
		rescount++;
		if ( (rescount > stoplow) && (rescount < starthi) ) continue;
		if ( (rescount > stophi) && (stophi !=0) ) break;
		Real distance = 0;
		core::Vector ca_i = pose.residue(rescount).atom("CA").xyz();
		core::Vector natca_i = native_->residue(i).atom("CA").xyz();
		distance = (ca_i-natca_i).length_squared();
		err += distance;
		core::Vector o_i = pose.residue(rescount).atom("O").xyz();
		core::Vector nato_i = native_->residue(i).atom("O").xyz();
		distance = (o_i-nato_i).length_squared();
		err += distance;
		core::Vector c_i = pose.residue(rescount).atom("C").xyz();
		core::Vector natc_i = native_->residue(i).atom("C").xyz();
		distance = (c_i-natc_i).length_squared();
		err += distance;
		if ( pose.residue(rescount).name3()!="GLY" ) {
			core::Vector cb_i = pose.residue(rescount).atom("CB").xyz();
			core::Vector natcb_i = native_->residue(i).atom("CB").xyz();
			distance = (cb_i-natcb_i).length_squared();
			err += distance;
			looplength+=1;
		}
		core::Vector n_i = pose.residue(rescount).atom("N").xyz();
		core::Vector natn_i = native_->residue(i).atom("N").xyz();
		distance = (n_i-natn_i).length_squared();
		err += distance;
		looplength +=4;
	}
	Real RMS = std::sqrt(err/looplength);
	return RMS;
}

Size
LoopGrower::check_coordinates(core::pose::Pose& pose, Size lower_pose, Size upper_pose, Size lower_fasta, Real radius){

	Size radiuscount = 0;
	std::string atomnames[1] = {"CA"};
	Size fastares = lower_fasta;
	for ( Size i=lower_pose; i<=upper_pose; i++ ) {
		bool inradius = false;
		for ( Size j=0; j<1; j++ ) {
			core::Vector atomcoords = pose.residue(i).atom(atomnames[j]).xyz();
			core::Size atomnumber = pose.residue(i).atom_index(atomnames[j]);
			std::pair< Size, Size > atomid = std::make_pair(fastares,atomnumber);
			for ( Size k=1; k<=scoringcoords_[atomid].size(); k++ ) {
				core::Vector storedcoords = scoringcoords_[atomid][k];
				core::Real distance = (atomcoords-storedcoords).length_squared();
				if ( distance < radius ) inradius = true;
			}
		}
		if ( inradius ) radiuscount++;
		fastares++;
	}
	return radiuscount;

}

void LoopGrower::update_to_stored( core::pose::Pose& growpose, core::pose::Pose& growpose_cen, const core::chemical::ResidueTypeCOPs& restypes_pose,
	const core::chemical::ResidueTypeCOPs& restypes_pose_cen, int & lower_pose, int & upper_pose, int & lower_fasta, int & upper_fasta, Size newreslow,
	Size newreshi, bool is_nterm, bool is_cterm){

	for ( int i=1; i<=(int)newreslow; ++i ) {
		core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[lower_fasta+i-1] ));
		growpose.conformation().safely_append_polymer_residue_after_seqpos(*newres, lower_pose+i-1, true);
		core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[lower_fasta+i-1] ));
		growpose_cen.conformation().safely_append_polymer_residue_after_seqpos(*newres_cen, lower_pose+i-1, true);
	}
	lower_pose +=newreslow;
	upper_pose +=newreslow;
	lower_fasta += newreslow;
	for ( int i=1; i<=(int)newreshi; ++i ) {
		core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[upper_fasta-i+1] ));
		growpose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, lower_pose+1, true);
		core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[upper_fasta-i+1] ));
		growpose_cen.conformation().safely_prepend_polymer_residue_before_seqpos(*newres_cen, lower_pose+1, true);
	}
	upper_pose += newreshi;
	upper_fasta -= newreshi;
	if ( upper_fasta<=lower_fasta ) {
		if ( !is_nterm && !is_cterm ) {
			core::conformation::remove_upper_terminus_type_from_conformation_residue( growpose.conformation(), lower_pose );
			core::conformation::remove_lower_terminus_type_from_conformation_residue( growpose.conformation(), lower_pose+1 );
			core::conformation::remove_upper_terminus_type_from_conformation_residue( growpose_cen.conformation(), lower_pose );
			core::conformation::remove_lower_terminus_type_from_conformation_residue( growpose_cen.conformation(), lower_pose+1 );

			growpose.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");
			growpose_cen.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");

			core::pose::remove_variant_type_from_pose_residue( growpose, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( growpose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			core::pose::remove_variant_type_from_pose_residue( growpose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( growpose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
		}
		growpose_cen.conformation().chains_from_termini();
		growpose.conformation().chains_from_termini();
		core::pose::PDBInfo newpdbinfocen(growpose_cen, false);
		core::pose::PDBInfo newpdbinfo(growpose, false);
		newpdbinfocen.attach_to(growpose_cen.conformation());
		newpdbinfo.attach_to(growpose.conformation());
	}
	add_user_csts(growpose);
	add_user_csts(growpose_cen);
}
void LoopGrower::rescoresolutionset( LoopPartialSolutionStore& solutionset, core::pose::Pose& fa_pose, core::pose::Pose& cen_pose, Size torsionrangelo, Size torsionrangehi ){
	LoopPartialSolutionStore rescoredset = solutionset;
	rescoredset.clear();
	core::pose::Pose scpose = fa_pose;
	core::pose::Pose sheetpose = cen_pose;
	for ( Size i=1; i<=solutionset.size(); i++ ) {
		LoopPartialSolution lps = solutionset[i];
		utility::vector1< std::pair<Size,Size> > history = lps.get_history();
		lps.apply(cen_pose, torsionrangelo, torsionrangehi);
		lps.apply(fa_pose, torsionrangelo, torsionrangehi);
		Real RMS = lps.get_rms();
		Real GDT = lps.get_gdt();
		Real bonus_score = lps.get_bonus_score();
		if ( samplesheets_ ) lps.apply_sheets(cen_pose);
		if ( checksymm_ ) transform_to_closest_symmunit(cen_pose,fa_pose,torsionrangelo);
		//turn minimize off if you want to rescore without repacking
		if ( (pack_min_cycles_ !=0 || cenrot_) ) {
			refine_cycle( fa_pose, cen_pose, torsionrangelo, torsionrangehi, false, i, i, 0, 0);
			//only apply sheets if they aren't already added in the refinement cycle
		} else if ( samplesheets_ ) {
			lps.apply_sheets(fa_pose);
		}
		Real score = modifiedscore( fa_pose, cen_pose, torsionrangelo, torsionrangehi );
		//turning this bonus off when sample sheets isn't active so I can do rescore testing
		if ( samplesheets_ ) {
			score += bonus_score;
		}

		//update the torsions after refinecycle
		if ( pack_min_cycles_ !=0 || cenrot_ ) {
			LoopPartialSolution falps(fa_pose, torsionrangelo, torsionrangehi, score);
			falps.set_rms(RMS);
			falps.set_gdt(GDT);
			if ( trackfragments_ ) falps.set_history(history);
			falps.set_score(score);
			if ( samplesheets_ ) {
				//Storing the sheets needs to come before strand_score which removes them to calculate the sheet score bonus
				store_sheets(fa_pose, falps);
				Real strand_score = sheetscore(fa_pose, cen_pose, torsionrangelo, torsionrangehi);
				if ( continuous_sheets_ ) {
					falps.set_bonus_score(bonus_score+strand_score);
				} else {
					falps.set_bonus_score(0);
				}
			}
			rescoredset.store(falps);
		} else {
			lps.set_score(score);
			rescoredset.store(lps);
		}

		TRACER << " score for " << lps.get_id() << " " << lps.get_rms() << " is " << score << std::endl;
		if ( debug_ && fafilter_ ) fa_pose.dump_pdb("afterfafilter_"+utility::to_string(lps.get_score())+"_"+utility::to_string(RMS)+".pdb");
		lps.set_score(score);
		if ( dumpbeam_ ) {
			if ( pack_min_cycles_ != 0 || cenrot_ ) {
				fa_pose.conformation().chains_from_termini();
				core::pose::PDBInfo newpdbinfo(fa_pose, false);
				newpdbinfo.attach_to(fa_pose.conformation());
				fa_pose.dump_pdb("afterrescore_"+utility::to_string(score)+"_"+utility::to_string(RMS)+".pdb");
			} else {
				cen_pose.conformation().chains_from_termini();
				core::pose::PDBInfo newpdbinfo(cen_pose, false);
				newpdbinfo.attach_to(cen_pose.conformation());
				cen_pose.dump_pdb("afterrescore_"+utility::to_string(score)+"_"+utility::to_string(RMS)+".pdb");
			}
		}
		fa_pose = scpose;
		cen_pose = sheetpose;
	}
	if ( dumperrors_ ) {
		Real bestrms = 99999;
		Size rank;
		Real bestscore;
		rescoredset.sort();
		bool reportednn = false;
		for ( Size i=1; i<=rescoredset.size(); i++ ) {
			LoopPartialSolution lps = rescoredset[i];
			lps = rescoredset[i];
			Real rms = lps.get_rms();
			if ( !reportednn && rms <=1.5 ) {
				TRACER << " the nnrms is " << rms << " scoring " << lps.get_score() << " rank " << i << std::endl;
				reportednn = true;
			}
			if ( rms < bestrms ) {
				bestrms = rms;
				rank = i;
				bestscore = lps.get_score();
			}
		}
		TRACER << "the best rms after rescore is " << bestrms << " and it's rank is " << rank << " with a score of " << bestscore << std::endl;
	}
	solutionset = rescoredset;
}

void
LoopGrower::store_sheets( core::pose::Pose& pose, LoopPartialSolution& lps ){

	if ( pose.total_residue() > total_residues_ ) {
		Size startpos = total_residues_+1;
		Size endpos = total_residues_+sheetsize_;
		Size baseres = 0;
		for ( Size i=startpos; i<=endpos; i++ ) {
			if ( pose.fold_tree().is_jump_point(i) ) {
				baseres = pose.fold_tree().get_parent_residue(i);
			}
		}
		runtime_assert( baseres != 0 );
		lps.store_sheetpositions( pose, startpos, endpos, numjumps_+1, baseres );
	}
	if ( pose.total_residue() > total_residues_+sheetsize_ ) {
		TR << " storing second " << std::endl;
		Size startpos = total_residues_+sheetsize_+1;
		Size endpos = startpos+sheetsize_-1;
		Size baseres = 0;
		for ( Size i=startpos; i<=endpos; i++ ) {
			if ( pose.fold_tree().is_jump_point(i) ) {
				baseres = pose.fold_tree().get_parent_residue(i);
			}
		}
		runtime_assert( baseres != 0 );
		lps.store_sheetpositions( pose, startpos, endpos, numjumps_+2, baseres );
		TR << " stored two sheets " << std::endl;
	}
}

Real
LoopGrower::modifieddensity( core::pose::Pose& pose, Size rangelo, Size rangehi, Real density_weight, Size & includesheets){

	//This vector scores the original scores of the backbones as well as the modified scores for the side chain to use in the modified dens score later.
	utility::vector1< utility::vector1< Real> > residueatomdensities;
	utility::vector1< Real > worstatoms;
	Size res_windowlength = 8; //it may be  worthwhile to make this something that can be set in the xml don't forget to also put this in the non window_dens score
	Real worst_atom = -9999;
	Real total_fast_dens = 0;
	Size atomcounter = 0;

	//this is just for a test
	//(*sf_)(pose);

	//Store the modified fast densities of the relative atoms
	Real sc_dens = 0;
	for ( core::Size i = rangelo; i <= pose.total_residue(); i++ ) {
		if ( i > rangehi && i <= total_residues_ ) continue;
		utility::vector1< Real > resatomdens;
		utility::vector1< Real > scatomdens;
		Real bbresiduedens = 0;
		Real residueworstatom = -999999;
		Real worstscatom = -999999;
		core::conformation::Residue currentres = pose.residue(i);
		Size resatomcount = currentres.nheavyatoms();
		for ( Size j=1; j<=resatomcount; j++ ) {
			std::string atomname = currentres.atom_name(j);
			if ( atomname.find("CEN") != std::string::npos ) {
				resatomcount-=1;
				continue;
			}
			atomcounter++;
			Real atomdens;// = density_weight*-core::scoring::electron_density::getDensityMap().matchAtomFast( i, j, currentres, pose, NULL );
			if ( skeleton_file_ != "" ) {
				atomdens = density_weight*-skeleton_.matchAtomFast( i, j, currentres, pose );
			} else {
				atomdens = density_weight*-core::scoring::electron_density::getDensityMap().matchAtomFast( i, j, currentres, pose, NULL );
			}
			total_fast_dens += atomdens;
			//track the backbone and sc density seperately
			if ( j<=4 ) {
				bbresiduedens += atomdens;
				resatomdens.push_back(atomdens);
			}
			if ( j>4 ) scatomdens.push_back(atomdens);
			//find the worst backbone atoms and check CB and CG for side chain to use for continuous scaling
			if ( j<=3 && pose.residue(i).name3() != "replacewith(GLY)toIgnoreGlycine" ) {
				if ( atomdens > worst_atom ) {
					worst_atom = atomdens;
				}
				if ( atomdens > residueworstatom ) {
					residueworstatom = atomdens;
				}
			} else if ( atomname == " CB " || atomname == " CG " ) {
				if ( atomdens > worstscatom ) {
					worstscatom = atomdens;
				}
			}
		}
		//add side chain density
		for ( Size j=1; j<= scatomdens.size(); j++ ) {
			if ( worstscatom > scatomdens[j] ) {
				sc_dens += worstscatom*sc_scale_;
				resatomdens.push_back(worstscatom*sc_scale_);
			} else {
				sc_dens += scatomdens[j]*sc_scale_;
				resatomdens.push_back(scatomdens[j]*sc_scale_);
			}
		}
		worstatoms.push_back( residueworstatom );
		residueatomdensities.push_back(resatomdens);
	}

	//Which if any sheets to keep
	Real mod_dens = 0;
	bool keepfirstsheet = true;
	bool keepsecondsheet = true;
	worst_atom = -99999;
	Size nonsheetrescount = rangehi-rangelo+2;
	for ( Size i=nonsheetrescount-sheetsize_+1; i<=residueatomdensities.size(); i++ ) {
		Real worst = worstatoms[i];
		if ( worst > worst_atom && i <= nonsheetrescount ) worst_atom = worst;
		if ( worst > worst_atom*sheet_tolerance_ && i > nonsheetrescount && i <= nonsheetrescount+sheetsize_ ) keepfirstsheet = false;
		if ( worst > worst_atom*sheet_tolerance_ && i > nonsheetrescount+sheetsize_ ) keepsecondsheet = false;
	}
	if ( keepfirstsheet || keepsecondsheet ) {
		Size nsheets = (pose.total_residue()-total_residues_)/sheetsize_;
		for ( Size i=1; i<=nsheets; i++ ) {
			Size lower = total_residues_+1;
			Size upper = total_residues_+sheetsize_;
			if ( i > 1 ) {
				lower = total_residues_+sheetsize_+1;
				upper = pose.total_residue();
			}
			Real hbond = get_resrange_hbond_energy(pose,lower,upper);
			if ( hbond == 0 ) {
				if ( i == 1 ) keepfirstsheet = false;
				if ( i == 2 ) keepsecondsheet = false;
			}
		}
	}
	if ( !keepfirstsheet && !keepsecondsheet ) includesheets=0;
	if ( keepfirstsheet && !keepsecondsheet ) includesheets=1;
	if ( !keepfirstsheet && keepsecondsheet ) includesheets=2;
	if ( keepfirstsheet && keepsecondsheet ) includesheets=3;

	//Calculate modified fast dens score
	for ( Size i=1; i<=residueatomdensities.size(); i++ ) {
		int lower_window = i-res_windowlength;
		int upper_window = i+res_windowlength;
		bool sheetres = false;
		if ( i > nonsheetrescount ) sheetres = true;
		if ( lower_window < 1 ) lower_window = 1;
		if ( upper_window > (int)rangehi-(int)rangelo ) upper_window = rangehi-rangelo+1;
		worst_atom = -99999;
		bool addresforfast = true;
		//The 4 in the next line is the number of sheet residues added
		if ( i > nonsheetrescount && i <= nonsheetrescount+sheetsize_ && !keepfirstsheet ) addresforfast = false;
		if ( i > nonsheetrescount+sheetsize_ && !keepsecondsheet ) addresforfast = false;
		for ( int j=lower_window; j<=upper_window; j++ ) {
			if ( j > (int)nonsheetrescount ) continue;
			Real worst = worstatoms[j];
			if ( worst > worst_atom ) worst_atom = worst;
		}
		if ( !sheetres ) {
			Real basescore = 0;
			Real fullscore = 0;
			for ( Size ii=1; ii<=residueatomdensities[i].size(); ii++ ) {
				Real atomscore = residueatomdensities[i][ii];
				fullscore+=atomscore;
				if ( atomscore > worst_atom ) {
					basescore += atomscore;
				} else {
					basescore += worst_atom;
				}
			}
			//TRACER << " basescore fullscore " << basescore << " " << fullscore << " " << std::endl;
			core::Real residue_dens = basescore + continuous_weight_*(fullscore-basescore);
			mod_dens += residue_dens;
		} else if ( addresforfast ) {
			Real basescore = 0;
			Real fullscore = 0;
			for ( Size ii=1; ii<=residueatomdensities[i].size(); ii++ ) {
				Real atomscore = residueatomdensities[i][ii];
				fullscore+= atomscore;
				if ( atomscore > worst_atom ) {
					basescore += atomscore;
				} else {
					basescore += worst_atom;
				}
			}
			core::Real residue_dens = (basescore + continuous_weight_*(fullscore-basescore))*sheetbonus_;
			mod_dens += residue_dens;
		}
	}
	TRACER.Debug << " mod dens is " << mod_dens << " total fast dens " << total_fast_dens << " sc dens is " << sc_dens << std::endl;

	Real adjusted_dens = 0;
	if ( windowdensweight_ != 0 ) {
		Real totalresdens = 0;
		for ( core::Size i = rangelo; i <= pose.total_residue(); i++ ) {
			if ( i > rangehi && i <= total_residues_ ) continue;
			if ( includesheets == 0 && i > rangehi ) continue;
			if ( includesheets == 1 && i > rangehi && i > total_residues_+sheetsize_ ) continue;
			if ( includesheets == 2 && i > rangehi && i <= total_residues_+5 ) continue;
			core::conformation::Residue currentres = pose.residue(i);
			core::scoring::electron_density::getDensityMap().setSCscaling( 0.2 );
			Real residuedens = windowdensweight_ * -core::scoring::electron_density::getDensityMap().matchRes(i, currentres, pose, 0, false);
			totalresdens += residuedens;
		}
		//this deals with really bad cases they may get a negative score due to the way continuous density is calculated. This code just gives them a poor score rather than using those equations
		if ( total_fast_dens >=-1 ) {
			adjusted_dens = 0;
		} else {
			adjusted_dens = totalresdens*(std::min(mod_dens,-0.0001)/total_fast_dens);
		}
	} else {
		adjusted_dens = mod_dens;
	}
	return adjusted_dens;

}
//Only takes centroid pose
Real LoopGrower::sheetscore( core::pose::Pose& fa_pose, core::pose::Pose& cen_pose, Size rangelo, Size rangehi ){

	Size irrelevant_number;
	Real density_weight;
	if ( cenrot_ ) {
		density_weight = cenrot_sf_->get_weight(scoring::elec_dens_fast);
	} else {
		density_weight = cen_sf_->get_weight(scoring::elec_dens_fast);
	}
	Real original_score = modifieddensity( cen_pose, rangelo, rangehi, density_weight, irrelevant_number );
	if ( cen_pose.total_residue() > total_residues_ ) {
		cen_pose.conformation().delete_residue_range_slow(total_residues_+1, cen_pose.total_residue());
		fa_pose.conformation().delete_residue_range_slow(total_residues_+1, fa_pose.total_residue());
	}
	Real nosheets = modifieddensity( cen_pose, rangelo, rangehi, density_weight, irrelevant_number );
	Real strandscore = original_score - nosheets;
	if ( strandscore > 0 ) strandscore = 0;
	return strandscore;
}

Real LoopGrower::modifiedscore( core::pose::Pose& fapose, core::pose::Pose& cen_pose, Size rangelo, Size rangehi){

	core::scoring::ScoreFunctionOP sf_nodens;
	if ( pack_min_cycles_ != 0 && famin_ ) {
		sf_nodens = sf_->clone();
		//sf_nodens = cen_sf_->clone();
	} else if ( cenrot_ ) {
		sf_nodens = cenrot_sf_->clone();
	} else {
		sf_nodens = cen_sf_->clone();
	}
	sf_nodens->set_weight( core::scoring::elec_dens_fast , 0.0 );
	Real density_weight;
	Real basescore = startingscore_;
	if ( ( pack_min_cycles_ == 0 || !famin_ ) && !cenrot_ ) {
		basescore = censtartingscore_;
	}
	//if(!cenrot_ ) basescore = censtartingscore_;
	if ( pack_min_cycles_ !=0 ) {
		density_weight = sf_->get_weight(scoring::elec_dens_fast);
	} else {
		if ( cenrot_ ) {
			density_weight = cenrot_sf_->get_weight(scoring::elec_dens_fast);
		} else {
			density_weight = cen_sf_->get_weight(scoring::elec_dens_fast);
		}
	}
	core::pose::Pose pose;
	if ( pack_min_cycles_ != 0 || cenrot_ ) {
		pose = fapose;
	} else {
		pose = cen_pose;
	}

	Size bestsheets = 0;
	Real adjusted_dens = 0;
	if ( windowdensweight_ != 0 || density_weight != 0 ) {
		adjusted_dens = modifieddensity( pose, rangelo, rangehi, density_weight, bestsheets );
	}
	//TR << "no sheet dens " << adjusted_dens << std::endl;
	if ( bestsheets == 0 ) {
		cen_pose.conformation().delete_residue_range_slow(total_residues_+1, pose.total_residue());
		fapose.conformation().delete_residue_range_slow(total_residues_+1, pose.total_residue());
	}
	if ( bestsheets == 1 ) {
		cen_pose.conformation().delete_residue_range_slow(total_residues_+sheetsize_+1, pose.total_residue());
		fapose.conformation().delete_residue_range_slow(total_residues_+sheetsize_+1, pose.total_residue());
	}
	if ( bestsheets == 2 ) {
		cen_pose.conformation().delete_residue_range_slow(total_residues_+1, total_residues_+sheetsize_);
		fapose.conformation().delete_residue_range_slow(total_residues_+1, total_residues_+sheetsize_);
	}

	Real score;
	if ( pack_min_cycles_ !=0 ) {
		if ( famin_ ) {
			score = (*sf_nodens)(fapose) + adjusted_dens;
		} else {
			Real nton3hbonds = nton3_hbond_score(cen_pose);
			score = (*sf_nodens)(cen_pose) + adjusted_dens;
			Real hbond_weight = cen_sf_->get_weight(scoring::hbond_sr_bb);
			score -= nton3hbonds*hbond_weight*.7;
		}

	} else {
		if ( cenrot_ ) {
			Real nton3hbonds = nton3_hbond_score(fapose);
			score = (*sf_nodens)(fapose) + adjusted_dens;
			Real hbond_weight = cenrot_sf_->get_weight(scoring::hbond_sr_bb);
			score -= nton3hbonds*hbond_weight*.7;
		} else {
			Real nton3hbonds = nton3_hbond_score(cen_pose);
			score = (*sf_nodens)(cen_pose) + adjusted_dens;
			Real hbond_weight = cen_sf_->get_weight(scoring::hbond_sr_bb);
			score -= nton3hbonds*hbond_weight*.7;
		}
	}
	TRACER.Debug << "pre subtraction score is " << score << " basescore is " << basescore << std::endl;
	score -= basescore;
	if ( !rescorebeams_ ) TRACER << " score is " << score << std::endl;
	return score;
}
void
LoopGrower::add_fragment_csts( core::pose::Pose &pose, Size startfasta, Size endfasta, Size natstoplow, Size natstarthi, Size startlower ) {

	core::fragment::FragSetOP frags = fragments_[1];

	// stolen from SecondaryStructure.cc
	if ( frags->global_offset() != 0 ) {
		TR.Error << "[ERROR] SecondaryStructure computations must be carried out with local coordinates (global offset of fragments must be 0)." << std::endl;
		runtime_assert( false );
	}

	// 1 - collect fragment statistics
	Size frag_nres = frags->max_pos();
	utility::vector1< utility::vector1< core::Real > > phi_distr( frag_nres, utility::vector1< core::Real >(36,0.0) );
	utility::vector1< utility::vector1< core::Real > > psi_distr( frag_nres, utility::vector1< core::Real >(36,0.0) );
	utility::vector1< int > N( frag_nres, 0 );

	for ( core::fragment::FragID_Iterator it=frags->begin(), eit=frags->end(); it!=eit; ++it ) { //carefully checked that I don't change FrameData
		core::Size loop_start = 1;
		core::Size loop_end = it->frame().length();
		for ( core::Size fpos = loop_start; fpos <= loop_end; ++fpos ) {
			core::fragment::BBTorsionSRFDCOP res_i =
				utility::pointer::dynamic_pointer_cast<const core::fragment::BBTorsionSRFD> (it->fragment().get_residue( fpos ) );

			Size pos = it->frame().seqpos( fpos );
			core::Real phi = std::fmod( res_i->torsion(1), 360.0 ); if ( phi<0 ) phi +=360.0;
			core::Real psi = std::fmod( res_i->torsion(2), 360.0 ); if ( psi<0 ) psi +=360.0;

			core::Size phibin = (core::Size) std::floor( phi/10.0 ); if ( phibin == 36 ) phibin=0;
			core::Size psibin = (core::Size) std::floor( psi/10.0 ); if ( psibin == 36 ) psibin=0;

			phi_distr[pos][phibin+1]+=1.0;
			psi_distr[pos][psibin+1]+=1.0;
			N[pos]++;
		}
	}

	// smoothing
	//core::Real smoothing[7] = {0.09, 0.14, 0.175, 0.19, 0.175, 0.14, 0.09};  // very soft
	core::Real smoothing[7] = {0.01, 0.05, 0.24, 0.40, 0.24, 0.05, 0.01}; // less soft
	for ( int i=1; i<=(int)frag_nres; ++i ) {
		utility::vector1< core::Real > phi_i = phi_distr[i];
		utility::vector1< core::Real > psi_i = psi_distr[i];
		for ( int j=1; j<=36; ++j ) {
			phi_distr[i][j] = smoothing[0]*phi_i[1+((j+32)%36)]
				+ smoothing[1]*phi_i[1+((j+33)%36)]
				+ smoothing[2]*phi_i[1+((j+34)%36)]
				+ smoothing[3]*phi_i[1+((j+35)%36)]
				+ smoothing[4]*phi_i[1+((j+0)%36)]
				+ smoothing[5]*phi_i[1+((j+1)%36)]
				+ smoothing[6]*phi_i[1+((j+2)%36)];

			psi_distr[i][j] = smoothing[0]*psi_i[1+((j+32)%36)]
				+ smoothing[1]*psi_i[1+((j+33)%36)]
				+ smoothing[2]*psi_i[1+((j+34)%36)]
				+ smoothing[3]*psi_i[1+((j+35)%36)]
				+ smoothing[4]*psi_i[1+((j+0)%36)]
				+ smoothing[5]*psi_i[1+((j+1)%36)]
				+ smoothing[6]*psi_i[1+((j+2)%36)];

			phi_distr[i][j] /= N[i];
			psi_distr[i][j] /= N[i];
		}
	}

	// convert to energy
	core::Real mest=0.1;
	for ( int i=1; i<=(int)frag_nres; ++i ) {
		for ( int j=1; j<=36; ++j ) {
			// shift so max is 0
			phi_distr[i][j] = -std::log( mest/36.0 + (1-mest)*phi_distr[i][j] ) + std::log( mest/36.0 ) ;
			psi_distr[i][j] = -std::log( mest/36.0 + (1-mest)*psi_distr[i][j] ) + std::log( mest/36.0 ) ;
		}

	}

	// finally .. add constraints
	Size pos = startlower;
	for ( int i=(int)startfasta; i<=(int)endfasta; ++i ) {
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;

		if ( i > (int)natstoplow && i < (int)natstarthi ) continue;
		char fastaseq = seq_->at(i);
		char poseseq = pose.residue(pos).name1();
		if ( fastaseq != poseseq ) {
			TRACER << " the fragment constraint you are trying to apply does not match the correct sequence" << std::endl;
			runtime_assert( fastaseq == poseseq );
		}

		if ( !pose.residue(pos).is_protein() ) continue;

		if ( pos>1 && pose.residue(pos-1).is_protein() ) {
			FuncOP phi_func( new CircularSplineFunc( 1.0, phi_distr[i] ) );
			ConstraintOP phi_cst( new DihedralConstraint(
				core::id::AtomID(3,pos-1),core::id::AtomID(1,pos),core::id::AtomID(2,pos),core::id::AtomID(3,pos), phi_func  ) );
			pose.add_constraint( scoring::constraints::ConstraintCOP( phi_cst ) );
		}

		if ( pos<frag_nres && pose.residue(pos+1).is_protein() ) {
			FuncOP psi_func( new CircularSplineFunc( 1.0, psi_distr[i] ) );
			ConstraintOP psi_cst( new DihedralConstraint(
				core::id::AtomID(1,pos),core::id::AtomID(2,pos),core::id::AtomID(3,pos),core::id::AtomID(1,pos+1), psi_func  ) );
			pose.add_constraint( scoring::constraints::ConstraintCOP( psi_cst ) );
		}
		pos++;
	}
}
void
LoopGrower::add_user_csts(core::pose::Pose & pose){
	if ( (basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user()) || (basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_file ].user()) ) {
		//map sequence
		TRACER << " applying csts " << std::endl;
		core::sequence::SequenceOP t_pdb_seq( new core::sequence::Sequence( pose.sequence(), "pose_seq" ));
		core::sequence::SWAligner sw_align;
		core::sequence::ScoringSchemeOP ss(  new core::sequence::SimpleScoringScheme(120, 0, -100, 0));
		core::sequence::SequenceAlignment fasta2template_;
		fasta2template_ = sw_align.align(seq_, t_pdb_seq, ss);
		core::id::SequenceMapping sequencemap = fasta2template_.sequence_mapping(1,2);
		core::id::SequenceMappingOP seqmapOP( new core::id::SequenceMapping( sequencemap ));
		core::scoring::constraints::ConstraintSetOP constraints( new core::scoring::constraints::ConstraintSet() );
		if ( pose.is_centroid() ) {
			core::scoring::constraints::add_constraints_from_cmdline_to_pose( cen_seqpose_ );
			constraints = cen_seqpose_.constraint_set()->remapped_clone( cen_seqpose_, pose, seqmapOP );
			pose.constraint_set( constraints );
		} else {
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( fa_seqpose_ );
			constraints = fa_seqpose_.constraint_set()->remapped_clone( fa_seqpose_, pose, seqmapOP );
			pose.constraint_set( constraints );
		}
	}
}


Real LoopGrower::single_grow( core::pose::Pose& growpose, core::pose::Pose& growpose_cen, LoopPartialSolutionStore& solutionset, const core::chemical::ResidueTypeCOPs& restypes_pose,
	const core::chemical::ResidueTypeCOPs& restypes_pose_cen, Size lower_pose, Size upper_pose, Size upper_term, int lower_fasta, int upper_fasta, Size torsionrangelo, Size torsionrangehi,
	bool insert_lower, Size initial_melt_left, bool is_cterm, bool is_nterm, Size cycle, int n_to_insert ){

	int insert_pose=0, insert_fasta=0;
	Size maxfrag = fragments_[1]->max_frag_length();
	for ( Size i=1; i<=fragments_.size(); i++ ) {
		core::fragment::FragSetOP fragset = fragments_[i];
		if ( maxfrag > fragset->max_frag_length() ) {
			maxfrag = fragset->max_frag_length();
		}
	}
	maxfrag_ = maxfrag;
	Real bestdensity = 99;
	Real bestrms = 9999;
	Real bestgdt = 0;
	Size is_lower;
	if ( insert_lower ) {
		is_lower = 1;
	} else {
		is_lower = 2;
	}

	LoopPartialSolutionStore solutionsetnew( beamwidth_, rmscutoff_, master_beam_width_, master_beam_cutoff_ );
	solutionsetnew.setfilterparams( fragmelt_, rmswindow_, parallelcount_, beamscorecutoff_, dumperrors_, dumpbeam_, writebeams_, clustercheck_, fafilter_, samplesheets_, filterprevious_,
		checksymm_, asymmdump_, dumpfinalbeam_ );
	if ( filterbeams_ != "" && filterprevious_ ) {
		utility::vector1<LoopPartialSolution> filtersolutions = solutionset.get_filteronly_solutions();
		solutionsetnew.store_filteronly_solutions(filtersolutions);
	}


	// copy side chains
	//if ( restoresc ) restore_sc->apply( growpose );


	// special logic if we are done
	if ( upper_fasta-lower_fasta+1 < n_to_insert ) {
		n_to_insert = upper_fasta-lower_fasta+1;
	}


	if ( insert_lower ) {
		for ( int i=1; i<=(int)n_to_insert; ++i ) {
			core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[lower_fasta+i-1] ));
			growpose.conformation().safely_append_polymer_residue_after_seqpos(*newres, lower_pose+i-1, true);
			core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[lower_fasta+i-1] ));
			growpose_cen.conformation().safely_append_polymer_residue_after_seqpos(*newres_cen, lower_pose+i-1, true);
		}

		insert_pose = lower_pose+n_to_insert-((int)maxfrag-1);
		TRACER << "insert " << n_to_insert << " lower at " << insert_pose << std::endl;
		insert_fasta = lower_fasta - initial_melt_left + ( n_to_insert -((int)maxfrag-(int)initial_melt_left));
		initial_melt_left = fragmelt_;
		TRACER << "lower_fasta = " << lower_fasta << std::endl;
		lower_pose += n_to_insert;
		lower_fasta += n_to_insert;
		upper_pose += n_to_insert;
	} else {
		for ( int i=1; i<=(int)n_to_insert; ++i ) {
			core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(*restypes_pose[upper_fasta-i+1] ));
			growpose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, lower_pose+1, true);
			core::conformation::ResidueOP newres_cen( core::conformation::ResidueFactory::create_residue(*restypes_pose_cen[upper_fasta-i+1] ));
			growpose_cen.conformation().safely_prepend_polymer_residue_before_seqpos(*newres_cen, lower_pose+1, true);
		}
		insert_pose = upper_term;
		insert_fasta = upper_fasta-n_to_insert+1;
		TRACER << "lower pose & upper pose " << lower_pose << " " << upper_pose << std::endl;
		upper_pose += n_to_insert;
		upper_fasta -= n_to_insert;
		TRACER << "insert " << n_to_insert << " upper at " << insert_pose << std::endl;
	}

	// loop is complete
	if ( upper_fasta<lower_fasta ) {
		// add chainbreak vars, turn on chainbreak score
		if ( !is_nterm && !is_cterm ) {
			core::conformation::remove_upper_terminus_type_from_conformation_residue( growpose.conformation(), lower_pose );
			core::conformation::remove_lower_terminus_type_from_conformation_residue( growpose.conformation(), lower_pose+1 );
			core::conformation::remove_upper_terminus_type_from_conformation_residue( growpose_cen.conformation(), lower_pose );
			core::conformation::remove_lower_terminus_type_from_conformation_residue( growpose_cen.conformation(), lower_pose+1 );

			growpose.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");
			growpose_cen.conformation().declare_chemical_bond(lower_pose, "C", lower_pose+1, "N");

			core::pose::add_variant_type_to_pose_residue( growpose, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::add_variant_type_to_pose_residue( growpose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			core::pose::add_variant_type_to_pose_residue( growpose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::add_variant_type_to_pose_residue( growpose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );

			sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
			cen_sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
			cenrot_sf_->set_weight( core::scoring::chainbreak, chainbreak_ );
		} else if ( is_nterm ) {
			core::conformation::add_lower_terminus_type_to_conformation_residue( growpose.conformation(), lower_pose+1 );
			core::conformation::add_lower_terminus_type_to_conformation_residue( growpose_cen.conformation(), lower_pose+1 );
		} else if ( is_cterm ) {
			core::conformation::add_upper_terminus_type_to_conformation_residue( growpose.conformation(), lower_pose );
			core::conformation::add_upper_terminus_type_to_conformation_residue( growpose_cen.conformation(), lower_pose );
		}
	}

	//store side chains.
	protocols::moves::MoverOP restore_sc;
	if ( pack_min_cycles_ != 0 || cenrot_ ) {
		restore_sc = protocols::simple_moves::ReturnSidechainMoverOP( new protocols::simple_moves::ReturnSidechainMover( growpose ));
	}

	core::Real beam_maximum = -99999;
	core::Real beam_minimum = 99999;

	// expand each previous rounds' hit
	Size ctr=1;
	torsionrangehi += n_to_insert;
	if ( is_cterm ) torsionrangehi  = upper_pose;
	TRACER << " the current solution set size is " << solutionset.size() << std::endl;
	core::Size total_lower;
	if ( !is_nterm ) {
		total_lower = lower_pose-torsionrangelo+1;
	} else {
		total_lower = 0;
	}
	Size omit_pose = insert_pose;
	insert_pose_ = insert_pose;

	Size sheetrangelo;
	Size sheetrangehi;
	if ( insert_lower ) {
		sheetrangelo = lower_pose-sheetsize_+1;
		sheetrangehi = lower_pose;
	} else {
		sheetrangelo = lower_pose+1;
		sheetrangehi = sheetrangelo+sheetsize_-1;
	}
	//disable sheet sampling if there aren't enough base residues (this shouldn't actually ever happen).
	if ( sheetrangelo < torsionrangelo || sheetrangehi > torsionrangehi ) samplesheets_ = false;
	TR << " sheets from " << sheetrangelo << " to " << sheetrangehi << std::endl;
	SheetSampler sheetsampler(sheetrangelo, sheetrangehi, 200);

	// add distance constraints if the loops are nearly closed"
	if ( !is_nterm && !is_cterm && (upper_fasta-lower_fasta+1 <= 10) ) {
		core::Real len = gapdist(upper_fasta-lower_fasta+1);
		//if symmetry is unclear check all transforms
		if ( checksymm_ ) {
			transform_to_closest_symmunit(growpose_cen,growpose,lower_pose);
		}


		core::scoring::constraints::BoundFuncOP myfunc(new core::scoring::constraints::BoundFunc( 0.0, len, 1.0, "x"));
		core::scoring::constraints::AtomPairConstraintOP close_cst(
			new core::scoring::constraints::AtomPairConstraint(
			core::id::AtomID(growpose.residue(lower_pose).atom_index("C"),lower_pose),
			core::id::AtomID(growpose.residue(lower_pose+1).atom_index("N"),lower_pose+1),
			myfunc
			));

		growpose.add_constraint( close_cst );
		growpose_cen.add_constraint( close_cst );
		Real constraint_weight;
		if ( pack_min_cycles_ != 0 ) {
			constraint_weight = sf_->get_weight(core::scoring::atom_pair_constraint);
		} else if ( cenrot_ ) {
			constraint_weight = cenrot_sf_->get_weight(core::scoring::atom_pair_constraint);
		} else {
			constraint_weight = cen_sf_->get_weight(core::scoring::atom_pair_constraint);
		}
		if ( constraint_weight  ==  0 ) {
			TRACER << "WARNING::Constraints on closure are turned off when you need them!" << std::endl;
			//cen_sf_->show(std::cout, growpose_cen);
		}
	}

	//applies user constraints only if they are set in the command line
	add_user_csts(growpose);
	add_user_csts(growpose_cen);

	for ( Size j=1; j<=solutionset.size(); j++ ) {
		LoopPartialSolution lps = solutionset[j];
		utility::vector1< std::pair<Size,Size> > fraghistory;
		if ( trackfragments_ ) {
			fraghistory = lps.get_history();
		}
		Real bonus_score = lps.get_bonus_score();
		int omit_low = insert_lower ? omit_pose:omit_pose-1;
		int omit_hi = insert_lower ? omit_pose+n_to_insert+1:omit_pose+n_to_insert;
		lps.apply (growpose_cen, torsionrangelo, omit_low, omit_hi, torsionrangehi);  // pose will get updated from pose_cen

		core::Real beam_score = lps.get_score();
		if ( beam_score > beam_maximum ) {
			beam_maximum = beam_score;
		}
		if ( beam_score < beam_minimum ) {
			beam_minimum = beam_score;
		}

		// random subset of fragments
		core::pose::Pose scpose = growpose;
		core::pose::Pose sheetsamplepose = growpose_cen;
		Size fragtracker = 1;
		if ( nativegrow_ ) fragtracker = 0;
		for ( Size fraglistitr=1; fraglistitr<=fragments_.size(); fraglistitr++ ) {
			core::fragment::FragSetOP fragset = fragments_[fraglistitr];
			int fraglength = fragset->max_frag_length();
			if ( insert_lower ) {
				insert_fasta = lower_fasta - fraglength;
				insert_pose = lower_pose - fraglength + 1;
			} else {
				insert_fasta = upper_fasta - fraglength + maxfrag +1;
			}
			if ( insert_fasta < 1 ) {
				continue;
			}

			//TRACER << "insert pose " << insert_pose << " torsion range lo and hi " << torsionrangelo << " " << torsionrangehi << " fraglength " << fraglength << std::endl;
			if ( (int)insert_pose < (int)torsionrangelo || (int)insert_pose + fraglength - 1 > (int)torsionrangehi ) {
				//TRACER << "iskipping " << fraglength << "mers for being outside the torsionranges " << std::endl;
				continue;
			}
			core::fragment::Frame frame=library_[fraglistitr][insert_fasta];
			frame.shift_to(insert_pose);
			utility::vector1< int > allfrags(frame.nr_frags());
			for ( int i=1; i<=(int)frame.nr_frags(); ++i ) allfrags[i] = i;
			if ( !trackfragments_ ) std::random_shuffle( allfrags.begin(), allfrags.end() );
			int ntrials = std::min( (int)frame.nr_frags(), (int)fragtrials_ );
			bool finished = false;

			//we need the original pose to reset sheets and sidechains for cenrot and fa and the sheets for centroid
			total_residues_ = growpose_cen.total_residue();
			for ( int i = nativegrow_ ? 0:1; i<=ntrials && !finished; ) {
				//reseting torsions
				lps.apply (growpose_cen, torsionrangelo, omit_low, omit_hi, torsionrangehi);
				std::string beamid = utility::to_string(cycle)+"."+utility::to_string(ctr)+"."+utility::to_string(i)+"."+utility::to_string(is_lower);
				// copy side chains
				// this logic is redundant now that we need to use it for sheet sampling fix it later
				if ( cenrot_ ) {
					growpose = scpose;
				}
				if ( pack_min_cycles_ != 0 ) {
					growpose = scpose;
				}
				if ( samplesheets_ ) {
					growpose_cen = sheetsamplepose;
				}

				//Steal Native Torsions for testing
				if ( i==0 ) {
					Real natscore = 0;

					//Set Initial Native Loop
					LoopPartialSolution nativeloop( *native_, insert_fasta, insert_fasta+maxfrag-1, natscore) ; // added for native
					nativeloop.apply (growpose_cen, insert_pose, insert_pose+maxfrag-1);

					if ( fragtrials_ == 0 ) {
						TRACER << "frag trials 0 " << std::endl;
						finished = true;
					}
					if ( debug_ ) growpose_cen.dump_pdb("afternative_beforesheets_"+utility::to_string( cycle )+"_"+utility::to_string( ctr )+"_"+utility::to_string( is_lower )+".pdb" );

					//If you are stealing torsions from the native and have the sheet sampler turned on it will always build them for now
					if ( samplesheets_ ) {
						sheetsampler.apply(growpose_cen);
					}

					if ( debug_ ) growpose_cen.dump_pdb("afternative_aftersheets_"+utility::to_string( cycle )+"_"+utility::to_string( ctr )+"_"+utility::to_string( is_lower )+".pdb" );
				} else {
					int idx = allfrags[i];
					frame.apply( idx, growpose_cen );  // pose will get updated from pose_cen

					if ( samplesheets_ ) {
						core::fragment::FragData tempfrag = frame.fragment(i);
						std::string secstructure = tempfrag.secstruct();
						bool meets_criteria = false;
						if ( sheetcriteria_ == 1 ) {
							if ( boost::algorithm::ends_with(secstructure, "EEE") ) {
								meets_criteria = true;
							}
						}
						if ( sheetcriteria_ == 2 ) {
							Size res1 = sheetrangelo + (sheetrangehi-sheetrangelo)/2;
							Size res2 = res1-1;
							if ( is_beta(growpose_cen.phi(res1),growpose_cen.psi(res1)) || is_beta(growpose_cen.phi(res2), growpose_cen.psi(res2)) ) {
								meets_criteria = true;
							}
						}
						if ( sheetcriteria_ == 3 ) meets_criteria = true;
						TR << " sheet criteria is " << sheetcriteria_ << " meets criteria " << meets_criteria << std::endl;
						if ( meets_criteria ) sheetsampler.apply(growpose_cen);
					}
				}
				if ( !is_nterm && !is_cterm && (upper_fasta-lower_fasta+1 <= 10) ) {
					//if symmetry is unclear check all transforms
					if ( checksymm_ ) {
						transform_to_closest_symmunit(growpose_cen,growpose,lower_pose);
					}
				}

				//calculate RMS of loop to native.
				if ( native_ ) {
					Real RMS = RMStonative(growpose_cen, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose-fragmelt_+1, lower_pose+fragmelt_-1, torsionrangehi);
					//TRACER << "The RMS before refinement ===== " << RMS << std::endl;
					if ( debug_ ) growpose_cen.dump_pdb("afterfrag_"+utility::to_string(RMS)+"_"+utility::to_string( cycle )+"."+utility::to_string( ctr )+"."+utility::to_string( i )+"_"+utility::to_string( is_lower )+".pdb" );
				}

				//output here also
				if ( debug_ && !native_ ) growpose_cen.dump_pdb("afterfrag_"+utility::to_string( cycle )+"."+utility::to_string( ctr )+"."+utility::to_string( i )+"_"+utility::to_string( is_lower )+".pdb" );



				//this tells refine cycle where to place the sheet residues for refinement
				int refine_lower = lower_pose-minmelt_;
				int refine_upper = lower_pose+1+minmelt_;
				if ( refine_lower < (int)torsionrangelo ) refine_lower = torsionrangelo;
				if ( refine_upper > (int)torsionrangehi ) refine_upper = torsionrangehi;
				refine_cycle( growpose, growpose_cen, torsionrangelo, torsionrangehi, false, cycle, ctr, i, is_lower);


				core::scoring::ScoreFunctionOP sf_densonly;
				if ( core::pose::symmetry::is_symmetric( growpose ) ) {
					//sf_densonly = core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction );
					core::scoring::symmetry::SymmetricScoreFunctionOP symmdens( new core::scoring::symmetry::SymmetricScoreFunction );
					sf_densonly = utility::pointer::dynamic_pointer_cast< core::scoring::symmetry::SymmetricScoreFunction >(symmdens);
				} else {
					sf_densonly = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
				}
				sf_densonly->set_weight( core::scoring::elec_dens_fast , sf_->get_weight(core::scoring::elec_dens_fast) );
				core::scoring::ScoreFunctionOP sf_nodens;
				sf_nodens = sf_->clone();
				sf_nodens->set_weight( core::scoring::elec_dens_fast , 0.0 );

				//Store Solutions

				Real originalscore;
				Real modscore;
				if ( pack_min_cycles_ != 0 ) {
					originalscore = (*sf_)(growpose);
					modscore = modifiedscore(growpose, growpose_cen, torsionrangelo, torsionrangehi);
				} else {
					if ( cenrot_ ) {
						originalscore = (*cenrot_sf_)(growpose);
					} else {
						originalscore = (*cen_sf_)(growpose_cen);
					}
					modscore = modifiedscore(growpose, growpose_cen, torsionrangelo, torsionrangehi);
				}
				modscore += bonus_score;

				//report and store beam
				Real dens_only;
				if ( pack_min_cycles_ !=0 ) {
					dens_only = (*sf_densonly)(growpose);
				} else {
					dens_only = (*sf_densonly)(growpose_cen);
				}
				LoopPartialSolution newlps = LoopPartialSolution(growpose, torsionrangelo, torsionrangehi, modscore);
				std::pair<Size,Size> historypair = std::make_pair( ctr, fragtracker );
				if ( trackfragments_ ) newlps.set_history( fraghistory, historypair );

				//calculate RMS of loop to native.
				if ( native_ ) {
					Real RMS;
					Real GDT;
					if ( pack_min_cycles_ !=0 || fafilter_ || cenrot_ ) {
						RMS = RMStonative(growpose, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
						GDT = GDThatonative(growpose, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
					} else {
						GDT = GDThatonative(growpose_cen, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
						RMS = RMStonative(growpose_cen, rmsrangelo_, rmsrangehi_, lower_fasta-1, upper_fasta+1, torsionrangelo, lower_pose, lower_pose, torsionrangehi);
					}
					if ( GDT > bestgdt ) bestgdt = GDT;
					if ( RMS < bestrms ) {
						bestrms = RMS;
					}
					TRACER << "The RMS after refinement ===== " << RMS << " GDT == " << GDT << std::endl;
					TRACER << "The best rms for this grow cycle = " << bestrms << " best GDT = " << bestgdt << std::endl;
					newlps.set_rms(RMS);
					newlps.set_gdt(GDT);
					std::string id = utility::to_string(cycle)+"."+ utility::to_string(ctr)+"."+utility::to_string(fragtracker);
					newlps.set_id(id);
				}
				if ( samplesheets_ ) {
					store_sheets(growpose_cen, newlps);
					Real strand_score = sheetscore(growpose, growpose_cen, torsionrangelo, torsionrangehi);

					//if sheets werent added or aren't good we reset the score bonus back to 0
					if ( strand_score == 0 || !continuous_sheets_ ) {
						bonus_score = 0;
					}

					if ( fafilter_ ) {
						newlps.set_bonus_score(bonus_score);
					} else {
						Real kept_bonus = strand_score+bonus_score;
						if ( !continuous_sheets_ ) kept_bonus = 0;
						newlps.set_bonus_score(kept_bonus);
					}
				}

				solutionsetnew.store(newlps);

				TRACER << "beamscoring " << cycle << "." << ctr << "." << i << "." << is_lower << " " << modscore << " dens only " << dens_only << " original score " << originalscore << std::endl;

				if ( debug_ ) {
					if ( pack_min_cycles_ !=0 ) growpose.dump_pdb("fullatomafterrefine_"+utility::to_string( cycle )+"_"+utility::to_string( ctr )+"_"+utility::to_string( i )+"_"+utility::to_string( is_lower )+".pdb" );
					if ( samplesheets_ !=0 ) growpose_cen.dump_pdb("aftersheetrefine_"+utility::to_string( cycle )+"_"+utility::to_string( ctr )+"_"+utility::to_string( i )+"_"+utility::to_string( is_lower )+".pdb" );
				}


				i++;
				fragtracker++;

				//best score for the direction
				if ( dens_only < bestdensity || bestdensity == 99 ) {
					bestdensity = dens_only;
				}
				//reset to the starting pose
				growpose = scpose;
				growpose_cen = sheetsamplepose;
			}
			ctr++;
		}
	}


	if ( coordfile_ != "" ) {
		coordinate_filter(solutionsetnew, growpose_cen, insert_lower, lower_fasta, upper_fasta, torsionrangelo, torsionrangehi);
	}

	if ( pack_min_cycles_ != 0 ) restore_sc->apply( growpose );
	//Size totalres = torsionrangehi-torsionrangelo;
	//Real beamscorecutoff = totalres * beamscorecutoff_ +8;
	if ( fafilter_ || cenrotfilter_ ) {
		fafilter(solutionsetnew, growpose, growpose_cen, total_lower, torsionrangelo, torsionrangehi, cycle, lower_fasta, upper_fasta, lower_pose);
	} else if ( pack_min_cycles_ != 0 ) {
		if ( parametercheck_ ) {
			solutionsetnew.parametercheck_filter(growpose, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
				lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
		} else {
			solutionsetnew.filter(growpose, beamwidth_, total_lower, torsionrangelo, torsionrangehi, cycle );
		}
	} else {
		if ( parametercheck_ ) {
			solutionsetnew.parametercheck_filter(growpose_cen, fragmelt_, total_lower, torsionrangelo, torsionrangehi, rmswindow_, dumpbeam_,
				lowest_ranked_native_, master_beam_width_, beamscorecutoff_, writebeams_, cycle, parallelcount_ );
		} else {
			solutionsetnew.filter(growpose_cen, beamwidth_, total_lower, torsionrangelo, torsionrangehi, cycle );
		}
	}
	if ( trackfragments_ ) {
		for ( Size i=1; i<=solutionsetnew.size(); i++ ) {
			LoopPartialSolution lps = solutionsetnew[i];
			utility::vector1< std::pair<Size,Size> > newhistory = lps.get_history();
			TRACER << " The fragments required to make " << lps.get_rms() << " were: ";
			for ( Size j=1; j<=newhistory.size(); j++ ) {
				TRACER << newhistory[j].first << "." << newhistory[j].second << " ";

			}
			TRACER << std::endl;
		}
	}
	// loop is complete
	if ( upper_fasta<=lower_fasta ) {
		if ( !is_nterm && !is_cterm ) {
			core::pose::remove_variant_type_from_pose_residue( growpose, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( growpose, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			core::pose::remove_variant_type_from_pose_residue( growpose_cen, core::chemical::CUTPOINT_LOWER, lower_pose   );
			core::pose::remove_variant_type_from_pose_residue( growpose_cen, core::chemical::CUTPOINT_UPPER, lower_pose+1 );
			sf_->set_weight( core::scoring::chainbreak, 0.0 );
			cen_sf_->set_weight( core::scoring::chainbreak, 0.0 );
			cenrot_sf_->set_weight( core::scoring::chainbreak, 0.0 );
		}
	}
	//core::Size num_missing = upper_fasta-lower_fasta+1;
	//if( skeleton_.has_density() && !is_nterm && !is_cterm && (num_missing <= 10) && (num_missing > 0) ) {
	//  solutionsetnew.skeleton_filter( growpose_cen, skeleton_, torsionrangelo, torsionrangehi, lower_pose, num_missing);
	//}
	solutionset = solutionsetnew;
	//TRACER << "the range of starting beams for cycle " << utility::to_string(cycle) << " was " << beam_maximum << " " << beam_minimum << std::endl;
	growpose.remove_constraints();
	growpose_cen.remove_constraints();
	return bestdensity;
}

}
}
