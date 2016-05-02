// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/components/VarLengthBuild.cc
/// @brief  Component that performs a simplified version of a protocol for
///         variable length remodeling of protein backbone segments.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/fragment_functions.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/symmetry/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/loops/LoopMoverFactory.hh>

// utility headers
#include <utility/vector1.hh>

#include <core/fragment/FrameIterator.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Auto Headers
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif


namespace protocols {
namespace forge {
namespace components {


// Tracer instance for this file
// Named after the original location of this code
static THREAD_LOCAL basic::Tracer TR( "protocols.forge.components.VarLengthBuild" );


/// @brief default constructor
VarLengthBuild::VarLengthBuild() :
	Super( "VLB" ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	recover_original_on_failure_( true ),
	loop_mover_str_( "RemodelLoopMover" ),
	cache_fragments_( true ),
	vall_memory_usage_( VLB_VallMemoryUsage::KEEP_IN_MEMORY ),
	num_fragpick_( 200 ),
	use_fullmer_( false ),
	max_linear_chainbreak_( 0.07 ),
	fragments_picked_( false ),
	user_provided_movers_apply_cycle_(3),
	loop_mover_fold_tree_constant_(false),
	restart_mode_( false ),
	ignore_cmdline_enzdes_cstfile_(false)
{}


/// @brief BuildManager constructor
VarLengthBuild::VarLengthBuild( BuildManager const & manager ) :
	Super( "VLB" ),
	manager_( manager ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	recover_original_on_failure_( true ),
	loop_mover_str_( "RemodelLoopMover" ),
	cache_fragments_( true ),
	vall_memory_usage_( VLB_VallMemoryUsage::KEEP_IN_MEMORY ),
	num_fragpick_( 200 ),
	use_fullmer_( false ),
	max_linear_chainbreak_( 0.07 ),
	fragments_picked_( false ),
	user_provided_movers_apply_cycle_(3),
	loop_mover_fold_tree_constant_(false),
	restart_mode_( false ),
	ignore_cmdline_enzdes_cstfile_(false)
{}

/// @brief BuildManager + remodelData constructor
VarLengthBuild::VarLengthBuild( BuildManager const & manager, RemodelData const & remodel_data ) :
	Super( "VLB" ),
	manager_( manager ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	remodel_data_( remodel_data ),
	recover_original_on_failure_( true ),
	loop_mover_str_( "RemodelLoopMover" ),
	cache_fragments_( true ),
	vall_memory_usage_( VLB_VallMemoryUsage::KEEP_IN_MEMORY ),
	num_fragpick_( 200 ),
	use_fullmer_( false ),
	max_linear_chainbreak_( 0.07 ),
	fragments_picked_( false ),
	user_provided_movers_apply_cycle_(3),
	loop_mover_fold_tree_constant_(false),
	restart_mode_( false ),
	ignore_cmdline_enzdes_cstfile_(false)
{}


/// @brief copy constructor
VarLengthBuild::VarLengthBuild( VarLengthBuild const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	manager_( rval.manager_ ),
	sfx_( rval.sfx_->clone() ),
	recover_original_on_failure_( rval.recover_original_on_failure_ ),
	loop_mover_str_( rval.loop_mover_str_ ),
	cache_fragments_( rval.cache_fragments_ ),
	vall_memory_usage_( rval.vall_memory_usage_ ),
	num_fragpick_( rval.num_fragpick_ ),
	use_fullmer_( rval.use_fullmer_ ),
	original_sequence_( rval.original_sequence_ ),
	new_secondary_structure_override_( rval.new_secondary_structure_override_ ),
	new_sequence_override_( rval.new_sequence_override_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	abego_( rval.abego_ ),
	rcgs_( rval.rcgs_ ),
	setup_movers_(rval.setup_movers_),
	user_provided_movers_(rval.user_provided_movers_),
	user_provided_movers_apply_cycle_(rval.user_provided_movers_apply_cycle_),
	loop_mover_fold_tree_constant_(rval.loop_mover_fold_tree_constant_),
	restart_mode_( rval.restart_mode_ ),
	ignore_cmdline_enzdes_cstfile_(rval.ignore_cmdline_enzdes_cstfile_)
{
	// assign fragments only if caching is active
	if ( cache_fragments_ ) {
		fragments_picked_ = rval.fragments_picked_;
		fragfull_ = rval.fragfull_;
		frag9_ = rval.frag9_;
		frag3_ = rval.frag3_;
		frag1_ = rval.frag1_;
	} else {
		fragments_picked_ = false;
	}
}


/// @brief default destructor
VarLengthBuild::~VarLengthBuild() {}


/// @brief clone this object
protocols::moves::MoverOP VarLengthBuild::clone() const {
	return protocols::moves::MoverOP( new VarLengthBuild( *this ) );
}


/// @brief create a new instance of this type of object
protocols::moves::MoverOP VarLengthBuild::fresh_instance() const {
	return protocols::moves::MoverOP( new VarLengthBuild() );
}


/// @brief set ScoreFunction used during build
void VarLengthBuild::scorefunction( ScoreFunction const & sfx ) {
	sfx_ = sfx.clone();
}


/// @brief set ScoreFunction used during build
void VarLengthBuild::scorefunction( ScoreFunctionOP const & sfx ) {
	sfx_ = sfx->clone();
}


/// @brief set build manager; also clears any cached fragments
void VarLengthBuild::manager( BuildManager const & manager ) {
	clear_fragments();
	manager_ = manager;
}


/// @brief get rid of all rcgs
void
VarLengthBuild::clear_rcgs() {
	rcgs_.clear();
}


/// @brief adding an rcg. this will set the rcg internal vlb backpointer
void
VarLengthBuild::add_rcg( RemodelConstraintGeneratorOP rcg )
{

	//we should probably also clone this
	rcg->set_vlb( VarLengthBuildAP( utility::pointer::static_pointer_cast< VarLengthBuild >( get_self_ptr() ) ) );

	rcgs_.push_back( rcg );
}

void
VarLengthBuild::clear_setup_movers() {
	setup_movers_.clear();
}

void
VarLengthBuild::add_setup_mover( moves::MoverOP mover_in )
{
	//maybe we should clone the mover?
	setup_movers_.push_back( mover_in );
}

void
VarLengthBuild::clear_user_provided_movers() {
	user_provided_movers_.clear();
}

void
VarLengthBuild::add_user_provided_mover( moves::MoverOP mover_in )
{
	//maybe we should clone the mover?
	user_provided_movers_.push_back( mover_in );
}


/// @brief clear any currently cached fragments
void VarLengthBuild::clear_fragments() {
	fragments_picked_ = false;
	fragfull_.reset(); // to NULL
	frag9_.reset(); // to NULL
	frag3_.reset(); // to NULL
	frag1_.reset(); // to NULL
}


/// @brief run protocol on given Pose
/// @return if procedure successful, return Pose with modifications and a
///  sealed fold tree, otherwise return Pose with modifications and the
///  in-progress cut fold tree
/// @remarks Before invoking this function it's best to make sure
///  the secondary structure in the Pose is marked via the method
///  that you would prefer, e.g. by Dssp (protocols::jumping::Dssp),
///  by the old Rosetta++ binning method (core::pose::set_ss_from_phipsi)
///  or by external method such as reading in a file.
void VarLengthBuild::apply( Pose & pose ) {
	using namespace core::chemical;
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::fold_tree_from_pose;
	using protocols::forge::methods::restore_residues;

	// reset status flags
	set_last_move_status( FAIL_DO_NOT_RETRY );

	// sanity check
	if ( !manager_.empty() ) {
		set_last_move_status( MS_SUCCESS );
	}

	// Cache name of the residue type set the BuildInstructions are
	// using -- currently assume they are all equivalent, as this is
	// enforced in the BuildInstruction compatibility check.
	String bi_rts_name = ( **manager_.begin() ).residue_type_set().name();

	// make backup copy for e.g. side-chain transferal later
	// REPEAT: also used for monomeric repeat
	Pose archive_pose = pose;

	utility::vector1<Real> cached_phi, cached_psi, cached_omega;

	// alter pose
	Original2Modified original2modified; // keep track of old -> new mapping
	if ( get_last_move_status() == MS_SUCCESS ) {
		// alter residue type set if necessary
		if ( pose.residue( 1 ).residue_type_set()->name() != bi_rts_name ) {
			core::util::switch_to_residue_type_set( pose, bi_rts_name );
		}

		// modify
		if ( restart_mode_ ) { // assume Pose has already been modified by the manager somewhere

			//manager_.dummy_modify( pose.n_residue() );
			//original2modified = manager_.original2modified();

			//since archive_pose is actually a copy of the modified version of the
			//pose it will have 1-to-1 mapping to the pose.  So need to make a fake
			//original2modified map if running in restart_mode
			Original2Modified dummy = manager_.original2modified();
			for ( Original2Modified::const_iterator it = dummy.begin(), end=dummy.end(); it != end; ++ it ) {
				//TR << "idx restart mode " << it->second << std::endl;
				if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
					if ( (*it).second <= (remodel_data_.sequence.length()*2) ) {
						original2modified[(*it).second] = (*it).second;
					}
				} else {
					original2modified[(*it).second] = (*it).second;
				}
			}

		} else {
			original2modified = manager_.modify( pose );
		}
	}
	if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {

		Size len_start = remodel_data_.sequence.length();

		// length should honor the blueprint length, and subsequently if the input
		// pose is longer than blueprint, allow copying of phi-psi beyond the last
		// residue. only for repeats with matching blueprint and pdb lengths.
		//cache the phi psi angles
		//if (len_start < pose.total_residue() && len_start * (basic::options::option[basic::options::OptionKeys::remodel::repeat_structure]) == pose.total_residue() ){
		if ( len_start < pose.total_residue() ) {
			for ( Size i = len_start+1; i <= pose.total_residue(); i++ ) {
				cached_phi.push_back( pose.phi( i ) );
				cached_psi.push_back( pose.psi( i ) );
				cached_omega.push_back( pose.omega( i ));
			}
			//std::cout << "HERE1-1" << std::endl;

			Size max_pdb_index = remodel_data_.blueprint.size()*2;
			//std::cout << "max_pdb index" << max_pdb_index <<  std::endl;

			while ( pose.total_residue() > max_pdb_index ) {
				//std::cout << "pose total" << pose.total_residue() <<  std::endl;
				pose.conformation().delete_residue_slow(pose.total_residue());
			}

			//similarly update archive pose in repeat cases
			while ( archive_pose.total_residue() > max_pdb_index ) {
				archive_pose.conformation().delete_residue_slow(archive_pose.total_residue());
			}
		}
	}

	// REPEAT: used for fragment picking and others
	repeat_tail_length_ =0;
	if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
		// adding a tail to the starter monomer pose
		if ( pose.total_residue() < (remodel_data_.sequence.length()*2) ) {

			Size len_diff = (2*remodel_data_.sequence.length()) - pose.total_residue();
			// append a tail of the same length
			for ( Size i = 1; i<= len_diff; i++ ) {
				//core::chemical::ResidueTypeSet const & rsd_set = (pose.residue(2).residue_type_set());
				core::chemical::ResidueTypeSetCOP const & rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
				core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set->name_map("VAL") ) );
				pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,pose.total_residue(), true);
				pose.conformation().insert_ideal_geometry_at_polymer_bond(pose.total_residue()-1);
				pose.set_omega(pose.total_residue()-1,180);
			}
		} else if ( pose.total_residue() > (remodel_data_.sequence.length()*2) ) {
			while ( pose.total_residue() != (2* remodel_data_.sequence.length()) ) {
				pose.conformation().delete_residue_slow(pose.total_residue());
			}
		}

		// similarly for archive_pose, process the same way for length
		// adding a tail to the starter monomer pose
		if ( archive_pose.total_residue() < (remodel_data_.sequence.length()*2) ) {

			Size len_diff = (2*remodel_data_.sequence.length()) - archive_pose.total_residue();
			// append a tail of the same length
			for ( Size i = 1; i<= len_diff; i++ ) {
				//core::chemical::ResidueTypeSet const & rsd_set = (archive_pose.residue(2).residue_type_set());
				core::chemical::ResidueTypeSetCOP const & rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
				core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set->name_map("VAL") ) );
				archive_pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,archive_pose.total_residue(), true);
				archive_pose.conformation().insert_ideal_geometry_at_polymer_bond(archive_pose.total_residue()-1);
				archive_pose.set_omega(archive_pose.total_residue()-1,180);
			}
		} else if ( archive_pose.total_residue() > (remodel_data_.sequence.length()*2) ) {
			while ( archive_pose.total_residue() != (2* remodel_data_.sequence.length()) ) {
				archive_pose.conformation().delete_residue_slow(archive_pose.total_residue());
			}
		}

		// take care of terminal types
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, pose.total_residue());
		core::pose::add_variant_type_to_pose_residue( archive_pose, core::chemical::UPPER_TERMINUS_VARIANT, archive_pose.total_residue());


		assert( pose.total_residue() == (2* remodel_data_.sequence.length()));
		repeat_tail_length_ = remodel_data_.sequence.length();

		//update the new lengthened pose with pose angles.
		for ( Size i = remodel_data_.sequence.length()+1, j=1; i<= pose.total_residue(); i++,j++ ) {
			if ( cached_phi.size() != 0 ) {
				pose.set_phi(i, cached_phi[j]);
				pose.set_psi(i, cached_psi[j]);
				pose.set_omega(i, cached_omega[j]);
			} else {
				pose.set_phi(i, 150);
				pose.set_psi(i, 150);
				pose.set_omega(i, 180);
			}
		}
	}

	//flo may '12, give user supplied setup movers
	//a chance to modify the pose before centroid building happens
	//pose.dump_pdb("vlb_bef_setup_movers.pdb");
	for (  utility::vector1< moves::MoverOP >::iterator mover_it( setup_movers_.begin() ); mover_it != setup_movers_.end(); ++mover_it ) {
		(*mover_it)->apply( pose );
	}
	//pose.dump_pdb("vlb_aft_setup_movers.pdb");

	// fix the abego string to deal with ligands that might have been added to the pose
	if ( abego_.size() > 0 ) {
		TR << "Abego size=" << abego_.size() << " and pose size=" << pose.total_residue() << std::endl;
		//runtime_assert( abego_.size() <= pose.total_residue() );
		// account for possible ligands in the pose
		for ( core::Size i=abego_.size(); i<pose.total_residue(); ++i ) {
			abego_.push_back("X");
		}
	}

	// centroid level protocol
	if ( get_last_move_status() == MS_SUCCESS ) {
		if ( pose.residue( 1 ).residue_type_set()->name() != core::chemical::CENTROID ) {
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
		}
		//pose.dump_pdb("vlb_bef_centroid_build.pdb");
		if ( centroid_build( pose ) ) {
			set_last_move_status( MS_SUCCESS );
		} else {
			pose=archive_pose;
			set_last_move_status( FAIL_DO_NOT_RETRY );
			return;
		}
	}
	//pose.dump_pdb("vlb_aft_centroid_build.pdb");
	//archive_pose.dump_pdb("arc_pose_vlb_aft_centroid_build.pdb");

	// flip back to prior residue type set if necessary
	if ( pose.residue( 1 ).residue_type_set()->name() != archive_pose.residue( 1 ).residue_type_set()->name() ) {
		core::util::switch_to_residue_type_set( pose, archive_pose.residue( 1 ).residue_type_set()->name() );
	}

	if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
		if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure] == 1 ) {
			// do nothing
		} else {
			//remove the added residue
			//pose.conformation().delete_residue_slow(pose.total_residue());
			//need to extend the archive pose, otherwise the connectivity is wrong

			for ( Size ii = 1; ii<=remodel_data_.sequence.length(); ++ii ) {
				ResidueType const & rsd_type(pose.residue_type(ii));
				Size res1 = ii;
				Size res2 = ii+remodel_data_.sequence.length();
				replace_pose_residue_copying_existing_coordinates(archive_pose,res1,rsd_type);
				replace_pose_residue_copying_existing_coordinates(archive_pose,res2,rsd_type);
			}

			using namespace protocols::loops;
			using protocols::forge::methods::intervals_to_loops;

			std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();

			LoopsOP loops( new Loops() );

			//Temporary fix: special case for denovo type, needed artificially
			//padding residues.  in reality all loops would be padded, but for
			//archive pose extension, it is only special for de novo case
			//std::cout << (*(loop_intervals.begin())).left << "left" << std::endl;
			//std::cout << (*(loop_intervals.begin())).right << "right" << std::endl;
			if ( loop_intervals.size() == 1 && (*(loop_intervals.begin())).left == 1 && (*(loop_intervals.begin())).right == remodel_data_.blueprint.size() ) {
				// std::cout << "man handle loop." << std::endl;
				loops->add_loop( Loop(1, remodel_data_.blueprint.size()+2, 0, 0, true) );
			} else {
				loops = LoopsOP( new Loops( intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) ) );
			}

			Pose bufferPose(archive_pose);
			protocols::forge::remodel::RemodelLoopMover RLM(loops);
			RLM.set_repeat_tail_length(repeat_tail_length_);
			RLM.repeat_generation_with_additional_residue( bufferPose, archive_pose );
			if ( option[ OptionKeys::symmetry::symmetry_definition].user() ) {
				//symmetrize if both rep+sym are used
				simple_moves::symmetry::SetupForSymmetryMover pre_mover;
				pre_mover.apply(archive_pose);
				archive_pose.pdb_info()->obsolete(true);
			}

		}
	}
	// recover side chains in fixed regions

	restore_residues( original2modified, archive_pose, pose );

	// finalize wrt to success/failure
	if ( get_last_move_status() == MS_SUCCESS ) {
		if ( !basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::bypass_closure]() ) {
			// seal the tree
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				pose.fold_tree ( core::pose::symmetry::sealed_symmetric_fold_tree( pose ) );
			} else {
				if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {}
				else {
					pose.fold_tree( fold_tree_from_pose( pose, pose.fold_tree().root(), MoveMap() ) );
				}
			}
		}
	} else if ( recover_original_on_failure_ ) {

		// user desires the original Pose ; unadressed issue with symmetry/repeat
		// propagation, as archive_pose has grown to full length, it's not correct
		// to reassign it to pose
		pose = archive_pose;

	}
	// pose.dump_pdb("vlb_aft_everything.pdb");
	// std::cout << "fold tree after everything: " << pose.fold_tree() << std::endl;
	// std::cout << "cutpoints: " << pose.fold_tree().cutpoints() << std::endl;
}


std::string
VarLengthBuild::get_name() const {
	return "VarLengthBuild";
}

/// @brief run centroid level protocol on given Pose
/// @return true if regions modeled within tolerances, false otherwise
bool VarLengthBuild::centroid_build( Pose & pose ) {

	using core::fragment::picking_old::FragmentLibraryManager;
	using core::kinematics::MoveMap;
	using protocols::loops::Loop;

	using protocols::forge::methods::count_cutpoints;
	using protocols::forge::methods::find_cutpoint;
	using protocols::forge::methods::linear_chainbreak;

	typedef utility::vector1< Interval > Intervals;

	// safety, clear the energies object
	pose.energies().clear();

	// grab new secondary structure
	String ss;
	if ( !new_secondary_structure_override_.empty() ) { // user has overridden the string auto-setup
		/*
		// first check to make sure length of the override string corresponds
		// to the working pose
		if ( new_secondary_structure_override_.length() != pose.n_residue() ) {
		// fail fast
		TR.Error << "ERROR: new secondary structure override string not equal in length to newly modified Pose."
		<< " Expected " << pose.n_residue() << " but found " << new_secondary_structure_override_.length()
		<< ".  Programmer mistake!" << std::endl;
		runtime_assert( false );
		}
		*/
		ss = new_secondary_structure_override_;

	} else { // auto-setup, take directly from the pose
		ss = pose.secstruct();
	}

	// construct new amino acid sequence, possibly empty
	String aa;
	if ( !new_sequence_override_.empty() ) { // user has overridden the string auto-setup

		// first check to make sure length of the override string corresponds
		// to the working pose
		//  if ( new_sequence_override_.length() != pose.n_residue() ) {
		//   // fail fast
		//   TR.Error << "ERROR: new sequence override string not equal in length to newly modified Pose."
		//    << " Expected " << pose.n_residue() << " but found " << new_sequence_override_.length()
		//    << ".  Programmer mistake!" << std::endl;
		//   runtime_assert( false );
		//  }

		aa = new_sequence_override_;

	} else if ( !original_sequence_.empty() ) { // auto-setup

		aa.append( pose.n_residue(), '.' ); // do dummy fill w/ '.' first

		// get mapping of positions that exist in both original sequence and
		// "new" modified pose and fill in the original sequence info
		Original2Modified o2m = manager_.original2modified();
		for ( Original2Modified::const_iterator i = o2m.begin(), ie = o2m.end(); i != ie; ++i ) {
			aa.at( i->second - 1 ) = original_sequence_.at( i->first - 1 );
		}

		// fill in the sequence of the new regions from the "new" modified Pose
		// that we're working on
		Positions np = manager_.new_positions();
		for ( Positions::const_iterator i = np.begin(), ie = np.end(); i != ie; ++i ) {
			aa.at( (*i) - 1 ) = pose.residue( *i ).name1();
		}
	}

	// data to feed to other parts of rosetta
	Intervals fragment_only_regions;
	loops::LoopsOP loops( new loops::Loops() );

	// identify regions to rebuild and pick fragments
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();

	if ( basic::options::option[basic::options::OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user() ) {
		//unfortunately hacky...  pre-process interval sets to make sure fragments
		//aren't picked for the insertion region;  should move this processing to
		//buildManager at some point

		//find insertion
		Size insertStartIndex = remodel_data_.dssp_updated_ss.find_first_of("I");
		Size insertEndIndex = remodel_data_.dssp_updated_ss.find_last_of("I");

		//loop over the interval set to find insertion and split it into two sections
		for ( std::set< Interval >::iterator i = loop_intervals.begin(), ie = loop_intervals.end(); i != ie; ++i ) {
			Interval interval = *i;

			if ( interval.left <= insertStartIndex && interval.right >= insertEndIndex && ((insertEndIndex-insertStartIndex) != 0) ) {
				//found insertion

				//create new left interval
				Interval split_interval_left(interval.left, insertStartIndex);
				//create new right interval
				Interval split_interval_right(insertEndIndex, interval.right);

				//insert the intervals to loop_intervals definition
				loop_intervals.insert(split_interval_left);
				loop_intervals.insert(split_interval_right);

				//delete the current interval
				loop_intervals.erase( i );

				break; //expect only one insertion, so can jump out if found one.
			}
		}

		// pick fragments for insertion case
		for ( std::set< Interval >::const_iterator i = loop_intervals.begin(), ie = loop_intervals.end(); i != ie; ++i ) {
			Interval interval = *i;
			if ( !( cache_fragments_ && fragments_picked_ ) ) {
				pick_all_fragments( ss, aa, abego_, interval, num_fragpick_ );
			}
		}
		fragments_picked_ = true;
	}

	//after finishing picking fragments, revert the interval, in case the
	//intervals are split by pdb insertion.  need to revert because the insertion
	//should be part of a single loop and not flanked by two connection rebuilt
	//loops
	loop_intervals = manager_.intervals_containing_undefined_positions();

	for ( std::set< Interval >::const_iterator i = loop_intervals.begin(), ie = loop_intervals.end(); i != ie; ++i ) {
		Interval interval = *i;

		Size n_cuts = count_cutpoints( pose, interval.left, interval.right );

		if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
			//if (interval.right == pose.total_residue()){
			if ( interval.right == remodel_data_.blueprint.size() ) {
				//interval.right = interval.right + repeat_tail_length_; // pad interval to include the extra shadow residue in pose
				interval.right = interval.right + 2; // pad interval to include the extra shadow residue in pose
			}
		}
		TR << "VLB count_cutpoints " << n_cuts << " interval.left " << interval.left << " interval.right " << interval.right << std::endl;

		// multi-cutpoint region handling not implemented yet
		runtime_assert( n_cuts < 2 );

		// setup regions
		if ( interval.left != 1 && interval.right != pose.n_residue() ) { //internal loop

			Size cutpoint = find_cutpoint( pose, interval.left, interval.right );

			loops->add_loop( Loop( interval.left, interval.right, cutpoint, 0.0, true ) );
			if ( cutpoint == 0 ) {
				loops->choose_cutpoints(pose);
			}

		} else if ( n_cuts == 0 ) { // fragment only region

			if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
				loops->add_loop( Loop( interval.left, interval.right, 0, 0.0, true ) );//pick additional frame for connection to repeats
			} else {
				loops->add_loop( Loop( interval.left, interval.right, 0, 0.0, true ) );
			}
		}

		// Pick fragments.  Choose num_fragpick_ ( default 200 ) per position per size.
		if ( !( cache_fragments_ && fragments_picked_ ) ) {
			//pick_all_fragments( ss, aa, interval, num_fragpick_ );
			pick_all_fragments( ss, aa, abego_, interval, num_fragpick_ );
		}
	}

	// we're done picking fragments; report status and do memory management
	fragments_picked_ = true;


	if ( use_fullmer_ ) {
		TR << "total full-mer fragments: " << fragfull_->size() << std::endl;
	}
	TR << "total 9-mer fragments: " << frag9_->size() << std::endl;
	TR << "total 3-mer fragments: " << frag3_->size() << std::endl;
	TR << "total 1-mer fragments: " << frag1_->size() << std::endl;

	switch ( vall_memory_usage_ ) {
	case VLB_VallMemoryUsage::KEEP_IN_MEMORY :
		break;
	case VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS :
		if ( cache_fragments_ ) {
			/// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
			FragmentLibraryManager::get_instance()->clear_Vall();
		}
		break;
	case VLB_VallMemoryUsage::ALWAYS_CLEAR :
		/// WARNING WARNING WARNING! THREAD UNSAFE!  WHY WOULD YOU THINK THIS IS A GOOD IDEA?
		FragmentLibraryManager::get_instance()->clear_Vall();
		break;
	default :
		break;
	}

	//setup eventual remodel constraints
	setup_remodel_constraints( pose );

	if ( (!ignore_cmdline_enzdes_cstfile_) && basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ) {

		protocols::forge::remodel::RemodelEnzdesCstModuleOP cstOP( new protocols::forge::remodel::RemodelEnzdesCstModule(remodel_data_) );

		//safety
		pose.remove_constraints();
		//wipe out cst_cache
		protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> set_cst_cache( NULL );
		//wipe out observer too
		pose.observer_cache().set( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER, NULL , false);

		//RemodelEnzdesCstModule cst(remodel_data);
		cstOP->use_backbone_only_blocks();
		cstOP->apply(pose);
		cstOP->enable_constraint_scoreterms(sfx_);
	}

	if ( basic::options::option[ basic::options::OptionKeys::constraints::cst_file ].user() && !basic::options::option[basic::options::OptionKeys::remodel::cst_fa_only].user() ) {
		//safety
		pose.remove_constraints();

		protocols::simple_moves::ConstraintSetMoverOP constraint( new protocols::simple_moves::ConstraintSetMover() );
		constraint->apply( pose );

		sfx_->set_weight(core::scoring::atom_pair_constraint, 1.0);
		sfx_->set_weight(core::scoring::dihedral_constraint, 10.0);
	}


	// setup loop building protocol
	MoverOP loop_mover = loop_mover_instance( loops, manager_.movemap() );

	// Run loop modeling.  The loop movers return the original fold tree
	// after apply().  Do we want that...?  There's also no good way to
	// check that the loop mover actually finished with a closed loop,
	// which is troubling.  For now we work around this by post-evaluating
	// the chainbreak at the original cutpoint.
	loop_mover->apply( pose );

	//remove the remodel constraints
	remove_remodel_constraints( pose );

	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;

	if ( basic::options::option[basic::options::OptionKeys::remodel::RemodelLoopMover::bypass_closure]() ) {
		return cbreaks_pass;
	}

	for ( Loops::const_iterator l = loops->begin(), le = loops->end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = linear_chainbreak( pose, l->cut() );
			TR << "centroid_build: loop: " << l->start() << "-" << l->stop() << ", final chainbreak = " << c << ", max_linear_chainbreak_: " << max_linear_chainbreak_ << std::endl;
			TR << "centroid_build: final chainbreak at " << l->cut() << " = " << c << " max tolerance " << max_linear_chainbreak_ <<  std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}

	return cbreaks_pass;
}


/// @param[in] loops The loops to model.
/// @param[in] false_mm Enforce False settings in this MoveMap.  Currently
///  only useful with the RemodelLoopMover.
VarLengthBuild::MoverOP VarLengthBuild::loop_mover_instance(
	loops::LoopsOP const loops,
	MoveMap const & false_mm
)
{
	using protocols::forge::remodel::RemodelLoopMover;
	using protocols::forge::remodel::RemodelLoopMoverOP;
	using protocols::loops::loop_mover::IndependentLoopMover;
	using protocols::loops::LoopMoverFactory;

	typedef utility::pointer::shared_ptr< IndependentLoopMover > IndependentLoopMoverOP;

	MoverOP lm;

	if ( loop_mover_str_ == "RemodelLoopMover" ) { // use RemodelLoopMover
		RemodelLoopMoverOP loop_mover( new RemodelLoopMover( loops ) );
		loop_mover->remodelData(remodel_data_);
		loop_mover->scorefunction( *sfx_ );

		if ( basic::options::option[basic::options::OptionKeys::remodel::repeat_structure].user() ) {
			loop_mover->set_repeat_tail_length( repeat_tail_length_ );
		}

		if ( use_fullmer_ ) {
			loop_mover->add_fragments( fragfull_ );
		}
		loop_mover->add_fragments( frag9_ );
		loop_mover->add_fragments( frag3_ );
		loop_mover->add_fragments( frag1_ );

		loop_mover->false_movemap( false_mm );

		if ( loop_mover_fold_tree_constant_ ) loop_mover->set_keep_input_foldtree(true);

		if ( user_provided_movers_.size() != 0 ) {
			loop_mover->set_user_provided_movers( user_provided_movers_ );
			loop_mover->user_provided_movers_apply_cycle( user_provided_movers_apply_cycle_ );
		}

		lm = loop_mover;

	} else { // use a mover from the IndependentLoopMover factory

		// protocols::loops really needs to get refactored; when this happens the
		// artificial cast and setup below will change.
		IndependentLoopMoverOP loop_mover( utility::pointer::static_pointer_cast< protocols::loops::loop_mover::IndependentLoopMover > ( loops::LoopMoverFactory::get_instance()->create_loop_mover( loop_mover_str_, loops ) ) );
		loop_mover->set_scorefxn( sfx_ );
		loop_mover->set_strict_loops( true ); // no sliding window

		if ( use_fullmer_ ) {
			loop_mover->add_fragments( fragfull_ );
		}
		loop_mover->add_fragments( frag9_ );
		loop_mover->add_fragments( frag3_ );
		loop_mover->add_fragments( frag1_ );

		lm = loop_mover;
	}

	return lm;
}


/// @brief pick fragments of size full, 9, 3, 1
/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
///            can be empty.  If empty, sequence bias is not used to pick fragments.
/// @param[in] interval The interval [left, right] to pick fragments from; Pose
///  numbering (i.e. 1-based indexing).
/// @param[in] n_frags The number of fragments to pick per position.
void VarLengthBuild::pick_all_fragments(
	String const & complete_ss,
	String const & complete_aa,
	utility::vector1< String > const & complete_abego,
	Interval const & interval,
	Size const n_frags
)
{
	using core::fragment::ConstantLengthFragSet;
	using core::fragment::OrderedFragSet;

	using protocols::forge::methods::smallmer_from_largemer;

	// pick full-mers
	if ( use_fullmer_ ) {
		if ( !fragfull_.get() ) {
			fragfull_ = OrderedFragSetOP( new OrderedFragSet() );
		}

		fragfull_->add( pick_fragments( complete_ss, complete_aa, complete_abego, interval, interval.length(), n_frags ) );
	}

	// pick 9-mers
	if ( !frag9_.get() ) {
		frag9_ = ConstantLengthFragSetOP( new ConstantLengthFragSet( 9 ) );
	}
	if ( basic::options::option[basic::options::OptionKeys::in::file::frag9].user() ) {
		frag9_->read_fragment_file( basic::options::option[basic::options::OptionKeys::in::file::frag9]());
	} else {
		frag9_->add( pick_fragments( complete_ss, complete_aa, complete_abego, interval, 9, n_frags ) );
	}

	// pick 3-mers
	if ( !frag3_.get() ) {
		frag3_ = ConstantLengthFragSetOP( new ConstantLengthFragSet( 3 ) );
	}
	ConstantLengthFragSetOP tmp_frag3( new ConstantLengthFragSet( 3 ) );
	if ( basic::options::option[basic::options::OptionKeys::in::file::frag3].user() ) {
		tmp_frag3->read_fragment_file( basic::options::option[basic::options::OptionKeys::in::file::frag3]());
		frag3_->read_fragment_file( basic::options::option[basic::options::OptionKeys::in::file::frag3]());
	} else {
		tmp_frag3->add( pick_fragments( complete_ss, complete_aa, complete_abego, interval, 3, n_frags ) );
		frag3_->add( *tmp_frag3 );
	}

	// make 1-mers from 3-mers
	if ( !frag1_.get() ) {
		frag1_ = ConstantLengthFragSetOP( new ConstantLengthFragSet( 1 ) );
	}
	frag1_->add( *smallmer_from_largemer( tmp_frag3->begin(), tmp_frag3->end(), 1 ) );
}


/// @brief pick fragments of a given length, padding when necessary
/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
///            can be empty.  If empty, sequence bias is not used to pick fragments.
/// @param[in] interval The interval [left, right] to pick fragments from; Pose
///  numbering (i.e. 1-based indexing).
/// @param[in] frag_length The desired length of the fragments
/// @param[in] n_frags The number of fragments to pick per position.
VarLengthBuild::FrameList VarLengthBuild::pick_fragments(
	String const & complete_ss,
	String const & complete_aa,
	utility::vector1< String > const & complete_abego,
	Interval const & interval,
	Size const frag_length,
	Size const n_frags
)
{
	using core::fragment::Frame;
	using core::fragment::FrameOP;
	using core::fragment::IndependentBBTorsionSRFD;

	using core::fragment::picking_old::vall::pick_fragments;
	using core::fragment::picking_old::vall::pick_fragments_by_ss;
	using core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa;

	FrameList frames;

	for ( Size j = 0, je = interval.length(); j < je; ++j ) {
		TR << "picking " << n_frags << " " << frag_length << "-mers for position " << ( interval.left + j ) << std::endl;

		String ss_sub = complete_ss.substr( interval.left + j - 1, frag_length );
		if ( ss_sub.length() < frag_length ) {
			ss_sub.append( frag_length - ss_sub.length(), 'D' );
		}

		String aa_sub;
		if ( !complete_aa.empty() ) {
			aa_sub = complete_aa.substr( interval.left + j - 1, frag_length );
			if ( aa_sub.length() < frag_length ) {
				aa_sub.append( frag_length - aa_sub.length(), '.' );
			}
		} else {
			aa_sub = "";
		}
		TR << "complete_ss length: " << complete_ss.length() << " complete_abego size: " << complete_abego.size() << std::endl;
		utility::vector1< String > abego_sub;
		if ( complete_abego.size() > 0 ) {
			//runtime_assert( complete_ss.length() == complete_abego.size() );
			//causing issues with symmetry
			Size pos( 1 );
			abego_sub.resize( frag_length );
			for ( Size ii = interval.left + j; ii <= interval.left + j + frag_length - 1; ++ii, ++pos ) {
				if ( ii > complete_abego.size() ) {
					abego_sub[ pos ] = "X";
				} else {
					abego_sub[ pos ] = complete_abego[ ii ];
				}
			}
		} else {
			abego_sub.clear(); // make sure it is empty
		}

		FrameOP frame( new Frame( interval.left + j, frag_length ) );

		frame->add_fragment( pick_fragments( ss_sub, aa_sub, abego_sub, n_frags, true, IndependentBBTorsionSRFD() ) );

		// pick wrt sec.struct and possibly sequence bias
		// if ( !complete_aa.empty() ) {
		//  String aa_sub = complete_aa.substr( interval.left + j - 1, frag_length );
		//  if ( aa_sub.length() < frag_length ) {
		//   aa_sub.append( frag_length - aa_sub.length(), '.' );
		//  }
		//  frame->add_fragment( pick_fragments_by_ss_plus_aa( ss_sub, aa_sub, n_frags, true, IndependentBBTorsionSRFD() ) );
		// } else {
		//  frame->add_fragment( pick_fragments_by_ss( ss_sub, n_frags, true, IndependentBBTorsionSRFD() ) );
		// }

		frames.push_back( frame );
	}

	return frames;
}


/// @brief telling all rcgs to setup and insert constraints into the pose
void
VarLengthBuild::setup_remodel_constraints(
	Pose & pose )
{
	if ( rcgs_.empty() ) {
		return; // nothing to do...
	}

	for ( utility::vector1< RemodelConstraintGeneratorOP >::iterator rcg_it = rcgs_.begin(); rcg_it != rcgs_.end(); ++rcg_it ) {
		(*rcg_it)->set_seqmap( this->manager().sequence_mapping());
		TR << "Adding remodel constraints to pose using " << (*rcg_it)->get_name() << std::endl;
		(*rcg_it)->add_remodel_constraints_to_pose( pose );
	}
}


/// @brief telling all rcgs to remove their constraints from the pose
void
VarLengthBuild::remove_remodel_constraints(
	Pose & pose )
{
	if ( rcgs_.empty() ) {
		return; // nothing to do...
	}

	for ( utility::vector1< RemodelConstraintGeneratorOP >::iterator rcg_it = rcgs_.begin(); rcg_it != rcgs_.end(); ++rcg_it ) {
		(*rcg_it)->remove_remodel_constraints_from_pose( pose );
	}
}


} // components
} // forge
} // protocols
