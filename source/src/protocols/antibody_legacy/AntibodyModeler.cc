// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @author Aroop Sircar ( aroopsircar@yahoo.com )


// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>

#include <core/chemical/VariantType.hh>
#include <core/fragment/FragData.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/fragment/FragID.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameList.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::format;

#include <protocols/jd2/ScoreMap.hh>
// AUTO-REMOVED #include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/antibody_legacy/AntibodyClass.hh>
#include <protocols/antibody_legacy/AntibodyModeler.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/rigid/RB_geometry.hh>
//#include <protocols/evaluation/PoseEvaluator.hh>
//#include <protocols/evaluation/RmsdEvaluator.hh>
// AUTO-REMOVED #include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/jd2/Job.hh>
// AUTO-REMOVED #include <protocols/jd2/JobOutputter.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody_legacy/CDRH3Modeler.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/antibody_legacy/GraftMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/TrialMover.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.AntibodyModeler" );

namespace protocols {
namespace antibody_legacy {
// Aroop's Magic number, do not change it
// (and dont try and use it anywhere else)

using namespace core;

// default constructor
AntibodyModeler::AntibodyModeler() : Mover(),
	init_for_input_yet_( false ) {
	Mover::type( "AntibodyModeler" );
	set_default();
	init_from_options();
}

// default destructor
AntibodyModeler::~AntibodyModeler() {}

//clone
protocols::moves::MoverOP AntibodyModeler::clone() const {
	return( protocols::moves::MoverOP( new AntibodyModeler() ) );
}

void AntibodyModeler::set_default() {
	TR <<  "Setting up default settings" << std::endl;
	model_h3_ = true;
	snugfit_ = true;
	native_present_ = false;
	graft_l1_ = true;
	graft_l2_ = true;
	graft_l3_ = true;
	graft_h1_ = true;
	graft_h2_ = true;
	graft_h3_ = false;
	benchmark_ = false;
	camelid_ = false;
	camelid_constraints_ = false;

}

void AntibodyModeler::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	TR <<  "Reading Options" << std::endl;
	model_h3_ = option[ OptionKeys::antibody::model_h3 ]();
	snugfit_ = option[ OptionKeys::antibody::snugfit ]();
	native_present_ = option[ OptionKeys::in::file::native ].user();
	if( native_present_ )
		native_filename_ = option[ OptionKeys::in::file::native ]();
	graft_l1_ = option[ OptionKeys::antibody::graft_l1 ]();
	graft_l2_ = option[ OptionKeys::antibody::graft_l2 ]();
	graft_l3_ = option[ OptionKeys::antibody::graft_l3 ]();
	graft_h1_ = option[ OptionKeys::antibody::graft_h1 ]();
	graft_h2_ = option[ OptionKeys::antibody::graft_h2 ]();
	graft_h3_ = option[ OptionKeys::antibody::graft_h3 ]();
	benchmark_ = option[ OptionKeys::run::benchmark ]();
	camelid_ = option[ OptionKeys::antibody::camelid ]();
	camelid_constraints_ = option[ OptionKeys::antibody::
	                               camelid_constraints ]();
	cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
	if( camelid_ ) {
		graft_l1_ = false;
		graft_l2_ = false;
		graft_l3_ = false;
		snugfit_ = false;
	}
	if( camelid_constraints_ )
		model_h3_ = false;

}

void AntibodyModeler::set_snugdock_foldtree( pose::Pose & pose_in ) {
	TR << "setting up Snug Dock fold tree" << std::endl;
	TR << "Snug Dock Fold Tree: " << std::endl;
	TR << pose_in.fold_tree() << std::endl;
	return;
}

void AntibodyModeler::init_on_new_input() {
	init_for_input_yet_ = true;

	// read native_pose
	if ( native_present_ ) {
		core::import_pose::pose_from_pdb(	native_pose_, native_filename_	);
		pose::set_ss_from_phipsi( native_pose_ );
	} else
		native_pose_ = start_pose_;

	antibody_in_.set_Fv( start_pose_, camelid_ );

	if( model_h3_ ) {
		// Read standard Rosetta fragments file
		read_and_store_fragments();

		// Read in CDR H3 C-terminal fragment file
		read_H3_cter_fragment( antibody_in_,
		                       H3_base_library_,
		                       camelid_);
	}

	return;
} // init_on_new_input


void AntibodyModeler::apply( pose::Pose & pose_in ) {
	using namespace core;
	using namespace chemical;
	using namespace id;
	using namespace fragment;
	using namespace scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;

	if( model_h3_ && ( cst_weight_ != 0.00 ) ) {
		protocols::simple_moves::ConstraintSetMoverOP cdr_constraint( new protocols::simple_moves::ConstraintSetMover() );
		cdr_constraint->apply( pose_in );
	}
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	start_pose_ = pose_in;

	if( !init_for_input_yet_ ) init_on_new_input();

	pose::set_ss_from_phipsi( pose_in );

	// display constraints and return
	if( camelid_constraints_ ) {
		display_constraint_residues();
		return;
	}

	// Junk
	/*
	Size tag = 0;
	for( Size i = antibody_in_.cdrh_[3][2] - 4;
			 i <= antibody_in_.cdrh_[3][2] + 1; i++ ) {

		Real phi = pose_in.phi( i );
		Real psi = pose_in.psi( i );
		Real omega = pose_in.omega( i );
		Real size = (antibody_in_.cdrh_[3][2] - antibody_in_.cdrh_[3][1])+1;

		if( i == (antibody_in_.cdrh_[3][2] + 1) ) {
			std::cout << "OUTPUT: n+1, " << antibody_in_.Fv.pdb_info()->number(i)
								<< antibody_in_.Fv.pdb_info()->icode(i) << ", "
								<< antibody_in_.Fv.residue(i).name1()
								<< ", " << omega << ", " << phi
								<< ", " << psi << ", " << size << std::endl;
		}
		else if( i == antibody_in_.cdrh_[3][2] ) {
			std::cout << "OUTPUT: n, " << antibody_in_.Fv.pdb_info()->number(i)
								<< antibody_in_.Fv.pdb_info()->icode(i)
								<< ", " << antibody_in_.Fv.residue(i).name1()
								<< ", " << omega << ", " << phi
								<< ", " << psi << ", " << size << std::endl;
		}
		else {
			std::cout << "OUTPUT: n-" << antibody_in_.cdrh_[3][2] - i
								<< ", " << antibody_in_.Fv.pdb_info()->number(i)
								<< antibody_in_.Fv.pdb_info()->icode(i)
								<< ", "<< antibody_in_.Fv.residue(i).name1()
								<< ", "	<< omega << ", " << phi
								<< ", " << psi << ", " << size << std::endl;
		}
	}
	return;
	*/
	// End Junk

	SequenceMoverOP model_antibody( new SequenceMover() );

	GraftMoverOP graft_move( new GraftMover() );
	graft_move->enable_graft_l1( graft_l1_ );
	graft_move->enable_graft_l2( graft_l2_ );
	graft_move->enable_graft_l3( graft_l3_ );
	graft_move->enable_graft_h1( graft_h1_ );
	graft_move->enable_graft_h2( graft_h2_ );
	graft_move->enable_graft_h3( graft_h3_ );
	graft_move->enable_benchmark_mode( benchmark_ );
	graft_move->set_camelid( camelid_ );
	graft_move->set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new pose::Pose ( native_pose_ ) ) ) );
	model_antibody->add_mover( graft_move );

	if ( model_h3_ ) {
		CDRH3ModelerOP model_cdrh3( new CDRH3Modeler( offset_frags_ ) );
		model_cdrh3->enable_benchmark_mode( benchmark_ );
		model_cdrh3->set_camelid( camelid_ );
		model_cdrh3->model_h3( model_h3_ );
		model_cdrh3->store_H3_cter_fragment( H3_base_library_ );
		model_cdrh3->set_native_pose( core::pose::PoseCOP( core::pose::PoseOP( new pose::Pose ( native_pose_ ) ) ) );
		model_cdrh3->set_centroid_loop_building( true );
		model_cdrh3->set_fullatom_loop_building( true );
		model_antibody->add_mover( model_cdrh3 );
	}

	model_antibody->apply( antibody_in_.Fv );

	if ( !camelid_ && snugfit_ ) {
		all_cdr_VL_VH_fold_tree( antibody_in_.Fv, antibody_in_.all_cdr_loops );
		relax_cdrs();
		repulsive_ramp ( antibody_in_.Fv, antibody_in_.all_cdr_loops );
		snugfit_mcm_protocol ( antibody_in_.Fv, antibody_in_.all_cdr_loops );

		// align Fv to native.Fv
		pose::Pose native_pose;
		if( get_native_pose() )
			native_pose = *get_native_pose();
		else
			native_pose = antibody_in_.Fv;
		Antibody native_ab( native_pose, camelid_ );
		antibody_in_.align_to_native( native_ab );
	} else if ( model_h3_ )
		relax_cdrs();

	pose_in = antibody_in_.Fv;

	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( pose_in, true );

	// Define CDR H3 loop
	Size frag_size = (antibody_in_.cdrh_[3][2]-antibody_in_.cdrh_[3][1]) + 3;
	Size cutpoint =  antibody_in_.cdrh_[3][1] + int( frag_size / 2 );
	loops::Loop cdr_h3( antibody_in_.cdrh_[3][1], antibody_in_.cdrh_[3][2],
	                    cutpoint,	0, false );

	// Fold Tree
	simple_one_loop_fold_tree( pose_in, cdr_h3 );

	// Redefining CDR H3 cutpoint variants
	loops::add_single_cutpoint_variant( pose_in, cdr_h3 );

	// score functions
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );
	scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.00 );

	// add scores to map for outputting constraint score
	( *scorefxn )( pose_in );
	protocols::jd2::ScoreMap::nonzero_energies(score_map_, scorefxn, pose_in);
	Real constraint_score = score_map_[ "atom_pair_constraint" ];

	// removing constraint score
	scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.00 );
	// add scores to map for output
	( *scorefxn )( pose_in );
	protocols::jd2::ScoreMap::nonzero_energies(score_map_, scorefxn, pose_in);

	score_map_[ "AA_H3" ] = global_loop_rmsd( pose_in, native_pose_, "h3" );
	score_map_[ "AB_H2" ] = global_loop_rmsd( pose_in, native_pose_, "h2" );
	score_map_[ "AB_H1" ] = global_loop_rmsd( pose_in, native_pose_, "h1" );
	if( !camelid_ ) {
		score_map_[ "AC_L3" ] = global_loop_rmsd( pose_in, native_pose_, "l3");
		score_map_[ "AD_L2" ] = global_loop_rmsd( pose_in, native_pose_, "l2");
		score_map_[ "AE_L1" ] = global_loop_rmsd( pose_in, native_pose_, "l1");
	}
	score_map_[ "AF_constraint" ] = constraint_score;

	//using pose::datacache::CacheableDataType::SCORE_MAP;
	using namespace basic::datacache;
	pose_in.data().set( core::pose::datacache::CacheableDataType::SCORE_MAP, DataCache_CacheableData::DataOP( new basic::datacache::DiagnosticData(score_map_) ));
}// end apply

std::string
AntibodyModeler::get_name() const {
	return "AntibodyModeler";
}

void
AntibodyModeler::read_and_store_fragments() {
	using namespace chemical;
	using namespace id;
	using namespace fragment;
	using namespace core::scoring;

	if ( !model_h3_ )
		return;

	// fragment initialization
	utility::vector1< FragSetOP > frag_libs;

	protocols::loops::read_loop_fragments( frag_libs );

	Size frag_size = (antibody_in_.cdrh_[3][2]-antibody_in_.cdrh_[3][1]) + 3;
	Size cutpoint =  antibody_in_.cdrh_[3][1] + int( frag_size / 2 );
	setup_simple_fold_tree(  antibody_in_.cdrh_[3][1] - 1, cutpoint,
	                         antibody_in_.cdrh_[3][2] + 1,
	                         antibody_in_.Fv.total_residue(),
	                         antibody_in_.Fv );

	FragSetOP offset_3mer_frags;
	// a fragset of same type should be able to handle everything
	offset_3mer_frags = frag_libs[2]->empty_clone();
	FrameList loop_3mer_frames;
	Size offset = 0;
	frag_libs[2]->region_simple( 1, frag_size, loop_3mer_frames );
	for ( FrameList::const_iterator it = loop_3mer_frames.begin(),
	        eit = loop_3mer_frames.end(); it!=eit; ++it ) {
		FrameOP short_frame = (*it)->clone_with_frags();
		offset++;
		short_frame->shift_to( ( antibody_in_.cdrh_[3][1] - 2 ) + offset  );
		offset_3mer_frags->add( short_frame );
	}

	FragSetOP offset_9mer_frags;
	// a fragset of same type should be able to handle everything
	offset_9mer_frags = frag_libs[1]->empty_clone();
	FrameList loop_9mer_frames;
	offset = 0;
	frag_libs[1]->region_simple( 1, frag_size, loop_9mer_frames );
	for ( FrameList::const_iterator it = loop_9mer_frames.begin(),
	        eit = loop_9mer_frames.end(); it!=eit; ++it ) {
		FrameOP short_frame = (*it)->clone_with_frags();
		offset++;
		short_frame->shift_to( ( antibody_in_.cdrh_[3][1] - 2 ) + offset  );
		offset_9mer_frags->add( short_frame );
	}

	offset_frags_.push_back( offset_9mer_frags );
	offset_frags_.push_back( offset_3mer_frags );

	return;
} // read_and_store_fragments

void
AntibodyModeler::setup_simple_fold_tree(
    Size jumppoint1,
    Size cutpoint,
    Size jumppoint2,
    Size nres,
    pose::Pose & pose_in ) {

	using namespace kinematics;

	TR << "ABM Setting up simple fold tree" << std::endl;

	FoldTree f;
	f.clear();

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
	f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR << "ABM Done: Setting up simple fold tree" << std::endl;

} // setup_simple_fold_tree

///////////////////////////////////////////////////////////////////////////
/// @begin relax_cdrs
///
/// @brief relaxes all cdrs simultaneously
///
/// @detailed based on the all_cdrs loop definiton, minimizes only
///           those regions. A standard dfpmin is utilized with the
///           given score function and chain -break and chain-overlap
///           set. The allow_bb/chi arrays are changed accordingly but
///           then are reset to their initial states before exiting
///           the routine. Similarly the fold tree and jump movements
///           are restored to their initial states
///
/// @param[out]
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors Aroop 02/15/2010
///
/// @last_modified 02/15/2010
///////////////////////////////////////////////////////////////////////////
void
AntibodyModeler::relax_cdrs() {
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::moves;
	// Storing initial fold tree
	kinematics::FoldTree const input_tree( antibody_in_.Fv.fold_tree() );

	// changing to all cdr fold tree
	antibody_in_.all_cdr_fold_tree();

	// adding cutpoint variants for chainbreak score computation
	loops::add_cutpoint_variants( antibody_in_.Fv );

	Size const nres = antibody_in_.Fv.total_residue();

	//setting MoveMap
	kinematics::MoveMapOP allcdr_map;
	allcdr_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	allcdr_map->clear();
	allcdr_map->set_chi( false );
	allcdr_map->set_bb( false );
	utility::vector1< bool> is_flexible( nres, false );
	bool include_neighbors( false );
	select_loop_residues( antibody_in_.Fv, antibody_in_.all_cdr_loops,
	                      include_neighbors, is_flexible );
	allcdr_map->set_bb( is_flexible );
	include_neighbors = true;
	antibody_in_.Fv.update_residue_neighbors(); // Need to have updated neighbors
	select_loop_residues( antibody_in_.Fv, antibody_in_.all_cdr_loops,
	                      include_neighbors, is_flexible );
	allcdr_map->set_chi( is_flexible );
	for( Size ii = 1; ii <= antibody_in_.all_cdr_loops.num_loop(); ii++ )
		allcdr_map->set_jump( ii, false );

	// score functions
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::get_score_function();
	scorefxn->set_weight( core::scoring::chainbreak, 10. / 3. );
	scorefxn->set_weight( core::scoring::overlap_chainbreak, 10. / 3. );

	Real min_tolerance = 0.1;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
	bool nb_list = true;
	protocols::simple_moves::MinMoverOP all_cdr_min_mover( new protocols::simple_moves::MinMover( allcdr_map,
	        scorefxn, min_type, min_tolerance, nb_list ) );
	all_cdr_min_mover->apply( antibody_in_.Fv );

	if( !benchmark_ ) {
		protocols::simple_moves::PackRotamersMoverOP repack( new protocols::simple_moves::PackRotamersMover( scorefxn ) );
		setup_packer_task( antibody_in_.Fv );
		( *scorefxn )( antibody_in_.Fv );
		tf_->push_back( TaskOperationCOP( new RestrictToInterface( is_flexible ) ) );
		repack->task_factory( tf_ );
		repack->apply( antibody_in_.Fv );

		protocols::simple_moves::RotamerTrialsMinMoverOP rtmin( new protocols::simple_moves::RotamerTrialsMinMover(
		    scorefxn, tf_ ) );
		rtmin->apply( antibody_in_.Fv );
	}

	// Restoring pose fold tree
	antibody_in_.Fv.fold_tree( input_tree );

	return;
} // relax_cdrs

///////////////////////////////////////////////////////////////////////////
/// @begin all_cdr_VL_VH_fold_tree
///
/// @brief change to all CDR and VL-VH dock fold tree
///
/// @detailed
///
/// @param[out]
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors Aroop 07/13/2010
///
/// @last_modified 07/13/2010
///////////////////////////////////////////////////////////////////////////
void
AntibodyModeler::all_cdr_VL_VH_fold_tree(
    pose::Pose & pose_in,
    const loops::Loops & loops_in ) {

	using namespace kinematics;

	Size nres = pose_in.total_residue();
	core::pose::PDBInfoCOP pdb_info = pose_in.pdb_info();
	char second_chain = 'H';
	Size rb_cutpoint(0);

	for ( Size i = 1; i <= nres; ++i ) {
		if( pdb_info->chain( i ) == second_chain) {
			rb_cutpoint = i-1;
			break;
		}
	}

	Size jump_pos1 ( core::pose::residue_center_of_mass( pose_in, 1,
	                 rb_cutpoint ) );
	Size jump_pos2 ( core::pose::residue_center_of_mass( pose_in,rb_cutpoint+1,
	                 nres ) );

	// make sure rb jumps do not reside in the loop region
	for( loops::Loops::const_iterator it= loops_in.begin(),
	        it_end = loops_in.end(); it != it_end; ++it ) {
		if ( jump_pos1 >= ( it->start() - 1 ) &&
		        jump_pos1 <= ( it->stop() + 1) )
			jump_pos1 = it->stop() + 2;
		if ( jump_pos2 >= ( it->start() - 1 ) &&
		        jump_pos2 <= ( it->stop() + 1) )
			jump_pos2 = it->start() - 2;
	}

	// make a simple rigid-body jump first
	setup_simple_fold_tree(jump_pos1,rb_cutpoint,jump_pos2,nres, pose_in );

	// add the loop jump into the current tree,
	// delete some old edge accordingly
	FoldTree f( pose_in.fold_tree() );

	for( loops::Loops::const_iterator it=loops_in.begin(),
	        it_end=loops_in.end(); it != it_end; ++it ) {
		Size const loop_start ( it->start() );
		Size const loop_stop ( it->stop() );
		Size const loop_cutpoint ( it->cut() );
		Size edge_start(0), edge_stop(0);
		//bool edge_found = false;
		const FoldTree & f_const = f;
		Size const num_jump = f_const.num_jump();
		for( FoldTree::const_iterator it2=f_const.begin(),
		        it2_end=f_const.end(); it2 !=it2_end; ++it2 ) {
			edge_start = std::min( it2->start(), it2->stop() );
			edge_stop = std::max( it2->start(), it2->stop() );
			if ( ! it2->is_jump() && loop_start > edge_start
			        && loop_stop < edge_stop ) {
				//edge_found = true;  // set but never used ~Labonte
				break;
			}
		}

		f.delete_unordered_edge( edge_start, edge_stop, Edge::PEPTIDE);
		f.add_edge( loop_start-1, loop_stop+1, num_jump+1 );
		f.add_edge( edge_start, loop_start-1, Edge::PEPTIDE );
		f.add_edge( loop_start-1, loop_cutpoint, Edge::PEPTIDE );
		f.add_edge( loop_cutpoint+1, loop_stop+1, Edge::PEPTIDE );
		f.add_edge( loop_stop+1, edge_stop, Edge::PEPTIDE );
	}

	f.reorder(1);
	pose_in.fold_tree(f);

	return;

} // all_cdr_VL_VH_fold_tree

///////////////////////////////////////////////////////////////////////////
/// @begin repulsive_ramp
///
/// @brief ramping up the fullatom repulsive weight slowly to allow the
///        partners to relieve clashes and make way for each other
///
/// @detailed This routine is specially targetted to the coupled
///           optimization of docking partners and the loop region.  The
///           loop modelling & all previous  steps  involve mainly
///           centroid  mode .On switching  on fullatom mode, one is bound
///           to end up with clashes.To relieve the clashes, it is
///           essential to slowly  dial up the  repulsive weight instead of
///           turning it on to the maximum  value in one single step
///
/// @param[in] input pose which is assumed to have a docking fold tree
///
/// @global_read fa_rep : fullatom repulsive weight
///
/// @global_write fa_rep ( It is reset to the original value at the end )
///
/// @remarks A particular portion is  commented out,which can be
///          uncommented if one  uses a  low  resolution  homology  model.
///          Check details in the beginning of the commented out region
///
/// @references
///
/// @authors Aroop 07/13/2010
///
/// @last_modified 07/13/2010
///////////////////////////////////////////////////////////////////////////
void
AntibodyModeler::repulsive_ramp(
    pose::Pose & pose_in,
    loops::Loops loops_in ) {

	Size nres = pose_in.total_residue();

	//setting MoveMap
	kinematics::MoveMapOP cdr_dock_map;
	cdr_dock_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	cdr_dock_map->clear();
	cdr_dock_map->set_chi( false );
	cdr_dock_map->set_bb( false );
	utility::vector1< bool> is_flexible( nres, false );
	bool include_neighbors( false );
	select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
	cdr_dock_map->set_bb( is_flexible );
	include_neighbors = true;
	antibody_in_.Fv.update_residue_neighbors(); // Need to have updated neighbors
	select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
	cdr_dock_map->set_chi( is_flexible );
	cdr_dock_map->set_jump( 1, true );
	for( Size ii = 2; ii <= loops_in.num_loop() + 1; ii++ )
		cdr_dock_map->set_jump( ii, false );


	// score functions
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::ScoreFunctionFactory::
	           create_score_function( "docking", "docking_min" );
	scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );

	// score functions
	core::scoring::ScoreFunctionOP pack_scorefxn;
	pack_scorefxn = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( pose_in, true );

	using namespace core::chemical;
	for ( loops::Loops::const_iterator it = loops_in.begin(),
	        it_end = loops_in.end();	it != it_end; ++it ) {
		core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_LOWER, it->cut() );
		core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_UPPER,it->cut()+1);
	}
	// add scores to map
	( *scorefxn )( pose_in );

	// dampen fa_rep weight
	Real rep_weight_max = scorefxn->get_weight( core::scoring::fa_rep );
	Size rep_ramp_cycles(3);
	Size cycles(4);
	Real minimization_threshold(15.0);
	//Real func_tol(1.0);
	//mjo commenting out 'nb_list' because it is unused and causes a warning
	//bool nb_list( true );
	if( benchmark_ ) {
		rep_ramp_cycles = 1;
		cycles = 1;
		minimization_threshold = 150.0;
		//func_tol = 10.0;  // set but never used ~Labonte
	}

	Real rep_ramp_step = (rep_weight_max - 0.02) / Real(rep_ramp_cycles-1);
	for ( Size i = 1; i <= rep_ramp_cycles; i++ ) {
		Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
		scorefxn->set_weight( core::scoring::fa_rep, rep_weight );
		snugfit_MC_min ( pose_in, cdr_dock_map, cycles, minimization_threshold,
		                 scorefxn, pack_scorefxn, is_flexible);
	}

	return;
} // repulsive_ramp

void
AntibodyModeler::snugfit_MC_min (
    pose::Pose & pose_in,
    kinematics::MoveMapOP cdr_dock_map,
    Size cycles,
    Real minimization_threshold,
    core::scoring::ScoreFunctionOP scorefxn,
    core::scoring::ScoreFunctionOP pack_scorefxn,
    utility::vector1< bool> is_flexible ) {

	using namespace moves;
	bool nb_list = true;
	Size nres = pose_in.total_residue();

	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( cdr_dock_map, scorefxn,
	        "dfpmin_armijo_nonmonotone", minimization_threshold, nb_list ) );

	//set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_perturb( new rigid::RigidBodyPerturbMover(pose_in,
	        *cdr_dock_map, 2.0, 0.1 , rigid::partner_downstream, true ) );

	setup_packer_task( pose_in );
	//set up sidechain movers for rigid body jump and loop & neighbors
	utility::vector1_size rb_jump;
	rb_jump.push_back( 1 );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	// selecting movable c-terminal residues
	ObjexxFCL::FArray1D_bool loop_residues( nres, false );
	for( Size i = 1; i <= nres; i++ )
		loop_residues( i ) = is_flexible[ i ]; // check mapping
	using namespace protocols::toolbox::task_operations;
	tf_->push_back( TaskOperationCOP( new RestrictToInterface( rb_jump, loop_residues ) ) );

	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
	    pack_scorefxn, tf_ ) );
	SequenceMoverOP rb_mover( new SequenceMover );
	rb_mover->add_mover( rb_perturb );
	rb_mover->add_mover( pack_rottrial );
	JumpOutMoverOP rb_mover_min( new JumpOutMover( rb_mover, min_mover,
	        scorefxn,	minimization_threshold ) );

	Real temperature = 0.8;
	MonteCarloOP mc( new MonteCarlo( pose_in, *scorefxn, temperature ) );
	TrialMoverOP rb_mover_min_trial( new TrialMover( rb_mover_min, mc) );
	RepeatMoverOP first_mcm_cycles( new RepeatMover( rb_mover_min_trial,
	        cycles ) );
	first_mcm_cycles->apply( pose_in );

	return;

} // snugfit_MC_min

void
AntibodyModeler::snugfit_mcm_protocol(
    pose::Pose & pose_in,
    loops::Loops loops_in ) {

	using namespace moves;
	bool nb_list = true;
	Size nres = pose_in.total_residue();

	//MC move
	Real trans_mag ( 0.1 );
	Real rot_mag ( 5.0 );

	// rb minimization
	std::string min_type = "dfpmin_armijo_nonmonotone";
	Real min_threshold ( 15.0 ); /* score unit */

	// score functions
	using namespace core::scoring;
	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = core::scoring::ScoreFunctionFactory::
	           create_score_function( "docking", "docking_min" );
	scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::overlap_chainbreak, 10./3. );

	// score functions
	core::scoring::ScoreFunctionOP pack_scorefxn;
	pack_scorefxn = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );

	// remove cutpoints variants for all cdrs
	// "true" forces removal of variants even from non-cutpoints
	loops::remove_cutpoint_variants( pose_in, true );

	using namespace core::chemical;
	for ( loops::Loops::const_iterator it = loops_in.begin(),
	        it_end = loops_in.end();	it != it_end; ++it ) {
		core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_LOWER, it->cut() );
		core::pose::add_variant_type_to_pose_residue( pose_in, CUTPOINT_UPPER,it->cut()+1);
	}

	//setting MoveMap
	kinematics::MoveMapOP cdr_dock_map;
	cdr_dock_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	cdr_dock_map->clear();
	cdr_dock_map->set_chi( false );
	cdr_dock_map->set_bb( false );
	utility::vector1< bool> is_flexible( nres, false );
	bool include_neighbors( false );
	select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
	cdr_dock_map->set_bb( is_flexible );
	include_neighbors = true;
	antibody_in_.Fv.update_residue_neighbors(); // Need to have updated neighbors
	select_loop_residues( pose_in, loops_in, include_neighbors, is_flexible);
	cdr_dock_map->set_chi( is_flexible );
	cdr_dock_map->set_jump( 1, true );
	for( Size ii = 2; ii <= loops_in.num_loop() + 1; ii++ )
		cdr_dock_map->set_jump( ii, false );


	//set up minimizer movers
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( cdr_dock_map, scorefxn, min_type,
	        min_threshold, nb_list ) );

	//set up rigid body movers
	rigid::RigidBodyPerturbMoverOP rb_perturb( new rigid::RigidBodyPerturbMover( pose_in,
	        *cdr_dock_map, rot_mag, trans_mag, rigid::partner_downstream, true ) );

	setup_packer_task( pose_in );
	//set up sidechain movers for rigid body jump and loop & neighbors
	utility::vector1_size rb_jump;
	rb_jump.push_back( 1 );
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	// selecting movable c-terminal residues
	ObjexxFCL::FArray1D_bool loop_residues( nres, false );
	for( Size i = 1; i <= nres; i++ )
		loop_residues( i ) = is_flexible[ i ]; // check mapping
	using namespace protocols::toolbox::task_operations;
	tf_->push_back( TaskOperationCOP( new RestrictToInterface( rb_jump, loop_residues ) ) );



	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
	    pack_scorefxn, tf_ ) );

	protocols::simple_moves::PackRotamersMoverOP pack_interface_repack( new protocols::simple_moves::PackRotamersMover(
	    pack_scorefxn ) );
	pack_interface_repack->task_factory(tf_);

	Real temperature = 0.8;
	MonteCarloOP mc( new MonteCarlo( pose_in, *scorefxn, temperature ) );

	TrialMoverOP pack_interface_trial( new TrialMover(pack_interface_repack,
	        mc ) );

	protocols::docking::SidechainMinMoverOP scmin_mover( new
	protocols::docking::SidechainMinMover( core::scoring::ScoreFunctionOP( pack_scorefxn ), core::pack::task::TaskFactoryCOP( tf_ ) ) );
	TrialMoverOP scmin_trial( new TrialMover( scmin_mover, mc ) );

	SequenceMoverOP rb_mover( new SequenceMover );
	rb_mover->add_mover( rb_perturb );
	rb_mover->add_mover( pack_rottrial );

	JumpOutMoverOP rb_mover_min( new JumpOutMover( rb_mover, min_mover,
	        scorefxn, min_threshold) );
	TrialMoverOP rb_mover_min_trial( new TrialMover( rb_mover_min, mc  ) );

	SequenceMoverOP repack_step( new SequenceMover );
	repack_step->add_mover( rb_mover_min_trial );
	repack_step->add_mover( pack_interface_trial );
	repack_step->add_mover( scmin_trial );

	CycleMoverOP rb_mover_min_trial_repack( new CycleMover );
	for ( Size i=1; i < 8; ++i )
		rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
	rb_mover_min_trial_repack->add_mover( repack_step );

	//set up initial repack mover
	SequenceMoverOP initial_repack( new SequenceMover );
	initial_repack->add_mover( pack_interface_trial );
	initial_repack->add_mover( scmin_trial );

	//set up initial and final min_trial movers for docking
	TrialMoverOP minimize_trial( new TrialMover( min_mover, mc ) );

	//set up mcm cycles and mcm_repack cycles
	RepeatMoverOP mcm_four_cycles( new RepeatMover( rb_mover_min_trial, 4 ) );

	Size cycles = 3;
	if ( benchmark_ )
		cycles = 1;
	RepeatMoverOP mcm_final_cycles( new RepeatMover(
	    rb_mover_min_trial_repack, cycles ) );

	SequenceMoverOP snugfit_mcm( new SequenceMover );
	snugfit_mcm->add_mover( initial_repack );
	snugfit_mcm->add_mover( minimize_trial );
	snugfit_mcm->add_mover( mcm_four_cycles );
	snugfit_mcm->add_mover( mcm_final_cycles );
	snugfit_mcm->add_mover( minimize_trial );

	snugfit_mcm->apply ( pose_in );

	return;
} // snugfit_mcm_protocol

void
AntibodyModeler::setup_packer_task(
    pose::Pose & pose_in ) {
	using namespace pack::task;
	using namespace pack::task::operation;

	if( init_task_factory_ ) {
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory( *init_task_factory_ ) );
		TR << "AbModeler Reinitializing Packer Task" << std::endl;
		return;
	} else
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory );

	TR << "AbModeler Setting Up Packer Task" << std::endl;

	tf_->push_back( TaskOperationCOP( new OperateOnCertainResidues( ResLvlTaskOperationOP( new PreventRepackingRLT ), ResFilterOP( new ResidueLacksProperty("PROTEIN") ) ) ) );
	tf_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf_->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf_->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	tf_->push_back( TaskOperationCOP( new NoRepackDisulfides ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	init_task_factory_ = tf_;

	TR << "AbModeler Done: Setting Up Packer Task" << std::endl;

} // setup_packer_task

Real
AntibodyModeler::global_loop_rmsd (
    const pose::Pose & pose_in,
    const pose::Pose & native_pose,
    std::string cdr_type
) {

	using namespace scoring;

	Size loop_start( 1 );
	Size loop_end( antibody_in_.Fv.total_residue() );
	if( cdr_type == "l1" ) {
		loop_start = antibody_in_.cdrl_[1][1];
		loop_end = antibody_in_.cdrl_[1][2];
	} else if( cdr_type == "l2" ) {
		loop_start = antibody_in_.cdrl_[2][1];
		loop_end = antibody_in_.cdrl_[2][2];
	} else if( cdr_type == "l3" ) {
		loop_start = antibody_in_.cdrl_[3][1];
		loop_end = antibody_in_.cdrl_[3][2];
	} else if( cdr_type == "h1" ) {
		loop_start = antibody_in_.cdrh_[1][1];
		loop_end = antibody_in_.cdrh_[1][2];
	} else if( cdr_type == "h2" ) {
		loop_start = antibody_in_.cdrh_[2][1];
		loop_end = antibody_in_.cdrh_[2][2];
	} else if( cdr_type == "h3" ) {
		loop_start = antibody_in_.cdrh_[3][1];
		loop_end = antibody_in_.cdrh_[3][2];
	}

	using ObjexxFCL::FArray1D_bool;
	FArray1D_bool superpos_partner ( pose_in.total_residue(), false );

	for ( Size i = loop_start; i <= loop_end; ++i )
		superpos_partner(i) = true;

	using namespace core::scoring;
	Real rmsG = rmsd_no_super_subset( native_pose, pose_in,
	                                  superpos_partner, is_protein_CA );
	return ( rmsG );
} // global_loop_rmsd

void
AntibodyModeler::display_constraint_residues() {

	// Detecting di-sulfide bond

	Size H1_Cys(0), H3_Cys(0);

	if( antibody_in_.Fv.residue( antibody_in_.Fv.pdb_info()->pdb2pose( 'H',
	                             32 ) ).name3() == "CYS" )
		H1_Cys = antibody_in_.Fv.pdb_info()->pdb2pose( 'H', 32 );
	else if( antibody_in_.Fv.residue( antibody_in_.Fv.pdb_info()->pdb2pose(
	                                      'H', 33 ) ).name3() == "CYS" )
		H1_Cys = antibody_in_.Fv.pdb_info()->pdb2pose( 'H', 33 );

	for( Size ii = antibody_in_.cdrh_[3][1]; ii <= antibody_in_.cdrh_[3][2];
	        ii++ )
		if( antibody_in_.Fv.residue(ii).name3() == "CYS" )
			H3_Cys = ii;

	if( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) )
		TR << "CONSTRAINTS: "
		   << "AtomPair CA " << H1_Cys << " CA " << H3_Cys
		   << " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;

	// Specifying extended kink

	Size hfr_46(0), h3_closest(0);
	hfr_46 = antibody_in_.Fv.pdb_info()->pdb2pose( 'H', 46 );
	if( antibody_in_.extended_ )
		h3_closest = antibody_in_.cdrh_[3][2] - 5;
	if( h3_closest != 0 )
		TR << "CONSTRAINTS: "
		   << "AtomPair CA " << hfr_46 << " CA " << h3_closest
		   << " BOUNDED 6.5 9.1 0.7 DISTANCE; mean 8.0 sd 0.7" << std::endl;

	return;
} // display_constraint_residues



} // end antibody
} // end protocols

