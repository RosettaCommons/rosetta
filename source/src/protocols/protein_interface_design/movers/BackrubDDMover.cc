// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-27 University of Washington
// (C) 199x-27 University of California Santa Cruz
// (C) 199x-27 University of California San Francisco
// (C) 199x-27 Johns Hopkins University
// (C) 199x-27 University of North Carolina, Chapel Hill
// (C) 199x-27 Vanderbilt University

/// @file protocols/protein_interface_design/movers/BackrubDDMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/BackrubDDMover.hh>
#include <protocols/protein_interface_design/movers/BackrubDDMoverCreator.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>

// Package headers
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/protein_interface_design/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
// AUTO-REMOVED #include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <basic/options/option.hh>
#include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <protocols/scoring/Interface.hh>
#include <core/kinematics/Edge.hh>
#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/kinematics/MoveMap.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
//#include <protocols/simple_moves/BackboneMover.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>

// C++ headers
#include <map>
#include <vector>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <boost/foreach.hpp>

#include <core/scoring/ScoreFunction.hh>
#include <utility/vector0.hh>
#include <utility/keys/Key3Vector.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

static thread_local basic::Tracer TR( "protocols.protein_interface_design.BackrubDDMover" );

typedef core::Real Real;
typedef core::pose::Pose Pose;

std::string
BackrubDDMoverCreator::keyname() const
{
	return BackrubDDMoverCreator::mover_name();
}

protocols::moves::MoverOP
BackrubDDMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BackrubDDMover() );
}

std::string
BackrubDDMoverCreator::mover_name()
{
	return "BackrubDD";
}

BackrubDDMover::BackrubDDMover() :
	simple_moves::DesignRepackMover( "BackrubDD" ),
	backrub_partner1_( false ),
	backrub_partner2_( true ),
	interface_distance_cutoff_( 8.0 ),
	backrub_moves_( 1000 ),
	mc_kt_( 0.6 ),
	sidechain_move_prob_( 0.25 ),
	small_move_prob_( 0.0 ),
	bbg_move_prob_( 0.25 )
{
	residues_.clear();
}

using namespace protocols::protein_interface_design;
using namespace core;
using namespace pack::task;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

// all constants for the backrub mover were taken from Colin's backrub.cc application
BackrubDDMover::BackrubDDMover
( ScoreFunctionCOP scorefxn,
	bool const backrub_partner1,
	bool const backrub_partner2,
	core::Real const interface_distance_cutoff,
	core::Size const backrub_moves,
	core::Real const mc_kt,
	core::Real const sidechain_move_prob,
	std::vector<core::Size> const & residues
)
		: simple_moves::DesignRepackMover( "Backrub" )
{
	core::Real const mm_bend_weight( 1.0 );

	using namespace core::scoring;

	backrub_partner1_ = backrub_partner1;
	backrub_partner2_ = backrub_partner2;
	backrub_moves_ = backrub_moves;
	mc_kt_ = mc_kt;
	interface_distance_cutoff_ = interface_distance_cutoff;
	sidechain_move_prob_ = sidechain_move_prob;

	runtime_assert( backrub_moves_ );
	scorefxn_repack_ = scorefxn->clone();
	scorefxn_repack_->set_weight( mm_bend, mm_bend_weight );

	// pivot atoms default to "CA" so that non-protein atoms are not considered during backrub scoring
	using namespace basic::options;
	methods::EnergyMethodOptions emo( scorefxn_repack_->energy_method_options() );
	emo.bond_angle_central_atoms_to_score( option[ OptionKeys::backrub::pivot_atoms ] );
	scorefxn_repack_->set_energy_method_options( emo );

	residues_ = residues;
}

BackrubDDMover::~BackrubDDMover() {}

protocols::moves::MoverOP
BackrubDDMover::clone() const {
	return( protocols::moves::MoverOP( new BackrubDDMover( *this ) ) );
}

void
BackrubDDMover::apply( Pose & pose )
{
	// following are the default values that Colin uses. Ugly way of setting it up, that should probably be replaced by
	// commandline options
	utility::vector1< std::string > const pivot_atoms( 1, "CA" );
	core::Size const min_atoms( 3 );
	core::Size const max_atoms( 34 );
	//core::Real const mc_kt( .6 ); // now setable
	core::Real const sc_prob_uniform( 0.1 );

	kinematics::FoldTree const saved_ft( pose.fold_tree() );
	bool make_new_ft( false );
	BOOST_FOREACH( kinematics::Edge const edge, saved_ft ){
		if( edge.start() > edge.stop() ){
			make_new_ft = true;
			break;
		}
	}

	if( make_new_ft )
		protocols::protein_interface_design::star_fold_tree( pose );

	// backrub setup based on collin's backrub.cc

	// set up the BackrubMover
	protocols::backrub::BackrubMover backrub_mover;
	protocols::simple_moves::sidechain_moves::SidechainMover sidechain_mover;
	protocols::simple_moves::SmallMover smallmover;
	protocols::simple_moves::BBG8T3AMover bbg8t3amover;

	smallmover.nmoves( 1 );
	bbg8t3amover.factorA( 0.5 ); // values suggested by Yuan
	bbg8t3amover.factorB( 10.0 );
	// read known and unknown optimization parameters from the database
	backrub_mover.branchopt().read_database();

	core::pose::PoseCOP pose_copy( new core::pose::Pose( pose ) );
	backrub_mover.set_input_pose( pose_copy );
	backrub_mover.set_native_pose( pose_copy );
	sidechain_mover.set_input_pose( pose_copy );
	sidechain_mover.set_native_pose( pose_copy );
	smallmover.set_input_pose( pose_copy );
	smallmover.set_native_pose( pose_copy );
	bbg8t3amover.set_input_pose( pose_copy );
	bbg8t3amover.set_native_pose( pose_copy );

	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	backrub_mover.clear_segments();
	TaskFactoryOP main_task_factory;
	TaskFactoryCOP ancestral_task( task_factory() );
	if( ancestral_task )
		main_task_factory = TaskFactoryOP( new TaskFactory( *ancestral_task ) );
	else
		main_task_factory = TaskFactoryOP( new TaskFactory );
	//RestrictToInterfaceOperationOP rtio = new RestrictToInterfaceOperation;
	//rtio->interface_cutoff( 8.0 );
	if( prevent_repacking().size() ){
		operation::OperateOnCertainResiduesOP prevent_repacking_on_certain_res( new operation::OperateOnCertainResidues );
		prevent_repacking_on_certain_res->residue_indices( prevent_repacking() );
		prevent_repacking_on_certain_res->op( ResLvlTaskOperationCOP( new PreventRepackingRLT ) );
		main_task_factory->push_back( prevent_repacking_on_certain_res );
	}
	main_task_factory->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	main_task_factory->push_back( TaskOperationCOP( new IncludeCurrent ) );
	//main_task_factory->push_back( rtio );
	main_task_factory->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	main_task_factory->push_back( TaskOperationCOP( new NoRepackDisulfides ) );
	if( !backrub_partner1_ )
		main_task_factory->push_back( TaskOperationCOP( new PreventChainFromRepackingOperation( 1 ) ) );
	if( !backrub_partner2_ )
		main_task_factory->push_back( TaskOperationCOP( new PreventChainFromRepackingOperation( 2 ) ) );
	if( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() )
	{
		main_task_factory->push_back( TaskOperationCOP( new ReadResfile ) );
	}

	using ObjexxFCL::FArray1D_bool;
	utility::vector1< core::Size > resnums;
	if( pose.conformation().num_chains() == 2 ){
		Size const rb_jump( 1 );

		FArray1D_bool partner1( pose.total_residue() );
		pose.fold_tree().partition_by_jump( rb_jump, partner1 ); // partner1 is true for all residues in partner1; false o/w
		Size const begin2( pose.conformation().chain_begin( 2 ) ); // the starting residue of partner2

		protocols::scoring::Interface interface( rb_jump );
		interface.distance( interface_distance_cutoff_ );
		interface.calculate( pose );
		// list of residues to backrub
		if( !residues_.size() ) {
			bool first( true ); bool last( false ); // mark all interface residues + 1 spanning residue on each side for backrub
			for (Size i = 1; i <= pose.total_residue(); i++) {
				if( !pose.residue( i ).is_protein() ) continue;
				if( (( partner1( i ) && backrub_partner1_ ) || ( !partner1(i) && backrub_partner2_ )) &&
					interface.is_interface( i ) && (!( i==begin2-1 || i==begin2) || (backrub_partner1_ && backrub_partner2_))) {
					if( first && i != 1 && (i!= begin2 || (backrub_partner1_ && backrub_partner2_) ) && pose.residue( i-1 ).is_protein() )
						resnums.push_back( i - 1 );
					first = false; last = true;
					resnums.push_back( i );
				}
				else {
					if( last ) { resnums.push_back( i ); }
					last = false; first = true;
				}
			}
		}
		else resnums = residues_;
	}
	else if ( !residues_.size() ) { // pose does not have 2 chains, backrub all (protein)
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			if ( !pose.residue( i ).is_protein() ) continue;
			resnums.push_back( i );
		}
	}
	if( residues_.size() ){// add user-defined residues to the list of backrubbable residues
		resnums.insert( resnums.begin(), residues_.begin(), residues_.end() );
		std::sort( resnums.begin(), resnums.end() );
		std::unique( resnums.begin(), resnums.end() );
	}

/// movemap is used by smallmoves
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->clear();
	movemap->set_bb( false );
	core::Size const sm_begin( backrub_partner1_ ? pose.conformation().chain_begin( 1 ) : pose.conformation().chain_begin( 2 ) );
	core::Size const sm_end( backrub_partner2_ ? pose.conformation().chain_end( 2 ) : pose.conformation().chain_end( 1 ) );

	for( core::Size resi( sm_begin ); resi <= sm_end; ++resi )
		movemap->set_bb( resi, true );
	smallmover.movemap( movemap );
	bbg8t3amover.movemap( movemap );

  //take out the residues from the resnum vector that are not allowed to be repacked
	//as well as 1 before and 1 residue after that one
	using namespace core::scoring;
	for( utility::vector1< core::Size > ::const_iterator prev_rep = prevent_repacking_.begin(); prev_rep!=prevent_repacking_.end(); ++prev_rep ){
  	utility::vector1< core::Size >::iterator it = std::find( resnums.begin(), resnums.end(),*prev_rep );
  	if( it != resnums.end() )
  	{
    	resnums.erase( it );
  		utility::vector1< core::Size >::iterator it_next = std::find( resnums.begin(), resnums.end(),*prev_rep+1 );
			if( it_next != resnums.end() ) resnums.erase( it_next );
  		utility::vector1< core::Size >::iterator it_previous = std::find( resnums.begin(), resnums.end(),*prev_rep-1 );
			if( it_previous != resnums.end() ) resnums.erase( it_previous );
  	}
	}

	// C-beta atoms should not be altered during packing because branching atoms are optimized
	main_task_factory->push_back( TaskOperationCOP( new PreserveCBeta ) );

	// set up the SidechainMover
	sidechain_mover.set_task_factory( main_task_factory );
	sidechain_mover.set_prob_uniform( sc_prob_uniform );

	backrub_mover.add_mainchain_segments( resnums, pivot_atoms, min_atoms, max_atoms );

	core::Size const br_segments( backrub_mover.num_segments() );
	TR << "Backrub Segments Added: " << br_segments <<"\n";
	if( br_segments == 0 ){
		TR<<"No segments to backrub. skipping backrub."<<std::endl;
		return;
	}

	TR << "Score After PDB Load:" << std::endl;
	scorefxn_repack_->show(TR, pose);

	backrub_mover.optimize_branch_angles( pose );
// SJF It doesn't make sense to idealize sidechains in docking
//	sidechain_mover.idealize_sidechains( pose );

	TR << "Score After Branch Angle Optimization/Side Chain Idealization:\n" ;
	scorefxn_repack_->show( TR, pose );

	protocols::moves::MonteCarlo mc( pose, *scorefxn_repack_, mc_kt_ );

	mc.reset( pose );

	TR << "Running " << backrub_moves_ << " trials..." << std::endl;

	for ( Size i = 1; i <= backrub_moves_; ++i )
	{
		std::string move_type;

		// could use random mover for this...
		core::Real const move_prob( numeric::random::rg().uniform() );
		if ( move_prob > small_move_prob_ + bbg_move_prob_ + sidechain_move_prob_ ) {
			backrub_mover.apply( pose );
			move_type = backrub_mover.type();
		} else if( move_prob > sidechain_move_prob_ + bbg_move_prob_ ){
			smallmover.apply( pose );
			move_type = smallmover.type();
		} else if( move_prob > sidechain_move_prob_ ) {
			bbg8t3amover.apply( pose );
			move_type = bbg8t3amover.type();
		}
		else {
			sidechain_mover.apply(pose);
			move_type = sidechain_mover.type();
		}

		mc.boltzmann( pose, move_type );
	}

	mc.show_counters();

	TR << "Last Score:" << std::endl;
	scorefxn_repack_->show(TR, pose);

	pose = mc.lowest_score_pose();

	TR << "Low Score:\n";
	scorefxn_repack_->show(TR, pose);
	TR.flush();
	pose = mc.lowest_score_pose();
	// write parameters for any sets of branching atoms for which there were not optimization coefficients
	backrub_mover.branchopt().write_database();

	pose.fold_tree( saved_ft );
}

std::string
BackrubDDMover::get_name() const {
	return BackrubDDMoverCreator::mover_name();
}

void BackrubDDMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose)
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ));
	backrub_partner1_ = tag->getOption<bool>( "partner1", 0 );
	backrub_partner2_ = tag->getOption<bool>( "partner2", 1 );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff", 8.0 );
	backrub_moves_ = tag->getOption<core::Size>( "moves", 1000 );
	sidechain_move_prob_ = tag->getOption<core::Real>( "sc_move_probability", 0.25 );
	small_move_prob_ = tag->getOption< core::Real >( "small_move_probability", 0.0 );
	bbg_move_prob_ = tag->getOption< core::Real >( "bbg_move_probability", 0.25 );
	runtime_assert( sidechain_move_prob_ + small_move_prob_ + bbg_move_prob_ <= 1.0 );
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	scorefxn_repack_->set_weight( mm_bend, 1.0 );
	// pivot atoms default to "CA" so that non-protein atoms are not considered during backrub scoring
	using namespace basic::options;
	methods::EnergyMethodOptions emo( scorefxn_repack_->energy_method_options() );
	emo.bond_angle_central_atoms_to_score( option[ OptionKeys::backrub::pivot_atoms ] );
	scorefxn_repack_->set_energy_method_options( emo );

	utility::vector0< TagCOP > const & backrub_tags( tag->getTags() );
	for( utility::vector0< TagCOP >::const_iterator br_it=backrub_tags.begin(); br_it!=backrub_tags.end(); ++br_it ) {
		TagCOP const br_tag_ptr = *br_it;
		if( br_tag_ptr->getName() == "residue" ) {
			core::Size const resnum( core::pose::get_resnum( br_tag_ptr, pose ) );
			residues_.push_back( resnum );
		}
		if( br_tag_ptr->getName() == "span" ) {
			string const begin_str( br_tag_ptr->getOption<string>( "begin" ) );
			string const end_str( br_tag_ptr->getOption<string>( "end" ) );
			core::Size const begin( core::pose::parse_resnum( begin_str, pose ) );
			core::Size const end( core::pose::parse_resnum( end_str, pose ) );
			runtime_assert( end > begin );
			runtime_assert( begin>=1);
			runtime_assert( end<=pose.total_residue() );
			for( core::Size i=begin; i<=end; ++i ) residues_.push_back( i );
		}
	}
	TR<<"backrub mover over residues: ";
	for( std::vector< core::Size >::const_iterator it=residues_.begin() ; it!=residues_.end(); ++it )
								TR<<*it<<" ";
	TR<<std::endl;
}


} //movers
} //protein_interface_design
} //protocols
