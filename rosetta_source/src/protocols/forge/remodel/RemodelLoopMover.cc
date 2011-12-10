// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.cc
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelLoopMoverCreator.hh>

// package headers
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/util.hh>

// project headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
//#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//#include <basic/options/keys/Remodel.OptionKeys.gen.hh>
#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// boost headers
#include <boost/format.hpp>

// C++ headers
#include <algorithm>
#include <set>
#include <sstream>
#include <utility>

#include <protocols/jd2/JobDistributor.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <boost/lexical_cast.hpp>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <utility/string_util.hh>


namespace protocols {
namespace forge {
namespace remodel {


// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.remodel.RemodelLoopMover" );

// RNG
static numeric::random::RandomGenerator RG( 9788221 ); // magic number, don't change

std::string
RemodelLoopMoverCreator::keyname() const
{
	return RemodelLoopMoverCreator::mover_name();
}

protocols::moves::MoverOP
RemodelLoopMoverCreator::create_mover() const {
	return new RemodelLoopMover;
}

std::string
RemodelLoopMoverCreator::mover_name()
{
	return "RemodelLoop";
}


/// @brief default constructor
RemodelLoopMover::RemodelLoopMover() :
	Super( "RemodelLoop"  ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	max_linear_chainbreak_( 0.07 ),
	randomize_loops_( true ),
	allowed_closure_attempts_( 1 ), //switched from 3 so no accumulation
	simultaneous_cycles_( 2 ),
	independent_cycles_( 8 ),
	boost_closure_cycles_( 30 ),
	temperature_( 2.0 )
{
	set_param_from_options();
}


/// @brief loops constructor
RemodelLoopMover::RemodelLoopMover( Loops const & loops ) :
	Super( "RemodelLoop" ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	loops_( loops ),
	max_linear_chainbreak_( 0.07 ),
	randomize_loops_( true ),
	allowed_closure_attempts_( 1 ),
	simultaneous_cycles_( 2 ),
	independent_cycles_( 8 ),
	boost_closure_cycles_( 30 ),
	temperature_( 2.0 )
{
	set_param_from_options();
}


/// @brief copy constructor
RemodelLoopMover::RemodelLoopMover( RemodelLoopMover const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	sfx_( rval.sfx_ ),
	false_movemap_( rval.false_movemap_ ),
	loops_( rval.loops_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	randomize_loops_( rval.randomize_loops_ ),
	allowed_closure_attempts_( rval.allowed_closure_attempts_ ),
	simultaneous_cycles_( rval.simultaneous_cycles_ ),
	independent_cycles_( rval.independent_cycles_ ),
	boost_closure_cycles_( rval.boost_closure_cycles_ ),
	temperature_( rval.temperature_ ),
	fragsets_( rval.fragsets_ )
{}


/// @brief default destructor
RemodelLoopMover::~RemodelLoopMover() {}


/// @brief clone this object
RemodelLoopMover::MoverOP RemodelLoopMover::clone() const {
	return new RemodelLoopMover( *this );
}


/// @brief create this type of object
RemodelLoopMover::MoverOP RemodelLoopMover::fresh_instance() const {
	return new RemodelLoopMover();
}


/// @brief set parameters from options
void RemodelLoopMover::set_param_from_options(){
	using namespace basic::options;
	using namespace OptionKeys::remodel;
	if( option[ OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak ].user() )    max_linear_chainbreak_ = option[ OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::randomize_loops ].user() )          randomize_loops_       = option[ OptionKeys::remodel::RemodelLoopMover::randomize_loops ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts ].user() ) allowed_closure_attempts_ = option[ OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::simultaneous_cycles ].user() )      simultaneous_cycles_      = option[ OptionKeys::remodel::RemodelLoopMover::simultaneous_cycles ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::independent_cycles ].user() )       independent_cycles_       = option[ OptionKeys::remodel::RemodelLoopMover::independent_cycles ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::boost_closure_cycles ].user() ) 	   boost_closure_cycles_     = option[ OptionKeys::remodel::RemodelLoopMover::boost_closure_cycles ].value();
	if( option[ OptionKeys::remodel::RemodelLoopMover::temperature ].user() ) 	           temperature_              = option[ OptionKeys::remodel::RemodelLoopMover::temperature ].value();
}

/// @registere options
void RemodelLoopMover::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
 	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::randomize_loops );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::simultaneous_cycles );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::independent_cycles );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::boost_closure_cycles );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::temperature );
}

/// @brief the ScoreFunction to use during modeling;
RemodelLoopMover::ScoreFunction const & RemodelLoopMover::scorefunction() const {
	return *sfx_;
}


/// @brief the ScoreFunction to use during modeling
void RemodelLoopMover::scorefunction( ScoreFunction const & sfx ) {
	//sfx_ = new ScoreFunction( sfx );
	sfx_ =  sfx.clone();
}


/// @brief add a fragment set
void RemodelLoopMover::add_fragments( FragSetCOP fragset ) {
	if ( fragset->size() > 0 ) {
		fragsets_.push_back( fragset->clone() );
	}
}


/// @brief clear all fragment sets
void RemodelLoopMover::clear_fragments() {
	fragsets_.clear();
}

void RemodelLoopMover::repeat_generation_with_additional_residue(Pose &pose, Pose & repeat_pose)
{
	using namespace core::pose;
	using core::Size;
	using namespace basic::options;

	//testing repeat units
	//two pose symmetry strategy

	core::pose::Pose non_terminal_pose(pose);

  if (option[ OptionKeys::remodel::repeat_structure].user()){

		//remove the extra tail from pose, but still use the pose with tail to build
		Size tail_count = repeat_tail_length_;
		while ( tail_count ){
			non_terminal_pose.conformation().delete_residue_slow(non_terminal_pose.total_residue());
			tail_count--;
		}

		repeat_pose=non_terminal_pose;
    core::pose::remove_lower_terminus_type_from_pose_residue(non_terminal_pose, 1);


    Size repeatFactor = option[ OptionKeys::remodel::repeat_structure];
    Size count = 1;
    while ( repeatFactor !=1){ // the argument should be total number of copies
      for (Size rsd = 1; rsd <= non_terminal_pose.total_residue(); rsd++){
        Size current_term = repeat_pose.total_residue();
				//intentially insert behind the last residue, this way blueprint definition will cover the junction with fragments
				if (rsd == non_terminal_pose.total_residue()){
        repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(), true);
				}else {
        repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(), false);
				/*
				for (int i =1; i<= repeat_pose.total_residue(); i++){
			std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) <<  " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
		} */
				}
      }
      Size junction = (non_terminal_pose.total_residue())* count;
			repeat_pose.conformation().insert_ideal_geometry_at_polymer_bond(junction);
			//std::cout << "junction " << junction << std::endl;
			repeat_pose.set_phi(junction+1,-150);
			repeat_pose.set_psi(junction,150);
			repeat_pose.set_omega(junction,180);
      //std::cout << "repeat Factor : " << repeatFactor << std::endl;
      count++;
      repeatFactor--;
    }
		//terminus residue check
		for (Size res=2; res< repeat_pose.total_residue(); res++){
			if (repeat_pose.residue(res).is_terminus()){
				std::cout<< "FIX TERMINUS " << res << std::endl;
				core::pose::remove_upper_terminus_type_from_pose_residue(repeat_pose, res);
				core::pose::remove_lower_terminus_type_from_pose_residue(repeat_pose, res);
			}
		}
		//take care of foldtree
		core::kinematics::FoldTree f;
		f.simple_tree(repeat_pose.total_residue());
		repeat_pose.fold_tree(f);

		//repeat_pose.dump_pdb("repeat.pdb");
		/*
		for (int i =1; i<= repeat_pose.total_residue(); i++){
			std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) << " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
		}*/
  //  pose = repeat_pose;
	//	std::cout << repeat_pose.fold_tree()<< std::endl;
  }
}


void RemodelLoopMover::repeat_generation(Pose &pose, Pose & repeat_pose)
{
	using namespace core::pose;
	using core::Size;
	using namespace basic::options;
	//testing repeat units
	//two pose symmetry strategy

	core::pose::Pose non_terminal_pose(pose);
	//core::pose::Pose repeat_pose;
  if (option[ OptionKeys::remodel::repeat_structure].user()){
		repeat_pose=non_terminal_pose;
  	core::pose::remove_lower_terminus_type_from_pose_residue(non_terminal_pose, 1);

    Size repeatFactor = option[ OptionKeys::remodel::repeat_structure];
    Size count = 1;
    while ( repeatFactor !=1){ // the argument should be total number of copies
      for (Size rsd = 1; rsd <= non_terminal_pose.total_residue(); rsd++){
        Size current_term = repeat_pose.total_residue();
				//intentially insert behind the last residue, this way blueprint definition will cover the junction with fragments
				if (rsd == non_terminal_pose.total_residue()){
        repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(),true);
				}else {
        repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(),false);
				}
      }
      Size junction = (non_terminal_pose.total_residue())* count;
        repeat_pose.conformation().insert_ideal_geometry_at_polymer_bond(junction);
        repeat_pose.set_omega(junction,180);
      //std::cout << "repeat Factor : " << repeatFactor << std::endl;
      count++;
      repeatFactor--;
    }
		//terminus residue check
		for (Size res=2; res< repeat_pose.total_residue(); res++){
			if (repeat_pose.residue(res).is_terminus()){
				std::cout<< "FIX TERMINUS " << res << std::endl;
				core::pose::remove_upper_terminus_type_from_pose_residue(repeat_pose, res);
				core::pose::remove_lower_terminus_type_from_pose_residue(repeat_pose, res);
			}
		}
		//take care of foldtree
		core::kinematics::FoldTree f;
		f.simple_tree(repeat_pose.total_residue());
		repeat_pose.fold_tree(f);
		//repeat_pose.dump_pdb("repeat.pdb");
  //  pose = repeat_pose;
	//	std::cout << repeat_pose.fold_tree()<< std::endl;
  }
}

void RemodelLoopMover::repeat_propagation( //utility function
	core::pose::Pose & pose,
	core::pose::Pose & repeat_pose,
	core::Size repeat_number
)
{
	using core::Size;
	Size segment_length = (repeat_pose.n_residue())/repeat_number;
	//std::cout << "DEBUG: segment lenght = " << segment_length << std::endl;

	for (Size rep = 0; rep < repeat_number; rep++){
		for (Size res = 1; res <= segment_length; res++){
				//std::cout << "DEBUG: res+segmentlength*rep = " << res+(segment_length*rep) << std::endl;
				Real loop_phi = 0;
				Real loop_psi = 0;
				if (res == 1 ){
					loop_phi = pose.phi(segment_length+1);
					loop_psi = pose.psi(segment_length+1);
				} else {
					loop_phi = pose.phi(res);
					loop_psi = pose.psi(res);
				}

				repeat_pose.set_phi(res+( segment_length*rep), loop_phi );
				repeat_pose.set_psi(res+( segment_length*rep), loop_psi );
				repeat_pose.set_omega( res+(segment_length*rep), pose.omega(res) );
		}
	}

	//loop over the tail fragment to the first fragment


  //pose = repeat_pose;
  //pose.dump_pdb("test_repeat1.pdb");
  //repeat_pose.dump_pdb("repeat2.pdb");
/*
		for (int i =1; i<= repeat_pose.total_residue(); i++){
			std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) << " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
		} */
}


/// @brief apply defined moves to given Pose
/// @remarks Sets protocols::moves::MS_SUCCESS upon successful closure of
///  all loops, otherwise sets protocols::moves::FAIL_RETRY.
void RemodelLoopMover::apply( Pose & pose ) {

  using namespace basic::options;
  using namespace core::scoring::constraints;
	using core::kinematics::FoldTree;
	using protocols::jd2::JobDistributor;
	using protocols::forge::methods::fold_tree_from_pose;

	// archive
	FoldTree const archive_ft = pose.fold_tree();

	FoldTree sealed_ft;
	if (core::pose::symmetry::is_symmetric(pose) ) {
		sealed_ft = core::pose::symmetry::sealed_symmetric_fold_tree( pose );
	} else {
			sealed_ft = fold_tree_from_pose( pose, pose.fold_tree().root(), MoveMap() ); // used during structure accumulation
	}
	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			repeat_generation_with_additional_residue(pose, repeat_pose_);

			//take care of constraints:

			//user defined constraint from file, make sure to do this first as it
			//replaces cst object in pose and that wipes out everything already set.

			//only use this type of cst file in this case
			if (basic::options::option[ OptionKeys::constraints::cst_file ].user()){

				protocols::simple_moves::ConstraintSetMoverOP repeat_constraint = new protocols::simple_moves::ConstraintSetMover();
				repeat_constraint->apply( repeat_pose_ );
			}

			// ResidueTypeLinkingConstraints
			Size repeat_number = basic::options::option[ OptionKeys::remodel::repeat_structure];
			Real bonus = 10;
			//std::cout << "RESIDUETYPELINKING CST" << std::endl;
		  Size segment_length = (repeat_pose_.n_residue())/repeat_number;
		  for (Size rep = 1; rep < repeat_number; rep++ ){ // from 1 since first segment don't need self-linking
			  for (Size res = 1; res <= segment_length; res++){
					 repeat_pose_.add_constraint( new ResidueTypeLinkingConstraint(repeat_pose_, res, res+(segment_length*rep), bonus));
	//				std::cout << res << " " << res+(segment_length*rep) << std::endl;
			  }
		  }
/*
			std::stringstream templateRangeSS;
			templateRangeSS << "1-" << segment_length;

			//Dihedral (NCS) Constraints
			//std::cout << "NCS CST" << std::endl;
			protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
		  for (Size rep = 1; rep < repeat_number; rep++){ // from 1 since first segment don't need self-linking
				std::stringstream targetSS;
				targetSS << 1+(segment_length*rep) << "-" << segment_length + (segment_length*rep);
		//		std::cout << templateRangeSS.str() << " " << targetSS.str() << std::endl;

					setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		  }
*/
	}

	// for accumulation of closed structures (only return the best)
	std::multimap< Real, PoseOP > accumulator;

	// setup parameters -- linearly scale the chainbreak weight
	Real const final_standard_cbreak_weight = 5.0;
	Real const cbreak_increment = final_standard_cbreak_weight / total_standard_cycles();

	// currently no scaling during boost_closure
	Real const final_boost_closure_cbreak_weight = 5.0;
	Real const boost_closure_cbreak_increment = ( final_boost_closure_cbreak_weight - final_standard_cbreak_weight ) / boost_closure_cycles();

	assert( final_boost_closure_cbreak_weight >= final_standard_cbreak_weight );

	// mark linear chainbreak in scoring function; this will be incremented
	// within simultaneous and independent stages
	ScoreFunctionOP sfxOP = sfx_;
	sfxOP->set_weight( core::scoring::linear_chainbreak, 0.0 );

	//REPEAT TEST
	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
		sfxOP->set_weight(core::scoring::atom_pair_constraint, 1.0 * basic::options::option[ OptionKeys::remodel::repeat_structure] );
		}

	// randomize loops
	if( randomize_loops_ ) {
		randomize_stage( pose );
		(*sfxOP)( pose );
	} else {
		TR << "Randomize stage was skipped " << std::endl;
	}

	// setup monte carlo
	Real const temp = temperature_;
	MonteCarlo mc( *sfxOP, temp ); // init without pose
	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
				(*sfxOP)(repeat_pose_);
	      mc.reset(repeat_pose_);
	}
	else {
			mc.reset(pose);
	}

	for ( Size attempt = 1; attempt <= allowed_closure_attempts_; ++attempt ) {
		TR << "* closure_attempt " << attempt << std::endl;

		// reset score function at the beginning of each attempt
		mc.score_function( *sfxOP );

		// simultaneous loop movements using fragment + ccd_move (default 20% of the time)
		simultaneous_stage( pose, mc, cbreak_increment );

		// closure is hard, so attempt closure for each loop independently using
		// fragment + ccd_move (default 80% of the time)
		independent_stage( pose, mc, cbreak_increment );

		// "boost": if any loops are not closed but within a chainbreak interval, attempt
		// to close them with 1-mer + ccd_move.
		boost_closure_stage( pose, mc, boost_closure_cbreak_increment );

		// check to see if all loops closed, if so rescore w/out chainbreak
		// and store in accumulator
		if ( check_closure_criteria( pose ) ) {
			//Pose temp_pose(pose);

			//make a pointer copy for storage, for REPEATs, store only the monomer
			//pose
			PoseOP pose_prime = new Pose( pose );
			pose_prime->fold_tree( sealed_ft );

			//this is for scoring in repeat context
			if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
					repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
					(*sfxOP)(repeat_pose_);

					//this has to be set because when copying to pose_prime, it lost the
					//first phi angle
					pose_prime->set_phi(1, repeat_pose_.phi(1));

					accumulator.insert( std::make_pair( repeat_pose_.energies().total_energy(), pose_prime ) );
			}	else {
				(*sfxOP)( *pose_prime );
				accumulator.insert( std::make_pair( pose_prime->energies().total_energy(), pose_prime ) );
			}

			//reset pose to monomer, if building repeats
			//pose=temp_pose;

			// now randomize the loops again for a new starting point
			if ( attempt < allowed_closure_attempts_ ) {
				randomize_stage( pose );
				if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
					repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
					(*sfxOP)(repeat_pose_);
					mc.reset( repeat_pose_ );
				} else{
						mc.reset( pose);
				}
			}
		} else {
			// Still broken, so perform a random smallest-mer insertion into each
			// loop before cycling again, otherwise too easy for the trajectory to
			// get trapped.
			insert_random_smallestmer_per_loop( pose, true );
			if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
				repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
					(*sfxOP)(repeat_pose_);
				mc.reset( repeat_pose_ );
			} else{
				mc.reset( pose);
			}
		}

//		std::ostringstream ss;
//		ss << "rlm." << attempt << ".";
//		JobDistributor::get_instance()->job_outputter()->other_pose(
//			JobDistributor::get_instance()->current_job(),
//			pose,
//			ss.str()
//		);
	}

	TR << "* " << accumulator.size() << " / " << allowed_closure_attempts_ << "   closed / attempts " << std::endl;

	// return the best structure if available, otherwise mark failure
	if ( !accumulator.empty() ) {
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
				pose = *( accumulator.begin()->second );
				repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
				(*sfxOP)(repeat_pose_);
				pose = repeat_pose_;
		} else {
			pose = *( accumulator.begin()->second );
		}

		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}

	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
		//do nothing?
	} else {
	// set original topology
		pose.fold_tree( archive_ft );
	}

}

std::string
RemodelLoopMover::get_name() const {
	return "RemodelLoopMover";
}

/// @brief randomize loops
void RemodelLoopMover::randomize_stage( Pose & pose ) {
	TR << "** randomize_stage" << std::endl;

	// simul movemap -- all loops moveable
	MoveMap movemap;
	mark_loops_moveable( loops_, movemap, true );
	enforce_false_movemap( movemap );

	// set appropriate topology
	if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
	}
	else {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, loops_, f_new );
			pose.fold_tree( f_new );
		} else {
			pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, loops_ ) );
		}
	}

	// init fragment mover for each fragment set
	FragmentMoverOPs frag_movers = create_fragment_movers( movemap );

	Size const n_moveable = count_moveable_residues( movemap, 1, pose.n_residue() );

	// insert random number of fragments = n_frag_movers * moveable_residues
	for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
		for ( Size j = 0; j < n_moveable; ++j ) {
			(*i)->apply( pose );
		}
	}

	check_closure_criteria( pose, true );
}


/// @brief find the smallest fragment size and insert a single such
///  smallmer into each loop; for breaking up trapped trajectories
/// @param[in,out] pose The pose to modify.
/// @param[in] only_broken_loop If true, only insert into broken loops,
///  otherwise insert into all. (default true)
void RemodelLoopMover::insert_random_smallestmer_per_loop(
	Pose & pose,
	bool const only_broken_loops
)
{
	using core::kinematics::FoldTree;
	using namespace basic::options;

	// determine the right set of loops to insert fragments
	Loops loops_to_model;

	if ( only_broken_loops ) {
		loops_to_model = determine_loops_to_model( pose );
	} else {
		loops_to_model = loops_;
	}

	// set appropriate topology
	if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
	}
	else {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, loops_to_model, f_new );
			pose.fold_tree( f_new );
		} else {
		pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, loops_to_model ) );
		}
	}
	// find the size of the smallest fragments
	Size smallestmer_size = ( *fragsets_.begin() )->max_frag_length();
	for ( FragSetOPs::const_iterator f = fragsets_.begin(), fe = fragsets_.end(); f != fe; ++f ) {
		smallestmer_size = std::min( smallestmer_size, (*f)->max_frag_length() );
	}

	// insert fragments
	FragmentMoverOPs frag_movers = create_per_loop_fragment_movers( loops_to_model, smallestmer_size );

	for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
		(*i)->apply( pose );
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			//Pose temp_pose(pose);
			repeat_propagation(pose, repeat_pose_,  basic::options::option[ OptionKeys::remodel::repeat_structure]);
			//pose = temp_pose;
		}
	}
}



/// @brief simultaneous stage: multiple loop movement prior to MC accept/reject
void RemodelLoopMover::simultaneous_stage(
	Pose & pose,
	MonteCarlo & mc,
	Real const cbreak_increment
)
{
	using core::kinematics::FoldTree;

 	using namespace basic::options;
  using namespace OptionKeys::remodel;
	using protocols::loops::add_cutpoint_variants;
	using protocols::loops::ccd_moves;
	using protocols::loops::remove_cutpoint_variants;
	using numeric::random::random_permutation;

	TR << "** simultaneous_stage" << std::endl;

	// Make a local copy of the Loops list.  At this stage all loops
	// are malleable so we don't use determine_loops_to_model().
	Loops loops_to_model = loops_;
	TR << "   n_loops = " << loops_to_model.size() << std::endl;

	if ( loops_to_model.size() == 0 ) { // nothing to do...
		return;
	}

	// Create fragment movers for each loop for each fragment set.  We want
	// to allow an equal probability of movement per-loop, rather than
	// per-residue.
	FragmentMoverOPs frag_movers = create_per_loop_fragment_movers( loops_to_model );
	assert( !frag_movers.empty() );

	// set appropriate topology
	if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
	}
	else {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, loops_to_model, f_new );
			pose.fold_tree( f_new );
		} else {
		pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, loops_to_model ) );
	}
}
	// add cutpoint variants
	add_cutpoint_variants( pose );
	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
		repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
		mc.reset(repeat_pose_);
	}	else{
	mc.reset( pose );
	}

	// setup master movemap covering all loops -- used only for tracking purposes
	MoveMap movemap;
	mark_loops_moveable( loops_to_model, movemap, true );
	enforce_false_movemap( movemap );

	// parameters
	Size const n_moveable = count_moveable_residues( movemap, 1, pose.n_residue() );
	Size const n_standard_cycles = total_standard_cycles();
	Size const max_outer_cycles = simultaneous_cycles();
	Size const max_inner_cycles = std::max( 50 * loops_to_model.size(), 10 * n_moveable );

	// reset counters
	mc.reset_counters();

	// simul frag + ccd_move
	for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {
		// increment the chainbreak weight
		ScoreFunctionOP sfxOP = mc.score_function().clone();
		sfxOP->set_weight(
			core::scoring::linear_chainbreak,
			sfxOP->get_weight( core::scoring::linear_chainbreak ) + cbreak_increment
		);

		if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
			sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
		}

		mc.score_function( *sfxOP );

		// recover low
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if (basic::options::option[ OptionKeys::remodel::repeat_structure] == 1){
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}

			for (Size res = 1; res<=copy_size; res++){
				pose.set_phi(res,mc.lowest_score_pose().phi(res));
				pose.set_psi(res,mc.lowest_score_pose().psi(res));
				pose.set_omega(res,mc.lowest_score_pose().omega(res));
			}
		}else{
			pose = mc.lowest_score_pose();
		}

		for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {

			if ( RG.uniform() * n_standard_cycles > outer || pose.fold_tree().num_cutpoint() == 0 ) {
				// fragments
				random_permutation( frag_movers.begin(), frag_movers.end(), RG );
				for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
					(*i)->apply( pose );
					if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
						//Pose temp_pose(pose);
						repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "simul_frag" );
						//pose=temp_pose;
					}else {
						mc.boltzmann( pose, "simul_frag" );
					}
				}
			} else {
				// per-loop ccd
				random_permutation( loops_to_model.v_begin(), loops_to_model.v_end(), RG );
				for ( Loops::const_iterator l = loops_to_model.begin(), le = loops_to_model.end(); l != le; ++l ) {
					if ( !l->is_terminal( pose ) ) {

						if (!option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
							ccd_moves( 10, pose, movemap, (int)l->start(), (int)l->stop(), (int)l->cut() );
						}
						if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
							//Pose temp_pose(pose);
							repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "ccd_move" );
						//	pose=temp_pose;
						}else {
							mc.boltzmann( pose, "ccd_move" );
						}
					}
				}
			}

		} // inner_cycles

	} // outer_cycles

	// recover low
	if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if (basic::options::option[ OptionKeys::remodel::repeat_structure] == 1){
				copy_size = pose.total_residue() - repeat_tail_length_;
			} else {
				copy_size = pose.total_residue();
			}
		for (Size res = 1; res<=copy_size; res++){
			pose.set_phi(res,mc.lowest_score_pose().phi(res));
			pose.set_psi(res,mc.lowest_score_pose().psi(res));
			pose.set_omega(res,mc.lowest_score_pose().omega(res));
		}
		//mc.lowest_score_pose().dump_pdb("simultaneous_stage.pdb");
	}
	else{
		pose = mc.lowest_score_pose();
	}

	// report status
//	mc.score_function().show_line_headers( TR );
//	TR << std::endl;

  if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
//		mc.score_function().show_line( TR, repeat_pose_ );
		mc.score_function().show( TR, repeat_pose_ );
//		repeat_pose_.dump_scored_pdb("checkRepeat1.pdb", mc.score_function());
	}	else {
		mc.score_function().show( TR, pose );
	}

	TR << std::endl;
	mc.show_state();

	check_closure_criteria( pose, true );
	remove_cutpoint_variants( pose );
}


/// @brief independent stage: single loop movement prior to MC accept/reject
void RemodelLoopMover::independent_stage(
	Pose & pose,
	MonteCarlo & mc,
	Real const cbreak_increment
)
{
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::add_cutpoint_variants;
	using namespace basic::options;
	using namespace OptionKeys::remodel;
	using protocols::loops::ccd_moves;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** independent_stage" << std::endl;

	// setup loops
	Loops loops_to_model = determine_loops_to_model( pose );
	TR << "   n_loops = " << loops_to_model.size() << std::endl;

	if ( loops_to_model.size() == 0 ) { // nothing to do...
		return;
	}

	// parameters
	Size const n_standard_cycles = total_standard_cycles();
	Size const max_outer_cycles = independent_cycles();

	// per-loop frag + ccd_move
	for ( Loops::iterator l = loops_to_model.v_begin(), le = loops_to_model.v_end(); l != le; ++l ) {
		Loop & loop = *l;

		// alter cutpoint to one before the end of the loop (either direction,
		// based on option) if closure is bypassed.  this is already set in VLB,
		// but reenforce here.
		if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
		//	if (option[OptionKeys::remodel::RemodelLoopMover::force_cutting_N].user()){
		//		loop.set_cut(loop.start()+1); //1 because remodel needs one residue flanking
		//	}
		//	else {
				//loop.set_cut(loop.stop()-1); //1 because remodel needs one residue flanking
		//	}
		}

		// movemap
		MoveMap movemap;
		mark_loop_moveable( loop, movemap, true );
		enforce_false_movemap( movemap );

		// fragment movers
		FragmentMoverOPs frag_movers = create_fragment_movers( movemap );
		assert( !frag_movers.empty() );

		// parameters
		Size const n_moveable = count_moveable_residues( movemap, loop.start(), loop.stop() );
		Size const max_inner_cycles = std::max( static_cast< Size >( 50 ), 10 * n_moveable );

		// set appropriate topology
		if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
		}
		else {
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				protocols::loops::set_single_loop_fold_tree( pose, loop );
			} else {
				protocols::forge::methods::set_single_loop_fold_tree( pose, loop );
			}
		}

		// add cutpoint variants
		add_cutpoint_variants( pose );
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
			mc.reset(repeat_pose_);
		}	else{
		mc.reset( pose );
		}
		// reset counters
		mc.reset_counters();

		// do closure
		for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {
			// increment the chainbreak weight
			ScoreFunctionOP sfxOP = mc.score_function().clone();
			sfxOP->set_weight(
				core::scoring::linear_chainbreak,
				sfxOP->get_weight( core::scoring::linear_chainbreak ) + cbreak_increment
			);

			if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			mc.score_function( *sfxOP );

			// recover low
			if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if (basic::options::option[ OptionKeys::remodel::repeat_structure] == 1){
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
				}
			}
			else{
				pose = mc.lowest_score_pose();
			}

			for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
				// fragments
				if ( loop.is_terminal( pose ) || RG.uniform() * n_standard_cycles > ( outer + simultaneous_cycles() ) || pose.fold_tree().num_cutpoint() == 0 ) {
					random_permutation( frag_movers.begin(), frag_movers.end(), RG );
					for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
						(*i)->apply( pose );
						if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
							//Pose temp_pose(pose);
							repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "frag" );
							//pose=temp_pose;
						}else {
							mc.boltzmann( pose, "frag" );
						}
					}
				} else { // ccd
					if (!option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
						ccd_moves( 10, pose, movemap, (int)loop.start(), (int)loop.stop(), (int)loop.cut() );
					}
					if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
					//	Pose temp_pose(pose);
						repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "ccd_move" );
					//	pose=temp_pose;
					}else {
					mc.boltzmann( pose, "ccd_move" );
					}
				}

			} // inner_cycles

		} // outer_cycles

		// recover low
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if (basic::options::option[ OptionKeys::remodel::repeat_structure] == 1){
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
					//mc.lowest_score_pose().dump_pdb("independent_stage.pdb");
				}
		} else{
				pose = mc.lowest_score_pose();
		}

		// report status
	//	mc.score_function().show_line_headers( TR );
//	TR << std::endl;

		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
	//		mc.score_function().show_line( TR, repeat_pose_ );
			mc.score_function().show( TR, repeat_pose_ );
	//	repeat_pose_.dump_scored_pdb("checkRepeat2.pdb", mc.score_function());
		}	else {
			mc.score_function().show( TR, pose );
		}

		TR << std::endl;
		mc.show_state();

		// remove cutpoints
		remove_cutpoint_variants( pose );
	}

	check_closure_criteria( pose, true );
}


/// @brief lock down stage: close loops within some threshold
///  w/ smallest-mer (typically 1-mer) + ccd_move only
/// @details This stage differs from simultaneous_stage() and independent_stage()
///  in that once a loop is closed, it breaks out of the closure cycle and goes
///  to the next one.  The rationale is that once we hit the boost_closure_stage()
///  we are desperate to close the loop, so we sacrifice diversity and instead
///  just seek a closed solution.
void RemodelLoopMover::boost_closure_stage(
	Pose & pose,
	MonteCarlo & mc,
	Real const cbreak_increment
)
{
	// Notes: the current implementation of boost_closure_stage() differs somewhat
	// from the equivalent, original "boost" stage within ++Remodel/EpiGraft.
	// In the old implementation, the boost cycles continued the
	// independent_stage().  As a result, the % time spent in fragment insertion
	// vs ccd continued to drop until the procedure was basically only doing
	// ccd when the total number of cycles reached 100.  In addition, the
	// chainbreak weight continued to increment until it reached 50 (at total
	// cycle 100).
	using protocols::simple_moves::ClassicFragmentMover;

	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::add_cutpoint_variants;
	using namespace basic::options;
	using namespace OptionKeys::remodel;
	using protocols::loops::ccd_moves;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** boost_closure_stage" << std::endl;

	// setup loops
	Loops pre_loops_to_model = determine_loops_to_model( pose );
	Loops loops_to_model;

	// filter for non-terminal loops that are within tolerance
	Real const cbreak_tolerance = 1.0;
	for ( Loops::const_iterator l = pre_loops_to_model.begin(), le = pre_loops_to_model.end(); l != le; ++l ) {
		if ( !l->is_terminal( pose ) ) {
			Real const cbreak = linear_chainbreak( pose, l->cut() );
			if ( cbreak < cbreak_tolerance ) {
				loops_to_model.add_loop( *l );
			}
		}
	}

	TR << "   n_loops = " << loops_to_model.size() << std::endl;
	if ( loops_to_model.size() == 0 ) { // nothing to do...
		return;
	}

	// find the size of the smallest fragments
	Size smallestmer_size = ( *fragsets_.begin() )->max_frag_length();
	for ( FragSetOPs::const_iterator f = fragsets_.begin(), fe = fragsets_.end(); f != fe; ++f ) {
		smallestmer_size = std::min( smallestmer_size, (*f)->max_frag_length() );
	}
	TR << "** boost_closure stage: using fragment size = " << smallestmer_size << std::endl;

	// Parameters.  Note that fragments often get rejected at this stage, so
	// recommend keeping the number of insertions low and the number of ccd_move
	// high.
	Size const max_outer_cycles = boost_closure_cycles();
	Real const frag_mover_probability = 0.25; // do 1-mer insertions only 25% of the time

	// 1-mer frag + ccd_move
	for ( Loops::const_iterator l = loops_to_model.begin(), le = loops_to_model.end(); l != le; ++l ) {
		Loop const & loop = *l;

		// movemap
		MoveMap movemap;
		mark_loop_moveable( loop, movemap, true );
		enforce_false_movemap( movemap );

		// prepare smallest-mer fragment movers
		FragmentMoverOPs frag1_movers = create_fragment_movers( movemap, smallestmer_size );

		// parameters
		Size const n_moveable = count_moveable_residues( movemap, loop.start(), loop.stop() );
		Size const max_inner_cycles = std::max( static_cast< Size >( 50 ), 10 * n_moveable );

		// set appropriate topology
		if (basic::options::option[basic::options::OptionKeys::remodel::no_jumps]){
		}
		else {
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				protocols::loops::set_single_loop_fold_tree( pose, loop );
			} else {
				protocols::forge::methods::set_single_loop_fold_tree( pose, loop );
			}
		}
		// add cutpoint variants
		add_cutpoint_variants( pose );
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
			mc.reset(repeat_pose_);
		}	else{
			mc.reset( pose );
		}

		// reset counters
		mc.reset_counters();

		// do closure
		for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {
			// increment the chainbreak weight
			ScoreFunctionOP sfxOP = mc.score_function().clone();
			sfxOP->set_weight(
				core::scoring::linear_chainbreak,
				sfxOP->get_weight( core::scoring::linear_chainbreak ) + cbreak_increment
			);

			if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			mc.score_function( *sfxOP );

			pose = mc.lowest_score_pose();

			// Going into the boost_closure stage implies we are "desperate" to close
			// the loop and don't care about diversity anymore, so continue to
			// cycle only until the loop is closed.
			if ( linear_chainbreak( pose, loop.cut() ) <= max_linear_chainbreak_ ) {
				break;
			}

			for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
				if ( (!frag1_movers.empty() && RG.uniform() < frag_mover_probability) || pose.fold_tree().num_cutpoint() == 0 ) { // 1-mer insertions

					random_permutation( frag1_movers.begin(), frag1_movers.end(), RG );
					for ( FragmentMoverOPs::iterator i = frag1_movers.begin(), ie = frag1_movers.end(); i != ie; ++i ) {
						(*i)->apply( pose );
						if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
						//	Pose temp_pose(pose);
							repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "frag1" );
						//	pose=temp_pose;
						}
						else{
							mc.boltzmann( pose, "frag1" );
						}
					}

				} else { // ccd_move
	      	if (!option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
						ccd_moves( 10, pose, movemap, (int)loop.start(), (int)loop.stop(), (int)loop.cut() );
					}
					if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
						//Pose temp_pose(pose);
						repeat_propagation(pose, repeat_pose_, basic::options::option[ OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "frag1" );
						//pose=temp_pose;
					} else {
						mc.boltzmann( pose, "ccd_move" );
					}
				}
			} // inner_cycles

		} // outer_cycles

		// recover low
		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if (basic::options::option[ OptionKeys::remodel::repeat_structure] == 1){
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}
					for (Size res = 1; res<=copy_size; res++){
          pose.set_phi(res,mc.lowest_score_pose().phi(res));
          pose.set_psi(res,mc.lowest_score_pose().psi(res));
          pose.set_omega(res,mc.lowest_score_pose().omega(res));
        }
				//	mc.lowest_score_pose().dump_pdb("boost_stage.pdb");
    } else{
        pose = mc.lowest_score_pose();
    }


		// report status
	//	mc.score_function().show_line_headers( TR );
//		TR << std::endl;

		if (basic::options::option[ OptionKeys::remodel::repeat_structure].user()){
	//		mc.score_function().show_line( TR, repeat_pose_ );
			mc.score_function().show( TR, repeat_pose_ );
	//	repeat_pose_.dump_scored_pdb("checkRepeat3.pdb", mc.score_function());
		}	else {
			mc.score_function().show( TR, pose );
		}

		TR << std::endl;
		mc.show_state();

		// remove cutpoint variants
		remove_cutpoint_variants( pose );
	}

	check_closure_criteria( pose, true );
}


/// @brief determine which loops need modeling wrt to given Pose
/// @remarks Skips closed loops and shuffles the order of the remaining
///  loops.
RemodelLoopMover::Loops RemodelLoopMover::determine_loops_to_model( Pose & pose ) {
	using protocols::forge::methods::linear_chainbreak;
	Loops loops_to_model;

	for ( Loops::const_iterator l = loops_.begin(), le = loops_.end(); l != le; ++l ) {
		bool skip_loop = false;

		if ( !l->is_terminal( pose ) ) {
			skip_loop |= ( linear_chainbreak( pose, l->cut() ) <= max_linear_chainbreak_ ); // loop already closed?
		}

		if ( !skip_loop ) {
			loops_to_model.add_loop( *l );
		}
	}

	// shuffle the order
	random_permutation( loops_to_model.v_begin(), loops_to_model.v_end(), RG );

	return loops_to_model;
}


/// @brief check all loops for closure criteria
/// @param[in] pose The pose being checked.
/// @param[in] show_in_tracer Output state of each loop to tracer?
/// @return true if all criteria pass, false otherwise
bool RemodelLoopMover::check_closure_criteria(
	Pose & pose,
	bool const show_in_tracer
)
{
	using protocols::forge::methods::linear_chainbreak;
	using namespace basic::options;
	using namespace OptionKeys::remodel;

	//breakout case if we don't care if the loops are closed
	if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
		return true;
	}

	// boost::format here does not appear to be doing what I want it to do...
	// The format string is probably borked.
	boost::format fmt( "%|5t|%1% %|5t|%2% %|5t|%3% %|8t|%4%" );
	if ( show_in_tracer ) {
		TR << fmt % "start" % "stop" % "cut" % "cbreak" << std::endl;
	}

	bool all_loops_pass = true;

	for ( Loops::const_iterator l = loops_.begin(), le = loops_.end(); l != le; ++l ) {
		Real cbreak = 0.0;
		if ( !l->is_terminal( pose ) ) {
			cbreak = linear_chainbreak( pose, l->cut() );
			all_loops_pass &= ( cbreak <= max_linear_chainbreak_ );
		}

		if ( show_in_tracer ) {
			TR << fmt % l->start() % l->stop() % l->cut() % cbreak << std::endl;
		}
	}

	return all_loops_pass;
}


/// @brief return fragment movers for the list of internally kept fragment sets
/// @param[in] movemap Use this movemap when initializing fragment movers.
/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
///  size is this number.  If zero, uses all fragment sets.
RemodelLoopMover::FragmentMoverOPs
RemodelLoopMover::create_fragment_movers(
	MoveMap const & movemap,
	Size const largest_frag_size
)
{
	using protocols::simple_moves::ClassicFragmentMover;

	FragmentMoverOPs frag_movers;
	for ( FragSetOPs::const_iterator f = fragsets_.begin(), fe = fragsets_.end(); f != fe; ++f ) {

		if ( largest_frag_size == 0 || (*f)->max_frag_length() <= largest_frag_size ) {
			ClassicFragmentMover * cfm = new ClassicFragmentMover( *f, movemap.clone() );
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			frag_movers.push_back( cfm );
		}

	}

	return frag_movers;
}


/// @brief append fragment movers for the list of internally kept fragment sets
/// @param[in] movemap Use this movemap when initializing fragment movers.
/// @param[out] frag_movers Append fragment movers to this list.
/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
///  size is this number.  If zero, uses all fragment sets.
void RemodelLoopMover::create_fragment_movers(
	MoveMap const & movemap,
	FragmentMoverOPs & frag_movers,
	Size const largest_frag_size
) {
	using protocols::simple_moves::ClassicFragmentMover;

	for ( FragSetOPs::const_iterator f = fragsets_.begin(), fe = fragsets_.end(); f != fe; ++f ) {

		if ( largest_frag_size == 0 || (*f)->max_frag_length() <= largest_frag_size ) {
			ClassicFragmentMover * cfm = new ClassicFragmentMover( *f, movemap.clone() );
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			frag_movers.push_back( cfm );
		}

	}
}


/// @brief create per-loop fragment movers: 1 fragment mover for each loop (uses
///  movemaps to lock down non-loop residues)
/// @param[in] loops The loops to use.
/// @param[in] largest_frag_size Only use fragment sets whose largest fragment
///  size is this number.  If zero, uses all fragment sets.
RemodelLoopMover::FragmentMoverOPs RemodelLoopMover::create_per_loop_fragment_movers(
	Loops const & loops,
	Size const largest_frag_size
)
{
	// Create fragment movers for each loop for each fragment set.  Here we
	// want to allow an equal probability of movement per-loop, rather than
	// per-residue.
	FragmentMoverOPs frag_movers;
	for ( Loops::const_iterator l = loops.begin(), le = loops.end(); l != le; ++l ) {
		MoveMap mm;
		mark_loop_moveable( *l, mm, true );
		enforce_false_movemap( mm );
		create_fragment_movers( mm, frag_movers, largest_frag_size );
	}

	return frag_movers;
}


/// @brief enforce settings in the false movemap
void RemodelLoopMover::enforce_false_movemap( MoveMap & movemap ) {
	// enforce everything in the false movemap
	movemap.import_false( false_movemap_ );
}


/// @brief mark bb/chi torsions of multiple loops moveable in a movemap
/// @param[in] loops The loops to use.
/// @param[out] movemap The movemap to modify.
/// @param[in] allow_omega Allow bb omega to move? (should be yes when
///  doing either fragment insertion or scoring function has omega
///  tether, otherwise should probably be no)
void RemodelLoopMover::mark_loops_moveable(
	Loops const & loops,
	MoveMap & movemap,
	bool const allow_omega
)
{
	for ( Loops::const_iterator l = loops.begin(), le = loops.end(); l != le; ++l ) {
		mark_loop_moveable( *l, movemap, allow_omega );
	}
}


/// @brief mark bb/chi torsion of a single loop moveable in movemap
/// @param[in] loops The loop to use.
/// @param[out] movemap The movemap to modify.
/// @param[in] allow_omega Allow bb omega to move? (should be yes when
///  doing either fragment insertion or scoring function has omega
///  tether, otherwise should probably be no)
void RemodelLoopMover::mark_loop_moveable(
	Loop const & loop,
	MoveMap & movemap,
	bool const allow_omega
)
{
	using core::id::BB;
	using core::id::omega_torsion;
	using core::id::TorsionID;

	for ( Size i = loop.start(), ie = loop.stop(); i <= ie; ++i ) {
		movemap.set_bb( i, true );
		movemap.set_chi( i, true );

		if ( !allow_omega ) {
			movemap.set( TorsionID( i, BB, omega_torsion ), false );
		}
	}
}


/// @brief count number of residues with moveable backbone torsions in the
///  given range [left, right]
RemodelLoopMover::Size RemodelLoopMover::count_moveable_residues(
	MoveMap const & movemap,
	Size const left,
	Size const right
)
{
	Size n_moveable = 0;

	// Count total number of residues w/ moveable backbone in movemap.
	// Depending on types of fragments (possible non-backbone?) consider
	// changing this in the future to check chi/other dof as well.
	for ( Size i = left; i <= right; ++i ) {
		if ( movemap.get_bb( i ) ) {
			++n_moveable;
		}
	}

	return n_moveable;
}

/// @brief parse xml
void
RemodelLoopMover::parse_my_tag(
	TagPtr const tag,
	DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	typedef utility::vector1< String > StringVec;
	using utility::string_split;
	using core::fragment::FragSet;

	// set score function
	String const sfxn ( tag->getOption<String>( "scorefxn", "" ) );
	if( sfxn != "" ) {
		sfx_ = data.get< ScoreFunction * >( "scorefxns", sfxn );
		TR << "score function, " << sfxn << ", is used. " << std::endl;
	}

	// set loops
	String const loops_string ( tag->getOption<String>( "loops", "" ) );
	runtime_assert( loops_string != "" );

	// Loop definition is given by begin-end.cutpoint, or begin-end.
	// For example, 10-23.15 indicates that the beginning(end) residue of a loop is 10(23) and
	// 15 is the cutpoint. If you omit a cutpoint input, it is given randomly between begin and end residues.
	Size num( 0 );
	StringVec loops ( string_split( loops_string, ',' ) );
	for ( StringVec::const_iterator loop( loops.begin() ), end( loops.end() ); loop!=end; ++loop ) {
		StringVec elements ( string_split( *loop, '.' ) );
		runtime_assert( elements.size() == 2 || elements.size() == 1 );
		Size cutpoint( 0 );
		StringVec interval ( string_split( elements[1], '-' ) );
		Size left = boost::lexical_cast<Size>( interval[1] );
		Size right = boost::lexical_cast<Size>( interval[2] );
		runtime_assert( left < right && interval.size() == 2 );

		// if begin or end residues are termini, there is no cutpoint
		if ( pose.residue( left ).is_lower_terminus() || pose.residue( right ).is_upper_terminus() ||
				 pose.residue( left ).is_upper_terminus() || pose.residue( right ).is_lower_terminus() ){
			cutpoint = 0;
		} else {
			if ( elements.size() == 2 ){ // if there is a input of cutpoint
				cutpoint = boost::lexical_cast<Size>( elements[2] );
			} else { // if there is no input of  cutpoint
				cutpoint = RG.random_range( 1, right - left ) + left;
			}
		}
		TR << "Modelling loop " << ++num << ": left=" << left << ",right=" << right << ",cutpoint=" << cutpoint << std::endl;
		// add a loop
		loops_.add_loop( Loop( left, right, cutpoint, 0.0, true ) );
	}

	// set fragsets
	FragSetOP fsop;
	String const fsets_string ( tag->getOption<String>( "fragsets", "" ) );
	runtime_assert( ! fsets_string.empty() );
	StringVec const fsets( string_split( fsets_string, ',' ) );
	for ( StringVec::const_iterator it( fsets.begin() ), end( fsets.end() ); it!= end; ++it )	{
		if ( data.has( "fragsets", *it ) ) {
			fsop = data.get< FragSet* >( "fragsets", *it );
		} else {
			utility_exit_with_message("fragsets " + *it + " not found in DataMap.");
		}
		add_fragments( fsop );
	}

	// set other parameters
	max_linear_chainbreak_    = tag->getOption<Real>( "max_linear_chainbreak", 0.07 );
	randomize_loops_          = tag->getOption<bool>( "randomize_loops", true );
	allowed_closure_attempts_ = tag->getOption<Size>( "allowed_closure_attempts", 1 );
	simultaneous_cycles_      = tag->getOption<Size>( "simultaneous_cycles", 2 );
	independent_cycles_       = tag->getOption<Size>( "independent_cycles", 8 );
	boost_closure_cycles_     = tag->getOption<Size>( "boost_closure_cycles", 30 );
	temperature_              = tag->getOption<Real>( "temperature", 2.0 );

}


} // remodel
} // forge
} // protocols

