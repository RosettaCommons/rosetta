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
#include <protocols/forge/remodel/RemodelData.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>

// package headers
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/methods/util.hh>

//datacache
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/io/pdb/file_data.hh>
#include <core/id/TorsionID.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/ConstantLengthFragSet.fwd.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackPotential.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackData.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <core/sequence/ABEGOManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/VallLookbackFragMover.hh>

#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/indexed_structure_store.OptionKeys.gen.hh>

//loophash
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashMap.hh>


//#include <basic/options/keys/Remodel.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

// only needed for intermediate output
//#include <protocols/jd2/JobDistributor.hh>
//#include <protocols/jd2/JobOutputter.hh>
//#include <protocols/jd2/Job.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// boost headers
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

// C++ headers
#include <algorithm>
#include <set>
#include <sstream>
#include <utility>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <boost/lexical_cast.hpp>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <utility/string_util.hh>

using namespace basic::options;

namespace protocols {
namespace forge {
namespace remodel {


// Tracer instance for this file
// Named after the original location of this code
static thread_local basic::Tracer TR( "protocols.forge.remodel.RemodelLoopMover" );

// RNG

/// @brief default constructor
RemodelLoopMover::RemodelLoopMover() :
	Super( "RemodelLoop"  ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	max_linear_chainbreak_( 0.07 ),
	randomize_loops_( true ),
	allowed_closure_attempts_( 1 ), //switched from 3 so no accumulation
	loophash_cycles_( 8 ),
	simultaneous_cycles_( 2 ),
	independent_cycles_( 8 ),
	user_provided_movers_apply_cycle_(3),
	boost_closure_cycles_( 30 ),
	temperature_( 2.0 ),
	keep_input_foldtree_( false ),
	seal_foldtree_( true )
{
	set_param_from_options();
}


/// @brief loops constructor
RemodelLoopMover::RemodelLoopMover( loops::LoopsOP const loops ) :
	Super( "RemodelLoop" ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "remodel_cen" ) ),
	loops_( loops ),
	max_linear_chainbreak_( 0.07 ),
	randomize_loops_( true ),
	allowed_closure_attempts_( 1 ),
	loophash_cycles_( 8 ),
	simultaneous_cycles_( 2 ),
	independent_cycles_( 8 ),
	user_provided_movers_apply_cycle_(3),
	boost_closure_cycles_( 30 ),
	temperature_( 2.0 ),
	keep_input_foldtree_( false ),
	seal_foldtree_( true )
{
	// Options that don't have default values.  Check before access.
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
	loophash_cycles_( rval.loophash_cycles_ ),
	simultaneous_cycles_( rval.simultaneous_cycles_ ),
	independent_cycles_( rval.independent_cycles_ ),
	user_provided_movers_(rval.user_provided_movers_),
	user_provided_movers_apply_cycle_(rval.user_provided_movers_apply_cycle_),
	boost_closure_cycles_( rval.boost_closure_cycles_ ),
	temperature_( rval.temperature_ ),
	fragsets_( rval.fragsets_ ),
	keep_input_foldtree_(rval.keep_input_foldtree_),
	seal_foldtree_(rval.seal_foldtree_)
{}


/// @brief default destructor
RemodelLoopMover::~RemodelLoopMover() {}


/// @brief clone this object
RemodelLoopMover::MoverOP RemodelLoopMover::clone() const {
	return RemodelLoopMover::MoverOP( new RemodelLoopMover( *this ) );
}


/// @brief create this type of object
RemodelLoopMover::MoverOP RemodelLoopMover::fresh_instance() const {
	return RemodelLoopMover::MoverOP( new RemodelLoopMover() );
}


/// @brief set parameters from options
void RemodelLoopMover::set_param_from_options(){

	if(option[OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak].user() )
		max_linear_chainbreak_ =option[OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak];
	if(option[OptionKeys::remodel::RemodelLoopMover::randomize_loops].user() )
		randomize_loops_       =option[OptionKeys::remodel::RemodelLoopMover::randomize_loops];
	if(option[OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts].user() )
		allowed_closure_attempts_ =option[OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts];
	if(option[OptionKeys::remodel::RemodelLoopMover::loophash_cycles].user() )
		loophash_cycles_       =option[OptionKeys::remodel::RemodelLoopMover::loophash_cycles];
	if(option[OptionKeys::remodel::RemodelLoopMover::simultaneous_cycles].user() )
		simultaneous_cycles_      =option[OptionKeys::remodel::RemodelLoopMover::simultaneous_cycles];
	if(option[OptionKeys::remodel::RemodelLoopMover::independent_cycles].user() )
		independent_cycles_       =option[OptionKeys::remodel::RemodelLoopMover::independent_cycles];
	if(option[OptionKeys::remodel::RemodelLoopMover::boost_closure_cycles].user() )
		boost_closure_cycles_     =option[OptionKeys::remodel::RemodelLoopMover::boost_closure_cycles];
	if(option[OptionKeys::remodel::RemodelLoopMover::temperature].user() )
		temperature_              =option[OptionKeys::remodel::RemodelLoopMover::temperature];
	// flo may '12 the below option is different bc it's
	// replacing the way it was used in the code until now,
	// where the check for the user didn't happen
	keep_input_foldtree_ =option[OptionKeys::remodel::no_jumps];
}

/// @registere options
void RemodelLoopMover::register_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
 	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::randomize_loops );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::allowed_closure_attempts );
	option.add_relevant( OptionKeys::remodel::RemodelLoopMover::loophash_cycles );
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

void RemodelLoopMover::remodelData(protocols::forge::remodel::RemodelData const remodel_data){
		remodel_data_ = remodel_data;
}

void
RemodelLoopMover::set_user_provided_movers( utility::vector1< moves::MoverOP > const & movers )
{
	user_provided_movers_.clear();
	user_provided_movers_ = movers;
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

	//this function currently don't extract asymmetric unit because for sym+repeat mode, it is called before symmetrizing structures

	using core::Size;
	using namespace basic::options;

	//testing repeat units
	//two pose symmetry strategy

	core::pose::Pose non_terminal_pose(pose);

  if(option[OptionKeys::remodel::repeat_structure].user()){

		//remove the extra tail from pose, but still use the pose with tail to build
		Size tail_count = repeat_tail_length_;
		while ( tail_count ){
			non_terminal_pose.conformation().delete_residue_slow(non_terminal_pose.total_residue());
			tail_count--;
		}

		Real jxnl_phi = pose.phi(non_terminal_pose.total_residue()-1);
		Real jxnl_psi = pose.psi(non_terminal_pose.total_residue()-1);
		Real jxnl_omega = pose.omega(non_terminal_pose.total_residue()-1);
		Real jxn_phi = pose.phi(non_terminal_pose.total_residue());
		Real jxn_psi = pose.psi(non_terminal_pose.total_residue());
		Real jxn_omega = pose.omega(non_terminal_pose.total_residue());
		Real jxnh_phi = pose.phi(non_terminal_pose.total_residue()+1);
		Real jxnh_psi = pose.psi(non_terminal_pose.total_residue()+1);
		Real jxnh_omega = pose.omega(non_terminal_pose.total_residue()+1);

		repeat_pose=non_terminal_pose;

		//reset foldtree and simply depend on coordinates for appendage
		core::kinematics::FoldTree ft = repeat_pose.fold_tree();
		ft.simple_tree(repeat_pose.total_residue());
		repeat_pose.fold_tree(ft);
    core::pose::remove_lower_terminus_type_from_pose_residue(non_terminal_pose, 1);

    Size repeatFactor =option[OptionKeys::remodel::repeat_structure];
    Size count = 1;

    while ( repeatFactor !=1){ // the argument should be total number of copies
      for (Size rsd = 1; rsd <= non_terminal_pose.total_residue(); rsd++){
				//intentially insert behind the last residue, this way blueprint definition will cover the junction with fragments
				if (rsd == non_terminal_pose.total_residue()){
					repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(), false);
				} else {
					repeat_pose.conformation().safely_append_polymer_residue_after_seqpos( non_terminal_pose.residue(rsd),repeat_pose.total_residue(), false);
				/*
				for (int i =1; i<= repeat_pose.total_residue(); i++){
				std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) <<  " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
				} */
				}
      }

      Size junction = (non_terminal_pose.total_residue())* count;
			repeat_pose.conformation().insert_ideal_geometry_at_polymer_bond(junction);
			/*
			//std::cout << "junction " << junction << std::endl;
			repeat_pose.set_phi(junction+1,-150);
			repeat_pose.set_psi(junction,150);
			repeat_pose.set_omega(junction,180);
			*/
			//std::cout << "junction " << junction << std::endl;
			repeat_pose.set_phi(junction-1, jxnl_phi);
			repeat_pose.set_psi(junction-1, jxnl_psi);
			repeat_pose.set_omega(junction-1, jxnl_omega);
			repeat_pose.set_phi(junction, jxn_phi);
			repeat_pose.set_psi(junction, jxn_psi);
			repeat_pose.set_omega(junction, jxn_omega);
			repeat_pose.set_phi(junction+1, jxnh_phi);
			repeat_pose.set_psi(junction+1, jxnh_psi);
			repeat_pose.set_omega(junction+1,jxnh_omega);
      //std::cout << "repeat Factor : " << repeatFactor << std::endl;
      count++;
      repeatFactor--;
    }
		//terminus residue check
		for (Size res=2; res< repeat_pose.total_residue(); res++){
			if (repeat_pose.residue(res).is_terminus()){
		//		std::cout<< "FIX TERMINUS " << res << std::endl;
				core::pose::remove_upper_terminus_type_from_pose_residue(repeat_pose, res);
				core::pose::remove_lower_terminus_type_from_pose_residue(repeat_pose, res);
			}
		}
	  if (! repeat_pose.residue(1).is_terminus()){
		  core::pose::add_variant_type_to_pose_residue( repeat_pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
	  }
	  if (! repeat_pose.residue( repeat_pose.total_residue() ).is_terminus()){
		  core::pose::add_variant_type_to_pose_residue( repeat_pose, core::chemical::UPPER_TERMINUS_VARIANT, repeat_pose.total_residue());
	  }

	using namespace protocols::loops;
	using namespace core::kinematics;

	  //set private for the class
	  repeat_length_=repeat_pose.total_residue();

	Size repeat_number =option[OptionKeys::remodel::repeat_structure];
	Size segment_length = repeat_length_/repeat_number;
	//std::cout << "DEBUG: segment lenght = " << segment_length << std::endl;

	//propagate jumps

  core::kinematics::FoldTree f;

	LoopsOP repeat_loops( new Loops() );

	std::set< Size > lower_termini;
	std::set< Size > upper_termini;
	for ( Size i = 1, ie = pose.conformation().num_chains(); i <= ie; ++i ) {
		lower_termini.insert( pose.conformation().chain_begin( i ) );
		upper_termini.insert( pose.conformation().chain_end( i ) );
	}

	for ( Loops::iterator l = loops_->v_begin(), le = loops_->v_end(); l != le; ++l ) {
    Loop & loop = *l;
		if (loop.start() == 1 && loop.stop() == segment_length+2){ // +2 for the padded shadow residues from repeat interval setup
			break;  //don't need to do anything about foldtree if fully de novo
		} else {

		repeat_loops->push_back(loop); //load the first guy
		//TR << loop.start() << " " << loop.stop() << " " << loop.cut() << std::endl;

		//math to figure out the different segments
			for (Size rep = 0; rep < repeat_number; rep++){

				//find equals iterator end means not found.  If the loops is internal, increment by repeats
				if ( lower_termini.find( loop.start() ) == lower_termini.end() && upper_termini.find( loop.stop() ) == upper_termini.end()){
					Size tempStart = loop.start()+(segment_length*rep);
					Size tempStop = loop.stop()+(segment_length*rep);
					Size tempCut = loop.cut()+(segment_length*rep);
					if (tempStop > repeat_pose.total_residue()){
							tempStop = repeat_pose.total_residue();
						if ( tempCut > tempStop ){
							tempCut = tempStop;
						}
					}
					LoopOP newLoop( new Loop( tempStart, tempStop , tempCut) );
					//TR << "adding loop in repeat propagation: " << tempStart << " " << tempStop <<  " " << tempCut << std::endl;
					repeat_loops->push_back(*newLoop);
				}
				else {
					LoopOP newLoop( new Loop(loop.start(), loop.stop(), loop.cut()) );
					repeat_loops->push_back(*newLoop);
				}
			}
		}
	}
	  if (!repeat_loops->empty()){
		  f = protocols::forge::methods::fold_tree_from_loops(repeat_pose, *repeat_loops);
	  } else {
		  f.simple_tree(repeat_pose.total_residue());
	  }
		FoldTree PFT = pose.fold_tree();
		PFT.reorder(1);
		pose.fold_tree(PFT);
		f.reorder(1);
    repeat_pose.fold_tree(f);
		TR << "REPEAT POSE FT: " <<  repeat_pose.fold_tree() << std::endl;

	for ( Loops::iterator l = repeat_loops->v_begin(), le = repeat_loops->v_end(); l != le; ++l ) {
		Loop & loop = *l;
    for (Size jxn=loop.start(); jxn <=loop.stop(); jxn++){
      //std::cout << "GEN fix junction at " << jxn << std::endl;
			if (jxn != repeat_pose.total_residue()){ // shouldn't happen because definition on repeat1, but for safety
				repeat_pose.conformation().insert_ideal_geometry_at_polymer_bond(jxn);
			}
    }
	}

		if(option[OptionKeys::remodel::no_jumps].user()){
			remove_cutpoint_variants( pose );
			remove_cutpoint_variants( repeat_pose );
			FoldTree FT;
			FT.simple_tree(repeat_pose.total_residue());
			repeat_pose.fold_tree(FT);

			FoldTree PFT;
			PFT.simple_tree(pose.total_residue());
			pose.fold_tree(PFT);

			return;
		}

		//take the jumps and set repeating RT
		//Size jump_offset = pose.num_jump(); //pose at this stage should at most have only one jump;
		//for (Size rep = 1; rep < repeat_number; rep++){
			for (Size i = 1; i<= repeat_pose.num_jump();){
				for (Size j = 1; j<= pose.num_jump(); j++,i++){
								if (i > repeat_pose.num_jump()){
									break;
								}
								FoldTree FT = repeat_pose.fold_tree();
								FT.renumber_jumps();

								numeric::xyzMatrix< Real > Rot = pose.jump(j).get_rotation();
								numeric::xyzVector< Real > Trx = pose.jump(j).get_translation();
								TR <<  "pose FT: " << pose.fold_tree() << std::endl;
								TR <<  "rpps FT: " << repeat_pose.fold_tree() << std::endl;
								TR <<  "set ROT-TRANS from " << j << " to " << i << std::endl;

								FT.set_jump_atoms( i , pose.fold_tree().jump_edge(j).upstream_atom(), pose.fold_tree().jump_edge(j).downstream_atom());
								repeat_pose.fold_tree(FT);

								Jump tempJump = repeat_pose.jump(i);
								tempJump.set_rotation( Rot );
								tempJump.set_translation( Trx );
								repeat_pose.conformation().set_jump( i, tempJump );
				}
			}

// all subsequent steps of repeat_structure gets RGF, so the initial build actually don't need it.

/*								TR <<  "preRGF " << std::endl;
		if(option[OptionKeys::remodel::helical_rise].user() &&
		   option[OptionKeys::remodel::helical_radius].user() &&
		   option[OptionKeys::remodel::helical_omega].user()){
					RGF_.restore_original_cst(repeat_pose);
					RGF_.setup_helical_constraint(repeat_pose);
		}
								TR <<  "postRGF " << std::endl;
*/
		//}
		//check tree:
		//TR << f << std::endl;
/*
		//debugging
	for ( Loops::iterator l = loops_->v_begin(), le = loops_->v_end(); l != le; ++l ) {
    Loop & loop = *l;
			for (int i = loop.start(); i<= loop.stop(); i++){
				for (Size rep = 0; rep < repeat_number; rep++){
					repeat_pose.set_phi( i+(segment_length*rep), numeric::random::rg().uniform());
					repeat_pose.set_psi( i+(segment_length*rep), numeric::random::rg().uniform());
			}
		}
	}
*/

/*

		//take care of foldtree
		core::kinematics::FoldTree f;
		f.simple_tree(repeat_pose.total_residue());
		repeat_pose.fold_tree(f);
*/
	/*
		for (int i =1; i<= repeat_pose.total_residue(); i++){
			std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) << " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
		}
  //  pose = repeat_pose;
	//	std::cout << repeat_pose.fold_tree()<< std::endl;
*/
  }
}

void RemodelLoopMover::repeat_sync( //utility function
	core::pose::Pose & repeat_pose,
	core::Size repeat_number
)
{
	using core::Size;
	using namespace protocols::loops;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;

	//repeat_pose.dump_pdb("repeatPose_in_rep_propagate.pdb");

	Size segment_length = repeat_length_/repeat_number;
	//std::cout << "DEBUG: segment lenght = " << segment_length << std::endl;

	//propagate jumps

  core::kinematics::FoldTree f;

	LoopsOP repeat_loops( new Loops() );

	utility::vector1<Size> linkPositions;

	for ( Loops::iterator l = loops_->v_begin(), le = loops_->v_end(); l != le; ++l ) {
    Loop & loop = *l;

		//collect the movable positions
		for (Size i = loop.start(); i<=loop.stop(); i++ ){
			linkPositions.push_back(i);
		}

		bool vallLookbackActive = repeat_pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA);
		VallLookbackDataOP vall_history_repeat_pose(0);
		if(vallLookbackActive){
			vall_history_repeat_pose = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >( repeat_pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
		}

		for ( Size i = 1; i<= linkPositions.size(); i++){

				Size res = linkPositions[i];
				Real loop_phi = 0;
				Real loop_psi = 0;
				Real loop_omega = 0;
				Real loop_rmsd_history = 0;
				Real loop_res_changed_history = 0;
				char loop_secstruct = 'H';
				if (res <= segment_length ){ // should already be, just to be sure

					if (res == 1){ //if the first and last positions are involved, loop around
						loop_phi = repeat_pose.phi(segment_length+1);
						loop_psi = repeat_pose.psi(segment_length+1);
						loop_omega = repeat_pose.omega(segment_length+1);
						loop_secstruct = repeat_pose.secstruct(segment_length+1);
						repeat_pose.set_phi( 1, loop_phi);
						repeat_pose.set_psi( 1, loop_psi);
						repeat_pose.set_omega( 1, loop_omega);
						repeat_pose.set_secstruct(1,loop_secstruct);
						if(vallLookbackActive){
							loop_rmsd_history = vall_history_repeat_pose->get_rmsd(1);
							loop_res_changed_history = vall_history_repeat_pose->get_res_changed(1);
						}

					} else {
						loop_phi = repeat_pose.phi(res);
						loop_psi = repeat_pose.psi(res);
						loop_omega = repeat_pose.omega(res);
						loop_secstruct = repeat_pose.secstruct(res);
						if(vallLookbackActive) {
							loop_rmsd_history = vall_history_repeat_pose->get_rmsd(res);
							loop_res_changed_history = vall_history_repeat_pose->get_res_changed(res);
						}

					}

					for (Size rep = 1; rep < repeat_number; rep++){
						repeat_pose.set_phi(res+( segment_length*rep), loop_phi );
						repeat_pose.set_psi(res+( segment_length*rep), loop_psi );
						repeat_pose.set_omega( res+(segment_length*rep),loop_omega );
						repeat_pose.set_secstruct( res+(segment_length*rep),loop_secstruct );
						if(vallLookbackActive){
							vall_history_repeat_pose->set_rmsd(res+(segment_length*rep),loop_rmsd_history);
							vall_history_repeat_pose->set_res_changed(res+(segment_length*rep),loop_res_changed_history);
						}

					}
				}
				else if (res > segment_length ){ //for spanning builds

					//spanning, update the equivalent copy in the first section with
					//new first
					repeat_pose.set_phi(res-segment_length, repeat_pose.phi(res));
					repeat_pose.set_psi(res-segment_length, repeat_pose.psi(res));
					repeat_pose.set_omega(res-segment_length, repeat_pose.omega(res));
					repeat_pose.set_secstruct(res-segment_length, repeat_pose.secstruct(res));
					//then propagate
					loop_phi = repeat_pose.phi(res-segment_length);
					loop_psi = repeat_pose.psi(res-segment_length);
					loop_omega = repeat_pose.omega(res-segment_length);
					loop_secstruct = repeat_pose.secstruct(res-segment_length);
					for (Size rep = 1; rep < repeat_number; rep++){
						if (res+( segment_length*rep)<= repeat_length_){
							repeat_pose.set_phi(res+( segment_length*rep), loop_phi );
							repeat_pose.set_psi(res+( segment_length*rep), loop_psi );
							repeat_pose.set_omega( res+(segment_length*rep),loop_omega );
							repeat_pose.set_secstruct( res+(segment_length*rep),loop_secstruct );

						}
					}

				}
				else {
					TR << "ERROR in collecting positions: res: " << res <<  std::endl;
				}
		}
	}
	if (option[OptionKeys::remodel::helical_rise].user() &&
	    option[OptionKeys::remodel::helical_radius].user() &&
	    option[OptionKeys::remodel::helical_omega].user()){
					RGF_.restore_original_cst(repeat_pose);
					RGF_.setup_CM_helical_constraint(repeat_pose);
	//				RGF_.align_segment(repeat_pose);
	}
	//repeat_pose.dump_pdb("rep_test.pdb");
}

void RemodelLoopMover::repeat_propagation( //utility function
	core::pose::Pose & pose,
	core::pose::Pose & repeat_pose,
	core::Size repeat_number
)
{
	using core::Size;
	using namespace core::chemical;
	using namespace protocols::loops;
	using namespace basic::options;
	using namespace core::kinematics;
	using namespace core::pose::symmetry;
	using namespace core::scoring::constraints;
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;



	//repeat_pose.dump_pdb("repeatPose_in_rep_propagate.pdb");

	Size segment_length = repeat_length_/repeat_number;
	//std::cout << "DEBUG: segment lenght = " << segment_length << std::endl;

	Pose junk_for_copy;
	Pose junk_for_copy_repeat;
	bool is_sym = false;

	//if symmetric, try use master copy only
	if (core::pose::symmetry::is_symmetric(pose) && core::pose::symmetry::is_symmetric(repeat_pose)){

		// save the constraints, as it would be lost
		ConstraintSetOP pose_cst_set( new ConstraintSet( *pose.constraint_set() ) );
		ConstraintSetOP repeat_pose_cst_set( new ConstraintSet( *repeat_pose_.constraint_set() ) );

		extract_asymmetric_unit( pose, junk_for_copy, false);
		extract_asymmetric_unit( repeat_pose, junk_for_copy_repeat, false);
		is_sym = true;
		pose=junk_for_copy;
		repeat_pose=junk_for_copy_repeat;

		//reset constraints
		pose.constraint_set(pose_cst_set);
		repeat_pose.constraint_set(repeat_pose_cst_set);

		pose.pdb_info()->obsolete(true);
		repeat_pose.pdb_info()->obsolete(true);
	}

	//propagate jumps

  core::kinematics::FoldTree f;

	bool build_across_jxn = false;
	Size residues_beyond_jxn = 0;

	LoopsOP repeat_loops( new Loops() );
  std::set< Size > lower_termini;
  std::set< Size > upper_termini;
  for ( Size i = 1, ie = pose.conformation().num_chains(); i <= ie; ++i ) {
    lower_termini.insert( pose.conformation().chain_begin( i ) );
    upper_termini.insert( pose.conformation().chain_end( i ) );
  }

  for ( Loops::iterator l = loops_->v_begin(), le = loops_->v_end(); l != le; ++l ) {
    Loop & loop = *l;
    if (loop.start() == 1 && loop.stop() == segment_length+2){ // padded 2 shadow residues
      break;  //don't need to do anything about foldtree if fully de novo
    } else {

			// if any of the loops build go beyond the junction, keep track of it and
			// update the first segment accordingly
			if (loop.start() <= segment_length && loop.stop() > segment_length){
				//std::cout << "build across jxn: " << loop.start() << " " << loop.stop() << std::endl;
				build_across_jxn = true;
				residues_beyond_jxn = loop.stop() - segment_length;
				//std::cout << "build across jxn leftover: " << residues_beyond_jxn << std::endl;
			}

    repeat_loops->push_back(loop); //load the first guy
    //TR << loop.start() << " " << loop.stop() << " " << loop.cut() << std::endl;

    //math to figure out the different segments
      for (Size rep = 0; rep < repeat_number; rep++){

        //find equals iterator end means not found.  If the loops is internal, increment by repeats
        if ( lower_termini.find( loop.start() ) == lower_termini.end() && upper_termini.find( loop.stop() ) == upper_termini.end()){
          Size tempStart = loop.start()+(segment_length*rep);
          Size tempStop = loop.stop()+(segment_length*rep);
          Size tempCut = loop.cut()+(segment_length*rep);
          if (tempStop > repeat_pose.total_residue()){
              tempStop = repeat_pose.total_residue();
						if ( tempCut > tempStop ){
							tempCut = tempStop;
						}
          }
          LoopOP newLoop( new Loop( tempStart, tempStop , tempCut) );
          //TR << "adding loop in repeat propagation: " << tempStart << " " << tempStop <<  " " << tempCut << std::endl;
          repeat_loops->push_back(*newLoop);
        }
        else {
          LoopOP newLoop( new Loop(loop.start(), loop.stop(), loop.cut()) );
          repeat_loops->push_back(*newLoop);
        }
      }
    }
  }
/*
	for ( Loops::iterator l = loops_->v_begin(), le = loops_->v_end(); l != le; ++l ) {
    Loop & loop = *l;
		if (loop.start() == 1 && loop.stop() == segment_length){
			break;  //don't need to do anything about foldtree if fully de novo
		} else {

		repeat_loops->push_back(loop); //load the first guy
	//	TR << loop.start() << " " << loop.stop() << " " << loop.cut() << std::endl;

		//math to figure out the different segments
			for (Size rep = 0; rep < repeat_number; rep++){
				LoopOP newLoop = new Loop(loop.start()+(segment_length*rep), loop.stop()+(segment_length*rep),loop.cut()+(segment_length*rep));
	//			TR << "adding loop in repeat propagation: " << loop.start()+(segment_length*rep) << std::endl;
				repeat_loops->push_back(*newLoop);
			}
		}
	}
*/
	if (!repeat_loops->empty() && !option[OptionKeys::remodel::no_jumps].user()){
		f = protocols::forge::methods::fold_tree_from_loops(repeat_pose, *repeat_loops);
	} else {
		f.simple_tree(repeat_pose.total_residue());
	}

    repeat_pose.fold_tree(f);

	for ( Loops::iterator l = repeat_loops->v_begin(), le = repeat_loops->v_end(); l != le; ++l ) {
		Loop & loop = *l;
    for (Size jxn=loop.start(); jxn <=loop.stop(); jxn++){
      //std::cout << "fix junction at " << jxn << std::endl;
			if (jxn != repeat_pose.total_residue()){
				repeat_pose.conformation().insert_ideal_geometry_at_polymer_bond(jxn);
			}
    }
	}

	//repeat_pose.dump_pdb("repeatPose_in_rep_propagate_setTree.pdb");
/* // repeat_generation sets up the foldtree correctly, even changin cutpoint
 * in foldtree doesn't require jump update
		//take the jumps and set repeating RT
		Size jump_offset = pose.num_jump();
		for (Size rep = 1; rep < repeat_number; rep++){
			for (Size i = 1; i<= pose.num_jump(); i++){
				numeric::xyzMatrix< Real > Rot = pose.jump(i).get_rotation();
				numeric::xyzVector< Real > Trx = pose.jump(i).get_translation();
			//	TR <<  "set ROT-TRANS from " << i << " to " << i+jump_offset*rep << std::endl;
				Jump tempJump = repeat_pose.jump(i+(jump_offset*rep));
				tempJump.set_rotation( Rot );
				tempJump.set_translation( Trx );
				repeat_pose.conformation().set_jump( i+(jump_offset*rep), tempJump );
			}
		}
		//check tree:
		//TR << f << std::endl;
*/
	//repeat_pose.dump_pdb("repeatPose_in_rep_propagate_setJump.pdb");

	//take care of the start if build across jxn
	if (build_across_jxn){
		while (residues_beyond_jxn){
			pose.set_secstruct(residues_beyond_jxn, pose.secstruct(residues_beyond_jxn+segment_length));
			pose.set_phi(residues_beyond_jxn, pose.phi(residues_beyond_jxn+segment_length));
			pose.set_psi(residues_beyond_jxn, pose.psi(residues_beyond_jxn+segment_length));
			pose.set_omega(residues_beyond_jxn, pose.omega(residues_beyond_jxn+segment_length));
			residues_beyond_jxn--;
		}
	}
	bool vallLookbackActive = pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA);
	VallLookbackDataOP vall_history_repeat_pose_(0);
	VallLookbackDataOP vall_history_pose_(0);
	if(vallLookbackActive){
		vall_history_pose_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >( pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
		vall_history_repeat_pose_ = utility::pointer::static_pointer_cast<core::scoring::methods::VallLookbackData >( repeat_pose.data().get_ptr( CacheableDataType::VALL_LOOKBACK_DATA ));
	}
	for (Size rep = 0; rep < repeat_number; rep++){
		for (Size res = 1; res <= segment_length; res++){
				//std::cout << "DEBUG: res+segmentlength*rep = " << res+(segment_length*rep) << std::endl;
				Real loop_phi = 0;
				Real loop_psi = 0;
				Size rsd_type_position = res;
				if (res == 1 ){
					rsd_type_position = segment_length+1;
					loop_phi = pose.phi(segment_length+1);
					loop_psi = pose.psi(segment_length+1);
				} else {
					loop_phi = pose.phi(res);
					loop_psi = pose.psi(res);
				}
				ResidueType const & rsd_type(pose.residue_type(rsd_type_position));
				replace_pose_residue_copying_existing_coordinates(repeat_pose,res+(segment_length*rep),rsd_type);
				repeat_pose.set_phi(res+( segment_length*rep), loop_phi );
				repeat_pose.set_psi(res+( segment_length*rep), loop_psi );
				repeat_pose.set_omega( res+(segment_length*rep), pose.omega(res) );
				repeat_pose.set_secstruct( res+(segment_length*rep),pose.secstruct(res) );
				if(vallLookbackActive){
					Real rmsd_tmp = vall_history_pose_->get_rmsd(res);
					bool res_changed_tmp= vall_history_pose_->get_res_changed(res);
					vall_history_repeat_pose_->set_rmsd(res+( segment_length*rep),rmsd_tmp);
					vall_history_repeat_pose_->set_res_changed(res+( segment_length*rep),res_changed_tmp);
				}
		}
	}

	if(option[OptionKeys::remodel::helical_rise].user() &&
			option[OptionKeys::remodel::helical_radius].user() &&
			option[OptionKeys::remodel::helical_omega].user()){
					RGF_.restore_original_cst(repeat_pose);
					RGF_.setup_CM_helical_constraint(repeat_pose);
					//RGF_.align_segment(repeat_pose);
	}
//repeat_pose.dump_pdb("post_align.pdb");
	//re-symmetrize
	if (is_sym){
		// save the constraints, as it would be lost
		ConstraintSetOP pose_cst_set( new ConstraintSet( *pose.constraint_set() ) );
		ConstraintSetOP repeat_pose_cst_set( new ConstraintSet( *repeat_pose_.constraint_set() ) );

	//unfortunately need separate movers for the operations
		simple_moves::symmetry::SetupForSymmetryMover pre_mover1;
		pre_mover1.apply( pose );
		simple_moves::symmetry::SetupForSymmetryMover pre_mover2;
		pre_mover2.apply( repeat_pose );

		//reset constraints
		pose.constraint_set(pose_cst_set);
		repeat_pose.constraint_set(repeat_pose_cst_set);

		pose.pdb_info()->obsolete(true);
		repeat_pose.pdb_info()->obsolete(true);
	}


//repeat_pose.dump_pdb("symmetrize.pdb");
//exit(0);
	//repeat_pose.dump_pdb("repeatPose_in_rep_propagate_setAngle.pdb");
	//loop over the tail fragment to the first fragment


  //pose = repeat_pose;
  //pose.dump_pdb("test_repeat1.pdb");
  //repeat_pose.dump_pdb("repeat2.pdb");
/*
		for (int i =1; i<= repeat_pose.total_residue(); i++){
			std::cout << "repeat_pose Phi: "<< repeat_pose.phi(i) << " psi: " << repeat_pose.psi(i) << " omega: " << repeat_pose.omega(i) << " at " << i << std::endl;
		} */
}


///
/// @begin RemodelLoopMover::apply
///
/// @remarks Sets protocols::moves::MS_SUCCESS upon successful closure of all loops, otherwise sets protocols::moves::FAIL_RETRY.
///
void RemodelLoopMover::apply( Pose & pose ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace OptionKeys::indexed_structure_store;
	using namespace core;
	using namespace core::scoring::methods;
	using namespace core::scoring;
	using namespace chemical;
	using core::kinematics::FoldTree;
	using protocols::forge::methods::fold_tree_from_pose;
	// archive
	FoldTree const archive_ft = pose.fold_tree();
	//std::cout << "archived foldtree " << archive_ft << std::endl;


		//initialize values that would be lost once the pose is turn symmetrical
		unit_length_ = pose.total_residue();

	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_generation_with_additional_residue(pose, repeat_pose_);

		//initialize values that would be lost once the pose is turn symmetrical
		repeat_length_ = repeat_pose_.total_residue();
		pose.pdb_info()->obsolete(true);
			repeat_pose_.pdb_info()->obsolete(true);
		if (option[OptionKeys::symmetry::symmetry_definition].user() ) {
			simple_moves::symmetry::SetupForSymmetryMover pre_mover1;
			pre_mover1.apply( pose );
			simple_moves::symmetry::SetupForSymmetryMover pre_mover2;
			pre_mover2.apply( repeat_pose_ );
			pose.pdb_info()->obsolete(true);
			repeat_pose_.pdb_info()->obsolete(true);
		}


/* not fully implemented yet, developmental
 if (option[OptionKeys::symmetry::symmetry_definition].user() )  {
    simple_moves::symmetry::SetupForSymmetryMover pre_mover;

		//try making both symmetrical
    //pre_mover.apply( repeat_pose_ );
    pre_mover.apply( pose );
    // Remodel assumes chain ID is ' '
    //pose::PDBInfoOP pdb_info ( repeat_pose_.pdb_info() );
    pose::PDBInfoOP pdb_info ( pose.pdb_info() );
    for ( Size i=1; i<= pdb_info->nres(); ++i ){
      pdb_info->chain(i,' ');
    }
    pose.pdb_info( pdb_info );
  }
*/

		// take care of constraints:

		// user defined constraint from file, make sure to do this first as it
		// replaces cst object in pose and that wipes out everything already set.

		// only use this type of cst file in this case
		if (option[OptionKeys::constraints::cst_file].user() ) {
			protocols::simple_moves::ConstraintSetMoverOP repeat_constraint( new protocols::simple_moves::ConstraintSetMover() );
			repeat_constraint->apply( repeat_pose_ );
		}

/* // ResidueTypeLinkingConstraints
		Size repeat_number =option[OptionKeys::remodel::repeat_structure];
		Real bonus = 10;
		//std::cout << "RESIDUETYPELINKING CST" << std::endl;
		Size segment_length = ( repeat_pose_.n_residue() ) / repeat_number;
		for ( Size rep = 1; rep < repeat_number; rep++ ) { // from 1 since first segment don't need self-linking
			for ( Size res = 1; res <= segment_length; res++ ) {
				repeat_pose_.add_constraint( new scoring::constraints::ResidueTypeLinkingConstraint( repeat_pose_, res, res+(segment_length*rep), bonus ) );
				//std::cout << res << " " << res+(segment_length*rep) << std::endl;
			}
		}
*/
			if(option[OptionKeys::remodel::helical_rise].user() &&
			   option[OptionKeys::remodel::helical_radius].user() &&
			   option[OptionKeys::remodel::helical_omega].user()){
							RGF_.set_native_cst_set( repeat_pose_ );
				TR.Debug << "repeat_pose_ length: " << repeat_length_ << std::endl;
							RGF_.set_segment_size( repeat_length_/option[OptionKeys::remodel::repeat_structure] );
							//RGF_.setup_CM_helical_constraint( repeat_pose_ );
			}
	}

	FoldTree sealed_ft;
	if ( pose::symmetry::is_symmetric(pose) ) {
		sealed_ft = pose::symmetry::sealed_symmetric_fold_tree( pose );
	} else if ( seal_foldtree_ ) {
		sealed_ft = fold_tree_from_pose( pose, pose.fold_tree().root(), MoveMap() ); // used during structure accumulation
	} else {
		sealed_ft = pose.fold_tree();
	}


/*
			// ResidueTypeLinkingConstraints
			Size repeat_number =option[OptionKeys::remodel::repeat_structure];
			Real bonus = 10;
			//std::cout << "RESIDUETYPELINKING CST" << std::endl;
		  Size segment_length = (repeat_pose_.n_residue())/repeat_number;
		  for (Size rep = 1; rep < repeat_number; rep++ ){ // from 1 since first segment don't need self-linking
			  for (Size res = 1; res <= segment_length; res++){
					 repeat_pose_.add_constraint( new ResidueTypeLinkingConstraint(repeat_pose_, res, res+(segment_length*rep), bonus));
	//				std::cout << res << " " << res+(segment_length*rep) << std::endl;
			  }
		  }

		std::stringstream templateRangeSS;
		templateRangeSS << "1-" << segment_length;

		// Dihedral (NCS) Constraints
		//std::cout << "NCS CST" << std::endl;
		protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
		for (Size rep = 1; rep < repeat_number; rep++){ // from 1 since first segment don't need self-linking
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep) << "-" << segment_length + (segment_length*rep);
			//std::cout << templateRangeSS.str() << " " << targetSS.str() << std::endl;

			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}
*/

/*   If want cyclize peptide in frag insertion, uncomment here.  But in
 *   practice, it seems to screw up more than helping.  So currently do
 *   fragment insertion without cyclizing pose, and only enforce this in
 *   refinement stage

		if(option[OptionKeys::RemodelLoopMover::cyclic_peptide].user()){
			protocols::forge::methods::cyclize_pose(repeat_pose_);
		}
*/


	// for accumulation of closed structures (only return the best)
	std::multimap< Real, PoseOP > accumulator;
// within simultaneous and independent stages
	ScoreFunctionOP sfxOP = sfx_;
	sfxOP->set_weight( scoring::linear_chainbreak, 0.0 );
	if(option[OptionKeys::remodel::staged_sampling::staged_sampling].user()){
		//initialize options
		std::cout << "************************INSIDE STAGED SAMPLING" << std::endl;
		if(option[OptionKeys::remodel::staged_sampling::starting_sequence].user())
			set_starting_sequence(pose);
		if(option[OptionKeys::remodel::staged_sampling::start_w_ideal_helices].user())
				set_ideal_helices(pose);
		if(option[OptionKeys::remodel::staged_sampling::starting_pdb].user())
			set_starting_pdb(pose);
		assert(!((option[OptionKeys::remodel::staged_sampling::starting_pdb].user()) && (option[OptionKeys::remodel::staged_sampling::start_w_ideal_helices].user()))); // starting pdb not compatible with ideal helices


		//setup score functions and movemap and where to sample--------------
		MoveMap movemap;
		MoveMap movemapAll;
		std::string ss = remodel_data_.ss;;
		Size singleRepeat = (pose.total_residue()/2);
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			movemap.set_bb( i, true );
			movemap.set_chi( i, true );
			movemapAll.set_bb(i, true);
			movemapAll.set_chi(i, true);
		  if(option[OptionKeys::remodel::staged_sampling::sample_over_loops].user()){
				Size ssIndex = i-1;
				if(ssIndex>=singleRepeat)
						ssIndex = i-singleRepeat-1;
				if(ss[ssIndex] == 'H' || ss[ssIndex] == 'E'){
						movemap.set_bb(i,false);
						movemap.set_chi(i,false);
						}
				}
		}
		ScoreFunctionOP sfxStage0_OP =  ( core::scoring::ScoreFunctionFactory::create_score_function( "abinitio_remodel_cen" ) );
		ScoreFunctionOP sfxStage1_OP =  ( core::scoring::ScoreFunctionFactory::create_score_function( "abinitio_remodel_cen" ) );
		sfxStage0_OP->set_weight(scoring::atom_pair_constraint, 0.0);
		sfxStage1_OP->set_weight(scoring::atom_pair_constraint, 1.0);
		if(option[OptionKeys::remodel::repeat_structure].user()){
			sfxStage1_OP->set_weight(scoring::atom_pair_constraint, 1.0 *option[OptionKeys::remodel::repeat_structure]);

		}
		sfxStage1_OP->show_pretty(TR);
		if(option[fragment_threshold_distance].user()){
			//prime the score function because of const issues involved with the score function.
			VallLookbackPotential const & potential_( ScoringManager::get_instance()->get_vallLookbackPotential());
			potential_.lookback(pose);
			potential_.lookback(repeat_pose_);
		}
		//setup fragments so they sample correctly-----------
		Real fragScoreThreshold = 0.99999;  //1.00XX indicates 1 ABEGO or HLE mismatch.  I chose to use the numbers for future finer control
		if(!option[OptionKeys::remodel::staged_sampling::require_frags_match_blueprint])
				fragScoreThreshold = 999.0;
		//setup locations to sample--------------------------
		std::set<Size> sampleAllResidues = generate_residues_to_sample(false,pose,9);
		std::set<Size> sampleSubsetResidues = generate_residues_to_sample(true,pose,9);
		//Initialize with any full length fragments-------------------------------
		if(option[OptionKeys::remodel::use_same_length_fragments])
				//999 allows for frags > 9 resiudes
				abinitio_stage( pose,999, movemap,sfxStage1_OP,1,100,sampleAllResidues,true,"full_length_frags",false,fragScoreThreshold);
		//Sample with 9mers in all positions------------------------------
		//This should be read in from staging file.
		abinitio_stage( pose, 9, movemap,sfxStage0_OP,1,100,sampleSubsetResidues ,true,"9mers_loops",false,fragScoreThreshold);
		abinitio_stage( pose, 3, movemap,sfxStage0_OP,1,100,sampleSubsetResidues ,true,"3mers_loops",false,fragScoreThreshold);
		abinitio_stage( pose, 9, movemapAll,sfxStage1_OP,3,500,sampleAllResidues,true,"9mers_allPos",false,fragScoreThreshold);
		abinitio_stage( pose, 3, movemapAll,sfxStage1_OP,3,500,sampleAllResidues,true,"3mers_allPos",false,fragScoreThreshold);
		if(option[OptionKeys::remodel::staged_sampling::fa_relax_moves].user()){
				fa_relax_stage(pose);
		}
		//cleanup to integrate with Possu------------------------------------
		PoseOP pose_prime( new Pose( pose ) );
		pose_prime->fold_tree( sealed_ft );
		// this is for scoring in repeat context
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
			(*sfxOP)(repeat_pose_);
			//this has to be set because when copying to pose_prime, it lost the
			//first phi angle
			pose_prime->set_phi(1, repeat_pose_.phi(1));
			accumulator.insert( std::make_pair( repeat_pose_.energies().total_energy(), pose_prime ) );
		} else {
			(*sfxOP)( *pose_prime );
			accumulator.insert( std::make_pair( pose_prime->energies().total_energy(), pose_prime ) );
		}
	}
	else {
		// setup parameters -- linearly scale the chainbreak weight
		Real const final_standard_cbreak_weight = 5.0;
		Real const cbreak_increment = final_standard_cbreak_weight / total_standard_cycles();

		// currently no scaling during boost_closure
		Real const final_boost_closure_cbreak_weight = 5.0;
		Real const boost_closure_cbreak_increment = ( final_boost_closure_cbreak_weight - final_standard_cbreak_weight ) / boost_closure_cycles();

		assert( final_boost_closure_cbreak_weight >= final_standard_cbreak_weight );

		sfxOP->set_weight( scoring::linear_chainbreak, 0.0 );

		//REPEAT TEST
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			sfxOP->set_weight( scoring::atom_pair_constraint, 1.0 *option[OptionKeys::remodel::repeat_structure] );
			sfxOP->set_weight( scoring::coordinate_constraint, 1.0 *option[OptionKeys::remodel::repeat_structure] );
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
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			(*sfxOP)(repeat_pose_);
			mc.score_function().show(TR, repeat_pose_);
			TR << std::endl;
			//std::cout << "check_constraint" << std::endl;
			//repeat_pose_.constraint_set()->show_definition(TR, repeat_pose_);
			//TR << std::endl;
			mc.reset(repeat_pose_);
		}
		else {
			mc.reset(pose);
		}
		for ( Size attempt = 1; attempt <= allowed_closure_attempts_; ++attempt ) {

			TR << "* closure_attempt " << attempt << std::endl;

			// reset score function at the beginning of each attempt
			mc.score_function( *sfxOP );
			if(option[OptionKeys::remodel::RemodelLoopMover::use_loop_hash].user()) {
				loophash_stage(pose, mc, cbreak_increment );
			}

			// simultaneous loop movements using fragment + ccd_move (default 20% of the time)
			simultaneous_stage( pose, mc, cbreak_increment );

			// closure is hard, so attempt closure for each loop independently using
			// fragment + ccd_move (default 80% of the time)
			independent_stage( pose, mc, cbreak_increment );

			// "boost": if any loops are not closed but within a chainbreak interval, attempt
			// to close them with 1-mer + ccd_move.
			boost_closure_stage( pose, mc, boost_closure_cbreak_increment );

			//pose.dump_pdb("test.pdb");

			// check to see if all loops closed, if so rescore w/out chainbreak
			// and store in accumulator
			if ( check_closure_criteria( pose ) ) {
				// make a pointer copy for storage, for REPEATs, store only the monomer pose
				PoseOP pose_prime( new Pose( pose ) );
				pose_prime->fold_tree( sealed_ft );

				// this is for scoring in repeat context
				if (option[OptionKeys::remodel::repeat_structure].user() ) {
					repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
					(*sfxOP)(repeat_pose_);

					//this has to be set because when copying to pose_prime, it lost the
					//first phi angle
					pose_prime->set_phi(1, repeat_pose_.phi(1));

					accumulator.insert( std::make_pair( repeat_pose_.energies().total_energy(), pose_prime ) );

				} else {
					(*sfxOP)( *pose_prime );
					accumulator.insert( std::make_pair( pose_prime->energies().total_energy(), pose_prime ) );
				}

				// now randomize the loops again for a new starting point
				if ( attempt < allowed_closure_attempts_ ) {
					randomize_stage( pose );
					if (option[OptionKeys::remodel::repeat_structure].user() ) {
						repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
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
				if (option[OptionKeys::remodel::repeat_structure].user() ) {
					repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
					(*sfxOP)(repeat_pose_);
					mc.reset( repeat_pose_ );
				} else{
					mc.reset( pose);
				}
			}

			//std::ostringstream ss;
			//ss << "rlm." << attempt << ".";
			//JobDistributor::get_instance()->job_outputter()->other_pose( JobDistributor::get_instance()->current_job(), pose, ss.str() );
		}
	}
	TR << "* " << accumulator.size() << " / " << allowed_closure_attempts_ << "   closed / attempts " << std::endl;

	// return the best structure if available, otherwise mark failure
	if ( !accumulator.empty() ) {
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			pose = *( accumulator.begin()->second );
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
			(*sfxOP)(repeat_pose_);
			// an issue with symmetry info in equal operator assignment, so to be safe explicitly make Pose again
			if(option[OptionKeys::symmetry::symmetry_definition].user() ){
				Pose junk_for_copy;
				core::scoring::constraints::ConstraintSetOP pose_cst_set( new core::scoring::constraints::ConstraintSet( *pose.constraint_set() ) );
				core::pose::symmetry::extract_asymmetric_unit( repeat_pose_, junk_for_copy, false);
				pose = junk_for_copy;
				pose.constraint_set(pose_cst_set);
				pose.pdb_info()->obsolete(true);
				//resymmetrize
				protocols::simple_moves::symmetry::SetupForSymmetryMover pre_mover;
				pre_mover.apply(pose);
				pose.pdb_info()->obsolete(true);
			}
			else {
				pose = repeat_pose_;
			}
		} else {
			pose = *( accumulator.begin()->second );
		}

		set_last_move_status( protocols::moves::MS_SUCCESS );

	} else {

		// if use bypass_closure flag, all failure is passed as success.  the
		// problem is that the subsequent restore_sidechain function won't know if
		// it's fail or success.  although it doesn't matter if simply
		// bypass_closure, it is a problem when we use lh_filter for testing the
		// likeliness of a backbone structure.  In this case, a failed criteria
		// should exit. otherwise the length would all be messed up.  Since masking pose
		// as built successfully, pass it on as a repeat, not a monomer in this
		// case

		//but hold on.. the filter failure should be the only reason monomer is
		//passed on.  so maybe we don't need to grow it afterall
		/*
			if(option[OptionKeys::remodel::repeat_structure].user()SACHKO.user() && SACHKO_remodel_bypass_closure_){
			pose = *( accumulator.begin()->second );
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			(*sfxOP)(repeat_pose_);
			pose = repeat_pose_;
			}
		*/
		set_last_move_status( protocols::moves::FAIL_RETRY );
		if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user() &&
				option[OptionKeys::remodel::lh_filter_string].user()){
			//activates lh_filter for plausible backbone, in this case, a failed
			//loop should exit.
			TR << "fail in loop building: EXIT " << std::endl;
			exit(0);
		}
	}

	if (option[OptionKeys::remodel::repeat_structure].user() ) {
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

///
/// @begin RemodelLoopMover::randomize_stage
///
/// @brief
/// randomize loops
///
void RemodelLoopMover::randomize_stage( Pose & pose ) {

	using namespace basic::options;
	using namespace core;
	using core::kinematics::FoldTree;

	TR << "** randomize_stage" << std::endl;

	// archive
	FoldTree const archive_ft = pose.fold_tree();

	// simul movemap -- all loops moveable
	MoveMap movemap;
	mark_loops_moveable( loops_, movemap, true );
	enforce_false_movemap( movemap );

	// set appropriate topology
	if ( keep_input_foldtree_ ) {
	}
	else {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, *loops_, f_new );
			pose.fold_tree( f_new );
		} else {
			pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, *loops_ ) );
		}
	}

	// init fragment mover for each fragment set
	FragmentMoverOPs frag_movers = create_fragment_movers( movemap );

	Size const n_moveable = count_moveable_residues( movemap, 1, unit_length_ );

	// insert random number of fragments = n_frag_movers * moveable_residues
	for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
		for ( Size j = 0; j < n_moveable; ++j ) {
			(*i)->apply( pose );
		}
	}

	check_closure_criteria( pose, true );

	//return original foldtree
	pose.fold_tree(archive_ft);
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
	loops::LoopsOP loops_to_model( new loops::Loops() );

	if ( only_broken_loops ) {
		loops_to_model = determine_loops_to_model( pose );
	} else {
		loops_to_model = loops_;
	}

	// set appropriate topology
	if ( keep_input_foldtree_ ){
	}
	else {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			core::kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, *loops_to_model, f_new );
			pose.fold_tree( f_new );
		} else {
		pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, *loops_to_model ) );
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
		if(option[OptionKeys::remodel::repeat_structure].user()){
			//Pose temp_pose(pose);
			repeat_propagation(pose, repeat_pose_, option[OptionKeys::remodel::repeat_structure]);
			//pose = temp_pose;
		}
	}
}


//check for clashes based on atom distances
//only check atoms in AtomID vector
//should be way faster than calculating entire score
bool
fast_clash_check(
  core::pose::Pose const & pose,
  utility::vector1< core::id::AtomID > const check_atids,
  core::Real const clash_dist_cut
)
{
	using namespace core;
  Real const clash_dist2_cut( clash_dist_cut * clash_dist_cut );
  for( Size iatid = 1; iatid <= check_atids.size(); ++iatid ){
    Vector const at1_xyz( pose.xyz( check_atids[ iatid ] ) );
    for( Size res2 = 1; res2 <= pose.total_residue(); ++res2 ){
      for( Size at2 = 1; at2 <= pose.residue( res2 ).natoms(); ++at2 ){
        //skip virtual atoms!
        if( pose.residue( res2 ).atom_type( at2 ).lj_wdepth() == 0.0 ) continue;
        id::AtomID atid2( at2, res2 );
        //skip if atid2 is in check_atids
        bool skip_at2( false );
        for( Size jatid = 1; jatid <= check_atids.size(); ++jatid ){
          if( atid2 == check_atids[ jatid ] ){ skip_at2 = true; break; }
        }
        if( skip_at2 ) continue;
        Real const dist2( at1_xyz.distance_squared( pose.xyz( atid2 ) ) );
        if( dist2 < clash_dist2_cut ){
          //TR_unsat << "CLASH!: " << check_atids[ iatid ] << " - " << atid2 <<
          //   " = " << dist2 << std::endl;
          return true;
        }
      }
    }
  }
  return false;
}




/// @brief independent stage: single loop movement prior to MC accept/reject
void RemodelLoopMover::loophash_stage(
	Pose & pose,
	MonteCarlo & mc,
	Real const cbreak_increment
)
{
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::add_cutpoint_variants;
	using namespace basic::options;
	using namespace OptionKeys::remodel;
	using namespace protocols::loophash;
	using namespace numeric::geometry::hashing;
	using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** LoopHash_stage" << std::endl;

	Pose const constantPose(pose); //needed this because get_rt function for loophash doesn't honor the cut positions needed to build the loop.

	// setup loops
	loops::LoopsOP loops_to_model( new loops::Loops() );

	//find terminal loops and skip them; loophash can't handle terminal loops
	for ( Loops::const_iterator l = loops_->begin(), le = loops_->end(); l != le; ++l ) {
		bool skip_loop = false;
		if ( l->is_terminal( pose ) || l->start() == 1 ) {
			skip_loop = true;
		}
		if ( !skip_loop ) {
			loops_to_model->add_loop( *l );
		}
	}

	TR << "   n_loops = " << loops_to_model->size() << std::endl;

	// if filter is used, make sure the number of strings specified agree with
	// num loops.
	utility::vector1< std::string > filter_target;
	if(option[OptionKeys::remodel::lh_filter_string].user()){
		filter_target =option[OptionKeys::remodel::lh_filter_string];
		runtime_assert(filter_target.size() == loops_to_model->size());
	}

	if ( loops_to_model->size() == 0 ) { // nothing to do...
		return;
	}

	utility::vector1<core::Size> loopsizes;

	//find the loopsize

	for ( Loops::iterator l = loops_to_model->v_begin(), le = loops_to_model->v_end(); l != le; ++l ) {
		Loop & loop = *l;
		Size db_to_use = loop.stop() - loop.start() + 2; // +2 for the strange way loophash RT is setup.
		loopsizes.push_back(db_to_use);
	}


	// parameters
	//Size const n_standard_cycles = total_standard_cycles();
	//	Size const n_standard_cycles = 3;
	Size const max_outer_cycles = loophash_cycles();
	//	Size const max_outer_cycles = 1;

	// per-loop frag + ccd_move

	Size loop_number = 1; // for lh_filter_string index

	for ( Loops::iterator l = loops_to_model->v_begin(), le = loops_to_model->v_end(); l != le; ++l, loop_number++ ) {
		Loop & loop = *l;

	utility::vector1<core::Size> local_loopsizes;
	local_loopsizes.push_back(loopsizes[loop_number]);

	//test loophashing
	LoopHashLibraryOP loop_hash_library( new LoopHashLibrary ( local_loopsizes , 1, 0 ) );

		// movemap
		MoveMap movemap;
		mark_loop_moveable( loop, movemap, true );
		enforce_false_movemap( movemap );

		// fragment movers
		FragmentMoverOPs frag_movers = create_fragment_movers( movemap );
		assert( !frag_movers.empty() );

		// parameters
		//Size const n_moveable = count_moveable_residues( movemap, loop.start(), loop.stop() );
		//currently looping over all the hashed loops
		//Size const max_inner_cycles = std::max( static_cast< Size >( 50 ), 10 * n_moveable );

		// set appropriate topology
		if ( keep_input_foldtree_ ){
		}
		else {
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				protocols::loops::set_single_loop_fold_tree( pose, loop );
			} else {
				if(option[OptionKeys::remodel::repeat_structure].user()){
					//do nothing. remake foldtree will destroy the pose when propagate across unconnected segments
				} else {
					protocols::forge::methods::set_single_loop_fold_tree( pose, loop );
				}
			}
		}

		TR << pose.fold_tree() << std::endl;

		//pose.dump_pdb("fix_junction.pdb");
		// add cutpoint variants
		add_cutpoint_variants( pose );
		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			add_cutpoint_variants( repeat_pose_ );
			mc.reset(repeat_pose_);
		}	else{
			mc.reset( pose );
		}
		// reset counters
		mc.reset_counters();

		Size loopsize = loopsizes[loop_number];

		loop_hash_library->load_mergeddb();

		//container for the actual fragments to use
		std::vector < BackboneSegment > bs_vec_;

		const BackboneDB & bbdb_ = loop_hash_library->backbone_database();

		Real6 loop_transform;

		PoseOP pose_for_rt( new Pose() );

		// for a pose carrying ligand, find the last peptide edge and only use the
		// peptide part
		using namespace core::chemical;

		Size max_res = 0;
		for (core::Size i = 1; i<= constantPose.total_residue(); i++){
			if (!constantPose.residue_type(i).is_ligand()){ //if not ligand, and assume ligand is always at the end!
				max_res = i;
			}
		}
		utility::vector1< core::Size > residue_indices;
		for(core::Size i = 1; i <= max_res; ++i){
				residue_indices.push_back(i);
		}
		ResidueTypeSetCOP residue_set(
				ChemicalManager::get_instance()->residue_type_set( CENTROID )
		);
		core::io::pdb::pose_from_pose(*pose_for_rt, constantPose, *residue_set, residue_indices);

		core::kinematics::FoldTree f;
		f.simple_tree(pose_for_rt->total_residue());
		pose_for_rt->fold_tree(f);

		//potentially throwing error if this step below doesn't pass
		//std::cout << "start " << loop.start() << " stop " << loop.stop() << std::endl;
		TR << "loophashing from " << loop.start()-1 << " to " << loop.stop()+1 << " using " << loop.stop()+1 - (loop.start()-1) << " residue loops." << std::endl;

		get_rt_over_leap_fast( *pose_for_rt, loop.start()-1, loop.stop()+1, loop_transform);

/*
		std::cout << sqrt((loop_transform[1]* loop_transform[1]) + (loop_transform[2]*loop_transform[2]) + (loop_transform[3]*loop_transform[3])) << std::endl;
		std::cout << loop_transform[1] << std::endl;
		std::cout << loop_transform[2] << std::endl;
		std::cout << loop_transform[3] << std::endl;
		std::cout << loop_transform[4] << std::endl;
		std::cout << loop_transform[5] << std::endl;
		std::cout << loop_transform[6] << std::endl;
*/


		BackboneSegment backbone_;
		LoopHashMap &hashmap = loop_hash_library->gethash(loopsize);

		Size lh_ex_limit =option[OptionKeys::remodel::lh_ex_limit];
		std::vector < core::Size > leap_index_list;

		TR << "radius = ";
		for (Size radius = 0; radius <= lh_ex_limit ; radius++ ){
			hashmap.radial_lookup( radius, loop_transform, leap_index_list );
			TR << radius << "... " << leap_index_list.size() << std::endl;
			if (leap_index_list.size() < 1000){ //making sure at least harvest one segment to build.
				continue;
			} else {
				break;
			}
		}
		TR << std::endl;

		//not doing shuffle for now
		//numeric::random::random_permutation( leap_index_list.begin(), leap_index_list.end(), numeric::random::rg() );

		Size lh_frag_count = leap_index_list.size();
		if (leap_index_list.size() == 0){
			TR.Warning << "No fragment found within radius=" << lh_ex_limit << "A.  Skip closure..." << std::endl;
			return;
		}
		else{
			TR << "Collected " << leap_index_list.size() << " fragments within radius=" << lh_ex_limit << "A." << std::endl;
		}

		// sanity check
		for( std::vector < core::Size >::const_iterator itx = leap_index_list.begin(); itx != leap_index_list.end(); ++itx ){
					core::Size bb_index = *itx;
					LeapIndex cp = hashmap.get_peptide( bb_index );
					bbdb_.get_backbone_segment( cp.index, cp.offset , loopsize , backbone_ );
					bs_vec_.push_back( backbone_ );
		}

		if( bs_vec_.size() == 0 ) {
				TR.Warning << "No fragment found for loop size=" << loopsize << ".  Skip closure..." << std::endl;
				return;
		}

		// do closure
		for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {

			// increment the chainbreak weight
			ScoreFunctionOP sfxOP = mc.score_function().clone();
			sfxOP->set_weight(
				core::scoring::linear_chainbreak,
				sfxOP->get_weight( core::scoring::linear_chainbreak ) + cbreak_increment
				//sfxOP->get_weight( core::scoring::linear_chainbreak ) + 1
			);

			if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			if(option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide].user()){
				sfxOP->set_weight(
					core::scoring::atom_pair_constraint,
					sfxOP->get_weight( core::scoring::atom_pair_constraint) + cbreak_increment
				);
			}

			if(option[OptionKeys::remodel::lh_cbreak_selection].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak,option[OptionKeys::remodel::lh_cbreak_selection]);
			}

			mc.score_function( *sfxOP );

			// recover low
			if(option[OptionKeys::remodel::repeat_structure].user()){
				Size copy_size =0;
				if(option[OptionKeys::remodel::repeat_structure] == 1){
					copy_size = unit_length_-1;
				} else {
					copy_size = unit_length_;
				}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
					pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
				}
				repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			}
			else{
				pose = mc.lowest_score_pose();
			}

			// currently it seems the collection of loops aren't too many.  So run through them all, or whatever cycle defined by n_standard_cycles
			//for ( Size inner = 1; inner <= 1; ++inner ) {

			// fragments
			if(option[OptionKeys::remodel::lh_filter_string].user() ) {
			for ( std::vector< BackboneSegment >::iterator i = bs_vec_.begin(), ie = bs_vec_.end(); i != ie; ++i) {
			//if ( loop.is_terminal( pose ) || numeric::random::rg().uniform() * n_standard_cycles > ( outer + simultaneous_cycles() ) || pose.fold_tree().num_cutpoint() == 0 ) {

					//again, no shuffling in early development stage
					//random_permutation( bs_vec_.begin(), bs_vec_.end(), numeric::random::rg() );

						std::vector<core::Real> phi = (*i).phi();
						std::vector<core::Real> psi = (*i).psi();
						std::vector<core::Real> omega = (*i).omega();
						Size seg_length = (*i).length();

						//check sec. struct at stub.
						//Size idxresStart = (int)loop.start()-1;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
						//Size idxresStop = (int)loop.start()-1+(seg_length-1);  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
						Size idxresStart = 0;  // 0 means starting from the jump position!
						Size idxresStop = seg_length-1;

						//special case for DB's test
					  if(option[OptionKeys::remodel::lh_filter_string].user()){
										core::sequence::ABEGOManager AM;
										std::string alphabet;
										std::string target = filter_target[loop_number];
										//turn string to same case
										boost::to_upper(target);

										for (Size idx = idxresStart; idx <= idxresStop; idx++){
											alphabet += AM.index2symbol( AM.torsion2index(phi[idx],psi[idx], omega[idx],1));
										}
										runtime_assert(alphabet.length() == target.length());

										// if X is contained in the filter string, don't filter
										if ( target.find('X') != std::string::npos){ // X found
										//do nothing
										} else if ( alphabet.compare( target ) != 0 && target.find('X') == std::string::npos ){ //No X in filter and if alphabet and target don't match, skip segment
											TR.Debug << "lh frag at " << idxresStart << " and " << idxresStop << " not as " << target <<  ": " <<  alphabet << std::endl;
											lh_frag_count--;
											continue;
										} else if ( alphabet.compare(target ) == 0 && target.find('X') == std::string::npos ){
										//found a match string, do nothing
										}
										else {
											TR.Debug << "logic error somewhere in lh filter string" << std::endl;
										}
						}

						for ( Size i = 0; i < seg_length; i++){
							Size ires = (int)loop.start()-1+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
							if (ires > unit_length_ ) break;
						//	std::cout << phi[i] << " " << psi[i] << " " << omega[i] << " " << ires << std::endl;
							pose.set_phi( ires, phi[i]);
							pose.set_psi( ires, psi[i]);
							pose.set_omega( ires, omega[i]);
						}

						CCDLoopClosureMover ccd_mover( loop, core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new MoveMap( movemap ) ) ) );
						ccd_mover.max_cycles( 50 );  // Used to be 10 moves, which would result in 50 "tries" in the old code. ~Labonte
						if(option[OptionKeys::remodel::repeat_structure].user()){
							for ( Size i = 0; i < seg_length; i++){
								Size ires = (int)loop.start()-1+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
								if (ires > repeat_length_ ) break;
							//	std::cout << phi[i] << " " << psi[i] << " " << omega[i] << " " << ires << std::endl;
								repeat_pose_.set_phi( ires, phi[i]);
								repeat_pose_.set_psi( ires, psi[i]);
								repeat_pose_.set_omega( ires, omega[i]);
							}
							repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);

							//pass every build through ccd for now
							if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
								ccd_mover.apply( repeat_pose_ );
							}
							repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "loop_hash-ccd");
						} else {
							for ( Size i = 0; i < seg_length; i++){
								Size ires = (int)loop.start()-1+i;  // this is terrible, due to the use of std:vector.  i has to start from 0, but positions offset by 1.
								if (ires > unit_length_ ) break;
							//	std::cout << phi[i] << " " << psi[i] << " " << omega[i] << " " << ires << std::endl;
								pose.set_phi( ires, phi[i]);
								pose.set_psi( ires, psi[i]);
								pose.set_omega( ires, omega[i]);
							}
							if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
								ccd_mover.apply( pose );
							}
							mc.boltzmann( pose, "loop_hash ccd");
						}

			} // inner_cycles hashed loops
			TR << "Sec struc filtered fragment count = " << lh_frag_count << std::endl;
			}
		} // outer_cycles

		// recover low
		if(option[OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if(option[OptionKeys::remodel::repeat_structure] == 1){
				copy_size = unit_length_-1;
			} else {
				copy_size = unit_length_;
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

		if(option[OptionKeys::remodel::repeat_structure].user()){
			mc.score_function().show( TR, repeat_pose_ );
			remove_cutpoint_variants( repeat_pose_ );
			remove_cutpoint_variants( pose );
		}	else {
			mc.score_function().show( TR, pose );
			remove_cutpoint_variants( pose );
		}

		TR << std::endl;
		mc.show_state();

	}// for each loop
	check_closure_criteria( pose, true );

}

/// @brief abinitio_stage::Assumes no loops need to be closed.
void RemodelLoopMover::abinitio_stage(
	Pose & pose,
	Size const fragmentSize,
	MoveMap const movemap,
	ScoreFunctionOP sfxOP,
	Size const max_outer_cycles,
	Size const max_inner_cycles,
	std::set<Size> const & disallowedPos,
	bool const recover_low,
	std::string stage_name,
	bool const smoothMoves,
	Real const fragScoreThreshold
	)
{
	using namespace basic::options;
	using namespace core;
	using core::kinematics::FoldTree;
	using namespace chemical;
	using namespace OptionKeys::remodel;
	using numeric::random::random_permutation;
	using namespace core::scoring;
	using protocols::loops::add_cutpoint_variants;
	using protocols::loops::remove_cutpoint_variants;
	using namespace core::pose::datacache;
	using namespace core::scoring::methods;


	TR << "** abinitio_stage_" << stage_name << std::endl;

	Real const temp = temperature_;
	MonteCarlo mc( *sfxOP, temp ); // init without pose
	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
		(*sfxOP)(repeat_pose_);
		mc.reset(repeat_pose_);
	}
	else {
		mc.reset(pose);
	}

	FragmentMoverOPs frag_movers = create_fragment_movers_limit_size(movemap, fragmentSize,disallowedPos,smoothMoves,fragScoreThreshold);
	assert( !frag_movers.empty() );

	// set appropriate topology
	//will need to change the foldtree when doing things with symmetry. For now I'm keeping the original foldtree

	// add cutpoint variants
	if (pose.num_jump() >0){
		add_cutpoint_variants( pose );
	}
	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
			if (repeat_pose_.num_jump()>0){
				add_cutpoint_variants( repeat_pose_ );
			}
		mc.reset(repeat_pose_);
	} else {
		mc.reset( pose );
	}
	// reset counters
	mc.reset_counters();

	// simul frag + ccd_move
	for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {
		// increment the chainbreak weight
		ScoreFunctionOP sfxOP = mc.score_function().clone();
		mc.score_function( *sfxOP );
		// recover low
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			Size copy_size =0;
			if (option[OptionKeys::remodel::repeat_structure] == 1 ) {
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}

			for (Size res = 1; res<=copy_size; res++){
				pose.set_phi(res,mc.lowest_score_pose().phi(res));
				pose.set_psi(res,mc.lowest_score_pose().psi(res));
				pose.set_omega(res,mc.lowest_score_pose().omega(res));
				pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
			}
		}else{
			pose = mc.lowest_score_pose();
		}
		if(pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA)){
			//Efficiency could be improved by copying rmsd and lookback. But I figure a fresh copy might be better for now
			VallLookbackPotential const & potential_( ScoringManager::get_instance()->get_vallLookbackPotential());
			potential_.lookback(pose);
		}
		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
		}
		for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
			// fragments
			random_permutation( frag_movers.begin(), frag_movers.end(), numeric::random::rg() );
			for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
				if(option[OptionKeys::remodel::repeat_structure].user()){
					(*i)->apply( repeat_pose_ );
					repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
					mc.boltzmann( repeat_pose_, stage_name );
				}else {
					(*i)->apply( pose );
					mc.boltzmann( pose, stage_name );
				}
			}
		}
		mc.show_state();
		TR<< "showing constraints on repeat pose" << std::endl;
		repeat_pose_.constraint_set()->show_definition(TR, repeat_pose_);
		mc.score_function().show( TR, repeat_pose_ );
		// recover low
		if(recover_low){
			if(option[OptionKeys::remodel::repeat_structure].user()){
				Size copy_size =0;
				if(option[OptionKeys::remodel::repeat_structure] == 1){
					copy_size = pose.total_residue() - repeat_tail_length_;
				} else {
					copy_size = pose.total_residue();
				}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
					pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
				}
			}
			else{
				pose = mc.lowest_score_pose();
			}
		}
		else{
				if(option[OptionKeys::remodel::repeat_structure].user()){
					Size copy_size =0;
					if(option[OptionKeys::remodel::repeat_structure] == 1){
						copy_size = pose.total_residue() - repeat_tail_length_;
					} else {
						copy_size = pose.total_residue();
					}
					for (Size res = 1; res<=copy_size; res++){
						pose.set_phi(res,repeat_pose_.phi(res));
						pose.set_psi(res,repeat_pose_.psi(res));
						pose.set_omega(res,repeat_pose_.omega(res));
						pose.set_secstruct(res,repeat_pose_.secstruct(res));
					}
				}
		}
	}
	if(pose.data().has( CacheableDataType::VALL_LOOKBACK_DATA)){
		//Efficiency could be improved by copying rmsd and lookback. But I figure a fresh copy might be better for now
		VallLookbackPotential const & potential_( ScoringManager::get_instance()->get_vallLookbackPotential());
		potential_.lookback(pose);
		potential_.lookback(repeat_pose_);
	}
}

/// @brief  relax stage
 void RemodelLoopMover::fa_relax_stage(
				 Pose & pose
				 ){
				using namespace basic::options;
				using namespace core;
				using namespace chemical;
				using namespace core::scoring;
				using namespace core::pose::symmetry;
				core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(TALARIS_2013 ));
				//scorefxn->set_weight(core::scoring::dihedral_constraint, 10.0 );
				//scorefxn->set_weight(core::scoring::coordinate_constraint, 5.0);
				Pose fa_pose = pose;
				if (option[OptionKeys::remodel::repeat_structure].user() ) {
						repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
						fa_pose = repeat_pose_;
						//protocols::simple_moves::symmetry::SetupNCSMover setup_ncs = generate_ncs_csts(pose);
						//setup_ncs.apply(fa_pose);
				}
				core::util::switch_to_residue_type_set( fa_pose, core::chemical::FA_STANDARD);
				protocols::relax::FastRelax frelax(scorefxn,1);//only 1 stage
				TR << "Relaxing pose" << std::endl;
				frelax.apply(fa_pose);
				TR << "Finished relaxing" << std::endl;
				core::util::switch_to_residue_type_set( fa_pose, core::chemical::CENTROID);
				//copy fa_pose to original_pose
				if (option[OptionKeys::remodel::repeat_structure].user() ) {
						Size copy_size =0;
						if (option[OptionKeys::remodel::repeat_structure] == 1 ) {
								copy_size = pose.total_residue()-1;
						} else {
								copy_size = pose.total_residue();
						}
						for (Size res = 1; res<=copy_size; res++){
								pose.set_phi(res,fa_pose.phi(res));
								pose.set_psi(res,fa_pose.psi(res));
								pose.set_omega(res,fa_pose.omega(res));
								pose.set_secstruct(res,fa_pose.secstruct(res));
						}
				}
				else{//non repeat case
						pose = fa_pose;
				}
 }

protocols::simple_moves::symmetry::SetupNCSMover RemodelLoopMover::generate_ncs_csts(Pose & pose){
		using namespace basic::options;
		using namespace core;
		using namespace core::pose::symmetry;
		using namespace protocols;
		protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
		Size asym_length = pose.total_residue();
		Size repeat_number =option[OptionKeys::remodel::repeat_structure];
		Size segment_length = asym_length/2;
		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		std::stringstream templateRangeSS;
		// take care of the terminal repeat, since the numbers are offset.
		templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
		std::stringstream targetSS;
		targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
		TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
		setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		return(setup_ncs);
}


/// @brief small_move_stage::Assumes no loops need to be closed.
void RemodelLoopMover::small_move_stage(
		Pose & pose,
		MoveMap const movemap,
		ScoreFunctionOP sfxOP,
		Size const max_outer_cycles,
		Size const max_inner_cycles,
		bool const recover_low,
		Real const h_range,
		Real const e_range,
		Real const l_range)
{
	using namespace basic::options;
	using namespace core;
	using core::kinematics::FoldTree;
	using namespace chemical;
	using namespace OptionKeys::remodel;
	using numeric::random::random_permutation;

	using protocols::loops::add_cutpoint_variants;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** small_move_stage_" << std::endl;

	Real temp = temperature_;
	MonteCarlo mc( *sfxOP, temp ); // init without pose
	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
		(*sfxOP)(repeat_pose_);
		mc.reset(repeat_pose_);
	}
	else {
		mc.reset(pose);
	}

	// add cutpoint variants
	if (pose.num_jump() >0){
		add_cutpoint_variants( pose );
	}
	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
			if (repeat_pose_.num_jump()>0){
				add_cutpoint_variants( repeat_pose_ );
			}
		mc.reset(repeat_pose_);
	} else {
		mc.reset( pose );
	}

	// reset counters
	mc.reset_counters();

		Size nmoves = 1;

		core::kinematics::MoveMapOP mm_temp( new core::kinematics::MoveMap( movemap ) );
simple_moves::SmallMoverOP small_mover( new simple_moves::SmallMover( mm_temp, temp, nmoves) );
		small_mover->angle_max( 'H', h_range );
		small_mover->angle_max( 'E', e_range );
		small_mover->angle_max( 'L', l_range );


	for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {
		// increment the chainbreak weight
		ScoreFunctionOP sfxOP = mc.score_function().clone();
		mc.score_function( *sfxOP );
		// recover low
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			Size copy_size =0;
			if (option[OptionKeys::remodel::repeat_structure] == 1 ) {
				copy_size = pose.total_residue()-1;
			} else {
				copy_size = pose.total_residue();
			}

			for (Size res = 1; res<=copy_size; res++){
				pose.set_phi(res,mc.lowest_score_pose().phi(res));
				pose.set_psi(res,mc.lowest_score_pose().psi(res));
				pose.set_omega(res,mc.lowest_score_pose().omega(res));
				pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
			}
		}else{
			pose = mc.lowest_score_pose();
		}
		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
		}
		for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
			// fragments
				if(option[OptionKeys::remodel::repeat_structure].user()){
					small_mover->apply( repeat_pose_ );
					repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
					mc.boltzmann( repeat_pose_, "small_moves" );
				}else {
					small_mover->apply( pose );
					mc.boltzmann( pose, "small_moves" );
				}
			}
		}
		mc.show_state();
		TR<< "showing constraints on repeat pose" << std::endl;
		repeat_pose_.constraint_set()->show_definition(TR, repeat_pose_);
		mc.score_function().show( TR, repeat_pose_ );
		// recover low
		if(recover_low){
			if(option[OptionKeys::remodel::repeat_structure].user()){
				Size copy_size =0;
				if(option[OptionKeys::remodel::repeat_structure] == 1){
					copy_size = pose.total_residue() - repeat_tail_length_;
				} else {
					copy_size = pose.total_residue();
				}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
					pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
				}
			}
			else{
				pose = mc.lowest_score_pose();
			}
		}
		else{
				if(option[OptionKeys::remodel::repeat_structure].user()){
					Size copy_size =0;
					if(option[OptionKeys::remodel::repeat_structure] == 1){
						copy_size = pose.total_residue() - repeat_tail_length_;
					} else {
						copy_size = pose.total_residue();
					}
					for (Size res = 1; res<=copy_size; res++){
						pose.set_phi(res,repeat_pose_.phi(res));
						pose.set_psi(res,repeat_pose_.psi(res));
						pose.set_omega(res,repeat_pose_.omega(res));
						pose.set_secstruct(res,repeat_pose_.secstruct(res));
					}
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
	using namespace basic::options;
	using namespace core;
	using core::kinematics::FoldTree;

	using namespace OptionKeys::remodel;
	using protocols::loops::add_cutpoint_variants;
	using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
	using protocols::loops::remove_cutpoint_variants;
	using numeric::random::random_permutation;

	TR << "** simultaneous_stage" << std::endl;
	loops::LoopsOP loops_to_model( new loops::Loops(*loops_) );
	TR << "   n_loops = " << loops_to_model->size() << std::endl;
	if ( loops_to_model->size() == 0 ) { // nothing to do...
		return;
	}
	//TR << "starting foldtree in SIMU stage " << pose.fold_tree() << std::endl;

	// Create fragment movers for each loop for each fragment set.  We want to allow an equal probability of movement per-loop, rather than per-residue.
	FragmentMoverOPs frag_movers = create_per_loop_fragment_movers( loops_to_model );
	assert( !frag_movers.empty() );

	// set appropriate topology
	if ( keep_input_foldtree_ ){
	}
	else {
		if ( pose::symmetry::is_symmetric( pose ) ) {
			kinematics::FoldTree f_new;
			protocols::loops::fold_tree_from_loops( pose, *loops_to_model, f_new );
			pose.fold_tree( f_new );
		} else {
			pose.fold_tree( protocols::forge::methods::fold_tree_from_loops( pose, *loops_to_model ) );
		}
	}

	// add cutpoint variants
	if (pose.num_jump() >0){
		add_cutpoint_variants( pose );
	}
	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure] );
			if (repeat_pose_.num_jump()>0){
				add_cutpoint_variants( repeat_pose_ );
			}
		mc.reset(repeat_pose_);
	} else {
		mc.reset( pose );
	}

	// setup master movemap covering all loops -- used only for tracking purposes
	MoveMap movemap;
	mark_loops_moveable( loops_to_model, movemap, true );
	enforce_false_movemap( movemap );

	// parameters
	Size const n_moveable = count_moveable_residues( movemap, 1, unit_length_ );
	Size const n_standard_cycles = total_standard_cycles();
	Size const max_outer_cycles = simultaneous_cycles();
	Size const max_inner_cycles = std::max( 50 * loops_to_model->size(), 10 * n_moveable );
	bool apply_user_provided_movers( user_provided_movers_.size() != 0 );

	// reset counters
	mc.reset_counters();

	// simul frag + ccd_move
	for ( Size outer = 1; outer <= max_outer_cycles; ++outer ) {

		// increment the chainbreak weight
		ScoreFunctionOP sfxOP = mc.score_function().clone();
		sfxOP->set_weight( scoring::linear_chainbreak, sfxOP->get_weight( scoring::linear_chainbreak ) + cbreak_increment );
		if (option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user() ) {
			sfxOP->set_weight( scoring::linear_chainbreak, 0 );
		}
		if(option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide].user()){
		//	sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			sfxOP->set_weight( core::scoring::atom_pair_constraint, 0);//ramping from 0
		}

		mc.score_function( *sfxOP );

		// recover low
		if (option[OptionKeys::remodel::repeat_structure].user() ) {
			Size copy_size =0;
			if (option[OptionKeys::remodel::repeat_structure] == 1 ) {
				copy_size = unit_length_-1;
			} else {
				copy_size = unit_length_;
			}

			for (Size res = 1; res<=copy_size; res++){
				pose.set_phi(res,mc.lowest_score_pose().phi(res));
				pose.set_psi(res,mc.lowest_score_pose().psi(res));
				pose.set_omega(res,mc.lowest_score_pose().omega(res));
				pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
			}
		}else{
			pose = mc.lowest_score_pose();
		}

		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation( pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
		}

		for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {

			if ( numeric::random::rg().uniform() * n_standard_cycles > outer || pose.fold_tree().num_cutpoint() == 0 ) {
				// fragments
				random_permutation( frag_movers.begin(), frag_movers.end(), numeric::random::rg() );
				for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
					if(option[OptionKeys::remodel::repeat_structure].user()){
						(*i)->apply( repeat_pose_ );
						repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "simul_frag" );
					}else {
						(*i)->apply( pose );
						mc.boltzmann( pose, "simul_frag" );
					}
				}
			} else {
				// per-loop ccd
				random_permutation( loops_to_model->v_begin(), loops_to_model->v_end(), numeric::random::rg() );
				for ( Loops::const_iterator l = loops_to_model->begin(), le = loops_to_model->end(); l != le; ++l ) {
					if ( !l->is_terminal( pose ) ) {
						CCDLoopClosureMover ccd_mover( *l, core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new MoveMap( movemap ) ) ) );
						ccd_mover.max_cycles( 50 );  // Used to be 10 moves, which would result in 50 "tries" in the old code. ~Labonte
						if(option[OptionKeys::remodel::repeat_structure].user()){
							if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
								ccd_mover.apply( repeat_pose_ );
							}
							repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "ccd_move" );
						}else {
							if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
								ccd_mover.apply( pose );
							}
							mc.boltzmann( pose, "ccd_move" );
						}
					}
				}
			}

			if( apply_user_provided_movers && ( inner % user_provided_movers_apply_cycle_ == 0 ) ){
				for( utility::vector1< moves::MoverOP >::iterator move_it( user_provided_movers_.begin() ); move_it != user_provided_movers_.end(); ++move_it ){
					(*move_it)->apply( pose );
					mc.boltzmann( pose, "user_provided_simul" );
					//if( inner % 50 == 0 ){
					//	static Size simulposecount = 0;
					//	simulposecount++;
					//	pose.dump_pdb("simulstage"+utility::to_string( simulposecount )+".pdb" );
					//}
				}
			}

		} // inner_cycles

	} // outer_cycles

	// recover low
	if(option[OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if(option[OptionKeys::remodel::repeat_structure] == 1){
				copy_size = unit_length_ - repeat_tail_length_;
			} else {
				copy_size = unit_length_;
			}
		for (Size res = 1; res<=copy_size; res++){
			pose.set_phi(res,mc.lowest_score_pose().phi(res));
			pose.set_psi(res,mc.lowest_score_pose().psi(res));
			pose.set_omega(res,mc.lowest_score_pose().omega(res));
		  pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
		}
	}
	else{
		pose = mc.lowest_score_pose();
	}

	// report status
//	mc.score_function().show_line_headers( TR );
//	TR << std::endl;

  if(option[OptionKeys::remodel::repeat_structure].user()){
//		mc.score_function().show_line( TR, repeat_pose_ );
		mc.score_function().show( TR, repeat_pose_ );
		check_closure_criteria( repeat_pose_, true );
		TR << std::endl;
		mc.show_state();
		remove_cutpoint_variants( repeat_pose_ );
		remove_cutpoint_variants( pose );
	}	else {
		TR << "\n";
		mc.score_function().show( TR, pose );
		TR << std::endl;
		mc.show_state();
		check_closure_criteria( pose, true );
		remove_cutpoint_variants( pose );
	}

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
	using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** independent_stage" << std::endl;

	// setup loops
	loops::LoopsOP pre_loops_to_model = determine_loops_to_model( pose );
	loops::LoopsOP loops_to_model( new loops::Loops() );

	// filter for non-terminal loops
	for ( Loops::const_iterator l = pre_loops_to_model->begin(), le = pre_loops_to_model->end(); l != le; ++l ) {
		if(option[OptionKeys::remodel::repeat_structure].user()){
			//take out the terminal loop in repeat cases
			if ( !l->is_terminal( pose ) ) {
					loops_to_model->add_loop( *l );
			}
			//however, need to address de novo building cases
			if (l->start() == 1 && l->stop() == unit_length_){
					loops_to_model->add_loop( *l );
			}
		}
		else {
				loops_to_model->add_loop( *l );
		}
	}

	if(option[OptionKeys::remodel::no_jumps].user()){
		// if using no_jumps, chainbreak based loop determination will skip this stage, but we obviously wants to build something...
		loops_to_model = loops_;
	}
/*
	//TEST 9/28/2012
  //find terminal loops and skip them; loophash can't handle terminal loops
	if(option[OptionKeys::remodel::repeat_structure].user()){
					for ( Loops::const_iterator l = loops_->begin(), le = loops_->end(); l != le; ++l ) {
						bool skip_loop = false;
						if ( l->is_terminal( pose ) || l->start() == 1) {
							skip_loop = true;
						}
						if ( !skip_loop ) {
							loops_to_model->add_loop( *l );
						}
					}
	}
*/
	TR << "   n_loops = " << loops_to_model->size() << std::endl;

	if ( loops_to_model->size() == 0 ) { // nothing to do...
		return;
	}

	// parameters
	Size const n_standard_cycles = total_standard_cycles();
	Size const max_outer_cycles = independent_cycles();

	// per-loop frag + ccd_move
	for ( Loops::iterator l = loops_to_model->v_begin(), le = loops_to_model->v_end(); l != le; ++l ) {
		Loop & loop = *l;

		// alter cutpoint to one before the end of the loop (either direction,
		// based on option) if closure is bypassed.  this is already set in VLB,
		// but reenforce here.
		if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
		//	if(option[OptionKeys::RemodelLoopMover::force_cutting_N].user()){
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
		bool apply_user_provided_movers( user_provided_movers_.size() != 0 );

		// set appropriate topology
		if ( keep_input_foldtree_ ){
		}
		else {
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				protocols::loops::set_single_loop_fold_tree( pose, loop );
			} else {
				if(option[OptionKeys::remodel::repeat_structure].user()){
					//do nothing. remake foldtree will destroy the pose when propagate across unconnected segments
				} else {
					protocols::forge::methods::set_single_loop_fold_tree( pose, loop );
				}
			}
		}

		// add cutpoint variants
		add_cutpoint_variants( pose );
		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			add_cutpoint_variants( repeat_pose_ );
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

			if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			if(option[OptionKeys::remodel::no_jumps].user() && option[OptionKeys::remodel::two_chain_tree].user() ){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			if(option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide].user()){
		//	sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			sfxOP->set_weight(
				core::scoring::atom_pair_constraint,
				sfxOP->get_weight( core::scoring::atom_pair_constraint) + cbreak_increment
			);
		}

			mc.score_function( *sfxOP );

			// recover low
			if(option[OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if(option[OptionKeys::remodel::repeat_structure] == 1){
				copy_size = unit_length_-1;
			} else {
				copy_size = unit_length_;
			}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
				  pose.set_secstruct(res, mc.lowest_score_pose().secstruct(res));
				}
			}
			else{
				pose = mc.lowest_score_pose();
			}

			if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);			}

			for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
				// fragments
				if ( loop.is_terminal( pose ) || numeric::random::rg().uniform() * n_standard_cycles > ( outer + simultaneous_cycles() ) || pose.fold_tree().num_cutpoint() == 0 ) {
					random_permutation( frag_movers.begin(), frag_movers.end(), numeric::random::rg() );
					for ( FragmentMoverOPs::iterator i = frag_movers.begin(), ie = frag_movers.end(); i != ie; ++i ) {
						if(option[OptionKeys::remodel::repeat_structure].user()){
							(*i)->apply( repeat_pose_ );
							repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "frag" );
						}else {
							(*i)->apply( pose );
							mc.boltzmann( pose, "frag" );
						}
					}
				} else { // ccd
					CCDLoopClosureMover ccd_mover( loop, core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new MoveMap( movemap ) ) ) );
					ccd_mover.max_cycles( 50 );  // Used to be 10 moves, which would result in 50 "tries" in the old code. ~Labonte
					if( option[ OptionKeys::remodel::repeat_structure ].user() ) {
						if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
							ccd_mover.apply( repeat_pose_ );
						}
						repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "ccd_move" );
					}else {
						if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
							ccd_mover.apply( pose );
						}
						mc.boltzmann( pose, "ccd_move" );
					}
				}

				if( apply_user_provided_movers && ( inner % user_provided_movers_apply_cycle_ == 0 ) ){
					for( utility::vector1< moves::MoverOP >::iterator move_it( user_provided_movers_.begin() ); move_it != user_provided_movers_.end(); ++move_it ){
						(*move_it)->apply( pose );
						mc.boltzmann( pose, "user_provided_indep" );
						//if( inner % 50 == 0 ){
						//	static Size indepposecount = 0;
						//	indepposecount++;
						//	pose.dump_pdb("indepstage"+utility::to_string( indepposecount )+".pdb" );
						//}
					}
				}

			} // inner_cycles

		} // outer_cycles

		// recover low
		if(option[OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if(option[OptionKeys::remodel::repeat_structure] == 1){
				copy_size = unit_length_-1;
			} else {
				copy_size = unit_length_;
			}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
				  pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
				}
		} else{
				pose = mc.lowest_score_pose();
		}

		// report status
	//	mc.score_function().show_line_headers( TR );
//	TR << std::endl;

		if(option[OptionKeys::remodel::repeat_structure].user()){
	//		mc.score_function().show_line( TR, repeat_pose_ );
			mc.score_function().show( TR, repeat_pose_ );
			remove_cutpoint_variants( repeat_pose_ );
	//	repeat_pose_.dump_scored_pdb("checkRepeat2.pdb", mc.score_function());
		}	else {
			TR << "\n";
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
	using protocols::loops::loop_closure::ccd::CCDLoopClosureMover;
	using protocols::loops::remove_cutpoint_variants;

	TR << "** boost_closure_stage" << std::endl;

	// setup loops
	loops::LoopsOP pre_loops_to_model = determine_loops_to_model( pose );
	loops::LoopsOP loops_to_model( new loops::Loops() );

	// filter for non-terminal loops that are within tolerance
	Real const cbreak_tolerance = 1.0;
	for ( Loops::const_iterator l = pre_loops_to_model->begin(), le = pre_loops_to_model->end(); l != le; ++l ) {
		if ( !l->is_terminal( pose ) || l->start() != 1 ) {
			Real const cbreak = linear_chainbreak( pose, l->cut() );
			if ( cbreak < cbreak_tolerance ) {
				loops_to_model->add_loop( *l );
			}
			else if (cbreak >= cbreak_tolerance ){
				TR << "at least one loop is beyond tolerance" << std::endl;
				return;
			}
		}
	}

	TR << "   n_loops = " << loops_to_model->size() << std::endl;
	if ( loops_to_model->size() == 0 ) { // nothing to do...
		TR << "   no loops to work on with boost_closure. returning." << std::endl;
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
	for ( Loops::const_iterator l = loops_to_model->begin(), le = loops_to_model->end(); l != le; ++l ) {
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
		if ( keep_input_foldtree_){
		}
		else {
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				protocols::loops::set_single_loop_fold_tree( pose, loop );
			} else {
				if(option[OptionKeys::remodel::repeat_structure].user()){
					//do nothing. remake foldtree will destroy the pose when propagate across unconnected segments
				} else {
					protocols::forge::methods::set_single_loop_fold_tree( pose, loop );
				}
			}
		}
		// add cutpoint variants
		add_cutpoint_variants( pose );
		if(option[OptionKeys::remodel::repeat_structure].user()){
			repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			add_cutpoint_variants( repeat_pose_ );
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

			if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
				sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
			}

			if(option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide].user()){
		//	sfxOP->set_weight( core::scoring::linear_chainbreak, 0);
				sfxOP->set_weight(
					core::scoring::atom_pair_constraint,
					sfxOP->get_weight( core::scoring::atom_pair_constraint) + cbreak_increment
				);
			}

			mc.score_function( *sfxOP );

      // recover low
      if(option[OptionKeys::remodel::repeat_structure].user()){
				Size copy_size =0;
				if(option[OptionKeys::remodel::repeat_structure] == 1){
					copy_size = unit_length_-1;
				} else {
					copy_size = unit_length_;
				}
				for (Size res = 1; res<=copy_size; res++){
					pose.set_phi(res,mc.lowest_score_pose().phi(res));
					pose.set_psi(res,mc.lowest_score_pose().psi(res));
					pose.set_omega(res,mc.lowest_score_pose().omega(res));
				  pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
				}
      } else {
				pose = mc.lowest_score_pose();
			}

			if(option[OptionKeys::remodel::repeat_structure].user()){
				repeat_propagation(pose, repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
			}

			// Going into the boost_closure stage implies we are "desperate" to close
			// the loop and don't care about diversity anymore, so continue to
			// cycle only until the loop is closed.
			if ( linear_chainbreak( pose, loop.cut() ) <= max_linear_chainbreak_ ) {
				break;
			}

			for ( Size inner = 1; inner <= max_inner_cycles; ++inner ) {
				if ( (!frag1_movers.empty() && numeric::random::rg().uniform() < frag_mover_probability) || pose.fold_tree().num_cutpoint() == 0 ) { // 1-mer insertions

					random_permutation( frag1_movers.begin(), frag1_movers.end(), numeric::random::rg() );
					for ( FragmentMoverOPs::iterator i = frag1_movers.begin(), ie = frag1_movers.end(); i != ie; ++i ) {
						if(option[OptionKeys::remodel::repeat_structure].user()){
							(*i)->apply( repeat_pose_ );
							repeat_sync(repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
							mc.boltzmann( repeat_pose_, "frag1" );
						}
						else{
							(*i)->apply( pose );
							mc.boltzmann( pose, "frag1" );
						}
					}

				} else { // ccd_move
					CCDLoopClosureMover ccd_mover( loop, core::kinematics::MoveMapCOP( core::kinematics::MoveMapOP( new MoveMap( movemap ) ) ) );
					ccd_mover.max_cycles( 50 );  // Used to be 10 moves, which would result in 50 "tries" in the old code. ~Labonte
					if(option[OptionKeys::remodel::repeat_structure].user()){
						if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
							ccd_mover.apply( repeat_pose_ );
						}
						repeat_sync( repeat_pose_,option[OptionKeys::remodel::repeat_structure]);
						mc.boltzmann( repeat_pose_, "frag1" );
					} else {
						if ( !option[ OptionKeys::remodel::RemodelLoopMover::bypass_closure ].user() ) {
							ccd_mover.apply( pose );
						}
						mc.boltzmann( pose, "ccd_move" );
					}
				}
			} // inner_cycles

		} // outer_cycles

		// recover low
		if(option[OptionKeys::remodel::repeat_structure].user()){
			Size copy_size =0;
		  if(option[OptionKeys::remodel::repeat_structure] == 1){
				copy_size = unit_length_-1;
			} else {
				copy_size = unit_length_;
			}
					for (Size res = 1; res<=copy_size; res++){
          pose.set_phi(res,mc.lowest_score_pose().phi(res));
          pose.set_psi(res,mc.lowest_score_pose().psi(res));
          pose.set_omega(res,mc.lowest_score_pose().omega(res));
				  pose.set_secstruct(res,mc.lowest_score_pose().secstruct(res));
					}
    } else{
        pose = mc.lowest_score_pose();
    }


		// report status
	//	mc.score_function().show_line_headers( TR );
//		TR << std::endl;

		if(option[OptionKeys::remodel::repeat_structure].user()){
	//		mc.score_function().show_line( TR, repeat_pose_ );
			mc.score_function().show( TR, repeat_pose_ );
		remove_cutpoint_variants( repeat_pose_ );
	//	repeat_pose_.dump_scored_pdb("checkRepeat3.pdb", mc.score_function());
		}	else {
			TR << "\n";
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
loops::LoopsOP RemodelLoopMover::determine_loops_to_model( Pose & pose ) {
	using protocols::forge::methods::linear_chainbreak;

	loops::LoopsOP loops_to_model( new loops::Loops() );

	for ( Loops::const_iterator l = loops_->begin(), le = loops_->end(); l != le; ++l ) {
		bool skip_loop = false;
		/*
		std::cout << "loop_ start " << l->start() << " loop_ end" << l->stop() <<  " loop_ cut " << l->cut() << std::endl;
		std::cout << "linear chainbreak " << linear_chainbreak( pose, l->cut() )  << std::endl;
		std::cout << "eval: linear_chainbreak( pose, l->cut() ) <= max_linear_chainbreak_ " << ( linear_chainbreak( pose, l->cut() ) <= max_linear_chainbreak_ )  << std::endl;
		*/
		if ( !l->is_terminal( pose ) ) {
			skip_loop |= ( linear_chainbreak( pose, l->cut() ) <= max_linear_chainbreak_ ); // loop already closed?
		}

		if ( !skip_loop ) {
			loops_to_model->add_loop( *l );
		}
	}

	// shuffle the order
	random_permutation( loops_to_model->v_begin(), loops_to_model->v_end(), numeric::random::rg() );

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

	boost::format format( "%|5t|%1% %|5t|%2% %|5t|%3% %|8t|%4%" );

	//breakout case if we don't care if the loops are closed
	if(option[OptionKeys::remodel::RemodelLoopMover::bypass_closure].user()){
		//special case for DB, filtering on close RMS with loophash
		//only trigger if lh_closure_filter is used before evaluation
		if(option[OptionKeys::remodel::lh_closure_filter].user()){
			TR << "using cbreak filter under bypass_closure." << std::endl;
		  bool all_loops_pass = true;
			for ( Loops::const_iterator l = loops_->begin(), le = loops_->end(); l != le; ++l ) {
				Real cbreak = 0.0;
				if ( !l->is_terminal( pose ) ) {
					cbreak = linear_chainbreak( pose, l->cut() );
					TR << "chain break " << cbreak << std::endl;
					all_loops_pass &= ( cbreak <= max_linear_chainbreak_ );
				}

				if ( show_in_tracer ) {
					TR << format % l->start() % l->stop() % l->cut() % cbreak << std::endl;
				}
			}
			return all_loops_pass;
		}
		else {
			return true;
		}
	}

	// boost::format here does not appear to be doing what I want it to do...
	// The format string is probably borked.
	if ( show_in_tracer ) {
		TR << format % "start" % "stop" % "cut" % "cbreak" << std::endl;
	}

	bool all_loops_pass = true;

	for ( Loops::const_iterator l = loops_->begin(), le = loops_->end(); l != le; ++l ) {
		Real cbreak = 0.0;
		if ( !l->is_terminal( pose ) ) {
			cbreak = linear_chainbreak( pose, l->cut() );
			all_loops_pass &= ( cbreak <= max_linear_chainbreak_ );
		}

		if ( show_in_tracer ) {
			TR << format % l->start() % l->stop() % l->cut() % cbreak << std::endl;
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
			ClassicFragmentMoverOP cfm( new ClassicFragmentMover( *f, movemap.clone() ) );
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			frag_movers.push_back( cfm );
		}

	}

	return frag_movers;
}

/// @brief return fragment movers for the list of internally kept fragment sets
/// @param[in] movemap Use this movemap when initializing fragment movers.
//  @param[in] note size 999 indicates any size > 9 fragments are allowed.
/// @param[in] limits size and position of fragments
RemodelLoopMover::FragmentMoverOPs
RemodelLoopMover::create_fragment_movers_limit_size(
	MoveMap const & movemap,
	Size const frag_size,
	std::set<Size> const & allowedPos,
	bool const smoothMoves,
	Real fragScoreThreshold
	)
{
	using namespace protocols::simple_moves;
	using namespace core::fragment;
	using namespace basic::options;
	using namespace OptionKeys::indexed_structure_store;
	FragmentMoverOPs frag_movers;
	for ( FragSetOPs::const_iterator f = fragsets_.begin(); f != fragsets_.end(); ++f ) {
		if((*f)->max_frag_length()==frag_size || ((frag_size == 999)&&((*f)->max_frag_length()>9))) {
			ConstantLengthFragSetOP tmp_frags( new ConstantLengthFragSet((*f)->max_frag_length()) );
			for( ConstFrameIterator frame_i = (*f)->begin(); frame_i != (*f)->end(); ++frame_i ){
				if(allowedPos.find((*frame_i)->start()) != allowedPos.end()){
						FrameOP tmp_frame = (*frame_i)->clone();
						Size frag_ct= 0;
						for(Size ii = 1; ii<=(*frame_i)->nr_frags(); ++ii){
								if((*frame_i)->fragment(ii).score()<fragScoreThreshold){
										frag_ct++;
										tmp_frame->add_fragment((*frame_i)->fragment_ptr(ii));
									}
						}
						assert(frag_ct != 0);//0 frags at this position. You have chosen a bad abego definition.
						tmp_frags->add(tmp_frame);
					}
			}
			ClassicFragmentMoverOP cfm;
			if(option[fragment_threshold_distance].user()){
				cfm = ClassicFragmentMoverOP(new VallLookbackFragMover(*f,movemap.clone()));
			}
			else{
				if(smoothMoves)
					cfm = ClassicFragmentMoverOP( new SmoothFragmentMover( *f, movemap.clone(), FragmentCostOP( new GunnCost ) ) );
				else
					cfm = ClassicFragmentMoverOP( new ClassicFragmentMover( *f, movemap.clone() ) );
				}
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
			ClassicFragmentMoverOP cfm( new ClassicFragmentMover( *f, movemap.clone() ) );
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
	loops::LoopsOP const loops,
	Size const largest_frag_size
)
{
	// Create fragment movers for each loop for each fragment set.  Here we
	// want to allow an equal probability of movement per-loop, rather than
	// per-residue.
	FragmentMoverOPs frag_movers;
	for ( Loops::const_iterator l = loops->begin(), le = loops->end(); l != le; ++l ) {
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
	loops::LoopsOP const loops,
	MoveMap & movemap,
	bool const allow_omega
)
{
	for ( Loops::const_iterator l = loops->begin(), le = loops->end(); l != le; ++l ) {
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

///@brief copies phi,psi,omega from a starting pose.
void RemodelLoopMover::set_starting_pdb(Pose & pose){
		using namespace basic::options;
		using namespace OptionKeys::remodel;
		using core::Size;
		PoseOP inputPose = core::import_pose::pose_from_pdb(option[OptionKeys::remodel::staged_sampling::starting_pdb]);
		assert(inputPose->total_residue() == pose.total_residue());
		for(Size ii=2; ii<pose.total_residue(); ++ii){
				pose.set_phi(ii,inputPose->phi(ii));
				pose.set_psi(ii,inputPose->psi(ii));
				pose.set_omega(ii,inputPose->omega(ii));
				pose.set_secstruct(ii,inputPose->secstruct(ii));
		}
		Size repeatRes = (pose.total_residue()/2)+1;
		pose.set_phi(1,inputPose->phi(repeatRes));
		pose.set_psi(1,inputPose->psi(repeatRes));
		pose.set_omega(1,inputPose->omega(repeatRes));
		pose.set_secstruct(1,inputPose->secstruct(repeatRes));
}

///@brief sets helices to there ideal value before any sampling begins.
void RemodelLoopMover::set_ideal_helices(Pose & pose){
		using namespace basic::options;
		using namespace OptionKeys::remodel;
		using core::Size;
		//At this point the pose is 2x. So we need to copy the helical residues twice.
		std::string ss = remodel_data_.ss;
		Size repeatRes = (pose.total_residue()/2);
		for(Size ii=1; ii<=ss.size(); ++ii){
				if(ss[ii-1] == 'H'){
						pose.set_phi(ii,-63.8);
						pose.set_phi(ii+repeatRes,-63.8);
						pose.set_psi(ii,-41.1);
						pose.set_psi(ii+repeatRes,-41.1);
						pose.set_omega(ii,180);
						pose.set_omega(ii+repeatRes,180);
						pose.set_secstruct(ii,'H');
						pose.set_secstruct(ii+repeatRes,'H');
				}
		}
}

void RemodelLoopMover::set_starting_sequence(Pose & pose){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using core::Size;
		std::string const & swap_sequence =option[OptionKeys::remodel::staged_sampling::starting_sequence];
		for(Size ii=1; ii<=swap_sequence.size(); ++ii){
		  char aa = swap_sequence[ii-1];
		  AA my_aa = aa_from_oneletter_code( aa );
			ResidueTypeSetCOP const &residue_set(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID ));
		  ResidueTypeCOPs const & rsd_type_list( residue_set->aa_map( my_aa ) );
		  ResidueType const & rsd_type( *(rsd_type_list[ 1 ]) );
		  if(option[OptionKeys::remodel::repeat_structure].user()){
				replace_pose_residue_copying_existing_coordinates(pose,ii,rsd_type);//pose has two coppies. This is the first
				replace_pose_residue_copying_existing_coordinates(pose,ii+swap_sequence.size(),rsd_type);
			}
			else
				replace_pose_residue_copying_existing_coordinates(pose,ii,rsd_type);
		}
}

std::set<core::Size> RemodelLoopMover::generate_residues_to_sample(bool chooseSubsetResidues, Pose & pose,Size fragmentSize){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using core::Size;
		std::set<Size> allowedRes;
		if(!chooseSubsetResidues || !option[OptionKeys::remodel::staged_sampling::residues_to_sample].user() || !option[OptionKeys::remodel::staged_sampling::sample_over_loops].user()){
			for(Size ii=1; ii<=pose.total_residue(); ++ii){
				allowedRes.insert(ii);
			}
		}
		else{
				if(option[OptionKeys::remodel::staged_sampling::residues_to_sample].user()){
						std::string const & allowedRes_str =option[OptionKeys::remodel::staged_sampling::residues_to_sample];
				utility::vector1< std::string > const res_keys( utility::string_split( allowedRes_str , ',' ) );
				BOOST_FOREACH( std::string const key, res_keys ){
						Size const res( utility::string2int( key ) );
						allowedRes.insert(res);
						}
				}
				if(option[OptionKeys::remodel::staged_sampling::sample_over_loops].user()){
				std::string ss = remodel_data_.ss;
				Size repeatRes = (pose.total_residue()/2)+1;
				char lastRes = ss[0];
				for(Size ii=1; ii<=ss.size(); ++ii){
						if((ss[ii-1] == 'H' || ss[ii-1] == 'E')&&(lastRes == 'L')){
								if(ii-(fragmentSize-1) > 0)
										allowedRes.insert(ii-(fragmentSize-1));
										allowedRes.insert(ii-(fragmentSize-1)+repeatRes);
					}
				}
				//edge case
				if(((ss[ss.size()-1]=='H')||(ss[ss.size()-1]=='E'))&&(ss[0] == 'L')){
						allowedRes.insert(ss.size());
				}
				}
		}
		return(allowedRes);
}

} // remodel
} // forge
} // protocols

