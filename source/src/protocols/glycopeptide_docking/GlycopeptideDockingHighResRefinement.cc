// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingHighResRefinement.cc
/// @details High resolution refinement of glycosyltransferase peptide/glycopeptide substrate complex.
/// This protocol applies MCM sampling in torsional and rigid body space with sidechain
/// repacking of the substrate and substrate-enzyme interface. The attractive and repulsive terms are
/// ramped down and up respectively to allow softer potential at the beginning of the run.
/// Protocol similar to high resolution refinement of flexpepdock.
/// @author Yashes Srinivasan (yashess@gmail.com), Sai Pooja Mahajan (saipooja@gmail.com)

// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingHighResRefinement.hh>

// Protocol headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh>
#include <protocols/glycopeptide_docking/utils.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSpanSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>

// Numeric headers

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>

// XSD Includes

#include <string>

static basic::Tracer TR( "protocols.glycopeptide_docking.GlycopeptideDockingHighResRefinement" );

namespace protocols {
namespace glycopeptide_docking {

using namespace std;
using namespace utility;
using namespace core;
using namespace protocols;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycopeptideDockingHighResRefinement::GlycopeptideDockingHighResRefinement():
	protocols::moves::Mover( GlycopeptideDockingHighResRefinement::mover_name() )
{
}

/// @brief Constructor with arguments
GlycopeptideDockingHighResRefinement::GlycopeptideDockingHighResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf, core::kinematics::FoldTreeOP ft_docking, core::kinematics::FoldTreeOP ft_substrate):
	protocols::moves::Mover( GlycopeptideDockingHighResRefinement::mover_name() ), flags_( flags ), sf_( sf ), ft_docking_(ft_docking), ft_substrate_(ft_substrate)
{
	setup();
}

GlycopeptideDockingHighResRefinement::GlycopeptideDockingHighResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf):
	protocols::moves::Mover( GlycopeptideDockingHighResRefinement::mover_name() ), flags_( flags ), sf_( sf )
{
	setup();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycopeptideDockingHighResRefinement::GlycopeptideDockingHighResRefinement( GlycopeptideDockingHighResRefinement const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycopeptideDockingHighResRefinement::~GlycopeptideDockingHighResRefinement(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////
void
GlycopeptideDockingHighResRefinement::setup(){
	using namespace minimization_packing;
	using namespace pack::task;
	using namespace rigid;
	using namespace kinematics;
	using namespace simple_moves;
	using namespace docking;
	using namespace scoring;
	using namespace select;


	//////////////////////////
	//////// Defaults ///////
	//////////////////////////

	target_atr_ = sf_->get_weight( fa_atr );
	target_rep_ = sf_->get_weight( fa_rep );
	kt_ = 0.8;

	/////////////////////////
	////// Constants ///////
	////////////////////////

	starting_ramp_down_factor_ = 3.25;
	starting_ramp_up_factor_ = 0.45455;

	//////////////////////////
	//////// MoveMaps ////////
	//////////////////////////

	std::cout<<"setting stuff up"<<std::endl;
	MoveMapOP jump_mm( utility::pointer::make_shared< MoveMap>() );
	jump_mm->set_jump( flags_->jump_num_substrate(), true );

	std::cout<<"setting up movemap torsion"<<std::endl;
	MoveMapOP torsion_mm( utility::pointer::make_shared< MoveMap >() );
	torsion_mm->set_bb_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	torsion_mm->set_chi_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	if ( flags_->get_allow_glycan_torsion_moves() ) {
		torsion_mm->set_bb( flags_->last_residue_substrate() + 1 , true); //should be applied to all sugar residues
		torsion_mm->set_chi( flags_->last_residue_substrate() + 1 , true);
	}
	torsion_mm->set_jump( flags_->jump_num_substrate(), true );
	torsion_mm->set_jump( flags_->jump_num_substrate(), true );
	torsion_mm->set_branches( true );


	//////////////////////////
	//////// Movers //////////
	//////////////////////////
	std::cout<<"setting up small mover"<<std::endl;
	small_mover_ = SmallMoverOP( utility::pointer::make_shared<  SmallMover>( torsion_mm, kt_, 3 /* nmoves? */ ) );
	shear_mover_ = ShearMoverOP( utility::pointer::make_shared< ShearMover>( torsion_mm, kt_, 3 /* nmoves? */ ) );
	//linkage_conformer_mover_ = LinkageConformerMoverOP( new LinkageConformerMover( ... ) );      /* TODO: Use for glycan elongation */
	perturber_ = RigidBodyPerturbMoverOP ( utility::pointer::make_shared< RigidBodyPerturbMover>( flags_->jump_num_substrate(), 3, 0.2 ) );
	// slider_ = FaDockingSlideIntoContactOP ( new FaDockingSlideIntoContact( flags_->jump_num_interface() ) );     /* TODO: Write a function that slides the substrate */
	std::cout<<"setting up minimizer"<<std::endl;
	torsion_minimizer_ = MinMoverOP( utility::pointer::make_shared< MinMover>( torsion_mm, sf_, "lbfgs_armijo_nonmonotone", 0.001, true ) );
	jump_minimizer_ = MinMoverOP( utility::pointer::make_shared< MinMover>( jump_mm, sf_, "lbfgs_armijo_nonmonotone", 0.001, true ) );

	// Residue selectors
	residue_selector::ResidueSpanSelectorOP substrate( utility::pointer::make_shared< residue_selector::ResidueSpanSelector >( flags_->first_residue_substrate(), flags_->last_residue_substrate() ) );
	residue_selector::ResidueSpanSelectorOP anchor( utility::pointer::make_shared< residue_selector::ResidueSpanSelector >( flags_->anchor_residue_substrate(), flags_->anchor_residue_substrate() ) );
	residue_selector::NotResidueSelectorOP not_substrate( utility::pointer::make_shared< residue_selector::NotResidueSelector>( substrate ) );
	residue_selector::GlycanResidueSelectorOP glycan( utility::pointer::make_shared< residue_selector::GlycanResidueSelector>() );

	// Task operations
	operation::RestrictToRepackingRLTOP rtrp_rlt( utility::pointer::make_shared< operation::RestrictToRepackingRLT>() );
	operation::PreventRepackingRLTOP prp_rlt( utility::pointer::make_shared< operation::PreventRepackingRLT>() );
	operation::RestrictToRepackingOP rtrp( utility::pointer::make_shared< operation::RestrictToRepacking>() );
	operation::IncludeCurrentOP ic( utility::pointer::make_shared< operation::IncludeCurrent>() );
	simple_task_operations::RestrictToInterfaceOP rti( utility::pointer::make_shared< simple_task_operations::RestrictToInterface>( flags_->jump_num_substrate(), flags_->get_interface_distance() ) );

	// Generate task factories and packers
	TaskFactoryOP tf_substrate( utility::pointer::make_shared<TaskFactory>() );
	tf_substrate->push_back( operation::OperateOnResidueSubsetOP( utility::pointer::make_shared< operation::OperateOnResidueSubset>( rtrp_rlt, substrate ) ) );
	tf_substrate->push_back( operation::OperateOnResidueSubsetOP( utility::pointer::make_shared< operation::OperateOnResidueSubset>( prp_rlt, not_substrate ) ) );
	if ( flags_->pack_anchor() ) {
		tf_substrate->push_back( operation::OperateOnResidueSubsetOP( utility::pointer::make_shared< operation::OperateOnResidueSubset>( prp_rlt, anchor ) ) );
		tf_substrate->push_back( operation::OperateOnResidueSubsetOP( utility::pointer::make_shared< operation::OperateOnResidueSubset>( prp_rlt, glycan ) ) );
	}
	packer_substrate_ = PackRotamersMoverOP( utility::pointer::make_shared< PackRotamersMover>( sf_ ) );
	packer_substrate_->task_factory( tf_substrate );

	TaskFactoryOP tf_interface( utility::pointer::make_shared< TaskFactory>() );
	tf_interface->push_back( rtrp );
	tf_interface->push_back( ic );
	tf_interface->push_back( rti );
	if ( flags_->pack_anchor() ) tf_interface->push_back( operation::OperateOnResidueSubsetOP( utility::pointer::make_shared< operation::OperateOnResidueSubset>( prp_rlt, anchor ) ) );
	// tf_interface->push_back( operation::OperateOnResidueSubsetOP( new operation::OperateOnResidueSubset( prp_rlt, glycan ) ) );  /* TODO: Run more tests to see if this is beneficial? */
	packer_interface_ = PackRotamersMoverOP( utility::pointer::make_shared< PackRotamersMover>( sf_ ) );
	packer_interface_->task_factory( tf_interface );
}

/// @details Applies a set of samplers to the pose in high resolution mode.
/// At the beginning, the fa_atr term is upweighted and the fa_rep term is downweighted.
/// The terms are ramped down and up respectively to their score function weights over
/// n_cycles_. For  each cycles, the pose is subjected to: \n
/// 1. N cycles of rigid body moves
/// 1.1. Rigid body perturbation
/// 1.2 Sidechain packing of substrate
/// 1.3 Sidechain packing of interface
/// 1.4 Minimization across the substrate-enzyme jump
/// 1.5 Metropolis criterion (MCM with 1.4)
/// 2. N cycles of torsion moves "induced fit"
/// 2.1 Apply small or shear mover
/// 2.2 Sidechain packing of substrate
/// 2.3 Sidechain packing of interface
/// 2.4 Minimize and apply Metropolis criterion (MCM)

/// TODO: add sugar bbsampler to sample torsions for
/// higher glycans. Protocol currently tested for monosaccharides.
/// TODO: Test sampling of "donor residue" in the torsional space.
void
GlycopeptideDockingHighResRefinement::apply( core::pose::Pose& pose ){

	if ( !ft_docking_ ) { //can this be done better?
		setup_glycosylation_foldtree(pose,flags_,ft_docking_);
	} else {
		pose.fold_tree(*ft_docking_);
	}
	TR << endl << " High res fold tree: " << endl;
	TR << " " << pose.fold_tree() << endl;
	show(cout);
	core::Real fa_atr_score_weight = sf_->get_weight( scoring::fa_atr);
	core::Real fa_rep_score_weight = sf_->get_weight( scoring::fa_rep);
	mc_ = moves::MonteCarloOP( utility::pointer::make_shared< moves::MonteCarlo>( pose, *sf_, kt_ ) );
	sf_->set_weight( scoring::fa_atr, target_atr_ * starting_ramp_down_factor_ );
	sf_->set_weight( scoring::fa_rep, target_atr_ * starting_ramp_up_factor_ );
	show(cout);

	for ( core::Size cycle = 1; cycle <= flags_->high_res_outer_cycles(); cycle++ ) {
		core::Size tenper = flags_->high_res_outer_cycles() / 10;
		if ( tenper == 0 ) { tenper = 1; }
		if ( cycle % tenper == 0 ) {  // Ramp every ~10% of n_cycles.
			Real fraction = Real( cycle ) / flags_->high_res_outer_cycles();
			ramp_score_weight( scoring::fa_atr, target_atr_, fraction );
			ramp_score_weight( scoring::fa_rep, target_rep_, fraction );
			mc_->reset( pose );
		}
		// Rigid body moves /* TODO: Make this its own function to expand */
		for ( core::Size rigid_body_moves = 1; rigid_body_moves <= flags_->get_backbone_moves(); rigid_body_moves++ ) {
			perturber_->apply( pose );
			packer_substrate_->apply( pose );
			if ( rigid_body_moves == flags_->get_backbone_moves() / 2 ) packer_interface_->apply( pose );
			jump_minimizer_->apply( pose );
			mc_->boltzmann( pose );
		}

		// Substrate torsion moves /* TODO: Make this its own function to expand */
		if ( flags_->enable_backbone_moves_substrate() ) {
			for ( core::Size backbone_moves = 1; backbone_moves <= flags_->get_backbone_moves(); backbone_moves++ ) {
				if ( backbone_moves % 2 == 0 ) {
					small_mover_->apply( pose );
					packer_substrate_->apply( pose );
				} else {
					shear_mover_->apply( pose );
					packer_substrate_->apply( pose );
				}
				/*   else if (backbone_moves % 4 == 2) {
				pose.fold_tree(*ft_substrate_);
				small_mover_->apply( pose );
				packer_substrate_->apply( pose );
				}
				else {
				std::cout<<"Backbone shear substrate "<<std::endl;
				pose.fold_tree(*ft_substrate_);
				shear_mover_->apply( pose );
				packer_substrate_->apply( pose );
				}
				*/
				//pymol_mover_->apply(pose);
				if ( flags_->debug_pdbs() ) {
					std::ostringstream ss;
					ss << "highresposedump_cycle"<<cycle<<"_bk"<< backbone_moves;
					jd2::JobOP job( jd2::JobDistributor::get_instance()->current_job() );
					TR.Debug<< "Dumping Low Res Pdb "<<ss.str()<<std::endl;
					write_debug_pdb(pose,job->nstruct_max(),job->nstruct_index(),ss.str());
				}

				if ( backbone_moves % flags_->get_interface_moves() == 0 ) {
					//pose.fold_tree(*ft_docking_);
					packer_interface_->apply( pose );
				}
				torsion_minimizer_->apply( pose );
				mc_->boltzmann( pose );
			}
		}
		core::Real sampled_distance(0.0);
		sampled_distance = calculate_sampled_distance(pose,flags_->glycosylation_residue_substrate(),flags_->get_sugar_donor());
		if ( sampled_distance > 10.0 ) {
			TR<< "Exiting High Res - sampled distance "<<sampled_distance<<std::endl;
			//Set weights to original weight before exiting
			sf_->set_weight( scoring::fa_atr, fa_atr_score_weight );
			sf_->set_weight( scoring::fa_rep, fa_rep_score_weight );
			break;
		}


	}//end of outer for loop
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycopeptideDockingHighResRefinement::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingHighResRefinement::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared<  GlycopeptideDockingHighResRefinement>() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingHighResRefinement::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingHighResRefinement>( *this ) );
}

std::string GlycopeptideDockingHighResRefinement::get_name() const {
	return mover_name();
}

std::string GlycopeptideDockingHighResRefinement::mover_name() {
	return "GlycopeptideDockingHighResRefinement";
}


////////////////////////////////////////////////////////////////////////////////
//// private methods ////
/////////////////////////

/// @details Function to ramp up and ramp down score weights.
/// fa_atr is ramped down from an upweighted weight.
/// fa_rep is ramped up from an downweighted weight.
/// The refinement protocol currently uses it to ramp
/// down attractive terms and ramp up repulsive terms.
void
GlycopeptideDockingHighResRefinement::ramp_score_weight( core::scoring::ScoreType const method,
	core::Real const target,
	core::Real const fraction_completion )
{
	Real factor;
	Real const current_weight( sf_->get_weight( method ) );

	if ( current_weight < target ) {
		Real const factor_range_size( 1 - starting_ramp_up_factor_ );
		factor = starting_ramp_up_factor_ + fraction_completion * factor_range_size;
	} else if ( current_weight > target ) {
		Real const factor_range_size( starting_ramp_down_factor_ - 1 );
		factor = starting_ramp_down_factor_ - fraction_completion * factor_range_size;
	} else {
		factor = 1;
	}
	Real const new_weight( target * factor );
	sf_->set_weight( method, new_weight );
}


std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingHighResRefinement const & mover )
{
	mover.show(os);
	return os;
}


} //glycopeptide_docking
} //protocols


