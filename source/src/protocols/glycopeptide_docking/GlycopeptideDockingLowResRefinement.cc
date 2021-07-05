// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingLowResRefinement.cc
/// @details Low resolution refinement of glycosyltransferase peptide/glycopeptide substrate complex.
/// This protocol allows enhanced sampling of the peptide/glycopeptide dofs in the low resolution mode with
/// simulated annealing. Current protocol was only tested for monosaccharides.
/// It does not treat the glycan-branch in a special manner. This behavior will change as the protocol is
/// developed in future versions for higher glycans.
/// Protocol similar to low resolution refinement of flexpepdock.
/// @author Yashes Srinivasan (yashess@gmail.com), Sai Pooja Mahajan (saipooja@gmail.com)

// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingLowResRefinement.hh>

// Protocol headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh>
#include <protocols/glycopeptide_docking/utils.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/FoldTree.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>

// XSD Includes

#include <string>
#include <cmath>

static basic::Tracer TR( "protocols.glycopeptide_docking.GlycopeptideDockingLowResRefinement" );

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
GlycopeptideDockingLowResRefinement::GlycopeptideDockingLowResRefinement():
	protocols::moves::Mover( GlycopeptideDockingLowResRefinement::mover_name() )
{

}

/// @brief Sets up movemap for torsion sampling.
/// Sets jumps. Torsions are sampled for backbone and chi angle
/// for the peptide substrate only (along with glycans that branch out).
/// Sets up minimizers for peptide torsion sampling and for enzyme-peptide interface jump.

GlycopeptideDockingLowResRefinement::GlycopeptideDockingLowResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf, core::kinematics::FoldTreeOP ft_docking, core::kinematics::FoldTreeOP ft_substrate ): protocols::moves::Mover( GlycopeptideDockingLowResRefinement::mover_name() ), flags_( flags ), sf_( sf ), ft_docking_(ft_docking), ft_substrate_(ft_substrate)
{
	setup_backbone_and_jump_movemap();
}

GlycopeptideDockingLowResRefinement::GlycopeptideDockingLowResRefinement( protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::scoring::ScoreFunctionOP sf):protocols::moves::Mover( GlycopeptideDockingLowResRefinement::mover_name() ), flags_( flags ), sf_(sf)
{
	setup_backbone_and_jump_movemap();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycopeptideDockingLowResRefinement::GlycopeptideDockingLowResRefinement( GlycopeptideDockingLowResRefinement const & src ):
	protocols::moves::Mover( src )
{

}

void
GlycopeptideDockingLowResRefinement::setup_backbone_and_jump_movemap()
{
	using namespace scoring;
	using namespace minimization_packing;
	using namespace kinematics;

	MoveMapOP torsion_mm( utility::pointer::make_shared< MoveMap >() );
	torsion_mm->set_bb_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	torsion_mm->set_chi_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	torsion_mm->set_branches_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	torsion_mm->set_jump( flags_->jump_num_substrate(), true );

	torsion_minimizer_ = MinMoverOP( utility::pointer::make_shared< MinMover>( torsion_mm, sf_, "lbfgs_armijo_nonmonotone", 0.01, true ) );

	MoveMapOP jump_mm( utility::pointer::make_shared< MoveMap>() );
	jump_mm->set_jump( flags_->jump_num_substrate(), true );

	jump_minimizer_ = MinMoverOP( utility::pointer::make_shared< MinMover >( jump_mm, sf_, "lbfgs_armijo_nonmonotone", 0.001, true ) );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycopeptideDockingLowResRefinement::~GlycopeptideDockingLowResRefinement(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @details Low resolution refinement primarly for sampling
/// peptide conformational space. We use simulated annealing
/// varying mc temperature from 2.0 to 0.6 (final) temperature.
/// The protocol applies the following steps over 30 cycles of simulated annealing:
/// 1. rigidbody perturber mover\n
/// 2. minimization across the jump\n
/// 3. Samples torsions of the peptide with small and shear mover\n
/// 4. Minimizes peptide\n
void
GlycopeptideDockingLowResRefinement::apply( core::pose::Pose& pose ){
	using namespace rigid;
	/*Ramp down temperature form high temperature (2.0) to Room Temperature 0.6
	* Higher temperature allows faster sampling at the beginning*/

	core::Real init_MC_temp = 2.0;
	core::Real final_MC_temp = 0.6;
	core::Real gamma = pow( (final_MC_temp / init_MC_temp) , 1.0 / (flags_->low_res_outer_cycles() - 1.0) );
	core::Real MC_temp = init_MC_temp;
	core::Real torsions_acceptance;
	core::Real distance_substrate_enzyme(0.0);

	if ( !ft_docking_ ) {
		setup_glycosylation_foldtree(pose,flags_,ft_docking_);
	} else {
		pose.fold_tree(*ft_docking_);
	}
	TR << endl << " Low res fold tree: " << endl;
	TR << " " << pose.fold_tree() << endl;

	mc_ = moves::MonteCarloOP( utility::pointer::make_shared< moves::MonteCarlo>( pose, *sf_, MC_temp ) );

	core::kinematics::MoveMapOP mm_upstream( utility::pointer::make_shared< core::kinematics::MoveMap >() );
	mm_upstream->set_bb_true_range( flags_->first_residue_substrate(), flags_->anchor_residue_substrate() - 1 );
	mm_upstream->set_chi_true_range( flags_->first_residue_substrate(), flags_->anchor_residue_substrate() - 1 );
	core::kinematics::MoveMapOP mm_downstream( utility::pointer::make_shared< core::kinematics::MoveMap>() );
	mm_downstream->set_bb_true_range( flags_->anchor_residue_substrate() + 1, flags_->last_residue_substrate() );
	mm_downstream->set_chi_true_range( flags_->anchor_residue_substrate() + 1, flags_->last_residue_substrate() );
	core::kinematics::MoveMapOP mm_full( utility::pointer::make_shared< core::kinematics::MoveMap>() );
	mm_full->set_bb_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );
	mm_full->set_chi_true_range( flags_->first_residue_substrate(), flags_->last_residue_substrate() );

	protocols::rigid::RigidBodyPerturbMoverOP perturber = protocols::rigid::RigidBodyPerturbMoverOP ( utility::pointer::make_shared< protocols::rigid::RigidBodyPerturbMover>( flags_->jump_num_substrate(), 2, 0.2 ) );
	distance_substrate_enzyme = calculate_sampled_distance(pose,flags_->glycosylation_residue_substrate(),flags_->get_sugar_donor());
	TR<< "Initial distance "<<distance_substrate_enzyme<<std::endl;

	for ( core::Size cycle_out = 1; cycle_out <= flags_->low_res_outer_cycles(); cycle_out++ ) {

		mc_->set_temperature( MC_temp );
		//pose.fold_tree(*ft_docking_);
		perturber->apply( pose );
		jump_minimizer_->apply( pose );
		sample_torsions( pose, flags_->low_res_inner_cycles(), torsions_acceptance,  mm_upstream );
		sample_torsions( pose, flags_->low_res_inner_cycles(), torsions_acceptance,  mm_downstream );
		torsion_minimizer_->apply( pose );

		MC_temp = MC_temp*gamma;
		if ( flags_->debug_pdbs() ) {
			jd2::JobOP job( jd2::JobDistributor::get_instance()->current_job() );
			std::ostringstream ss;
			ss << "lowresposedump_cycle"<< cycle_out;
			TR << "Dumping Low Res Pdb "<<ss.str()<<std::endl;
			write_debug_pdb(pose,job->nstruct_max(),job->nstruct_index(),ss.str());
		}
		distance_substrate_enzyme = calculate_sampled_distance(pose,flags_->glycosylation_residue_substrate(),flags_->get_sugar_donor());
		if ( distance_substrate_enzyme > 8.0 ) {
			TR<< "Exiting Low Res - sampled distance "<<distance_substrate_enzyme<<std::endl;
			// pymol_mover_->apply(pose);
			//reset to initial MC temperature before exiting.
			mc_->set_temperature( init_MC_temp );
			return;
		}
	}

	mc_->set_temperature( final_MC_temp );
	//Additional sampling at the final MC_temperature
	for ( core::Size cycle = 1; cycle <= flags_->low_res_outer_cycles(); cycle++ ) {
		perturber->apply( pose );
		jump_minimizer_->apply( pose );
		sample_torsions( pose, flags_->low_res_inner_cycles(), torsions_acceptance,  mm_upstream );
		sample_torsions( pose, flags_->low_res_inner_cycles(), torsions_acceptance,  mm_downstream );
		torsion_minimizer_->apply( pose );
		// pymol_mover_->apply(pose);
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycopeptideDockingLowResRefinement::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingLowResRefinement::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingLowResRefinement>() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingLowResRefinement::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingLowResRefinement>( *this ) );
}

std::string GlycopeptideDockingLowResRefinement::get_name() const {
	return mover_name();
}

std::string GlycopeptideDockingLowResRefinement::mover_name() {
	return "GlycopeptideDockingLowResRefinement";
}


////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

/// @details Samples torsional degrees of freedom with small and shear mover.
void
GlycopeptideDockingLowResRefinement::sample_torsions( core::pose::Pose & pose, core::Size cycles, core::Real & acceptance_rate, core::kinematics::MoveMapOP mm )
{
	using namespace protocols::moves;
	using namespace protocols::simple_moves;

	// setup sub-moves
	SmallMoverOP small_mover( utility::pointer::make_shared< protocols::simple_moves::SmallMover>( mm, mc_->temperature(), 10 ) );
	ShearMoverOP shear_mover( utility::pointer::make_shared< protocols::simple_moves::ShearMover>( mm, mc_->temperature(), 10 ) );

	RandomMoverOP random_mover( utility::pointer::make_shared< protocols::moves::RandomMover>() );
	random_mover->add_mover( small_mover, 1.0 );
	random_mover->add_mover( shear_mover, 1.0 );

	// wrap with monte-carlo trial mover
	TrialMoverOP mc_trial( utility::pointer::make_shared< TrialMover>( random_mover, mc_ ) );
	mc_trial->keep_stats_type( accept_reject ); // track stats (for acceptance rate)

	for ( core::Size cycle = 1; cycle <= cycles; cycle++ ) {
		mc_trial->apply( pose );
	}

	pose = mc_->lowest_score_pose();
	mc_->reset( pose );
	acceptance_rate = mc_trial->acceptance_rate();
}

std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingLowResRefinement const & mover )
{
	mover.show(os);
	return os;
}


} //glycopeptide_docking
} //protocols
