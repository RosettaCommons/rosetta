// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/scoring_util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/scoring_util.hh>
#include <protocols/stepwise/modeler/constraint_util.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

using namespace core;

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.scoring_util" );

namespace protocols {
namespace stepwise {
namespace modeler {

/////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
get_minimize_scorefxn( core::pose::Pose const & pose,
	core::scoring::ScoreFunctionCOP scorefxn,
	options::StepWiseModelerOptionsCOP options ){
	using namespace core::scoring;
	ScoreFunctionOP minimize_scorefxn = scorefxn->clone();
	if ( minimize_scorefxn->get_weight( atom_pair_constraint ) == 0.0 ) minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 ); //go ahead and turn these on
	if ( minimize_scorefxn->get_weight( coordinate_constraint ) == 0.0 && options->coordinate_constraints_during_minimize() ) {
		minimize_scorefxn->set_weight( coordinate_constraint, 1.0 ); // go ahead and turn these on
	}
	if ( options->coordinate_constraints_during_minimize() ) {
		// this will exit with error if coordinate_constraints are not turned on
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, minimize_scorefxn );
	}
	// minimize_scorefxn->set_weight( linear_chainbreak, 150.0 ); // original SWA protein value.
	minimize_scorefxn->set_weight( linear_chainbreak, 5.0 ); // unify with RNA.
	if ( options->turn_off_rna_chem_map_during_optimize() )  minimize_scorefxn->set_weight( rna_chem_map, 0.0 ); // Just for now...
	if ( options->cart_min() && ( minimize_scorefxn->get_weight( cart_bonded ) == 0.0 ) ) minimize_scorefxn->set_weight( cart_bonded, 1.0 );
	if ( options->mapfile_activated() && minimize_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) minimize_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
	return minimize_scorefxn;
}


////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
initialize_sample_scorefxn( core::scoring::ScoreFunctionCOP scorefxn,
	pose::Pose const & pose,
	options::StepWiseModelerOptionsCOP options ){

	using namespace core::scoring;

	ScoreFunctionOP sample_scorefxn;
	std::string const & sample_weights = options->pack_weights();

	if ( sample_weights.size() == 0 ) {
		sample_scorefxn = scorefxn->clone();
		/////////////////////////////////////////////////////
		sample_scorefxn->set_weight( fa_rep, 0.12 ); // from RNA.
		sample_scorefxn->set_weight( linear_chainbreak, 0.0 ); // from RNA.
		sample_scorefxn->set_weight( chainbreak, 0.0 ); // from RNA.
		/////////////////////////////////////////////////////
	} else { // this used to happen for proteins -- may decide to deprecate.
		sample_scorefxn = ScoreFunctionFactory::create_score_function( sample_weights );
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, sample_scorefxn );
		sample_scorefxn->set_weight( linear_chainbreak, 0.2 /*arbitrary*/ ); // will cause problem with RNA?!
		if ( options->mapfile_activated()  && sample_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) {
			sample_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
		}
	}
	if ( options->turn_off_rna_chem_map_during_optimize() )  sample_scorefxn->set_weight( rna_chem_map, 0.0 ); // Just for now...

	return sample_scorefxn;
}

////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionCOP
initialize_pack_scorefxn( core::scoring::ScoreFunctionCOP sample_scorefxn,
	pose::Pose const & /* put back in if we want to try contains_protein hack*/ ){


	// NOTE:
	// Modern solvation terms (i.e., geom_sol_fast, occ_sol_fitted, etc.) are not necessary for RNA packing.
	// Although important for modeling DNA/RNA base stacking, they are too slow right now.
	// In the case of RNA, packing only involves O2' hydroxyls; fa_sol is sufficient for such cases.
	// The fa_sol term contains both lk_polar and lk_nonpolar terms.
	// - Caleb Geniesse 12.08.2014

	using namespace core::scoring;

	// if ( !contains_protein( pose ) ) return initialize_o2prime_pack_scorefxn( sample_scorefxn );

	ScoreFunctionOP pack_scorefxn = sample_scorefxn->clone();
	methods::EnergyMethodOptions const & energy_method_options = sample_scorefxn->energy_method_options();

	// hack for speed -- geom_sol & lk_nonpolar are too slow right now.
	// [see also: O2PrimePacker]
	if ( sample_scorefxn->has_nonzero_weight( geom_sol ) ||
			 sample_scorefxn->has_nonzero_weight( geom_sol_fast ) ||
			 sample_scorefxn->has_nonzero_weight( lk_ball ) ||
			sample_scorefxn->has_nonzero_weight( occ_sol_fitted )||
			sample_scorefxn->has_nonzero_weight( lk_polar )      ) {

		Real lk_weight(  sample_scorefxn->get_weight( lk_nonpolar ) );
		if ( sample_scorefxn->has_nonzero_weight( lk_nonpolar ) ) {
			runtime_assert( !energy_method_options.etable_options().no_lk_polar_desolvation );
			runtime_assert( !sample_scorefxn->has_nonzero_weight( fa_sol ) );
		} else { // new -- not using lk_nonpolar -- instead using fa_sol with no_lk_polar_desolvation flags
			runtime_assert( energy_method_options.etable_options().no_lk_polar_desolvation );
			runtime_assert( sample_scorefxn->has_nonzero_weight( fa_sol ) );
			lk_weight = sample_scorefxn->get_weight( fa_sol );

			// get fa_sol to handle polars again -- for speed.
			etable::EtableOptionsOP pack_etable_options( new etable::EtableOptions( energy_method_options.etable_options() ) );
			pack_etable_options->no_lk_polar_desolvation  = true;
			methods::EnergyMethodOptionsOP pack_energy_method_options = energy_method_options.clone();
			pack_energy_method_options->etable_options( *pack_etable_options);
			pack_scorefxn->set_energy_method_options( *pack_energy_method_options );
		}

		pack_scorefxn->set_weight( geom_sol,      0.0 );
		pack_scorefxn->set_weight( geom_sol_fast, 0.0 );
		pack_scorefxn->set_weight( lk_ball,   0.0 );
		pack_scorefxn->set_weight( occ_sol_fitted, 0.0 );
		pack_scorefxn->set_weight( lk_polar,       0.0 );
		pack_scorefxn->set_weight( lk_nonpolar,    0.0 );
		pack_scorefxn->set_weight( fa_sol,         lk_weight );

	} else {
		runtime_assert( !pack_scorefxn->has_nonzero_weight( lk_nonpolar ) ); // only allow lk_nonpolar if with geom_sol or occ_sol_fitted.
	}

	if ( pack_scorefxn->has_nonzero_weight( free_dof ) ) {
		Real const free_dof_weight = pack_scorefxn->get_weight( free_dof );
		pack_scorefxn->set_weight( free_suite, free_dof_weight );
		pack_scorefxn->set_weight( free_2HOprime, free_dof_weight );
		pack_scorefxn->set_weight( free_side_chain, free_dof_weight );
		pack_scorefxn->set_weight( free_dof, 0.0 ); // do not calculate free_base
	}
	pack_scorefxn->set_weight( free_base, 0.0 );

	// these also take too long to compute.
	pack_scorefxn->set_weight( fa_stack, 0.0 );
	pack_scorefxn->set_weight( ch_bond, 0.0 );
	pack_scorefxn->set_weight( rna_chem_map, 0.0 ); // Just for now...

	return pack_scorefxn;
}

////////////////////////////////////////////////////////////////////////
// may be deprecated in favor of initialize_pack_scorefxn above.
core::scoring::ScoreFunctionCOP
initialize_o2prime_pack_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){

	using namespace core::scoring;

	ScoreFunctionOP o2prime_pack_scorefxn( new ScoreFunction );

	// Each of the following terms have been pretty optimized for the packer (trie, etc.)
	o2prime_pack_scorefxn->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
	o2prime_pack_scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
	o2prime_pack_scorefxn->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
	o2prime_pack_scorefxn->set_weight( free_2HOprime, scorefxn->get_weight( free_2HOprime ) );
	///Warning, don't include hbond_intra, since hbond_intra HAS NOT been been optimized for packing!


	if ( scorefxn->has_nonzero_weight( lk_nonpolar ) ) {
		//// note that geom_sol is not optimized well --> replace with lk_sol for now.
		o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) );
	} else {
		runtime_assert( scorefxn->has_nonzero_weight( fa_sol ) );
		o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( fa_sol ) );
	}
	//This sets NO_HB_ENV_DEP, INCLUDE_INTRA_RES_RNA_HB and etcs.
	o2prime_pack_scorefxn->set_energy_method_options( scorefxn->energy_method_options() );

	return o2prime_pack_scorefxn;
}

} //modeler
} //stepwise
} //protocols
