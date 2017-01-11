// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/PartitionContactScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/PartitionContactScreener.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.PartitionContactScreener" );

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// In StepWise Assembly enumeration, a massive acceleration come from filtering out poses that either have too many clashes
//  or are not making any contacts between partitions whose relative orientations are being sampled.
//
// This screener replaces ProteinAtrRepScreener and RNA_AtrRepScreener, which used Rosetta's scorefunction and comparison
//  of atr/rep after a move, compared to values for a pose 'split' at the moving DOFs.  Now, we just explicitly
//  calculate atr/rep between residues in the two partitions separated by the moving DOFS.
//
// -- rhiju, 2014
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
PartitionContactScreener::PartitionContactScreener( pose::Pose const & pose,
	modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
	bool const use_loose_rep_cutoff,
	core::scoring::methods::EnergyMethodOptions const & options /* how to setup etable */ ):
	pose_( pose ),
	moving_res_list_( working_parameters->working_moving_res_list() ),
	fa_atr_weight_( 0.23 ), // historical -- used to use Rosetta scoring
	fa_rep_weight_( 0.12 ), // historical -- used to use Rosetta scoring
	rep_cutoff_(  4.0 ), // 10.0 for proteins?
	atr_cutoff_( -1.0 ), // -0.4 for proteins? [allow a single H-bond to count]
	use_loose_rep_cutoff_( use_loose_rep_cutoff ),
	close_chain_( modeler::figure_out_moving_cutpoints_closed_from_moving_res( pose, moving_res_list_ ).size() > 0 ),
	actual_rep_cutoff_( rep_cutoff_ )
{
	runtime_assert( moving_res_list_.size() > 0 );
	initialize_actual_rep_cutoff();
	initialize_evaluator( options );
}

//Destructor
PartitionContactScreener::~PartitionContactScreener()
{}

//////////////////////////////////////////////////////////////////////////////
void
PartitionContactScreener::initialize_actual_rep_cutoff(){
	actual_rep_cutoff_ = rep_cutoff_;
	if ( use_loose_rep_cutoff_ ) {
		actual_rep_cutoff_ = 200.0; // KIC needs a much higher cutoff -- atoms can get really close
	} else if ( close_chain_ ) {
		actual_rep_cutoff_ = 10.0;  //Parin's old parameter
	}
}

//////////////////////////////////////////////////////////////////////////////
void
PartitionContactScreener::initialize_evaluator( core::scoring::methods::EnergyMethodOptions const & options ){
	core::scoring::etable::EtableCOP etable(core::scoring::ScoringManager::get_instance()->etable( options ) );
	eval_ = core::scoring::etable::AnalyticEtableEvaluatorOP( new core::scoring::etable::AnalyticEtableEvaluator( *etable ) );
}

//////////////////////////////////////////////////////////////////////////////
bool
PartitionContactScreener::check_screen(){
	bool makes_contact( false ), atr_ok( false ), rep_ok( false );
	for ( Size const res : moving_res_list_ ) {
		check_screen( res, atr_ok, rep_ok );
		if ( !rep_ok ) return false;
		if (  atr_ok ) makes_contact = true;
	}
	return makes_contact;
}

//////////////////////////////////////////////////////////////////////////////
//
// Look at residue pair energy for residues in partitions on either side of
//  moving connection (as defined by moving_res).
//
// Could be sped up a little if we precalculate neighbor graph and/or cache
// energies between different moving_res.
//
//////////////////////////////////////////////////////////////////////////////
void
PartitionContactScreener::check_screen( Size const moving_res, bool & atr_ok, bool & rep_ok ) const {

	using namespace core::scoring;
	using namespace core::scoring::etable;
	using namespace core::scoring::etable::count_pair;

	atr_ok = false;
	rep_ok = false;

	utility::vector1< Size > root_partition_res, moving_partition_res;
	modeler::figure_out_root_and_moving_partition_res( pose_, moving_res, root_partition_res, moving_partition_res );

	EnergyMap emap;
	//Distance const contact_dist( 4.0 );
	for ( Size const i : root_partition_res ) {
		core::conformation::Residue const & rsd_i = pose_.residue( i );

		for ( Size const j : moving_partition_res ) {
			core::conformation::Residue const & rsd_j = pose_.residue( j );

			CPCrossoverBehavior crossover = CP_CROSSOVER_3;
			if ( rsd_i.polymeric_sequence_distance( rsd_j ) == 1 ) crossover = CP_CROSSOVER_4;
			CountPairFunctionCOP count_pair =  CountPairFactory::create_count_pair_function( rsd_i, rsd_j, crossover );

			eval_->residue_atom_pair_energy( rsd_i, rsd_j, *count_pair, emap );
			if ( ( fa_rep_weight_ * emap[ fa_rep ] ) > actual_rep_cutoff_ ) {
				return;
			}
		}
	}
	rep_ok = true;

	// AMW TODO: instead of multiplying by weights then comparing to x
	// just compare to a precalculated value, x/weight...
	Real const weighted_atr_score = fa_atr_weight_ * emap[ fa_atr ];
	Real const weighted_rep_score = fa_rep_weight_ * emap[ fa_rep ];
	Real const weighted_contact_score = weighted_atr_score + weighted_rep_score;
	if ( close_chain_ ) { // no check on attractive.
		atr_ok = true;
	} else if ( use_loose_rep_cutoff_ ) {
		if ( weighted_atr_score < atr_cutoff_ &&
				weighted_contact_score < (actual_rep_cutoff_ - rep_cutoff_)  )  atr_ok = true;
	} else {
		if ( weighted_atr_score < atr_cutoff_ &&
				weighted_contact_score < 0.0 )  atr_ok = true;
	}

	return;

}


} //screener
} //stepwise
} //protocols
