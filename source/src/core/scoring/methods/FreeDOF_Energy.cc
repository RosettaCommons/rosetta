// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FreeDOF_Energy.cc
/// @brief  FreeMoiety energy method implementation
/// @author Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <core/scoring/methods/FreeDOF_Energy.hh>
#include <core/scoring/methods/FreeDOF_EnergyCreator.hh>
#include <core/scoring/methods/FreeDOF_Options.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Package Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

using namespace core::pose::full_model_info;

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.FreeDOF_Energy", basic::t_info );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the FreeDOF_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
FreeDOF_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new FreeDOF_Energy( options ) );
}

ScoreTypes
FreeDOF_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( free_suite );
	sts.push_back( free_2HOprime );
	sts.push_back( free_side_chain );
	sts.push_back( free_base );
	sts.push_back( free_res );
	sts.push_back( free_dof );
	return sts;
}

/// ctor
FreeDOF_Energy::FreeDOF_Energy( methods::EnergyMethodOptions const & energy_method_options ) :
	parent( methods::EnergyMethodCreatorOP( new FreeDOF_EnergyCreator ) ),
	energy_method_options_( energy_method_options ),
	options_( energy_method_options.free_dof_options() ),
	free_res_weights_( energy_method_options_.method_weights( free_res ) )
{
}

FreeDOF_Energy::~FreeDOF_Energy() {}

/// clone
core::scoring::methods::EnergyMethodOP
FreeDOF_Energy::clone() const
{
	return core::scoring::methods::EnergyMethodOP( new FreeDOF_Energy( energy_method_options_ ) );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

void
FreeDOF_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const{
	nonconst_full_model_info( pose ); // does the setup.
}

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
FreeDOF_Energy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();
	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();

	Real free_suite_energy( 0.0 ), free_2HOprime_energy( 0.0 ), free_side_chain_energy( 0.0 );

	if ( rsd.has_variant_type( chemical::VIRTUAL_PHOSPHATE ) ) {
		Size const seqpos_in_full_model = res_list[ rsd.seqpos() ];
		if ( seqpos_in_full_model > 1 && !cutpoint_open_in_full_model.has_value( seqpos_in_full_model - 1 ) ) {
			free_suite_energy += options_.free_suite_bonus();
		}
	}

	// Note: In following, bonus/penalties for phosphates & sugars were folded into free_suite_energy
	// to have concise outfiles. That historical choice was confusing. Its now better to fold all into free_dof, and
	// allowing separate bonuses to be set through tags like FREE_SUGAR_BONUS in .wts file, or set from command-line.
	if ( rsd.has_variant_type( chemical::FIVE_PRIME_PACKABLE_PHOSPHATE ) ) { // this always comes with a virtual phosphate
		free_suite_energy += options_.pack_phosphate_penalty();
	}
	if ( rsd.has_variant_type( chemical::THREE_PRIME_PACKABLE_PHOSPHATE ) ) {
		free_suite_energy += options_.pack_phosphate_penalty();
	}
	if ( rsd.has_variant_type( chemical::FIVE_PRIME_PHOSPHATE ) ) {
		// weird, I know. these variants should go on as intermediates in stepwise assembly, replacing
		// virtual phosphates. Ideally the bonus & the penalty should cancel, but somehow that is too stringent
		// and would end up disallowing phosphates from ever being instantiated.
		free_suite_energy += options_.free_suite_bonus() + options_.pack_phosphate_penalty();
	}
	if ( rsd.has_variant_type( chemical::THREE_PRIME_PHOSPHATE ) ) {
		free_suite_energy += options_.pack_phosphate_penalty();
	}
	if ( rsd.has_variant_type( chemical::VIRTUAL_RIBOSE ) ) {
		free_suite_energy += options_.free_sugar_bonus();
	}
	// end block of terms that historically were folded into free_suite score term.

	if ( rsd.has_variant_type( chemical::VIRTUAL_O2PRIME_HYDROGEN ) ) free_2HOprime_energy += options_.free_2HOprime_bonus();

	// Following could be more fine-grained if we allow for virtualization of 'tips' of side chains. Currently
	//  all-or-nothing, though.
	if ( rsd.has_variant_type( chemical::VIRTUAL_SIDE_CHAIN ) )  free_side_chain_energy += options_.free_side_chain_bonus() * rsd.nchi() ;

	emap[ free_suite      ] += free_suite_energy;
	emap[ free_2HOprime   ] += free_2HOprime_energy;
	emap[ free_side_chain ] += free_side_chain_energy;
	emap[ free_dof        ] += free_suite_energy + free_2HOprime_energy + free_side_chain_energy;
}

//////////////////////////////////////////////////////////////////////////////////////////
void
FreeDOF_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & totals
) const
{

	/////////////////
	// free_res
	/////////////////
	// this should only be computed once for a full_model run, and not for the other poses.
	// (it looks over all the other_pose's). The good news is that when the other poses are scored
	// by OtherPoseEnergy, the weight on free_res is set to zero, so it won't be included here.
	if ( scorefxn.has_nonzero_weight( free_res ) ) {
		using namespace core::chemical;
		utility::vector1< char > missing_residues;
		pose::full_model_info::get_number_missing_residues_and_connections( pose, missing_residues );
		Real free_res_energy( 0.0 );
		for ( Size n = 1; n <= missing_residues.size(); n++ ) {
			if ( !oneletter_code_specifies_aa( missing_residues[ n ] ) ) continue;
			AA const aa = aa_from_oneletter_code( missing_residues[n] );
			if ( Size( aa ) > free_res_weights_.size() ) continue;
			free_res_energy += free_res_weights_[ aa  ];
		}
		totals[ free_res ] += free_res_energy;
	}

	///////////////////////////////////////////////////////////////////////
	// entropy terms -- cost of base stacking & hydrogen bonding.
	// set as constants, chosen initially as guesses.
	///////////////////////////////////////////////////////////////////////
	if ( !scorefxn.has_nonzero_weight( free_base ) && !scorefxn.has_nonzero_weight( free_dof ) ) return;

	utility::vector1< Real > base_stack_energy, base_hbond_energy, sugar_hbond_energy;
	accumulate_stack_energy( pose, scorefxn, base_stack_energy );
	get_hbond_energy( pose, scorefxn, base_hbond_energy, sugar_hbond_energy );

	static Real const base_stack_energy_cutoff_( -0.25 );
	static Real const base_hbond_energy_cutoff_( -0.25 );
	static Real const sugar_hbond_energy_cutoff_( -0.20 );

	Real total_free_base_score( 0.0 );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !pose.residue( n ).is_RNA() ) continue;
		Real free_base_score( 0.0 );
		if ( base_hbond_energy[ n ] <= base_hbond_energy_cutoff_ ) {
			// assume that a base hydrogen bond locks the residue into place.
			continue;
		} else {
			if ( base_stack_energy[ n ] <= base_stack_energy_cutoff_ ) {
				free_base_score += -1.0; // make this an option
			} else {
				// no stacks - big bonus!
				free_base_score += -2.0; // make this an option
			}
			if ( sugar_hbond_energy[ n ] <= sugar_hbond_energy_cutoff_ ) {
				free_base_score += 1.0;
			}
		}
		total_free_base_score += free_base_score;
	}

	totals[ free_base ] += total_free_base_score;
}

///////////////////////////////////////////////////////////////////////
void
FreeDOF_Energy::do_fa_stack_scorefunction_checks( ScoreFunction const & scorefxn ) const {
	EnergyMap const & weights( scorefxn.weights() );
	runtime_assert( scorefxn.has_nonzero_weight( fa_stack ) ||
		( scorefxn.has_nonzero_weight( fa_stack_lower ) &&
		weights[ fa_stack_upper ] == weights[ fa_stack_lower ] ) );

	if ( !scorefxn.check_methods_in_right_order( fa_stack, free_base ) ) {
		utility_exit_with_message( "fa_stack should be before free_base/free_dof in .wts file!" );
	}

}


///////////////////////////////////////////////////////////////////////
// makes use of stack_energy.
///////////////////////////////////////////////////////////////////////
void
FreeDOF_Energy::accumulate_stack_energy(
	pose::Pose & pose,
	ScoreFunction const & scorefxn,
	utility::vector1< Real > & stack_energy
) const
{
	using namespace core::scoring;

	stack_energy = utility::vector1< Real >( pose.size(), 0.0 );
	do_fa_stack_scorefunction_checks( scorefxn );

	EnergyGraph const & energy_graph = pose.energies().energy_graph();
	EnergyMap const & weights( scorefxn.weights() );

	Real fa_stack_weight( weights[ fa_stack ] );
	if ( fa_stack_weight == 0.0 ) fa_stack_weight = weights[ fa_stack_upper ];

	// accumulate energies, but keep track of which residue supplied the base.
	for ( Size i=1, i_end = pose.size(); i<= i_end; ++i ) {
		if ( !pose.residue( i ).is_RNA() ) continue;
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const j( edge.get_second_node_ind() );
			runtime_assert( i < j );
			stack_energy[ i ] += fa_stack_weight * edge[ fa_stack_lower ];
			stack_energy[ j ] += fa_stack_weight * edge[ fa_stack_upper ];
		}
	}
}


///////////////////////////////////////////////////////////////////////
void
FreeDOF_Energy::get_hbond_energy(
	pose::Pose & pose,
	ScoreFunction const & scorefxn,
	utility::vector1< Real > & base_hbond_energy,
	utility::vector1< Real > & sugar_hbond_energy
) const
{
	using namespace core::scoring::hbonds;

	EnergyMap const & weights( scorefxn.weights() );

	// unclear whether we can avoid repeating work elsewhere through caching this.
	HBondSet hbond_set;
	fill_hbond_set( pose, false /*calculate derivative*/, hbond_set );

	base_hbond_energy= utility::vector1< Real > ( pose.size(), 0.0 );
	sugar_hbond_energy= utility::vector1< Real > ( pose.size(), 0.0 );

	for ( Size n = 1; n <= hbond_set.nhbonds(); ++n ) {

		if ( !hbond_set.allow_hbond(n) ) continue;
		HBond const & hbond_(hbond_set.hbond(n));
		Size const i = hbond_.don_res();
		Size const j = hbond_.acc_res();

		HBEvalType const hbe_type = hbond_.eval_type();
		Real const scorefxn_weight = hb_eval_type_weight( hbe_type, weights, (i == j), hbond_set.hbond_options().put_intra_into_total() );
		Real hbE = hbond_.energy() /*raw energy*/ * hbond_.weight() /*env-dep-wt*/ * scorefxn_weight;

		// note double-counting here if hbond is base/base.
		if ( pose.residue_type( i ).is_RNA() && pose.residue_type( i ).RNA_type().is_RNA_base_atom( hbond_.don_hatm() ) ) {
			base_hbond_energy[ i ] += hbE;
		}
		if ( pose.residue_type( j ).is_RNA() && pose.residue_type( j ).RNA_type().is_RNA_base_atom( hbond_.acc_atm() ) ) {
			base_hbond_energy[ j ] += hbE;
		}

		if ( pose.residue_type( i ).is_RNA() && ( hbond_.don_hatm() == pose.residue_type( i ).RNA_type().ho2prime_index() ) ) {
			sugar_hbond_energy[ i ] += hbE;
		}
		if ( pose.residue_type( j ).is_RNA() && ( hbond_.acc_atm() == pose.residue_type( j ).RNA_type().o2prime_index() ) ) {
			sugar_hbond_energy[ j ] += hbE;
		}
	}
}

/// @brief FreeDOF_Energy is context independent; indicates that no context graphs are required
void
FreeDOF_Energy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
FreeDOF_Energy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

