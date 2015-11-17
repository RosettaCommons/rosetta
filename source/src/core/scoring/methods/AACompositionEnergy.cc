// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/AACompositionEnergy.cc
/// @brief An EnergyMethod that penalizes deviation from a desired amino acid composition.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/methods/AACompositionEnergy.hh>
#include <core/scoring/methods/AACompositionEnergyCreator.hh>
#include <core/scoring/methods/AACompositionEnergySetup.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <utility/numbers.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR("core.scoring.methods.AACompositionEnergy");

/// @brief This must return a fresh instance of the AACompositionEnergy class, never an instance already in use.
///
methods::EnergyMethodOP
AACompositionEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return methods::EnergyMethodOP( new AACompositionEnergy( options ) );
}

/// @brief Defines the score types that this energy method calculates.
///
ScoreTypes
AACompositionEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( aa_composition );
	return sts;
}

/// @brief Options constructor.
///
AACompositionEnergy::AACompositionEnergy ( core::scoring::methods::EnergyMethodOptions const &options ) :
	parent1( methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) ),
	parent2( ),
	setup_helpers_()
{
	core::Size const n_setup_helpers( options.aa_composition_setup_file_count() );
	for ( core::Size i=1; i<=n_setup_helpers; ++i ) {
		setup_helpers_.push_back( AACompositionEnergySetupOP( new AACompositionEnergySetup ) );
		setup_helpers_[i]->initialize_from_file( options.aa_composition_setup_file(i) );
	}
	if ( TR.Debug.visible() ) report();
}

/// @brief Copy constructor.
///
AACompositionEnergy::AACompositionEnergy( AACompositionEnergy const &src ) :
	parent1( methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) ),
	parent2( src ),
	setup_helpers_() //CLONE the helper data below; don't copy them.
{
	for ( core::Size i=1, imax=src.setup_helpers_.size(); i<=imax; ++i ) {
		setup_helpers_.push_back( src.setup_helpers_[i]->clone() );
	}
}

/// @brief Default destructor.
///
AACompositionEnergy::~AACompositionEnergy() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
EnergyMethodOP AACompositionEnergy::clone() const {
	return EnergyMethodOP( new AACompositionEnergy(*this) );
}

/// @brief AACompositionEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void AACompositionEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief AACompositionEnergy is version 1.0 right now.
///
core::Size AACompositionEnergy::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
void AACompositionEnergy::finalize_total_energy( core::pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const
{
	//Number of residues:
	core::Size const nres( pose.n_residue() );

	//Vector of residue owning pointers:
	utility::vector1< core::conformation::ResidueCOP > resvector;
	resvector.reserve(nres);

	//Populate the vector with const owning pointers to the residues:
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		resvector.push_back( pose.residue(ir).get_self_ptr() );
	}

	totals[ aa_composition ] += calculate_energy( resvector ); //Using the vector of residue owning pointers, calculate the repeat energy (unweighted) and set the aa_composition to this value.

	return;
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called by finalize_total_energy().  Calls calculate_energy_by_restype() and
/// calculate_energy_by_properties().
core::Real AACompositionEnergy::calculate_energy( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const
{
	return calculate_energy_by_restype( resvect ) + calculate_energy_by_properties( resvect );
}

/// @brief Calculate the total energy based on residue types, given a vector of const owning pointers to residues.
/// @details Called by calculate_energy().
core::Real AACompositionEnergy::calculate_energy_by_restype( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const {

	core::Real outer_accumulator(0.0);

	for ( core::Size ihelper=1, ihelpermax=setup_helper_count(); ihelper<=ihelpermax; ++ihelper ) { //loop through all setup_helper objects

		// Const owning pointer to the setup helper:
		AACompositionEnergySetupCOP helper( setup_helper_cop( ihelper ) );
		if ( helper->n_residue_types() == 0 ) return 0.0; //Return 0 if we're not tracking any residue types.

		// Number of residues:
		core::Size const nres( resvect.size() );

		// Figure out expected number of residues of each type:
		utility::vector1 < signed long > expected;
		expected.resize( helper->n_residue_types(), 0 );
		for ( core::Size i=1, imax=expected.size(); i<=imax; ++i ) {
			if ( helper->expected_by_type_absolute(i) >= 0 ) expected[i] = helper->expected_by_type_absolute(i);
			else expected[i] = static_cast< signed long >( utility::round( static_cast<double>( nres ) * static_cast<double>( helper->expected_by_type_fraction(i) ) ) );
		}

		// Counts of the various residue types that we're counting.  The indicies are determined by the res_type_index_mappings_ map in the setup_helper object.
		utility::vector1 < signed long > counts;
		counts.resize( helper->n_residue_types(), 0 ); //Make sure that we have a counter for each residue type that we're counting

		// Loop through all residues and count residue types.
		for ( core::Size ir=1; ir<=nres; ++ir ) {
			if ( helper->has_type( resvect[ir]->type().name3() ) ) { //If this residue's name3 is in the map, it's a type that we're counting.
				++counts[ helper->res_type_index( resvect[ir]->type().name3() ) ]; //So add 1 to the count for that residue type.
			}
		}

		// Accumulator for penalty function
		core::Real accumulator(0.0);
		for ( core::Size i=1, imax=counts.size(); i<=imax; ++i ) { //Loop through the counts and accumulate the appropriate penalty
			counts[i] -= expected[i]; // Calculate the DELTA count.
			accumulator += helper->type_penalty( counts[i], i );
		}
		outer_accumulator += accumulator;

	} //looping through all helpers

	return outer_accumulator;
}

/// @brief Calculate the total energy based on residue properties, given a vector of const owning pointers to residues.
/// @details Called by calculate_energy().
core::Real AACompositionEnergy::calculate_energy_by_properties( utility::vector1< core::conformation::ResidueCOP > const &resvect ) const {

	core::Real outer_accumulator(0.0);

	for ( core::Size ihelper=1, ihelpermax=setup_helper_count(); ihelper<=ihelpermax; ++ihelper ) { //loop through all setup_helper objects

		// Const owning pointer to the setup helper:
		AACompositionEnergySetupCOP helper( setup_helper_cop(ihelper) );
		if ( helper->n_property_sets() == 0 ) return 0.0; //Return 0 if we're not tracking any properties.

		// Number of residues:
		core::Size const nres( resvect.size() );

		// Figure out expected number of residues of each property:
		utility::vector1 < signed long > expected;
		expected.resize( helper->n_property_sets(), 0 );
		for ( core::Size i=1, imax=expected.size(); i<=imax; ++i ) {
			if ( helper->expected_by_properties_absolute(i) >= 0 ) expected[i] = helper->expected_by_properties_absolute(i);
			else {
				expected[i] = static_cast< signed long >( utility::round( static_cast<double>( nres ) * static_cast<double>( helper->expected_by_properties_fraction(i) ) ) );
			}
		}

		// Counts of the residues corresponding to the various property sets that we're counting.  The indicies are property set indices that are internal to the setup_helper object.
		utility::vector1 < signed long > counts;
		counts.resize( helper->n_property_sets(), 0 ); //Make sure that we have a counter for each residue type that we're counting

		// A storage container for the property sets that each residue is in.
		utility::vector1< core::Size > this_rsd_property_sets;

		// Loop through all residues and count residues in each property set.
		for ( core::Size ir=1; ir<=nres; ++ir ) {
			helper->property_set_indices_matching_residue(resvect[ir]->type(), this_rsd_property_sets); //Get the property sets that this residue belongs to.
			for ( core::Size j=1, jmax=this_rsd_property_sets.size(); j<=jmax; ++j ) {
				++counts[this_rsd_property_sets[j]]; //Increment the counter for EVERY property set that this residue belongs to.
			}
		}

		// Accumulator for penalty function
		core::Real accumulator(0.0);
		for ( core::Size i=1, imax=counts.size(); i<=imax; ++i ) { //Loop through the counts and accumulate the appropriate penalty
			counts[i] -= expected[i]; // Calculate the DELTA count.
			accumulator += helper->property_penalty( counts[i], i );
		}

		outer_accumulator += accumulator;
	} //End loop through helper objects

	return outer_accumulator;
}

/// @brief Get a summary of all loaded data.
///
void AACompositionEnergy::report() const {
	if ( !TR.Debug.visible() ) return; //Do nothing if I don't have a tracer.

	TR.Debug << std::endl << "Summary of data loaded by AACompositionEnergy object:" << std::endl;

	for ( core::Size i=1, imax=setup_helper_count(); i<=imax; ++i ) {
		TR.Debug << "AACompositionEnergySetup #" << i << ":" << std::endl;
		TR.Debug << setup_helper_cop(i)->report();
	}

	TR.Debug << std::endl;

	TR.Debug.flush();

	return;
}


} // methods
} // scoring
} // core
