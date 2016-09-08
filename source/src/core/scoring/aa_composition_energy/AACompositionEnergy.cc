// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/aa_composition_energy/AACompositionEnergy.cc
/// @brief An EnergyMethod that penalizes deviation from a desired amino acid composition.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.  It has also been modified to permit sub-regions of a pose to be scored.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/scoring/aa_composition_energy/AACompositionEnergy.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergyCreator.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <utility/numbers.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/scoring/aa_composition_energy/AACompositionConstraint.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

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
namespace aa_composition_energy {

static THREAD_LOCAL basic::Tracer TR("core.scoring.aa_composition_energy.AACompositionEnergy");

/// @brief This must return a fresh instance of the AACompositionEnergy class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
AACompositionEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return core::scoring::methods::EnergyMethodOP( new AACompositionEnergy( options ) );
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
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) ),
	parent2( ),
	setup_helpers_(),
	setup_helpers_for_packing_(),
	setup_helper_masks_for_packing_()
{
	//The following reads from disk the first time only, and caches the data in memory:
	setup_helpers_ = core::scoring::ScoringManager::get_instance()->get_cloned_aa_comp_setup_helpers( options );
	if ( TR.Debug.visible() ) report();
}

/// @brief Copy constructor.
///
AACompositionEnergy::AACompositionEnergy( AACompositionEnergy const &src ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) ),
	parent2( src ),
	setup_helpers_(), //CLONE the helper data below; don't copy them.
	setup_helpers_for_packing_(), //CLONE these below, too -- don't copy them.
	setup_helper_masks_for_packing_(  src.setup_helper_masks_for_packing_ )
{
	for ( core::Size i=1, imax=src.setup_helpers_.size(); i<=imax; ++i ) {
		setup_helpers_.push_back( src.setup_helpers_[i]->clone() );
	}
	for ( core::Size i=1, imax=src.setup_helpers_for_packing_.size(); i<=imax; ++i ) {
		setup_helpers_for_packing_.push_back( src.setup_helpers_for_packing_[i]->clone() );
	}
}

/// @brief Default destructor.
///
AACompositionEnergy::~AACompositionEnergy() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP AACompositionEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new AACompositionEnergy(*this) );
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
	core::Size const nres( pose.size() );

	//Vector of residue owning pointers:
	utility::vector1< core::conformation::ResidueCOP > resvector;
	resvector.reserve(nres);

	//Populate the vector with const owning pointers to the residues:
	for ( core::Size ir=1; ir<=nres; ++ir ) {
		resvector.push_back( pose.residue(ir).get_self_ptr() );
	}

	//Get the AACompositionEnergySetup objects from the pose and append them to the setup_helpers_ list, making a new setup_helpers list:
	utility::vector1< AACompositionEnergySetupCOP > setup_helpers;
	utility::vector1< core::select::residue_selector::ResidueSubset > masks;
	get_helpers_from_pose( pose, setup_helpers, masks ); //Pulls AACompositionEnergySetupCOPs from pose; generates masks from ResidueSelectors simultaneously.
	runtime_assert( masks.size() == setup_helpers.size() ); //Should be guaranteed to be true.

	totals[ aa_composition ] += calculate_energy( resvector, setup_helpers, masks ); //Using the vector of residue owning pointers, calculate the repeat energy (unweighted) and set the aa_composition to this value.

	return;
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.  Requires
/// that setup_residuearrayannealablenergy_for_packing() be called first.
core::Real
AACompositionEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect
) const {
	//for(core::Size i=1, imax=resvect.size(); i<=imax; ++i) { TR << i << ":\t" << resvect[i]->name() << "\t" << resvect[i]->seqpos() << std::endl; } //DELETE ME
	return calculate_energy( resvect, setup_helpers_for_packing_, setup_helper_masks_for_packing_);
}


/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called by finalize_total_energy() and during packer runs.
core::Real
AACompositionEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	utility::vector1< AACompositionEnergySetupCOP > const &setup_helpers,
	utility::vector1< core::select::residue_selector::ResidueSubset > const &masks
) const
{
	runtime_assert( setup_helpers.size() == masks.size() ); //Should always be true.

	core::Real outer_accumulator(0.0);

	for ( core::Size ihelper=1, ihelpermax=setup_helpers.size(); ihelper<=ihelpermax; ++ihelper ) { //loop through all setup_helper objects

		// Const owning pointer to the setup helper:
		AACompositionEnergySetupCOP helper( setup_helpers[ihelper] );
		if ( helper->n_property_sets() == 0 ) return 0.0; //Return 0 if we're not tracking any properties.

		// Number of residues:
		core::Size const nres( resvect.size() );

		// Number of unmasked residues for the current helper:
		core::Size n_res_total(0); //Only computed if fractional counts are used.
		if ( helper->use_fract_expected_values_any() || helper->use_fract_ranges_any() ) {
			for ( core::Size j=1, jmax=masks[ihelper].size(); j<=jmax; ++j ) { if ( masks[ihelper][j] ) ++n_res_total; } //Count the number of unmasked residues for the current AACompositionEnergySetup object.
		}

		// Figure out expected number of residues of each property:
		utility::vector1 < signed long > expected;
		expected.resize( helper->n_property_sets(), 0 );
		for ( core::Size i=1, imax=expected.size(); i<=imax; ++i ) {
			if ( helper->expected_by_properties_absolute(i) >= 0 ) expected[i] = helper->expected_by_properties_absolute(i);
			else {
				expected[i] = static_cast< signed long >( utility::round( static_cast<double>( n_res_total ) * static_cast<double>( helper->expected_by_properties_fraction(i) ) ) );
			}
		}

		// Counts of the residues corresponding to the various property sets that we're counting.  The indicies are property set indices that are internal to the setup_helper object.
		utility::vector1 < signed long > counts;
		utility::vector1 < core::Real > fract_counts;
		counts.resize( helper->n_property_sets(), 0 ); //Make sure that we have a counter for each property set that we're counting.
		fract_counts.resize( helper->n_property_sets(), 0.0 ); //Make sure that we have a counter for each property set that we're counting.

		// A storage container for the property sets that each residue is in.
		utility::vector1< core::Size > this_rsd_property_sets;

		// Loop through all residues and count residues in each property set.
		for ( core::Size ir=1; ir<=nres; ++ir ) {
			if ( !masks[ihelper][ir] ) continue; //Skip residues that are masked.
			helper->property_set_indices_matching_residue(resvect[ir]->type(), this_rsd_property_sets); //Get the property sets that this residue belongs to.
			for ( core::Size j=1, jmax=this_rsd_property_sets.size(); j<=jmax; ++j ) {
				++counts[this_rsd_property_sets[j]]; //Increment the counter for EVERY property set that this residue belongs to.
			}
		}

		// Accumulator for penalty function
		core::Real accumulator(0.0);
		for ( core::Size i=1, imax=counts.size(); i<=imax; ++i ) { //Loop through the counts and accumulate the appropriate penalty
			counts[i] -= expected[i]; // Calculate the DELTA count.
			if ( helper->use_fract_ranges(i) ) {
				fract_counts[i] = ( n_res_total > 0 ? static_cast<core::Real>(counts[i]) / static_cast<core::Real>(n_res_total) : 0 );
				accumulator += helper->fract_property_penalty( fract_counts[i], i );
			} else {
				accumulator += helper->property_penalty( counts[i], i );
			}
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

/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
AACompositionEnergy::setup_residuearrayannealableenergy_for_packing (
	core::pose::Pose const &pose,
	core::scoring::ScoreFunction const &/*sfxn*/
) {
	if ( TR.Debug.visible() ) TR.Debug << "Setting up the AACompostionEnergy for packing." << std::endl;

	get_helpers_from_pose( pose, setup_helpers_for_packing_, setup_helper_masks_for_packing_ );

	runtime_assert( setup_helper_masks_for_packing_.size() == setup_helpers_for_packing_.size() ); //Should be guaranteed to be true.

	if ( TR.Debug.visible() ) {
		for ( core::Size i=1, imax=setup_helper_masks_for_packing_.size(); i<=imax; ++i ) {
			TR.Debug << "Setup helper " << i << ":" << std::endl;
			TR.Debug << "Mask: ";
			for ( core::Size j=1, jmax=setup_helper_masks_for_packing_[i].size(); j<=jmax; ++j ) {
				TR.Debug << (setup_helper_masks_for_packing_[i][j] ? "1" : "0");
			}
			TR.Debug << std::endl;
			TR.Debug << setup_helpers_for_packing_[i]->report();
		}
		TR.Debug.flush();
	}
	return;
}

/// @brief Given a pose, pull out the AACompositionEnergySetup objects stored in SequenceConstraints in the pose and
/// append them to the setup_helpers_ vector, returning a new vector.  This also generates a vector of masks simultaneously.
/// @param [in] pose The pose from which the AACompositionEnergySetupCOPs will be extracted.
/// @param [out] setup_helpers The output vector of AACompositionEnergySetupCOPs that is the concatenation of those stored in setup_helpers_ and those from the pose.
/// @param [out] masks The output vector of ResidueSubsets, which will be equal in size to the helpers vector.
/// @details The output vectors are first cleared by this operation.
void
AACompositionEnergy::get_helpers_from_pose(
	core::pose::Pose const &pose,
	utility::vector1< AACompositionEnergySetupCOP > &setup_helpers,
	utility::vector1< core::select::residue_selector::ResidueSubset > &masks
) const {
	setup_helpers.clear();
	masks.clear();
	if ( setup_helpers_.size() > 0 ) {
		setup_helpers = setup_helpers_; //Copy the setup_helpers_ list.
		masks.resize( setup_helpers_.size(), core::select::residue_selector::ResidueSubset( pose.size(), true ) ); //All of the helpers in the setup_helpers_ list should be applied globally.
	}

	core::Size const n_sequence_constraints( pose.constraint_set()->n_sequence_constraints() );
	if ( n_sequence_constraints > 0 ) {
		for ( core::Size i=1; i<=n_sequence_constraints; ++i ) {
			AACompositionConstraintCOP cur_cst( utility::pointer::dynamic_pointer_cast<AACompositionConstraint const>( pose.constraint_set()->sequence_constraint(i) ) );
			if ( !cur_cst ) continue; //Continue if this isn't an AACompositionConstraint.
			setup_helpers.push_back( cur_cst->aa_composition_energy_setup() ); //Append the AACompositionEnergySetup object stored in the current sequence constraint to the list to be used.
			core::select::residue_selector::ResidueSelectorCOP selector( cur_cst->selector() ); //Get the ResidueSelector in the current sequence constraint object, if there is one.  (May be NULL).
			if ( selector ) { //If we have a ResidueSelector, generate a mask from the pose and store it in the masks list.
				masks.push_back( selector->apply( pose ) );
			} else { //If not, add an all-true mask
				masks.push_back( core::select::residue_selector::ResidueSubset( pose.size(), true ) );
			}
		}
	}

	runtime_assert( setup_helpers.size() == masks.size() ); //Should be guaranteed true.

	return;
}


} // aa_composition_energy
} // scoring
} // core
