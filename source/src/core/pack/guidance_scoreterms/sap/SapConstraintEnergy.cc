// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/sap/SapConstraintEnergy.cc
/// @brief Energy method that enforces the sap_constraint
/// @details sap_constraint constrains your protein to be soluble.
/// @author Brian Coventry (bcov@uw.edu)

// Unit headers
#include <core/pack/guidance_scoreterms/sap/SapConstraintEnergy.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintHelper.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraintOptions.hh>
#include <core/pack/guidance_scoreterms/sap/SapConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/SapMathConstraint.hh>
#include <core/pack/guidance_scoreterms/sap/util.hh>

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

using namespace core::scoring;

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace sap {

static basic::Tracer TR("core.pack.guidance_scoreterms.sap.SapConstraintEnergy");


core::scoring::methods::EnergyMethodOP
SapConstraintEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return utility::pointer::make_shared<SapConstraintEnergy>( options );
}

ScoreTypes
SapConstraintEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( scoring::sap_constraint );
	return sts;
}

SapConstraintEnergy::SapConstraintEnergy ( core::scoring::methods::EnergyMethodOptions const & ) :
	parent1( utility::pointer::make_shared<SapConstraintEnergyCreator>() ),
	parent2( ),
	disabled_(false),
	helpers_()
{}


core::scoring::methods::EnergyMethodOP SapConstraintEnergy::clone() const {
	return utility::pointer::make_shared< SapConstraintEnergy >( *this );
}

void SapConstraintEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	return;
}

core::Size SapConstraintEnergy::version() const
{
	return 1; // Initial versioning
}

void
SapConstraintEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( disabled_ ) return; //Do nothing when this energy is disabled.

	utility::vector1<SapConstraintHelperOP> helpers = get_helpers_from_pose( pose );

	utility::vector1< std::pair< utility::vector1< std::pair< Real, SapConstraintHelperCOP > >, SapMathConstraintCOP > > math_csts =
		get_math_csts_from_pose( pose, helpers );

	utility::vector1< core::conformation::ResidueCOP > res_vector;
	pack::rotamer_set::RotamerSetsOP rotsets = rotamer_sets_from_pose( pose, res_vector );

	// We're using the helpers here to calculate sap correctly. The reason we don't use calculate_sap() here is because
	//  we need to helpers to have their values calculated correctly for the math constraints
	for ( SapConstraintHelperOP const & helper : helpers ) {
		helper->init_with_pose( pose, *rotsets );
		Real sap = helper->calculate_energy( res_vector, 0 );
		if ( helper->options()->full_accuracy_when_scoring() ) {
			sap = helper->set_accurate_sasa_and_recalc( pose );
		}
		emap[ sap_constraint ] += helper->options()->transform_sap_to_score( sap, false );
	}

	for ( auto const & math_cst : math_csts ) {
		emap[ sap_constraint ] += math_cst.second->get_score( math_cst.first );
	}

}

core::Real
SapConstraintEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	utility::vector1< core::Size > const &,
	core::Size const substitution_position
) const {
	if ( disabled_ ) return 0.0;

	Real score = 0;

	for ( SapConstraintHelperOP const & helper : helpers_ ) {
		Real sap = helper->calculate_energy( resvect, substitution_position );
		score += helper->options()->transform_sap_to_score( sap, true );
	}

	for ( auto const & math_cst : math_csts_ ) {
		score += math_cst.second->get_score( math_cst.first );
	}

	return score;
}

void
SapConstraintEnergy::commit_considered_substitution() {
	for ( SapConstraintHelperOP const & helper : helpers_ ) {
		helper->commit_considered_substitution();
	}
}



/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
SapConstraintEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose &pose,
	core::pack::rotamer_set::RotamerSets const & rotamer_sets,
	core::scoring::ScoreFunction const &/*sfxn*/
) {
	disabled_ = false; //Enabled for packing.

	helpers_ = get_helpers_from_pose( pose );

	math_csts_ = get_math_csts_from_pose( pose, helpers_ );

	for ( SapConstraintHelperOP const & helper : helpers_ ) {
		helper->init_with_pose( pose, rotamer_sets );
	}
}

/// @brief Clear the cached data from the pose after packing.
///
void
SapConstraintEnergy::clean_up_residuearrayannealableenergy_after_packing(
	core::pose::Pose & pose
) {
	for ( auto & helper : helpers_ ) {
		Real actual_sap = calculate_sap( pose, helper->options()->score_selector(), helper->options()->sap_calculate_selector(),
			helper->options()->sasa_selector() );
		helper->report_final_score( actual_sap );
	}
	helpers_.clear();
	math_csts_.clear();
}

/// @brief Disable this energy during minimization.
void
SapConstraintEnergy::setup_for_minimizing(
	pose::Pose & /*pose*/, ScoreFunction const & /*sfxn*/,
	kinematics::MinimizerMapBase const & /*minmap*/
) const {
	TR << "Disabling SapConstraintEnergy during minimization." << std::endl;
	disabled_ = true;
}

/// @brief Re-enable this energy after minimization.
void
SapConstraintEnergy::finalize_after_minimizing( pose::Pose & /*pose*/ ) const {
	TR << "Re-enabling SapConstraintEnergy following minimization." << std::endl;
	disabled_ = false;
}

utility::vector1<SapConstraintHelperOP>
SapConstraintEnergy::get_helpers_from_pose(
	core::pose::Pose const &pose
) const {
	utility::vector1<SapConstraintHelperOP> helpers;

	core::Size const n_sequence_constraints( pose.constraint_set()->n_sequence_constraints() );
	if ( n_sequence_constraints > 0 ) {
		for ( core::Size i=1; i<=n_sequence_constraints; ++i ) {
			SapConstraintCOP cur_cst( utility::pointer::dynamic_pointer_cast<SapConstraint const>(
				pose.constraint_set()->sequence_constraint(i) ) );
			if ( !cur_cst ) continue; //Continue if this isn't a SapConstraint.

			SapConstraintHelperOP helper = utility::pointer::make_shared<SapConstraintHelper>( cur_cst->get_const_options() );

			helpers.push_back( helper );
		}
	}
	return helpers;
}

utility::vector1< std::pair< utility::vector1< std::pair< Real, SapConstraintHelperCOP > >, SapMathConstraintCOP > >
SapConstraintEnergy::get_math_csts_from_pose(
	core::pose::Pose const & pose,
	utility::vector1<SapConstraintHelperCOP> const & helpers
) const {
	utility::vector1< std::pair< utility::vector1< std::pair< Real, SapConstraintHelperCOP > >, SapMathConstraintCOP > > math_csts;

	core::Size const n_sequence_constraints( pose.constraint_set()->n_sequence_constraints() );
	if ( n_sequence_constraints > 0 ) {
		for ( core::Size i=1; i<=n_sequence_constraints; ++i ) {
			SapMathConstraintCOP cur_cst( utility::pointer::dynamic_pointer_cast<SapMathConstraint const>(
				pose.constraint_set()->sequence_constraint(i) ) );
			if ( !cur_cst ) continue; //Continue if this isn't a SapMathConstraint.

			math_csts.emplace_back( cur_cst->parse_helpers( helpers ), cur_cst );
		}
	}
	return math_csts;
}


} //sap
} //guidance_scoreterms
} //pack
} //core
