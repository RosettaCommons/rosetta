// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/RamaMutationSelector.hh
/// @brief  Selects positions that would have a rama_prepro score below a given threshold IF mutated to a given residue type.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/RamaMutationSelector.hh>
#include <protocols/cyclic_peptide/RamaMutationSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/conformation/Residue.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "protocols.cyclic_peptide.RamaMutationSelector" );

namespace protocols {
namespace cyclic_peptide {

/// @brief Constructor.
///
RamaMutationSelector::RamaMutationSelector():
	ResidueSelector(),
	target_type_(""),
	score_threshold_(0.0),
	rama_prepro_multiplier_(0.45)
{
}

/// @brief Destructor.
///
RamaMutationSelector::~RamaMutationSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
RamaMutationSelector::RamaMutationSelector(RamaMutationSelector const & src):
	ResidueSelector( src ),
	target_type_(src.target_type_),
	score_threshold_(src.score_threshold_),
	rama_prepro_multiplier_(src.rama_prepro_multiplier_)
{
}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
RamaMutationSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new RamaMutationSelector(*this) );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
RamaMutationSelector::ResidueSubset
RamaMutationSelector::apply(
	core::pose::Pose const &pose
) const {
	core::Size const nres( pose.total_residue() ); //Total residue count in pose.
	bool const has_target_type( !target_type().empty() ); //Is there a target type to which we're mutating?
	ResidueSubset selection( nres, false );

	core::chemical::ResidueTypeSetCOP rsd_type_set( pose.residue_type_set_for_pose( /*any type -- centroid or fullatom*/ ) );
	core::scoring::ScoringManager const &score_man( *(core::scoring::ScoringManager::get_instance()) );
	core::scoring::RamaPrePro const &rama( score_man.get_RamaPrePro() );
	utility::vector1< core::Real > gradient; //Dummmy var needed for RamaPrePro evaluations.

	for ( core::Size ir(1); ir<=nres; ++ir ) { //Loop through all residues in the pose.
		if ( !pose.residue(ir).is_polymer() ) {
			TR.Warning << "residue " << ir << " is not polymeric!  Skipping and not selecting." << std::endl;
			continue; //Skip non-polymeric residues.
		}
		if ( is_terminus(pose.residue(ir)) ) {
			TR.Warning << "residue " << ir << " is a terminus!  Skipping and not selecting." << std::endl;
			continue; //Skip terminal residues.
		}

		core::chemical::ResidueTypeCOP restype( has_target_type ? nullptr /*set below*/ : pose.residue_type_ptr(ir) );
		if ( has_target_type ) {
			core::chemical::ResidueTypeFinder finder( *rsd_type_set );
			finder.residue_base_name( target_type() ).variants( pose.residue_type(ir).variant_type_enums(), pose.residue_type(ir).custom_variant_types() );
			restype = finder.get_representative_type();
			if ( restype == nullptr ) {
				TR.Warning << "No representative type found for " << target_type() << " equivalent to " << pose.residue(ir).name() << "." << std::endl;
				continue; //No sensible type.
			}
		}
		TR << "Set residue type for position " << ir << " to " << restype->name() << "." << std::endl;

		core::chemical::ResidueTypeCOP other_restype;
		if ( pose.residue(ir).has_upper_connect() ) {
			core::Size const other_res( pose.residue(ir).connected_residue_at_resconn( pose.residue_type(ir).upper_connect_id() ) );
			if ( other_res > 0 && other_res <= nres ) {
				other_restype = pose.residue_type_ptr( other_res );
			}
		}
		debug_assert(other_restype != nullptr); //Should be true.

		utility::vector1< core::Real > mainchain_tors( pose.residue(ir).mainchain_torsions().size() - 1 );
		for ( core::Size i(1), imax( pose.residue(ir).mainchain_torsions().size() ); i<imax; ++i ) mainchain_tors[i] = numeric::nonnegative_principal_angle_degrees( pose.residue(ir).mainchain_torsions()[i] );

		core::Real scoreval(0.0);
		rama.eval_rpp_rama_score( pose.conformation(), restype, other_restype, mainchain_tors, scoreval, gradient, false);
		core::Real const scoreval_multiplied( scoreval * rama_prepro_multiplier() );
		TR << "The rama_prepro energy for position " << ir << " is " << scoreval << " (" << scoreval_multiplied << " when multiplied by the weight coefficient), which is " << ( scoreval_multiplied <= score_threshold() ? "below" : "above" ) << " the threshold for selection." << std::endl;
		if ( scoreval_multiplied <= score_threshold() ) selection[ir] = true;

	}

	return selection;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
RamaMutationSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	set_target_type( tag->getOption<std::string>( "target_type", target_type() ) );
	set_score_threshold( tag->getOption<core::Real>( "score_threshold", score_threshold() ) );
	set_rama_prepro_multiplier( tag->getOption<core::Real>( "rama_prepro_multiplier", rama_prepro_multiplier() ) );
}

std::string RamaMutationSelector::get_name() const
{
	return RamaMutationSelector::class_name();
}

std::string RamaMutationSelector::class_name()
{
	return "RamaMutationSelector";
}

void RamaMutationSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "target_type", xs_string,
		"The target residue type (full name, not three-letter code) that we are considering.  This residue selector "
		"will select positions which, when mutated to this residue type, have a rama_prepro score below the score "
		"threshold.  If the type is set to an empty string (\"\"), the current residue type at each position is "
		"considered instead.",
		"" )
		+ XMLSchemaAttribute::attribute_w_default( "score_threshold", xsct_real,
		"The rama_prepro score threshold above which positions are not selected.",
		"0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rama_prepro_multiplier", xsct_real,
		"A weighting factor for rama_prepro scoring.  The calculated rama_prepro value is multiplied by this factor "
		"before being compared to the threshold value.  This is for convenience: rama_prepro is generally given a weight "
		"different than 1.0 in the score function.  This defaults to 0.45 to match the weight in the beta_nov15 score "
		"function.",
		"0.45" );

	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(),
		"The RamaMutationSelector selects positions which, when mutated to a desired residue type, would have a rama_prepro "
		"score below a user-defined threshold.  This is useful for selecting those positions that ought to be mutated to a "
		"particular conformationally-constrained residue type, such as proline or alpha-aminoisobutyric acid.",
		attributes );

}

/// @brief Set the residue type (full name) to which we are considering
/// mutations.  If set to an empty string, selection is based on rama_prepro scoring
/// of the residue type at each position.
void
RamaMutationSelector::set_target_type(
	std::string const &type_in
) {
	target_type_ = type_in;
	TR << "Set target type to " << target_type_ << "." << std::endl;
}

/// @brief Set the score threshold.
/// @details Positions which, when mutated to the target type, have a rama_prepro
/// score above this threshold are not selected.
void
RamaMutationSelector::set_score_threshold(
	core::Real const &threshold_in
) {
	score_threshold_ = threshold_in;
	TR << "Set rama_prepro score threshold to " << score_threshold_ << "." << std::endl;
}

/// @brief Set the weight multiplier of the rama_prepro term.
/// @details Defaults to 0.45 to match beta_nov15.
void
RamaMutationSelector::set_rama_prepro_multiplier(
	core::Real const &multiplier_in
) {
	runtime_assert_string_msg( multiplier_in > 0.0, "Error in protocols::cyclic_peptide::RamaMutationSelector::set_rama_prepro_multiplier(): The weight must be greater than zero!" );
	rama_prepro_multiplier_ = multiplier_in;
	TR << "Set the rama_prepro multiplier to " << rama_prepro_multiplier_ << "." << std::endl;
}

core::select::residue_selector::ResidueSelectorOP
RamaMutationSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new RamaMutationSelector );
}

/// @brief Is a position in a pose a terminal position?
/// @details Returns true if (a) it has a terminal type, (b) its lower_connect is missing,
/// (c) its lower_connect is unconnected, (d) its upper_connect is missing, or (e) its
/// upper_connect is unconnected.
bool
RamaMutationSelector::is_terminus(
	core::conformation::Residue const &rsd
) const {

	return (
		rsd.is_terminus() ||
		!rsd.has_lower_connect() ||
		!rsd.has_upper_connect() ||
		rsd.connected_residue_at_lower() == 0 ||
		rsd.connected_residue_at_upper() == 0
	);
}

std::string
RamaMutationSelectorCreator::keyname() const {
	return RamaMutationSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
RamaMutationSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	RamaMutationSelector::provide_xml_schema( xsd );
}

} //protocols
} //cyclic_peptide

#ifdef    SERIALIZATION
template< class Archive >
void
protocols::cyclic_peptide::RamaMutationSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( target_type_ ) ); // std::string
	arc( CEREAL_NVP( score_threshold_ ) ); // core::Real
	arc( CEREAL_NVP( rama_prepro_multiplier_ ) ); // core::Real
}

template< class Archive >
void
protocols::cyclic_peptide::RamaMutationSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( target_type_ ); // std::string
	arc( score_threshold_ ); // core::Real
	arc( rama_prepro_multiplier_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::cyclic_peptide::RamaMutationSelector );
CEREAL_REGISTER_TYPE( protocols::cyclic_peptide::RamaMutationSelector )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_cyclic_peptide_RamaMutationSelector )
#endif // SERIALIZATION
