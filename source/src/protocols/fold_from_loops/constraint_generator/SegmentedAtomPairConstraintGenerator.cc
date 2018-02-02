// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/SegmentedAtomPairConstraintGenerator.cc
/// @brief Given a set of non-continuous selected segments, generates differently scored atom pair constraints
///       for the resides in each segment and between segments.
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Unit headers
#include <protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGenerator.hh>
#include <protocols/fold_from_loops/constraint_generator/SegmentedAtomPairConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.fold_from_loops.SegmentedAtomPairConstraintGenerator" );

namespace protocols {
namespace fold_from_loops {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
SegmentedAtomPairConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new SegmentedAtomPairConstraintGenerator );
}

std::string
SegmentedAtomPairConstraintGeneratorCreator::keyname() const
{
	return SegmentedAtomPairConstraintGenerator::class_name();
}

SegmentedAtomPairConstraintGenerator::SegmentedAtomPairConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( SegmentedAtomPairConstraintGenerator::class_name() ),
	selector_( new core::select::residue_selector::TrueResidueSelector ),
	reference_pose_( new core::pose::Pose ),
	inner_( default_inner() ),
	outer_( default_outer() ),
	do_inner_( default_do_inner() ),
	do_outer_( default_do_outer() )
{}

SegmentedAtomPairConstraintGenerator::~SegmentedAtomPairConstraintGenerator() = default;

protocols::constraint_generator::ConstraintGeneratorOP
SegmentedAtomPairConstraintGenerator::clone() const
{
	return SegmentedAtomPairConstraintGeneratorOP( new SegmentedAtomPairConstraintGenerator( *this ) );
}

void
SegmentedAtomPairConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	bool const use_native = tag->getOption< bool >( "native", false );
	if ( use_native ) {
		set_reference_pose( protocols::constraint_generator::get_native_pose() );
		if ( ! reference_pose_->empty() ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'native' option for SegmentedAtomPairConstraintGenerator specified, but no native pose is availible.");
		}
	}

	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( *selector );

	if ( !tag->hasTag( "Inner" ) ) {
		do_inner_ = false;
	} else {
		utility::tag::TagCOP const& innertag = tag->getTag( "Inner" );
		do_inner_ = true;
		ConstraintConditions inner = default_inner();
		set_inner_sd( innertag->getOption< core::Real >( "sd", inner.sd ) );
		set_inner_weight( innertag->getOption< core::Real >( "weight", inner.weight ) );
		set_inner_ca_only( innertag->getOption< bool >( "ca_only", inner.ca_only ) );
		set_inner_min_seq_sep( innertag->getOption< core::Size >( "min_seq_sep", inner.min_seq_sep ) );
		set_inner_use_harmonic_function( innertag->getOption< bool >( "use_harmonic", inner.use_harmonic ) );
		set_inner_unweighted_function( innertag->getOption< bool >( "unweighted", inner.unweighted ) );
	}
	if ( !tag->hasTag( "Outer" ) ) {
		do_outer_ = false;
	} else {
		utility::tag::TagCOP const& outertag = tag->getTag( "Outer" );
		do_outer_ = true;
		ConstraintConditions outer = default_outer();
		set_outer_sd( outertag->getOption< core::Real >( "sd", outer.sd ) );
		set_outer_weight( outertag->getOption< core::Real >( "weight", outer.weight ) );
		set_outer_ca_only( outertag->getOption< bool >( "ca_only", outer.ca_only ) );
		set_outer_max_distance( outertag->getOption< core::Real >( "max_distance", outer.max_distance ) );
		set_outer_use_harmonic_function( outertag->getOption< bool >( "use_harmonic", outer.use_harmonic ) );
		set_outer_unweighted_function( outertag->getOption< bool >( "unweighted", outer.unweighted ) );
	}

	if ( !selector_ ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "SegmentedAtomPairConstraintGenerator requires a residue selector, but one is not set.\n" );
	}
}

core::scoring::constraints::ConstraintCOPs
SegmentedAtomPairConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	using namespace protocols::constraint_generator;
	debug_assert( selector_ );
	core::select::residue_selector::ResidueRanges ranges;
	ranges.from_subset( selector_->apply( pose ) );

	core::scoring::constraints::ConstraintCOPs constraints;
	// Inner Constraints Params
	AtomPairConstraintGenerator innerapc;
	if ( !reference_pose_->empty() ) innerapc.set_reference_pose( reference_pose_ );
	innerapc.set_id( "INNER" );
	innerapc.set_sd( inner_.sd );
	innerapc.set_weight( inner_.weight );
	innerapc.set_ca_only( inner_.ca_only );
	innerapc.set_min_seq_sep( inner_.min_seq_sep );
	innerapc.set_max_distance( inner_.max_distance );
	innerapc.set_use_harmonic_function( inner_.use_harmonic );
	innerapc.set_unweighted_function( inner_.unweighted );
	// Outer Constraints Params
	AtomPairConstraintGenerator outerapc;
	if ( !reference_pose_->empty() ) outerapc.set_reference_pose( reference_pose_ );
	outerapc.set_id( "OUTER" );
	outerapc.set_sd( outer_.sd );
	outerapc.set_weight( outer_.weight );
	outerapc.set_ca_only( outer_.ca_only );
	outerapc.set_min_seq_sep( outer_.min_seq_sep );
	outerapc.set_max_distance( outer_.max_distance );
	outerapc.set_use_harmonic_function( outer_.use_harmonic );
	outerapc.set_unweighted_function( outer_.unweighted );


	for ( core::Size i = 1; i<= ranges.size(); ++i ) {
		core::select::residue_selector::ResidueIndexSelector inner_selector( ranges[i].to_string() );
		if ( do_inner_ ) {
			innerapc.set_residue_selector( inner_selector );
			constraints.append( innerapc.apply( pose ) );
		}
		for ( core::Size j = i + 1; j<= ranges.size(); ++j ) {
			if ( do_outer_ ) {
				core::select::residue_selector::ResidueIndexSelector outer_selector( ranges[j].to_string() );
				outerapc.set_residue_selector( inner_selector );
				outerapc.set_secondary_residue_selector( outer_selector );
				constraints.append( outerapc.apply( pose ) );
			}
		}
	}
	return constraints;
}

void
SegmentedAtomPairConstraintGenerator::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}

void
SegmentedAtomPairConstraintGenerator::set_inner_sd( core::Real const sd )
{
	inner_.sd = sd;
}

void
SegmentedAtomPairConstraintGenerator::set_inner_ca_only( bool const ca_only )
{
	inner_.ca_only = ca_only;
}

void
SegmentedAtomPairConstraintGenerator::set_inner_use_harmonic_function( bool const use_harmonic )
{
	inner_.use_harmonic = use_harmonic;
}

void
SegmentedAtomPairConstraintGenerator::set_inner_unweighted_function( bool const unweighted )
{
	inner_.unweighted = unweighted;
}

void
SegmentedAtomPairConstraintGenerator::set_inner_min_seq_sep( core::Size const min_seq_sep )
{
	inner_.min_seq_sep = min_seq_sep;
}

void
SegmentedAtomPairConstraintGenerator::set_inner_weight( core::Real const weight )
{
	inner_.weight = weight;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_sd( core::Real const sd )
{
	outer_.sd = sd;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_ca_only( bool const ca_only )
{
	outer_.ca_only = ca_only;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_use_harmonic_function( bool const use_harmonic )
{
	outer_.use_harmonic = use_harmonic;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_unweighted_function( bool const unweighted )
{
	outer_.unweighted = unweighted;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_max_distance( core::Real const max_dist )
{
	outer_.max_distance = max_dist;
}

void
SegmentedAtomPairConstraintGenerator::set_outer_weight( core::Real const weight )
{
	outer_.weight = weight;
}

void
SegmentedAtomPairConstraintGenerator::set_reference_pose( core::pose::PoseCOP ref_pose )
{
	reference_pose_->detached_copy( *ref_pose );
}

void
SegmentedAtomPairConstraintGenerator::set_do_inner( bool pick )
{
	do_inner_ = pick;
}

void
SegmentedAtomPairConstraintGenerator::set_do_outer( bool pick )
{
	do_outer_ = pick;
}

void
SegmentedAtomPairConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist, innerlist, outerlist;
	attlist
		+ XMLSchemaAttribute( "name", xs_string, "The name given to this instance." )
		+ XMLSchemaAttribute::attribute_w_default( "native", xsct_rosetta_bool, "Restrain to native distance?", "false" );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "Selector specifying residues to be constrained. When not provided, all residues are selected" );

	ConstraintConditions inner = default_inner();
	innerlist
		+ XMLSchemaAttribute::attribute_w_default( "sd", xsct_real, "Standard deviation for distance constraint", std::to_string(inner.sd) )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Weight of distance constraint", std::to_string(inner.weight) )
		+ XMLSchemaAttribute::attribute_w_default( "ca_only", xsct_rosetta_bool, "Only make constraints between alpha carbons", std::to_string(inner.ca_only) )
		+ XMLSchemaAttribute::attribute_w_default( "use_harmonic", xsct_rosetta_bool, "If true, use harmonic function instead of SOG function", std::to_string(inner.use_harmonic) )
		+ XMLSchemaAttribute::attribute_w_default( "unweighted", xsct_rosetta_bool, "If true, SCALARWEIGHTEDFUNC is not added to the constraint definition", std::to_string(inner.unweighted) )
		+ XMLSchemaAttribute::attribute_w_default( "min_seq_sep", xsct_non_negative_integer, "Minimum sequence separation between constrained residues", std::to_string(inner.min_seq_sep) );

	ConstraintConditions  outer = default_outer();
	outerlist
		+ XMLSchemaAttribute::attribute_w_default( "sd", xsct_real, "Standard deviation for distance constraint", std::to_string(outer.sd) )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Weight of distance constraint", std::to_string(outer.weight) )
		+ XMLSchemaAttribute::attribute_w_default( "ca_only", xsct_rosetta_bool, "Only make constraints between alpha carbons", std::to_string(outer.ca_only) )
		+ XMLSchemaAttribute::attribute_w_default( "use_harmonic", xsct_rosetta_bool, "If true, use harmonic function instead of SOG function", std::to_string(outer.use_harmonic) )
		+ XMLSchemaAttribute::attribute_w_default( "unweighted", xsct_rosetta_bool, "If true, SCALARWEIGHTEDFUNC is not added to the constraint definition", std::to_string(outer.unweighted)  )
		+ XMLSchemaAttribute::attribute_w_default( "max_distance", xsct_real, "Do not add constraints if atoms are farther apart than this", std::to_string(outer.max_distance) );

	XMLSchemaRepeatableCTNodeOP innernode( new XMLSchemaRepeatableCTNode );
	innernode->set_element_w_attributes( "Inner", innerlist, "Describes the pair constraint properties between residues in the same selection contiguous segment" );
	XMLSchemaRepeatableCTNodeOP outernode( new XMLSchemaRepeatableCTNode );
	outernode->set_element_w_attributes( "Outer", outerlist, "Describes the pair constraint properties between residues in different selection segments" );

	XMLSchemaRepeatableCTNodeOP root_node(new XMLSchemaRepeatableCTNode);
	root_node->set_element_w_attributes( class_name(), attlist, "Generates different atom pair constraints between residues in the same contiguous segment and between segments." );
	root_node->set_root_node_naming_func( & protocols::constraint_generator::ConstraintGeneratorFactory::complex_type_name_for_constraint_generator );
	root_node->set_kids_naming_func( & add_sub_ct_name );
	root_node->add_child( innernode );
	root_node->add_child( outernode );
	root_node->recursively_write_ct_to_schema( xsd );
}

void
SegmentedAtomPairConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SegmentedAtomPairConstraintGenerator::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
