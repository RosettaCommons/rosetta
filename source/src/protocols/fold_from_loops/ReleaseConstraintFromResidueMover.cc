// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/RelaseConstraintFromResidueMoverCreator.hh
/// @brief Remove constraints from the residues specified by a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/ReleaseConstraintFromResidueMover.hh>
#include <protocols/fold_from_loops/ReleaseConstraintFromResidueMoverCreator.hh>

// Protocol headers


// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.fold_from_loops.ReleaseConstraintFromResidueMover" );

namespace protocols {
namespace fold_from_loops {

ReleaseConstraintFromResidueMover::ReleaseConstraintFromResidueMover():
	protocols::moves::Mover( ReleaseConstraintFromResidueMover::mover_name() )
{
	core::select::residue_selector::ResidueSelectorCOP true_selector( new core::select::residue_selector::TrueResidueSelector );
	core::select::residue_selector::NotResidueSelector false_selector( true_selector );
	set_residue_selector( false_selector );
}

ReleaseConstraintFromResidueMover::ReleaseConstraintFromResidueMover( core::select::residue_selector::ResidueSelectorCOP const & selector ):
	protocols::moves::Mover( ReleaseConstraintFromResidueMover::mover_name() ),
	selector_( selector )
{
}

ReleaseConstraintFromResidueMover::~ReleaseConstraintFromResidueMover()= default;

void
ReleaseConstraintFromResidueMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
	if ( selector ) set_residue_selector( *selector );
}

protocols::moves::MoverOP
ReleaseConstraintFromResidueMover::clone() const
{
	return protocols::moves::MoverOP( new ReleaseConstraintFromResidueMover( *this ) );
}

protocols::moves::MoverOP
ReleaseConstraintFromResidueMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ReleaseConstraintFromResidueMover );
}

void
ReleaseConstraintFromResidueMover::apply( core::pose::Pose & pose )
{
	core::select::residue_selector::ResidueSubset const subset = selector_->apply( pose );
	if ( core::select::residue_selector::all_false_selection( subset ) ) return;
	core::scoring::constraints::ConstraintCOPs csts = pose.constraint_set()->get_all_constraints();
	core::scoring::constraints::ConstraintCOPs todel;
	for ( auto const & constraint : csts ) {
		utility::vector1< core::Size > residues = constraint->residues();
		for ( auto const & res : residues ) {
			if ( subset[res] ) {
				todel.push_back( constraint );
				break;
			}
		}
	}
	pose.remove_constraints( todel );
}

void
ReleaseConstraintFromResidueMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP const & selector )
{
	selector_ = selector;
}

void
ReleaseConstraintFromResidueMover::set_residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}

std::string ReleaseConstraintFromResidueMover::get_name() const {
	return mover_name();
}

std::string ReleaseConstraintFromResidueMover::mover_name() {
	return "ReleaseConstraintFromResidue";
}

void ReleaseConstraintFromResidueMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "Selector specifying residues to be released." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Release constraints from the specified residues.", attlist );
}

std::string ReleaseConstraintFromResidueMoverCreator::keyname() const {
	return ReleaseConstraintFromResidueMover::mover_name();
}

protocols::moves::MoverOP
ReleaseConstraintFromResidueMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReleaseConstraintFromResidueMover );
}

void ReleaseConstraintFromResidueMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReleaseConstraintFromResidueMover::provide_xml_schema( xsd );
}


} //protocols
} //fold_from_loops
