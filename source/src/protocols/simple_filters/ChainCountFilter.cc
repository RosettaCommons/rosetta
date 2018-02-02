// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/ChainCountFilter.cc
/// @brief Count the amount of chains in a pose
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#include <protocols/simple_filters/ChainCountFilter.hh>
#include <protocols/simple_filters/ChainCountFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility>
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.filters.ChainCountFilter" );

/// @brief default ctor
ChainCountFilter::ChainCountFilter() :
	parent( "ChainCount" ),
	threshold_( 1 )
{}

/// @return Whether the atom pair is within the cutoff
bool ChainCountFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	report( TR.Debug, pose );
	if ( dist <= core::Real(threshold_) ) return true;
	return false;
}

core::Real
ChainCountFilter::compute( core::pose::Pose const & pose ) const
{
	return pose.num_chains();
}

core::Real
ChainCountFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	return( dist );
}

void ChainCountFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	out << "Total number of chains: " << dist << std::endl;
}

void ChainCountFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	threshold( tag->getOption< core::Size >( "max_chains", 1 ) );
}

std::string ChainCountFilter::name() const {
	return class_name();
}

std::string ChainCountFilter::class_name() {
	return "ChainCountFilter";
}

void ChainCountFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "max_chains", xs_integer, "Max allowed number of chains (filter fails if more) ", "1" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Counts the number of chains in the Pose.", attlist );
}

std::string ChainCountFilterCreator::keyname() const {
	return ChainCountFilter::class_name();
}

protocols::filters::FilterOP
ChainCountFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ChainCountFilter );
}

void ChainCountFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ChainCountFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
