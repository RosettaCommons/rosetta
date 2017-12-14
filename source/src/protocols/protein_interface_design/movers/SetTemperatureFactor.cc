// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/SetTemperatureFactor.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/SetTemperatureFactor.hh>
#include <protocols/protein_interface_design/movers/SetTemperatureFactorCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/ResId.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <protocols/filters/Filter.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;
using protocols::moves::ResId;

static basic::Tracer TR( "protocols.protein_interface_design.movers.SetTemperatureFactor" );

// XRW TEMP std::string
// XRW TEMP SetTemperatureFactorCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SetTemperatureFactor::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetTemperatureFactorCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetTemperatureFactor );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetTemperatureFactor::mover_name()
// XRW TEMP {
// XRW TEMP  return "SetTemperatureFactor";
// XRW TEMP }

SetTemperatureFactor::SetTemperatureFactor() :
	Mover( SetTemperatureFactor::mover_name() ),
	filter_( /* NULL */ ),
	scaling_( 1.0 )
{
}


SetTemperatureFactor::~SetTemperatureFactor() = default;

void
SetTemperatureFactor::apply( core::pose::Pose & pose )
{
	auto * resid = dynamic_cast< ResId *>( filter_.get() );
	runtime_assert( resid );

	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		resid->set_resid( resi );
		core::Real const filter_value( filter_->report_sm( pose ) );
		for ( core::Size j = 1; j <= pose.residue( resi ).natoms(); ++j ) {
			pose.pdb_info()->temperature( resi, j, filter_value * scaling() );
		}
	}
}

// XRW TEMP std::string
// XRW TEMP SetTemperatureFactor::get_name() const {
// XRW TEMP  return SetTemperatureFactor::mover_name();
// XRW TEMP }

void
SetTemperatureFactor::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const & filters, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	std::string const filter_name( tag->getOption<std::string>( "filter" ) );
	auto find_ap_filter( filters.find( filter_name));
	bool const ap_filter_found( find_ap_filter != filters.end() );
	runtime_assert( ap_filter_found );
	filter( find_ap_filter->second->clone() );
	scaling( tag->getOption< core::Real >( "scaling", 1.0 ) );
}

protocols::moves::MoverOP
SetTemperatureFactor::clone() const {
	return( protocols::moves::MoverOP( new SetTemperatureFactor( *this ) ));
}

void
SetTemperatureFactor::filter( protocols::filters::FilterOP filter )
{
	auto * resid = dynamic_cast< ResId *>( filter.get() );
	if ( !resid ) {
		utility_exit_with_message( "filter must be derived from ResId class for this to work." );
	}
	filter_ = filter;
}

void
SetTemperatureFactor::scaling( core::Real const scaling ){
	scaling_ = scaling;
}

core::Real
SetTemperatureFactor::scaling() const{
	return scaling_;
}

std::string SetTemperatureFactor::get_name() const {
	return mover_name();
}

std::string SetTemperatureFactor::mover_name() {
	return "SetTemperatureFactor";
}

void SetTemperatureFactor::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "filter", xs_string, "Filter to apply" )
		+ XMLSchemaAttribute::attribute_w_default( "scaling", xsct_real, "Scale filter value by this amount prior to storing in temperature factor", "1.0" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Co-opt the PDB field for temperature factors to instead store filter values", attlist );
}

std::string SetTemperatureFactorCreator::keyname() const {
	return SetTemperatureFactor::mover_name();
}

protocols::moves::MoverOP
SetTemperatureFactorCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetTemperatureFactor );
}

void SetTemperatureFactorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetTemperatureFactor::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
