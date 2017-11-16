// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PatchdockTransform.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PatchdockTransform.hh>
#include <protocols/protein_interface_design/movers/PatchdockTransformCreator.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/protein_interface_design/read_patchdock.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.PatchdockTransform" );

// XRW TEMP std::string
// XRW TEMP PatchdockTransformCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return PatchdockTransform::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP PatchdockTransformCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new PatchdockTransform );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP PatchdockTransform::mover_name()
// XRW TEMP {
// XRW TEMP  return "PatchdockTransform";
// XRW TEMP }

PatchdockTransform::PatchdockTransform() :
	Mover()
{
	pd_reader_ = PatchdockReaderOP( new protocols::protein_interface_design::PatchdockReader ); // initialize the patchdock reader object
}

PatchdockTransform::~PatchdockTransform() {}

void
PatchdockTransform::apply( core::pose::Pose & pose )
{
	if ( pd_reader()->random_entry() ) { //randomize the entry number
		core::Size const actual_last_entry( std::min( pd_reader()->to_entry(), pd_reader()->number_of_patchdock_entries() ) );
		TR<<"sampling a number between "<<pd_reader()->from_entry()<<" and "<<actual_last_entry<<std::endl;

		pd_reader()->patchdock_entry_num( ( core::Size ) floor( numeric::random::rg().uniform() * ( actual_last_entry - pd_reader()->from_entry() + 1 ) ) + pd_reader()->from_entry() );
		TR<<"Patchdock entry: "<<pd_reader()->patchdock_entry_num()<<std::endl;
	}
	protocols::protein_interface_design::Transformation t( pd_reader()->read_patchdock_entry() );
	pd_reader()->transform_pose( pose, 2/*chain*/, t );
}

// XRW TEMP std::string
// XRW TEMP PatchdockTransform::get_name() const {
// XRW TEMP  return PatchdockTransform::mover_name();
// XRW TEMP }

void
PatchdockTransform::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &,
	protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	pd_reader()->patchdock_fname( tag->getOption< std::string >( "fname", pd_reader()->patchdock_fname() ) );
	pd_reader()->from_entry( tag->getOption< core::Size >( "from_entry", pd_reader()->from_entry() ) );
	pd_reader()->to_entry( tag->getOption< core::Size >( "to_entry", pd_reader()->to_entry() ) );
	pd_reader()->random_entry( tag->getOption< bool > ( "random_entry", pd_reader()->random_entry() ) );

	TR << "PatchdockTransform parsed with parameters: fname " << pd_reader()->patchdock_fname() << " random_entry: " <<
		pd_reader()->random_entry() << " from_entry: " << pd_reader()->from_entry() << " to entry: " <<
		pd_reader()->to_entry() << std::endl;
	runtime_assert( pd_reader()->from_entry() <= pd_reader()->to_entry() );
	runtime_assert( pd_reader()->from_entry() >= 1 );
	runtime_assert( pd_reader()->patchdock_fname() != "" );
}

protocols::moves::MoverOP
PatchdockTransform::clone() const {
	return( protocols::moves::MoverOP( new PatchdockTransform( *this ) ));
}

PatchdockReaderOP
PatchdockTransform::pd_reader() const{
	return pd_reader_;
}

void
PatchdockTransform::pd_reader( PatchdockReaderOP p ){
	pd_reader_ = p;
}

std::string PatchdockTransform::get_name() const {
	return mover_name();
}

std::string PatchdockTransform::mover_name() {
	return "PatchdockTransform";
}

void PatchdockTransform::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute( "fname", xs_string, "Patchdock filename to be read" )
		+ XMLSchemaAttribute( "from_entry", xs_string, "Entry from which to read" )
		+ XMLSchemaAttribute( "to_entry", xs_string, "Entry to which to write" )
		+ XMLSchemaAttribute( "random_entry", xsct_rosetta_bool, "Use random from and to entries" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string PatchdockTransformCreator::keyname() const {
	return PatchdockTransform::mover_name();
}

protocols::moves::MoverOP
PatchdockTransformCreator::create_mover() const {
	return protocols::moves::MoverOP( new PatchdockTransform );
}

void PatchdockTransformCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PatchdockTransform::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
