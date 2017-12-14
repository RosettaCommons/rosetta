// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fldsgn/MatchResiduesMover.cc
/// @brief
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit Headers
#include <protocols/fldsgn/MatchResiduesMover.hh>
#include <protocols/fldsgn/MatchResiduesMoverCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/rms_util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <utility/exit.hh>

// Parser headers
#include <protocols/moves/Mover.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>


// Boost Headers
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/lexical_cast.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fldsgn.MatchResiduesMover" );

namespace protocols {
namespace fldsgn {

// -------------  Mover Creator -------------

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MatchResiduesMoverCreator::create_mover() const { return protocols::moves::MoverOP( new MatchResiduesMover ); }

// XRW TEMP std::string
// XRW TEMP MatchResiduesMoverCreator::keyname() const { return "MatchResiduesMover"; }

// -------------  Mover Creator -------------

// @brief default constructor
MatchResiduesMover::MatchResiduesMover() = default;


// @brief destructor
MatchResiduesMover::~MatchResiduesMover() = default;

protocols::moves::MoverOP
MatchResiduesMover::fresh_instance() const{
	return protocols::moves::MoverOP( new MatchResiduesMover() );
}

protocols::moves::MoverOP
MatchResiduesMover::clone() const{
	return protocols::moves::MoverOP( new MatchResiduesMover( *this ) );
}


void
MatchResiduesMover::apply( core::pose::Pose & pose )
{
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;

	VecSize matched_pos;
	core::Real rms = compute(pose, matched_pos);
	bool pass( rms < threshold() );
	if ( pass ) {
		set_last_move_status( MS_SUCCESS );
	} else {
		set_last_move_status( FAIL_DO_NOT_RETRY );
	}

	if ( superimpose_ && pass ) {
		superimpose_comb(pose, matched_pos);
	}
}

void
MatchResiduesMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	MatchResidues::parse_my_tag(tag, data, filters, movers, pose);
	superimpose_ = tag->getOption< bool >("superimpose", false);
}

std::string MatchResiduesMover::get_name() const {
	return mover_name();
}

std::string MatchResiduesMover::mover_name() {
	return "MatchResiduesMover";
}

void MatchResiduesMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	XMLSchemaSimpleSubelementList ssl;
	MatchResidues::provide_attributes_and_subelements( attlist, ssl );
	attlist + XMLSchemaAttribute::attribute_w_default( "superimpose", xsct_rosetta_bool, "Superimpose the input pose by its matched position", "false" );
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string MatchResiduesMoverCreator::keyname() const {
	return MatchResiduesMover::mover_name();
}

protocols::moves::MoverOP
MatchResiduesMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MatchResiduesMover );
}

void MatchResiduesMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MatchResiduesMover::provide_xml_schema( xsd );
}


} // fldsgn
} // protocols
