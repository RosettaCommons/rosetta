// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/StorePoseSnapshot.cc
/// @brief Stores current residue indices in the pose as a ReferencePose.  As residues are added or subtracted,
/// this permits the user to set up movers based on the current residue indices rather than the modified indices
/// of a future state.
/// @author Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory.

// Unit headers
#include <protocols/simple_moves/StorePoseSnapshot.hh>
#include <protocols/simple_moves/StorePoseSnapshotCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


// Utility Headers

// Unit Headers

// C++ headers


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace core::chemical;
using namespace std;
using namespace numeric::conversions;

using core::pose::Pose;
using core::conformation::Residue;

static basic::Tracer TR( "protocols.simple_moves.StorePoseSnapshot" );

// XRW TEMP std::string
// XRW TEMP StorePoseSnapshotCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return StorePoseSnapshot::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StorePoseSnapshotCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new StorePoseSnapshot );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP StorePoseSnapshot::mover_name()
// XRW TEMP {
// XRW TEMP  return "StorePoseSnapshot";
// XRW TEMP }

StorePoseSnapshot::~StorePoseSnapshot() = default;

/// @brief Default constructor
///
StorePoseSnapshot::StorePoseSnapshot() :
	parent(),
	reference_pose_name_( "" ),
	overwrite_current_refpose_( false )
{}

/// @brief Copy constructor
///
StorePoseSnapshot::StorePoseSnapshot( StorePoseSnapshot const &) = default;

/// @brief Apply function -- actually apply this mover to the pose, modifying the pose.
///
void StorePoseSnapshot::apply( Pose & pose ) {
	runtime_assert_string_msg(
		reference_pose_name()!="",
		"Error in protocols::simple_moves::StorePoseSnapshot::apply(): The reference pose name is currently blank.  This must be set before the apply() function is called."
	);
	pose.reference_pose_from_current( reference_pose_name(), overwrite_current_refpose_); //Yes, that's all this mover does.  It's that simple.
	if ( TR.visible() ) {
		TR << "Stored pose snapshot " << reference_pose_name() << "." << std::endl;
		TR.flush();
	}
	return;
}

/// @brief Get the mover name.
///
// XRW TEMP std::string
// XRW TEMP StorePoseSnapshot::get_name() const {
// XRW TEMP  return StorePoseSnapshot::mover_name();
// XRW TEMP }

/// @brief Parse RosettaScripts XML to set up this mover.
/// @details This is called at script initialization, long before the apply()
/// function is called.
void StorePoseSnapshot::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & //pose
)
{
	runtime_assert_string_msg( tag->hasOption("reference_pose_name"), "Error in protocols::simple_moves::StorePoseSnapshot::parse_my_tag():  When parsing options for the StorePoseSnapshot mover, no \"reference_pose_name\" option was found.  This is required." );

	set_reference_pose_name( tag->getOption< std::string >( "reference_pose_name", "" ) );
	if ( TR.visible() ) TR << "Set reference pose name to " << reference_pose_name() << "." << std::endl;

	overwrite_current_refpose_ = tag->getOption< bool >("override_current", overwrite_current_refpose_);

	if ( TR.visible() ) TR.flush();
	return;
}

/// @brief Set the name of the reference pose object that will be created and stored in the pose.
///
void StorePoseSnapshot::set_reference_pose_name( std::string const &name_in ) {
	runtime_assert_string_msg( name_in!="", "Error in protocols::simple_moves::StorePoseSnapshot::set_reference_pose_name(): The name cannot be an empty string." );
	reference_pose_name_ = name_in;
	return;
}

/// @brief Return the name of the reference pose object that will be created and stored in the pose.
///
std::string StorePoseSnapshot::reference_pose_name( ) const { return reference_pose_name_; }

std::string StorePoseSnapshot::get_name() const {
	return mover_name();
}

std::string StorePoseSnapshot::mover_name() {
	return "StorePoseSnapshot";
}

void StorePoseSnapshot::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute(
		"reference_pose_name", xs_string,
		"The name of the snapshot or reference pose object. Many different "
		"snapshots may be stored (by applying different StorePoseSnapshot movers "
		"at different points in the protocol), and may subsequently be "
		"referred to by different names")
		+ XMLSchemaAttribute(
		"override_current", xsct_rosetta_bool,
		"If there already exists a reference pose, overwrite it." );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Rosetta's sequential residue numbering can create headaches in protocols "
		"in which residues are to be inserted or deleted, and in which one wishes "
		"subsequently to refer to residues past the insertion or deletion point. "
		"The StorePoseSnapshot mover is intended as a means of permitting a user "
		"to store a named snapshot or \"reference pose\" at a particular point in a "
		"protocol, then use the residue numbering of the pose at that point in the "
		"protocol for movers and filters that might be applied sometime later, "
		"after pose indices may have been altered by insertions or deletions.",
		attlist );
}

std::string StorePoseSnapshotCreator::keyname() const {
	return StorePoseSnapshot::mover_name();
}

protocols::moves::MoverOP
StorePoseSnapshotCreator::create_mover() const {
	return protocols::moves::MoverOP( new StorePoseSnapshot );
}

void StorePoseSnapshotCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StorePoseSnapshot::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
