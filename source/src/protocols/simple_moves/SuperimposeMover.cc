// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   SuperimposeMover.cc
///
/// @brief
/// @author Ingemar Andre

// unit headers
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/simple_moves/SuperimposeMoverCreator.hh>

// type headers
#include <core/types.hh>
#include <core/id/types.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>

// utility header
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <utility/tag/Tag.hh>

//option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SuperimposeMover" );

// XRW TEMP std::string
// XRW TEMP SuperimposeMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SuperimposeMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SuperimposeMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SuperimposeMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SuperimposeMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "Superimpose";
// XRW TEMP }

SuperimposeMover::SuperimposeMover() :
	protocols::moves::Mover("SuperimposeMover"),
	ref_pose_(/* 0 */)
{}

SuperimposeMover::SuperimposeMover( Pose const & ref_pose ) :
	protocols::moves::Mover("SuperimposeMover"),
	ref_pose_(core::pose::PoseOP( new Pose(ref_pose) ))
{}

SuperimposeMover::SuperimposeMover(Pose const & ref_pose, core::Size ref_start, core::Size ref_end, core::Size target_start, core::Size target_end, bool CA_only):
	protocols::moves::Mover("SuperimposeMover"),
	ref_pose_(core::pose::PoseOP( new Pose(ref_pose) )),
	ref_start_(ref_start),
	ref_end_(ref_end),
	target_start_(target_start),
	target_end_(target_end),
	CA_only_(CA_only)
{
}

SuperimposeMover::~SuperimposeMover() = default;

protocols::moves::MoverOP
SuperimposeMover::clone() const
{
	return protocols::moves::MoverOP( new SuperimposeMover( *this ) );
}

protocols::moves::MoverOP
SuperimposeMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SuperimposeMover() );
}

void
SuperimposeMover::set_reference_pose( Pose const & pose,Size start, Size end ) {
	ref_pose_ = core::pose::PoseOP( new Pose(pose) );
	ref_start_ = start;
	ref_end_ = (end == 0) ? pose.size() : end;
	runtime_assert(ref_start_ > 0 && ref_start_ < ref_end_ && ref_end_ <= pose.size());
}

void
SuperimposeMover::set_target_range( Size start, Size end ) {
	target_start_ = start;
	target_end_ = end;
	runtime_assert(target_start_ > 0 && target_start_ < target_end_);
}

void
SuperimposeMover::set_ca_only(bool setting){
	CA_only_ = setting;
}

/// @details copied and modified from calpha_superimpose_pose
core::Real
SuperimposeMover::superimpose(
	core::pose::Pose & mod_pose,
	core::pose::Pose const & ref_pose,
	Size ref_start,
	Size ref_end,
	Size target_start,
	Size /*target_end*/
)
{
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	std::map< core::id::AtomID, core::id::AtomID> atom_id_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::AtomID::BOGUS_ATOM_ID() );
	for ( Size i_target = target_start, i_ref = ref_start; i_ref <= ref_end; ++i_ref, ++i_target ) {
		if ( ! mod_pose.residue(i_target).has("CA") ) continue;
		if ( ! ref_pose.residue(i_ref).has("CA") ) continue;

		core::id::AtomID const id1( mod_pose.residue(i_target).atom_index("CA"), i_target );
		core::id::AtomID const id2( ref_pose.residue(i_ref).atom_index("CA"), i_ref );
		atom_map.set( id1, id2 );
		atom_id_map.insert( std::make_pair(id1, id2) );

	}
	return core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}

/// @details copied and modified from calpha_superimpose_pose
core::Real
SuperimposeMover::superimposebb(
	core::pose::Pose & mod_pose,
	core::pose::Pose const & ref_pose,
	Size ref_start,
	Size ref_end,
	Size target_start,
	Size /*target_end*/
)
{
	core::id::AtomID_Map< core::id::AtomID > atom_map;
	std::map< core::id::AtomID, core::id::AtomID> atom_id_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, core::id::AtomID::BOGUS_ATOM_ID() );
	for ( Size i_target = target_start, i_ref = ref_start; i_ref <= ref_end; ++i_ref, ++i_target ) {

		if ( ! mod_pose.residue(i_target).has("N") ) continue;
		if ( ! ref_pose.residue(i_ref).has("N") ) continue;
		core::id::AtomID const id1( mod_pose.residue(i_target).atom_index("N"), i_target );
		core::id::AtomID const id2( ref_pose.residue(i_ref).atom_index("N"), i_ref );
		atom_map.set( id1, id2 );
		atom_id_map.insert( std::make_pair(id1, id2) );

		if ( ! mod_pose.residue(i_target).has("CA") ) continue;
		if ( ! ref_pose.residue(i_ref).has("CA") ) continue;
		core::id::AtomID const id3( mod_pose.residue(i_target).atom_index("CA"), i_target );
		core::id::AtomID const id4( ref_pose.residue(i_ref).atom_index("CA"), i_ref );
		atom_map.set( id3, id4 );
		atom_id_map.insert( std::make_pair(id3, id4) );

		if ( ! mod_pose.residue(i_target).has("C") ) continue;
		if ( ! ref_pose.residue(i_ref).has("C") ) continue;
		core::id::AtomID const id5( mod_pose.residue(i_target).atom_index("C"), i_target );
		core::id::AtomID const id6( ref_pose.residue(i_ref).atom_index("C"), i_ref );
		atom_map.set( id5, id6 );
		atom_id_map.insert( std::make_pair(id5, id6) );

		if ( ! mod_pose.residue(i_target).has("O") ) continue;
		if ( ! ref_pose.residue(i_ref).has("O") ) continue;
		core::id::AtomID const id7( mod_pose.residue(i_target).atom_index("O"), i_target );
		core::id::AtomID const id8( ref_pose.residue(i_ref).atom_index("O"), i_ref );
		atom_map.set( id7, id8 );
		atom_id_map.insert( std::make_pair(id7, id8) );

	}
	return core::scoring::superimpose_pose( mod_pose, ref_pose, atom_map );
}


void
SuperimposeMover::apply( Pose & pose ) {
	using namespace basic::options;

	if ( ref_pose_ == nullptr ) {
		TR << "using -in:file:native as the reference pose " <<  std::endl;
		ref_pose_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ].value() , core::import_pose::PDB_file);
	}

	const Size ref_start = ref_start_;
	const Size target_start = target_start_;
	const Size ref_end = (ref_end_ == 0) ? ref_pose_->size() : ref_end_;
	const Size target_end = (target_end_ == 0) ? pose.size() : target_end_;

	TR << "ref_start: "<< ref_start << " ref_end " << ref_end <<std::endl;
	TR << "target_start: "<< target_start << " target_end " << target_end <<std::endl;
	runtime_assert(ref_start > 0 && ref_start < ref_end && ref_end <= pose.size());
	runtime_assert_msg(ref_end - ref_start == target_end - target_start, "segments to superimpose have different lengths!");

	if ( CA_only_ ) {
		core::Real rms  = superimpose( pose, *ref_pose_, ref_start, ref_end, target_start, target_end );
		TR << "CA RMS to reference: " << rms << std::endl;
	} else {
		core::Real rms  = superimposebb( pose, *ref_pose_, ref_start, ref_end, target_start, target_end );
		TR << "Backbone RMS to reference: " << rms << std::endl;
	}
}

// XRW TEMP std::string
// XRW TEMP SuperimposeMover::get_name() const {
// XRW TEMP  return "SuperimposeMover";
// XRW TEMP }

void
SuperimposeMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	ref_start_ = tag->getOption< Size >("ref_start",1);
	ref_end_ = tag->getOption< Size >("ref_end",0);
	target_start_ = tag->getOption< Size >("target_start",1);
	target_end_ = tag->getOption< Size >("target_end",0);
	CA_only_ = tag->getOption< bool >("CA_only",1);
	if ( tag->hasOption("ref_pose") ) ref_pose_ = core::import_pose::pose_from_file(tag->getOption< std::string >("ref_pose"), core::import_pose::PDB_file);
	else if ( tag->hasOption("spm_reference_name") ) {
		ref_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data, "spm_reference_name");
	}
}

std::string SuperimposeMover::get_name() const {
	return mover_name();
}

std::string SuperimposeMover::mover_name() {
	return "Superimpose";
}

void SuperimposeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// Note: target_start and target_end are formally not required, but you
	// may get a no-op as a result if they retain their meaningless default values.
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "ref_start", xsct_non_negative_integer, "Starting residue for superposition, in reference", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "ref_end", xsct_non_negative_integer, "Ending residue for superposition, in reference", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "target_start", xsct_non_negative_integer, "Starting residue for superposition", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "target_end", xsct_non_negative_integer, "Ending residue for superposition", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "CA_only", xsct_rosetta_bool, "Superimpose by CA coordinates only", "1" )
		+ XMLSchemaAttribute( "ref_pose", xs_string, "Reference pose file name" );
	//+ XMLSchemaAttribute( "spm_reference_name", xs_string, "Reference pose name " );
	rosetta_scripts::attributes_for_saved_reference_pose( attlist, "spm_reference_name" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string SuperimposeMoverCreator::keyname() const {
	return SuperimposeMover::mover_name();
}

protocols::moves::MoverOP
SuperimposeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SuperimposeMover );
}

void SuperimposeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SuperimposeMover::provide_xml_schema( xsd );
}


} // moves
} // protocols
