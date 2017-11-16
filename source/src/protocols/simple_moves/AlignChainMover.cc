// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AlignChainMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/AlignChainMover.hh>
#include <protocols/simple_moves/AlignChainMoverCreator.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.simple_moves.AlignChainMover" );

#include <utility/tag/Tag.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/xyzVector.hh>
#include <protocols/toolbox/superimpose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace simple_moves {

// XRW TEMP std::string
// XRW TEMP AlignChainMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AlignChainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AlignChainMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AlignChainMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AlignChainMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AlignChain";
// XRW TEMP }

AlignChainMover::AlignChainMover()
: moves::Mover("AlignChain"),
	pose_( /* NULL */ ),
	source_chain_( 0 ),
	target_chain_( 0 ),
	align_to_com_( 0 )
{
}

core::pose::PoseOP
AlignChainMover::pose() const{ return pose_; }

void
AlignChainMover::pose( core::pose::PoseOP pose ){ pose_ = pose; }

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coord( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size const pos : positions ) {
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

void
AlignChainMover::apply( Pose & in_pose )
{
	//default behavior
	if ( ! align_to_com_ ) {
		utility::vector1< core::Size > in_pose_positions, target_positions;
		in_pose_positions.clear(); target_positions.clear();
		runtime_assert( pose() != nullptr );

		core::Size const in_pose_chain_begin( source_chain() == 0 ? 1 : in_pose.conformation().chain_begin( source_chain() ) );
		core::Size const in_pose_chain_end  ( source_chain() == 0 ? in_pose.size() : in_pose.conformation().chain_end( source_chain() ) );
		core::Size const target_pose_chain_begin( target_chain() == 0 ? 1 : pose()->conformation().chain_begin( target_chain() ) );
		core::Size const target_pose_chain_end  ( target_chain() == 0 ? pose()->size() : pose()->conformation().chain_end( target_chain() ) );
		TR<<"In_pose from residue: "<<in_pose_chain_begin<<" to_residue: "<<in_pose_chain_end<<"\ntarget_pose from residue: "<<target_pose_chain_begin<<" to_residue: "<<target_pose_chain_end<<std::endl;
		for ( core::Size i = in_pose_chain_begin; i<=in_pose_chain_end; ++i ) {
			in_pose_positions.push_back( i );
		}
		for ( core::Size i = target_pose_chain_begin; i<=target_pose_chain_end; ++i ) {
			target_positions.push_back( i );
		}

		utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coord( in_pose, in_pose_positions ) ), ref_coords( Ca_coord( *pose(), target_positions ) );

		numeric::xyzMatrix< core::Real > rotation;
		numeric::xyzVector< core::Real > to_init_center, to_fit_center;

		using namespace protocols::toolbox;

		superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

		apply_superposition_transform( in_pose, rotation, to_init_center, to_fit_center );
	} else {
		in_pose.center();
	}
}

// XRW TEMP std::string
// XRW TEMP AlignChainMover::get_name() const {
// XRW TEMP  return AlignChainMover::mover_name();
// XRW TEMP }

moves::MoverOP
AlignChainMover::clone() const
{
	return moves::MoverOP( new AlignChainMover( *this ) );
}

moves::MoverOP
AlignChainMover::fresh_instance() const
{
	return moves::MoverOP( new AlignChainMover );
}

void
AlignChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	source_chain( tag->getOption< core::Size >( "source_chain", 0 ) );
	target_chain( tag->getOption< core::Size >( "target_chain", 0 ) );
	std::string const fname( tag->getOption< std::string >( "target_name" ) );

	align_to_com_ = tag->getOption< bool >( "align_to_com", "0" );
	if ( ! align_to_com_ ) {
		pose( core::import_pose::pose_from_file( fname ) );
	}

	TR << "source_chain: " << source_chain() << " target_chain: " << target_chain() << " pdb name: " << fname << " align_to_com: " << align_to_com_ << std::endl;
}

std::string AlignChainMover::get_name() const {
	return mover_name();
}

std::string AlignChainMover::mover_name() {
	return "AlignChain";
}

void AlignChainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"target_name", xs_string,
		"file name of the target pose on disk. "
		"The pose will be read just once at the start of the run and saved in memory to save I/O");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"source_chain", xsct_non_negative_integer,
		"the chain number in the working pose. 0 means the entire pose", "0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"target_chain", xsct_non_negative_integer,
		"the chain number in the target pose. 0 means the entire pose", "0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"align_to_com", xsct_rosetta_bool,
		"If true, will align pose to COM instead of target", "0" );
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Align a chain in the working pose to a chain in a pose on disk (CA superposition). "
		"All chains in the moving pose are rotated into the new coordinate frame, "
		"but the rotation is calculated on the specified chain. Specifying the 0th chain "
		"results in a whole-Pose alignment. target and source chains must have the same length of residues.",
		attlist );
}

std::string AlignChainMoverCreator::keyname() const {
	return AlignChainMover::mover_name();
}

protocols::moves::MoverOP
AlignChainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AlignChainMover );
}

void AlignChainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignChainMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
