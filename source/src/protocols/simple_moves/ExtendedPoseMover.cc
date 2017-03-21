// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ExtendedPoseMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/simple_moves/ExtendedPoseMover.hh>

// Utility headers
#include <utility>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// C/C++ headers
#include <string>

#include <utility/exit.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// Option Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/excn/Exceptions.hh>
#include <core/kinematics/Jump.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/simple_moves/ExtendedPoseMoverCreator.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.simple_moves.ExtendedPoseMover" );

namespace protocols {
namespace simple_moves {

ExtendedPoseMover::ExtendedPoseMover(std::string  sequence,
	std::string  residue_type_set)
: sequence_(std::move(sequence)), residue_type_set_(std::move(residue_type_set)) {}

bool ExtendedPoseMover::valid() const {
	return sequence() != "";
}

void ExtendedPoseMover::apply(core::pose::Pose& pose) {
	// Ensure that this instance is in a valid state
	debug_assert(valid());
	pose.clear();

	protocols::loops::Loops loops;
	core::pose::make_pose_from_sequence(pose, sequence_, residue_type_set_);
	protocols::loops::set_extended_torsions_and_idealize_loops(pose, loops);

	if ( chain() != "" ) {
		core::pose::PDBInfoOP info( new core::pose::PDBInfo( pose, true ) );
		info->set_chains( (char) chain()[0] );
		pose.pdb_info( info );
	}
}

// XRW TEMP std::string ExtendedPoseMover::get_name() const {
// XRW TEMP  return "ExtendedPoseMover";
// XRW TEMP }

const std::string& ExtendedPoseMover::sequence() const {
	return sequence_;
}

const std::string& ExtendedPoseMover::residue_type_set() const {
	return residue_type_set_;
}

void ExtendedPoseMover::sequence(const std::string& sequence) {
	sequence_ = sequence;
}

void ExtendedPoseMover::residue_type_set(const std::string& residue_type_set) {
	residue_type_set_ = residue_type_set;
}

protocols::moves::MoverOP ExtendedPoseMover::clone() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::ExtendedPoseMover(*this) );
}

protocols::moves::MoverOP ExtendedPoseMover::fresh_instance() const {
	return protocols::moves::MoverOP( new protocols::simple_moves::ExtendedPoseMover() );
}

void ExtendedPoseMover::parse_my_tag(const utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	const protocols::filters::Filters_map&,
	const protocols::moves::Movers_map&,
	const core::pose::Pose&) {
	using namespace core::sequence;

	// required options
	if ( tag->hasOption( "fasta" ) ) {
		utility::vector1< SequenceOP > sequences = read_fasta_file( tag->getOption< string >( "fasta" ) );
		sequence( sequences.front()->sequence() );
		if ( sequences.size() > 1 ) {
			std::ostringstream ss;
			ss << "In " << tag->getOption< string >( "name" ) << ": The fasta file " << tag->getOption< string >( "fasta" ) << " contained >1 sequence; using the first one only. Use multiple ExtendedPoseMovers for multiple chains." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
		}
	} else if ( tag->hasOption( "sequence" ) ) {
		sequence( tag->getOption<string>( "sequence" ) );
	} else if ( tag->getOption<bool>("use_fasta", false ) ) {
		string filename = basic::options::option[ basic::options::OptionKeys::in::file::fasta ]()[1];
		sequence( read_fasta_file( filename ).front()->sequence() );
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption("Failed to specify required option `sequence` or fasta file");
	}

	// additional options
	if ( tag->hasOption("residue_type_set") ) {
		residue_type_set(tag->getOption<string>("residue_type_set"));
	}

	chain( tag->getOption< string >( "chain", "" ) );
	if ( chain().length() != 1 ) {
		throw utility::excn::EXCN_RosettaScriptsOption( chain()+" is an invalid chain code in "+tag->getOption< string >( "name" )+"." );
	}
}

std::string ExtendedPoseMover::get_name() const {
	return mover_name();
}

std::string ExtendedPoseMover::mover_name() {
	return "ExtendedPoseMover";
}

void ExtendedPoseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;
	AttributeList attlist;
	attlist
		+ Attr( "fasta", xs_string, "Fasta file name (complete sequence?)" )
		+ Attr( "sequence", xs_string, "Sequence string" )
		+ Attr::attribute_w_default( "use_fasta", xsct_rosetta_bool, "Take sequence from fasta file as specified in comandline option \"-fasta\"", "false" )
		+ Attr( "residue_type_set", xs_string, "Residue representation, centroid, fa_standard, centroid_rot ...." )
		+ Attr( "chain", xsct_char, "Chain ID" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Creates an extended, idealized pose from the sequence and"
		"residue type set specified in the constructor.", attlist );
}

std::string ExtendedPoseMoverCreator::keyname() const {
	return ExtendedPoseMover::mover_name();
}

protocols::moves::MoverOP
ExtendedPoseMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ExtendedPoseMover );
}

void ExtendedPoseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtendedPoseMover::provide_xml_schema( xsd );
}


}  // namespace simple_moves
}  // namespace protocols
