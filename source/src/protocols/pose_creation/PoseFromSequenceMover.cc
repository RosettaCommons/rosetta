// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_creation/PoseFromSequenceMover.cc
/// @brief A class for generating a pose from a sequence/fasta
/// @author Dan Farrell (danpf@uw.edu)

// Unit headers
#include <protocols/pose_creation/PoseFromSequenceMover.hh>
#include <protocols/pose_creation/PoseFromSequenceMoverCreator.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "protocols.pose_creation.PoseFromSequenceMover" );

namespace protocols {
namespace pose_creation {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PoseFromSequenceMover::PoseFromSequenceMover():
	protocols::moves::Mover( PoseFromSequenceMover::mover_name() ),
	sequence_(""), residue_type_set_("fa_standard"), extended_(false)
{
}

/// @brief Default constructor
PoseFromSequenceMover::PoseFromSequenceMover(
	std::string const & sequence, std::string const & residue_type_set /*= "fa_standard"*/,
	bool const extended /*= false*/) :
	protocols::moves::Mover( PoseFromSequenceMover::mover_name() ),
	sequence_(sequence), residue_type_set_(residue_type_set), extended_(extended)
{
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PoseFromSequenceMover::PoseFromSequenceMover( PoseFromSequenceMover const & src ):
	protocols::moves::Mover(src),
	sequence_(src.get_sequence()),
	residue_type_set_(src.get_residue_type_set()),
	extended_(src.get_extended())
{
}

PoseFromSequenceMover & PoseFromSequenceMover::operator=( PoseFromSequenceMover const & src ) {
	set_sequence(src.get_sequence());
	set_residue_type_set(src.get_residue_type_set());
	set_extended(src.get_extended());
	return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PoseFromSequenceMover::~PoseFromSequenceMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
PoseFromSequenceMover::set_sequence(std::string const & sequence) {sequence_ = sequence;}

std::string const &
PoseFromSequenceMover::get_sequence() const {return sequence_;}

void
PoseFromSequenceMover::set_residue_type_set(std::string const & residue_type_set) {residue_type_set_ = residue_type_set;}

std::string const &
PoseFromSequenceMover::get_residue_type_set() const {return residue_type_set_;}

void
PoseFromSequenceMover::set_extended(bool const extended) {extended_ = extended;}

bool
PoseFromSequenceMover::get_extended() const {return extended_;}

/// @brief Apply the mover
void
PoseFromSequenceMover::apply( core::pose::Pose& pose ) {
	show(TR);
	pose.clear();
	if ( get_sequence().length() == 0 ) {
		TR.Warning << "Making a pose with sequence len 0 will just clear the pose!" << std::endl;
	}
	core::pose::make_pose_from_sequence(pose, get_sequence(), get_residue_type_set());

	if ( get_extended() ) {
		protocols::loops::Loops loops;
		protocols::loops::set_extended_torsions_and_idealize_loops(pose, loops);
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
PoseFromSequenceMover::show(std::ostream & output) const
{
	output << "extended: " << get_extended()
		<< ", residue_type_set: " << get_residue_type_set()
		<< ", seq: " << get_sequence() << ", ";
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PoseFromSequenceMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap&
) {
	std::string sequence("");
	if ( tag->hasOption( "fasta" ) && tag->hasOption( "sequence" ) ) {
		throw CREATE_EXCEPTION(
			utility::excn::RosettaScriptsOptionError,
			"PoseFromSequenceMover cannot create sequence if BOTH \"fasta\" and \"sequence\" tags are present");
	}
	if ( tag->hasOption( "fasta" ) ) {
		utility::vector1< core::sequence::SequenceOP > const sequences(core::sequence::read_fasta_file( tag->getOption< std::string >( "fasta" ) ));
		if ( tag->getOption<bool>( "use_all_in_fasta", false ) ) {
			for ( auto const & seqOP : sequences ) {
				sequence += seqOP->sequence() + "/";
			}
			sequence.pop_back();
		} else {
			sequence = sequences.front()->sequence();
		}
	} else if ( tag->hasOption( "sequence" ) ) {
		sequence = tag->getOption<std::string>( "sequence" );
	} else {
		throw CREATE_EXCEPTION(
			utility::excn::RosettaScriptsOptionError,
			"PoseFromSequenceMover requires \"fasta\" OR \"sequence\" tags to be set");
	}
	set_sequence(sequence);
	TR << "Sequence to build: " << get_sequence() << std::endl;

	if ( tag->hasOption("residue_type_set") ) {
		set_residue_type_set((tag->getOption<std::string>("residue_type_set")));
	}

	set_extended(tag->getOption<bool>("extended", false));

	TR << "Setting residue_type_set to: " << get_residue_type_set() << ", extended: " << get_extended() << std::endl;
}

void PoseFromSequenceMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using Attr = XMLSchemaAttribute;
	AttributeList attlist;
	attlist
		+ Attr( "fasta", xs_string, "Fasta file name" )
		+ Attr( "sequence", xs_string, "sequence as a string" )
		+ Attr::attribute_w_default( "use_all_in_fasta", xsct_rosetta_bool, "instead of using only first sequence in a fasta, use them all (and join them via '/')", "false" )
		+ Attr::attribute_w_default( "extended", xsct_rosetta_bool, "extend the pose", "false" )
		+ Attr::attribute_w_default( "residue_type_set", xs_string, "Residue representation, centroid, fa_standard", "fa_standard");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"This is a mover that can create a pose from a sequence, or fasta. if you"
		" use a fasta that has multiple sequences in it, and use the 'use_all_in_fasta'"
		" tag (set to true) it will use all of the sequences in the fasta and separate"
		" them by '/'.  Also you can build the pose extended if you set the 'extended'"
		" tag.",
		attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PoseFromSequenceMover::fresh_instance() const
{
	return utility::pointer::make_shared< PoseFromSequenceMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PoseFromSequenceMover::clone() const
{
	return utility::pointer::make_shared< PoseFromSequenceMover >( *this );
}

std::string PoseFromSequenceMover::get_name() const {
	return mover_name();
}

std::string PoseFromSequenceMover::mover_name() {
	return "PoseFromSequenceMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
PoseFromSequenceMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< PoseFromSequenceMover >();
}

std::string
PoseFromSequenceMoverCreator::keyname() const
{
	return PoseFromSequenceMover::mover_name();
}

void PoseFromSequenceMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PoseFromSequenceMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, PoseFromSequenceMover const & mover )
{
	mover.show(os);
	return os;
}

} //pose_creation
} //protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pose_creation::PoseFromSequenceMover::save( Archive & arc ) const {
	arc( CEREAL_NVP( sequence_ ) );
	arc( CEREAL_NVP( residue_type_set_ ) );
	arc( CEREAL_NVP( extended_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pose_creation::PoseFromSequenceMover::load( Archive & arc ) {
	arc( sequence_ );
	arc( residue_type_set_ );
	arc( extended_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pose_creation::PoseFromSequenceMover );
CEREAL_REGISTER_TYPE( protocols::pose_creation::PoseFromSequenceMover )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_pose_creation_PoseFromSequenceMover )
#endif // SERIALIZATION
