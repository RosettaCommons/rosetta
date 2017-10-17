// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/silent/SilentFileOptions.hh
/// @brief  Options for constructing a pose from a silent file
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <core/io/silent/SilentFileOptions.hh>

// Package Headers
#include <core/io/silent/SilentStructFactory.hh>

// Plaform Headers
#include <core/types.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace io {
namespace silent {

SilentFileOptions::SilentFileOptions()
{
	read_from_global_options();
}

SilentFileOptions::SilentFileOptions( utility::options::OptionCollection const & options )
{
	read_from_options( options );
}

void SilentFileOptions::read_from_global_options()
{
	read_from_options( basic::options::option );
}

void SilentFileOptions::read_from_options( utility::options::OptionCollection const & options )
{
	// NOTE: namespace basic::options is intentionally not using'ed here
	// to avoid accidental reads from the (evil) global option collection.
	using namespace basic::options::OptionKeys;

	keep_input_scores_ = options[ in::file::keep_input_scores ];
	out_user_tag_set_ = options[ out::user_tag ].user();
	out_user_tag_ = options[ out::user_tag ];
	out_weight_silent_scores_ = options[ out::file::weight_silent_scores ];
	in_silent_score_prefix_ = options[ in::file::silent_score_prefix ];
	in_silent_scores_wanted_set_ = options[ in::file::silent_scores_wanted ].user();
	if ( in_silent_scores_wanted_set_ ) in_silent_scores_wanted_ = options[ in::file::silent_scores_wanted ] ;

	in_fullatom_ = options[ in::file::fullatom ];
	write_failures_set_ = options[ run::write_failures ].user();
	write_failures_ = options[ run::write_failures ];

	in_silent_struct_type_set_ = options[ in::file::silent_struct_type ].user();
	in_silent_struct_type_ = options[ in::file::silent_struct_type ];
	out_silent_struct_type_set_ = options[ out::file::silent_struct_type ].user();
	out_silent_struct_type_ = options[ out::file::silent_struct_type ];

	read_through_errors_ = options[ in::file::silent_read_through_errors ];
	select_random_set_ = options[ in::file::silent_select_random ].user();
	select_random_ = options[ in::file::silent_select_random ];
	select_range_start_ = options[ in::file::silent_select_range_start ];
	select_range_len_ = options[ in::file::silent_select_range_len ];
	select_range_mul_ = options[ in::file::silent_select_range_mul ];
	force_silent_bitflip_on_read_set_ = options[ in::file::force_silent_bitflip_on_read ].user();
	force_silent_bitflip_on_read_ = options[ in::file::force_silent_bitflip_on_read ];
	print_all_score_headers_ = options[ out::file::silent_print_all_score_headers ];

}

void SilentFileOptions::read_from_tag( utility::tag::TagCOP tag )
{
	keep_input_scores_ = tag->getOption< bool >( "keep_input_scores", false );
	out_user_tag_set_ = tag->hasOption( "out_user_tag" );
	out_user_tag_ = tag->getOption< std::string >( "out_user_tag", "" );
	out_weight_silent_scores_ = tag->getOption< bool >( "out_weight_silent_scores", true );
	in_silent_score_prefix_ = tag->getOption< std::string >( "in_silent_score_prefix", "" );

	in_silent_scores_wanted_set_ = tag->hasOption( "silent_scores_wanted" );
	std::string cs_in_silent_scores_wanted = tag->getOption< std::string >( "silent_scores_wanted", "" );
	in_silent_scores_wanted_ = utility::string_split( cs_in_silent_scores_wanted, ',', std::string() );

	in_fullatom_ = tag->getOption< bool >( "in_fullatom", false );
	write_failures_set_ = tag->hasOption( "write_failures" );
	write_failures_ = tag->getOption< bool >( "write_failures", true );

	in_silent_struct_type_set_ = tag->hasOption( "in_silent_struct_type" );
	in_silent_struct_type_ = tag->getOption< std::string >( "in_silent_struct_type", "protein" );
	out_silent_struct_type_set_ = tag->hasOption( "out_silent_struct_type" );
	out_silent_struct_type_ = tag->getOption< std::string >( "out_silent_struct_type", "protein" );

	read_through_errors_ = tag->getOption< bool >( "read_through_errors", false );
	select_random_set_ = tag->hasOption( "select_random" );
	select_random_ = tag->getOption< int >( "select_random", 0 );
	select_range_start_ = tag->getOption< int >( "select_range_start", -1 );
	select_range_len_ = tag->getOption< int >( "select_range_len", 1 );
	select_range_mul_ = tag->getOption< int >( "select_range_mul", 1 );
	force_silent_bitflip_on_read_set_ = tag->hasOption( "force_silent_bitflip_on_read" );
	force_silent_bitflip_on_read_ = tag->getOption< bool >( "force_silent_bitflip_on_read", false );
	print_all_score_headers_ = tag->getOption< bool >( "print_all_score_headers", false );

}

void SilentFileOptions::list_read_options( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;

	read_options
		+ in::file::force_silent_bitflip_on_read
		+ in::file::fullatom
		+ in::file::keep_input_scores
		+ in::file::silent_read_through_errors
		+ in::file::silent_score_prefix
		+ in::file::silent_scores_wanted
		+ in::file::silent_select_random
		+ in::file::silent_select_range_len
		+ in::file::silent_select_range_mul
		+ in::file::silent_select_range_start
		+ in::file::silent_struct_type
		+ out::file::silent_print_all_score_headers
		+ out::file::silent_struct_type
		+ out::file::weight_silent_scores
		+ out::user_tag
		+ run::write_failures;
}

void SilentFileOptions::append_attributes_for_tag_parsing(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::AttributeList & attributes
)
{
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;

	// Create a restriction that lists all of the available silent-file types
	utility::vector1< std::string > ss_names = SilentStructFactory::get_instance()->get_ss_names();
	XMLSchemaRestriction silent_struct_type;
	silent_struct_type.name( "silent_struct_type" );
	silent_struct_type.base_type( xs_string );
	for ( auto const & ss : ss_names ) silent_struct_type.add_restriction( xsr_enumeration, ss );
	xsd.add_top_level_element( silent_struct_type );


	attributes
		+ Attr::attribute_w_default( "keep_input_scores", xsct_rosetta_bool, "Keep/Don't keep scores from input file in Pose", "false" )
		+ Attr( "out_user_tag", xs_string , "add this tag to structure tags: e.g., a process id" )
		+ Attr::attribute_w_default( "out_weight_silent_scores", xsct_rosetta_bool , "Weight scores in silent-file output; if false, raw, unweighted scores are output", "true" )
		+ Attr::attribute_w_default( "in_silent_score_prefix", xs_string , "Prefix that is appended to all scores read in from a silent-file; i.e."
		" allows you to rename existing scores to avoid name collisions when writing out the structures and their scores later", "" )
		+ Attr::attribute_w_default( "silent_scores_wanted", xs_string, "Comma separated list of the set of score terms to put into the Pose when reading it in", "false" )
		+ Attr::attribute_w_default( "in_fullatom", xsct_rosetta_bool, "Enable full-atom input of PDB or centroid structures; only used if"
		" the input silent file is ambiguous, otherwise, the information in the silent file itself is used to determine whether the input structures"
		" should be read as fullatom or centroid", "false" )
		+ Attr::attribute_w_default( "write_failures", xsct_rosetta_bool , "write failed structures to output", "true" )
		+ Attr::attribute_w_default( "in_silent_struct_type", xs_string, "Type of SilentStruct object to use in silent-file input. Only used if"
		" the input silent struct file does not describe its type in the REMARK lines at the top of the file (most do).", "protein" )
		+ Attr::attribute_w_default( "out_silent_struct_type", xs_string , "Type of SilentStruct object to use in silent-file output. The default"
		" 'protein' silent struct is only used if the structure is ideal; non-ideal structures are written out using a binary silent file;"
		" binary silent files are probably best in any case.", "protein" )
		+ Attr::attribute_w_default( "read_through_errors", xsct_rosetta_bool , "will ignore decoys with errors and continue reading", "false" )
		+ Attr::attribute_w_default( "select_random", xs_integer, "Instead of processing all of the structures in the input silent file, select a random"
		" subset of this number of decoys from every silent-file read; ignored if this value is zero", "0" )
		+ Attr::attribute_w_default( "select_range_start", xs_integer, "When given a postive value, then when jobs are being created based on the input"
		" silent struct file, instead of processing all of the Poses in the file, the first job will be taken as the nth Pose in the file,"
		" and then 'select_range_len' poses from there will be run. A negative value means consider the whole set", "-1" )
		+ Attr::attribute_w_default( "select_range_len", xs_integer, "For use with 'select_range_start': this says how large the block should"
		" be which should be run; a negative value means 'consider all decoys from the start of the"
		" block to the end'. A value of '1' (the default) means run a single job", "1" )
		+ Attr::attribute_w_default( "select_range_mul", xs_integer, "For use with 'select_range_start': Select a blocksize multiplier. This param"
		" basically multiplies 'select_range_start' when defining which structure to start on", "1" )
		+ Attr::attribute_w_default( "force_silent_bitflip_on_read", xsct_rosetta_bool, "Force bit-flipping when reading binary silent files."
		"  This is useful if files are produced on a little-endian system and read on a big-endian system.", "false" )
		+ Attr::attribute_w_default( "print_all_score_headers", xsct_rosetta_bool, "Print a SCORE header for every SilentStruct in a silent-file", "false" );
}

bool
SilentFileOptions::keep_input_scores() const
{
	return keep_input_scores_;
}

bool
SilentFileOptions::out_user_tag_set() const
{
	return out_user_tag_set_;
}

std::string
SilentFileOptions::out_user_tag() const
{
	return out_user_tag_;
}

bool
SilentFileOptions::out_weight_silent_scores() const
{
	return out_weight_silent_scores_;
}

std::string
SilentFileOptions::in_silent_score_prefix() const
{
	return in_silent_score_prefix_;
}

bool
SilentFileOptions::in_silent_scores_wanted_set() const
{
	return in_silent_scores_wanted_set_;
}

utility::vector1< std::string >
SilentFileOptions::in_silent_scores_wanted() const
{
	return in_silent_scores_wanted_;
}

bool
SilentFileOptions::in_fullatom() const
{
	return in_fullatom_;
}

bool
SilentFileOptions::write_failures_set() const
{
	return write_failures_set_;
}

bool
SilentFileOptions::write_failures() const
{
	return write_failures_;
}

bool
SilentFileOptions::in_silent_struct_type_set() const
{
	return in_silent_struct_type_set_;
}

std::string
SilentFileOptions::in_silent_struct_type() const
{
	return in_silent_struct_type_;
}

bool
SilentFileOptions::out_silent_struct_type_set() const
{
	return out_silent_struct_type_set_;
}

std::string
SilentFileOptions::out_silent_struct_type() const
{
	return out_silent_struct_type_;
}

bool
SilentFileOptions::read_through_errors() const
{
	return read_through_errors_;
}

bool
SilentFileOptions::select_random_set() const
{
	return select_random_set_;
}

int
SilentFileOptions::select_random() const
{
	return select_random_;
}

int
SilentFileOptions::select_range_start() const
{
	return select_range_start_;
}

int
SilentFileOptions::select_range_len() const
{
	return select_range_len_;
}

int
SilentFileOptions::select_range_mul() const
{
	return select_range_mul_;
}

bool SilentFileOptions::force_silent_bitflip_on_read_set() const
{
	return force_silent_bitflip_on_read_set_;
}


bool
SilentFileOptions::force_silent_bitflip_on_read() const
{
	return force_silent_bitflip_on_read_;
}

bool
SilentFileOptions::print_all_score_headers() const
{
	return print_all_score_headers_;
}

void
SilentFileOptions::keep_input_scores( bool setting )
{
	keep_input_scores_ = setting;
}

void
SilentFileOptions::out_user_tag( std::string setting )
{
	out_user_tag_set_ = true;
	out_user_tag_ = setting;
}

void
SilentFileOptions::out_weight_silent_scores( bool setting )
{
	out_weight_silent_scores_ = setting;
}

void
SilentFileOptions::in_silent_score_prefix( std::string setting )
{
	in_silent_score_prefix_ = setting;
}

void
SilentFileOptions::in_silent_scores_wanted( utility::vector1< std::string > const & setting )
{
	in_silent_scores_wanted_set_ = true;
	in_silent_scores_wanted_ = setting;
}

void
SilentFileOptions::in_fullatom( bool setting )
{
	in_fullatom_ = setting;
}

void
SilentFileOptions::write_failures( bool setting )
{
	write_failures_set_ = true;
	write_failures_ = setting;
}

void
SilentFileOptions::in_silent_struct_type( std::string setting )
{
	in_silent_struct_type_ = setting;
}

void
SilentFileOptions::out_silent_struct_type( std::string setting )
{
	out_silent_struct_type_ = setting;
}

void
SilentFileOptions::read_through_errors( bool setting )
{
	read_through_errors_ = setting;
}

void
SilentFileOptions::select_random( int setting )
{
	select_random_ = setting;
}

void
SilentFileOptions::select_range_start( int setting )
{
	select_range_start_ = setting;
}

void
SilentFileOptions::select_range_len( int setting )
{
	select_range_len_ = setting;
}

void
SilentFileOptions::select_range_mul( int setting )
{
	select_range_mul_ = setting;
}

void
SilentFileOptions::force_silent_bitflip_on_read( bool setting )
{
	force_silent_bitflip_on_read_set_ = true;
	force_silent_bitflip_on_read_ = setting;
}

void
SilentFileOptions::print_all_score_headers( bool setting )
{
	print_all_score_headers_ = setting;
}


} // namespace
} // namespace
} // namespace

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::silent::SilentFileOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( keep_input_scores_ ) ); // _Bool
	arc( CEREAL_NVP( out_user_tag_set_ ) ); // _Bool
	arc( CEREAL_NVP( out_user_tag_ ) ); // std::string
	arc( CEREAL_NVP( out_weight_silent_scores_ ) ); // _Bool
	arc( CEREAL_NVP( in_silent_score_prefix_ ) ); // std::string
	arc( CEREAL_NVP( in_silent_scores_wanted_set_ ) ); // _Bool
	arc( CEREAL_NVP( in_silent_scores_wanted_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( in_fullatom_ ) ); // _Bool
	arc( CEREAL_NVP( write_failures_set_ ) ); // _Bool
	arc( CEREAL_NVP( write_failures_ ) ); // _Bool
	arc( CEREAL_NVP( in_silent_struct_type_set_ ) ); // _Bool
	arc( CEREAL_NVP( in_silent_struct_type_ ) ); // std::string
	arc( CEREAL_NVP( out_silent_struct_type_set_ ) ); // _Bool
	arc( CEREAL_NVP( out_silent_struct_type_ ) ); // std::string
	arc( CEREAL_NVP( read_through_errors_ ) ); // _Bool
	arc( CEREAL_NVP( select_random_set_ ) ); // _Bool
	arc( CEREAL_NVP( select_random_ ) ); // int
	arc( CEREAL_NVP( select_range_start_ ) ); // int
	arc( CEREAL_NVP( select_range_len_ ) ); // int
	arc( CEREAL_NVP( select_range_mul_ ) ); // int
	arc( CEREAL_NVP( force_silent_bitflip_on_read_set_ ) ); // _Bool
	arc( CEREAL_NVP( force_silent_bitflip_on_read_ ) ); // _Bool
	arc( CEREAL_NVP( print_all_score_headers_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::silent::SilentFileOptions::load( Archive & arc ) {
	arc( keep_input_scores_ ); // _Bool
	arc( out_user_tag_set_ ); // _Bool
	arc( out_user_tag_ ); // std::string
	arc( out_weight_silent_scores_ ); // _Bool
	arc( in_silent_score_prefix_ ); // std::string
	arc( in_silent_scores_wanted_set_ ); // _Bool
	arc( in_silent_scores_wanted_ ); // utility::vector1<std::string>
	arc( in_fullatom_ ); // _Bool
	arc( write_failures_set_ ); // _Bool
	arc( write_failures_ ); // _Bool
	arc( in_silent_struct_type_set_ ); // _Bool
	arc( in_silent_struct_type_ ); // std::string
	arc( out_silent_struct_type_set_ ); // _Bool
	arc( out_silent_struct_type_ ); // std::string
	arc( read_through_errors_ ); // _Bool
	arc( select_random_set_ ); // _Bool
	arc( select_random_ ); // int
	arc( select_range_start_ ); // int
	arc( select_range_len_ ); // int
	arc( select_range_mul_ ); // int
	arc( force_silent_bitflip_on_read_set_ ); // _Bool
	arc( force_silent_bitflip_on_read_ ); // _Bool
	arc( print_all_score_headers_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::silent::SilentFileOptions );
CEREAL_REGISTER_TYPE( core::io::silent::SilentFileOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_io_silent_SilentFileOptions )
#endif // SERIALIZATION
