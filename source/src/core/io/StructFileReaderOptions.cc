// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/StructFileReaderOptions.cc
/// @brief  Definitions for StructFileReaderOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/io/StructFileReaderOptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace io {

StructFileReaderOptions::StructFileReaderOptions() :
	StructFileRepOptions()
{
	init_from_options( basic::options::option );
}

StructFileReaderOptions::StructFileReaderOptions( utility::options::OptionCollection const & options ) :
	StructFileRepOptions( options )
{
	init_from_options( options );
}

StructFileReaderOptions::~StructFileReaderOptions() = default;

void StructFileReaderOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	StructFileRepOptions::parse_my_tag( tag );

	set_new_chain_order(  tag->getOption< bool >( "new_chain_order", false ) );
	set_obey_ENDMDL(  tag->getOption< bool >( "obey_ENDMDL", false ) );
	set_read_pdb_header( tag->getOption< bool >( "preserve_header", false ) );
	set_glycam_pdb_format( tag->getOption< bool >( "glycam_pdb_format", false ) );
}

std::string StructFileReaderOptions::type() const { return "StructFileReaderOptions"; }

// accessors
bool StructFileReaderOptions::new_chain_order() const { return new_chain_order_; }
bool StructFileReaderOptions::obey_ENDMDL() const { return obey_ENDMDL_; }
bool StructFileReaderOptions::read_pdb_header() const { return read_pdb_header_; }
bool StructFileReaderOptions::glycam_pdb_format() const { return glycam_pdb_format_; }

// mutators
void StructFileReaderOptions::set_new_chain_order( bool setting ) { new_chain_order_ = setting; }
void StructFileReaderOptions::set_obey_ENDMDL( bool setting ) { obey_ENDMDL_ = setting; }
void StructFileReaderOptions::set_read_pdb_header( bool setting ) { read_pdb_header_ = setting; }
void StructFileReaderOptions::set_glycam_pdb_format( bool setting ) { glycam_pdb_format_ = setting; }

void
StructFileReaderOptions::list_options_read( utility::options::OptionKeyList & read_options )
{
	using namespace basic::options::OptionKeys;
	StructFileRepOptions::list_options_read( read_options );
	read_options
		+ in::file::new_chain_order
		+ in::file::obey_ENDMDL
		+ run::preserve_header
		+ carbohydrates::glycam_pdb_format;
}

void
StructFileReaderOptions::append_schema_attributes( utility::tag::AttributeList & attributes )
{
	StructFileRepOptions::append_schema_attributes( attributes );
	using namespace utility::tag;
	typedef XMLSchemaAttribute Attr;

	attributes
		+ Attr::attribute_w_default( "new_chain_order", xsct_rosetta_bool, "TO DO", "false" )
		+ Attr::attribute_w_default( "obey_ENDMDL", xsct_rosetta_bool, "TO DO", "false" )
		// + Attr::attribute_w_default( "preserve_header", xsct_rosetta_bool, "TO DO", "false" ) // already done by base class
		+ Attr::attribute_w_default( "glycam_pdb_format", xsct_rosetta_bool, "TO DO", "false" );

}


void StructFileReaderOptions::init_from_options( utility::options::OptionCollection const & options )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_new_chain_order( options[ in::file::new_chain_order ]() );
	set_obey_ENDMDL( options[ in::file::obey_ENDMDL ].value() );
	set_read_pdb_header( options[ run::preserve_header ]() );
	set_glycam_pdb_format( options[ carbohydrates::glycam_pdb_format ]() );
}

bool
StructFileReaderOptions::operator == ( StructFileReaderOptions const & other ) const
{
	if ( ! StructFileRepOptions::operator==( other ) ) return false;

	if ( new_chain_order_ != other.new_chain_order_ ) return false;
	if ( obey_ENDMDL_ != other.obey_ENDMDL_  ) return false;
	if ( read_pdb_header_ != other.read_pdb_header_  ) return false;
	if ( glycam_pdb_format_ != other.glycam_pdb_format_  ) return false;

	return true;
}

bool
StructFileReaderOptions::operator < ( StructFileReaderOptions const & other ) const
{
	if ( StructFileRepOptions::operator< ( other ) ) return true;
	if ( StructFileRepOptions::operator== ( other ) ) return false;

	if ( new_chain_order_ <  other.new_chain_order_ ) return true;
	if ( new_chain_order_ == other.new_chain_order_ ) return false;
	if ( obey_ENDMDL_ <  other.obey_ENDMDL_  ) return true;
	if ( obey_ENDMDL_ == other.obey_ENDMDL_  ) return false;
	if ( read_pdb_header_ <  other.read_pdb_header_  ) return true;
	if ( read_pdb_header_ == other.read_pdb_header_  ) return false;
	if ( glycam_pdb_format_ <  other.glycam_pdb_format_  ) return true;
	//if ( glycam_pdb_format_ == other.glycam_pdb_format_  ) return false;
	return false;

}


} // namespace io
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::io::StructFileReaderOptions::save( Archive & arc ) const {
	arc( cereal::base_class< core::io::StructFileRepOptions >( this ) );
	arc( CEREAL_NVP( new_chain_order_ ) ); // _Bool
	arc( CEREAL_NVP( obey_ENDMDL_ ) ); // _Bool
	arc( CEREAL_NVP( read_pdb_header_ ) ); // _Bool
	arc( CEREAL_NVP( glycam_pdb_format_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::io::StructFileReaderOptions::load( Archive & arc ) {
	arc( cereal::base_class< core::io::StructFileRepOptions >( this ) );
	arc( new_chain_order_ ); // _Bool
	arc( obey_ENDMDL_ ); // _Bool
	arc( read_pdb_header_ ); // _Bool
	arc( glycam_pdb_format_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::io::StructFileReaderOptions );
CEREAL_REGISTER_TYPE( core::io::StructFileReaderOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_io_StructFileReaderOptions )
#endif // SERIALIZATION
