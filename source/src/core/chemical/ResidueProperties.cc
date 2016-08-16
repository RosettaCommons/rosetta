// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/ResidueProperties.cc
/// @brief   Method definitions for ResidueProperties.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit header
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

// C++ headers
#include <map>
#include <iostream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.chemical.ResidueProperties" );


namespace core {
namespace chemical {

using namespace core;

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////
// Constructor with owning ResidueType
ResidueProperties::ResidueProperties( ResidueType const * residue_type ) : utility::pointer::ReferenceCount()
{
	init( residue_type );
}

// "Copy constructor"
ResidueProperties::ResidueProperties( ResidueProperties const & object_to_copy, ResidueType const * new_owner ) :
	utility::pointer::ReferenceCount( object_to_copy )
{
	residue_type_ = new_owner;
	copy_data( *this, object_to_copy );
}

// Destructor
ResidueProperties::~ResidueProperties() {}


// Standard Rosetta methods ///////////////////////////////////////////////////
// General methods
void
ResidueProperties::show( std::ostream & output ) const
{
	using namespace std;
	using namespace utility;

	output << " Properties:";
	vector1< string > const & properties( get_list_of_properties() );
	Size const n_properties( properties.size() );
	for ( core::uint i = 1; i <= n_properties; ++i ) {
		output << ' ' << properties[ i ];
	}
	output << endl;

	output << " Variant types:";
	vector1< string > const & variants( get_list_of_variants() );
	Size const n_variants( variants.size() );
	for ( core::uint i = 1; i <= n_variants; ++i ) {
		output << ' ' << variants[ i ];
	}
	output << endl;
}


// Accessors/Mutators
bool
ResidueProperties::has_property( std::string const & property ) const
{
	using namespace std;
	using namespace utility::excn;

	// should we use NO_PROPERTY return value instead of exception catching? Would aid in gdb/lldb debugging to not use exceptions for this. See has_variant_type. -- rhiju
	try {
		return general_property_status_[ get_property_from_string( property ) ];
	} catch ( EXCN_Msg_Exception const & e ) {
		TR.Warning << e.msg() << endl;
	}
	return false;
}

void
ResidueProperties::set_property( std::string const & property, bool const setting )
{
	using namespace std;
	using namespace utility::excn;

	// should we use NO_PROPERTY return value instead of exception catching? Would aid in gdb/lldb debugging to not use exceptions for this. See set_variant_type. -- rhiju
	try {
		general_property_status_[ get_property_from_string( property ) ] = setting;
	} catch ( EXCN_Msg_Exception const & e ) {
		utility_exit_with_message( e.msg() );
	}
}


bool
ResidueProperties::is_variant_type( std::string const & variant_type ) const
{
	using namespace std;
	using namespace utility::excn;

	VariantType vtype = get_variant_from_string( variant_type );
	if ( vtype != NO_VARIANT ) {
		return variant_type_status_[ vtype ];
	} else {
		if ( ! has_custom_variant_types_ ) {
			TR.Debug << "Rosetta does not recognize the variant: " << variant_type << "; has it been added to variant_types.list?" << std::endl;
		} else {
			return custom_variant_types_.has_value( variant_type );
		}
	}
	return false;
}

void
ResidueProperties::set_variant_type( std::string const & variant_type, bool const setting )
{
	using namespace std;
	using namespace utility;
	using namespace utility::excn;


	VariantType vtype = get_variant_from_string( variant_type );
	if ( vtype != NO_VARIANT ) {
		variant_type_status_[ vtype ] = setting;
	} else {
		if ( ! has_custom_variant_types_ ) {
			utility_exit_with_message( "Rosetta does not recognize the variant: " + variant_type + "; has it been added to variant_types.list?" );
		} else {
			if ( setting /* == true */ ) {
				if ( custom_variant_types_.has_value( variant_type ) ) {
					TR.Trace << "Custom variant " << variant_type <<
						" already exists in " << residue_type_->name() << endl;
				} else {
					TR.Trace << "Adding the custom variant " << variant_type <<
						" to " << residue_type_->name() << endl;
					custom_variant_types_.push_back( variant_type );
				}
			} else /* setting == false */ {
				vector1< string >::iterator i =
					find( custom_variant_types_.begin(), custom_variant_types_.end(), variant_type );
				if ( i == custom_variant_types_.end() ) {
					utility_exit_with_message( "Rosetta does not recognize the custom variant " + variant_type +
						" in " + residue_type_->name() );
				} else {
					TR.Trace << "Removing the custom variant " << variant_type <<
						" from " << residue_type_->name() << endl;
					custom_variant_types_.erase( i );
				}
			}
		}
	}
}


/// @note  This function was copied from its old location in ResidueType.
void
ResidueProperties::add_numeric_property( std::string const & tag, core::Real const value )
{
	using namespace std;

	numeric_properties_.insert( make_pair( tag, value ) );
}

/// @note  This function was copied from its old location in ResidueType.
void
ResidueProperties::add_string_property( std::string const & tag, std::string const & value )
{
	using namespace std;

	string_properties_.insert( make_pair( tag, value ) );
}


// Generate and return a list of strings representing the properties of this ResidueType.
utility::vector1< std::string >
ResidueProperties::get_list_of_properties() const
{
	using namespace std;
	using namespace utility;

	vector1< string > list;

	for ( ResidueProperty property = FIRST_PROPERTY; property <= N_PROPERTIES; ++property ) {
		if ( general_property_status_[ property ] ) {
			list.push_back( get_string_from_property( property ) );
		}
	}

	return list;
}

/// @brief Return a list of VariantType enums for this ResidueType.
/// @details This will not include custom, string-based variant types generated on the fly.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< VariantType >
ResidueProperties::get_list_of_variant_enums() const {
	utility::vector1< VariantType > output_list;
	for ( VariantType variant = FIRST_VARIANT; variant <= N_VARIANTS; ++variant ) {
		if ( variant_type_status_[ variant ] ) output_list.push_back(variant);
	}
	return output_list;
}

/// @brief Get a const-access reference to the list of custom VariantType strings for this ResidueType.
/// @details This will not include enum-based standard variants.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
utility::vector1< std::string > const &
ResidueProperties::get_list_of_custom_variants_by_reference() const {
	return custom_variant_types_;
}


/// @brief Generate and return a list of strings representing the VariantTypes of this ResidueType.
utility::vector1< std::string >
ResidueProperties::get_list_of_variants() const
{
	using namespace std;
	using namespace utility;

	vector1< string > list;

	for ( VariantType variant = FIRST_VARIANT; variant <= N_VARIANTS; ++variant ) {
		if ( variant_type_status_[ variant ] ) {
			list.push_back( get_string_from_variant( variant ) );
		}
	}

	// Consider "custom" variants if necessary.
	if ( has_custom_variant_types_ ) {
		Size const n_custom_variants( custom_variant_types_.size() );
		for ( core::uint i = 1; i <= n_custom_variants; ++i ) {
			list.push_back( custom_variant_types_[ i ] );
		}
	}

	return list;
}


// Private methods ////////////////////////////////////////////////////////////
// Initialize data members.
void
ResidueProperties::init( ResidueType const * residue_type )
{
	residue_type_ = residue_type;
	general_property_status_.resize( N_PROPERTIES, false );
	variant_type_status_.resize( N_VARIANTS, false );
	has_custom_variant_types_ = false;
}

// Copy all data members from <from> to <to>.
void
ResidueProperties::copy_data( ResidueProperties & to, ResidueProperties const & from )
{
	to.general_property_status_ = from.general_property_status_;
	to.variant_type_status_ = from.variant_type_status_;
	to.has_custom_variant_types_ = from.has_custom_variant_types_;
	to.custom_variant_types_ = from.custom_variant_types_;
	to.numeric_properties_ = from.numeric_properties_;
	to.string_properties_ = from.string_properties_;
}


// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that ResidueProperties can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, ResidueProperties const & object_to_output )
{
	object_to_output.show(output);
	return output;
}

// This allows one to use a for loop with ResidueProperty enum values.
// This is safe, because the ResidueProperty values are set automatically before compiling.
ResidueProperty &
operator++( ResidueProperty & property )
{
	property = static_cast< ResidueProperty >( static_cast< int >( property ) + 1 );
	return property;
}

// This allows one to use a for loop with VariantType enum values.
// This is safe, because the VariantType values are set automatically before compiling.
VariantType &
operator++( VariantType & variant )
{
	variant = static_cast< VariantType >( static_cast< int >( variant ) + 1 );
	return variant;
}

}  // namespace chemical
}  // namespace core


