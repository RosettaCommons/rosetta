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
#include <core/chemical/ResidueTypeBase.hh>
#include <core/chemical/ResidueProperties.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

// C++ headers
#include <map>
#include <iostream>
#include <algorithm>

// Construct tracer.
static basic::Tracer TR( "core.chemical.ResidueProperties" );


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

using namespace core;

// Public methods /////////////////////////////////////////////////////////////
// Standard methods ///////////////////////////////////////////////////////////

ResidueProperties::ResidueProperties() :
	utility::VirtualBase(),
	general_property_status_( N_PROPERTIES, false ),
	variant_type_status_( N_VARIANTS, false )
{}


// Constructor with owning ResidueType
ResidueProperties::ResidueProperties( ResidueTypeBase const & residue_type ) :
	utility::VirtualBase(),
	parent_residue_type_( residue_type.name() ),
	general_property_status_( N_PROPERTIES, false ),
	variant_type_status_( N_VARIANTS, false )
{}

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
	// YES, we should. Done.  --Vikram

	ResidueProperty const prop( get_property_from_string( property ) );
	if ( prop == NO_PROPERTY ) {
		if ( TR.Warning.visible() ) TR.Warning << "Error in core::chemical::ResidueProperties::set_property(): Rosetta does not recognize the property \"" << property << "\".  Has it been added to the \"general_properties.list\" file?" << std::endl;
		return false;
	}
	return general_property_status_[ prop ];
}

void
ResidueProperties::set_property( std::string const & property, bool const setting )
{
	using namespace std;
	using namespace utility::excn;

	// should we use NO_PROPERTY return value instead of exception catching? Would aid in gdb/lldb debugging to not use exceptions for this. See set_variant_type. -- rhiju
	// YES, we should. Done.  --Vikram

	ResidueProperty const prop( get_property_from_string( property ) );
	runtime_assert_string_msg( prop != NO_PROPERTY, "Error in core::chemical::ResidueProperties::set_property(): Rosetta does not recognize the property \"" + property + "\".  Has it been added to the \"general_properties.list\" file?" );
	general_property_status_[ prop ] = setting;
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
						" already exists in " << parent_residue_type_ << endl;
				} else {
					TR.Trace << "Adding the custom variant " << variant_type <<
						" to " << parent_residue_type_ << endl;
					custom_variant_types_.push_back( variant_type );
				}
			} else /* setting == false */ {
				auto i =
					find( custom_variant_types_.begin(), custom_variant_types_.end(), variant_type );
				if ( i == custom_variant_types_.end() ) {
					utility_exit_with_message( "Rosetta does not recognize the custom variant " + variant_type +
						" in " + parent_residue_type_ );
				} else {
					TR.Trace << "Removing the custom variant " << variant_type <<
						" from " << parent_residue_type_ << endl;
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

bool
ResidueProperties::operator==( ResidueProperties const & other ) const
{
	return parent_residue_type_ == other.parent_residue_type_ &&
		general_property_status_ == other.general_property_status_ &&
		variant_type_status_ == other.variant_type_status_ &&
		has_custom_variant_types_ == other.has_custom_variant_types_ &&
		numeric_properties_ == other.numeric_properties_ &&
		string_properties_ == other.string_properties_ &&
		custom_variant_types_.size() == other.custom_variant_types_.size() &&
		std::is_permutation(custom_variant_types_.begin(), custom_variant_types_.end(), other.custom_variant_types_.begin()); // Need to compare in order-independent fashion.
}

// Helper methods /////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that ResidueProperties can be "printed" in PyRosetta).
std::ostream &
operator<<( std::ostream & output, ResidueProperties const & object_to_output )
{
	object_to_output.show(output);
	return output;
}

std::ostream &
operator<<( std::ostream & output, ResidueProperty const & object_to_output )
{
	output << ResidueProperties::get_string_from_property( object_to_output );
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

ResiduePropertiesOP
deep_copy( ResidueProperties const & source) {
	return utility::pointer::make_shared< ResidueProperties >( source );
}

}  // namespace chemical
}  // namespace core



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::ResidueProperties::save( Archive & arc ) const {
	arc( CEREAL_NVP( parent_residue_type_ ) );
	arc( CEREAL_NVP( general_property_status_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( variant_type_status_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( has_custom_variant_types_ ) ); // _Bool
	arc( CEREAL_NVP( custom_variant_types_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( numeric_properties_ ) ); // std::map<std::string, core::Real>
	arc( CEREAL_NVP( string_properties_ ) ); // std::map<std::string, std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueProperties::load( Archive & arc ) {
	arc( parent_residue_type_ ); // std::string
	arc( general_property_status_ ); // utility::vector1<_Bool>
	arc( variant_type_status_ ); // utility::vector1<_Bool>
	arc( has_custom_variant_types_ ); // _Bool
	arc( custom_variant_types_ ); // utility::vector1<std::string>
	arc( numeric_properties_ ); // std::map<std::string, core::Real>
	arc( string_properties_ ); // std::map<std::string, std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueProperties );
CEREAL_REGISTER_TYPE( core::chemical::ResidueProperties )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueProperties )
#endif // SERIALIZATION
