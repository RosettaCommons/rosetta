// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/ResidueProperties.hh
/// @brief   Declarations and simple accessor/mutator definitions for ResidueProperties.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_core_chemical_ResidueProperties_HH
#define INCLUDED_core_chemical_ResidueProperties_HH

// Unit headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueProperties.fwd.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <iostream>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

/// @details This is a container class for the large assortment of properties associated with ResidueTypes.
/// It prevents ResidueType from becoming cluttered with an over-abundance of properties and related methods.
/// @remarks This is the first step in a major refactor of how properties are handled in Rosetta.
/// For now, I have just gathered all the properties related code into one place, so that changes to the system can be
/// more readily made.  Previous behavior has been maintained.  In the future, I have several ideas for further
/// improvements. ~Labonte
class ResidueProperties : public utility::pointer::ReferenceCount {
public:  // Standard methods //////////////////////////////////////////////////
	/// @brief  Constructor with owning ResidueType
	// Note: Can't use ResidueTypeCAP here because ResidueProperties is created in the ResidueType c'tor,
	// where the self weak pointer is not available yet. residue_type_ here is used only for Tracer output.
	ResidueProperties( ResidueType const * residue_type );

	/// @brief  Copy constructor
	ResidueProperties( ResidueProperties const & object_to_copy, ResidueType const * new_owner );

	// Destructor
	~ResidueProperties() override;

private:
	// Empty constructor -- should not ever be called; not implemented
	// An instance of ResidueProperties must have a pointer back to its owning ResidueType.
	ResidueProperties();

	// Copy constructor -- should not ever be called; not implemented
	// An instance of ResidueProperties must have a pointer back to its owning ResidueType.
	ResidueProperties( ResidueProperties const & object_to_copy );

	// Assignment operator -- should not ever be called; not implemented
	// An instance of ResidueProperties must have a pointer back to its owning ResidueType.
	ResidueProperties & operator=( ResidueProperties const & object_to_copy );


public:  // Standard Rosetta methods //////////////////////////////////////////
	/// @brief  Generate string representation of ResidueProperties for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;


public:  // Accessors/Mutators ////////////////////////////////////////////////
	/// @brief  Get whether or not this ResidueType has the requested property.
	inline bool
	has_property( ResidueProperty const property ) const
	{
		return general_property_status_[ property ];
	}

	/// @brief  Get whether or not this ResidueType has the requested property by string.
	bool has_property( std::string const & property ) const;


	/// @brief  Set the status of the given property for this ResidueType.
	void
	set_property( ResidueProperty const property, bool const setting )
	{
		general_property_status_[ property ] = setting;
	}

	/// @brief  Set the status of the given property for this ResidueType by string.
	void set_property( std::string const & property, bool const setting );


	/// @brief  Get whether or not this ResidueType is of the requested VariantType.
	bool
	is_variant_type( VariantType const variant_type ) const
	{
		// If this is ever called with NO_VARIANT, problems!
		if ( variant_type == NO_VARIANT ) {
			for ( Size ii = 1; ii <= variant_type_status_.size(); ++ii ) {
				if ( variant_type_status_[ ii ] ) return false;
			}
			return true;
		}
		return variant_type_status_[ variant_type ];
	}

	/// @brief  Get whether or not this ResidueType is of the requested VariantType by string.
	bool is_variant_type( std::string const & variant_type ) const;

	/// @brief Set the status of a given VariantType for this ResidueType.
	void
	set_variant_type( VariantType const variant_type, bool const setting )
	{
		variant_type_status_[ variant_type ] = setting;
	}

	/// @brief Set the status of a given VariantType for this ResidueType by string.
	void set_variant_type( std::string const & variant_type, bool const setting );

	/// @brief Does this ResidueType contain additional VariantTypes than the standard list?
	bool
	has_custom_variant_types() const
	{
		return has_custom_variant_types_;
	}

	/// @brief   Turn on the ability to create VariantTypes "on-the-fly".
	/// @details Custom" VariantTypes as strings are permitted for the enzdes and metalloproteins cases.
	/// Do not enable unless you have a good reason to, as string look-ups are less efficient and more error-prone.
	void
	enable_custom_variant_types()
	{
		has_custom_variant_types_ = true;
	}

	/// @brief  Add a numeric property.
	void add_numeric_property( std::string const & tag, core::Real const value );

	/// @brief  Add a string property.
	void add_string_property( std::string const & tag, std::string const & value );

	std::map< std::string, core::Real > const &
	numeric_properties() const
	{
		return numeric_properties_;
	}

	std::map< std::string, std::string > const &
	string_properties() const
	{
		return string_properties_;
	}


public:  // Other public methods //////////////////////////////////////////////
	/// @brief  Generate and return a list of strings representing the properties of this ResidueType.
	utility::vector1< std::string > get_list_of_properties() const;

	/// @brief Return a list of VariantType enums for this ResidueType.
	/// @details This will not include custom, string-based variant types generated on the fly.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	utility::vector1< VariantType > get_list_of_variant_enums() const;

	/// @brief Get a const-access reference to the list of custom VariantType strings for this ResidueType.
	/// @details This will not include enum-based standard variants.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	utility::vector1< std::string > const & get_list_of_custom_variants_by_reference() const;

	/// @brief  Generate and return a list of strings representing the VariantTypes of this ResidueType.
	/// @details This will include both custom, string-based variants made on-the-fly AND the standard
	/// variants that are enumerated.
	utility::vector1< std::string > get_list_of_variants() const;

	/// @brief  Return a list of custom VariantTypes only for this ResidueType.
	utility::vector1< std::string >
	get_list_of_custom_variants() const
	{
		return custom_variant_types_;
	}

	// TODO: Add n_properties() and size() and search for get_list_of_properties().size() calls.
	// TODO: Add n_variants() and search for get_list_of_variants().size() calls.
	// TODO: Add is_not_a_variant() and search for get_list_of_variants().empty()  calls.

	/// @brief Generate the map of (string->VariantType).
	/// @details Used to initialize the const static map.  Since c++11 initialization of const static singletons
	/// is inherently threadsafe, we don't subsequently have to worry about mutexes for access to this.
	/// @note This is implemented in core/chemical/VariantType_mappings.cc, an automatically-generated file.  DO NOT EDIT
	/// THAT FILE DIRECTLY.  Intead, edit the generating script.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	static std::map< std::string, VariantType > const * generate_string_to_variant_map();

	/// @brief Generate the map of (string->ResidueProperty).
	/// @details Used to initialize the const static map.  Since c++11 initialization of const static singletons
	/// is inherently threadsafe, we don't subsequently have to worry about mutexes for access to this.
	/// @note This is implemented in core/chemical/ResidueProperty_mappings.cc, an automatically-generated file.  DO NOT EDIT
	/// THAT FILE DIRECTLY.  Intead, edit the generating script.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	static std::map< std::string, ResidueProperty > const * generate_string_to_property_map();

	/// @brief Generate the map of (VariantType->string).  (This is actually a vector since we're going from an enum -- which can be
	/// a vector index -- to a string).
	/// @details Used to initialize the const static vector.  Since c++11 initialization of const static singletons
	/// is inherently threadsafe, we don't subsequently have to worry about mutexes for access to this.
	/// @note This is implemented in core/chemical/VariantType_mappings.cc, an automatically-generated file.  DO NOT EDIT
	/// THAT FILE DIRECTLY.  Intead, edit the generating script.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	static utility::vector1< std::string > const * generate_variant_to_string_vector();

	/// @brief Generate the map of (ResidueProperty->string).  (This is actually a vector since we're going from an enum -- which can be
	/// a vector index -- to a string).
	/// @details Used to initialize the const static vector.  Since c++11 initialization of const static singletons
	/// is inherently threadsafe, we don't subsequently have to worry about mutexes for access to this.
	/// @note This is implemented in core/chemical/ResidueProperty_mappings.cc, an automatically-generated file.  DO NOT EDIT
	/// THAT FILE DIRECTLY.  Intead, edit the generating script.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	static utility::vector1< std::string > const * generate_property_to_string_vector();

private:  // Private methods //////////////////////////////////////////////////
	// Initialize data members.
	void init( ResidueType const * residue_type );

	// Copy all data members from <from> to <to>.
	void copy_data( ResidueProperties & to, ResidueProperties const & from );


public:

	/// @brief Static constant data access
	/// @details Get the ResidueProperty enum value from the corresponding string.
	/// This private static class method is defined in ResidueProperty_mappings.cc,
	/// which is auto-generated by the add_ResidueType_enum_files.py script.
	/// Note -- made public because it's a handy function for other code to use (VKM, 22 July 2015).
	static ResidueProperty const & get_property_from_string( std::string const & property );

	/// @brief Get a string from the corresponding ResidueProperty enum value.
	/// @details This private static class method is defined in ResidueProperty_mappings.cc,
	/// which is auto-generated by the add_ResidueType_enum_files.py script.
	/// Note -- made public because it's a handy function for other code to use (VKM, 22 July 2015).
	static std::string const & get_string_from_property( ResidueProperty const property );


public:  // FIXME: Reset this to private; I made it public to work as a temp fix for the following files:
	// - ZnHash.cc
	// - EnzdesBaseProtocol.cc
	// - SecondaryMatchProtocol.cc
	// - SaveAndRetrieveSidechains.cc
	// - pose_manipulation.cc
	// - ModifyVariantTypeMover.cc
	// - ReturnSidechainMover.cc
	// - EnzConstraintParameters.cc
	// - FoldTreeHybridize.cc
	// - protocols/stepwise/modeler/rna/util.cc
	// Get the VariantType enum value from the corresponding string.
	// This private static class method is defined in VariantType_mappings.cc,
	// which is auto-generated by the add_ResidueType_enum_files.py script.
	static VariantType const & get_variant_from_string( std::string const & variant );

public:  // FIXME: Reset this to private; I made it public to work as a temp fix for ResidueTypeSet.cc.
	// Get a string from the corresponding VariantType enum value.
	// This private static class method is defined in VariantType_mappings.cc,
	// which is auto-generated by the add_ResidueType_enum_files.py script.
	static std::string const & get_string_from_variant( VariantType const variant );


private:  // Private data /////////////////////////////////////////////////////
	/// @brief Back pointer to the owning ResidueType
	/// @details (Luki made this a const * during the switch to smart pointers, because get_self_weak_ptr() is not available
	/// within ResidueType's constructors.)
	ResidueType const * residue_type_;  // useful for tracer output

	/// @brief Storage of general properties.
	utility::vector1< bool > general_property_status_;  // indexed by ResidueProperty

	/// @brief The patch operations/variant types that describe this residue.
	utility::vector1< bool > variant_type_status_;  // indexed by VariantType

	/// @brief "Custom" variants as strings are permitted for the enzdes and metalloproteins cases.
	bool has_custom_variant_types_;
	utility::vector1< std::string > custom_variant_types_;

	/// @brief Arbitrary numeric properties with string names.
	std::map<std::string,core::Real> numeric_properties_;

	/// @brief Arbitrary string properties with string names.
	std::map<std::string,std::string> string_properties_;
#ifdef    SERIALIZATION
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< ResidueProperties > & construct );
#endif // SERIALIZATION

};


// Insertion operator (overloaded so that ResidueProperties can be "printed" in PyRosetta).
std::ostream & operator<<( std::ostream & output, ResidueProperties const & object_to_output );

// This allows one to use a for loop with ResidueProperty enum values.
ResidueProperty & operator++( ResidueProperty & property );

// This allows one to use a for loop with VariantType enum values.
VariantType & operator++( VariantType & variant );

}  // namespace chemical
}  // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueProperties )
#endif // SERIALIZATION


#endif  // INCLUDED_core_chemical_ResidueProperties_HH
