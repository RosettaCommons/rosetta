// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/serialization/serialization.hh
/// @brief  Commons serlialization routines
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_serialization_serialization_HH
#define INCLUDED_utility_serialization_serialization_HH

#ifdef SERIALIZATION

//#include <utility/vector1.hh>

// Archives
#include <cereal/archives/binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>

/////// Macros (*sigh*) for working with serialization.  A necessary evil /////////

// This macro should be used in the .srlz.cc file for classes that implement
// the serialize method
#define SERIALIZABLE(T) \
template void T::serialize(cereal::BinaryOutputArchive&); \
template void T::serialize(cereal::BinaryInputArchive&); \
template void T::serialize(cereal::XMLOutputArchive&); \
template void T::serialize(cereal::XMLInputArchive&); \
template void T::serialize(cereal::JSONOutputArchive&); \
template void T::serialize(cereal::JSONInputArchive&)

// This macro should be used in the .srlz.cc file for classes that implement
// the split save and load methods.
#define SAVE_AND_LOAD_SERIALIZABLE(T) \
template void T::save( typename cereal::BinaryOutputArchive &) const; \
template void T::save( typename cereal::XMLOutputArchive &) const; \
template void T::save( typename cereal::JSONOutputArchive &) const; \
template void T::load( typename cereal::BinaryInputArchive &); \
template void T::load( typename cereal::XMLInputArchive &); \
template void T::load( typename cereal::JSONInputArchive &)


// PyRosetta macro to forward declare template instantiation for the split save and load methods.
#define EXTERN_SAVE_AND_LOAD_SERIALIZABLE(T) \
extern template void T::save( typename cereal::BinaryOutputArchive &) const; \
extern template void T::save( typename cereal::XMLOutputArchive &) const; \
extern template void T::save( typename cereal::JSONOutputArchive &) const; \
extern template void T::load( typename cereal::BinaryInputArchive &); \
extern template void T::load( typename cereal::XMLInputArchive &); \
extern template void T::load( typename cereal::JSONInputArchive &)


#define EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( T )						\
template void save( typename cereal::BinaryOutputArchive &, T const & ); \
template void save( typename cereal::XMLOutputArchive &, T const & ); \
template void save( typename cereal::JSONOutputArchive &, T const & ); \
template void load( typename cereal::BinaryInputArchive &, T & ); \
template void load( typename cereal::XMLInputArchive &, T & ); \
template void load( typename cereal::JSONInputArchive &, T & )

#define SAVE_AND_LOAD_AND_CONSTRUCT_SERIALIZABLE( T ) \
template void T::save( typename cereal::BinaryOutputArchive &) const; \
template void T::save( typename cereal::XMLOutputArchive &) const; \
template void T::save( typename cereal::JSONOutputArchive &) const; \
template void T::load_and_construct( typename cereal::BinaryInputArchive &, cereal::construct< T > & construct ); \
template void T::load_and_construct( typename cereal::XMLInputArchive &, cereal::construct< T > & construct ); \
template void T::load_and_construct( typename cereal::JSONInputArchive &, cereal::construct< T > & construct )

// This macro is to force the instantiation of a templated funtion for each
// of the active output archive types -- in the future, different compiler
// options will cause the instantiation of different archive types.  This
// function should return a type named "return_type", be named "func", and
// take two parameters: an archive, and a second parameter, whose type should
// be "param."  Perhaps alternative flavors of this macro will be necessary.
// The idea would be for the given function to somehow serialize the second
// parameter into the archive.
#define INSTANTIATE_FOR_OUTPUT_ARCHIVES( return_type, func, param )	\
template return_type func( typename cereal::BinaryOutputArchive &, param ); \
template return_type func( typename cereal::XMLOutputArchive &, param ); \
template return_type func( typename cereal::JSONOutputArchive &, param )

// This macro is to force the instantiation of a templated funtion for each
// of the active output archive types -- in the future, different compiler
// options will cause the instantiation of different archive types.  This
// function should return a type named "return_type", be named "func", and
// take two parameters: an archive, and a second parameter, whose type should
// be "param."  Perhaps alternative flavors of this macro will be necessary.
// The idea would be for the given function to somehow deserialize an object into
// the second parameter from the archive.
#define INSTANTIATE_FOR_INPUT_ARCHIVES( return_type, func, param )	\
template return_type func( typename cereal::BinaryInputArchive &, param ); \
template return_type func( typename cereal::XMLInputArchive &, param ); \
template return_type func( typename cereal::JSONInputArchive &, param  )

// This macro is for specializing handling of COPs
// for example, when you want to point some to a Global (non-serialized) resource.
// Place this macro in the header file of the type, so that the redefinitions are
// found by every class which uses/serializes the class.
// - type_name is a *fully namespace qualified* name of the base type to specialize
// - save_function is the (fully namespace qualified) name of the function which takes a
//   templated Archive by reference and the COP by *value* and saves to the archive.
// - load_function is the (fully namespace qualified) name of the function which takes a
//   templated Archive by reference and the COP by *reference* and loads from the archive.
// IMPORTANT!! - If in your save_function/load_function you archive the object (as opposed
// to a parent/subclass OP) *you* are responsible for handling deduplication of pointers.
// This can easily be accomplished with the arc.registerSharedPointer() facility.
// (See core/chemical/ResidueType.srlz.cc for example.)
template<typename T> struct LOAD_DEPENDENT_FALSE: std::false_type {};

#define SPECIAL_COP_SERIALIZATION_HANDLING( type_name, save_function, load_function ) \
namespace cereal { \
template < class Archive > \
void save( Archive & arc, utility::pointer::shared_ptr< type_name const > const & ptr ) \
{ save_function( arc, ptr ); } \
template < class Archive > \
void save( Archive & arc, utility::pointer::shared_ptr< type_name > const & ptr ) \
{ save_function( arc, ptr ); } \
template < class Archive > \
void load( Archive & arc, utility::pointer::shared_ptr< type_name const > & ptr ) \
{ load_function( arc, ptr ); } \
template < class Archive > \
void load( Archive &, utility::pointer::shared_ptr< type_name > & ) \
{ static_assert( LOAD_DEPENDENT_FALSE<Archive>::value , "You cannot load into a " #type_name "OP - must use a " #type_name "COP" ); } \
}

namespace utility {
namespace serialization {

}
}

#endif // SERIALIZATION

#endif // INCLUDED_utility_serialization_serialization_HH
