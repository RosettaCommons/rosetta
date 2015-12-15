// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/serialization/serialization.hh
/// @brief  Commons serlialization routines
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_serialization_serialization_HH
#define INCLUDED_utility_serialization_serialization_HH

#ifdef SERIALIZATION

#include <utility/vector1.hh>

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

#define EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( T ) \
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




namespace utility {
namespace serialization {

}
}

#endif // SERIALIZATION

#endif // INCLUDED_utility_serialization_serialization_HH
