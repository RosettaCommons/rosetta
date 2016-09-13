// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/SingletonBase.hh
/// @brief  A base class for all singltons using CRTP and managing the complexity of safely
///         initializing singltons in a thread-safe way.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_SingletonBase_HH
#define INCLUDED_utility_SingletonBase_HH

namespace utility {

/// @brief SingletonBase is meant to serve as a base class for singleton classes in Rosetta
/// handling the initialization of the singleton in a thread-safe way.
///
/// @details The derived class must a) implement a private, static function:
/// T * create_singleton_instance()
/// so that the SingletonBase class can invoke this function, and b) declare the
/// SingletonBase class to be a friend, so that it can invoke this function
template < class T >
class SingletonBase
{
public:
	/// @brief public constructor (the derived class must have a private constructor, of course).
	SingletonBase() {}

	/// @brief Safely instantiate a singleton class in a (possibly)
	/// multithreaded context.
	static
	T *
	get_instance() {
		// The C++11 memory model ensures function-scope static variables are initialized in a thread-safe manner
		// http://stackoverflow.com/questions/11711920/how-to-implement-multithread-safe-singleton-in-c11-without-using-mutex
		static T* instance_{ T::create_singleton_instance() };
		return instance_; // Return pointer, to keep current interface
	}

private:
	/// @brief Private, unimplemented copy constructor -- uncopyable.
	SingletonBase( SingletonBase< T > const & );

	/// @brief Private, unimplemented assignment operator -- uncopyable.
	SingletonBase< T > const & operator = ( SingletonBase< T > const & rhs );

private:

};


} // utility

#endif
