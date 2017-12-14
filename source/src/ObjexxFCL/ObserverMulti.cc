// ObserverMulti: Combined Subject + Multi Observer Abstract Base Class
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/SetWrapper.hh>

// C++ Headers
#include <cassert>


namespace ObjexxFCL {


// ObserverMulti: Combined Subject + Multi Observer Abstract Base Class


/// @brief Destructor
ObserverMulti::~ObserverMulti()
{
	notify_destructed();
	delete observers_p_;
}


/// @brief Insert an Observer
void
ObserverMulti::insert_observer( Observer & observer ) const
{
	assert( this != &observer );
	assert( acyclic( observer ) );
	if ( ! observers_p_ ) observers_p_ = new Observers;
	(*observers_p_)().insert( &observer );
}


/// @brief Remove an Observer
void
ObserverMulti::do_remove_observer( Observer & observer ) const
{
	assert( observers_p_ );
	(*observers_p_)().erase( &observer );
}


/// @brief Has At Least One Observer?
bool
ObserverMulti::do_has_observer() const
{
	assert( observers_p_ );
	return !(*observers_p_)().empty();
}


/// @brief Notify Observers That This Subject is Being Destructed
void
ObserverMulti::do_notify_destructed() const
{
	assert( observers_p_ );
	for ( auto io : (*observers_p_)() ) {
		assert( io );
		io->destructed( *this );
	}
}


// ObserverMulti


} // namespace ObjexxFCL
