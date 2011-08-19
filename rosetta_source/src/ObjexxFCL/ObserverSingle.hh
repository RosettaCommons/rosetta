#ifndef INCLUDED_ObjexxFCL_ObserverSingle_hh
#define INCLUDED_ObjexxFCL_ObserverSingle_hh


// ObserverSingle: Combined Subject + Single Observer Abstract Base Class
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
#include <ObjexxFCL/Observer.hh>

// C++ Headers
#include <cassert>


namespace ObjexxFCL {


/// @brief ObserverSingle: Combined Subject + Single Observer Abstract Base Class
class ObserverSingle :
	public Observer
{


protected: // Creation


	/// @brief Default Constructor
	inline
	ObserverSingle() :
		observer_p_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	ObserverSingle( ObserverSingle const & ) :
		Observer(),
		observer_p_( 0 )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~ObserverSingle()
	{
		notify_destructed();
	}


protected: // Assignment


	/// @brief Copy Assignment
	inline
	ObserverSingle &
	operator =( ObserverSingle const & )
	{
		return *this;
	}


public: // Subject Inspector


	/// @brief Insert an Observer
	inline
	void
	insert_observer( Observer & observer ) const
	{
		assert( acyclic( observer ) );
		assert( ( ! observer_p_ ) || ( observer_p_ == &observer ) );
		observer_p_ = &observer;
	}


	/// @brief Remove an Observer
	inline
	void
	remove_observer( Observer & observer ) const
	{
		if ( observer_p_ == &observer ) observer_p_ = 0;
	}


	/// @brief Has At Least One Observer?
	inline
	bool
	has_observer() const
	{
		return observer_p_;
	}


	/// @brief Notify Observers That This Subject is Being Destructed
	inline
	void
	notify_destructed() const
	{
		if ( observer_p_ ) observer_p_->destructed( *this );
	}


	/// @brief Observer
	inline
	Observer *
	observer_p() const
	{
		return observer_p_;
	}


private: // Data


	/// @brief Observer of this Subject (non-owning pointer)
	mutable Observer * observer_p_;


}; // ObserverSingle


// Types
typedef  ObserverSingle  SubjectSingle;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ObserverSingle_HH
