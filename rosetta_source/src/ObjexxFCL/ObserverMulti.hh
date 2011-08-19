#ifndef INCLUDED_ObjexxFCL_ObserverMulti_hh
#define INCLUDED_ObjexxFCL_ObserverMulti_hh


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
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>

// C++ Headers
#include <cassert>


namespace ObjexxFCL {


/// @brief ObserverMulti: Combined Subject + Multi Observer Abstract Base Class
class ObserverMulti :
	public Observer
{


public: // Types


	typedef  SetWrapper< Observer * >  Observers;


protected: // Creation


	/// @brief Default Constructor
	inline
	ObserverMulti() :
		observers_p_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	ObserverMulti( ObserverMulti const & ) :
		Observer(),
		observers_p_( 0 )
	{}


public: // Creation


	/// @brief Destructor
	virtual
	~ObserverMulti();


protected: // Assignment


	/// @brief Copy Assignment
	inline
	ObserverMulti &
	operator =( ObserverMulti const & )
	{
		return *this;
	}


public: // Subject Inspector


	/// @brief Insert an Observer
	void
	insert_observer( Observer & observer ) const;


	/// @brief Remove an Observer
	inline
	void
	remove_observer( Observer & observer ) const
	{
		if ( observers_p_ ) do_remove_observer( observer );
	}


	/// @brief Has At Least One Observer?
	inline
	bool
	has_observer() const
	{
		return ( observers_p_ ? do_has_observer() : false );
	}


	/// @brief Observers Pointer
	inline
	Observers const *
	observers_p() const
	{
		return observers_p_;
	}


	/// @brief Observers
	inline
	Observers const &
	observers() const
	{
		assert( observers_p_ );
		return *observers_p_;
	}


	/// @brief Notify Observers That This Subject is Being Destructed
	inline
	void
	notify_destructed() const
	{
		if ( observers_p_ ) do_notify_destructed();
	}


private: // Functions


	/// @brief Remove an Observer
	void
	do_remove_observer( Observer & observer ) const;


	/// @brief Has At Least One Observer?
	bool
	do_has_observer() const;


	/// @brief Notify Observers That This Subject is Being Destructed
	void
	do_notify_destructed() const;


private: // Data


	/// @brief Observers of this Subject
	mutable Observers * observers_p_;


}; // ObserverMulti


// Types
typedef  ObserverMulti  SubjectMulti;


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_ObserverMulti_HH
