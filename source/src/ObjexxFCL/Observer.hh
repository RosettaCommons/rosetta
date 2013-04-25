#ifndef INCLUDED_ObjexxFCL_Observer_hh
#define INCLUDED_ObjexxFCL_Observer_hh


// Observer: Combined Subject + Observer Abstract Base Class
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
#include <ObjexxFCL/Observer.fwd.hh>


namespace ObjexxFCL {


/// @brief Observer: Combined Subject + Observer Abstract Base Class
class Observer
{


protected: // Creation


	/// @brief Default Constructor
	inline
	Observer()
	{}


	/// @brief Copy Constructor
	inline
	Observer( Observer const & )
	{}


public: // Creation


	/// @brief Destructor
	inline
	virtual
	~Observer()
	{}


protected: // Assignment


	/// @brief Copy Assignment
	inline
	Observer &
	operator =( Observer const & )
	{
		return *this;
	}


public: // Subject Inspector


	/// @brief Insert an Observer
	virtual
	void
	insert_observer( Observer & ) const = 0;


	/// @brief Remove an Observer
	virtual
	void
	remove_observer( Observer & ) const = 0;


	/// @brief Has At Least One Observer?
	virtual
	bool
	has_observer() const = 0;


	/// @brief Notify Observers That This Subject Has Changed
	void
	notify() const;


	/// @brief Acyclic After Adding an Observer of This Subject?
	bool
	acyclic( Observer & ) const;


public: // Observer Modifier


	/// @brief Update
	virtual
	void
	update() = 0;


	/// @brief Update for Destruction of a Subject
	virtual
	void
	destructed( Subject const & ) = 0;


}; // Observer


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_Observer_HH
