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
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMediator.hh>


namespace ObjexxFCL {


// Observer: Combined Subject + Observer Abstract Base Class


	/// @brief Notify Observers That This Subject Has Changed
	void
	Observer::notify() const
	{
		if ( has_observer() ) internal::ObserverMediator::notify( *this );
	}


	/// @brief Acyclic After Adding an Observer of This Subject?
	bool
	Observer::acyclic( Observer & observer ) const
	{
		return internal::ObserverMediator::acyclic( *this, observer );
	}


// Observer


} // namespace ObjexxFCL
