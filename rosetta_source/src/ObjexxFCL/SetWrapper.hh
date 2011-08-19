#ifndef INCLUDED_ObjexxFCL_SetWrapper_hh
#define INCLUDED_ObjexxFCL_SetWrapper_hh


// SetWrapper: Insulating Wrapper of std::set That Can Be Forward Declared
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
#include <ObjexxFCL/SetWrapper.fwd.hh>

// C++ Headers
#include <set>


namespace ObjexxFCL {


/// @brief SetWrapper: Insulating Wrapper of std::set That Can Be Forward Declared
/// @note For objects that manage their own memory not pointers to owned objects
template< typename T >
class SetWrapper
{


public: // Types


	typedef  std::set< T >  Container;

	// STL style
	typedef  T  value_type;
	typedef  typename Container::iterator  iterator;
	typedef  typename Container::const_iterator  const_iterator;

	// C++ style
	typedef  T  Value;
	typedef  typename Container::iterator  Iterator;
	typedef  typename Container::const_iterator  ConstIterator;


public: // Creation


	/// @brief Default Constructor
	inline
	SetWrapper()
	{}


	/// @brief Destructor
	inline
	~SetWrapper()
	{}


public: // Inspector


	/// @brief set Accessor
	inline
	Container const &
	operator ()() const
	{
		return container_;
	}


public: // Modifier


	/// @brief set Accessor
	inline
	Container &
	operator ()()
	{
		return container_;
	}


private: // Data


	/// @brief std::set being wrapped
	Container container_;


}; // SetWrapper


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_SetWrapper_HH
