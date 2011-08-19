#ifndef INCLUDED_ObjexxFCL_DimensionExpression_hh
#define INCLUDED_ObjexxFCL_DimensionExpression_hh


// DimensionExpression: DimensionExpression Interface Class
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
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Observer.fwd.hh>

// C++ Headers
#include <cassert>
#include <iosfwd>


namespace ObjexxFCL {


/// @brief DimensionExpression: DimensionExpression Interface Class
class DimensionExpression
{


protected: // Creation


	/// @brief Default Constructor
	inline
	DimensionExpression()
	{}


	/// @brief Copy Constructor
	inline
	DimensionExpression( DimensionExpression const & )
	{}


public: // Creation


	/// @brief Clone
	virtual
	DimensionExpression *
	clone() const = 0;


	/// @brief Clone with Dimension Substitution
	virtual
	DimensionExpression *
	clone( Dimension const & ) const = 0;


	/// @brief Destructor
	inline
	virtual
	~DimensionExpression()
	{}


public: // Conversion


	/// @brief int Conversion
	inline
	operator int() const
	{
		assert( initialized() );
		return static_cast< int >( value() );
	}


	/// @brief double Conversion
	inline
	operator double() const
	{
		assert( initialized() );
		return value();
	}


private: // Assignment


	/// @brief Copy Assignment
	DimensionExpression &
	operator =( DimensionExpression const & ); // Unimplemented: Dimension handles assignment


public: // Inspector


	/// @brief Initialized?
	virtual
	bool
	initialized() const = 0;


	/// @brief Integer?
	virtual
	bool
	integer() const = 0;


	/// @brief Constant?
	virtual
	bool
	constant() const = 0;


	/// @brief Reference?
	virtual
	bool
	reference() const = 0;


	/// @brief Reducible?
	virtual
	bool
	reducible() const = 0;


	/// @brief Value
	virtual
	double
	operator ()() const = 0;


	/// @brief Value
	virtual
	double
	value() const = 0;


	/// @brief Integer Value
	inline
	virtual
	int
	ivalue() const
	{
		return static_cast< int >( value() );
	}


	/// @brief Integer Value: Zero if Uninitialized
	inline
	virtual
	int
	zvalue() const
	{
		return ( initialized() ? static_cast< int >( value() ) : 0 );
	}


	/// @brief Insert an Observer
	virtual
	void
	insert_observer( Observer & ) const = 0;


	/// @brief Remove an Observer
	virtual
	void
	remove_observer( Observer & ) const = 0;


public: // Modifier


	/// @brief Update for Destruction of a Subject
	virtual
	void
	destructed( Subject const & ) = 0;


}; // DimensionExpression


// Comparison


/// @brief DimensionExpression == DimensionExpression
inline
bool
operator ==( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return ( ( exp1.initialized() ) && ( exp2.initialized() ) && ( exp1.value() == exp2.value() ) );
}


/// @brief DimensionExpression != DimensionExpression
inline
bool
operator !=( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return !( exp1 == exp2 );
}


/// @brief DimensionExpression < DimensionExpression
inline
bool
operator <( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return ( ( exp1.initialized() ) && ( exp2.initialized() ) && ( exp1.value() < exp2.value() ) );
}


/// @brief DimensionExpression <= DimensionExpression
inline
bool
operator <=( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return ( ( exp1.initialized() ) && ( exp2.initialized() ) && ( exp1.value() <= exp2.value() ) );
}


/// @brief DimensionExpression > DimensionExpression
inline
bool
operator >( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return ( ( exp1.initialized() ) && ( exp2.initialized() ) && ( exp1.value() > exp2.value() ) );
}


/// @brief DimensionExpression >= DimensionExpression
inline
bool
operator >=( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return ( ( exp1.initialized() ) && ( exp2.initialized() ) && ( exp1.value() >= exp2.value() ) );
}


/// @brief int == DimensionExpression
inline
bool
operator ==( int const i, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( i == exp.value() ) );
}


/// @brief int != DimensionExpression
inline
bool
operator !=( int const i, DimensionExpression const & exp )
{
	return !( i == exp );
}


/// @brief int < DimensionExpression
inline
bool
operator <( int const i, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( i < exp.value() ) );
}


/// @brief int <= DimensionExpression
inline
bool
operator <=( int const i, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( i <= exp.value() ) );
}


/// @brief int > DimensionExpression
inline
bool
operator >( int const i, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( i > exp.value() ) );
}


/// @brief int >= DimensionExpression
inline
bool
operator >=( int const i, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( i >= exp.value() ) );
}


/// @brief DimensionExpression == int
inline
bool
operator ==( DimensionExpression const & exp, int const i )
{
	return ( ( exp.initialized() ) && ( exp.value() == i ) );
}


/// @brief DimensionExpression != int
inline
bool
operator !=( DimensionExpression const & exp, int const i )
{
	return !( exp == i );
}


/// @brief DimensionExpression < int
inline
bool
operator <( DimensionExpression const & exp, int const i )
{
	return ( ( exp.initialized() ) && ( exp.value() < i ) );
}


/// @brief DimensionExpression <= int
inline
bool
operator <=( DimensionExpression const & exp, int const i )
{
	return ( ( exp.initialized() ) && ( exp.value() <= i ) );
}


/// @brief DimensionExpression > int
inline
bool
operator >( DimensionExpression const & exp, int const i )
{
	return ( ( exp.initialized() ) && ( exp.value() > i ) );
}


/// @brief DimensionExpression >= int
inline
bool
operator >=( DimensionExpression const & exp, int const i )
{
	return ( ( exp.initialized() ) && ( exp.value() >= i ) );
}


/// @brief double == DimensionExpression
inline
bool
operator ==( double const d, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( d == exp.value() ) );
}


/// @brief double != DimensionExpression
inline
bool
operator !=( double const d, DimensionExpression const & exp )
{
	return !( d == exp );
}


/// @brief double < DimensionExpression
inline
bool
operator <( double const d, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( d < exp.value() ) );
}


/// @brief double <= DimensionExpression
inline
bool
operator <=( double const d, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( d <= exp.value() ) );
}


/// @brief double > DimensionExpression
inline
bool
operator >( double const d, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( d > exp.value() ) );
}


/// @brief double >= DimensionExpression
inline
bool
operator >=( double const d, DimensionExpression const & exp )
{
	return ( ( exp.initialized() ) && ( d >= exp.value() ) );
}


/// @brief DimensionExpression == double
inline
bool
operator ==( DimensionExpression const & exp, double const d )
{
	return ( ( exp.initialized() ) && ( exp.value() == d ) );
}


/// @brief DimensionExpression != double
inline
bool
operator !=( DimensionExpression const & exp, double const d )
{
	return !( exp == d );
}


/// @brief DimensionExpression < double
inline
bool
operator <( DimensionExpression const & exp, double const d )
{
	return ( ( exp.initialized() ) && ( exp.value() < d ) );
}


/// @brief DimensionExpression <= double
inline
bool
operator <=( DimensionExpression const & exp, double const d )
{
	return ( ( exp.initialized() ) && ( exp.value() <= d ) );
}


/// @brief DimensionExpression > double
inline
bool
operator >( DimensionExpression const & exp, double const d )
{
	return ( ( exp.initialized() ) && ( exp.value() > d ) );
}


/// @brief DimensionExpression >= double
inline
bool
operator >=( DimensionExpression const & exp, double const d )
{
	return ( ( exp.initialized() ) && ( exp.value() >= d ) );
}


// I/O


/// @brief Stream Output
std::ostream &
operator <<( std::ostream & stream, DimensionExpression const & exp );


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_DimensionExpression_HH
