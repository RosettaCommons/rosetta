#ifndef INCLUDED_ObjexxFCL_DimensionExpressions_hh
#define INCLUDED_ObjexxFCL_DimensionExpressions_hh


// DimensionExpressions: DimensionExpression Headers for Sources that Use Dimension Expressions
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
#include <ObjexxFCL/DimensionExpressionCon.hh>
#include <ObjexxFCL/DimensionExpressionRef.hh>
#include <ObjexxFCL/DimensionExpressionSum.hh>
#include <ObjexxFCL/DimensionExpressionSub.hh>
#include <ObjexxFCL/DimensionExpressionMul.hh>
#include <ObjexxFCL/DimensionExpressionDiv.hh>
#include <ObjexxFCL/DimensionExpressionMin.hh>
#include <ObjexxFCL/DimensionExpressionMax.hh>
#include <ObjexxFCL/DimensionExpressionPow.hh>
#include <ObjexxFCL/DimensionExpressionSquare.hh>
#include <ObjexxFCL/DimensionExpressionCube.hh>


namespace ObjexxFCL {


// Dimension Math


/// @brief +Dimension
inline
DimensionExpressionRef
operator +( Dimension const & dim )
{
	return DimensionExpressionRef( dim );
}


/// @brief -Dimension
inline
DimensionExpressionMul
operator -( Dimension const & dim )
{
	return DimensionExpressionMul( new DimensionExpressionCon( -1 ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension + Dimension
inline
DimensionExpressionSum
operator +( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionSum( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief Dimension + int
inline
DimensionExpressionSum
operator +( Dimension const & dim, int const value )
{
	return DimensionExpressionSum( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief int + Dimension
inline
DimensionExpressionSum
operator +( int const value, Dimension const & dim )
{
	return DimensionExpressionSum( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension + double
inline
DimensionExpressionSum
operator +( Dimension const & dim, double const value )
{
	return DimensionExpressionSum( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief double + Dimension
inline
DimensionExpressionSum
operator +( double const value, Dimension const & dim )
{
	return DimensionExpressionSum( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension - Dimension
inline
DimensionExpressionSub
operator -( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionSub( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief Dimension - int
inline
DimensionExpressionSub
operator -( Dimension const & dim, int const value )
{
	return DimensionExpressionSub( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief int - Dimension
inline
DimensionExpressionSub
operator -( int const value, Dimension const & dim )
{
	return DimensionExpressionSub( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension - double
inline
DimensionExpressionSub
operator -( Dimension const & dim, double const value )
{
	return DimensionExpressionSub( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief double - Dimension
inline
DimensionExpressionSub
operator -( double const value, Dimension const & dim )
{
	return DimensionExpressionSub( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension * Dimension
inline
DimensionExpressionMul
operator *( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionMul( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief Dimension * int
inline
DimensionExpressionMul
operator *( Dimension const & dim, int const value )
{
	return DimensionExpressionMul( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief int * Dimension
inline
DimensionExpressionMul
operator *( int const value, Dimension const & dim )
{
	return DimensionExpressionMul( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension * double
inline
DimensionExpressionMul
operator *( Dimension const & dim, double const value )
{
	return DimensionExpressionMul( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief double * Dimension
inline
DimensionExpressionMul
operator *( double const value, Dimension const & dim )
{
	return DimensionExpressionMul( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension / Dimension
inline
DimensionExpressionDiv
operator /( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionDiv( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief Dimension / int
inline
DimensionExpressionDiv
operator /( Dimension const & dim, int const value )
{
	return DimensionExpressionDiv( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief int / Dimension
inline
DimensionExpressionDiv
operator /( int const value, Dimension const & dim )
{
	return DimensionExpressionDiv( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension / double
inline
DimensionExpressionDiv
operator /( Dimension const & dim, double const value )
{
	return DimensionExpressionDiv( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief double / Dimension
inline
DimensionExpressionDiv
operator /( double const value, Dimension const & dim )
{
	return DimensionExpressionDiv( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


// DimensionExpression Math


/// @brief +DimensionExpression
inline
DimensionExpressionMul
operator +( DimensionExpression const & exp )
{
	return DimensionExpressionMul( new DimensionExpressionCon( 1 ), exp.clone() );
}


/// @brief -DimensionExpression
inline
DimensionExpressionMul
operator -( DimensionExpression const & exp )
{
	return DimensionExpressionMul( new DimensionExpressionCon( -1 ), exp.clone() );
}


/// @brief DimensionExpression + DimensionExpression
inline
DimensionExpressionSum
operator +( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionSum( exp1.clone(), exp2.clone() );
}


/// @brief DimensionExpression + Dimension
inline
DimensionExpressionSum
operator +( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionSum( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension + DimensionExpression
inline
DimensionExpressionSum
operator +( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionSum( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief DimensionExpression + int
inline
DimensionExpressionSum
operator +( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionSum( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief int + DimensionExpression
inline
DimensionExpressionSum
operator +( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionSum( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression + double
inline
DimensionExpressionSum
operator +( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionSum( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief double + DimensionExpression
inline
DimensionExpressionSum
operator +( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionSum( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression - DimensionExpression
inline
DimensionExpressionSub
operator -( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionSub( exp1.clone(), exp2.clone() );
}


/// @brief DimensionExpression - Dimension
inline
DimensionExpressionSub
operator -( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionSub( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension - DimensionExpression
inline
DimensionExpressionSub
operator -( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionSub( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief DimensionExpression - int
inline
DimensionExpressionSub
operator -( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionSub( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief int - DimensionExpression
inline
DimensionExpressionSub
operator -( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionSub( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression - double
inline
DimensionExpressionSub
operator -( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionSub( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief double - DimensionExpression
inline
DimensionExpressionSub
operator -( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionSub( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression * DimensionExpression
inline
DimensionExpressionMul
operator *( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionMul( exp1.clone(), exp2.clone() );
}


/// @brief DimensionExpression * Dimension
inline
DimensionExpressionMul
operator *( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionMul( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension * DimensionExpression
inline
DimensionExpressionMul
operator *( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionMul( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief DimensionExpression * int
inline
DimensionExpressionMul
operator *( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionMul( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief int * DimensionExpression
inline
DimensionExpressionMul
operator *( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionMul( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression * double
inline
DimensionExpressionMul
operator *( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionMul( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief double * DimensionExpression
inline
DimensionExpressionMul
operator *( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionMul( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression / DimensionExpression
inline
DimensionExpressionDiv
operator /( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionDiv( exp1.clone(), exp2.clone() );
}


/// @brief DimensionExpression / Dimension
inline
DimensionExpressionDiv
operator /( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionDiv( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief Dimension / DimensionExpression
inline
DimensionExpressionDiv
operator /( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionDiv( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief DimensionExpression / int
inline
DimensionExpressionDiv
operator /( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionDiv( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief int / DimensionExpression
inline
DimensionExpressionDiv
operator /( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionDiv( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief DimensionExpression / double
inline
DimensionExpressionDiv
operator /( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionDiv( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief double / DimensionExpression
inline
DimensionExpressionDiv
operator /( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionDiv( new DimensionExpressionCon( value ), exp.clone() );
}


// Min


/// @brief min( Dimension, Dimension )
inline
DimensionExpressionMin
min( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionMin( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief min( Dimension, int )
inline
DimensionExpressionMin
min( Dimension const & dim, int const value )
{
	return DimensionExpressionMin( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief min( int, Dimension )
inline
DimensionExpressionMin
min( int const value, Dimension const & dim )
{
	return DimensionExpressionMin( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief min( Dimension, double )
inline
DimensionExpressionMin
min( Dimension const & dim, double const value )
{
	return DimensionExpressionMin( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief min( double, Dimension )
inline
DimensionExpressionMin
min( double const value, Dimension const & dim )
{
	return DimensionExpressionMin( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief min( Dimension, DimensionExpression )
inline
DimensionExpressionMin
min( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionMin( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief min( DimensionExpression, Dimension )
inline
DimensionExpressionMin
min( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionMin( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief min( DimensionExpression, DimensionExpression )
inline
DimensionExpressionMin
min( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionMin( exp1.clone(), exp2.clone() );
}


/// @brief min( DimensionExpression, int )
inline
DimensionExpressionMin
min( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionMin( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief min( int, DimensionExpression )
inline
DimensionExpressionMin
min( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionMin( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief min( DimensionExpression, double )
inline
DimensionExpressionMin
min( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionMin( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief min( double, DimensionExpression )
inline
DimensionExpressionMin
min( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionMin( new DimensionExpressionCon( value ), exp.clone() );
}


// Max


/// @brief max( Dimension, Dimension )
inline
DimensionExpressionMax
max( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionMax( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief max( Dimension, int )
inline
DimensionExpressionMax
max( Dimension const & dim, int const value )
{
	return DimensionExpressionMax( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief max( int, Dimension )
inline
DimensionExpressionMax
max( int const value, Dimension const & dim )
{
	return DimensionExpressionMax( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief max( Dimension, double )
inline
DimensionExpressionMax
max( Dimension const & dim, double const value )
{
	return DimensionExpressionMax( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief max( double, Dimension )
inline
DimensionExpressionMax
max( double const value, Dimension const & dim )
{
	return DimensionExpressionMax( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief max( Dimension, DimensionExpression )
inline
DimensionExpressionMax
max( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionMax( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief max( DimensionExpression, Dimension )
inline
DimensionExpressionMax
max( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionMax( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief max( DimensionExpression, DimensionExpression )
inline
DimensionExpressionMax
max( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionMax( exp1.clone(), exp2.clone() );
}


/// @brief max( DimensionExpression, int )
inline
DimensionExpressionMax
max( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionMax( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief max( int, DimensionExpression )
inline
DimensionExpressionMax
max( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionMax( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief max( DimensionExpression, double )
inline
DimensionExpressionMax
max( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionMax( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief max( double, DimensionExpression )
inline
DimensionExpressionMax
max( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionMax( new DimensionExpressionCon( value ), exp.clone() );
}


// Pow


/// @brief pow( Dimension, Dimension )
inline
DimensionExpressionPow
pow( Dimension const & dim1, Dimension const & dim2 )
{
	return DimensionExpressionPow( new DimensionExpressionRef( dim1 ), new DimensionExpressionRef( dim2 ) );
}


/// @brief pow( Dimension, int )
inline
DimensionExpressionPow
pow( Dimension const & dim, int const value )
{
	return DimensionExpressionPow( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief pow( int, Dimension )
inline
DimensionExpressionPow
pow( int const value, Dimension const & dim )
{
	return DimensionExpressionPow( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief pow( Dimension, double )
inline
DimensionExpressionPow
pow( Dimension const & dim, double const value )
{
	return DimensionExpressionPow( new DimensionExpressionRef( dim ), new DimensionExpressionCon( value ) );
}


/// @brief pow( double, Dimension )
inline
DimensionExpressionPow
pow( double const value, Dimension const & dim )
{
	return DimensionExpressionPow( new DimensionExpressionCon( value ), new DimensionExpressionRef( dim ) );
}


/// @brief pow( Dimension, DimensionExpression )
inline
DimensionExpressionPow
pow( Dimension const & dim, DimensionExpression const & exp )
{
	return DimensionExpressionPow( new DimensionExpressionRef( dim ), exp.clone() );
}


/// @brief pow( DimensionExpression, Dimension )
inline
DimensionExpressionPow
pow( DimensionExpression const & exp, Dimension const & dim )
{
	return DimensionExpressionPow( exp.clone(), new DimensionExpressionRef( dim ) );
}


/// @brief pow( DimensionExpression, DimensionExpression )
inline
DimensionExpressionPow
pow( DimensionExpression const & exp1, DimensionExpression const & exp2 )
{
	return DimensionExpressionPow( exp1.clone(), exp2.clone() );
}


/// @brief pow( DimensionExpression, int )
inline
DimensionExpressionPow
pow( DimensionExpression const & exp, int const value )
{
	return DimensionExpressionPow( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief pow( int, DimensionExpression )
inline
DimensionExpressionPow
pow( int const value, DimensionExpression const & exp )
{
	return DimensionExpressionPow( new DimensionExpressionCon( value ), exp.clone() );
}


/// @brief pow( DimensionExpression, double )
inline
DimensionExpressionPow
pow( DimensionExpression const & exp, double const value )
{
	return DimensionExpressionPow( exp.clone(), new DimensionExpressionCon( value ) );
}


/// @brief pow( double, DimensionExpression )
inline
DimensionExpressionPow
pow( double const value, DimensionExpression const & exp )
{
	return DimensionExpressionPow( new DimensionExpressionCon( value ), exp.clone() );
}


// Square


/// @brief square( Dimension )
inline
DimensionExpressionSquare
square( Dimension const & dim )
{
	return DimensionExpressionSquare( new DimensionExpressionRef( dim ) );
}


/// @brief square( DimensionExpression )
inline
DimensionExpressionSquare
square( DimensionExpression const & exp )
{
	return DimensionExpressionSquare( exp.clone() );
}


/// @brief square( int )
inline
DimensionExpressionCon
square( int const value )
{
	return DimensionExpressionCon( value * value );
}


/// @brief square( double )
inline
DimensionExpressionCon
square( double const value )
{
	return DimensionExpressionCon( value * value );
}


// Cube


/// @brief cube( Dimension )
inline
DimensionExpressionCube
cube( Dimension const & dim )
{
	return DimensionExpressionCube( new DimensionExpressionRef( dim ) );
}


/// @brief cube( DimensionExpression )
inline
DimensionExpressionCube
cube( DimensionExpression const & exp )
{
	return DimensionExpressionCube( exp.clone() );
}


/// @brief cube( int )
inline
DimensionExpressionCon
cube( int const value )
{
	return DimensionExpressionCon( value * value * value );
}


/// @brief cube( double )
inline
DimensionExpressionCon
cube( double const value )
{
	return DimensionExpressionCon( value * value * value );
}


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_DimensionExpressions_HH
