// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/xyzMatrix.hh
/// @brief  Fast 3x3 matrix
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class


#ifndef INCLUDED_numeric_xyzMatrix_hh
#define INCLUDED_numeric_xyzMatrix_hh


// Unit headers
#include <numeric/xyzMatrix.fwd.hh>

// Package headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/internal/RowVectors.hh>

// C++ headers
#include <utility/assert.hh>


namespace numeric {


/// @brief xyzMatrix: Fast 3x3 xyz matrix template

template< typename T >
class xyzMatrix
{


private: // Friends


	template< typename > friend class xyzMatrix;

	// Friend Functions (for speed of non-inlining debug builds)
	friend xyzVector< T > operator *<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > product<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > & inplace_product<>( xyzMatrix< T > const & m, xyzVector< T > & v );
	friend xyzVector< T > transpose_product<>( xyzMatrix< T > const & m, xyzVector< T > const & v );
	friend xyzVector< T > & inplace_transpose_product<>( xyzMatrix< T > const & m, xyzVector< T > & v );
	friend xyzMatrix< T > outer_product<>( xyzVector< T > const & a, xyzVector< T > const & b );
	friend xyzMatrix< T > projection_matrix<>( xyzVector< T > const & v );
	friend xyzMatrix< T > rotation_matrix<>( xyzVector< T > const & axis, T const & theta );
	friend xyzVector< T > rotation_axis<>( xyzMatrix< T > const & R, T & theta );
	friend xyzVector< T > eigenvalue_jacobi<>( xyzMatrix< T > const & a, T const & tol );
	friend xyzVector< T > eigenvector_jacobi<>( xyzMatrix< T > const & a, T const & tol, xyzMatrix< T > & J );
	friend xyzMatrix< T > inverse<>( xyzMatrix< T > const & a );


public: // Types


	// Project style
	typedef  T          Value;
	typedef  T &        Reference;
	typedef  T const &  ConstReference;
	typedef  T *        Pointer;
	typedef  T const *  ConstPointer;
	typedef  xyzVector< T >  Vector;

	// STL/boost style
	typedef  T          value_type;
	typedef  T &        reference;
	typedef  T const &  const_reference;
	typedef  T *        pointer;
	typedef  T const *  const_pointer;


public: // Creation


	/// @brief Default constructor
	/// @note  Values are uninitialized for efficiency
	inline
	xyzMatrix()
	{}


	/// @brief Copy constructor
	inline
	xyzMatrix( xyzMatrix const & m ) :
		xx_( m.xx_ ), xy_( m.xy_ ), xz_( m.xz_ ),
		yx_( m.yx_ ), yy_( m.yy_ ), yz_( m.yz_ ),
		zx_( m.zx_ ), zy_( m.zy_ ), zz_( m.zz_ )
	{}


	/// @brief Copy constructor
	template< typename U >
	inline
	xyzMatrix( xyzMatrix< U > const & m ) :
		xx_( m.xx_ ), xy_( m.xy_ ), xz_( m.xz_ ),
		yx_( m.yx_ ), yy_( m.yy_ ), yz_( m.yz_ ),
		zx_( m.zx_ ), zy_( m.zy_ ), zz_( m.zz_ )
	{}


	/// @brief Uniform value constructor
	inline
	explicit
	xyzMatrix( Value const & t ) :
		xx_( t ), xy_( t ), xz_( t ),
		yx_( t ), yy_( t ), yz_( t ),
		zx_( t ), zy_( t ), zz_( t )
	{}


	/// @brief Destructor
	inline
	~xyzMatrix()
	{}


public: // Creation: column values


	/// @brief Column-ordered value named constructor
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::cols( xx_a, yx_a, zx_a, xy_a, yy_a, zy_a, xz_a, yz_a, zz_a ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::cols( xx_a, yx_a, zx_a, xy_a, yy_a, zy_a, xz_a, yz_a, zz_a )
	inline
	static
	xyzMatrix
	cols(
		Value const & xx_a, Value const & yx_a, Value const & zx_a,
		Value const & xy_a, Value const & yy_a, Value const & zy_a,
		Value const & xz_a, Value const & yz_a, Value const & zz_a
	)
	{
		return xyzMatrix(
			xx_a, xy_a, xz_a,
			yx_a, yy_a, yz_a,
			zx_a, zy_a, zz_a
		);
	}


public: // Creation: column pointers


	/// @brief Pointer to contiguous column-ordered values constructor
	/// @see   cols
	/// @note  The cols named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage: xyzMatrix m( xyzMatrix::cols( cols_pointer ) )
	/// @note  Assignment usage: m = xyzMatrix::cols( cols_pointer )
	template< typename U >
	inline
	xyzMatrix( ColsPointer< U > const & c ) :
		xx_( *(c.p_) ),   xy_( *(c.p_+3) ), xz_( *(c.p_+6) ),
		yx_( *(c.p_+1) ), yy_( *(c.p_+4) ), yz_( *(c.p_+7) ),
		zx_( *(c.p_+2) ), zy_( *(c.p_+5) ), zz_( *(c.p_+8) )
	{}


	/// @brief Pointer to contiguous column-ordered values named constructor/assignment
	/// @warning No way to check that argument points to nine contiguous values
	/// @note  Constructor usage: xyzMatrix m( xyzMatrix::cols( cols_pointer ) )
	/// @note  Assignment usage: m = xyzMatrix::cols( cols_pointer )
	template< typename U >
	inline
	static
	ColsPointer< U >
	cols( U const * cp )
	{
		return ColsPointer< U >( cp );
	}


	/// @brief Pointers to contiguous columns constructor
	/// @see   cols
	/// @note  The cols named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::cols( x_col_pointer, y_col_pointer, z_col_pointer ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::cols( x_col_pointer, y_col_pointer, z_col_pointer )
	template< typename U >
	inline
	xyzMatrix( ColPointers< U > const & c ) :
		xx_( *c.xp_ ),     xy_( *c.yp_ ),     xz_( *c.zp_ ),
		yx_( *(c.xp_+1) ), yy_( *(c.yp_+1) ), yz_( *(c.zp_+1) ),
		zx_( *(c.xp_+2) ), zy_( *(c.yp_+2) ), zz_( *(c.zp_+2) )
	{}


	/// @brief Pointers to contiguous columns named constructor/assignment
	/// @warning No way to check that arguments each point to three contiguous values
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::cols( x_col_pointer, y_col_pointer, z_col_pointer ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::cols( x_col_pointer, y_col_pointer, z_col_pointer )
	template< typename U >
	inline
	static
	ColPointers< U >
	cols(
		U const * xp,
		U const * yp,
		U const * zp
	)
	{
		return ColPointers< U >( xp, yp, zp );
	}


	/// @brief Pointers to contiguous columns named constructor
	/// @warning No way to check that arguments each point to three contiguous values
	/// @note  Can be faster than cols() construction due to return value optimization
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::cols_constructor( x_col_pointer, y_col_pointer, z_col_pointer ) )
	template< typename U >
	inline
	static
	xyzMatrix
	cols_constructor(
		U const * xp,
		U const * yp,
		U const * zp
	)
	{
		return xyzMatrix(
			*xp,     *yp,     *zp,
			*(xp+1), *(yp+1), *(zp+1),
			*(xp+2), *(yp+2), *(zp+2)
		);
	}


public: // Creation: column xyzVectors


	/// @brief Columns constructor
	/// @see   cols
	/// @note  The cols named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::cols( x_col_vector, y_col_vector, z_col_vector ) )
	template< typename U >
	inline
	xyzMatrix( ColVectors< U > const & c ) :
		xx_( c.x_.x_ ), xy_( c.y_.x_ ), xz_( c.z_.x_ ),
		yx_( c.x_.y_ ), yy_( c.y_.y_ ), yz_( c.z_.y_ ),
		zx_( c.x_.z_ ), zy_( c.y_.z_ ), zz_( c.z_.z_ )
	{}


	/// @brief Column vectors named constructor/assignment
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::cols( x_col_vector, y_col_vector, z_col_vector ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::cols( x_col_vector, y_col_vector, z_col_vector )
	template< typename U >
	inline
	static
	ColVectors< U >
	cols(
		xyzVector< U > const & col_x,
		xyzVector< U > const & col_y,
		xyzVector< U > const & col_z
	)
	{
		return ColVectors< U >( col_x, col_y, col_z );
	}


	/// @brief xyzVector columns named constructor
	/// @note  Can be faster than cols() construction due to return value optimization
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::cols_constructor( x_col_vector, y_col_vector, z_col_vector ) )
	template< typename U >
	inline
	static
	xyzMatrix
	cols_constructor(
		xyzVector< U > const & col_x,
		xyzVector< U > const & col_y,
		xyzVector< U > const & col_z
	)
	{
		return xyzMatrix(
			col_x.x_, col_y.x_, col_z.x_,
			col_x.y_, col_y.y_, col_z.y_,
			col_x.z_, col_y.z_, col_z.z_
		);
	}


public: // Creation: row values


	/// @brief Row-ordered value named constructor
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::rows( xx_a, xy_a, xz_a, yx_a, yy_a, yz_a, zx_a, zy_a, zz_a ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::rows( xx_a, xy_a, xz_a, yx_a, yy_a, yz_a, zx_a, zy_a, zz_a )
	inline
	static
	xyzMatrix
	rows(
		Value const & xx_a, Value const & xy_a, Value const & xz_a,
		Value const & yx_a, Value const & yy_a, Value const & yz_a,
		Value const & zx_a, Value const & zy_a, Value const & zz_a
	)
	{
		return xyzMatrix(
			xx_a, xy_a, xz_a,
			yx_a, yy_a, yz_a,
			zx_a, zy_a, zz_a
		);
	}


public: // Creation: row pointers


	/// @brief Pointer to contiguous row-ordered values constructor
	/// @see   rows
	/// @note  The rows named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage: xyzMatrix m( xyzMatrix::rows( rows_pointer ) )
	/// @note  Assignment usage: m = xyzMatrix::rows( rows_pointer )
	template< typename U >
	inline
	xyzMatrix( RowsPointer< U > const & r ) :
		xx_( *(r.p_) ),   xy_( *(r.p_+1) ), xz_( *(r.p_+2) ),
		yx_( *(r.p_+3) ), yy_( *(r.p_+4) ), yz_( *(r.p_+5) ),
		zx_( *(r.p_+6) ), zy_( *(r.p_+7) ), zz_( *(r.p_+8) )
	{}


	/// @brief Pointer to contiguous row-ordered values named constructor/assignment
	/// @warning No way to check that argument points to nine contiguous values
	/// @note  Constructor usage: xyzMatrix m( xyzMatrix::rows( rows_pointer ) )
	/// @note  Assignment usage: m = xyzMatrix::rows( rows_pointer )
	template< typename U >
	inline
	static
	RowsPointer< U >
	rows( U const * rp )
	{
		return RowsPointer< U >( rp );
	}


	/// @brief Pointers to contiguous rows constructor
	/// @see   rows
	/// @note  The rows named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::rows( x_row_pointer, y_row_pointer, z_row_pointer ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::rows( x_row_pointer, y_row_pointer, z_row_pointer )
	template< typename U >
	inline
	xyzMatrix( RowPointers< U > const & r ) :
		xx_( *r.xp_ ), xy_( *(r.xp_+1) ), xz_( *(r.xp_+2) ),
		yx_( *r.yp_ ), yy_( *(r.yp_+1) ), yz_( *(r.yp_+2) ),
		zx_( *r.zp_ ), zy_( *(r.zp_+1) ), zz_( *(r.zp_+2) )
	{}


	/// @brief Pointers to contiguous rows named constructor/assignment
	/// @warning No way to check that arguments each point to three contiguous values
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::rows( x_row_pointer, y_row_pointer, z_row_pointer ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::rows( x_row_pointer, y_row_pointer, z_row_pointer )
	template< typename U >
	inline
	static
	RowPointers< U >
	rows(
		U const * xp,
		U const * yp,
		U const * zp
	)
	{
		return RowPointers< U >( xp, yp, zp );
	}


	/// @brief Pointers to contiguous rows named constructor
	/// @warning No way to check that arguments each point to three contiguous values
	/// @note  Can be faster than rows() construction due to return value optimization
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::rows_constructor( x_row_pointer, y_row_pointer, z_row_pointer ) )
	template< typename U >
	inline
	static
	xyzMatrix
	rows_constructor(
		U const * xp,
		U const * yp,
		U const * zp
	)
	{
		return xyzMatrix(
			*xp, *(xp+1), *(xp+2),
			*yp, *(yp+1), *(yp+2),
			*zp, *(zp+1), *(zp+2)
		);
	}


public: // Creation: row xyzVectors


	/// @brief Rows constructor
	/// @see   rows
	/// @note  The rows named constructor/assignment provides a simple wrapper interface
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::rows( x_row_vector, y_row_vector, z_row_vector ) )
	template< typename U >
	inline
	xyzMatrix( RowVectors< U > const & r ) :
		xx_( r.x_.x_ ), xy_( r.x_.y_ ), xz_( r.x_.z_ ),
		yx_( r.y_.x_ ), yy_( r.y_.y_ ), yz_( r.y_.z_ ),
		zx_( r.z_.x_ ), zy_( r.z_.y_ ), zz_( r.z_.z_ )
	{}


	/// @brief Row vectors named constructor/assignment
	/// @note  Constructor usage:
	///  xyzMatrix m( xyzMatrix::rows( x_row_vector, y_row_vector, z_row_vector ) )
	/// @note  Assignment usage:
	///  m = xyzMatrix::rows( x_row_vector, y_row_vector, z_row_vector )
	template< typename U >
	inline
	static
	RowVectors< U >
	rows(
		xyzVector< U > const & row_x,
		xyzVector< U > const & row_y,
		xyzVector< U > const & row_z
	)
	{
		return RowVectors< U >( row_x, row_y, row_z );
	}


	/// @brief xyzVector rows named constructor
	/// @note  Can be faster than rows() construction due to return value optimization
	/// @note  Constructor usage:
	///  xyzMatrix m(
	///   xyzMatrix::rows_constructor( x_row_vector, y_row_vector, z_row_vector ) )
	template< typename U >
	inline
	static
	xyzMatrix
	rows_constructor(
		xyzVector< U > const & row_x,
		xyzVector< U > const & row_y,
		xyzVector< U > const & row_z
	)
	{
		return xyzMatrix(
			row_x.x_, row_x.y_, row_x.z_,
			row_y.x_, row_y.y_, row_y.z_,
			row_z.x_, row_z.y_, row_z.z_
		);
	}


public: // Creation: special constructors


	/// @brief Diagonal value named constructor
	/// @note  diag refers to a diagonal matrix
	inline
	static
	xyzMatrix
	diag(
		Value const & xx_a,
		Value const & yy_a,
		Value const & zz_a
	)
	{
		return xyzMatrix(
			xx_a,        Value( 0 ),  Value( 0 ),
			Value( 0 ),  yy_a,        Value( 0 ),
			Value( 0 ),  Value( 0 ),  zz_a
		);
	}


	/// @brief Diagonal xyzVector named constructor
	/// @note  diag refers to a diagonal matrix
	template< typename U >
	inline
	static
	xyzMatrix
	diag( xyzVector< U > const & diag_a )
	{
		return xyzMatrix(
			diag_a.x_,  Value( 0 ), Value( 0 ),
			Value( 0 ), diag_a.y_,  Value( 0 ),
			Value( 0 ), Value( 0 ), diag_a.z_
		);
	}


	/// @brief Identity xyzMatrix named constructor
	/// @note  Can be faster than I() for construction due to return value optimization
	inline
	static
	xyzMatrix
	identity()
	{
		return xyzMatrix(
			Value( 1 ), Value( 0 ), Value( 0 ),
			Value( 0 ), Value( 1 ), Value( 0 ),
			Value( 0 ), Value( 0 ), Value( 1 )
		);
	}


private: // Creation


	/// @brief Row-ordered value constructor
	/// @note  Client code uses named constructors that specify the value order explicitly
	inline
	xyzMatrix(
		Value const & xx_a, Value const & xy_a, Value const & xz_a,
		Value const & yx_a, Value const & yy_a, Value const & yz_a,
		Value const & zx_a, Value const & zy_a, Value const & zz_a
	) :
		xx_( xx_a ), xy_( xy_a ), xz_( xz_a ),
		yx_( yx_a ), yy_( yy_a ), yz_( yz_a ),
		zx_( zx_a ), zy_( zy_a ), zz_( zz_a )
	{}


public: // Assignment: xyzMatrix


	/// @brief Copy assignment
	inline
	xyzMatrix &
	operator =( xyzMatrix const & m )
	{
		if ( this != &m ) {
			xx_ = m.xx_; xy_ = m.xy_; xz_ = m.xz_;
			yx_ = m.yx_; yy_ = m.yy_; yz_ = m.yz_;
			zx_ = m.zx_; zy_ = m.zy_; zz_ = m.zz_;
		}
		return *this;
	}


	/// @brief Copy assignment
	template< typename U >
	inline
	xyzMatrix &
	operator =( xyzMatrix< U > const & m )
	{
		xx_ = m.xx_; xy_ = m.xy_; xz_ = m.xz_;
		yx_ = m.yx_; yy_ = m.yy_; yz_ = m.yz_;
		zx_ = m.zx_; zy_ = m.zy_; zz_ = m.zz_;
		return *this;
	}


	/// @brief += xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	operator +=( xyzMatrix< U > const & m )
	{
		xx_ += m.xx_; xy_ += m.xy_; xz_ += m.xz_;
		yx_ += m.yx_; yy_ += m.yy_; yz_ += m.yz_;
		zx_ += m.zx_; zy_ += m.zy_; zz_ += m.zz_;
		return *this;
	}


	/// @brief -= xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	operator -=( xyzMatrix< U > const & m )
	{
		xx_ -= m.xx_; xy_ -= m.xy_; xz_ -= m.xz_;
		yx_ -= m.yx_; yy_ -= m.yy_; yz_ -= m.yz_;
		zx_ -= m.zx_; zy_ -= m.zy_; zz_ -= m.zz_;
		return *this;
	}


	/// @brief *= xyzMatrix
	/// @note  Same as right_multiply_by( xyzMatrix )
	template< typename U >
	inline
	xyzMatrix &
	operator *=( xyzMatrix< U > const & m )
	{
		Value x, y, z; // Temporaries

		// First row
		x = ( xx_ * m.xx_ ) + ( xy_ * m.yx_ ) + ( xz_ * m.zx_ );
		y = ( xx_ * m.xy_ ) + ( xy_ * m.yy_ ) + ( xz_ * m.zy_ );
		z = ( xx_ * m.xz_ ) + ( xy_ * m.yz_ ) + ( xz_ * m.zz_ );
		xx_ = x; xy_ = y; xz_ = z;

		// Second row
		x = ( yx_ * m.xx_ ) + ( yy_ * m.yx_ ) + ( yz_ * m.zx_ );
		y = ( yx_ * m.xy_ ) + ( yy_ * m.yy_ ) + ( yz_ * m.zy_ );
		z = ( yx_ * m.xz_ ) + ( yy_ * m.yz_ ) + ( yz_ * m.zz_ );
		yx_ = x; yy_ = y; yz_ = z;

		// Third row
		x = ( zx_ * m.xx_ ) + ( zy_ * m.yx_ ) + ( zz_ * m.zx_ );
		y = ( zx_ * m.xy_ ) + ( zy_ * m.yy_ ) + ( zz_ * m.zy_ );
		z = ( zx_ * m.xz_ ) + ( zy_ * m.yz_ ) + ( zz_ * m.zz_ );
		zx_ = x; zy_ = y; zz_ = z;

		return *this;
	}


public: // Assignment: pointer


	/// @brief Assignment from pointer to contiguous column-ordered values
	/// @note  Use via named assignment wrapper: m = xyzMatrix::cols( pointer )
	template< typename U >
	inline
	xyzMatrix &
	operator =( ColsPointer< U > const & c )
	{
		xx_ = *(c.p_);   xy_ = *(c.p_+3); xz_ = *(c.p_+6);
		yx_ = *(c.p_+1); yy_ = *(c.p_+4); yz_ = *(c.p_+7);
		zx_ = *(c.p_+2); zy_ = *(c.p_+5); zz_ = *(c.p_+8);
		return *this;
	}


	/// @brief Assignment from pointer to contiguous row-ordered values
	/// @note  Use via named assignment wrapper: m = xyzMatrix::rows( pointer )
	template< typename U >
	inline
	xyzMatrix &
	operator =( RowsPointer< U > const & r )
	{
		xx_ = *(r.p_);   xy_ = *(r.p_+1); xz_ = *(r.p_+2);
		yx_ = *(r.p_+3); yy_ = *(r.p_+4); yz_ = *(r.p_+5);
		zx_ = *(r.p_+6); zy_ = *(r.p_+7); zz_ = *(r.p_+8);
		return *this;
	}


	/// @brief Assignment from pointers to contiguous columns
	/// @note  Use via named assignment wrapper:
	///  m = xyzMatrix::cols( pointer, pointer, pointer )
	template< typename U >
	inline
	xyzMatrix &
	operator =( ColPointers< U > const & c )
	{
		xx_ = *c.xp_;     xy_ = *c.yp_;     xz_ = *c.zp_;
		yx_ = *(c.xp_+1); yy_ = *(c.yp_+1); yz_ = *(c.zp_+1);
		zx_ = *(c.xp_+2); zy_ = *(c.yp_+2); zz_ = *(c.zp_+2);
		return *this;
	}


	/// @brief Assignment from pointers to contiguous rows
	/// @note  Use via named assignment wrapper:
	///  m = xyzMatrix::rows( pointer, pointer, pointer )
	template< typename U >
	inline
	xyzMatrix &
	operator =( RowPointers< U > const & r )
	{
		xx_ = *r.xp_; xy_ = *(r.xp_+1); xz_ = *(r.xp_+2);
		yx_ = *r.yp_; yy_ = *(r.yp_+1); yz_ = *(r.yp_+2);
		zx_ = *r.zp_; zy_ = *(r.zp_+1); zz_ = *(r.zp_+2);
		return *this;
	}


public: // Assignment: xyzVector


	/// @brief xyzVector columns assignment
	/// @note  Use via named assignment wrapper:
	///  m = xyzMatrix::cols( xyzVector, xyzVector, xyzVector )
	template< typename U >
	inline
	xyzMatrix &
	operator =( ColVectors< U > const & c )
	{
		xx_ = c.x_.x_; xy_ = c.y_.x_; xz_ = c.z_.x_;
		yx_ = c.x_.y_; yy_ = c.y_.y_; yz_ = c.z_.y_;
		zx_ = c.x_.z_; zy_ = c.y_.z_; zz_ = c.z_.z_;
		return *this;
	}


	/// @brief xyzVector rows assignment
	/// @note  Use via named assignment wrapper:
	///  m = xyzMatrix::rows( xyzVector, xyzVector, xyzVector )
	template< typename U >
	inline
	xyzMatrix &
	operator =( RowVectors< U > const & r )
	{
		xx_ = r.x_.x_; xy_ = r.x_.y_; xz_ = r.x_.z_;
		yx_ = r.y_.x_; yy_ = r.y_.y_; yz_ = r.y_.z_;
		zx_ = r.z_.x_; zy_ = r.z_.y_; zz_ = r.z_.z_;
		return *this;
	}


public: // Assignment: scalar


	/// @brief = Value
	inline
	xyzMatrix &
	operator =( Value const & t )
	{
		xx_ = xy_ = xz_ = t;
		yx_ = yy_ = yz_ = t;
		zx_ = zy_ = zz_ = t;
		return *this;
	}


	/// @brief += Value
	inline
	xyzMatrix &
	operator +=( Value const & t )
	{
		xx_ += t; xy_ += t; xz_ += t;
		yx_ += t; yy_ += t; yz_ += t;
		zx_ += t; zy_ += t; zz_ += t;
		return *this;
	}


	/// @brief -= Value
	inline
	xyzMatrix &
	operator -=( Value const & t )
	{
		xx_ -= t; xy_ -= t; xz_ -= t;
		yx_ -= t; yy_ -= t; yz_ -= t;
		zx_ -= t; zy_ -= t; zz_ -= t;
		return *this;
	}


	/// @brief *= Value
	inline
	xyzMatrix &
	operator *=( Value const & t )
	{
		xx_ *= t; xy_ *= t; xz_ *= t;
		yx_ *= t; yy_ *= t; yz_ *= t;
		zx_ *= t; zy_ *= t; zz_ *= t;
		return *this;
	}


	/// @brief /= Value
	inline
	xyzMatrix &
	operator /=( Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value ( 1 ) / t );
		xx_ *= inv_t; xy_ *= inv_t; xz_ *= inv_t;
		yx_ *= inv_t; yy_ *= inv_t; yz_ *= inv_t;
		zx_ *= inv_t; zy_ *= inv_t; zz_ *= inv_t;
		return *this;
	}


public: // Methods: basic mathematical


	/// @brief xyzMatrix * xyzVector
	/// @note  Same as product( xyzMatrix, xyzVector )
	inline
	xyzVector<T>
	operator *(xyzVector<T> const & v ) const
	{
		return xyzVector<T>(
			xx_ * v.x_ + xy_ * v.y_ + xz_ * v.z_,
			yx_ * v.x_ + yy_ * v.y_ + yz_ * v.z_,
			zx_ * v.x_ + zy_ * v.y_ + zz_ * v.z_
						 );
	}


	/// @brief xyzMatrix + xyzMatrix
	friend
	inline
	xyzMatrix
	operator +( xyzMatrix const & a, xyzMatrix const & b )
	{
		return xyzMatrix(
			a.xx_ + b.xx_, a.xy_ + b.xy_, a.xz_ + b.xz_,
			a.yx_ + b.yx_, a.yy_ + b.yy_, a.yz_ + b.yz_,
			a.zx_ + b.zx_, a.zy_ + b.zy_, a.zz_ + b.zz_
		);
	}


	/// @brief xyzMatrix + Value
	friend
	inline
	xyzMatrix
	operator +( xyzMatrix const & m, Value const & t )
	{
		return xyzMatrix(
			m.xx_ + t, m.xy_ + t, m.xz_ + t,
			m.yx_ + t, m.yy_ + t, m.yz_ + t,
			m.zx_ + t, m.zy_ + t, m.zz_ + t
		);
	}


	/// @brief Value + xyzMatrix
	friend
	inline
	xyzMatrix
	operator +( Value const & t, xyzMatrix const & m )
	{
		return xyzMatrix(
			t + m.xx_, t + m.xy_, t + m.xz_,
			t + m.yx_, t + m.yy_, t + m.yz_,
			t + m.zx_, t + m.zy_, t + m.zz_
		);
	}


	/// @brief xyzMatrix - xyzMatrix
	friend
	inline
	xyzMatrix
	operator -( xyzMatrix const & a, xyzMatrix const & b )
	{
		return xyzMatrix(
			a.xx_ - b.xx_, a.xy_ - b.xy_, a.xz_ - b.xz_,
			a.yx_ - b.yx_, a.yy_ - b.yy_, a.yz_ - b.yz_,
			a.zx_ - b.zx_, a.zy_ - b.zy_, a.zz_ - b.zz_
		);
	}


	/// @brief xyzMatrix - Value
	friend
	inline
	xyzMatrix
	operator -( xyzMatrix const & m, Value const & t )
	{
		return xyzMatrix(
			m.xx_ - t, m.xy_ - t, m.xz_ - t,
			m.yx_ - t, m.yy_ - t, m.yz_ - t,
			m.zx_ - t, m.zy_ - t, m.zz_ - t
		);
	}


	/// @brief Value - xyzMatrix
	friend
	inline
	xyzMatrix
	operator -( Value const & t, xyzMatrix const & m )
	{
		return xyzMatrix(
			t - m.xx_, t - m.xy_, t - m.xz_,
			t - m.yx_, t - m.yy_, t - m.yz_,
			t - m.zx_, t - m.zy_, t - m.zz_
		);
	}


	/// @brief xyzMatrix * xyzMatrix
	friend
	inline
	xyzMatrix
	operator *( xyzMatrix const & a, xyzMatrix const & b )
	{
		return xyzMatrix(
			// First row
			( a.xx_ * b.xx_ ) + ( a.xy_ * b.yx_ ) + ( a.xz_ * b.zx_ ),
			( a.xx_ * b.xy_ ) + ( a.xy_ * b.yy_ ) + ( a.xz_ * b.zy_ ),
			( a.xx_ * b.xz_ ) + ( a.xy_ * b.yz_ ) + ( a.xz_ * b.zz_ ),

			// Second row
			( a.yx_ * b.xx_ ) + ( a.yy_ * b.yx_ ) + ( a.yz_ * b.zx_ ),
			( a.yx_ * b.xy_ ) + ( a.yy_ * b.yy_ ) + ( a.yz_ * b.zy_ ),
			( a.yx_ * b.xz_ ) + ( a.yy_ * b.yz_ ) + ( a.yz_ * b.zz_ ),

			// Third row
			( a.zx_ * b.xx_ ) + ( a.zy_ * b.yx_ ) + ( a.zz_ * b.zx_ ),
			( a.zx_ * b.xy_ ) + ( a.zy_ * b.yy_ ) + ( a.zz_ * b.zy_ ),
			( a.zx_ * b.xz_ ) + ( a.zy_ * b.yz_ ) + ( a.zz_ * b.zz_ )
		);
	}


	/// @brief xyzMatrix * Value
	friend
	inline
	xyzMatrix
	operator *( xyzMatrix const & m, Value const & t )
	{
		return xyzMatrix(
			m.xx_ * t, m.xy_ * t, m.xz_ * t,
			m.yx_ * t, m.yy_ * t, m.yz_ * t,
			m.zx_ * t, m.zy_ * t, m.zz_ * t
		);
	}


	/// @brief Value * xyzMatrix
	friend
	inline
	xyzMatrix
	operator *( Value const & t, xyzMatrix const & m )
	{
		return xyzMatrix(
			t * m.xx_, t * m.xy_, t * m.xz_,
			t * m.yx_, t * m.yy_, t * m.yz_,
			t * m.zx_, t * m.zy_, t * m.zz_
		);
	}


	/// @brief xyzMatrix / Value
	friend
	inline
	xyzMatrix
	operator /( xyzMatrix const & m, Value const & t )
	{
		assert( t != Value( 0 ) );
		Value const inv_t( Value( 1 ) / t );
		return xyzMatrix(
			m.xx_ * inv_t, m.xy_ * inv_t, m.xz_ * inv_t,
			m.yx_ * inv_t, m.yy_ * inv_t, m.yz_ * inv_t,
			m.zx_ * inv_t, m.zy_ * inv_t, m.zz_ * inv_t
		);
	}


public: // Methods: complex mathematical


	/// @brief Clear
	inline
	xyzMatrix &
	clear()
	{
		xx_ = xy_ = xz_ =
			yx_ = yy_ = yz_ =
			zx_ = zy_ = zz_ = Value( 0 );
		return *this;
	}


	/// @brief Set to the zero xyzMatrix
	inline
	xyzMatrix &
	zero()
	{
		xx_ = xy_ = xz_ =
			yx_ = yy_ = yz_ =
			zx_ = zy_ = zz_ = Value( 0 );
		return *this;
	}


	/// @brief Set to the identity xyzMatrix
	inline
	xyzMatrix &
	to_identity()
	{
		xx_ = yy_ = zz_ = Value( 1 );
		xy_ = xz_ = yx_ = yz_ = zx_ = zy_ = Value( 0 );
		return *this;
	}


	/// @brief Set to diagonal xyzMatrix from value
	/// @note  Resets the entire matrix (diag refers to a diagonal matrix)
	inline
	xyzMatrix &
	to_diag(
		Value const & xx_a,
		Value const & yy_a,
		Value const & zz_a
	)
	{
		xx_ = xx_a;
		yy_ = yy_a;
		zz_ = zz_a;
		xy_ = xz_ = yx_ = yz_ = zx_ = zy_ = Value( 0 );
		return *this;
	}


	/// @brief Set to diagonal xyzMatrix from xyzVector
	/// @note  Resets the entire matrix (diag refers to a diagonal matrix)
	template< typename U >
	inline
	xyzMatrix &
	to_diag( xyzVector< U > const & diag_a )
	{
		xx_ = diag_a.x_;
		yy_ = diag_a.y_;
		zz_ = diag_a.z_;
		xy_ = xz_ = yx_ = yz_ = zx_ = zy_ = Value( 0 );
		return *this;
	}


	/// @brief set diagonal of xyzMatrix from value
	/// @note  Resets the diagonal of the matrix only (diagonal refers to diagonal
	///   entries of a matrix)
	inline
	xyzMatrix &
	set_diagonal(
		Value const & xx_a,
		Value const & yy_a,
		Value const & zz_a
	)
	{
		xx_ = xx_a;
		yy_ = yy_a;
		zz_ = zz_a;
		return *this;
	}


	/// @brief Set diagonal of xyzMatrix from xyzVector
	/// @note  Resets the diagonal of the matrix only (diagonal refers to diagonal
	//    entries of a matrix)
	template< typename U >
	inline
	xyzMatrix &
	set_diagonal( xyzVector< U > const & diag_a )
	{
		xx_ = diag_a.x_;
		yy_ = diag_a.y_;
		zz_ = diag_a.z_;
		return *this;
	}


	/// @brief Add values to diagonal of xyzMatrix
	inline
	xyzMatrix &
	add_diagonal(
		Value const & xx_a,
		Value const & yy_a,
		Value const & zz_a
	)
	{
		xx_ += xx_a;
		yy_ += yy_a;
		zz_ += zz_a;
		return *this;
	}


	/// @brief Add xyzVector to diagonal of xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	add_diagonal( xyzVector< U > const & diag_a )
	{
		xx_ += diag_a.x_;
		yy_ += diag_a.y_;
		zz_ += diag_a.z_;
		return *this;
	}


	/// @brief Subtract values from diagonal of xyzMatrix
	inline
	xyzMatrix &
	subtract_diagonal(
		Value const & xx_a,
		Value const & yy_a,
		Value const & zz_a
	)
	{
		xx_ -= xx_a;
		yy_ -= yy_a;
		zz_ -= zz_a;
		return *this;
	}


	/// @brief Subtract xyzVector from diagonal of xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	subtract_diagonal( xyzVector< U > const & diag_a )
	{
		xx_ -= diag_a.x_;
		yy_ -= diag_a.y_;
		zz_ -= diag_a.z_;
		return *this;
	}


	/// @brief Transpose
	inline
	xyzMatrix &
	transpose()
	{
		Value temp = xy_;
		xy_ = yx_;
		yx_ = temp;

		temp = xz_;
		xz_ = zx_;
		zx_ = temp;

		temp = yz_;
		yz_ = zy_;
		zy_ = temp;

		return *this;
	}


	/// @brief Right multiply by xyzMatrix
	/// @note  Same as *= xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	right_multiply_by( xyzMatrix< U > const & m )
	{
		Value x, y, z; // Temporaries

		// First row
		x = ( xx_ * m.xx_ ) + ( xy_ * m.yx_ ) + ( xz_ * m.zx_ );
		y = ( xx_ * m.xy_ ) + ( xy_ * m.yy_ ) + ( xz_ * m.zy_ );
		z = ( xx_ * m.xz_ ) + ( xy_ * m.yz_ ) + ( xz_ * m.zz_ );
		xx_ = x; xy_ = y; xz_ = z;

		// Second row
		x = ( yx_ * m.xx_ ) + ( yy_ * m.yx_ ) + ( yz_ * m.zx_ );
		y = ( yx_ * m.xy_ ) + ( yy_ * m.yy_ ) + ( yz_ * m.zy_ );
		z = ( yx_ * m.xz_ ) + ( yy_ * m.yz_ ) + ( yz_ * m.zz_ );
		yx_ = x; yy_ = y; yz_ = z;

		// Third row
		x = ( zx_ * m.xx_ ) + ( zy_ * m.yx_ ) + ( zz_ * m.zx_ );
		y = ( zx_ * m.xy_ ) + ( zy_ * m.yy_ ) + ( zz_ * m.zy_ );
		z = ( zx_ * m.xz_ ) + ( zy_ * m.yz_ ) + ( zz_ * m.zz_ );
		zx_ = x; zy_ = y; zz_ = z;

		return *this;
	}


	/// @brief Right multiply by transpose xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	right_multiply_by_transpose( xyzMatrix< U > const & m )
	{
		Value x, y, z; // Temporaries

		// First row
		x = ( xx_ * m.xx_ ) + ( xy_ * m.xy_ ) + ( xz_ * m.xz_ );
		y = ( xx_ * m.yx_ ) + ( xy_ * m.yy_ ) + ( xz_ * m.yz_ );
		z = ( xx_ * m.zx_ ) + ( xy_ * m.zy_ ) + ( xz_ * m.zz_ );
		xx_ = x; xy_ = y; xz_ = z;

		// Second row
		x = ( yx_ * m.xx_ ) + ( yy_ * m.xy_ ) + ( yz_ * m.xz_ );
		y = ( yx_ * m.yx_ ) + ( yy_ * m.yy_ ) + ( yz_ * m.yz_ );
		z = ( yx_ * m.zx_ ) + ( yy_ * m.zy_ ) + ( yz_ * m.zz_ );
		yx_ = x; yy_ = y; yz_ = z;

		// Third row
		x = ( zx_ * m.xx_ ) + ( zy_ * m.xy_ ) + ( zz_ * m.xz_ );
		y = ( zx_ * m.yx_ ) + ( zy_ * m.yy_ ) + ( zz_ * m.yz_ );
		z = ( zx_ * m.zx_ ) + ( zy_ * m.zy_ ) + ( zz_ * m.zz_ );
		zx_ = x; zy_ = y; zz_ = z;

		return *this;
	}


	/// @brief Left multiply by xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	left_multiply_by( xyzMatrix< U > const & m )
	{
		Value x, y, z; // Temporaries

		// First column
		x = ( m.xx_ * xx_ ) + ( m.xy_ * yx_ ) + ( m.xz_ * zx_ );
		y = ( m.yx_ * xx_ ) + ( m.yy_ * yx_ ) + ( m.yz_ * zx_ );
		z = ( m.zx_ * xx_ ) + ( m.zy_ * yx_ ) + ( m.zz_ * zx_ );
		xx_ = x; yx_ = y; zx_ = z;

		// Second column
		x = ( m.xx_ * xy_ ) + ( m.xy_ * yy_ ) + ( m.xz_ * zy_ );
		y = ( m.yx_ * xy_ ) + ( m.yy_ * yy_ ) + ( m.yz_ * zy_ );
		z = ( m.zx_ * xy_ ) + ( m.zy_ * yy_ ) + ( m.zz_ * zy_ );
		xy_ = x; yy_ = y; zy_ = z;

		// Third column
		x = ( m.xx_ * xz_ ) + ( m.xy_ * yz_ ) + ( m.xz_ * zz_ );
		y = ( m.yx_ * xz_ ) + ( m.yy_ * yz_ ) + ( m.yz_ * zz_ );
		z = ( m.zx_ * xz_ ) + ( m.zy_ * yz_ ) + ( m.zz_ * zz_ );
		xz_ = x; yz_ = y; zz_ = z;

		return *this;
	}


	/// @brief Left multiply by transpose xyzMatrix
	template< typename U >
	inline
	xyzMatrix &
	left_multiply_by_transpose( xyzMatrix< U > const & m )
	{
		Value x, y, z; // Temporaries

		// First column
		x = ( m.xx_ * xx_ ) + ( m.yx_ * yx_ ) + ( m.zx_ * zx_ );
		y = ( m.xy_ * xx_ ) + ( m.yy_ * yx_ ) + ( m.zy_ * zx_ );
		z = ( m.xz_ * xx_ ) + ( m.yz_ * yx_ ) + ( m.zz_ * zx_ );
		xx_ = x; yx_ = y; zx_ = z;

		// Second column
		x = ( m.xx_ * xy_ ) + ( m.yx_ * yy_ ) + ( m.zx_ * zy_ );
		y = ( m.xy_ * xy_ ) + ( m.yy_ * yy_ ) + ( m.zy_ * zy_ );
		z = ( m.xz_ * xy_ ) + ( m.yz_ * yy_ ) + ( m.zz_ * zy_ );
		xy_ = x; yy_ = y; zy_ = z;

		// Third column
		x = ( m.xx_ * xz_ ) + ( m.yx_ * yz_ ) + ( m.zx_ * zz_ );
		y = ( m.xy_ * xz_ ) + ( m.yy_ * yz_ ) + ( m.zy_ * zz_ );
		z = ( m.xz_ * xz_ ) + ( m.yz_ * yz_ ) + ( m.zz_ * zz_ );
		xz_ = x; yz_ = y; zz_ = z;

		return *this;
	}


	/// @brief Identity xyzMatrix for expressions
	/// @note  identity() named constructor can be faster for construction
	/// @note  Returns reference to function-local static object so can be used in
	///  construction of global objects
	inline
	static
	xyzMatrix const &
	I()
	{
		static xyzMatrix const I_(
			Value( 1 ), Value( 0 ), Value( 0 ),
			Value( 0 ), Value( 1 ), Value( 0 ),
			Value( 0 ), Value( 0 ), Value( 1 )
		);

		return I_;
	}


public: // Properties: columns


	/// @brief Column x
	inline
	Vector
	col_x() const
	{
		return Vector( xx_, yx_, zx_ );
	}


	/// @brief Column x assignment
	inline
	xyzMatrix &
	col_x( Vector const & v )
	{
		xx_ = v.x_; yx_ = v.y_; zx_ = v.z_;
		return *this;
	}


	/// @brief Column y
	inline
	Vector
	col_y() const
	{
		return Vector( xy_, yy_, zy_ );
	}


	/// @brief Column y assignment
	inline
	xyzMatrix &
	col_y( Vector const & v )
	{
		xy_ = v.x_; yy_ = v.y_; zy_ = v.z_;
		return *this;
	}


	/// @brief Column z
	inline
	Vector
	col_z() const
	{
		return Vector( xz_, yz_, zz_ );
	}


	/// @brief Column z assignment
	inline
	xyzMatrix &
	col_z( Vector const & v)
	{
		xz_ = v.x_; yz_ = v.y_; zz_ = v.z_;
		return *this;
	}


	/// @brief Column( i ): 1-based index
	inline
	Vector
	col( int const i ) const
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		switch ( i ) {
		case 1 :
			return Vector( xx_, yx_, zx_ );
		case 2 :
			return Vector( xy_, yy_, zy_ );
		default : // Assume i == 3
			return Vector( xz_, yz_, zz_ );
		}
	}


	/// @brief Column( i, xyzVector ) assignment: 1-base index
	inline
	xyzMatrix &
	col( int const i, Vector const & v )
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		switch ( i ) {
		case 1 :
			xx_ = v.x_; yx_ = v.y_; zx_ = v.z_;
			break;
		case 2 :
			xy_ = v.x_; yy_ = v.y_; zy_ = v.z_;
			break;
		default : // Assume i == 3
			xz_ = v.x_; yz_ = v.y_; zz_ = v.z_;
		}
		return *this;
	}


public: // Properties: rows


	/// @brief Row x
	inline
	Vector
	row_x() const
	{
		return Vector( xx_, xy_, xz_ );
	}


	/// @brief Row x assignment
	inline
	xyzMatrix &
	row_x( Vector const & v )
	{
		xx_ = v.x_; xy_ = v.y_; xz_ = v.z_;
		return *this;
	}


	/// @brief Row y
	inline
	Vector
	row_y() const
	{
		return Vector( yx_, yy_, yz_ );
	}


	/// @brief Row y assignment
	inline
	xyzMatrix &
	row_y( Vector const & v )
	{
		yx_ = v.x_; yy_ = v.y_; yz_ = v.z_;
		return *this;
	}


	/// @brief Row z
	inline
	Vector
	row_z() const
	{
		return Vector( zx_, zy_, zz_ );
	}


	/// @brief Row z assignment
	inline
	xyzMatrix &
	row_z( Vector const & v )
	{
		zx_ = v.x_; zy_ = v.y_; zz_ = v.z_;
		return *this;
	}


	/// @brief Row ( i ): 1-based index
	inline
	Vector
	row( int const i ) const
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		switch ( i ) {
		case 1 :
			return Vector( xx_, xy_, xz_ );
		case 2 :
			return Vector( yx_, yy_, yz_ );
		default : // Assume i == 3
			return Vector( zx_, zy_, zz_ );
		}
	}


	/// @brief Row ( i, xyzVector ) assignment: 1-based index
	inline
	xyzMatrix &
	row( int const i, Vector const & v )
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		switch ( i ) {
		case 1 :
			xx_ = v.x_; xy_ = v.y_; xz_ = v.z_;
		case 2 :
			yx_ = v.x_; yy_ = v.y_; yz_ = v.z_;
		default : // Assume i == 3
			zx_ = v.x_; zy_ = v.y_; zz_ = v.z_;
		}
		return *this;
	}


public: // Properties: scalars


	/// @brief Value xx const
	inline
	Value const &
	xx() const
	{
		return xx_;
	}


	/// @brief Value xx
	inline
	Value &
	xx()
	{
		return xx_;
	}


	/// @brief Value xy const
	inline
	Value const &
	xy() const
	{
		return xy_;
	}


	/// @brief Value xy
	inline
	Value &
	xy()
	{
		return xy_;
	}


	/// @brief Value xz const
	inline
	Value const &
	xz() const
	{
		return xz_;
	}


	/// @brief Value xz
	inline
	Value &
	xz()
	{
		return xz_;
	}


	/// @brief Value yx const
	inline
	Value const &
	yx() const
	{
		return yx_;
	}


	/// @brief Value yx
	inline
	Value &
	yx()
	{
		return yx_;
	}


	/// @brief Value yy const
	inline
	Value const &
	yy() const
	{
		return yy_;
	}


	/// @brief Value yy
	inline
	Value &
	yy()
	{
		return yy_;
	}


	/// @brief Value yz const
	inline
	Value const &
	yz() const
	{
		return yz_;
	}


	/// @brief Value yz
	inline
	Value &
	yz()
	{
		return yz_;
	}


	/// @brief Value zx const
	inline
	Value const &
	zx() const
	{
		return zx_;
	}


	/// @brief Value zx
	inline
	Value &
	zx()
	{
		return zx_;
	}


	/// @brief Value zy const
	inline
	Value const &
	zy() const
	{
		return zy_;
	}


	/// @brief Value zy
	inline
	Value &
	zy()
	{
		return zy_;
	}


	/// @brief Value zz const
	inline
	Value const &
	zz() const
	{
		return zz_;
	}


	/// @brief Value zz
	inline
	Value &
	zz()
	{
		return zz_;
	}


	/// @brief xyzMatrix( i, j ) const: 1-based index
	inline
	Value const &
	operator ()( int const i, int const j ) const
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		assert( ( j > 0 ) && ( j <= 3 ) );
		switch ( i ) {
		case 1 :
			return ( j == 1 ? xx_ : ( j == 2 ? xy_ : xz_ ) );
		case 2 :
			return ( j == 1 ? yx_ : ( j == 2 ? yy_ : yz_ ) );
		default : // Assume i == 3
			return ( j == 1 ? zx_ : ( j == 2 ? zy_ : zz_ ) );
		}
	}


	/// @brief xyzMatrix( i, j ): 1-based index
	inline
	Value &
	operator ()( int const i, int const j )
	{
		assert( ( i > 0 ) && ( i <= 3 ) );
		assert( ( j > 0 ) && ( j <= 3 ) );
		switch ( i ) {
		case 1 :
			return ( j == 1 ? xx_ : ( j == 2 ? xy_ : xz_ ) );
		case 2 :
			return ( j == 1 ? yx_ : ( j == 2 ? yy_ : yz_ ) );
		default : // Assume i == 3
			return ( j == 1 ? zx_ : ( j == 2 ? zy_ : zz_ ) );
		}
	}


public: // Properties: value assignment


	/// @brief xx assignment
	inline
	void
	xx( Value const & xx_a )
	{
		xx_ = xx_a;
	}


	/// @brief xy assignment
	inline
	void
	xy( Value const & xy_a )
	{
		xy_ = xy_a;
	}


	/// @brief xz assignment
	inline
	void
	xz( Value const & xz_a )
	{
		xz_ = xz_a;
	}


	/// @brief yx assignment
	inline
	void
	yx( Value const & yx_a )
	{
		yx_ = yx_a;
	}


	/// @brief yy assignment
	inline
	void
	yy( Value const & yy_a )
	{
		yy_ = yy_a;
	}


	/// @brief yz assignment
	inline
	void
	yz( Value const & yz_a )
	{
		yz_ = yz_a;
	}


	/// @brief zx assignment
	inline
	void
	zx( Value const & zx_a )
	{
		zx_ = zx_a;
	}


	/// @brief zy assignment
	inline
	void
	zy( Value const & zy_a )
	{
		zy_ = zy_a;
	}


	/// @brief zz assignment
	inline
	void
	zz( Value const & zz_a )
	{
		zz_ = zz_a;
	}


public: // Properties: predicates


	/// @brief Is zero?
	inline
	bool
	is_zero() const
	{
		static Value const ZERO( 0 );
		return
			( xx_ == ZERO ) && ( xy_ == ZERO ) && ( xz_ == ZERO ) &&
			( yx_ == ZERO ) && ( yy_ == ZERO ) && ( yz_ == ZERO ) &&
			( zx_ == ZERO ) && ( zy_ == ZERO ) && ( zz_ == ZERO );
	}


	/// @brief Is identity?
	inline
	bool
	is_identity() const
	{
		static Value const ZERO( 0 );
		static Value const ONE( 1 );
		return
			( xx_ == ONE  ) && ( xy_ == ZERO ) && ( xz_ == ZERO ) &&
			( yx_ == ZERO ) && ( yy_ == ONE  ) && ( yz_ == ZERO ) &&
			( zx_ == ZERO ) && ( zy_ == ZERO ) && ( zz_ == ONE  );
	}


public: // Properties: calculated


	/// @brief Determinant
	inline
	Value
	det() const
	{
		return
			xx_ * ( ( yy_ * zz_ ) - ( zy_ * yz_ ) ) -
			xy_ * ( ( yx_ * zz_ ) - ( zx_ * yz_ ) ) +
			xz_ * ( ( yx_ * zy_ ) - ( zx_ * yy_ ) );
	}

	/// @brief Trace
	inline
	Value
	trace() const
	{
		return xx_ + yy_ + zz_;
	}


	/// @brief Transposed copy
	inline
	xyzMatrix
	transposed() const
	{
		return xyzMatrix(
			xx_, yx_, zx_,
			xy_, yy_, zy_,
			xz_, yz_, zz_
		);
	}

	// Adding member version so we can use M.inverse in PyRosetta
	inline
	xyzMatrix < T >
	inverse() const {
		xyzMatrix const & a(*this);
		T D = a.det();
		return xyzMatrix< T >(
			(a.yy_*a.zz_-a.yz_*a.zy_)/D, -(a.xy_*a.zz_-a.xz_*a.zy_)/D,  (a.xy_*a.yz_-a.xz_*a.yy_)/D,
			-(a.yx_*a.zz_-a.zx_*a.yz_)/D,  (a.xx_*a.zz_-a.xz_*a.zx_)/D, -(a.xx_*a.yz_-a.xz_*a.yx_)/D,
			(a.yx_*a.zy_-a.zx_*a.yy_)/D, -(a.xx_*a.zy_-a.xy_*a.zx_)/D,  (a.xx_*a.yy_-a.xy_*a.yx_)/D
		);
	}

public: // Comparison


	/// @brief xyzMatrix == xyzMatrix
	friend
	inline
	bool
	operator ==( xyzMatrix const & a, xyzMatrix const & b )
	{
		return
			( a.xx_ == b.xx_ ) && ( a.xy_ == b.xy_ ) && ( a.xz_ == b.xz_ ) &&
			( a.yx_ == b.yx_ ) && ( a.yy_ == b.yy_ ) && ( a.yz_ == b.yz_ ) &&
			( a.zx_ == b.zx_ ) && ( a.zy_ == b.zy_ ) && ( a.zz_ == b.zz_ );
	}


	/// @brief xyzMatrix != xyzMatrix
	friend
	inline
	bool
	operator !=( xyzMatrix const & a, xyzMatrix const & b )
	{
		return !( a == b );
	}


	/// @brief xyzMatrix < xyzMatrix
	friend
	inline
	bool
	operator <( xyzMatrix const & a, xyzMatrix const & b )
	{
		return
			( a.xx_ < b.xx_ ) && ( a.xy_ < b.xy_ ) && ( a.xz_ < b.xz_ ) &&
			( a.yx_ < b.yx_ ) && ( a.yy_ < b.yy_ ) && ( a.yz_ < b.yz_ ) &&
			( a.zx_ < b.zx_ ) && ( a.zy_ < b.zy_ ) && ( a.zz_ < b.zz_ );
	}


	/// @brief xyzMatrix <= xyzMatrix
	friend
	inline
	bool
	operator <=( xyzMatrix const & a, xyzMatrix const & b )
	{
		return
			( a.xx_ <= b.xx_ ) && ( a.xy_ <= b.xy_ ) && ( a.xz_ <= b.xz_ ) &&
			( a.yx_ <= b.yx_ ) && ( a.yy_ <= b.yy_ ) && ( a.yz_ <= b.yz_ ) &&
			( a.zx_ <= b.zx_ ) && ( a.zy_ <= b.zy_ ) && ( a.zz_ <= b.zz_ );
	}


	/// @brief xyzMatrix >= xyzMatrix
	friend
	inline
	bool
	operator >=( xyzMatrix const & a, xyzMatrix const & b )
	{
		return
			( a.xx_ >= b.xx_ ) && ( a.xy_ >= b.xy_ ) && ( a.xz_ >= b.xz_ ) &&
			( a.yx_ >= b.yx_ ) && ( a.yy_ >= b.yy_ ) && ( a.yz_ >= b.yz_ ) &&
			( a.zx_ >= b.zx_ ) && ( a.zy_ >= b.zy_ ) && ( a.zz_ >= b.zz_ );
	}


	/// @brief xyzMatrix > xyzMatrix
	friend
	inline
	bool
	operator >( xyzMatrix const & a, xyzMatrix const & b )
	{
		return
			( a.xx_ > b.xx_ ) && ( a.xy_ > b.xy_ ) && ( a.xz_ > b.xz_ ) &&
			( a.yx_ > b.yx_ ) && ( a.yy_ > b.yy_ ) && ( a.yz_ > b.yz_ ) &&
			( a.zx_ > b.zx_ ) && ( a.zy_ > b.zy_ ) && ( a.zz_ > b.zz_ );
	}


	/// @brief xyzMatrix == Value
	friend
	inline
	bool
	operator ==( xyzMatrix const & m, Value const & t )
	{
		return
			( m.xx_ == t ) && ( m.xy_ == t ) && ( m.xz_ == t ) &&
			( m.yx_ == t ) && ( m.yy_ == t ) && ( m.yz_ == t ) &&
			( m.zx_ == t ) && ( m.zy_ == t ) && ( m.zz_ == t );
	}


	/// @brief xyzMatrix != Value
	friend
	inline
	bool
	operator !=( xyzMatrix const & m, Value const & t )
	{
		return !( m == t );
	}


	/// @brief xyzMatrix < Value
	friend
	inline
	bool
	operator <( xyzMatrix const & m, Value const & t )
	{
		return
			( m.xx_ < t ) && ( m.xy_ < t ) && ( m.xz_ < t ) &&
			( m.yx_ < t ) && ( m.yy_ < t ) && ( m.yz_ < t ) &&
			( m.zx_ < t ) && ( m.zy_ < t ) && ( m.zz_ < t );
	}


	/// @brief xyzMatrix <= Value
	friend
	inline
	bool
	operator <=( xyzMatrix const & m, Value const & t )
	{
		return
			( m.xx_ <= t ) && ( m.xy_ <= t ) && ( m.xz_ <= t ) &&
			( m.yx_ <= t ) && ( m.yy_ <= t ) && ( m.yz_ <= t ) &&
			( m.zx_ <= t ) && ( m.zy_ <= t ) && ( m.zz_ <= t );
	}


	/// @brief xyzMatrix >= Value
	friend
	inline
	bool
	operator >=( xyzMatrix const & m, Value const & t )
	{
		return
			( m.xx_ >= t ) && ( m.xy_ >= t ) && ( m.xz_ >= t ) &&
			( m.yx_ >= t ) && ( m.yy_ >= t ) && ( m.yz_ >= t ) &&
			( m.zx_ >= t ) && ( m.zy_ >= t ) && ( m.zz_ >= t );
	}


	/// @brief xyzMatrix > Value
	friend
	inline
	bool
	operator >( xyzMatrix const & m, Value const & t )
	{
		return
			( m.xx_ > t ) && ( m.xy_ > t ) && ( m.xz_ > t ) &&
			( m.yx_ > t ) && ( m.yy_ > t ) && ( m.yz_ > t ) &&
			( m.zx_ > t ) && ( m.zy_ > t ) && ( m.zz_ > t );
	}


	/// @brief Value == xyzMatrix
	friend
	inline
	bool
	operator ==( Value const & t, xyzMatrix const & m )
	{
		return
			( t == m.xx_ ) && ( t == m.xy_ ) && ( t == m.xz_ ) &&
			( t == m.yx_ ) && ( t == m.yy_ ) && ( t == m.yz_ ) &&
			( t == m.zx_ ) && ( t == m.zy_ ) && ( t == m.zz_ );
	}


	/// @brief Value != xyzMatrix
	friend
	inline
	bool
	operator !=( Value const & t, xyzMatrix const & m )
	{
		return !( t == m );
	}


	/// @brief Value < xyzMatrix
	friend
	inline
	bool
	operator <( Value const & t, xyzMatrix const & m )
	{
		return
			( t < m.xx_ ) && ( t < m.xy_ ) && ( t < m.xz_ ) &&
			( t < m.yx_ ) && ( t < m.yy_ ) && ( t < m.yz_ ) &&
			( t < m.zx_ ) && ( t < m.zy_ ) && ( t < m.zz_ );
	}


	/// @brief Value <= xyzMatrix
	friend
	inline
	bool
	operator <=( Value const & t, xyzMatrix const & m )
	{
		return
			( t <= m.xx_ ) && ( t <= m.xy_ ) && ( t <= m.xz_ ) &&
			( t <= m.yx_ ) && ( t <= m.yy_ ) && ( t <= m.yz_ ) &&
			( t <= m.zx_ ) && ( t <= m.zy_ ) && ( t <= m.zz_ );
	}


	/// @brief Value >= xyzMatrix
	friend
	inline
	bool
	operator >=( Value const & t, xyzMatrix const & m )
	{
		return
			( t >= m.xx_ ) && ( t >= m.xy_ ) && ( t >= m.xz_ ) &&
			( t >= m.yx_ ) && ( t >= m.yy_ ) && ( t >= m.yz_ ) &&
			( t >= m.zx_ ) && ( t >= m.zy_ ) && ( t >= m.zz_ );
	}


	/// @brief Value > xyzMatrix
	friend
	inline
	bool
	operator >( Value const & t, xyzMatrix const & m )
	{
		return
			( t > m.xx_ ) && ( t > m.xy_ ) && ( t > m.xz_ ) &&
			( t > m.yx_ ) && ( t > m.yy_ ) && ( t > m.yz_ ) &&
			( t > m.zx_ ) && ( t > m.zy_ ) && ( t > m.zz_ );
	}

	/// @brief Show
	inline
	void
	show(std::ostream & output=std::cout) const
	{
		output << "(" << xx_ << ", " << xy_ << ", " << xz_ << ")" << std::endl;
		output << "(" << yx_ << ", " << yy_ << ", " << yz_ << ")" << std::endl;
		output << "(" << zx_ << ", " << zy_ << ", " << zz_ << ")" << std::endl;
	}


private: // Fields


	/// @brief Elements of the 3x3 matrix
	Value xx_, xy_, xz_;
	Value yx_, yy_, yz_;
	Value zx_, zy_, zz_;


}; // xyzMatrix


/// @brief xyzMatrix + xyzMatrix
template< typename T >
xyzMatrix< T >
operator +( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix + T
template< typename T >
xyzMatrix< T >
operator +( xyzMatrix< T > const & m, T const & t );


/// @brief T + xyzMatrix
template< typename T >
xyzMatrix< T >
operator +( T const & t, xyzMatrix< T > const & m );


/// @brief xyzMatrix - xyzMatrix
template< typename T >
xyzMatrix< T >
operator -( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix - T
template< typename T >
xyzMatrix< T >
operator -( xyzMatrix< T > const & m, T const & t );


/// @brief T - xyzMatrix
template< typename T >
xyzMatrix< T >
operator -( T const & t, xyzMatrix< T > const & m );


/// @brief xyzMatrix * xyzMatrix
template< typename T >
xyzMatrix< T >
operator *( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix * T
template< typename T >
xyzMatrix< T >
operator *( xyzMatrix< T > const & m, T const & t );


/// @brief T * xyzMatrix
template< typename T >
xyzMatrix< T >
operator *( T const & t, xyzMatrix< T > const & m );


/// @brief xyzMatrix / T
template< typename T >
xyzMatrix< T >
operator /( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix == xyzMatrix
template< typename T >
bool
operator ==( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix != xyzMatrix
template< typename T >
bool
operator !=( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix < xyzMatrix
template< typename T >
bool
operator <( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix <= xyzMatrix
template< typename T >
bool
operator <=( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix >= xyzMatrix
template< typename T >
bool
operator >=( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix > xyzMatrix
template< typename T >
bool
operator >( xyzMatrix< T > const & a, xyzMatrix< T > const & b );


/// @brief xyzMatrix == T
template< typename T >
bool
operator ==( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix != T
template< typename T >
bool
operator !=( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix < T
template< typename T >
bool
operator <( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix <= T
template< typename T >
bool
operator <=( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix >= T
template< typename T >
bool
operator >=( xyzMatrix< T > const & m, T const & t );


/// @brief xyzMatrix > T
template< typename T >
bool
operator >( xyzMatrix< T > const & m, T const & t );


/// @brief T == xyzMatrix
template< typename T >
bool
operator ==( T const & t, xyzMatrix< T > const & m );


/// @brief T != xyzMatrix
template< typename T >
bool
operator !=( T const & t, xyzMatrix< T > const & m );


/// @brief T < xyzMatrix
template< typename T >
bool
operator <( T const & t, xyzMatrix< T > const & m );


/// @brief T <= xyzMatrix
template< typename T >
bool
operator <=( T const & t, xyzMatrix< T > const & m );


/// @brief T >= xyzMatrix
template< typename T >
bool
operator >=( T const & t, xyzMatrix< T > const & m );


/// @brief T > xyzMatrix
template< typename T >
bool
operator >( T const & t, xyzMatrix< T > const & m );


/// @brief xyzMatrix * xyzVector
/// @note  Same as product( xyzMatrix, xyzVector )
// NOTE: this moved to be a member operator instead so it could work in PyRosetta
// template< typename T >
// inline
// xyzVector< T >
// operator *( xyzMatrix< T > const & m, xyzVector< T > const & v )
// {
// 	return xyzVector< T >(
// 		m.xx_ * v.x_ + m.xy_ * v.y_ + m.xz_ * v.z_,
// 		m.yx_ * v.x_ + m.yy_ * v.y_ + m.yz_ * v.z_,
// 		m.zx_ * v.x_ + m.zy_ * v.y_ + m.zz_ * v.z_
// 	);
// }


// PyRosetta work around for templates classes
//class xyzMatrix_Double : public xyzMatrix< double >
//{};

} // namespace numeric


#endif // INCLUDED_numeric_xyzMatrix_HH
