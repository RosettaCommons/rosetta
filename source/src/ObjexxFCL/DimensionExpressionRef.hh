#ifndef INCLUDED_ObjexxFCL_DimensionExpressionRef_hh
#define INCLUDED_ObjexxFCL_DimensionExpressionRef_hh


// DimensionExpressionRef: Dimension Reference DimensionExpression
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
#include <ObjexxFCL/DimensionExpression.hh>

// C++ Headers
#include <cassert>


namespace ObjexxFCL {


/// @brief DimensionExpressionRef: Dimension Reference DimensionExpression
class DimensionExpressionRef :
	public DimensionExpression
{


private: // Types


	typedef  DimensionExpression  Super;


public: // Creation


	/// @brief Copy Constructor
	inline
	DimensionExpressionRef( DimensionExpressionRef const & exp ) :
		Super(),
		dim_p_( exp.dim_p_ )
	{
		assert( dim_p_ );
	}


	/// @brief Dimension Constructor
	inline
	explicit
	DimensionExpressionRef( Dimension const & dim ) :
		dim_p_( &dim )
	{}


	/// @brief Clone
	inline
	DimensionExpressionRef *
	clone() const
	{
		return new DimensionExpressionRef( *this );
	}


	/// @brief Clone with Dimension Substitution
	DimensionExpression *
	clone( Dimension const & dim ) const;


	/// @brief Destructor
	inline
	virtual
	~DimensionExpressionRef()
	{}


public: // Inspector


	/// @brief Initialized?
	bool
	initialized() const;


	/// @brief Integer?
	inline
	bool
	integer() const
	{
		return true;
	}


	/// @brief Constant?
	inline
	bool
	constant() const
	{
		return false;
	}


	/// @brief Reference?
	inline
	bool
	reference() const
	{
		return true;
	}


	/// @brief Reducible?
	inline
	bool
	reducible() const
	{
		return false;
	}


	/// @brief Value
	double
	operator ()() const;


	/// @brief Value
	double
	value() const;


	/// @brief Insert an Observer
	void
	insert_observer( Observer & observer ) const;


	/// @brief Remove an Observer
	void
	remove_observer( Observer & observer ) const;


public: // Modifier


	/// @brief Update for Destruction of a Subject
	void
	destructed( Subject const & subject );


private: // Data


	/// @brief Pointer (non-owning) to Dimension referenced
	Dimension const * dim_p_;


}; // DimensionExpressionRef


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_DimensionExpressionRef_HH
