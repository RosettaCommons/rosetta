#ifndef INCLUDED_ObjexxFCL_DimensionExpressionCube_hh
#define INCLUDED_ObjexxFCL_DimensionExpressionCube_hh


// DimensionExpressionCube: DimensionExpression Cube Function
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
#include <ObjexxFCL/DimensionExpressionCon.hh>


namespace ObjexxFCL {


/// @brief DimensionExpressionCube: DimensionExpression Cube Function
class DimensionExpressionCube :
	public DimensionExpression
{


private: // Types


	typedef  DimensionExpression  Super;


public: // Creation


	/// @brief Copy Constructor
	inline
	DimensionExpressionCube( DimensionExpressionCube const & exp ) :
		Super(),
		exp_p_( exp.exp_p_ ? exp.exp_p_->clone() : static_cast< DimensionExpression * >( 0 ) )
	{
		assert( exp_p_ );
	}


	/// @brief Expression Constructor
	inline
	DimensionExpressionCube( DimensionExpression const & exp ) :
		exp_p_( exp.clone() )
	{
		assert( exp_p_ );
	}


	/// @brief Expression Pointer Constructor (Ownership Transfer)
	inline
	DimensionExpressionCube( DimensionExpression * exp_p_a ) :
		exp_p_( exp_p_a )
	{
		assert( exp_p_ );
	}


	/// @brief Clone
	inline
	DimensionExpression *
	clone() const
	{
		assert( exp_p_ );
		if ( constant() ) {
			if ( integer() ) {
				return new DimensionExpressionCon( exp_p_->ivalue() * exp_p_->ivalue() * exp_p_->ivalue() );
			} else {
				return new DimensionExpressionCon( exp_p_->value() * exp_p_->value() * exp_p_->value() );
			}
		} else {
			return new DimensionExpressionCube( exp_p_->clone() );
		}
	}


	/// @brief Clone with Dimension Substitution
	inline
	DimensionExpression *
	clone( Dimension const & dim ) const
	{
		assert( exp_p_ );
		if ( constant() ) {
			if ( integer() ) {
				return new DimensionExpressionCon( exp_p_->ivalue() * exp_p_->ivalue() * exp_p_->ivalue() );
			} else {
				return new DimensionExpressionCon( exp_p_->value() * exp_p_->value() * exp_p_->value() );
			}
		} else {
			return new DimensionExpressionCube( exp_p_->clone( dim ) );
		}
	}


	/// @brief Destructor
	inline
	virtual
	~DimensionExpressionCube()
	{
		assert( exp_p_ );
		delete exp_p_;
	}


public: // Inspector


	/// @brief Initialized?
	inline
	bool
	initialized() const
	{
		assert( exp_p_ );
		return ( exp_p_->initialized() );
	}


	/// @brief Integer?
	inline
	bool
	integer() const
	{
		assert( exp_p_ );
		return ( exp_p_->integer() );
	}


	/// @brief Constant?
	inline
	bool
	constant() const
	{
		assert( exp_p_ );
		return ( exp_p_->constant() );
	}


	/// @brief Reference?
	inline
	bool
	reference() const
	{
		assert( exp_p_ );
		return ( exp_p_->reference() );
	}


	/// @brief Reducible?
	inline
	bool
	reducible() const
	{
		assert( exp_p_ );
		return ( ( constant() ) || ( exp_p_->reducible() ) );
	}


	/// @brief Value
	inline
	double
	operator ()() const
	{
		assert( exp_p_ );
		return ( exp_p_->operator ()() * exp_p_->operator ()() * exp_p_->operator ()() );
	}


	/// @brief Value
	inline
	double
	value() const
	{
		assert( exp_p_ );
		return ( exp_p_->value() * exp_p_->value() * exp_p_->value() );
	}


	/// @brief Insert an Observer
	inline
	void
	insert_observer( Observer & observer ) const
	{
		assert( exp_p_ );
		exp_p_->insert_observer( observer );
	}


	/// @brief Remove an Observer
	inline
	void
	remove_observer( Observer & observer ) const
	{
		assert( exp_p_ );
		exp_p_->remove_observer( observer );
	}


public: // Modifier


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & subject )
	{
		assert( exp_p_ );
		exp_p_->destructed( subject );
	}


private: // Data


	/// @brief Pointer to expression
	DimensionExpression * exp_p_;


}; // DimensionExpressionCube


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_DimensionExpressionCube_HH
