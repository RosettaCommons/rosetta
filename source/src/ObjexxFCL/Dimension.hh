#ifndef INCLUDED_ObjexxFCL_Dimension_hh
#define INCLUDED_ObjexxFCL_Dimension_hh


// Dimension: Dynamic Dimension
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
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DimensionExpression.hh>

// C++ Headers
#include <cassert>
#include <iosfwd>

#ifdef CXX11
#include <utility>
#include <type_traits>  // for swap
#else
#include <algorithm>
#endif

namespace ObjexxFCL {


/// @brief Dimension: Dynamic Dimension
class Dimension :
	public ObserverMulti
{


private: // Friend


	friend class DynamicIndexRange;


public: // Types


	typedef  DimensionExpression  Expression;


public: // Creation


	/// @brief Default Constructor
	inline
	Dimension() :
		exp_p_( 0 ),
		initialized_( false ),
		value_( 0 )
	{}


	/// @brief Copy Constructor
	/// @note  Copy creates a reference to the passed Dimension: Not intended for pass-by-value
	explicit
	Dimension( Dimension const & dim );


	/// @brief int Constructor
	explicit
	Dimension( int const i );


	/// @brief double Constructor
	explicit
	Dimension( double const d );


	/// @brief Expression Constructor
	inline
	explicit
	Dimension( Expression const & exp ) :
		exp_p_( exp.clone() ),
		initialized_( exp_p_->initialized() ),
		value_( initialized_ ? exp_p_->ivalue() : 0 )
	{
		insert_as_observer();
	}


	/// @brief Expression Pointer Constructor (Ownership Transfer)
	inline
	explicit
	Dimension( Expression * exp_p_a ) :
		exp_p_( exp_p_a ),
		initialized_( exp_p_ ? exp_p_->initialized() : false ),
		value_( initialized_ ? exp_p_->ivalue() : 0 )
	{
		reduce_expression();
		insert_as_observer();
	}


	/// @brief Clone
	inline
	Dimension *
	clone() const
	{
		return new Dimension( exp_clone() );
	}


	/// @brief Reference Copy
	inline
	Dimension *
	reference_copy() const
	{
		return new Dimension( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~Dimension()
	{
		remove_as_observer();
		delete exp_p_;
	}


public: // Conversion


	/// @brief int Conversion
	inline
	operator int() const
	{
		assert( initialized_ );
		return value_;
	}


	/// @brief double Conversion
	inline
	operator double() const
	{
		assert( initialized_ );
		return static_cast< double >( value_ );
	}


public: // Assignment


	/// @brief Copy Assignment: Creates a reference to the passed Dimension
	Dimension &
	operator =( Dimension const & dim );


	/// @brief Expression Assignment
	inline
	Dimension &
	operator =( Expression const & exp )
	{
		assert( exp_p_ != &exp );
		remove_as_observer();
		delete exp_p_; exp_p_ = exp.clone( *this );
		insert_as_observer();
		update_notify();
		return *this;
	}


	/// @brief int Assignment
	Dimension &
	operator =( int const i );


	/// @brief double Assignment
	Dimension &
	operator =( double const d );


	/// @brief += Dimension
	Dimension &
	operator +=( Dimension const & dim );


	/// @brief += Expression
	Dimension &
	operator +=( Expression const & exp );


	/// @brief += int
	Dimension &
	operator +=( int const i );


	/// @brief += double
	Dimension &
	operator +=( double const d );


	/// @brief -= Dimension
	Dimension &
	operator -=( Dimension const & dim );


	/// @brief -= Expression
	Dimension &
	operator -=( Expression const & exp );


	/// @brief -= int
	Dimension &
	operator -=( int const i );


	/// @brief -= double
	Dimension &
	operator -=( double const d );


	/// @brief *= Dimension
	Dimension &
	operator *=( Dimension const & dim );


	/// @brief *= Expression
	Dimension &
	operator *=( Expression const & exp );


	/// @brief *= int
	Dimension &
	operator *=( int const i );


	/// @brief *= double
	Dimension &
	operator *=( double const d );


	/// @brief /= Dimension
	Dimension &
	operator /=( Dimension const & dim );


	/// @brief /= Expression
	Dimension &
	operator /=( Expression const & exp );


	/// @brief /= int
	Dimension &
	operator /=( int const i );


	/// @brief /= double
	Dimension &
	operator /=( double const d );


	/// @brief Dimension Value-Semantics Assignment
	inline
	Dimension &
	assign_value_of( Dimension const & dim )
	{
		if ( this != &dim ) {
			remove_as_observer();
			delete exp_p_; exp_p_ = dim.exp_clone();
			insert_as_observer();
			update();
		}
		notify();
		return *this;
	}


	/// @brief int Assignment if Bigger than Value or Smaller than Multiplier * Value
	Dimension &
	assign_if( int const i, double const m = 1.0 );


	/// @brief double Assignment if Bigger than Value or Smaller than Multiplier * Value
	Dimension &
	assign_if( double const d, double const m = 1.0 );


	/// @brief int Assignment if Bigger than Value or Smaller than Half Value
	Dimension &
	assign_if_half( int const i );


	/// @brief double Assignment if Bigger than Value or Smaller than Half Value
	Dimension &
	assign_if_half( double const d );


	/// @brief int Assignment if Bigger than Value
	Dimension &
	assign_if_bigger( int const i );


	/// @brief double Assignment if Bigger than Value
	Dimension &
	assign_if_bigger( double const d );


	/// @brief int Assignment if Bigger than Value or Smaller than Multiplier * Value: Notify if Changed
	Dimension &
	assign_if_nic( int const i, double const m = 1.0 );


	/// @brief double Assignment if Bigger than Value or Smaller than Multiplier * Value: Notify if Changed
	Dimension &
	assign_if_nic( double const d, double const m = 1.0 );


	/// @brief int Assignment if Bigger than Value or Smaller than Half Value: Notify if Changed
	Dimension &
	assign_if_half_nic( int const i );


	/// @brief double Assignment if Bigger than Value or Smaller than Half Value: Notify if Changed
	Dimension &
	assign_if_half_nic( double const d );


	/// @brief int Assignment if Bigger than Value: Notify if Changed
	Dimension &
	assign_if_bigger_nic( int const i );


	/// @brief double Assignment if Bigger than Value: Notify if Changed
	Dimension &
	assign_if_bigger_nic( double const d );


public: // Incrememt/Decrement


	/// @brief ++Dimension
	Dimension &
	operator ++();


	/// @brief Dimension++
	Dimension const
	operator ++( int );


	/// @brief --Dimension
	Dimension &
	operator --();


	/// @brief Dimension--
	Dimension const
	operator --( int );


public: // Inspector


	/// @brief Initialized?
	inline
	bool
	initialized() const
	{
		return initialized_;
	}


	/// @brief Constant?
	inline
	bool
	constant() const
	{
		return ( exp_p_ ? exp_p_->constant() : false );
	}


	/// @brief Reference?
	inline
	bool
	reference() const
	{
		return ( exp_p_ ? exp_p_->reference() : false );
	}


	/// @brief Reducible?
	inline
	bool
	reducible() const
	{
		return ( exp_p_ ? exp_p_->reducible() : false );
	}


	/// @brief Value
	inline
	int
	operator ()() const
	{
		assert( initialized_ );
		return value_;
	}


	/// @brief Value
	inline
	int
	value() const
	{
		assert( initialized_ );
		return value_;
	}


	/// @brief Value: Zero if Uninitialized
	inline
	int
	zvalue() const
	{
		return value_;
	}


	/// @brief Expression Pointer
	inline
	Expression const *
	exp_p() const
	{
		return exp_p_;
	}


	/// @brief Expression
	inline
	Expression const &
	exp() const
	{
		assert( exp_p_ );
		return *exp_p_;
	}


	/// @brief Expression Clone
	inline
	Expression *
	exp_clone() const
	{
		return ( exp_p_ ? exp_p_->clone() : static_cast< Expression * >( 0 ) );
	}


public: // Modifier


	/// @brief Clear the Dimension
	inline
	Dimension &
	clear()
	{
		if ( exp_p_ ) {
			remove_as_observer();
			delete exp_p_; exp_p_ = 0;
			initialized_ = false;
			value_ = 0;
		}
		notify();
		return *this;
	}


	/// @brief Clear the Dimension Without Notification
	inline
	Dimension &
	clear_no_notify()
	{
		if ( exp_p_ ) {
			remove_as_observer();
			delete exp_p_; exp_p_ = 0;
			initialized_ = false;
			value_ = 0;
		}
		return *this;
	}


	/// @brief Swap
	inline
	Dimension &
	swap( Dimension & dim )
	{
		if ( this != &dim ) {
			remove_as_observer();
			dim.remove_as_observer();
			std::swap( exp_p_, dim.exp_p_ );
			std::swap( initialized_, dim.initialized_ );
			std::swap( value_, dim.value_ );
			insert_as_observer();
			dim.insert_as_observer();
			notify();
		}
		return *this;
	}


	/// @brief Swap Without Notification
	inline
	Dimension &
	swap_no_notify( Dimension & dim )
	{
		if ( this != &dim ) {
			remove_as_observer();
			dim.remove_as_observer();
			std::swap( exp_p_, dim.exp_p_ );
			std::swap( initialized_, dim.initialized_ );
			std::swap( value_, dim.value_ );
			insert_as_observer();
			dim.insert_as_observer();
		}
		return *this;
	}


public: // Observer Modifier


	/// @brief Update
	inline
	void
	update()
	{
		initialized_ = ( exp_p_ ? exp_p_->initialized() : false );
		value_ = ( initialized_ ? exp_p_->ivalue() : 0 );
	}


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & subject )
	{
		if ( exp_p_ ) exp_p_->destructed( subject );
	}


public: // Friend


	/// @brief Swap
	friend
	inline
	void
	swap( Dimension & a, Dimension & b )
	{
		a.swap( b );
	}


	/// @brief Swap
	friend
	inline
	void
	swap_no_notify( Dimension & a, Dimension & b )
	{
		a.swap_no_notify( b );
	}


private: // Functions


	/// @brief Reduce Expression
	inline
	void
	reduce_expression()
	{
		if ( reducible() ) {
			Expression * reduced_exp_p = exp_p_->clone();
			delete exp_p_; exp_p_ = reduced_exp_p;
		}
	}


	/// @brief Insert as Observer of an Expression's Referenced Dimensions
	inline
	void
	insert_as_observer_of( Dimension const & dim )
	{
		dim.insert_observer( *this );
	}


	/// @brief Insert as Observer of an Expression's Referenced Dimensions
	inline
	void
	insert_as_observer_of( Expression const & exp )
	{
		exp.insert_observer( *this );
	}


	/// @brief Insert as Observer of the Expression's Referenced Dimensions
	inline
	void
	insert_as_observer()
	{
		if ( exp_p_ ) exp_p_->insert_observer( *this );
	}


	/// @brief Remove as Observer of the Expression's Referenced Dimensions
	inline
	void
	remove_as_observer()
	{
		if ( exp_p_ ) exp_p_->remove_observer( *this );
	}


	/// @brief Update and Notify
	inline
	void
	update_notify()
	{
		initialized_ = ( exp_p_ ? exp_p_->initialized() : false );
		value_ = ( initialized_ ? exp_p_->ivalue() : 0 );
		notify();
	}


	/// @brief Update and Notify if External State Changed
	inline
	void
	update_notify_if_changed()
	{
		bool const now_initialized = ( exp_p_ ? exp_p_->initialized() : false );
		int const now_value = ( now_initialized ? exp_p_->ivalue() : 0 );
		if ( ( initialized_ != now_initialized ) || ( value_ != now_value ) ) {
			initialized_ = now_initialized;
			value_ = now_value;
			notify();
		}
	}


private: // Data


	/// @brief Expression pointer (owned)
	Expression * exp_p_;

	/// @brief Cached initialization state
	bool initialized_;

	/// @brief Cached value: Kept in synch with expression value (0 if uninitialized)
	int value_;


}; // Dimension


/// @brief Swap
void
swap( Dimension & a, Dimension & b );


/// @brief Swap
void
swap_no_notify( Dimension & a, Dimension & b );


/// @brief Dimension == Dimension
inline
bool
operator ==( Dimension const & dim1, Dimension const & dim2 )
{
	return ( ( dim1.initialized() ) && ( dim2.initialized() ) && ( dim1.value() == dim2.value() ) );
}


/// @brief Dimension != Dimension
inline
bool
operator !=( Dimension const & dim1, Dimension const & dim2 )
{
	return !( dim1 == dim2 );
}


/// @brief Dimension < Dimension
inline
bool
operator <( Dimension const & dim1, Dimension const & dim2 )
{
	return ( ( dim1.initialized() ) && ( dim2.initialized() ) && ( dim1.value() < dim2.value() ) );
}


/// @brief Dimension <= Dimension
inline
bool
operator <=( Dimension const & dim1, Dimension const & dim2 )
{
	return ( ( dim1.initialized() ) && ( dim2.initialized() ) && ( dim1.value() <= dim2.value() ) );
}


/// @brief Dimension > Dimension
inline
bool
operator >( Dimension const & dim1, Dimension const & dim2 )
{
	return ( ( dim1.initialized() ) && ( dim2.initialized() ) && ( dim1.value() > dim2.value() ) );
}


/// @brief Dimension >= Dimension
inline
bool
operator >=( Dimension const & dim1, Dimension const & dim2 )
{
	return ( ( dim1.initialized() ) && ( dim2.initialized() ) && ( dim1.value() >= dim2.value() ) );
}


/// @brief int == Dimension
inline
bool
operator ==( int const i, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( i == dim.value() ) );
}


/// @brief int != Dimension
inline
bool
operator !=( int const i, Dimension const & dim )
{
	return !( i == dim );
}


/// @brief int < Dimension
inline
bool
operator <( int const i, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( i < dim.value() ) );
}


/// @brief int <= Dimension
inline
bool
operator <=( int const i, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( i <= dim.value() ) );
}


/// @brief int > Dimension
inline
bool
operator >( int const i, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( i > dim.value() ) );
}


/// @brief int >= Dimension
inline
bool
operator >=( int const i, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( i >= dim.value() ) );
}


/// @brief Dimension == int
inline
bool
operator ==( Dimension const & dim, int const i )
{
	return ( ( dim.initialized() ) && ( dim.value() == i ) );
}


/// @brief Dimension != int
inline
bool
operator !=( Dimension const & dim, int const i )
{
	return !( dim == i );
}


/// @brief Dimension < int
inline
bool
operator <( Dimension const & dim, int const i )
{
	return ( ( dim.initialized() ) && ( dim.value() < i ) );
}


/// @brief Dimension <= int
inline
bool
operator <=( Dimension const & dim, int const i )
{
	return ( ( dim.initialized() ) && ( dim.value() <= i ) );
}


/// @brief Dimension > int
inline
bool
operator >( Dimension const & dim, int const i )
{
	return ( ( dim.initialized() ) && ( dim.value() > i ) );
}


/// @brief Dimension >= int
inline
bool
operator >=( Dimension const & dim, int const i )
{
	return ( ( dim.initialized() ) && ( dim.value() >= i ) );
}


/// @brief double == Dimension
inline
bool
operator ==( double const d, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( d == dim.value() ) );
}


/// @brief double != Dimension
inline
bool
operator !=( double const d, Dimension const & dim )
{
	return !( d == dim );
}


/// @brief double < Dimension
inline
bool
operator <( double const d, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( d < dim.value() ) );
}


/// @brief double <= Dimension
inline
bool
operator <=( double const d, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( d <= dim.value() ) );
}


/// @brief double > Dimension
inline
bool
operator >( double const d, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( d > dim.value() ) );
}


/// @brief double >= Dimension
inline
bool
operator >=( double const d, Dimension const & dim )
{
	return ( ( dim.initialized() ) && ( d >= dim.value() ) );
}


/// @brief Dimension == double
inline
bool
operator ==( Dimension const & dim, double const d )
{
	return ( ( dim.initialized() ) && ( dim.value() == d ) );
}


/// @brief Dimension != double
inline
bool
operator !=( Dimension const & dim, double const d )
{
	return !( dim == d );
}


/// @brief Dimension < double
inline
bool
operator <( Dimension const & dim, double const d )
{
	return ( ( dim.initialized() ) && ( dim.value() < d ) );
}


/// @brief Dimension <= double
inline
bool
operator <=( Dimension const & dim, double const d )
{
	return ( ( dim.initialized() ) && ( dim.value() <= d ) );
}


/// @brief Dimension > double
inline
bool
operator >( Dimension const & dim, double const d )
{
	return ( ( dim.initialized() ) && ( dim.value() > d ) );
}


/// @brief Dimension >= double
inline
bool
operator >=( Dimension const & dim, double const d )
{
	return ( ( dim.initialized() ) && ( dim.value() >= d ) );
}


/// @brief Stream Input
std::istream &
operator >>( std::istream & stream, Dimension & dim );


/// @brief Stream Output
std::ostream &
operator <<( std::ostream & stream, Dimension const & dim );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.


namespace std {


/// @brief std::swap( Dimension, Dimension )
inline
void
swap( ObjexxFCL::Dimension & a, ObjexxFCL::Dimension & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_Dimension_HH
