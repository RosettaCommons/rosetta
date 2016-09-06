#ifndef INCLUDED_ObjexxFCL_DynamicIndexRange_hh
#define INCLUDED_ObjexxFCL_DynamicIndexRange_hh


// DynamicIndexRange: Dynamic Index Range
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
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/Dimension.hh>

// C++ Headers
//#include <type_traits> // for swap
#include <algorithm> // for swap
#include <utility> // also for swap??

namespace ObjexxFCL {


/// @brief DynamicIndexRange: Dynamic Index Range
///
/// @remarks
///  @li Initialized unless an active Dimension is uninitialized
///  @li Uninitialized range has ( size == 0 ) and its range values should not be accessed
///  @li Zero-size range is indicated by ( l - 1 == u ) and ( size == 0 )
///  @li Upper-unbounded range is indicated by ( l - 2 == u ) and ( size == npos )
///  @li Legal ranges have ( l - 2 <= u ) with l and u in their allowed ranges
class DynamicIndexRange :
	public IndexRange,
	public ObserverSingle // DynamicIndexRange only observed at most by the FArray that contains it
{


private: // Types


	typedef  IndexRange  Super;


public: // Types


	using Super::l;
	using Super::u;

	typedef  DimensionExpression  Expression;


public: // Creation


	/// @brief Default Constructor
	inline
	DynamicIndexRange() :
		l_dim_p_( 0 ),
		u_dim_p_( 0 )
	{}


	/// @brief Copy Constructor
	inline
	DynamicIndexRange( DynamicIndexRange const & I ) :
		Super( I ),
		ObserverSingle( I ),
		l_dim_p_( I.l_dim_clone() ),
		u_dim_p_( I.u_dim_clone() )
	{
		assert( legal_dynamic() );
		insert_as_observer();
	}


	/// @brief IndexRange Constructor
	inline
	DynamicIndexRange( IndexRange const & I ) :
		Super( I ),
		l_dim_p_( I.l_dim_clone() ),
		u_dim_p_( I.u_dim_clone() )
	{
		assert( legal_dynamic() );
		insert_as_observer();
	}


	/// @brief Upper Index Constructor
	inline
	DynamicIndexRange( int const u_a ) :
		Super( u_a ),
		l_dim_p_( 0 ),
		u_dim_p_( 0 )
	{
		assert( legal_static() );
	}


	/// @brief Unbounded Upper Index Constructor
	inline
	DynamicIndexRange( Star const & star ) :
		Super( star ),
		l_dim_p_( 0 ),
		u_dim_p_( 0 )
	{}


	/// @brief Upper Dimension Constructor
	inline
	DynamicIndexRange( Dimension const & u_dim_a ) :
		Super( u_dim_a.zvalue() ),
		l_dim_p_( 0 ),
		u_dim_p_( u_dim_a.reference_copy() )
	{
		assert( legal_dynamic() );
		size_dynamic();
		u_insert_as_observer();
	}


	/// @brief Upper Expression Constructor
	inline
	DynamicIndexRange( Expression const & u_exp_a ) :
		Super( u_exp_a.zvalue() ),
		l_dim_p_( 0 ),
		u_dim_p_( new Dimension( u_exp_a ) )
	{
		assert( legal_dynamic() );
		size_dynamic();
		u_insert_as_observer();
	}


	/// @brief Index Range Constructor
	inline
	DynamicIndexRange( int const l_a, int const u_a ) :
		Super( l_a, u_a ),
		l_dim_p_( 0 ),
		u_dim_p_( 0 )
	{
		assert( legal_static() );
	}


	/// @brief Dimension Range Constructor
	inline
	DynamicIndexRange( Dimension const & l_dim_a, Dimension const & u_dim_a ) :
		Super( l_dim_a.zvalue(), u_dim_a.zvalue() ),
		l_dim_p_( l_dim_a.reference_copy() ),
		u_dim_p_( u_dim_a.reference_copy() )
	{
		assert( legal_dynamic() );
		size_dynamic();
		insert_as_observer();
	}


	/// @brief Expression Range Constructor
	inline
	DynamicIndexRange( Expression const & l_exp_a, Expression const & u_exp_a ) :
		Super( l_exp_a.zvalue(), u_exp_a.zvalue() ),
		l_dim_p_( new Dimension( l_exp_a ) ),
		u_dim_p_( new Dimension( u_exp_a ) )
	{
		assert( legal_dynamic() );
		size_dynamic();
		insert_as_observer();
	}


	/// @brief Index and Dimension Constructor
	inline
	DynamicIndexRange( int const l_a, Dimension const & u_dim_a ) :
		Super( l_a, u_dim_a.zvalue() ),
		l_dim_p_( 0 ),
		u_dim_p_( u_dim_a.reference_copy() )
	{
		assert( legal_dynamic() );
		size_dynamic();
		u_insert_as_observer();
	}


	/// @brief Dimension and Index Constructor
	inline
	DynamicIndexRange( Dimension const & l_dim_a, int const u_a ) :
		Super( l_dim_a.zvalue(), u_a ),
		l_dim_p_( l_dim_a.reference_copy() ),
		u_dim_p_( 0 )
	{
		assert( legal_dynamic() );
		size_dynamic();
		l_insert_as_observer();
	}


	/// @brief Index and Expression Constructor
	inline
	DynamicIndexRange( int const l_a, Expression const & u_exp_a ) :
		Super( l_a, u_exp_a.zvalue() ),
		l_dim_p_( 0 ),
		u_dim_p_( new Dimension( u_exp_a ) )
	{
		assert( legal_dynamic() );
		size_dynamic();
		u_insert_as_observer();
	}


	/// @brief Expression and Index Constructor
	inline
	DynamicIndexRange( Expression const & l_exp_a, int const u_a ) :
		Super( l_exp_a.zvalue(), u_a ),
		l_dim_p_( new Dimension( l_exp_a ) ),
		u_dim_p_( 0 )
	{
		assert( legal_dynamic() );
		size_dynamic();
		l_insert_as_observer();
	}


	/// @brief Dimension and Expression Constructor
	inline
	DynamicIndexRange( Dimension const & l_dim_a, Expression const & u_exp_a ) :
		Super( l_dim_a.zvalue(), u_exp_a.zvalue() ),
		l_dim_p_( l_dim_a.reference_copy() ),
		u_dim_p_( new Dimension( u_exp_a ) )
	{
		assert( legal_dynamic() );
		size_dynamic();
		insert_as_observer();
	}


	/// @brief Expression and Dimension Constructor
	inline
	DynamicIndexRange( Expression const & l_exp_a, Dimension const & u_dim_a ) :
		Super( l_exp_a.zvalue(), u_dim_a.zvalue() ),
		l_dim_p_( new Dimension( l_exp_a ) ),
		u_dim_p_( u_dim_a.reference_copy() )
	{
		assert( legal_dynamic() );
		size_dynamic();
		insert_as_observer();
	}


	/// @brief Index and Unbounded Upper Index Constructor
	inline
	DynamicIndexRange( int const l_a, Star const & star ) :
		Super( l_a, star ),
		l_dim_p_( 0 ),
		u_dim_p_( 0 )
	{
		assert( legal_static() );
	}


	/// @brief Dimension and Unbounded Upper Index Constructor
	DynamicIndexRange( Dimension const & l_dim_a, Star const & star );


	/// @brief Expression and Unbounded Upper Index Constructor
	DynamicIndexRange( Expression const & l_exp_a, Star const & star );


	/// @brief Destructor
	inline
	virtual
	~DynamicIndexRange()
	{
		delete l_dim_p_;
		delete u_dim_p_;
	}


public: // Assignment


	/// @brief Copy Assignment
	inline
	DynamicIndexRange &
	operator =( DynamicIndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		notify();
		return *this;
	}


	/// @brief IndexRange Assignment
	inline
	DynamicIndexRange &
	operator =( IndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		notify();
		return *this;
	}


	/// @brief Upper Index Assignment
	inline
	DynamicIndexRange &
	operator =( int const u_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::operator =( u_a );
		assert( legal_static() );
		notify();
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	DynamicIndexRange &
	operator =( Star const & star )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::operator =( star );
		notify();
		return *this;
	}


	/// @brief Upper Dimension Assignment
	inline
	DynamicIndexRange &
	operator =( Dimension const & u_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::operator =( u_dim_a.zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Upper Expression Assignment
	inline
	DynamicIndexRange &
	operator =( Expression const & u_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::operator =( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief DynamicIndexRange Assignment
	inline
	DynamicIndexRange &
	assign( DynamicIndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		notify();
		return *this;
	}


	/// @brief IndexRange Assignment
	inline
	DynamicIndexRange &
	assign( IndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		notify();
		return *this;
	}


	/// @brief Upper Index Assignment
	inline
	DynamicIndexRange &
	assign( int const u_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::operator =( u_a );
		assert( legal_static() );
		notify();
		return *this;
	}


	/// @brief Unbounded Upper Index Assignment
	inline
	DynamicIndexRange &
	assign( Star const & star )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::operator =( star );
		notify();
		return *this;
	}


	/// @brief Upper Dimension Assignment
	inline
	DynamicIndexRange &
	assign( Dimension const & u_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::operator =( u_dim_a.zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Upper Expression Assignment
	inline
	DynamicIndexRange &
	assign( Expression const & u_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::operator =( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Index Range Assignment
	inline
	DynamicIndexRange &
	assign( int const l_a, int const u_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::assign( l_a, u_a );
		assert( legal_static() );
		notify();
		return *this;
	}


	/// @brief Dimension Range Assignment
	inline
	DynamicIndexRange &
	assign( Dimension const & l_dim_a, Dimension const & u_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::assign( l_dim_a.zvalue(), u_dim_a.zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Expression Range Assignment
	inline
	DynamicIndexRange &
	assign( Expression const & l_exp_a, Expression const & u_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::assign( l_dim_p_->zvalue(), u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Index and Dimension Assignment
	inline
	DynamicIndexRange &
	assign( int const l_a, Dimension const & u_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::assign( l_a, u_dim_a.zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Dimension and Index Assignment
	inline
	DynamicIndexRange &
	assign( Dimension const & l_dim_a, int const u_a )
	{
		delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = 0;
		Super::assign( l_dim_a.zvalue(), u_a );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Index and Expression Assignment
	inline
	DynamicIndexRange &
	assign( int const l_a, Expression const & u_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::assign( l_a, u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Expression and Index Assignment
	inline
	DynamicIndexRange &
	assign( Expression const & l_exp_a, int const u_a )
	{
		delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = 0;
		Super::assign( l_dim_p_->zvalue(), u_a );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Dimension and Expression Assignment
	inline
	DynamicIndexRange &
	assign( Dimension const & l_dim_a, Expression const & u_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::assign( l_dim_a.zvalue(), u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Expression and Dimension Assignment
	inline
	DynamicIndexRange &
	assign( Expression const & l_exp_a, Dimension const & u_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::assign( l_dim_p_->zvalue(), u_dim_a.zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Index and Unbounded Upper Index Assignment
	inline
	DynamicIndexRange &
	assign( int const l_a, Star const & star )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::assign( l_a, star );
		assert( legal_static() );
		notify();
		return *this;
	}


	/// @brief Dimension and Unbounded Upper Index Assignment
	DynamicIndexRange &
	assign( Dimension const & l_dim_a, Star const & star );


	/// @brief Expression and Unbounded Upper Index Assignment
	DynamicIndexRange &
	assign( Expression const & l_exp_a, Star const & star );


	/// @brief DynamicIndexRange Assignment Without Notification
	inline
	DynamicIndexRange &
	assign_no_notify( DynamicIndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		return *this;
	}


	/// @brief IndexRange Assignment Without Notification
	inline
	DynamicIndexRange &
	assign_no_notify( IndexRange const & I )
	{
		if ( this != &I ) {
			delete l_dim_p_; l_dim_p_ = I.l_dim_clone(); l_insert_as_observer();
			delete u_dim_p_; u_dim_p_ = I.u_dim_clone(); u_insert_as_observer();
			Super::operator =( I );
			assert( legal_dynamic() );
		}
		return *this;
	}


	/// @brief Index and Unbounded Upper Index Assignment Without Notification
	inline
	DynamicIndexRange &
	assign_no_notify( int const l_a, Star const & star )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		Super::assign( l_a, star );
		assert( legal_static() );
		return *this;
	}


public: // Predicate


	/// @brief Initialized?
	inline
	bool
	initialized() const
	{
		return ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	}


	/// @brief Lower Initialized?
	inline
	bool
	l_initialized() const
	{
		return ( l_dim_p_ ? l_dim_p_->initialized_ : true );
	}


	/// @brief Upper Initialized?
	inline
	bool
	u_initialized() const
	{
		return ( u_dim_p_ ? u_dim_p_->initialized_ : true );
	}


	/// @brief Legal?
	inline
	bool
	legal() const
	{
		return ( ( ( l_ >= l_min ) && ( u_ <= u_max ) && ( l_ - 2 <= u_ ) ) ||
		 ( ! ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) ) );
	}


	/// @brief Bounded?
	inline
	bool
	bounded() const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( Super::bounded() ) );
	}


	/// @brief Unbounded?
	inline
	bool
	unbounded() const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( Super::unbounded() ) );
	}


	/// @brief Not Unbounded?
	inline
	bool
	not_unbounded() const
	{
		return ( ( ! ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) ) || ( Super::not_unbounded() ) );
	}


	/// @brief Bounded with Positive Size?
	inline
	bool
	positive() const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( Super::positive() ) );
	}


	/// @brief Contains an Index?
	inline
	bool
	contains( int const i ) const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( ( l_ <= i ) && ( ( i <= u_ ) || ( size_ == npos ) ) ) );
	}


	/// @brief Contains Another IndexRange?
	inline
	bool
	contains( IndexRange const & I ) const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( Super::contains( I ) ) );
	}


	/// @brief Intersects Another IndexRange?
	inline
	bool
	intersects( IndexRange const & I ) const
	{
		return ( ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) && ( Super::intersects( I ) ) );
	}


public: // Inspector


public: // Modifier


	/// @brief Lower Index Set
	inline
	DynamicIndexRange &
	l( int const l_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		Super::l( l_a );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Lower Dimension Set
	inline
	DynamicIndexRange &
	l( Dimension const & l_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
		Super::l( l_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Lower Expression Set
	inline
	DynamicIndexRange &
	l( Expression const & l_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
		Super::l( l_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Lower Index Set Without Notification
	inline
	DynamicIndexRange &
	l_no_notify( int const l_a )
	{
		delete l_dim_p_; l_dim_p_ = 0;
		Super::l( l_a );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Lower Dimension Set Without Notification
	inline
	DynamicIndexRange &
	l_no_notify( Dimension const & l_dim_a )
	{
		delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
		Super::l( l_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Lower Expression Set Without Notification
	inline
	DynamicIndexRange &
	l_no_notify( Expression const & l_exp_a )
	{
		delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
		Super::l( l_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Upper Index Set
	inline
	DynamicIndexRange &
	u( int const u_a )
	{
		delete u_dim_p_; u_dim_p_ = 0;
		Super::u( u_a );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Unbounded Upper Index Set
	DynamicIndexRange &
	u( Star const & star );


	/// @brief Upper Dimension Set
	inline
	DynamicIndexRange &
	u( Dimension const & u_dim_a )
	{
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::u( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Upper Expression Set
	inline
	DynamicIndexRange &
	u( Expression const & u_exp_a )
	{
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::u( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		notify();
		return *this;
	}


	/// @brief Upper Index Set Without Notification
	inline
	DynamicIndexRange &
	u_no_notify( int const u_a )
	{
		delete u_dim_p_; u_dim_p_ = 0;
		Super::u( u_a );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Unbounded Upper Index Set Without Notification
	DynamicIndexRange &
	u_no_notify( Star const & star );


	/// @brief Upper Dimension Set Without Notification
	inline
	DynamicIndexRange &
	u_no_notify( Dimension const & u_dim_a )
	{
		delete u_dim_p_; u_dim_p_ = u_dim_a.reference_copy(); u_insert_as_observer();
		Super::u( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Upper Expression Set Without Notification
	inline
	DynamicIndexRange &
	u_no_notify( Expression const & u_exp_a )
	{
		delete u_dim_p_; u_dim_p_ = new Dimension( u_exp_a ); u_insert_as_observer();
		Super::u( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
		return *this;
	}


	/// @brief Expand to Contain an Index
	DynamicIndexRange &
	contain( int const i );


	/// @brief Expand to Contain an Index and Notify If Changed
	DynamicIndexRange &
	contain_nic( int const i );


	/// @brief Expand to Contain Another IndexRange
	DynamicIndexRange &
	contain( IndexRange const & I );


	/// @brief Expand to Contain Another IndexRange and Notify If Changed
	DynamicIndexRange &
	contain_nic( IndexRange const & I );


	/// @brief Intersect With Another IndexRange
	DynamicIndexRange &
	intersect( IndexRange const & I );


	/// @brief Intersect With Another IndexRange and Notify If Changed
	DynamicIndexRange &
	intersect_nic( IndexRange const & I );


	/// @brief Clear
	inline
	DynamicIndexRange &
	clear()
	{
		Super::clear();
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		notify();
		return *this;
	}


	/// @brief Clear Without Notification
	inline
	DynamicIndexRange &
	clear_no_notify()
	{
		Super::clear();
		delete l_dim_p_; l_dim_p_ = 0;
		delete u_dim_p_; u_dim_p_ = 0;
		return *this;
	}


	/// @brief Swap
	inline
	DynamicIndexRange &
	swap( DynamicIndexRange & I )
	{
		if ( this != &I ) {
			remove_as_observer();
			I.remove_as_observer();
			std::swap( l_dim_p_, I.l_dim_p_ );
			std::swap( u_dim_p_, I.u_dim_p_ );
			insert_as_observer();
			I.insert_as_observer();
			Super::swap( I );
			assert( legal_dynamic() );
			notify();
		}
		return *this;
	}


	/// @brief Swap Without Notification
	inline
	DynamicIndexRange &
	swap_no_notify( DynamicIndexRange & I )
	{
		if ( this != &I ) {
			remove_as_observer();
			I.remove_as_observer();
			std::swap( l_dim_p_, I.l_dim_p_ );
			std::swap( u_dim_p_, I.u_dim_p_ );
			insert_as_observer();
			I.insert_as_observer();
			Super::swap( I );
			assert( legal_dynamic() );
		}
		return *this;
	}


public: // Observer Modifier


	/// @brief Update
	inline
	void
	update()
	{
		if ( l_dim_p_ ) Super::l( l_dim_p_->zvalue() );
		if ( u_dim_p_ ) Super::u( u_dim_p_->zvalue() );
		assert( legal_dynamic() );
		size_dynamic();
	}


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & )
	{}


public: // Friend


	/// @brief Swap
	friend
	inline
	void
	swap( DynamicIndexRange & a, DynamicIndexRange & b )
	{
		a.swap( b );
	}


	/// @brief Swap
	friend
	inline
	void
	swap_no_notify( DynamicIndexRange & a, DynamicIndexRange & b )
	{
		a.swap_no_notify( b );
	}


private: // Functions


	/// @brief Legal DynamicIndexRange?
	inline
	bool
	legal_dynamic() const
	{
		return ( ( ( l_ >= l_min ) && ( u_ <= u_max ) && ( l_ - 2 <= u_ ) ) ||
		 ( ! ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) ) );
	}


	/// @brief Set Size to Zero if Uninitialized
	inline
	void
	size_dynamic()
	{
		if ( ! ( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) ) ) size_ = 0;
	}


	/// @brief Lower Dimension Clone
	inline
	Dimension *
	l_dim_clone() const
	{
		return ( l_dim_p_ ? l_dim_p_->clone() : static_cast< Dimension * >( 0 ) );
	}


	/// @brief Upper Dimension Clone
	inline
	Dimension *
	u_dim_clone() const
	{
		return ( u_dim_p_ ? u_dim_p_->clone() : static_cast< Dimension * >( 0 ) );
	}


	/// @brief Insert as Observer of the Dimensions
	inline
	void
	insert_as_observer()
	{
		if ( l_dim_p_ ) l_dim_p_->insert_observer( *this );
		if ( u_dim_p_ ) u_dim_p_->insert_observer( *this );
	}


	/// @brief Insert as Observer of the Lower Dimension
	inline
	void
	l_insert_as_observer()
	{
		if ( l_dim_p_ ) l_dim_p_->insert_observer( *this );
	}


	/// @brief Insert as Observer of the Upper Dimension
	inline
	void
	u_insert_as_observer()
	{
		if ( u_dim_p_ ) u_dim_p_->insert_observer( *this );
	}


	/// @brief Remove as Observer of the Dimensions
	inline
	void
	remove_as_observer()
	{
		if ( l_dim_p_ ) l_dim_p_->remove_observer( *this );
		if ( u_dim_p_ ) u_dim_p_->remove_observer( *this );
	}


private: // Data


	/// @brief Lower Dimension pointer (0 iff no Dimension)
	Dimension * l_dim_p_;

	/// @brief Upper Dimension pointer (0 iff no Dimension)
	Dimension * u_dim_p_;


}; // DynamicIndexRange


/// @brief Swap
void
swap( DynamicIndexRange & a, DynamicIndexRange & b );


/// @brief Swap
void
swap_no_notify( DynamicIndexRange & a, DynamicIndexRange & b );


} // namespace ObjexxFCL


#ifndef NO_STD_SWAP_OVERLOADS


// std::swap Overloads for Efficiency
//
// Technically you cannot add functions overloads to namespace std
// but this works with most compilers and makes it much faster if someone uses
// std::swap instead of swap or ObjexxFCL::swap.


namespace std {


/// @brief std::swap( DynamicIndexRange, DynamicIndexRange )
inline
void
swap( ObjexxFCL::DynamicIndexRange & a, ObjexxFCL::DynamicIndexRange & b )
{
	a.swap( b );
}


} // namespace std


#endif // NO_STD_SWAP_OVERLOADS


#endif // INCLUDED_ObjexxFCL_DynamicIndexRange_HH
