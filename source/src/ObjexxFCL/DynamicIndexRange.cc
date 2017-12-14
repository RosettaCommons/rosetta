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
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/DimensionExpressions.hh>


namespace ObjexxFCL {


// DynamicIndexRange: Dynamic Index Range


/// @brief Dimension and Unbounded Upper Index Constructor
DynamicIndexRange::DynamicIndexRange( Dimension const & l_dim_a, Star const & star ) :
	IndexRange( l_dim_a.zvalue(), star ),
	l_dim_p_( l_dim_a.reference_copy() ),
	u_dim_p_( new Dimension( *l_dim_p_ - 2 ) ) // Stays unbounded until upper index changed
{
	assert( legal_dynamic() );
	size_dynamic();
	insert_as_observer();
}


/// @brief Expression and Unbounded Upper Index Constructor
DynamicIndexRange::DynamicIndexRange( Expression const & l_exp_a, Star const & star ) :
	IndexRange( l_exp_a.zvalue(), star ),
	l_dim_p_( new Dimension( l_exp_a ) ),
	u_dim_p_( new Dimension( *l_dim_p_ - 2 ) ) // Stays unbounded until upper index changed
{
	assert( legal_dynamic() );
	size_dynamic();
	insert_as_observer();
}


/// @brief Dimension and Unbounded Upper Index Assignment
DynamicIndexRange &
DynamicIndexRange::assign( Dimension const & l_dim_a, Star const & star )
{
	delete l_dim_p_; l_dim_p_ = l_dim_a.reference_copy(); l_insert_as_observer();
	delete u_dim_p_; u_dim_p_ = new Dimension( *l_dim_p_ - 2 ); u_insert_as_observer();
	Super::assign( l_dim_a.zvalue(), star );
	assert( legal_dynamic() );
	size_dynamic();
	notify();
	return *this;
}


/// @brief Expression and Unbounded Upper Index Assignment
DynamicIndexRange &
DynamicIndexRange::assign( Expression const & l_exp_a, Star const & star )
{
	delete l_dim_p_; l_dim_p_ = new Dimension( l_exp_a ); l_insert_as_observer();
	delete u_dim_p_; u_dim_p_ = new Dimension( *l_dim_p_ - 2 ); u_insert_as_observer();
	Super::assign( l_dim_p_->zvalue(), star );
	assert( legal_dynamic() );
	size_dynamic();
	notify();
	return *this;
}


/// @brief Unbounded Upper Index Set
DynamicIndexRange &
DynamicIndexRange::u( Star const & star )
{
	delete u_dim_p_;
	if ( l_dim_p_ ) {
		u_dim_p_ = new Dimension( *l_dim_p_ - 2 ); u_insert_as_observer();
	} else {
		u_dim_p_ = nullptr;
	}
	Super::u( star );
	assert( legal_dynamic() );
	size_dynamic();
	notify();
	return *this;
}


/// @brief Unbounded Upper Index Set Without Notification
DynamicIndexRange &
DynamicIndexRange::u_no_notify( Star const & star )
{
	delete u_dim_p_;
	if ( l_dim_p_ ) {
		u_dim_p_ = new Dimension( *l_dim_p_ - 2 ); u_insert_as_observer();
	} else {
		u_dim_p_ = nullptr;
	}
	Super::u( star );
	assert( legal_dynamic() );
	size_dynamic();
	return *this;
}


/// @brief Expand to Contain an Index
DynamicIndexRange &
DynamicIndexRange::contain( int const i )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	if ( Super::bounded() ) { // Bounded
		if ( l_ > i ) l_no_notify( i );
		if ( u_ < i ) u_no_notify( i );
	} else { // Unbounded
		if ( l_ > i ) assign_no_notify( i, star ); // Reset u_ to maintain unbounded state
	}
	assert( legal_dynamic() );
	notify();
	return *this;
}


/// @brief Expand to Contain an Index and Notify If Changed
DynamicIndexRange &
DynamicIndexRange::contain_nic( int const i )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	if ( Super::bounded() ) { // Bounded
		bool changed( false );
		if ( l_ > i ) {
			l_no_notify( i );
			changed = true;
		}
		if ( u_ < i ) {
			u_no_notify( i );
			changed = true;
		}
		if ( changed ) notify();
	} else { // Unbounded
		if ( l_ > i ) assign( i, star ); // Reset u_ to maintain unbounded state
	}
	assert( legal_dynamic() );
	return *this;
}


/// @brief Expand to Contain Another IndexRange
DynamicIndexRange &
DynamicIndexRange::contain( IndexRange const & I )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	assert( I.initialized() );
	if ( I.IndexRange::positive() ) {
		if ( Super::bounded() ) { // Bounded
			if ( l_ > I.l_ ) l_no_notify( I.l_ );
			if ( I.IndexRange::bounded() ) { // I bounded
				if ( u_ < I.u_ ) u_no_notify( I.u_ );
			} else { // I unbounded: Make this IndexRange unbounded
				u_no_notify( star );
			}
		} else { // Unbounded
			if ( l_ > I.l_ ) assign_no_notify( I.l_, star ); // Reset u_ to maintain unbounded state
		}
		assert( legal_dynamic() );
	}
	notify();
	return *this;
}


/// @brief Expand to Contain Another IndexRange and Notify If Changed
DynamicIndexRange &
DynamicIndexRange::contain_nic( IndexRange const & I )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	assert( I.initialized() );
	if ( I.IndexRange::positive() ) {
		if ( Super::bounded() ) { // Bounded
			bool changed( false );
			if ( l_ > I.l_ ) {
				l_no_notify( I.l_ );
				changed = true;
			}
			if ( I.IndexRange::bounded() ) { // I bounded
				if ( u_ < I.u_ ) {
					u_no_notify( I.u_ );
					changed = true;
				}
			} else { // I unbounded: Make this IndexRange unbounded
				u_no_notify( star );
				changed = true;
			}
			if ( changed ) notify();
		} else { // Unbounded
			if ( l_ > I.l_ ) assign( I.l_, star ); // Reset u_ to maintain unbounded state
		}
		assert( legal_dynamic() );
	}
	return *this;
}


/// @brief Intersect With Another IndexRange
DynamicIndexRange &
DynamicIndexRange::intersect( IndexRange const & I )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	assert( I.initialized() );
	if ( Super::intersects( I ) ) { // I and this DynamicIndexRange have positive size
		if ( l_ <= u_ ) { // Bounded with positive size
			if ( l_ < I.l_ ) l_no_notify( I.l_ );
			if ( ( I.l_ <= I.u_ ) && ( u_ > I.u_ ) ) u_no_notify( I.u_ );
		} else { // Unbounded
			if ( l_ < I.l_ ) {
				l_no_notify( I.l_ );
				if ( I.l_ <= I.u_ ) {
					u_no_notify( I.u_ );
				} else { // I is unbounded
					u_no_notify( star ); // Reset u_ to maintain unbounded state
				}
			} else if ( I.l_ <= I.u_ ) {
				u_no_notify( I.u_ );
			}
		}
	} else { // Empty intersection: Set zero size
		u_no_notify( l_ - 1 );
	}
	assert( legal_dynamic() );
	notify();
	return *this;
}


/// @brief Intersect With Another IndexRange and Notify If Changed
DynamicIndexRange &
DynamicIndexRange::intersect_nic( IndexRange const & I )
{
	assert( ( l_dim_p_ ? l_dim_p_->initialized_ : true ) && ( u_dim_p_ ? u_dim_p_->initialized_ : true ) );
	assert( I.initialized() );
	if ( Super::intersects( I ) ) { // I and this DynamicIndexRange have positive size
		if ( l_ <= u_ ) { // Bounded with positive size
			bool changed( false );
			if ( l_ < I.l_ ) {
				l_no_notify( I.l_ );
				changed = true;
			}
			if ( ( I.l_ <= I.u_ ) && ( u_ > I.u_ ) ) {
				u_no_notify( I.u_ );
				changed = true;
			}
			if ( changed ) notify();
		} else { // Unbounded
			bool changed( false );
			if ( l_ < I.l_ ) {
				l_no_notify( I.l_ );
				if ( I.l_ <= I.u_ ) {
					u_no_notify( I.u_ );
				} else { // I is unbounded
					u_no_notify( star ); // Reset u_ to maintain unbounded state
				}
				changed = true;
			} else if ( I.l_ <= I.u_ ) {
				u_no_notify( I.u_ );
				changed = true;
			}
			if ( changed ) notify();
		}
	} else { // Empty intersection: Set zero size
		u( l_ - 1 );
	}
	assert( legal_dynamic() );
	return *this;
}


// DynamicIndexRange


} // namespace ObjexxFCL
