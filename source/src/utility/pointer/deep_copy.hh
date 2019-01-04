// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/deep_copy.hh
/// @brief  This is a utility to assist in making deep copies of OPs
///
/// The rationale is that if you have 100 member variables which are a straight copy,
/// but 1 member variable which needs a deep copy, then you have to manually make a copy constructor.
/// This is suboptimal, not only because you have to manually copy all the member variables,
/// but because when someone adds the 102st member variable, chances are they'll forget the copy constructor.
/// Also, having a user-defined copy constructor turns off the automatic move constructor/assignment.
///
/// By marking your member variable a DeepCopy, you can use it like a normal OP, but it will be deep copied
/// in the default copy constructor of your class.
///
/// IMPORTANT: This shouldn't be used for OPs inside of containers. Because the standard containers routinely copy their objects,
/// using a DeepCopyOP inside a container (versus directly in the object) will result in too much copying.
/// This class should only be used for "top level" members in classes.
/// (Also, don't use this class for parameters/return values - use a regular OP instead.)
///
/// Usage in class is like the following:
///
///     utility::pointer::DeepCopyOP< MyObject > my_object_;
///
/// which should be sufficient to make a deep copy of my_object_ whenever the containing object is copied.
/// The resulting DeepCopyOP object should (ideally) function identically to a regular OP (or COP).
///
/// You'll also need to make sure that there's a free deep_copy() object which implements the deep copy.
///
///     MyObjectOP deep_copy( MyObject const & source);
///
/// This function should be a free function in the same namespace as MyObject (for Koenig lookup purposes).
/// This should be declared in the MyObject *forward header*, and defined either in the MyObject header or cc.
/// (It needs to be declared in the forward header to allow rule-of-zero deletion of the copy constructor.)
///
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_utility_pointer_deep_copy_hh
#define INCLUDED_utility_pointer_deep_copy_hh

#include <utility/pointer/owning_ptr.fwd.hh>

#include <cstddef>
#include <type_traits>

namespace utility {
namespace pointer {

template< class T >
class DeepCopyOP {
public:

	typedef utility::pointer::shared_ptr< T > pointer_type;

	constexpr DeepCopyOP() noexcept: val_(nullptr) {}
	constexpr DeepCopyOP( std::nullptr_t const & null) noexcept: val_(null) {}
	DeepCopyOP( DeepCopyOP::pointer_type const & val ): val_(val) {}
	explicit DeepCopyOP( T * val ): val_(val) {}

	/// @brief Copy constructor -- This is what the effort is about.
	DeepCopyOP( DeepCopyOP const & r ) noexcept: val_( (r.val_ == nullptr) ? nullptr : deep_copy(*r.val_) ) {}

	/// @brief Move constructor
	/// @details Because the rvalue is disappearing, there's no need to actually make a deep copy
	/// -- there's no reference which should be sharing this value
	DeepCopyOP( DeepCopyOP const && r ) noexcept: val_( r.val_ ) {}

	~DeepCopyOP() = default;

	/// @brief Standard assignment operator, needed for implicit assignment operator of containing class
	// (Which should clone to prevent from sharing values with the RHS of the object.)
	DeepCopyOP & operator=( DeepCopyOP const & r ) {
		val_ = (r.val_ == nullptr) ? nullptr : deep_copy(*r.val_);
		return *this;
	}

	/// @brief Move assignment operator.
	/// @details The RHS is disappearing, so we can steal the value
	DeepCopyOP & operator=( DeepCopyOP && r ) {
		std::swap(val_,r.val_);
		return *this;
	}

	/// @brief Assignment from OP.
	/// @details This is what's used to manually set the value within the class.
	/// We don't need to clone.
	DeepCopyOP & operator=( DeepCopyOP::pointer_type const & val ) {
		val_ = val;
		return *this;
	}

	/// @brief Move assignement operator for OPs
	DeepCopyOP & operator=( DeepCopyOP::pointer_type && val ) {
		std::swap(val_,val);
		return *this;
	}

	/// @details We need to overload the nullptr_t otherwise we get compiler issues
	DeepCopyOP & operator=( std::nullptr_t const & val) {
		val_ = val;
		return *this;
	}

	/// @brief Conversion operator to allow a DeepCopyOP to be passed as a plain OP
	/// @details While this actually won't allow you to convert a COP to an OP (you'll get a obtuse-but-interpretable compiler error if you try)
	/// the remove_const is needed to avoid an unconditional compiler error for just using DeepCopyOP with a COP
	operator shared_ptr< typename std::remove_const< T >::type >() const { return val_; }

	/// @brief Conversion operator which allows an OP -> COP change.
	operator shared_ptr< typename std::add_const< T >::type >() const { return val_; }

	void swap( DeepCopyOP & r ) noexcept { swap(val_, r.val_); }
	void swap( DeepCopyOP::pointer_type & r ) noexcept { swap(val_, r); }

	T* get() const noexcept { return val_.get(); }

	T & operator*() const noexcept { return *val_; }
	T * operator->() const noexcept { return val_.get(); }

	typename DeepCopyOP::pointer_type get_op() const { return val_; }

	explicit operator bool() const noexcept { return bool(val_); }

	bool operator ==( DeepCopyOP const & r ) const {
		return val_ == r.val_;
	}

	bool operator ==( DeepCopyOP::pointer_type const & val ) const {
		return val_ == val;
	}

	bool operator ==( std::nullptr_t const & null ) const {
		return val_ == null;
	}

	bool operator !=( DeepCopyOP const & r ) const {
		return val_ != r.val_;
	}

	bool operator !=( DeepCopyOP::pointer_type const & val ) const {
		return val_ != val;
	}

	bool operator !=( std::nullptr_t const & null ) const {
		return val_ != null;
	}

private:
	DeepCopyOP::pointer_type val_;

};

/// @details It looks like dynamic_point_cast<>() can't recognize the implicit conversion ... correct for this.
template< class OUT, class T >
utility::pointer::shared_ptr< OUT >
dynamic_pointer_cast( DeepCopyOP< T > const & param ) {
	return dynamic_pointer_cast< OUT >( param.get_op() );
}

} // namespace pointer
} // namespace utility

#endif
