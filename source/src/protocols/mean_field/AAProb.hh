// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    AAProb.hh

/// @brief   Declarations and simple accessor/mutator definitions for AAProb.
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_AAProb_HH
#define INCLUDED_protocols_mean_field_AAProb_HH

// Unit header
#include <protocols/mean_field/AAProb.fwd.hh>

// Package headers
#include <protocols/mean_field/RotProb.fwd.hh>

// Project headers
#include <core/chemical/AA.hh>

// Utility headers
#include <core/types.hh>

// Numeric headers



// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {


/// @details  encapsulates information about an amino acid at a position and its probability
/// @remarks  used as the base unit in the specificity profile (AAMatrix) class
class AAProb {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Default constructor
	/// @details only necessary because will be used to initialize matrix of AAProb
	AAProb();

	/// @brief  standard constructor
	AAProb( core::Real prob, core::chemical::AA aa_ind, core::Size pos, core::Size nrot );

	/// @brief constructor from RotProb
	AAProb( protocols::mean_field::RotProb const & rp );

	/// @brief  Copy constructor
	AAProb( AAProb const & object_to_copy );

	/// @brief Assignment operator
	AAProb & operator=( AAProb const & object_to_copy );

	/// @brief Destructor
	~AAProb();

	/// @brief overloaded equality operators to fix instantiation issues with PyRosetta
	inline
	bool operator==( AAProb const & obj ) const
	{
		bool is_equal = false;
		if ( probability() == obj.probability() &&
				aa_ind() == obj.aa_ind() &&
				pos() == obj.pos() &&
				nrot() == obj.nrot() ) {
			is_equal = true;
		}
		return is_equal;
	}

	/// @brief overloaded equality operators to fix instantiation issues with PyRosetta
	inline
	bool operator!=( AAProb const & obj ) const
	{
		return !(*this == obj);
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details adds probability of obj to its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator+=( AAProb const & obj )
	{
		probability( probability() + obj.probability() );
		nrot( nrot() + obj.nrot() );
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details subtracts probability of obj from its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator-=( AAProb const & obj )
	{
		probability( probability() - obj.probability() );
		nrot( nrot() - obj.nrot() );
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details multiplies probability of obj by its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator*=( AAProb const & obj )
	{
		probability( probability() * obj.probability() );
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details divides its own probability by probability of obj
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator/=( AAProb const & obj )
	{
		probability( probability() / obj.probability() );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details adds prob to its own probability
	inline
	void operator+=( core::Real prob )
	{
		probability( probability() + prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details subtracts prob from its own probability
	inline
	void operator-=( core::Real prob )
	{
		probability( probability() - prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details multiplies its own probability by prob
	inline
	void operator*=( core::Real prob )
	{
		probability( probability() * prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details divides its own probability by prob
	inline
	void operator/=( core::Real prob )
	{
		probability( probability() / prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details adds prob to probability
	inline
	AAProb operator+( core::Real prob )
	{
		AAProb ap( (*this) );
		ap += prob;
		return ap;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details subtracts prob from probability
	inline
	AAProb operator-( core::Real prob )
	{
		AAProb ap( (*this) );
		ap -= prob;
		return ap;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details multiplies probability by prob
	inline
	AAProb operator*( core::Real prob )
	{
		AAProb ap( (*this) );
		ap *= prob;
		return ap;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details divides probability by prob
	inline
	AAProb operator/( core::Real prob )
	{
		AAProb ap( (*this) );
		ap /= prob;
		return ap;
	}

	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of AAProb for debugging purposes.
	void show( std::ostream & output=std::cout ) const;

	/// @brief Insertion operator (overloaded so that AAProb can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, AAProb const & object_to_output );


	// Accessors/Mutators
	/// @brief returns probability
	inline
	core::Real probability() const
	{
		return probability_;
	}

	/// @brief sets probability
	inline
	void probability( core::Real prob )
	{
		probability_ = prob;
	}

	/// @brief returns AA identity
	inline
	core::chemical::AA aa_ind() const
	{
		return aa_ind_;
	}

	/// @brief sets AA identity
	inline
	void aa_ind( core::chemical::AA aa )
	{
		aa_ind_ = aa;
	}

	/// @brief returns position in pose numbering
	inline
	core::Size pos() const
	{
		return pos_;
	}

	/// @brief sets position in pose numbering
	inline
	void pos( core::Size p )
	{
		pos_ = p;
	}

	/// @brief returns number of rotamers corresponding to AA identity at this position
	inline
	core::Size nrot()  const {
		return nrot_;
	}

	/// @brief sets number of rotamers corresponding to AA identity at this position
	inline
	void nrot( core::Size n ) {
		nrot_ = n;
	}

private:
	// Private methods /////////////////////////////////////////////////////////
	/// @brief Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data( AAProb & object_to_copy_to, AAProb const & object_to_copy_from );

	// Private data ////////////////////////////////////////////////////////////
	core::Real probability_;
	core::chemical::AA aa_ind_;
	core::Size pos_; //position in pose numbering
	core::Size nrot_; //number of rotamers

};  // class AAProb

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_AAProb_HH
