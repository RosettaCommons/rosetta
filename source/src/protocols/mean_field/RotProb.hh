// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RotProb.hh

/// @brief   Declarations and simple accessor/mutator definitions for RotProb, simple class to encapsulate rotamer data and probability of rotamer occurring
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_RotProb_HH
#define INCLUDED_protocols_mean_field_RotProb_HH

// Unit header
#include <protocols/mean_field/RotProb.fwd.hh>

// Package headers


// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

// Utility headers
//#include <utility/pointer/ReferenceCount.fwd.hh>
#include <core/types.hh>

// Numeric headers


// C++ headers
#include <iostream>

namespace protocols {
namespace mean_field {

/// @details  encapsulates information about a rotamer at a position and its probability
/// @remarks  used as the base unit in the conformational matrix (RotMatrix) class
class RotProb {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Default constructor
	/// @remarks only necessary because will be used to initialize matrix of RotProb
	RotProb();

	/// @brief main constructor to initialize RotProb according to parameters
	RotProb( core::Real prob, core::Size rot_ind, core::Size pos, core::conformation::ResidueCOP res );

	/// @brief  Copy constructor
	RotProb(RotProb const & object_to_copy);

	/// @brief  Assignment operator
	RotProb & operator=(RotProb const & object_to_copy);

	/// @brief Destructor
	~RotProb();

	/// @brief overloaded equality operators to fix instantiation issues with PyRosetta
	inline
	bool operator==( RotProb const & obj ) const
	{
		bool is_equal = false;
		if ( probability() == obj.probability() &&
				rot_ind() == obj.rot_ind() &&
				pos() == obj.pos() &&
				res() == obj.res() ) {
			is_equal = true;
		}
		return is_equal;
	}

	/// @brief overloaded equality operators to fix instantiation issues with PyRosetta
	inline
	bool operator!=( RotProb const & obj ) const
	{
		return !(*this == obj);
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details adds probability of obj to its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator+=( RotProb const & obj )
	{
		probability( probability() + obj.probability() );
	}

	/// @brief overloaded arithmetic operators to allow for easier matrix arithmetic
	/// @details subtracts probability of obj from its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator-=( RotProb const & obj )
	{
		probability( probability() - obj.probability() );
	}

	/// @brief overloaded arithmetic operator not necessary but implemented due to its being expected
	/// @details multiplies probability of obj by its own probability
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator*=( RotProb const & obj )
	{
		probability( probability() * obj.probability() );
	}

	/// @brief overloaded arithmetic operator not necessary but implemented due to its being expected
	/// @details divides its probability by probability of obj
	/// @remarks only overloading operators that modify object so as to avoid confusion about identity of non-probability data members
	inline
	void operator/=( RotProb const & obj )
	{
		probability( probability() / obj.probability() );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details adds prob to probability
	inline
	void operator+=( core::Real prob )
	{
		probability( probability() + prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details subtracts prob from probability
	inline
	void operator-=( core::Real prob)
	{
		probability( probability() - prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details multiplies probability by prob
	inline
	void operator*=( core::Real prob )
	{
		probability( probability() * prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details divides probability by prob
	inline
	void operator/=( core::Real prob)
	{
		probability( probability() / prob );
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details adds prob to probability
	inline
	RotProb operator+( core::Real prob )
	{
		RotProb rp( (*this) );
		rp += prob;
		return rp;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details subtracts prob from probability
	inline
	RotProb operator-( core::Real prob )
	{
		RotProb rp( (*this) );
		rp -= prob;
		return rp;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details multiplies probability by prob
	inline
	RotProb operator*( core::Real prob )
	{
		RotProb rp( (*this) );
		rp *= prob;
		return rp;
	}

	/// @brief overloaded arithmetic operator for arithmetic with Real probability
	/// @details divides probability by prob
	inline
	RotProb operator/( core::Real prob )
	{
		RotProb rp( (*this) );
		rp /= prob;
		return rp;
	}

	// Standard Rosetta methods ////////////////////////////////////////////////

	/// @brief  Generate string representation of RotProb for debugging and reporting purposes.
	void show( std::ostream & output = std::cout ) const;

	/// @brief Insertion operator (overloaded so that RotProb can be "printed" in PyRosetta).
	friend std::ostream & operator<<( std::ostream & output, RotProb const & object_to_output );


	// Accessors/Mutators
	/// @brief returns probability
	inline
	core::Real probability() const
	{
		return probability_;
	}

	/// @brief set probability
	inline
	void probability( core::Real prob )
	{
		probability_ = prob;
	}

	/// @brief returns rotamer index
	inline
	core::Size rot_ind() const
	{
		return rot_ind_;
	}

	/// @brief set rotamer index
	inline
	void rot_ind( core::Size aa )
	{
		rot_ind_ = aa;
	}

	/// @brief returns position
	inline
	core::Size pos() const
	{
		return pos_;
	}

	/// @brief set position
	inline
	void pos( core::Size p )
	{
		pos_ = p;
	}

	/// @brief returns COP to residue
	core::conformation::ResidueCOP res() const;

	/// @brief set ResidueCOP
	void res( core::conformation::ResidueCOP r );

	/// @brief returns AA type based on residue
	core::chemical::AA aa_ind() const;


private:

	// Private methods /////////////////////////////////////////////////////////
	/// @brief Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(RotProb & object_to_copy_to, RotProb const & object_to_copy_from);

	// Private data ////////////////////////////////////////////////////////////
	core::Real probability_;
	core::Size rot_ind_;
	core::Size pos_; //position in pose numbering
	core::conformation::ResidueCOP res_;

};  // class RotProb

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_RotProb_HH
