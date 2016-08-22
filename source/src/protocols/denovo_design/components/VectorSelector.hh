// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/VectorSelector.hh
/// @brief Selects items from a vector using different methods
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_VectorSelector_hh
#define INCLUDED_protocols_denovo_design_components_VectorSelector_hh

#include <protocols/denovo_design/components/VectorSelector.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Selects items from a vector using different methods
template< class T >
class VectorSelector : public utility::vector1< T > {
public:
	virtual T
	select() = 0;
};

/// @brief Selects items from a vector randomly
template< class T >
class RandomVectorSelector : public VectorSelector< T > {
public:
	RandomVectorSelector( core::Size const max_count ):
		VectorSelector< T >(),
		max_count_( max_count ),
		count_( 1 )
	{}

	virtual T
	select()
	{
		core::Real random = numeric::random::rg().uniform();
		core::Size const idx = extract_int( random, 1, this->size() );
		++count_;
		return (*this)[ idx ];
	}

private:
	RandomVectorSelector() {};

private:
	core::Size const max_count_;
	core::Size count_;
};

/// @brief Selects items sequentially from a vector
template< class T >
class EnumeratedVectorSelector : public VectorSelector< T > {
public:
	EnumeratedVectorSelector():
		VectorSelector< T >(),
		count_( 1 )
	{};

	virtual T
	select()
	{
		if ( count_ > this->size() ) reset();
		std::cout << "Count is now " << count_ << " Size is " << this->size() << std::endl;
		T selected = (*this)[ count_ ];
		++count_;
		return selected;
	}

	/// @brief To reset, randomly shuffle and reset counter
	void
	reset()
	{
		count_ = 1;
		numeric::random::random_permutation( *this, numeric::random::rg() );
	}

private:
	core::Size count_;
};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_VectorSelector_hh

