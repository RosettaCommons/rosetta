// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/DiscreteIntervalEncodingTree.srlz.hh
/// @brief  serialization routines for DIETs
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


#ifdef    SERIALIZATION

namespace numeric {

template<>
template< class Archive >
void
DietNode<platform::Size>::save( Archive & arc ) const {
	arc(lower_);
	arc(upper_);
	arc(left_);
	arc(right_);
}

template<>
template< class Archive >
void
DietNode<platform::Size>::load( Archive & arc ) {
	arc(lower_);
	arc(upper_);
	arc(left_);
	arc(right_);
}

template<>
template< class Archive >
void
DiscreteIntervalEncodingTree<platform::Size>::save( Archive & arc ) const {
	arc(root_);
}

template<>
template< class Archive >
void
DiscreteIntervalEncodingTree<platform::Size>::load( Archive & arc ) {
	arc(root_);
}

}

SAVE_AND_LOAD_SERIALIZABLE( numeric::DietNode<platform::Size> );
//CEREAL_REGISTER_TYPE( numeric::DietNode<platform::Size> )

SAVE_AND_LOAD_SERIALIZABLE( numeric::DiscreteIntervalEncodingTree<platform::Size> );
//CEREAL_REGISTER_TYPE( numeric::DiscreteIntervalEncodingTree<platform::Size> )

//CEREAL_REGISTER_DYNAMIC_INIT( numeric_diet )

#endif // SERIALIZATION

