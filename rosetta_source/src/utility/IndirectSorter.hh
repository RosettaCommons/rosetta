// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/IndirectSorter.hh
/// @brief  Indirect sort : sorts a container with some data without altering it
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_utility_IndirectSorter_hh
#define INCLUDED_utility_IndirectSorter_hh

// #include <utility/IndirectSorter.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <core/types.hh>

#include <algorithm>

namespace utility {

using namespace core;

/// @brief Sorts a container of things indirectly.
/// @detailed The original vector remains untouched. The class provides a vector
///	with indexes that defines the order user wants
template <class DataContainerType,class ComparatorType>
class IndirectSorter: public utility::pointer::ReferenceCount {
public:

    /// @brief Creates a new sorter for a given data container and sorting rule
    /// @detailed The given comparator must be a valid way to sort the given container,
    ///		e.g. when DataContainerType = vector1<Real> then ComparatorType
    ///		must be an object that implements bool operator()(Real,Real) method
    IndirectSorter(DataContainerType& data,ComparatorType& comp) : underlying_data_(data), comparator_(comp) {}

    ~IndirectSorter() {}

    /// @brief sorts a given container and returns a vector that defines the right order of the sorted data
    inline void sort( vector1<Size>& order) {

      order.clear();
      for(Size i=1;i<=underlying_data_.size();++i) {
        order.push_back(i);
      }
      std::sort( order.begin(),order.end(),*(this) );
    }

    /// @brief compares two indexes by comparing the two referenced objects
    /// @detailed returns the result of comparison between data[index_1] and data[index_2]
    /// 	The result of a comparison depends on a comparator definition provided by a user
    inline bool operator()(int index_1,int index_2) {

	return ( comparator_(underlying_data_[index_1],underlying_data_[index_2]) );
    }

private:
    DataContainerType &underlying_data_;
    ComparatorType &comparator_;
};

} // utility

#endif // INCLUDED_core_fragment_picking_IndirectSorter_HH
