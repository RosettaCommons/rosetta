// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Sarel Fleishman

#ifndef INCLUDED_basic_datacache_DataMapObj_hh
#define INCLUDED_basic_datacache_DataMapObj_hh

// Project headers
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace datacache {

/// @brief templated class to make any data type compatible with ReferenceCounts and OPs.
/// e.g., utility::pointer::owning_ptr< basic::datacache::DataMapObj< bool > > stop;
/// You can then place such constructs on the basic::datacache::DataMap
///
/// This class can now be depricated with the conversion to the (non-intrusive)
/// C++11 / boost shared_ptr system: arbitrary objects can now be put into smart
/// pointers.
template < class Ty >
class DataMapObj : public utility::pointer::ReferenceCount {
public:
	Ty obj;
};

#ifdef PYROSETTA
//template class DataMapObj<bool>;
class DataMapObj_bool : public DataMapObj<bool> {};
#endif

} // datacache
} // basic

#endif
