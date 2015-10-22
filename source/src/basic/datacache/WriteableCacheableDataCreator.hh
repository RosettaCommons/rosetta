// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/WriteableCacheableDataCreator.hh
/// @brief  Base class for WriteableCacheableDataCreators for the WriteableCacheableData load-time factory registration scheme
/// @author Justin Porter

#ifndef INCLUDED_basic_datacache_WriteableCacheableDataCreator_hh
#define INCLUDED_basic_datacache_WriteableCacheableDataCreator_hh

// Unit Headers
#include <basic/datacache/WriteableCacheableData.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <iosfwd>

namespace basic {
namespace datacache {

/// @brief Abstract base class for a Mover factory; the Creator class is responsible for
/// creating a particular mover class.
class WriteableCacheableDataCreator : public utility::pointer::ReferenceCount
{
public:
	WriteableCacheableDataCreator();
	virtual ~WriteableCacheableDataCreator();

	virtual WriteableCacheableDataOP create_data( std::istream &in ) const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< WriteableCacheableDataCreator > WriteableCacheableDataCreatorOP;
typedef utility::pointer::shared_ptr< WriteableCacheableDataCreator const > WriteableCacheableDataCreatorCOP;

} //namespace datacache
} //namespace basic

#endif
