// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/CacheableString.hh
/// @brief
/// @author Phil Bradley


#ifndef INCLUDED_basic_datacache_CacheableString_hh
#define INCLUDED_basic_datacache_CacheableString_hh

// unit headers
#include <basic/datacache/CacheableString.fwd.hh>

// package headers
#include <basic/datacache/CacheableData.hh>

// C++ headers
#include <string>


namespace basic {
namespace datacache {


///@brief Wrapper for std::string
class CacheableString : public CacheableData
{
public:
	CacheableString( std::string str ) : CacheableData(), str_(str) {}
	virtual ~CacheableString(){};
	virtual CacheableDataOP clone() const { return new CacheableString(*this); }
	virtual std::string const & str() const { return str_; }
private:
	std::string str_;
};


} // namespace datacache
} // namespace basic

#endif /* INCLUDED_basic_datacache_CacheableString_HH */
