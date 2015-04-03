// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /protocols/filters/filters.fwd.hh
/// @brief  forward declaration for Filter base class
/// @author Florian Richter, floric@u.washington.edu, Sarel Fleishman sarelf@u.washington.edu


#ifndef INCLUDED_protocols_filters_Filter_fwd_hh
#define INCLUDED_protocols_filters_Filter_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>
#include <string>
#include <map>
#include <boost/shared_ptr.hpp>

#include <utility/vector1.fwd.hh>


namespace protocols {
namespace filters {

// Forward
class Filter;
class FilterCollection;

// Types
typedef utility::pointer::shared_ptr< Filter >  FilterOP;
typedef utility::pointer::shared_ptr< Filter const >  FilterCOP;
typedef boost::shared_ptr < Filter > FilterSP;

typedef utility::pointer::shared_ptr< FilterCollection >  FilterCollectionOP;
typedef utility::pointer::shared_ptr< FilterCollection const >  FilterCollectionCOP;

typedef utility::vector1< FilterOP > FilterOPs;
typedef utility::vector1< FilterCOP > FilterCOPs;

typedef std::map< std::string const, FilterOP > Filters_map;

} // namespace protocols
} // namespace filters

#endif
