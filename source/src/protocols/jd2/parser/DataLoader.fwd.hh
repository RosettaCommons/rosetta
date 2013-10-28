// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.fwd.hh
/// @brief  Forward declartion of the XML parser's DataLoader class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd2_parser_DataLoader_fwd_hh
#define INCLUDED_protocols_jd2_parser_DataLoader_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {
namespace parser {

/// @brief A class for loading arbitrary data into the XML parser's basic::datacache::DataMap.
class DataLoader;

typedef utility::pointer::owning_ptr< DataLoader > DataLoaderOP;
typedef utility::pointer::owning_ptr< DataLoader const > DataLoaderCOP;

} //namespace parser
} //namespace jd2
} //namespace protocols

#endif
