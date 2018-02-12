// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/TagManager.fwd.hh
/// @brief
/// @detailed
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_TagManager_FWD_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_TagManager_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>
#include <map>

namespace protocols {
namespace multistage_rosetta_scripts {

class TagManager;
typedef utility::pointer::shared_ptr< TagManager > TagManagerOP;
typedef utility::pointer::shared_ptr< TagManager const > TagManagerCOP;

class NoFailDataMap;
typedef utility::pointer::shared_ptr< NoFailDataMap > NoFailDataMapOP;

typedef std::list< utility::tag::TagCOP > TagList;
typedef utility::pointer::shared_ptr< std::list< utility::tag::TagCOP > > TagListOP;
typedef utility::pointer::shared_ptr< std::list< utility::tag::TagCOP > const > TagListCOP;

typedef std::map< std::string, utility::tag::TagCOP > TagMap;
typedef utility::pointer::shared_ptr< std::map< std::string, utility::tag::TagCOP > > TagMapOP;
typedef utility::pointer::shared_ptr< std::map< std::string, utility::tag::TagCOP > const > TagMapCOP;

struct ParsedTagCache;
typedef utility::pointer::shared_ptr< ParsedTagCache > ParsedTagCacheOP;

} //multistage_rosetta_scripts
} //protocols

#endif
