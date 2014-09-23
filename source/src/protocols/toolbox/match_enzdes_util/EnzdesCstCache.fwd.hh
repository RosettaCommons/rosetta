// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file fwd.hh file for enzdes constraint cache
/// @brief
/// @author Florian Richter, floric@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCstCache_fwd_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_EnzdesCstCache_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class EnzdesCstCache;
class EnzdesCstParamCache;
class EnzCstTemplateResCache;

typedef utility::pointer::shared_ptr< EnzdesCstCache > EnzdesCstCacheOP;
typedef utility::pointer::shared_ptr< EnzdesCstCache const > EnzdesCstCacheCOP;

typedef utility::pointer::shared_ptr< EnzdesCstParamCache > EnzdesCstParamCacheOP;
typedef utility::pointer::shared_ptr< EnzdesCstParamCache const > EnzdesCstParamCacheCOP;

typedef utility::pointer::shared_ptr< EnzCstTemplateResCache > EnzCstTemplateResCacheOP;
typedef utility::pointer::shared_ptr< EnzCstTemplateResCache const > EnzCstTemplateResCacheCOP;

}
} // enzdes
} //protocols


#endif
