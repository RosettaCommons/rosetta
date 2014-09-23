// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilters.fwd.hh
/// @brief  core-level (very general) classes that take a pose and a residue index, and returns true or false
/// @author ashworth

#ifndef INCLUDED_core_pack_task_operation_ResFilters_fwd_hh
#define INCLUDED_core_pack_task_operation_ResFilters_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

class ResidueHasProperty;
class ResiduePDBInfoHasLabel;
class ResidueLacksProperty;
class ResidueName3Is;
class ResidueName3Isnt;
class ResidueIndexIs;
class ResidueIndexIsnt;
class ResiduePDBIndexIs;
class ResiduePDBIndexIsnt;
class ChainIs;
class ChainIsnt;

typedef utility::pointer::shared_ptr< ResidueHasProperty > ResidueHasPropertyOP;
typedef utility::pointer::shared_ptr< ResiduePDBInfoHasLabel >  ResiduePDBInfoHasLabelOP;
typedef utility::pointer::shared_ptr< ResidueLacksProperty > ResidueLacksPropertyOP;
typedef utility::pointer::shared_ptr< ResidueName3Is > ResidueName3IsOP;
typedef utility::pointer::shared_ptr< ResidueName3Isnt > ResidueName3IsntOP;
typedef utility::pointer::shared_ptr< ResidueIndexIs > ResidueIndexIsOP;
typedef utility::pointer::shared_ptr< ResidueIndexIsnt > ResidueIndexIsntOP;
typedef utility::pointer::shared_ptr< ResiduePDBIndexIs > ResiduePDBIndexIsOP;
typedef utility::pointer::shared_ptr< ResiduePDBIndexIsnt > ResiduePDBIndexIsntOP;
typedef utility::pointer::shared_ptr< ChainIs > ChainIsOP;
typedef utility::pointer::shared_ptr< ChainIsnt > ChainIsntOP;

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
