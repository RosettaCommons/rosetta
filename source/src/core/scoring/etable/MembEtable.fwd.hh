// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/memb_etable/MembEtable.fwd.hh
///
/// @brief		Generate the table for fa_atr/rep and fa_sol with membrane additions
/// @details	Used by the scoring manager. becasue computing LJ potentials is time
///				consuming, precomputes and discritizes the potential (broken down into bins).
///				Once bins are created, will smooth bins for better interpolation.
///				Last Modified: 5/13/14
///
/// @author		Patrick Barth
/// @author		(Updates) Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_etable_MembEtable_fwd_hh
#define INCLUDED_core_scoring_etable_MembEtable_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace etable {

class MembEtable;

typedef utility::pointer::owning_ptr< MembEtable > MembEtableOP;
typedef utility::pointer::access_ptr< MembEtable const > MembEtableCAP;

} // etable
} // scoring
} // core

#endif
