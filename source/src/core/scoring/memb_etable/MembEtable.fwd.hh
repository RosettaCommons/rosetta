// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/memb_etable/MembEtable.fwd.hh
/// @brief Table of pre-computed LK membrane solvation energies
/// @author Patrick Barth (original)
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_memb_etable_MembEtable_fwd_hh
#define INCLUDED_core_scoring_memb_etable_MembEtable_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace etable {

class MembEtable;

typedef utility::pointer::shared_ptr< MembEtable > MembEtableOP;
typedef utility::pointer::shared_ptr< MembEtable const > MembEtableCOP;
typedef utility::pointer::weak_ptr< MembEtable const > MembEtableCAP;

} // etable
} // scoring
} // core

#endif
