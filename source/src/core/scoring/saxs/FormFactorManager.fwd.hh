// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/FormFactorManager.fwd.hh
/// @brief Forward declaration of a form factor manager
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_saxs_FormFactorManager_fwd_hh
#define INCLUDED_core_scoring_saxs_FormFactorManager_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace saxs {

class FormFactorManager;

typedef utility::pointer::shared_ptr<FormFactorManager>
		FormFactorManagerOP;
typedef utility::pointer::shared_ptr<FormFactorManager const>
		FormFactorManagerCOP;

} // saxs
} // scoring
} // core


#endif /* INCLUDED_core_scoring_saxs_FormFactorManager_FWD_HH */
