// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/NV/NVlookup.fwd.hh
/// @brief  Typedefs and forward declarations for class NVlookup
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nv_NVlookup_fwd_hh
#define INCLUDED_core_scoring_nv_NVlookup_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace nv {

class NVlookup;
typedef utility::pointer::shared_ptr< NVlookup > NVlookupOP;
typedef utility::pointer::shared_ptr< NVlookup const > NVlookupCOP;

} //NV
} //scoring
} //core

#endif /*NVLOOKUP_FWD_HH_ */
