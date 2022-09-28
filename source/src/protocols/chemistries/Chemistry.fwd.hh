// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/chemistries/Chemistry.fwd.hh
///
/// @brief Forward headers for Chemistry
/// @author Steven Combs


#ifndef INCLUDED_protocols_chemistries_Chemistry_fwd_hh
#define INCLUDED_protocols_chemistries_Chemistry_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace chemistries {

class Chemistry;

typedef utility::pointer::shared_ptr<Chemistry> ChemistryOP;
typedef utility::pointer::shared_ptr<Chemistry const> ChemistryCOP;

}
}



#endif
