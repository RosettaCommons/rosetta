// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/wwPDBJobOutputter.fwd.hh
/// @brief  header file for wwPDBJobOutputter class, interstitial between FileJobOutputter and PDBJobOutputter or mmCIFJobOutputter to share functionality.
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_wwPDBJobOutputter_fwd_hh
#define INCLUDED_protocols_jd2_wwPDBJobOutputter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class wwPDBJobOutputter;
typedef utility::pointer::shared_ptr< wwPDBJobOutputter > wwPDBJobOutputterOP;

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_wwPDBJobOutputter_FWD_HH
