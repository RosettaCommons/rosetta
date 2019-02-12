// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/tcr/TCRmodel.fwd.hh
/// @brief  TCRmodel class forward declarations header
/// @author Ragul Gowthaman (ragul@umd.com)


#ifndef INCLUDED_protocols_tcr_TCRmodel_fwd_hh
#define INCLUDED_protocols_tcr_TCRmodel_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// C++ Headers
namespace protocols {
namespace tcr {

// Forward
class TCRmodel;

typedef utility::pointer::shared_ptr< TCRmodel > TCRmodelOP;
typedef utility::pointer::shared_ptr< TCRmodel const > TCRmodelCOP;
typedef  utility::pointer::weak_ptr< TCRmodel >  TCRmodelAP;
typedef  utility::pointer::weak_ptr< TCRmodel const >  TCRmodelCAP;


} //namespace tcr
} //namespace protocols

#endif //INCLUDED_protocols_TCRmodel_FWD_HH
