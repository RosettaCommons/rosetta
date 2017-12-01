// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/MPSpanAngleFilter.fwd.hh
/// @brief calculates the angle between the TM span and the membrane normal
/// @details naturla TM span angles are generally distributed between 0-60 degrees. this filter claculates the angle
/// of the speicifed TM span and either returns the angle, or a score corresponding to how different that angle is
/// from the natural distribution
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)

#ifndef INCLUDED_protocols_simple_filters_MPSpanAngleFilter_fwd_hh
#define INCLUDED_protocols_simple_filters_MPSpanAngleFilter_fwd_hh


namespace protocols {
namespace simple_filters {

class MPSpanAngleFilter;

}
}

#endif
