// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/roc_curve.fwd.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_roc_curve_FWD_HH
#define INCLUDED_numeric_roc_curve_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace numeric {


enum RocStatus {
	true_positive,
	true_negative,
	false_positive,
	false_negative
};

class RocPoint;
class RocCurve;

typedef utility::pointer::shared_ptr<RocPoint> RocPointOP;
typedef utility::pointer::shared_ptr<RocCurve> RocCurveOP;

}


#endif /*  INCLUDED_numeric_roc_curve_FWD_HH */
