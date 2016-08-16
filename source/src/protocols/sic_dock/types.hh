// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_sic_dock_types_hh
#define INCLUDED_protocols_sic_dock_types_hh

#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace sic_dock {

struct Vec3 { numeric::xyzVector<platform::Real> a,b,c; };
typedef utility::vector1<std::pair<platform::Size,Vec3> > TermInfo;

typedef numeric::xyzVector<platform::Real> Vec;
typedef numeric::xyzMatrix<platform::Real> Mat;
typedef numeric::xyzTransform<platform::Real> Xform;
typedef utility::vector1<Vec> Vecs;
typedef utility::vector1<Mat> Mats;
using numeric::Xforms;

}
}

#endif
