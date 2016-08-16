// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author


#ifndef INCLUDED_protocols_viewer_marchingCubes_hh
#define INCLUDED_protocols_viewer_marchingCubes_hh

namespace protocols {
namespace viewer {


// This here is the numeric information needed to do marching cubes.

extern const int POLY_CASES[][21];

extern const int VERTEX_OFF[][3];

extern const int EDGE_NGHBRS[][3];

}
}

#endif
