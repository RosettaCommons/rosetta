// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_viewer_viewers_hh
#define INCLUDED_protocols_viewer_viewers_hh

// Unit headers

// Package headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <protocols/viewer/GraphicsState.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/viewer/triangle.hh>

// Project headers
#include <utility/vector1.hh>

// C++ Headers
#include <string>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <utility/assert.hh>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <vector>


namespace protocols {
namespace viewer {

typedef void* (*VoidFunc)(void*);

static std::string empty_string("");

enum spheremode { SPHERE_MODE_BB, SPHERE_MODE_SC, SPHERE_MODE_LIGAND };

#ifndef GL_GRAPHICS ///////////////////////////////////////////////////////

inline
void
add_conformation_viewer(
	core::conformation::Conformation &,
	std::string const & = empty_string,
	int const = 900,
	int const = 900,
	bool const  = false,
	bool const  = false,
	core::Vector const = core::Vector( 0.0 ) ){
}

inline
void
add_monte_carlo_viewer(
	moves::MonteCarlo &,
	std::string const & = empty_string,
	int const = 900,
	int const = 900,
	bool = false
)
{
}

inline
int
viewer_main( VoidFunc worker_main ){ worker_main(NULL); return 0; }

#else // GL_GRAPHICS ////////////////////////////////////////////////////

void
add_conformation_viewer(
	core::conformation::Conformation & conformation,
	std::string const & = empty_string,
	int const length = 900,
	int const width = 900,
	bool const debug_pause = false,
	bool const set_center_vector = false,
	core::Vector const = core::Vector( 0.0 ) );

void
add_monte_carlo_viewer(
	moves::MonteCarlo & mc,
	std::string const & name_in = empty_string,
	int const length = 900,
	int const width  = 900,
	bool debug_pause=false
);


int
viewer_main( VoidFunc worker_main );

#endif

#if defined GL_GRAPHICS || BOINC_GRAPHICS

void
display_residues(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	core::id::AtomID const & anchor_id
);

//void
//display_residues_wireframe(
// utility::vector1< core::conformation::ResidueCOP > const & residues,
// core::id::AtomID const & anchor_id
//);

void
display_residues_wireframe(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	core::Vector const & center
);

void set_bg_color( core::Vector new_bg_color );

/// @brief Clear the background and fill it with the background colour.
///
void clear_bg();

/// @brief Draw a frame for a window.
///
void draw_frame();

/// @brief Fill window background with black.
///
void draw_black_bg();

/// @brief Draw a gradient for the window background.
///
void draw_gradient_bg();

void draw_pose(const core::pose::Pose & pose,
	GraphicsState & gs);

void draw_conformation( utility::vector1< core::conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss,
	GraphicsState & gs, core::Vector const & center);

void
render_density(
	GraphicsState &gs,
	utility::vector1< triangle > &triangles );

void
display_density(
	GraphicsState &gs,
	utility::vector1< triangle > &triangles ) ;

void
draw_conformation_and_density(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss,
	utility::vector1< triangle > &triangles,
	GraphicsState & gs,
	core::Vector const & center);

void
draw_sphere(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	spheremode const &sphere_mode
);

core::Vector
get_center( utility::vector1< core::conformation::ResidueCOP > const & residues );

void
add_monte_carlo_silent_viewer(
	moves::MonteCarlo & mc,
	std::string const name_in,
	bool fullatom
);

#endif ////////////////////////////////////////////////////////////////


/// @brief Allows for graceful exit of graphics viewers.
void
clear_conformation_viewers();

void
clear_conformation_viewers();


} // viewer
} // protocols


#endif
