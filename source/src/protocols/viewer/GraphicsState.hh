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

#ifndef INCLUDED_protocols_viewer_GraphicsState_hh
#define INCLUDED_protocols_viewer_GraphicsState_hh

// Unit headers

// Package headers
#include <core/types.hh>
#include <numeric/xyzVector.hh>

// Project headers

// C++ Headers

namespace protocols {
namespace viewer {
namespace graphics_states_param {
//lin backbone display state
const std::size_t Num_BBdisplayState = 4;
enum BBdisplayState { SHOW_NOBB, SHOW_CARTOON, SHOW_BACKBONE, SHOW_BBSPHERES };

//lin sidechain display state
const std::size_t Num_SCdisplayState = 4;
enum SCdisplayState { SHOW_NOSC, SHOW_STICK, SHOW_WIREFRAME, SHOW_SCSPHERES};

//lin color state
const std::size_t Num_ColorModes = 9;
enum ColorMode { RAINBOW_COLOR, CPK_COLOR, RESIDUE_COLOR, CHAIN_COLOR, RAINBOW_CPK_COLOR, RESIDUE_CPK_COLOR, RHIJU_COLOR, GHOST_COLOR };

//lin trajectory state
const std::size_t Num_TrajectoryState = 5;
enum TrajectoryState { SHOW_LOW, SHOW_BEST, SHOW_MC_TRIALS, SHOW_ALL_TRIALS };

// H state
const std::size_t Num_ShowHState = 2;
enum ShowHState { SHOW_NO_H, SHOW_H };

}
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace graphics_states_param;

//lin define the graphics state
class GraphicsState {
public:
	BBdisplayState BBdisplay_state;
	SCdisplayState SCdisplay_state;
	ColorMode Color_mode;
	TrajectoryState Trajectory_state;
	ShowHState show_H_state;
	core::Vector previous_vertex1, previous_vertex2, previous_width_vector;

	core::Real density_sigma;  // contour level of density
	bool density_redraw;
	std::size_t nres_for_graphics;

	GraphicsState() :
		BBdisplay_state (SHOW_CARTOON),//default
		SCdisplay_state (SHOW_STICK),//default
		Color_mode (RAINBOW_COLOR),//default
		Trajectory_state (SHOW_ALL_TRIALS),//default
		show_H_state (SHOW_NO_H), //default
		previous_vertex1( 0.0 ), previous_vertex2( 0.0 ), previous_width_vector( 0.0 ),
		density_sigma( 2.0 ), density_redraw(true),  //default
		nres_for_graphics( 0 )
	{}

	GraphicsState(
		BBdisplayState BBdisplay_state_in,
		SCdisplayState SCdisplay_state_in,
		ColorMode Color_mode_in,
		TrajectoryState Trajectory_state_in,
		ShowHState show_H_state_in
	) :
		BBdisplay_state (BBdisplay_state_in),
		SCdisplay_state (SCdisplay_state_in),
		Color_mode (Color_mode_in),
		Trajectory_state (Trajectory_state_in),
		show_H_state (show_H_state_in),
		previous_vertex1( 0.0 ), previous_vertex2( 0.0 ), previous_width_vector( 0.0 ),
		density_sigma( 2.0 ), density_redraw(true),  // default
		nres_for_graphics( 1 )
	{}
};

} // viewer
} // protocols


#endif
