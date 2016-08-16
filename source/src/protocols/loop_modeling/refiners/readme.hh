// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @namespace protocols::loop_modeling::refiners
///
/// @brief Algorithms for lowering loop region scores during sampling.
///
/// @details the classes in the namespace are broadly supposed to be
/// responsible for lowering the score in loop regions sampled using the
/// algorithms in the samplers namespace.  The primary refinement algorithms
/// are RepackingRefiner, RotamerTrialsRefiner, and MinimizationRefiner.
/// These algorithms together account for the lion's share of the total running
/// time in most loop modeling simulations.
