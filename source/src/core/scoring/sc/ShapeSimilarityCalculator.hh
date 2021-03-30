// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ShapeSimilarityCalculator.hh
/// @brief    Headers for the ShapeSimilarityCalculator
/// @details  The code modifies Luki Goldschmidt's
///           implementation of Lawrence & Coleman shape complementarity calculator
///           to allow for the comparison of similar surface shapes.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)

#ifndef INCLUDED_core_scoring_sc_ShapeSimilarityCalculator_hh
#define INCLUDED_core_scoring_sc_ShapeSimilarityCalculator_hh

// This code contains support for GPU acceleration using OpenCL.
// Build with scons option extras=opencl to enable GPU support.
//
// Note: Original Fotran code used floats rather than doubles. This code works
// correctly with either floats or core::Real, though core::Real is about 50%
// slower (tested on 64-bit machines). The code can be compiled with double (Real)
// precision floats by defining SC_PRECISION_REAL.

// #define SC_PRECISION_REAL
// #define USEOPENCL

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/sc/ShapeSimilarityCalculator.fwd.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <numeric/xyzVector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueVector.hh>



//// C++ headers
#include <vector>
#include <string>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////
// core::scoring::sc namespace
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Types


////////////////////////////////////////////////////////////
// Shape Similarity Calculator class definition
////////////////////////////////////////////////////////////

class ShapeSimilarityCalculator : public MolecularSurfaceCalculator {

public:
	ShapeSimilarityCalculator();
	~ShapeSimilarityCalculator() override;

	int Calc() override;


	core::Real CalcSs(
		core::pose::Pose const & pose,
		core::pose::Pose const & native,
		core::select::residue_selector::ResidueSelectorCOP selector_pose,
		core::select::residue_selector::ResidueSelectorCOP selector_native );

	bool median;

protected:
	// Surface generation configuration
	int AssignAttentionNumbers(std::vector<Atom>& atom) override;

private:

	// Dot trimming
	ScValue TrimPeripheralBand(std::vector<DOT> const &sdots, std::vector<const DOT*> &trimmed_dots);
	int TrimPeripheralBandCheckDot(DOT const &dot, std::vector<DOT> const &sdots);

	// Shape Similarity Kernel functions for parallalization
	core::Real CalcNeighborDistance(int const molecule, std::vector<const DOT*> const &my_dots, std::vector<const DOT*> const &their_dots);
	DOT const *CalcNeighborDistanceFindClosestNeighbor(DOT const &dot1, std::vector<const DOT*> const &their_dots);

#ifdef USEOPENCL
	protected:
	basic::gpu::GPU gpu;
	float gpuTrimPeripheralBand(const std::vector<DOT> &dots, std::vector<const DOT*> &trimmed_dots);
	int gpuFindClosestNeighbors(const std::vector<const DOT*> &my_dots, const std::vector<const DOT*> &their_dots, std::vector<const DOT*> &neighbors);
	core::Real GetTimerMs(clock_t &start);

	public:
	void gpuInit();
#endif



};

} //namespace sc
} //namespace filters
} //namespace protocols

#endif
