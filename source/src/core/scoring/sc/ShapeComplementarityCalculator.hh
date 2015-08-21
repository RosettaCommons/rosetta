// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/sc/ShapeComplementarityCalculator.hh
/// @brief  Headers for the Shape Complementarity Calculator
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifndef INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_hh
#define INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_hh

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
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <numeric/xyzVector.hh>

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
// Shape Complementarity Calculator class definition
////////////////////////////////////////////////////////////

class ShapeComplementarityCalculator : public MolecularSurfaceCalculator {

public:

	ShapeComplementarityCalculator();
	virtual ~ShapeComplementarityCalculator();

	virtual int Calc();
	virtual int Calc(core::pose::Pose const & pose, core::Size jump_id);

	// Static function that can be called without object instantiation
	static core::Real CalcSc(core::pose::Pose const & pose, core::Size jump_id =1, int quick =0);

protected:
	// Surface generation configuration
	virtual int AssignAttentionNumbers(std::vector<Atom>& atom);

private:

	// Dot trimming
	ScValue TrimPeripheralBand(std::vector<DOT> const &sdots, std::vector<const DOT*> &trimmed_dots);
	int TrimPeripheralBandCheckDot(DOT const &dot, std::vector<DOT> const &sdots);

	// Shape Complementarity Kernel functions for parallalization
	int CalcNeighborDistance(int const molecule, std::vector<const DOT*> const &my_dots, std::vector<const DOT*> const &their_dots);
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

