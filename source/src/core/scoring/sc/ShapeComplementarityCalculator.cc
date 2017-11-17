// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ShapeComplementarityCalculator.cc
/// @brief    Headers for the Shape Complementarity Calculator
/// @details Lawrence & Coleman shape complementarity calculator (based on CCP4's sc)
/// @author   Luki Goldschmidt <luki@mbi.ucla.edu>

/// This code was ported from the original Fortran code found in CCP4:
/// Sc (Version 2.0): A program for determining Shape Complementarity
/// Copyright Michael Lawrence, Biomolecular Research Institute
/// 343 Royal Parade Parkville Victoria Australia
///
/// This version contains support for GPU-acceleration OpenCL-
/// capable devices, which provides a 10-25x speed up over the CPU-only code
/// using a regular desktop video card with 4 processors (32 cores).
/// Build with scons option extras=opencl to enable GPU support.

#ifndef INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_cc
#define INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_cc

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>

#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// C headers
#include <stdio.h>

// C++ headers
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

#define UPPER_MULTIPLE(n,d) (((n)%(d)) ? (((n)/(d)+1)*(d)) : (n))

static basic::Tracer TR( "core.scoring.sc.ShapeComplementarityCalculator" );

using namespace core;

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////////////////////
// Public class functions
////////////////////////////////////////////////////////////////////////////

/// @brief
/// ShapeComplementarityCalculator constructor, initializes default settings

ShapeComplementarityCalculator::ShapeComplementarityCalculator() :
	MolecularSurfaceCalculator()
{
}

ShapeComplementarityCalculator::~ShapeComplementarityCalculator()
{
}

/// @brief
/// Run the SC calculation on Pose and return just the sc statistic or -1 on error
/// @details
/// This is a static function and can be called without instantiating ShapeComplementarityCalculator.
/// The jump_id is used to partition the pose into two molecular surfaces; the first jump (1)
/// is used is no jump_id is explicity specified. Those desiring more control as to what residues
/// make up either surface should use the AddResidue() or even add_atom() function instead.
/// Setting quick to true will perform a much faster calculation (~5-10 times faster) at the expense
/// of accuracy (about 0.05 units).
///
/// Example:
/// core::Real sc = core::scoring::sc::ShapeComplementarityCalculator( pose );
core::Real ShapeComplementarityCalculator::CalcSc(core::pose::Pose const & pose, core::Size jump_id, int quick)
{
	ShapeComplementarityCalculator sc;

	if ( quick ) sc.settings.density = 5;

	if ( sc.Calc(pose, jump_id) ) {
		return sc.GetResults().sc;
	} else {
		return -1;
	}
}

/// @brief Run the SC calculation on a Pose, partitionied by jump_id
/// @details
/// This non-static function requires an instance of the ShapeComplementarityCalculator class.
/// The jump_id is used to partition the pose into two molecular surfaces. To control what
/// residues make up either surface, use the AddResidue() or even add_atom() function instead.
/// Returns true on success. Results are retrieved with GetResults().
///
/// Example:
/// core::scoring::sc::ShapeComplementarityCalculator calc;
/// core::Real sc;
/// if(calc.Calc( pose ))
///   sc = calc.GetResults().sc;
int ShapeComplementarityCalculator::Calc(core::pose::Pose const & pose, core::Size jump_id)
{
	if ( jump_id > pose.num_jump() || jump_id <= 0 ) {
		TR.Error << "Jump ID out of bounds (pose has " << pose.num_jump() << " jumps)" << std::endl;
		return 0;
	}

	return MolecularSurfaceCalculator::Calc(pose, jump_id);
}

/// @brief Run the SC calculation for previously defined molecules (via AddResidue or add_atom calls)
/// @details
/// This function should be called the residues / atoms making up the two molecular surfaces
/// have been explicitly defined.
/// Returns true on success.

int ShapeComplementarityCalculator::Calc()
{
#ifdef USEOPENCL
	gpuInit();
#endif

	try {
		basic::gpu::Timer timer(TR.Debug);

		run_.results.valid = 0;

		if ( run_.atoms.empty() ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined");
		}
		if ( !run_.results.surface[0].nAtoms ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined for molecule 1");
		}
		if ( !run_.results.surface[1].nAtoms ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No atoms defined for molecule 2");
		}

		// Determine and assign the attention numbers for each atom
		AssignAttentionNumbers(run_.atoms);

		GenerateMolecularSurfaces();

		if ( !run_.dots[0].size() || !run_.dots[1].size() ) {
			throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "No molecular dots generated!");
		}

		// Cut away the periphery of each surface
		TR.Debug << "Trimming peripheral band, " << settings.band << "A range" << std::endl;

		std::vector<DOT const *> trimmed_dots[2];
		for ( int i = 0; i < 2; ++i ) {
			run_.results.surface[i].trimmedArea = TrimPeripheralBand(run_.dots[i], trimmed_dots[i]);
			if ( !trimmed_dots[i].size() ) {
				throw CREATE_EXCEPTION(utility::excn::Exception, "No molecular dots for surface " + utility::to_string( i ));
			}
			run_.results.surface[i].nTrimmedDots = trimmed_dots[i].size();
			run_.results.surface[i].nAllDots = run_.dots[i].size();
		}

		// Compute distance arrays and histograms for each surface
		TR.Debug << "Computing surface separation and vectors" << std::endl;

		CalcNeighborDistance(0, trimmed_dots[0], trimmed_dots[1]);
		CalcNeighborDistance(1, trimmed_dots[1], trimmed_dots[0]);

		run_.results.surface[2].d_mean = (run_.results.surface[0].d_mean + run_.results.surface[1].d_mean) / 2;
		run_.results.surface[2].d_median = (run_.results.surface[0].d_median + run_.results.surface[1].d_median) / 2;
		run_.results.surface[2].s_mean = (run_.results.surface[0].s_mean + run_.results.surface[1].s_mean) / 2;
		run_.results.surface[2].s_median = (run_.results.surface[0].s_median + run_.results.surface[1].s_median) / 2;

		run_.results.surface[2].nAtoms = (run_.results.surface[0].nAtoms + run_.results.surface[1].nAtoms);
		run_.results.surface[2].nBuriedAtoms = (run_.results.surface[0].nBuriedAtoms + run_.results.surface[1].nBlockedAtoms);
		run_.results.surface[2].nBlockedAtoms = (run_.results.surface[0].nBuriedAtoms + run_.results.surface[1].nBuriedAtoms);
		run_.results.surface[2].nAllDots = (run_.results.surface[0].nAllDots + run_.results.surface[1].nAllDots);
		run_.results.surface[2].nTrimmedDots = (run_.results.surface[0].nTrimmedDots + run_.results.surface[1].nTrimmedDots);
		run_.results.surface[2].trimmedArea = (run_.results.surface[0].trimmedArea + run_.results.surface[1].trimmedArea);

		run_.results.sc = run_.results.surface[2].s_median;
		run_.results.distance = run_.results.surface[2].d_median;
		run_.results.area = run_.results.surface[2].trimmedArea;
		run_.results.valid = 1;

		TR.Debug <<
			"Done. Atoms: " << run_.results.surface[0].nAtoms << " + " << run_.results.surface[1].nAtoms <<
			"; sc = " << run_.results.sc <<
			", area = " << run_.results.area <<
			", distance " << run_.results.distance <<
			std::endl;

		return 1;
	} catch ( ShapeComplementarityCalculatorException & e ) {
		TR.Error << "Failed: " << e.error << std::endl;
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////
// Protected class functions
////////////////////////////////////////////////////////////////////////////

// Determine assign the attention numbers for each atom
int ShapeComplementarityCalculator::AssignAttentionNumbers(std::vector<Atom> & )
{
	std::vector<Atom>::iterator pAtom1, pAtom2;

	for ( pAtom1 = run_.atoms.begin(); pAtom1 < run_.atoms.end(); ++pAtom1 ) {
		// find nearest neighbour in other molecule
		ScValue dist_min = 99999.0, r;
		for ( pAtom2 = run_.atoms.begin(); pAtom2 < run_.atoms.end(); ++pAtom2 ) {
			if ( pAtom1->molecule == pAtom2->molecule ) {
				continue;
			}
			r = pAtom1->distance(*pAtom2);
			if ( r < dist_min ) {
				dist_min = r;
			}
		}

		// check if within separator distance
		if ( dist_min >= settings.sep ) {
			// TR.Debug << "Atom ATTEN_BLOCKER: " << pAtom1->natom << std::endl;
			// too _far_ away from other molecule, blocker atom only
			pAtom1->atten = ATTEN_BLOCKER;
			++run_.results.surface[pAtom1->molecule].nBlockedAtoms;
		} else {
			// potential interface or neighbouring atom
			pAtom1->atten = ATTEN_BURIED_FLAGGED;
			++run_.results.surface[pAtom1->molecule].nBuriedAtoms;
		}
	}

	return 1;
}

// SC molecular dot trimming, vector dot product calculation and statistics
// Trim dots and retain only the peripheral band

ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::TrimPeripheralBand(
	std::vector<DOT> const &sdots,
	std::vector<DOT const *> &trimmed_dots)
{
	ScValue area = 0;

	if ( sdots.empty() ) return 0.0;

#ifdef USEOPENCL
	if(settings.gpu) {
		area = gpuTrimPeripheralBand(sdots, trimmed_dots);
	} else {
#endif

	// Loop over one surface
	// If a point is buried then see if there is an accessible point within distance band

	for ( std::vector<DOT>::const_iterator idot = sdots.begin(); idot < sdots.end(); ++idot ) {
		DOT const &dot = *idot;
		// Paralelleizable kernel function
		if ( dot.buried && TrimPeripheralBandCheckDot(dot, sdots) ) {
			area += dot.area;
			trimmed_dots.push_back(&dot);
		}
	}

#ifdef USEOPENCL
	}
#endif

	return area;
}

// Test a dot against a set of dots for collision
// NOTE: ~75% of time is spent in this function
int ShapeComplementarityCalculator::TrimPeripheralBandCheckDot(
	DOT const &dot,
	std::vector<DOT> const &sdots)
{
	// Caching of r2 only brings 0.5% speed boost
	ScValue r2 = pow(settings.band, 2);

	for ( std::vector<DOT>::const_iterator idot2 = sdots.begin(); idot2 < sdots.end(); ++idot2 ) {
		DOT const &dot2 = *idot2;
		if ( &dot == &dot2 ) continue;
		if ( dot2.buried ) continue;
		if ( dot.coor.distance_squared(dot2.coor) <= r2 ) return 0;
	}
	return 1;
}

////////////////////////////////////////////////////////////////////////////
// Calculate separation distance and and normal vector dot product (shape)
// distributions, mean and median of molecular surface
////////////////////////////////////////////////////////////////////////////

int ShapeComplementarityCalculator::CalcNeighborDistance(
	int const molecule,
	std::vector<DOT const*> const &my_dots,
	std::vector<DOT const*> const &their_dots)
{
	std::map<int,int> dbins; // Distance bins
	std::map<int,int> sbins; // Vector dot product bins (sc)
	ScValue norm_sum = 0.0, distmin_sum = 0.0;
	int ibin;
	ScValue total = 0.0;

	if ( my_dots.empty() || their_dots.empty() ) return 0;

	for ( std::vector<DOT const*>::const_iterator idot = my_dots.begin();
			idot < my_dots.end(); ++idot ) {
		total += (*idot)->area;
	}

#ifdef USEOPENCL
	std::vector<DOT const*> neighbors;
	std::vector<DOT const*>::const_iterator iNeighbor;

	if(settings.gpu) {
		gpuFindClosestNeighbors(my_dots, their_dots, neighbors);
		iNeighbor = neighbors.begin();
	}
#endif

	for ( std::vector<DOT const*>::const_iterator idot = my_dots.begin();
			idot < my_dots.end(); ++idot ) {
		DOT const &dot1 = **idot;

		ScValue distmin, r;
		DOT const *neighbor = NULL;

#ifdef USEOPENCL
		if ( settings.gpu ) {
			neighbor = *iNeighbor++;
		} else {
			neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);
		}
#else
		neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);
#endif

		if ( !neighbor ) continue;

		// having looked at all possible neighbours now accumulate stats
		distmin = neighbor->coor.distance(dot1.coor);
		distmin_sum += distmin;
		// decide which bin to put it into and then add to distance histogram
		ibin = (int)(distmin / settings.binwidth_dist);
		++dbins[ibin];

		//work out dot product
		r = dot1.outnml.dot(neighbor->outnml);

		// weight dot product
		// cpjx I think the weighting factor is the denominator 2 below?
		// cpjx    r = r * exp( - (distmin**2) / 2.)
		r = r * exp( - pow(distmin, 2) * settings.weight );
		// rounding errors a problem, so ensure std::abs(r) <1
		r = MIN(0.999, MAX(r, -0.999));
		norm_sum += r;

		// left_trunc ScValue to int ibin
		// otherwise: (int)-0.9 = 0.
		r /= settings.binwidth_norm;
		if ( r >= 0 ) {
			ibin = (int)r;
		} else {
			ibin = (int)r -1;
		}
		++sbins[ibin];
	}

	// Determine the last distance bin that has anything in it
	// Accumulate percentages and area from all filled distance bins
	ScValue abin, cumarea =0, cumperc = 0, perc, c;
	ScValue rleft =0, rmedian =0;
	std::map<int,int>::const_iterator it;

	TR.Trace << std::endl;
	TR.Trace << "Distance between surfaces D(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl;
	TR.Trace << "From - To\tArea\tCum. Area\t%\tCum. %" << std::endl;

	for ( it = dbins.begin(); it != dbins.end(); ++it ) {
		abin = total * (it->second) / my_dots.size();
		cumarea += abin;
		perc = abin * 100 / total;
		c = cumperc + perc;
		if ( cumperc <= 50 && c >= 50 ) {
			rleft = (it->first) * settings.binwidth_dist;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_dist / ( c - cumperc );
		}
		cumperc = c;

#ifndef WIN32
		if ( TR.Trace.visible() ) {
			char buf[128];

			snprintf(buf, sizeof(buf),
				"%.2f - %.2f\t%.1f\t%.1f\t%.1f\t%.1f",
				(ScValue)it->first * settings.binwidth_dist,
				(ScValue)it->first * settings.binwidth_dist + settings.binwidth_dist,
				abin, cumarea,
				perc, cumperc);

			TR.Trace << buf << std::endl;
		}
#endif
	}

	run_.results.surface[molecule].d_mean = distmin_sum / my_dots.size();
	run_.results.surface[molecule].d_median = rmedian;

	TR.Trace << std::endl;
	TR.Trace << "Surface complementarity S(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl;
	TR.Trace << "From - To\tNumber\t%\tCumm. %" << std::endl;

	cumperc = 0;
	for ( it = sbins.begin(); it != sbins.end(); ++it ) {
		perc = (ScValue)(it->second) * 100 / my_dots.size();
		c = cumperc + perc;
		if ( cumperc <= 50 && c >= 50 ) {
			rleft = (ScValue)(it->first) * settings.binwidth_norm;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_norm / ( c - cumperc );
		}
		cumperc = c;

#ifndef WIN32
		if ( TR.Trace.visible() ) {
			char buf[128];
			snprintf(buf, sizeof(buf),
				"%.2f - %.2f\t%d\t%.1f\t%.1f",
				(ScValue)-it->first * settings.binwidth_norm - settings.binwidth_norm,
				(ScValue)-it->first * settings.binwidth_norm,
				it->second, perc, cumperc);
			TR.Trace << buf << std::endl;
		}
#endif
	}
	run_.results.surface[molecule].s_mean= -norm_sum / my_dots.size();
	run_.results.surface[molecule].s_median = -rmedian;

	return 1;
}

// Find closest neighbor dot for a given dot
// NOTE: ~20% of time is spent in this function
DOT const *ShapeComplementarityCalculator::CalcNeighborDistanceFindClosestNeighbor(
	DOT const &dot1,
	std::vector<DOT const*> const &their_dots
) {
	ScValue distmin = 999999.0, d;
	DOT const *neighbor = NULL;

	// Loop over the entire surface: find and flag neighbour of each point
	// that we're interested in and store nearest neighbour pointer

	for ( std::vector<DOT const*>::const_iterator idot2 = their_dots.begin();
			idot2 < their_dots.end(); ++idot2 ) {
		DOT const &dot2 = **idot2;
		if ( !dot2.buried ) continue;
		d = dot2.coor.distance_squared(dot1.coor);
		if ( d <= distmin ) {
			distmin = d;
			neighbor = *idot2;
		}
	}
	return neighbor;
}

////////////////////////////////////////////////////////////////////////
// GPU SUPPORT FUNCTIONS

#ifdef USEOPENCL

core::Real inline ShapeComplementarityCalculator::GetTimerMs(clock_t &start)
{
	clock_t now = clock();
	core::Real d = (now - start)/(CLOCKS_PER_SEC/1000);
	return d;
}

void ShapeComplementarityCalculator::gpuInit()
{
	if(gpu.use()) {
	  if(TR.Debug.visible())
		  gpu.profiling(1);
		if(gpu.Init())
			settings.gpu = 1;
	}
	if(settings.gpu_threads < 32)
		settings.gpu_threads = gpu.device().threads;
	gpu.RegisterProgram("gpu/sc.cl");
}

ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::gpuTrimPeripheralBand(
		std::vector<DOT> const &dots,
		std::vector<DOT const*> &trimmed_dots)
{
	using namespace basic::gpu;

	int n, nBur, nAcc;
	int threads;
	ScValue area = 0;
	clock_t timer;

	threads = MIN(512, settings.gpu_threads);
	n = dots.size();
	timer = clock();

	// Host and device (GPU) memory pointers for dot coordinates and results
	float4 *hAccDotCoords, *phAccDotCoords;
	float4 *hBurDotCoords, *phBurDotCoords;
	char *hDotColl;

	hAccDotCoords = new float4[UPPER_MULTIPLE(n, threads)];
	hBurDotCoords = new float4[UPPER_MULTIPLE(n, threads)];
	hDotColl = new char[UPPER_MULTIPLE(n, threads)];

	if(!hAccDotCoords || !hBurDotCoords || !hDotColl) {
		throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "Out of host memory!");
	}

	// Make GPU copy of (x, y, z) buried and accessible coordinates
	phAccDotCoords = hAccDotCoords;
	phBurDotCoords = hBurDotCoords;
	for(std::vector<DOT>::const_iterator idot = dots.begin();
			idot < dots.end(); ++idot) {
#ifdef SC_PRECISION_REAL
		if(idot->buried) {
			phBurDotCoords->x = idot->coor.x();
			phBurDotCoords->y = idot->coor.y();
			phBurDotCoords->z = idot->coor.z();
			++phBurDotCoords;
		} else {
			phAccDotCoords->x = idot->coor.x();
			phAccDotCoords->y = idot->coor.y();
			phAccDotCoords->z = idot->coor.z();
			++phAccDotCoords++;
		}
#else
		// Quick copy
		if(idot->buried)
			*phBurDotCoords++ = *((float4*)&idot->coor.x());
		else
			*phAccDotCoords++ = *((float4*)&idot->coor.x());
#endif
	}
	nBur = phBurDotCoords - hBurDotCoords;
	nAcc = phAccDotCoords - hAccDotCoords;

	// Run kernel on GPU
	float r2 = pow(settings.band, 2);
	if(!gpu.ExecuteKernel("TrimPeripheralBand", nBur, threads, 32,
		GPU_IN, UPPER_MULTIPLE(nAcc, threads) * sizeof(*hAccDotCoords), hAccDotCoords,
		GPU_INT, nAcc,
		GPU_IN, UPPER_MULTIPLE(nBur, threads) * sizeof(*hBurDotCoords), hBurDotCoords,
		GPU_OUT, UPPER_MULTIPLE(nBur, threads) * sizeof(*hDotColl), hDotColl,
		GPU_FLOAT, r2,
		NULL)) {
		throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "Failed to launch GPU kernel TrimPeripheralBand!");
	}

	// Make a new list of dots that have no collisions
	char *p = hDotColl;
	for(std::vector<DOT>::const_iterator idot = dots.begin();
			idot < dots.end(); ++idot) {
		DOT const &dot1 = *idot;
		if(!idot->buried)
			continue;
		if(!*p++) {
			area += dot1.area;
			trimmed_dots.push_back(&dot1);
		}
	}

	delete [] hAccDotCoords;
	delete [] hBurDotCoords;
	delete [] hDotColl;

	TR.Debug << "Peripheral trimming GPU processing time: " << gpu.lastKernelRuntime() << " ms kernel, " << GetTimerMs(timer) << " ms total" << std::endl;

	return area;
}

int ShapeComplementarityCalculator::gpuFindClosestNeighbors(
	std::vector<DOT const*> const &my_dots,
	std::vector<DOT const*> const &their_dots,
	std::vector<DOT const*> &neighbors)
{
	using namespace basic::gpu;

	int nMyDots, nTheirDots, nNeighbors;
	int threads;
	clock_t timer;

	timer = clock();
	threads = MIN(512, settings.gpu_threads);

	// Memory pointers for my and their dot coordinate arrays, CPU and GPU
	float4 *hMyDotCoords, *phMyDotCoords;
	float4 *hTheirDotCoords, *phTheirDotCoords;

	// Dot point pointer map
	DOT const **hTheirDots, **phTheirDots;

	// Neighbor ID memory pointers
	::uint *hNeighbors;

	nMyDots = my_dots.size();
	nTheirDots = their_dots.size();
	nNeighbors = nMyDots;

	hMyDotCoords = new float4 [UPPER_MULTIPLE(nMyDots, threads)];
	hTheirDotCoords = new float4 [UPPER_MULTIPLE(nTheirDots, threads)];
	hTheirDots = new DOT const * [UPPER_MULTIPLE(nTheirDots, threads)];
	hNeighbors = new ::uint[UPPER_MULTIPLE(nNeighbors, threads)];

	// Make GPU copy of (x, y, z) dot coordinates for my dots
	phMyDotCoords = hMyDotCoords;
	for(std::vector<DOT const*>::const_iterator idot = my_dots.begin(); idot < my_dots.end(); ++idot) {
		phMyDotCoords->x = (*idot)->coor.x();
		phMyDotCoords->y = (*idot)->coor.y();
		phMyDotCoords->z = (*idot)->coor.z();
		++phMyDotCoords;
	}
	nMyDots = phMyDotCoords - hMyDotCoords;

	// Make GPU copy of (x, y, z) dot coordinates for their dots and keep a map
	phTheirDotCoords = hTheirDotCoords;
	phTheirDots = hTheirDots;
	for(std::vector<DOT const*>::const_iterator idot = their_dots.begin(); idot < their_dots.end(); ++idot) {
		if(!(*idot)->buried)
			continue;
		phTheirDotCoords->x = (*idot)->coor.x();
		phTheirDotCoords->y = (*idot)->coor.y();
		phTheirDotCoords->z = (*idot)->coor.z();
		++phTheirDotCoords;
		*phTheirDots++ = *idot;
	}
	nTheirDots = phTheirDotCoords - hTheirDotCoords;

	// Run kernel on GPU
	if(!gpu.ExecuteKernel("FindClosestNeighbor", nMyDots, threads, 32,
		GPU_IN, UPPER_MULTIPLE(nMyDots, threads) * sizeof(*hMyDotCoords), hMyDotCoords,
		GPU_IN, UPPER_MULTIPLE(nTheirDots, threads) * sizeof(*hTheirDotCoords), hTheirDotCoords,
		GPU_INT, nTheirDots,
		GPU_OUT, UPPER_MULTIPLE(nNeighbors, threads) * sizeof(*hNeighbors), hNeighbors,
		NULL)) {
		throw CREATE_EXCEPTION(ShapeComplementarityCalculatorException, "Failed to launch GPU kernel FindClosestNeighbor!");
	}

	for(int i = 0; i < nNeighbors; ++i)
		neighbors.push_back( hTheirDots[hNeighbors[i]] );

	delete [] hMyDotCoords;
	delete [] hTheirDotCoords;
	delete [] hTheirDots;
	delete [] hNeighbors;

	TR.Debug << "Find Neighbors GPU processing time: " << gpu.lastKernelRuntime() << " ms kernel, " << GetTimerMs(timer) << " ms total" << std::endl;

	return 1;
}

#endif // USEOPENCL

// The End
////////////////////////////////////////////////////////////////////////////

} // namespace sc
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_cc

// END //
