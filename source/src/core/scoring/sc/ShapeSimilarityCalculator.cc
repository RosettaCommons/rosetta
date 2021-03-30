// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/scoring/sc/ShapeSimilarityCalculator.cc
/// @brief    Headers for the Shape Similarity Calculator
/// @details  The code modifies Luki Goldschmidt's
///           implementation of Lawrence & Coleman shape complementarity calculator
///           to allow for the comparison of similar surface shapes.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)


#ifndef INCLUDED_core_scoring_sc_ShapeSimilarityCalculator_cc
#define INCLUDED_core_scoring_sc_ShapeSimilarityCalculator_cc

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/scoring/sc/ShapeSimilarityCalculator.hh>
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
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueVector.hh>



// C headers
#include <cstdio>

// C++ headers
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip> // for std::setprecision()

#define UPPER_MULTIPLE(n,d) (((n)%(d)) ? (((n)/(d)+1)*(d)) : (n))

static basic::Tracer TR( "core.scoring.sc.ShapeSimilarityCalculator" );

using namespace core;

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////////////////////
// Public class functions
////////////////////////////////////////////////////////////////////////////

/// @brief
/// ShapeSimilarityCalculator constructor, initializes default settings

ShapeSimilarityCalculator::ShapeSimilarityCalculator() :
	MolecularSurfaceCalculator(),
	median( false )
{}

ShapeSimilarityCalculator::~ShapeSimilarityCalculator() = default;

/// @brief
/// Run the SS calculation on Pose and Reference and return the ss statistic or -1 on error
/// @details
/// This is a static function and can be called without instantiating ShapeSimilarityCalculator.
/// The residue selectors for the pose and reference are used to select two surfaces to compare.
/// Setting quick to true will perform a much faster calculation (~5-10 times faster) at the expense
/// of accuracy (about 0.05 units).
///
/// Example:
/// core::Real sc = core::scoring::sc::ShapeSimilarityCalculator( pose );
core::Real ShapeSimilarityCalculator::CalcSs(
	core::pose::Pose const & pose,
	core::pose::Pose const & native,
	core::select::residue_selector::ResidueSelectorCOP selector_pose,
	core::select::residue_selector::ResidueSelectorCOP selector_native ) {

	using core::select::residue_selector::ResidueVector;

	ShapeSimilarityCalculator ssc;
	ShapeSimilarityCalculator ssc_full;
	ShapeSimilarityCalculator ssc_ref;
	ShapeSimilarityCalculator ssc_ref_full;


	ResidueVector const residues_pose( selector_pose->apply( pose ) );
	ResidueVector const residues_native( selector_native->apply( native ) );

	// this assumes all residues are from the same chain;
	Size pose_chain = pose.residue(residues_pose[1]).chain();
	Size native_chain = native.residue(residues_native[1]).chain();

	TR << "chain pose: " << pose_chain << " and chain native: " << native_chain << std::endl;


	// Dump information about residues
	TR << "Using residues for molecule surface (rosetta numbering):" << std::endl;
	TR << "  Scaffold: ";
	for ( auto r = residues_pose.begin(); r != residues_pose.end(); ++r ) {
		TR << (r == residues_pose.begin() ? "" : ", ") << *r;
	}
	TR << std::endl;
	TR << "  Native: ";
	for ( auto r = residues_native.begin(); r != residues_native.end(); ++r ) {
		TR << (r == residues_native.begin() ? "" : ", ") << *r;
	}
	TR << std::endl;


	for ( core::Size r : residues_pose ) {
		ssc.AddResidue( 0, pose.residue( r ) );
	}

	core::pose::PoseOP pose_subset = pose.split_by_chain(pose_chain);
	for ( auto r = pose_subset->begin(); r != pose_subset->end(); ++r ) {
		ssc_full.AddResidue( 0, *r);
	}

	for ( core::Size r : residues_native ) {
		ssc_ref.AddResidue( 0, native.residue( r ) );
	}

	core::pose::PoseOP native_subset = native.split_by_chain(native_chain);
	for ( auto r = native_subset->begin(); r != native_subset->end(); ++r ) {
		ssc_ref_full.AddResidue( 0, *r);
	}


	if ( !ssc.Calc() ) {
		return -1;
	}
	if ( !ssc_full.Calc() ) {
		return -1;
	}
	if ( !ssc_ref.Calc() ) {
		return -1;
	}
	if ( !ssc_ref_full.Calc() ) {
		return -1;
	}


	// To get all dots from the selected surface, we first compute the overall surface of the protein
	// followed by the surface of the selected residues, for the target and reference respectively.
	// The two dot clouds of the full surface and selected residues's surface are compared and only
	// dots that intersect are kept. This is to avoid dots that are an artefact of the surface
	// generation of selected residues and are actually buried in the protein core.
	std::vector<DOT const *> intersecting_dots[1];
	std::vector<DOT const *> intersecting_dots_native[1];

	for ( auto r1 = ssc.GetDots(0).begin(); r1 != ssc.GetDots(0).end(); ++r1 ) {
		for ( auto r2 = ssc_full.GetDots(0).begin(); r2 != ssc_full.GetDots(0).end(); ++r2 ) {
			core::Real dist = r1->coor.distance(r2->coor);
			if ( dist == 0 ) {
				DOT const &dot1 = *r1;
				intersecting_dots[0].push_back(&dot1);
				break;
			}
		}
	}

	for ( auto r1 = ssc_ref.GetDots(0).begin(); r1 != ssc_ref.GetDots(0).end(); ++r1 ) {
		for ( auto r2 = ssc_ref_full.GetDots(0).begin(); r2 != ssc_ref_full.GetDots(0).end(); ++r2 ) {
			core::Real dist = r1->coor.distance(r2->coor);
			if ( dist == 0 ) {
				DOT const &dot2 = *r1;
				intersecting_dots_native[0].push_back(&dot2);
				break;
			}
		}
	}

	sort( intersecting_dots[0].begin(), intersecting_dots[0].end() );
	sort( intersecting_dots_native[0].begin(), intersecting_dots_native[0].end() );

	intersecting_dots[0].erase( std::unique( intersecting_dots[0].begin(), intersecting_dots[0].end() ), intersecting_dots[0].end() );
	intersecting_dots_native[0].erase( std::unique( intersecting_dots_native[0].begin(), intersecting_dots_native[0].end() ), intersecting_dots_native[0].end() );

	core::Real shape_score1 = CalcNeighborDistance(0, intersecting_dots[0], intersecting_dots_native[0]);
	core::Real shape_score2 = CalcNeighborDistance(1, intersecting_dots_native[0], intersecting_dots[0]);

	core::Real ss_value = (shape_score1 + shape_score2)/2;

	TR << "Finishing SS calculations!" << std::endl;
	TR << shape_score1 << ", " << shape_score2 << std::endl;
	TR << "ss = " << ss_value << std::endl;
	return -ss_value;
}


/// @brief Run the SS calculation for previously defined molecules (via AddResidue or add_atom calls)
/// @details
/// This function is called for residues / atoms making up the two molecular surfaces
/// have been explicitly defined.
/// Returns true on success.

int ShapeSimilarityCalculator::Calc()
{
#ifdef USEOPENCL
	gpuInit();
#endif

	try {
		basic::gpu::Timer timer(TR.Debug);

		run_.results.valid = 0;

		if ( run_.atoms.empty() ) {

			throw CREATE_EXCEPTION(utility::excn::Exception, "Error in ShapeSimilarityCalculator::Calc(): No atoms defined");
		}
		if ( !run_.results.surface[0].nAtoms ) {
			throw CREATE_EXCEPTION(utility::excn::Exception, "Error in ShapeSimilarityCalculator::Calc(): No atoms defined for molecule 1");
		}

		// Determine and assign the attention numbers for each atom
		AssignAttentionNumbers(run_.atoms); // repalce this with a call to MolecularSurfaceCalculator::AssignAttentinoNumbers

		GenerateMolecularSurfaces();
		for ( auto r = run_.dots[0].begin(); r != run_.dots[0].end(); ++r ) {
			r->buried = 1;
		}

		return 1;

	} catch ( utility::excn::Exception & e ) {
		TR.Error << "Failed: " << e.msg() << std::endl;
	}

	return -1;
}

////////////////////////////////////////////////////////////////////////////
// Protected class functions
////////////////////////////////////////////////////////////////////////////

int ShapeSimilarityCalculator::AssignAttentionNumbers(std::vector<Atom> &)
{
	std::vector<Atom>::iterator pAtom;

	for ( pAtom = run_.atoms.begin(); pAtom < run_.atoms.end(); ++pAtom ) {
		pAtom->atten = ATTEN_BURIED_FLAGGED;
		++run_.results.surface[pAtom->molecule].nBuriedAtoms;
	}

	return 1;
}


/// @brief SS molecular dot trimming, vector dot product calculation and statistics
/// @details Trim dots and retain only the peripheral band (currently not used during evaluation)
ShapeSimilarityCalculator::ScValue ShapeSimilarityCalculator::TrimPeripheralBand(
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

	for ( auto idot = sdots.begin(); idot < sdots.end(); ++idot ) {
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

/// @brief Test a dot against a set of dots for collision
/// NOTE: ~75% of time is spent in this function
int ShapeSimilarityCalculator::TrimPeripheralBandCheckDot(
	DOT const &dot,
	std::vector<DOT> const &sdots)
{
	// Caching of r2 only brings 0.5% speed boost
	ScValue r2 = pow(settings.band, 2);
	//ScValue r2 = pow(0.5, 2);

	for ( auto idot2 = sdots.begin(); idot2 < sdots.end(); ++idot2 ) {
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

/// @brief Evaluate similarity of two points of the surface
/// @details
/// Identify the closes points of the two provided surfaces and compute the distance between the points and their normal vectors.
/// The normal vectors are used to evaluate the shape of the points via the normal vector dot product.
/// This operation is performed for all points of the specified surfaces and the mean of all values represent the overall shape similarity.
core::Real ShapeSimilarityCalculator::CalcNeighborDistance(
	int const molecule,
	std::vector<DOT const*> const &my_dots,
	std::vector<DOT const*> const &their_dots)
{
	std::map<int,int> dbins; // Distance bins
	std::map<int,int> sbins; // Vector dot product bins (ss)
	ScValue norm_sum = 0.0, distmin_sum = 0.0;
	int ibin;
	ScValue total = 0.0;

	// Count the dots that have been ejected due to having r = 0 and distmin = 0
	Size ejected = 0;
	if ( my_dots.empty() || their_dots.empty() ) return 0;

	for ( auto idot = my_dots.begin(); idot < my_dots.end(); ++idot ) {
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

	for ( auto idot = my_dots.begin(); idot < my_dots.end(); ++idot ) {
		DOT const &dot1 = **idot;

		ScValue distmin, r;

#ifdef USEOPENCL
		DOT const *neighbor = nullptr;
		if ( settings.gpu ) {
			neighbor = *iNeighbor++;
		} else {
			neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);
		}
#else
		DOT const *neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);
#endif

		if ( !neighbor ) continue;

		// having looked at all possible neighbours now accumulate stats
		distmin = neighbor->coor.distance(dot1.coor);

		//work out dot product
		r = dot1.outnml.dot(neighbor->outnml);

		// weight dot product
		// cpjx I think the weighting factor is the denominator 2 below?
		// cpjx    r = r * exp( - (distmin**2) / 2.)
		r = r * exp( - pow(distmin, 2) * settings.weight );
		// rounding errors a problem, so ensure std::abs(r) <1
		r = MIN(0.999, MAX(r, -0.999));

		// workaround to cancel out low scoring duplicates
		// There are some dots that have distances of 0 but get a score of 0, which are also not filtered when duplicated dots are removed.
		if ( r == 0.0 && distmin == 0.0 ) {
			ejected += 1;
			continue;
		}
		distmin_sum += distmin;
		// decide which bin to put it into and then add to distance histogram
		ibin = (int)(distmin / settings.binwidth_dist);
		++dbins[ibin];
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
		abin = total * (it->second) / (my_dots.size() - ejected);
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

	run_.results.surface[molecule].d_mean = distmin_sum / (my_dots.size() - ejected);
	run_.results.surface[molecule].d_median = rmedian;

	TR.Trace << std::endl;
	TR.Trace << "Surface similarity S(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl;
	TR.Trace << "From - To\tNumber\t%\tCumm. %" << std::endl;

	cumperc = 0;
	for ( it = sbins.begin(); it != sbins.end(); ++it ) {
		perc = (ScValue)(it->second) * 100 / (my_dots.size() - ejected);
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
	run_.results.surface[molecule].s_mean= -norm_sum / (my_dots.size() - ejected);
	run_.results.surface[molecule].s_median = -rmedian;
	TR << "ss mean: " << -norm_sum / (my_dots.size() - ejected) << std::endl;
	TR << "ss median: " << -rmedian << std::endl;

	if ( median ) {
		return -rmedian;
	} else {
		return -norm_sum / (my_dots.size() - ejected);
	}
}

/// @brief Find closest neighbor dot for a given dot
/// NOTE: ~20% of time is spent in this function
DOT const *ShapeSimilarityCalculator::CalcNeighborDistanceFindClosestNeighbor(
	DOT const &dot1,
	std::vector<DOT const*> const &their_dots
) {
	ScValue distmin = 999999.0, d;
	DOT const *neighbor = nullptr;

	// Loop over the entire surface: find and flag neighbour of each point
	// that we're interested in and store nearest neighbour pointer

	for ( auto idot2 = their_dots.begin();
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

core::Real inline ShapeSimilarityCalculator::GetTimerMs(clock_t &start)
{
	clock_t now = clock();
	core::Real d = (now - start)/(CLOCKS_PER_SEC/1000);
	return d;
}

void ShapeSimilarityCalculator::gpuInit()
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

ShapeSimilarityCalculator::ScValue ShapeSimilarityCalculator::gpuTrimPeripheralBand(
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
		throw CREATE_EXCEPTION(utility::excn::Exception, "Error in ShapeSimilarityCalculator::gpuTrimPeripheralBand(): Out of host memory!");
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
		throw CREATE_EXCEPTION(utility::excn::Exception, "Error in ShapeSimilarityCalculator::gpuTrimPeripheralBand: Failed to launch GPU kernel TrimPeripheralBand!");
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

int ShapeSimilarityCalculator::gpuFindClosestNeighbors(
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
		throw CREATE_EXCEPTION(utility::excn::Exception, "Error in ShapeSimilarityCalculator::gpuFindClosestNeighbor(): Failed to launch GPU kernel FindClosestNeighbor!");
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

#endif // INCLUDED_core_scoring_sc_ShapeSimilarityCalculator_cc
