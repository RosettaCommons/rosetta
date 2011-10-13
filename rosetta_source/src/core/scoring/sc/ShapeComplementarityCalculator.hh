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

// This code contains support for GPU acceleration using CUDA.
// Basic steps to enable GPU acceleration support:
//  1) define USECUDA (-DUSECUDA)
//  2) compile ShapeComplementarityCalculator_GPUKernels.cu with nvcc
//  3) link libcore.3.so file against libcudart.so (-lcudart)
// You will see a line like this when GPU is utilized:
//  core.scoring.sc.ShapeComplementarityCalculator: GPU support enabled:
//  GeForce GTX 580 [1594 MHz, capability 2.0] with 16 multi processors, 1024 threads.
//
// Note: Original Fotran code used floats rather than doubles. This code works
// correctly with either floats or core::Real, though core::Real is about 50%
// slower (tested on 64-bit machines). The code can be compiled with double (Real)
// precision floats by defining SC_PRECISION_REAL.

// #define SC_PRECISION_REAL
// #define USECUDA

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>
#include <numeric/xyzVector.hh>

//// C++ headers
#include <vector>
#include <map>
#include <string>

#ifdef USECUDA
#include <time.h>
#endif

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////
// core::scoring::sc namespace
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Types

// Calculation results
typedef struct _RESULTS {
	core::Real sc;                 // Interface shape complementarity statistic
	core::Real area;               // Area of interface (A^2)
	core::Size nAtoms;             // Number of Atoms used in calculation
	struct {                       // These detailed counts are for for both molecules (0, 1) and sum/average for both (2)
		core::Real d_mean;           // Mean seperation of molecules
		core::Real d_median;         // Median separation of molecules
		core::Real s_mean;           // Mean shape complementarity
		core::Real s_median;         // Median shape complementarity (this is sc above)
		core::Size nAtoms;           // Number of atoms in molecule
		core::Size nBuriedAtoms;     // Number of buried atoms in molecule
		core::Size nBlockedAtoms;    // Number of blocked (covered) atoms in molecule
		core::Size nAllDots;         // Number of molecule surface all dots
		core::Size nTrimmedDots;     // Number of molecule surface dots after trimming
		core::Size nBuriedDots;      // Number of buried dots (not on surface)
		core::Size nAccessibleDots;  // Number of accessible dots (on surface)
		core::Real trimmedArea;      // Trimmed area in Angstrom^2
	} surface[3];
	struct {			// Surface dot counts by type (from Connolly algorithm)
		core::Size convex;           // Number of convex surface dots
		core::Size concave;          // Number of concace surface dots
		core::Size toroidal;         // NUmber of toroidal surfac dots
	} dots;				// True if computed results are valid
	int valid;			
} RESULTS;

// Atom radius definition
typedef struct _ATOM_RADIUS {
	char residue[4];
	char atom[5];
#ifdef SC_PRECISION_REAL
	core::Real radius;
#else
	float radius;
#endif
} ATOM_RADIUS;

// Molecular dot
class Atom;
typedef struct _DOT {
#ifdef SC_PRECISION_REAL
	numeric::xyzVector < core::Real > coor, outnml;
	core::Real area;
#else
	numeric::xyzVector < float > coor, outnml;
	float area;
#endif
	int buried;
	int type;
	Atom const *atom;
} DOT;

// Molecular probe
typedef struct _PROBE {
	Atom const *pAtoms[3];
#ifdef SC_PRECISION_REAL
	core::Real height;
	numeric::xyzVector < core::Real > point, alt;
#else
	float height;
	numeric::xyzVector < float > point, alt;
#endif
} PROBE;

// PDB Atom
class Atom :
#ifdef SC_PRECISION_REAL
	public numeric::xyzVector < core::Real >
#else
	public numeric::xyzVector < float >
#endif
{
protected:
	core::Size natom;
	core::Size nresidue;
	char atom[4];
	char residue[4];

	friend class ShapeComplementarityCalculator;

	core::Size molecule;
#ifdef SC_PRECISION_REAL
	core::Real radius;
	core::Real density;
#else
	float radius;
	float density;
#endif
	int atten;
	int access;
	std::vector<Atom*> neighbors;
	std::vector<Atom*> buried;

public:
	Atom();
	~Atom();

	int operator ==(Atom const &atom2) {
		// <EVIL>
		// For speed reasons, rather than comparing the coordinates,
		// we'll compare pointers as all atoms are part of the same atoms
		// vector, so the same atom will have the same address
		// </EVIL>
		return this == &atom2;
	}
	int operator <=(Atom &atom2) {
		return this <= &atom2;
	}
};

////////////////////////////////////////////////////////////
// Shape Complementarity Calculator class definition
////////////////////////////////////////////////////////////

class ShapeComplementarityCalculator {

public:
#ifdef SC_PRECISION_REAL
	typedef core::Real ScValue;
	typedef numeric::xyzVector < core::Real > Vec3;
#else
	typedef float ScValue;
	typedef numeric::xyzVector < float > Vec3;
#endif

	struct {
		// From sc source
		core::Real rp;
		core::Real density;
		core::Real band;
		core::Real sep;
		core::Real weight;
		core::Real binwidth_dist;
		core::Real binwidth_norm;

		core::Size verbose;

#ifdef USECUDA
		core::SSize gpu;
		core::Size gpu_threads;
		core::Size gpu_proc;
#endif
	} settings;

	ShapeComplementarityCalculator();
	~ShapeComplementarityCalculator();
	int Init();
	void Reset();

	int AddAtom(int molecule, Atom &atom);
	core::Size AddResidue(int molecule, core::conformation::Residue const &residue);

	int Calc();
	int Calc(core::pose::Pose const & pose, core::Size jump_id);

	// Static function that can be called without object instantiation
	static core::Real CalcSc(core::pose::Pose const & pose, core::Size jump_id =1, int quick =0);

	RESULTS const & GetResults() { return run_.results; }
	std::vector<Atom> const & GetAtoms() { return run_.atoms; }

protected:
	// This is a constant list of atom radii; declared static to avoid re-loading
	// between computations
	static std::vector<ATOM_RADIUS> radii_;

	int AssignAtomRadius(Atom &atom);
	int WildcardMatch(char const *r, char const *pattern, int const l);
	int ReadScRadii();
	void AddDot(int const molecule, int const type, Vec3 const coor, ScValue const area, Vec3 const pcen, Atom const &atom);

private:
	struct {
		ScValue radmax;
		RESULTS results;
		std::vector<Atom> atoms;
		std::vector<DOT> dots[2];
		std::vector<PROBE> probes;
		Vec3 prevp;
		int prevburied;
	} run_;

	// Molecular surface generation
	int AssignAttentionNumbers(std::vector<Atom>& atom);
	int CalcDotsForAllAtoms(std::vector<Atom>& atoms);
	int CalcDotsForAtoms(std::vector<Atom>& atoms);
	int FindNeighbordsAndBuriedAtoms(Atom& atom);
	int FindNeighborsForAtom(Atom& atom1);

	int GenerateToroidalSurface(Atom& atom1, Atom& atom2, Vec3 const uij, Vec3 const tij, ScValue rij, int between);
	int GenerateConvexSurface(Atom const & atom1);
	int GenerateConcaveSurface();

	// Function names similar to original source
	int SecondLoop(Atom &pAtom1);
	int ThirdLoop(Atom &pAtom1, Atom &pAtom, Vec3 const &uij, Vec3 const &tij, ScValue const rij);
	int CheckAtomCollision2(Vec3 const &pijk, Atom const &atom1, Atom const &atom2, std::vector<Atom*> const &atoms);
	int CheckPointCollision(Vec3 const &pcen, std::vector<Atom*> const &atoms);
	int CheckProbeCollision(Vec3 const &point, std::vector<const PROBE*> const nears, ScValue const r2);

	// Dot trimming
	ScValue TrimPeripheralBand(std::vector<DOT> const &sdots, std::vector<const DOT*> &trimmed_dots);
	int TrimPeripheralBandCheckDot(DOT const &dot, std::vector<DOT> const &sdots);

	// Elementary functions
	ScValue DistancePointToLine(Vec3 const &cen, Vec3 const &axis, Vec3 const &pnt);
	ScValue SubArc(Vec3 const &cen, ScValue const rad, Vec3 const &axis, ScValue const density,	Vec3 const &x, Vec3 const &v, std::vector<Vec3> &points);
	ScValue SubDiv(Vec3 const &cen, ScValue const rad, Vec3 const &x, Vec3 const &y, ScValue angle, ScValue density, std::vector<Vec3> &points);
	ScValue SubCir(Vec3 const &cen, ScValue const rad, Vec3 const &north, ScValue const density, std::vector<Vec3> &points);

	// Kernel functions for parallalization
	int CalcNeighborDistance(int const molecule, std::vector<const DOT*> const &my_dots, std::vector<const DOT*> const &their_dots);
	DOT const *CalcNeighborDistanceFindClosestNeighbor(DOT const &dot1, std::vector<const DOT*> const &their_dots);

  // Sort callback
	static int _atom_distance_cb(void *a1, void *a2);

#ifdef USECUDA
public:
	void GPUInit();

private:
	void cudaThrowException(int err, char const *fn, int line);
	ScValue CudaTrimPeripheralBand(std::vector<DOT> const &dots, std::vector<DOT const *> &trimmed_dots);
	int CudaFindClosestNeighbors(std::vector<const DOT*> const &my_dots, std::vector<const DOT*> const &their_dots, std::vector<DOT const*> &neighbors);
	core::Real GetTimerMs(clock_t &start);

#endif

};

} //namespace sc
} //namespace filters
} //namespace protocols

#endif

