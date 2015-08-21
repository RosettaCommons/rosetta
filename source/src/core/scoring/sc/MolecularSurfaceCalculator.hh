// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/sc/ShapeComplementarityCalculator.hh
/// @brief  Headers for the Molecular Surface Calculator
/// @author Luki Goldschmidt (luki@mbi.ucla.edu), refactored by Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_scoring_sc_MolecularSurfaceCalculator_hh
#define INCLUDED_core_scoring_sc_MolecularSurfaceCalculator_hh

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
#include <core/scoring/sc/MolecularSurfaceCalculator.fwd.hh>
#include <numeric/xyzVector.hh>

/// utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <vector>
#include <string>
#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

// OpenCL headers
#include <basic/gpu/Timer.hh>
#ifdef USEOPENCL
#include <basic/gpu/GPU.hh>
#endif

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////
// core::scoring::sc namespace
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// Types
//
//TODO Need to refactor to remove sc-specific results from this typedef
// Calculation results
typedef struct _RESULTS {
	core::Real sc;                 // Interface shape complementarity statistic
	core::Real area;               // Area of interface (A^2)
	core::Real distance;           // Interface separation (A)
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
	struct {   // Surface dot counts by type (from Connolly algorithm)
		core::Size convex;           // Number of convex surface dots
		core::Size concave;          // Number of concace surface dots
		core::Size toroidal;         // NUmber of toroidal surfac dots
	} dots;    // True if computed results are valid
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
	Atom const * atom;
} DOT;

// Molecular probe
typedef struct _PROBE {
	Atom const * pAtoms[3];
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
public:
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
// Molecular Surface Calculator class definition
////////////////////////////////////////////////////////////

class MolecularSurfaceCalculator : public utility::pointer::ReferenceCount {

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

#ifdef USEOPENCL
		core::SSize gpu;
		core::Size gpu_threads;
#endif
	} settings;

	MolecularSurfaceCalculator();
	virtual ~MolecularSurfaceCalculator();
	virtual int Init();
	virtual void Reset();
#if defined(WIN32) && !defined(WIN_PYROSETTA)
	int AddAtomWIN32(int molecule, Atom &atom);
#else
	int AddAtom(int molecule, Atom &atom);
#endif
	core::Size AddResidue(int molecule, core::conformation::Residue const &residue);

	/// @brief Generate molecular surfaces for the given pose.
	///// @details
	// This function initializes the calculator, adds all residues in the given pose, and generates molecular surfaces.
	//
	// The pose is partitioned into separate molecules across the given jump. If the given jump is 0, the entire pose is
	// loaded as molecule 1.
	virtual int Calc(core::pose::Pose const & pose, core::Size jump_id = 0);

	/// @brief Generate molecular surfaces for loaded atoms.
	///// @details
	// This function generates molecular surfaces for atoms added via AddAtom and AddResidue.
	//
	// Init() must be called before this function.
	virtual int Calc();

	std::vector<Atom> const & GetAtoms() { return run_.atoms; }
	std::vector<DOT> const & GetDots(int const moleculeid) { return run_.dots[moleculeid]; }
	RESULTS const & GetResults() { return run_.results; }

protected:
	/// @brief Generate untrimmed surfaces for the defined molecules.
	///// @details
	/// This function should be called within a try/catch block for ShapeComplementarityCalculatorException.
	/// Raises exception on error.
	void GenerateMolecularSurfaces();

	// This is a constant list of atom radii; declared static to avoid re-loading
	// between computations
	static std::vector<ATOM_RADIUS> radii_;

	int AssignAtomRadius(Atom &atom);
	int WildcardMatch(char const *query, char const *pattern, int const l);
	int ReadScRadii();
	void AddDot(int const molecule, int const type, Vec3 const coor, ScValue const area, Vec3 const pcen, Atom const &atom);

	struct {
		ScValue radmax;
		RESULTS results;
		std::vector<Atom> atoms;
		std::vector<DOT> dots[2];
		std::vector<PROBE> probes;
		Vec3 prevp;
		int prevburied;

	} run_;

	// Surface generation configuration
	virtual int AssignAttentionNumbers(std::vector<Atom>& atom);

#ifdef MULTI_THREADED
private:
	utility::thread::ReadWriteMutex sc_radii_mutex_;
#endif

private:

	// Molecular surface generation
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


	// Elementary functions
	ScValue DistancePointToLine(Vec3 const &cen, Vec3 const &axis, Vec3 const &pnt);
	ScValue SubArc(Vec3 const &cen, ScValue const rad, Vec3 const &axis, ScValue const density, Vec3 const &x, Vec3 const &v, std::vector<Vec3> &points);
	ScValue SubDiv(Vec3 const &cen, ScValue const rad, Vec3 const &x, Vec3 const &y, ScValue angle, ScValue density, std::vector<Vec3> &points);
	ScValue SubCir(Vec3 const &cen, ScValue const rad, Vec3 const &north, ScValue const density, std::vector<Vec3> &points);

	// Sort callback
	// static int _atom_distance_cb(void *a1, void *a2);
};

} //namespace sc
} //namespace filters
} //namespace protocols

#endif
