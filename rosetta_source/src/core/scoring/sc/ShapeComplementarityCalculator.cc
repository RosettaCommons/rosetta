// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     core/scoring/sc/ShapeComplementarityCalculator.fwd.hh
/// @brief    Headers for the Shape Complementarity Calculator
/// @detailed Lawrence & Coleman shape complementarity calculator (based on CCP4's sc)
/// @author   Luki Goldschmidt <luki@mbi.ucla.edu>

/// This code was ported from the original Fortran code found in CCP4:
/// Sc (Version 2.0): A program for determining Shape Complementarity
/// Copyright Michael Lawrence, Biomolecular Research Institute
/// 343 Royal Parade Parkville Victoria Australia
///
/// This version contains support for GPU-acceleration by CUDA-capable devices,
/// which provides a 10-25x speed up over the CPU-only code using a regular desktop
/// video card with 4 processors (32 cores). Define USECUDA and compile with
/// NVIDIA's nvcc. Also see ShapeComplementarityCalculator.hh.

#ifndef INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_cc
#define INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_cc

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/Residue.functions.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.fwd.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <utility/exit.hh>

// C++ headers
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
// AUTO-REMOVED #include <math.h>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>



// Cuda headers
#ifdef USECUDA
#include <cuda.h>
#endif


static basic::Tracer tr("core.scoring.sc.ShapeComplementarityCalculator");

using namespace core;
using namespace core::scoring::sc;

namespace core {
namespace scoring {
namespace sc {

std::vector<ATOM_RADIUS> ShapeComplementarityCalculator::radii_;  // static const

////////////////////////////////////////////////////////////////////////////
// Public class functions
////////////////////////////////////////////////////////////////////////////

/// @begin ShapeComplementarityCalculator::ShapeComplementarityCalculator()
/// @brief
/// ShapeComplementarityCalculator constructor, initializes default settings

ShapeComplementarityCalculator::ShapeComplementarityCalculator()
{
	memset(&settings, 0, sizeof(settings));

	// Defaults from sc source
	settings.rp = 1.7;
	settings.density = 15.0;
	settings.band = 1.5;
	settings.sep = 8.0;
	settings.weight = 0.5;
	settings.binwidth_dist = 0.02;
	settings.binwidth_norm = 0.02;

#ifdef USECUDA
	settings.gpu = -1;
#endif
}

/// @begin ShapeComplementarityCalculator::Init()
/// @brief
/// Initializes calculation and GPU (if used)
/// Init() is also called implicitly by the static CalcSc() function.

int ShapeComplementarityCalculator::Init()
{
	if(radii_.empty()) {
		Reset();
		ReadScRadii();
	}
	if(radii_.empty())
		return 0;

	return 1;
}

ShapeComplementarityCalculator::~ShapeComplementarityCalculator()
{
}

/// @begin ShapeComplementarityCalculator::Reset()
/// @brief
/// Reset calculator for another calculation.
/// Must be used when the ShapeComplementarityCalculator instance is re-used.
/// @detailed
/// Atom, probe and surface dot vectors are reset here. We don't clear them
/// after the calculation is finished in case the caller would like to use those
/// data elsewhere.

void ShapeComplementarityCalculator::Reset()
{
	// Free data
	run_.atoms.clear();
	run_.probes.clear();
	run_.dots[0].clear();
	run_.dots[1].clear();
	// Clear structures
	memset(&run_, 0, sizeof(run_));
}

/// @begin ShapeComplementarityCalculator::CalcSc()
/// @brief
/// Run the SC calculation on Pose and return just the sc statistic or -1 on error
/// @detailed
/// This is a static function and can be called without instantiating ShapeComplementarityCalculator.
/// The jump_id is used to partition the pose into two molecular surfaces; the first jump (1)
/// is used is no jump_id is explicity specified. Those desiring more control as to what residues
/// make up either surface should use the AddResidue() or even AddAtom() function instead.
/// Setting quick to true will perform a much faster calculation (~5-10 times faster) at the expense
/// of accuracy (about 0.05 units).
///
/// Example:
/// core::Real sc = core::scoring::sc::ShapeComplementarityCalculator( pose );

core::Real ShapeComplementarityCalculator::CalcSc(core::pose::Pose const & pose, core::Size jump_id, int quick)
{
	ShapeComplementarityCalculator sc;

	if(quick)
		sc.settings.density = 5;

	if(sc.Calc(pose, jump_id))
		return sc.GetResults().sc;
	else
		return -1;
}

/// @begin ShapeComplementarityCalculator::Calc()
/// @brief Run the SC calculation on a Pose, partitionied by jump_id
/// @detailed
/// This non-static function requires an instance of the ShapeComplementarityCalculator class.
/// The jump_id is used to partition the pose into two molecular surfaces. To control what
/// residues make up either surface, use the AddResidue() or even AddAtom() function instead.
/// Returns true on success. Results are retrieved with GetResults().
///
/// Example:
/// core::scoring::sc::ShapeComplementarityCalculator calc;
/// core::Real sc;
/// if(calc.Calc( pose ))
///   sc = calc.GetResults().sc;

int ShapeComplementarityCalculator::Calc(core::pose::Pose const & pose, core::Size jump_id)
{
	if(!Init())
		return 0;

	if( jump_id > pose.num_jump() ) {
		tr << "Jump ID out of bounds (pose has " << pose.num_jump() << " jumps)" << std::endl;
		return 0;
	}

	// Partition pose by jump_id
	ObjexxFCL::FArray1D_bool is_upstream ( pose.total_residue(), false );
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );

	for(Size i = 1; i <= pose.n_residue(); ++i) {
		core::conformation::Residue const & residue = pose.residue(i);
		if(residue.type().name() == "VRT")
			continue;
		AddResidue(is_upstream[i] ? 0 : 1, residue);
	}

	return Calc();
}

/// @begin ShapeComplementarityCalculator::Calc
/// @brief Run the SC calculation for previously defined molecules (via AddResidue or AddAtom calls)
/// @detailed
/// This function should be called the residues / atoms making up the two molecular surfaces
/// have been explicitly defined.
/// Returns true on success.

int ShapeComplementarityCalculator::Calc()
{
#ifdef USECUDA
	GPUInit();
#endif

	try {

	run_.results.valid = 0;

	if(run_.atoms.empty())
		throw ShapeComplementarityCalculatorException("No atoms defined");
	if(!run_.results.surface[0].nAtoms)
		throw ShapeComplementarityCalculatorException("No atoms defined for molecule 1");
	if(!run_.results.surface[1].nAtoms)
		throw ShapeComplementarityCalculatorException("No atoms defined for molecule 2");

	// Determine assign the attention numbers for each atom
	AssignAttentionNumbers(run_.atoms);

	// Now compute the surface for the atoms in the interface and its neighbours
	VERBOSE("Generating molecular surface, " << settings.density << " dots/A^2" << std::endl);
	CalcDotsForAllAtoms(run_.atoms);

	if(!run_.dots[0].size() || !run_.dots[1].size())
		  throw ShapeComplementarityCalculatorException("No molecular dots generated!");

	VERBOSE("           Convex dots: " << run_.results.dots.convex << std::endl);
	VERBOSE("         Toroidal dots: " << run_.results.dots.toroidal << std::endl);
	VERBOSE("          Concave dots: " << run_.results.dots.concave << std::endl);
	VERBOSE("Total surface dots (1): " << run_.dots[0].size() << std::endl);
	VERBOSE("Total surface dots (2): " << run_.dots[1].size() << std::endl);
	VERBOSE("    Total surface dots: " << (run_.dots[0].size()+run_.dots[1].size()) << std::endl);

	// Cut away the periphery of each surface
	VERBOSE("Trimming peripheral band, " << settings.band << "A range" << std::endl);

	std::vector<DOT const *> trimmed_dots[2];
	for(int i = 0; i < 2; ++i) {
		run_.results.surface[i].trimmedArea = TrimPeripheralBand(run_.dots[i], trimmed_dots[i]);
		if(!trimmed_dots[i].size())
			throw ShapeComplementarityCalculatorException("No molecular dots for surface %d", i);
		run_.results.surface[i].nTrimmedDots = trimmed_dots[i].size();
		run_.results.surface[i].nAllDots = run_.dots[i].size();
	}

	// Compute distance arrays and histograms for each surface
	VERBOSE("Computing surface separation and vectors" << std::endl);

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
	//run_.results.surface[2].nBuriedDots = (run_.results.surface[0].nBuriedDots + run_.results.surface[1].nBuriedDots);
	//run_.results.surface[2].nAccessibleDots = (run_.results.surface[0].nAccessibleDots + run_.results.surface[1].nAccessibleDots);
	run_.results.surface[2].trimmedArea = (run_.results.surface[0].trimmedArea + run_.results.surface[1].trimmedArea);

	run_.results.sc = run_.results.surface[2].s_median;
	run_.results.distance = run_.results.surface[2].d_median;
	run_.results.area = run_.results.surface[2].trimmedArea;
	run_.results.valid = 1;

	return 1;

	} catch(ShapeComplementarityCalculatorException e) {
		tr << "Failed: " << e.error << std::endl;
	}

	return 0;
}

/// @begin ShapeComplementarityCalculator::ReadScRadii()
/// @brief Read atom radius definitions from file
/// @defailed
/// This function is implicitly called, but can be overloaded or
/// called explicitly for custom handling of the atom radii library.
/// Returns true on success

int ShapeComplementarityCalculator::ReadScRadii()
{
	char const *fn = "sc_radii.lib";
	ATOM_RADIUS radius;
	utility::io::izstream in;

	if(!basic::database::open(in, fn)) {
		tr << "Failed to read " << fn << std::endl;
		return 0;
	}

	radii_.clear();

	while( in.good() ) {
		memset(&radius, 0, sizeof(radius));
		in >> radius.residue >> radius.atom >> radius.radius;
		//VERBOSE("Atom Radius: " << radius.residue << ", " << radius.atom << ", " << radius.radius << std::endl);
		if(*radius.residue && *radius.atom && radius.radius > 0)
			radii_.push_back(radius);
	}

	VERBOSE("Atom radii read: " << radii_.size() << std::endl);

	return !radii_.empty();
}

/// @begin ShapeComplementarityCalculator::AddResidue()
/// @brief Add a rosetta residue to a specific molecule
/// @detailed
/// Call this function when explicitly defining which residues belong to
/// which the molecular surface. If partitioning by jump_id is sufficient
/// for your application, you may use the Calc() or CalcSc() functions
/// instead.
/// Returns number of atoms added for the specified residue.
///
/// Example:
/// core::scoring::sc::ShapeComplementarityCalculator calc;
/// core::Real sc;
/// calc.Init();
/// calc.Reset(); // Only needed when re-using the calculator
/// for(core::Size i = 1; i <= pose.n_residue(); i++)
///   calc.AddResidue((i < 100), pose.residue(i));
/// if(calc.Calc())
///   sc = calc.GetResults().sc;

core::Size ShapeComplementarityCalculator::AddResidue(
	int molecule,
	core::conformation::Residue const & residue)
{
	Atom scatom;
	int n =0;

	if(!Init())
		return 0;

	// Only use heavy atoms for SC calculation
	for(Size i = 1; i <= residue.nheavyatoms(); ++i) {
		// Skip virtual atoms
		if(residue.is_virtual(i))
			continue;
		numeric::xyzVector<Real> xyz = residue.xyz(i);
		scatom.x(xyz.x());
		scatom.y(xyz.y());
		scatom.z(xyz.z());
		scatom.nresidue = 0;
		strncpy(scatom.residue, residue.name3().c_str(), sizeof(scatom.residue)-1);
		strncpy(scatom.atom, residue.atom_name(i).c_str()+1, sizeof(scatom.atom)-1);
		if(AddAtom(molecule, scatom))
			++n;
	}

	return n;
}

/// @begin ShapeComplementarityCalculator::AddAtom()
/// @brief Add an atom to a molecule for computation.
/// @detailed
/// Add an core::scoring::sc::Atom to the molecule.
/// Normally this is called by AddResidue(). Explicit addition
/// of atoms via this function is rarely needed.
/// This function also looks-up the atom radius and density.
/// Returns true on success.

int ShapeComplementarityCalculator::AddAtom(
		int molecule,
		Atom &atom)
{
	if(AssignAtomRadius(atom)) {
		molecule = (molecule == 1);
		atom.density = settings.density;
		atom.molecule = molecule;
		atom.natom = ++run_.results.nAtoms;
		atom.access = 0;
		run_.atoms.push_back(atom);
		++run_.results.surface[molecule].nAtoms;
		/*
		printf("AddAtom[%d] %d: %s:%s (%10.4f, %10.4f, %10.4f) = %.4f\n", molecule, run_.results.surface[molecule].nAtoms, atom.residue, atom.atom, atom.x(), atom.y(), atom.z(), atom.radius);
		*/
		return 1;

	} else {
		tr << "Failed to assign atom radius for residue "
			<< atom.residue << ":" << atom.atom
			<< ". Skipping atom!" << std::endl;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
// Protected class functions
////////////////////////////////////////////////////////////////////////////

// Look up the atom radius for an atom
int ShapeComplementarityCalculator::AssignAtomRadius(Atom &atom)
{
	std::vector<ATOM_RADIUS>::const_iterator radius;

	// Assign radius with wildcard matching
	for(radius = radii_.begin(); radius != radii_.end(); ++radius) {
		if(WildcardMatch(atom.residue, radius->residue, sizeof(atom.residue)) &&
			WildcardMatch(atom.atom, radius->atom, sizeof(atom.atom))) {
				atom.radius = radius->radius;
				return 1;
		}
	}

	return 0;
}

// Inline residue and atom name matching function
int ShapeComplementarityCalculator::WildcardMatch(
		char const *r,
		char const *pattern,
		int l)
{
	while(--l > 0) {
		if((*pattern != '*') && (*r != *pattern) && !(*r == ' ' && !*pattern))
			return 0;
		++r;
		++pattern;
	}
	return 1;
}

// Determine assign the attention numbers for each atom
int ShapeComplementarityCalculator::AssignAttentionNumbers(std::vector<Atom> &atoms)
{
	std::vector<Atom>::iterator pAtom1, pAtom2;

	for(pAtom1 = run_.atoms.begin(); pAtom1 < run_.atoms.end(); ++pAtom1) {
		// find nearest neighbour in other molecule
		ScValue dist_min = 99999.0, r;
		for(pAtom2 = run_.atoms.begin(); pAtom2 < run_.atoms.end(); ++pAtom2) {
			if(pAtom1->molecule == pAtom2->molecule)
				continue;
			r = pAtom1->distance(*pAtom2);
			if(r < dist_min)
				dist_min = r;
		}

		// check if within separator distance
		if(dist_min >= settings.sep) {
			// too far away from other molecule, blocker atom only
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

////////////////////////////////////////////////////////////////////////////
// Molecular surface calculation
////////////////////////////////////////////////////////////////////////////
// M. L. Connolly J. Appl. Crystallogr., 16, p548 - p558 (1983)
// Compute the surface for the atoms in the interface and its neighbours
////////////////////////////////////////////////////////////////////////////

int ShapeComplementarityCalculator::CalcDotsForAllAtoms(std::vector<Atom> &atoms)
{
	// Calc maximum radius for atoms list
	run_.radmax = 0.0;
	for(std::vector<Atom>::const_iterator pAtom1 = run_.atoms.begin(); pAtom1 < run_.atoms.end(); ++pAtom1) {
		if(pAtom1->radius > run_.radmax)
			run_.radmax = pAtom1->radius;
	}

	// Add dots for each atom in the list
	for(std::vector<Atom>::iterator pAtom1 = run_.atoms.begin(); pAtom1 < run_.atoms.end(); ++pAtom1) {
		Atom &atom1 = *pAtom1;
		if(atom1.atten <= 0)
			continue;

		// Find neighbor
		if(!FindNeighbordsAndBuriedAtoms(atom1))
			continue;

		if(!atom1.access)
			continue;
		if(atom1.atten <= ATTEN_BLOCKER)
			continue;
		if(atom1.atten == ATTEN_6 && atom1.buried.empty())
			continue;

		// Generate convex surface
		GenerateConvexSurface(atom1);
	}

	// Concave surface generation
	if(settings.rp > 0)
		GenerateConcaveSurface();

	return 1;
}

// Private atom distance sort callback used below
Atom *_atom_distance_ref = NULL;
int ShapeComplementarityCalculator::_atom_distance_cb(void *a1, void *a2)
{
	ScValue d1 = _atom_distance_ref->distance(*((Atom*)a1));
	ScValue d2 = _atom_distance_ref->distance(*((Atom*)a2));
	return d1 < d2 ? -1: (d2 > d1 ? 1 : 0);
}

// Calculate surface dots around a single atom (main loop in original code)
int ShapeComplementarityCalculator::FindNeighbordsAndBuriedAtoms(Atom &atom1)
{
	if(!FindNeighborsForAtom(atom1))
		return 0;

	// sort neighbors by distance from atom1
	_atom_distance_ref = &atom1;
	std::sort(atom1.neighbors.begin(), atom1.neighbors.end(), ShapeComplementarityCalculator::_atom_distance_cb);
	_atom_distance_ref = NULL;

	SecondLoop(atom1);

	return atom1.neighbors.size();
}

// Make a list of neighboring atoms from atom1
// (loop to label 100 in original code)
int ShapeComplementarityCalculator::FindNeighborsForAtom(Atom &atom1)
{
	std::vector<Atom>::iterator iAtom2;
	std::vector<Atom*> &neighbors = atom1.neighbors;
	ScValue d2;
	ScValue bridge;
	ScValue bb2 = pow(4 * run_.radmax + 4 * settings.rp, 2);
	int nbb = 0;

	for(iAtom2 = run_.atoms.begin(); iAtom2 < run_.atoms.end(); ++iAtom2) {
		Atom &atom2 = *iAtom2;
		if(atom1 == atom2 || atom2.atten <= 0)
			continue;

		if(atom1.molecule == atom2.molecule) {

			d2 = atom1.distance_squared(atom2);

			if(d2 <= 0.0001)
				throw ShapeComplementarityCalculatorException("Coincident atoms: %d:%s:%s @ (%.3f, %.3f, %.3f) == %d:%s:%s @ (%.3f, %.3f, %.3f)",
					atom1.natom, atom1.residue, atom1.atom, atom1.x(), atom1.y(), atom1.z(),
					atom2.natom, atom2.residue, atom2.atom, atom2.x(), atom2.y(), atom2.z());

			bridge = atom1.radius + atom2.radius + 2 * settings.rp;
			if (d2 >= bridge * bridge)
				continue;

			neighbors.push_back(&atom2);

		} else {

			if(atom2.atten < ATTEN_BURIED_FLAGGED)
				continue;

			d2 = atom1.distance_squared(atom2);
			if (d2 < bb2)
				++nbb;

			bridge = atom1.radius + atom2.radius + 2 * settings.rp;
			if (d2 >= bridge * bridge)
				continue;

			atom1.buried.push_back(&atom2);
		}
	}

	if(atom1.atten == ATTEN_6 && !nbb)
		return 0;

	if(neighbors.empty()) {
		// no neighbors
		atom1.access = 1;
		// no convex surface generation
		return 0;
	}

	return neighbors.size();
}

// second loop (per original code)
int ShapeComplementarityCalculator::SecondLoop(Atom &atom1)
{
	Vec3 uij, tij;
	ScValue erj, eri, rij, density, dij, asymm, far, contain;
	int between;
	std::vector<Atom*> &neighbors = atom1.neighbors;

	eri = atom1.radius + settings.rp;

	for(std::vector<Atom*>::iterator iAtom2 = neighbors.begin(); iAtom2 < neighbors.end(); ++iAtom2) {
		Atom &atom2 = **iAtom2;

		if(atom2 <= atom1)
			continue;

		erj = atom2.radius + settings.rp;
		density = (atom1.density + atom2.density) / 2;
		dij = atom1.distance(atom2);

		uij = (atom2 - atom1) / dij;
		asymm = (eri * eri - erj * erj) / dij;
		between = (ABS(asymm) < dij);

		tij = ((atom1 + atom2) * 0.5f) + (uij * (asymm * 0.5f));

		far = (eri + erj) * (eri + erj) - dij * dij;
		if (far <= 0.0)
			continue;

		far = sqrt(far);

		contain = dij * dij - ((atom1.radius - atom2.radius) * (atom1.radius - atom2.radius));
		if (contain <= 0.0)
			continue;

		contain = sqrt(contain);
		rij = 0.5 * far * contain / dij;

		if (neighbors.size() <= 1) {
			atom1.access = 1;
			atom2.access = 1;
			break;
		}

		ThirdLoop(atom1, atom2, uij, tij, rij);

		if(atom1.atten > ATTEN_BLOCKER || (atom2.atten > ATTEN_BLOCKER && settings.rp > 0.0))
			GenerateToroidalSurface(atom1, atom2, uij, tij, rij, between);
	}

	return 1;
}

// third loop (per original code)
int ShapeComplementarityCalculator::ThirdLoop(
		Atom& atom1,
		Atom& atom2,
		Vec3 const &uij,
		Vec3 const &tij,
		ScValue rij)
{
	std::vector<Atom*> &neighbors = atom1.neighbors;
	ScValue eri, erj, erk, djk, dik;
	ScValue asymm, dt, dtijk2, hijk, isign, wijk, swijk, rkp2;
	int is0;
	Vec3 uik, dijk, pijk, utb, iujk, tv, tik, bijk, uijk;

	eri = atom1.radius + settings.rp;
	erj = atom2.radius + settings.rp;

	for(std::vector<Atom*>::iterator iAtom3 = neighbors.begin();
			iAtom3 < neighbors.end(); ++iAtom3) {
		Atom &atom3 = **iAtom3;
		if(atom3 <= atom2)
			continue;

		erk = atom3.radius + settings.rp;
		djk = atom2.distance(atom3);
		if(djk >= erj+erk)
			continue;

		dik = atom1.distance(atom3);
		if(dik >= eri+erk)
			continue;

		if(atom1.atten <= ATTEN_BLOCKER && atom2.atten <= ATTEN_BLOCKER && atom3.atten <= ATTEN_BLOCKER)
			continue;

		uik = (atom3 - atom1) / dik;
		dt = uij.dot(uik);
		wijk = acosf(dt);
		swijk = sinf(wijk);

		if(dt >= 1.0 || dt <= -1.0 || wijk <= 0.0 || swijk <= 0.0) {
			// collinear and other
			dtijk2 = tij.distance(atom3);
			rkp2 = erk * erk - rij * rij;
			if(dtijk2 < rkp2)
				// 600
				return 0;
			continue;
		}

		uijk = uij.cross(uik) / swijk;
		utb = uijk.cross(uij);
		asymm = (eri*eri - erk*erk) / dik;
		tik = (atom1 + atom3)*0.5f + uik*asymm*0.5f;

		// Was: tv = uik * (tik - tij);
		tv = (tik - tij);
		tv.x(uik.x() * tv.x());
		tv.y(uik.y() * tv.y());
		tv.z(uik.z() * tv.z());

		dt = tv.x() + tv.y() + tv.z();
		bijk = tij + utb * dt / swijk;
		hijk = eri*eri - bijk.distance_squared(atom1);
		if(hijk <= 0.0)
			// no height, skip
			continue;

		hijk = sqrt(hijk);
		for(is0 = 1; is0 <= 2; ++is0) {
			isign = 3 - 2 * is0;
			pijk = bijk + uijk * hijk * isign;

			// check for collision
			if(CheckAtomCollision2(pijk, atom2, atom3, neighbors))
				continue;

			// new probe position
			PROBE probe;
			if(isign > 0) {
				probe.pAtoms[0] = &atom1;
				probe.pAtoms[1] = &atom2;
				probe.pAtoms[2] = &atom3;
			} else {
				probe.pAtoms[0] = &atom2;
				probe.pAtoms[1] = &atom1;
				probe.pAtoms[2] = &atom3;
			}
			probe.height = hijk;
			probe.point = pijk;
			probe.alt = uijk * isign;
			run_.probes.push_back(probe);

			atom1.access = 1;
			atom2.access = 1;
			atom3.access = 1;
		}
	}

	return 1;
}

// Check two atoms against a list of neighbors for collision
int ShapeComplementarityCalculator::CheckAtomCollision2(
		Vec3 const &pijk,
		Atom const &atom1,
		Atom const &atom2,
		std::vector<Atom*> const &atoms)
{
	for(std::vector<Atom*>::const_iterator ineighbor = atoms.begin();
			ineighbor < atoms.end(); ++ineighbor) {
		Atom const &neighbor = **ineighbor;
		if(&atom1 == &neighbor || &atom2 == &neighbor)
			continue;
		if(pijk.distance_squared(neighbor) <= pow(neighbor.radius + settings.rp, 2))
			// collision detected
			return 1;
	}
	return 0;
}

// Generate convex surface for a specific atom
int ShapeComplementarityCalculator::GenerateConvexSurface(Atom const &atom1)
{
	std::vector<Atom*> const &neighbors = atom1.neighbors;
	Atom const *neighbor;
	Vec3 north(0, 0, 1);
	Vec3 south(0, 0, -1);
	Vec3 eqvec(1, 0, 0);
	Vec3 vtemp, vql, uij, tij, pij, cen, pcen;
	ScValue dt, erj, dij, eri, far, contain;
	ScValue area, asymm, cs, ps, rad, ri, rij, rj;

	ri = atom1.radius;
	eri = (atom1.radius + settings.rp);

	if(!neighbors.empty()) {
		// use first neighbor
		neighbor = neighbors[0];

		north = atom1 - *neighbor;
		north.normalize();

		vtemp.x(north.y()*north.y() + north.z()*north.z());
		vtemp.y(north.x()*north.x() + north.z()*north.z());
		vtemp.z(north.x()*north.x() + north.y()*north.y());
		vtemp.normalize();

		dt = vtemp.dot(north);
		if(ABS(dt) > 0.99)
			vtemp = Vec3(1, 0, 0);

		eqvec = north.cross(vtemp);
		eqvec.normalize();
		vql = eqvec.cross(north);

		rj = neighbor->radius;
		erj = neighbor->radius + settings.rp;
		dij = atom1.distance(*neighbor);
		uij = (*neighbor - atom1) / dij;

		asymm = (eri*eri - erj*erj) / dij;
		tij = ((atom1 + *neighbor) * 0.5f) + (uij * (asymm * 0.5f));
		far = pow(eri + erj, 2) - dij*dij;
		if(far <= 0.0)
			throw ShapeComplementarityCalculatorException("Imaginary far for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		far = sqrt(far);

		contain = pow(dij, 2) - pow(ri - rj, 2);
		if(contain <= 0.0)
			throw ShapeComplementarityCalculatorException("Imaginary contain for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		contain = sqrt(contain);
		rij = 0.5 * far * contain / dij;
		pij = tij + (vql * rij);
		south = (pij - atom1) / eri;

		if(north.cross(south).dot(eqvec) <= 0.0)
			throw ShapeComplementarityCalculatorException("Non-positive frame for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
	}

	// Generate subdivided arc
	std::vector<Vec3> lats;
	Vec3 o(0);
	cs = SubArc(o, ri, eqvec, atom1.density, north, south, lats);

	if(lats.empty())
		return 0;

	// Project onto north vector
	std::vector<Vec3> points;
	for(std::vector<Vec3>::iterator ilat = lats.begin(); ilat < lats.end(); ++ilat) {
		dt = ilat->dot(north);
		cen = atom1 + (north*dt);
		rad = ri*ri - dt*dt;
		if(rad <= 0.0)
			continue;
		rad = sqrt(rad);

		points.clear();
		ps = SubCir(cen, rad, north, atom1.density, points);
		if(points.empty())
			continue;
		area = ps * cs;

		for(std::vector<Vec3>::iterator point = points.begin();
				point < points.end(); ++point) {

			pcen = atom1 + ((*point - atom1) * (eri/ri));

			// Check for collision
			if(CheckPointCollision(pcen, neighbors))
				continue;

			// No collision, put point
			++run_.results.dots.convex;
			AddDot(atom1.molecule, 1, *point, area, pcen, atom1);
		}
	}
	return 1;
}

// Check a point for collision against a list of atoms
int ShapeComplementarityCalculator::CheckPointCollision(
		Vec3 const &pcen,
		std::vector<Atom*> const &atoms)
{
	for(std::vector<Atom*>::const_iterator ineighbor = atoms.begin()+1;
		ineighbor < atoms.end(); ++ineighbor) {
		if(pcen.distance(**ineighbor) <= ((*ineighbor)->radius + settings.rp))
			// collision detected
			return 1;
	}
	return 0;
}

// Generate toroidal surface between two atoms
int ShapeComplementarityCalculator::GenerateToroidalSurface(
		Atom &atom1,
		Atom &atom2,
		Vec3 const uij,
		Vec3 const tij,
		ScValue rij,
		int between)
{
	std::vector<Atom*> &neighbors = atom1.neighbors;
	ScValue density, ri, rj, rb, rci, rcj, rs, e, edens, eri, erj, erl, dtq, pcusp, anglei, anglej, dt, ts, ps, area;
	Vec3 pi, pj, axis, dij, pqi, pqj, qij, qjk, qj;

	std::vector<Vec3> subs;

	// following Fortran original
	// will be optimized by compiler
	density = (atom1.density + atom2.density) / 2;
	ri = atom1.radius;
	rj = atom2.radius;
	eri = (atom1.radius + settings.rp);
	erj = (atom2.radius + settings.rp);
	rci = rij * atom1.radius / eri;
	rcj = rij * atom2.radius / erj;
	rb = rij - settings.rp;

	if(rb <= 0.0)
		rb = 0.0;

	rs = (rci + 2 * rb + rcj) / 4;
	e = rs / rij;
	edens = e * e * density;

	ts = SubCir(tij, rij, uij, edens, subs);
	if(subs.empty())
		return 0;

	for(std::vector<Vec3>::iterator isub = subs.begin(); isub < subs.end(); ++isub) {
		Vec3 &sub = *isub;

		// check for collision
		int tooclose = 0;
		ScValue d2 = 0;
		for(std::vector<Atom*>::iterator ineighbor = neighbors.begin();
				!tooclose && ineighbor < neighbors.end() && !tooclose;
				++ineighbor) {
			Atom const &neighbor = **ineighbor;	// for readability
			if(atom2 == neighbor)
				continue;
			erl = neighbor.radius + settings.rp;
			d2 = sub.distance_squared(neighbor);
			tooclose = d2 < (erl * erl);
		}
		if(tooclose)
			continue;

		// no collision, toroidal arc generation
		Vec3 &pij = sub;
		atom1.access = 1;
		atom2.access = 1;

		if(atom1.atten == ATTEN_6 && atom2.atten == ATTEN_6&& atom1.buried.empty())
			continue;

		pi = (atom1 - pij) / eri;
		pj = (atom2 - pij) / erj;
		axis = pi.cross(pj);
		axis.normalize();

		dtq = pow(settings.rp, 2) - pow(rij, 2);
		pcusp = dtq > 0 && between;
		if(pcusp) {
			// point cusp -- two shortened arcs
			dtq = sqrt(dtq);
			qij = tij - uij * dtq;
			qjk = tij + uij * dtq;
			pqi = (qij - pij) / (ScValue)settings.rp;
			pqj = Vec3(0.0);

		} else {
			// no cusp
			pqi = pi + pj;
			pqi.normalize();
			pqj = pqi;
		}

		dt = pqi.dot(pi);
		if(dt >= 1.0 || dt <= -1.0)
			return 0;
		anglei = acosf(dt);

		dt = pqj.dot(pj);
		if(dt >= 1.0 || dt <= -1.0)
			return 0;
		anglej = acosf(dt);

		// convert two arcs to points
		if(atom1.atten >= ATTEN_2) {
			std::vector<Vec3> points;
			ps = SubArc(pij, settings.rp, axis, density, pi, pqi, points);
			for(std::vector<Vec3>::iterator point = points.begin(); point < points.end(); ++point) {
				area = ps * ts * DistancePointToLine(tij, uij, *point) / rij;
				++run_.results.dots.toroidal;
				AddDot(atom1.molecule, 2, *point, area, pij, atom1);
			}
		}

		if(atom2.atten >= ATTEN_2) {
			std::vector<Vec3> points;
			ps = SubArc(pij, settings.rp, axis, density, pqj, pj, points);
			for(std::vector<Vec3>::iterator point = points.begin(); point < points.end(); ++point) {
				area = ps * ts * DistancePointToLine(tij, uij, *point) / rij;
				++run_.results.dots.toroidal;
				AddDot(atom1.molecule, 2, *point, area, pij, atom2);
			}
		}
	}
	return 1;
}

// Generate concave surface for all probes
int ShapeComplementarityCalculator::GenerateConcaveSurface()
{
	std::vector<PROBE const *> lowprobs, nears;

	// collect low probes
	for(std::vector<PROBE>::iterator probe = run_.probes.begin();
			probe < run_.probes.end(); ++probe) {
		if(probe->height < settings.rp)
			lowprobs.push_back(&(*probe));
	}

	for(std::vector<PROBE>::iterator probe = run_.probes.begin();
			probe < run_.probes.end(); ++probe) {

		if(	probe->pAtoms[0]->atten == ATTEN_6 &&
			probe->pAtoms[1]->atten == ATTEN_6 &&
			probe->pAtoms[2]->atten == ATTEN_6) {
			continue;
		}

		Vec3 &pijk = probe->point, &uijk = probe->alt;
		ScValue hijk = probe->height;
		ScValue density = (
				probe->pAtoms[0]->density +
				probe->pAtoms[1]->density +
				probe->pAtoms[2]->density ) / 3;

		// gather nearby low probes
		nears.clear();
		for(std::vector<PROBE const *>::const_iterator lprobe = lowprobs.begin();
				lprobe < lowprobs.end(); ++lprobe) {
			if(&(*probe) == *lprobe)
				continue;

			ScValue d2 = pijk.distance_squared((*lprobe)->point);
			if(d2 > 4 * pow(settings.rp, 2))
				continue;

			nears.push_back(*lprobe);
		}

		// set up vectors from probe center to three atoms
		Vec3 vp[3], vectors[3];
		for(int i = 0; i < 3; ++i) {
			vp[i] = *(probe->pAtoms[i]) - pijk;
			vp[i].normalize();
		}

		// set up vectors to three cutting planes
		vectors[0] = vp[0].cross(vp[1]);
		vectors[1] = vp[1].cross(vp[2]);
		vectors[2] = vp[2].cross(vp[0]);
		vectors[0].normalize();
		vectors[1].normalize();
		vectors[2].normalize();

		// find latitude of highest vertex of triangle
		ScValue dm = -1.0;
		int mm = 0;
		for(int i = 0; i < 3; ++i) {
			ScValue dt = uijk.dot(vp[i]);
			if(dt > dm) {
				dm = dt;
				mm = i;
			}
		}

		// create arc for selecting latitudes
		Vec3 south = -uijk;
		Vec3 axis = vp[mm].cross(south);
		axis.normalize();

		std::vector<Vec3> lats;
		Vec3 o(0);
		ScValue cs;

		cs = SubArc(o, settings.rp, axis, density, vp[mm], south, lats);
		if(lats.empty())
			continue;

		std::vector<Vec3> points;
		for(std::vector<Vec3>::iterator ilat = lats.begin();
				ilat < lats.end(); ++ilat) {
			ScValue dt, area, rad, ps;
			Vec3 cen;

			dt = ilat->dot(south);
			cen = south * dt;
			rad = pow(settings.rp, 2) - pow(dt, 2);
			if(rad <= 0.0)
				continue;
			rad = sqrtf(rad);

			points.clear();
			ps = SubCir(cen, rad, south, density, points);
			if(points.empty())
				continue;

			area = ps * cs;

			for(std::vector<Vec3>::iterator point = points.begin();
					point < points.end(); ++point) {
				// check against 3 planes
				int bail = 0;
				for(int i = 0; i < 3; ++i) {
					ScValue dt = point->dot(vectors[i]);
					if(dt >= 0.0) {
						bail = 1;
						break;
					}
				}
				if(bail)
					continue;

				*point += pijk;

				if((hijk < settings.rp && !nears.empty()) &&
					CheckProbeCollision(*point, nears, pow(settings.rp, 2)))
						continue;

				// determine which atom the surface point is closest to
				int mc = 0;
				ScValue dmin = 2 * settings.rp;
				for(int i = 0; i < 3; ++i) {
					ScValue d = point->distance(*(probe->pAtoms[i])) -
							probe->pAtoms[i]->radius;
					if(d < dmin) {
						dmin = d;
						mc = i;
					}
				}

				// No collision, put point
				++run_.results.dots.concave;
				AddDot(probe->pAtoms[mc]->molecule, 3, *point, area, pijk, *probe->pAtoms[mc]);
			}
		}
	}
	return 1;
}

// Check a point against a set of probes for collision within radius^2
int ShapeComplementarityCalculator::CheckProbeCollision(
		Vec3 const &point,
		std::vector<PROBE const *> const nears,
		ScValue const r2)
{
	for(std::vector<const PROBE*>::const_iterator near = nears.begin();
			near < nears.end(); ++near) {
		if(point.distance_squared((*near)->point) < r2)
			// Collision
			return 1;
	}
	return 0;
}

// Add a molecular dot
void ShapeComplementarityCalculator::AddDot(
		int const molecule,
		int const type,
		Vec3 const coor,
		ScValue const area,
		Vec3 const pcen,
		Atom const &atom)
{
	DOT dot = { coor, Vec3(), area, 0, type, &atom };
	ScValue pradius = settings.rp, erl;

	// calculate outward pointing unit normal vector
	if(pradius <= 0)
		dot.outnml = coor - atom;
	else
		dot.outnml = (pcen - coor) / pradius;

	// determine whether buried

	// first check whether probe changed
	if(pcen.distance_squared(run_.prevp) <= 0.0) {
		dot.buried = run_.prevburied;

	} else {

		// check for collision with neighbors in other molecules
		dot.buried = 0;
		for(std::vector<Atom*>::const_iterator iNeighbor = atom.buried.begin();
				iNeighbor < atom.buried.end();
				++iNeighbor) {
			erl = (*iNeighbor)->radius + pradius;
			ScValue d = pcen.distance_squared(**iNeighbor);
			if(d <= erl*erl) {
				dot.buried = 1;
				break;
			}

		}

		run_.prevp = pcen;
		run_.prevburied = dot.buried;
	}

	run_.dots[molecule].push_back(dot);
}

//
// Calculate distance from point to line
ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::DistancePointToLine(
		Vec3 const &cen,
		Vec3 const &axis,
		Vec3 const &pnt)
{
	Vec3 vec = pnt - cen;
	ScValue dt = vec.dot(axis);
	ScValue d2 = vec.magnitude_squared() - pow(dt, 2);
	return d2 < 0.0 ? 0.0 : sqrt(d2);
}

// Generate sub arc of molecular dots centered around a defined point
ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::SubArc(
		Vec3 const &cen,
		ScValue const rad,
		Vec3 const &axis,
		ScValue const density,
		Vec3 const &x,
		Vec3 const &v,
		std::vector<Vec3> &points)
{
	Vec3 y;
	ScValue angle;
	ScValue dt1, dt2;

	y = axis.cross(x);
	dt1 = v.dot(x);
	dt2 = v.dot(y);
	angle = atan2(dt2, dt1);

	if(angle < 0.0)
		angle = angle + 2*PI;

	return SubDiv(cen, rad, x, y, angle, density, points);
}

// Subdivide defined arc and generate molecular dots
ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::SubDiv(
		Vec3 const &cen,
		ScValue const rad,
		Vec3 const &x,
		Vec3 const &y,
		ScValue const angle,
		ScValue const density,
		std::vector<Vec3> &points)
{
	ScValue delta, a, c, s, ps;
	int i;

	delta = 1.0 / (sqrt(density) * rad);
	a = - delta / 2;

	for(i = 0; i < MAX_SUBDIV; ++i) {
		a = a + delta;
		if(a > angle)
			break;

		c = rad * cosf(a);
		s = rad * sinf(a);
		points.push_back(Vec3(cen + x*c + y*s));
	}

	if (a + delta < angle)
		throw ShapeComplementarityCalculatorException("Too many subdivisions");

	if (!points.empty())
		ps = rad * angle / points.size();
	else
		ps = 0.0;

	return ps;
}

// Generate an arbitrary unit vector perpendicular to axis
ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::SubCir(
		Vec3 const &cen,
		ScValue const rad,
		Vec3 const &axis,
		ScValue const density,
		std::vector<Vec3> &points)
{
	Vec3 v1, v2, x, y;
	ScValue dt;

	v1.x(pow(axis.y(), 2) + pow(axis.z(), 2));
	v1.y(pow(axis.x(), 2) + pow(axis.z(), 2));
	v1.z(pow(axis.x(), 2) + pow(axis.y(), 2));
	v1.normalize();
	dt = v1.dot(axis);

	if(ABS(dt) > 0.99) {
		v1.x(1.0);
		v1.y(0.0);
		v1.z(0.0);
	}

	v2 = axis.cross(v1);
	v2.normalize();
	x = axis.cross(v2);
	x.normalize();
	y = axis.cross(x);

	return SubDiv(cen, rad, x, y, 2.0 * PI, density, points);
}

// SC molecular dot trimming, vector dot product calculation and statistics
// Trim dots and retain only the peripheral band

ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::TrimPeripheralBand(
		std::vector<DOT> const &sdots,
		std::vector<DOT const *> &trimmed_dots)
{
	ScValue area = 0;

	if(sdots.empty())
		return 0.0;

#ifdef USECUDA
	if(settings.gpu) {
		area = CudaTrimPeripheralBand(sdots, trimmed_dots);
	} else {
#endif

	// Loop over one surface
	// If a point is buried then see if there is an accessible point within distance band

	for(std::vector<DOT>::const_iterator idot = sdots.begin(); idot < sdots.end(); ++idot) {
		DOT const &dot = *idot;
		// Paralelleizable kernel function
		if(dot.buried && TrimPeripheralBandCheckDot(dot, sdots)) {
			area += dot.area;
			trimmed_dots.push_back(&dot);
		}
	}

#ifdef USECUDA
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

	for(std::vector<DOT>::const_iterator idot2 = sdots.begin(); idot2 < sdots.end(); ++idot2) {
		DOT const &dot2 = *idot2;
		if(&dot == &dot2)
			continue;
		if(dot2.buried)
			continue;
		if(dot.coor.distance_squared(dot2.coor) <= r2)
			return 0;
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
	std::map<int,int> dbins;	// Distance bins
	std::map<int,int> sbins;	// Vector dot product bins (sc)
	ScValue norm_sum = 0.0, distmin_sum = 0.0;
	int ibin;
	ScValue total = 0.0;

	if(my_dots.empty() || their_dots.empty())
		return 0;

	for(std::vector<DOT const*>::const_iterator idot = my_dots.begin();
			idot < my_dots.end(); ++idot) {
		//if((*idot)->buried)
		//	run_.results.surface[molecule].nBuriedDots++;
		//else
		//	run_.results.surface[molecule].nAccessibleDots++;
		total += (*idot)->area;
	}

#ifdef USECUDA
	std::vector<DOT const*> neighbors;
	std::vector<DOT const*>::const_iterator iNeighbor;

	if(settings.gpu) {
		CudaFindClosestNeighbors(my_dots, their_dots, neighbors);
		iNeighbor = neighbors.begin();
	}
#endif

	for(std::vector<DOT const*>::const_iterator idot = my_dots.begin();
			idot < my_dots.end(); ++idot) {
		DOT const &dot1 = **idot;

		ScValue distmin, r;
		DOT const *neighbor = NULL;

#ifdef USECUDA
		if(settings.gpu)
			neighbor = *iNeighbor++;
		else
#endif
		neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);

		if(!neighbor)
			continue;

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
		// cpjx		  r = r * exp( - (distmin**2) / 2.)
		r = r * exp( - pow(distmin, 2) * settings.weight );
		// rounding errors a problem, so ensure abs(r) <1
		r = MIN(0.999, MAX(r, -0.999));
		norm_sum += r;

		// left_trunc ScValue to int ibin
		// otherwise: (int)-0.9 = 0.
		r /= settings.binwidth_norm;
		if(r >= 0)
			ibin = (int)r;
		else
			ibin = (int)r -1;
		++sbins[ibin];
	}

	// Determine the last distance bin that has anything in it
	// Accumulate percentages and area from all filled distance bins
	ScValue abin, cumarea =0, cumperc = 0, perc, c;
	ScValue rleft =0, rmedian =0;
	std::map<int,int>::const_iterator it;

	VERBOSE(std::endl);
	VERBOSE("Distance between surfaces D(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl);
	VERBOSE("From - To\tArea\tCum. Area\t%\tCum. %" << std::endl);

	for(it = dbins.begin(); it != dbins.end(); ++it) {
		abin = total * (it->second) / my_dots.size();
		cumarea += abin;
		perc = abin * 100 / total;
		c = cumperc + perc;
		if(cumperc <= 50 && c >= 50) {
			rleft = (it->first) * settings.binwidth_dist;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_dist / ( c - cumperc );
		}
		cumperc = c;

		#ifndef WIN_PYROSETTA
			if(settings.verbose) {
				char buf[128];
				snprintf(buf, sizeof(buf),
					"%.2f - %.2f\t%.1f\t%.1f\t%.1f\t%.1f",
					(ScValue)it->first * settings.binwidth_dist,
					(ScValue)it->first * settings.binwidth_dist + settings.binwidth_dist,
					abin, cumarea,
					perc, cumperc);
				tr << buf << std::endl;
			}
		#endif
	}

	run_.results.surface[molecule].d_mean = distmin_sum / my_dots.size();
	run_.results.surface[molecule].d_median = rmedian;

	VERBOSE(std::endl);
	VERBOSE("Surface complementarity S(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl);
	VERBOSE("From - To\tNumber\t%\tCumm. %" << std::endl);

	cumperc = 0;
	for(it = sbins.begin(); it != sbins.end(); ++it) {
		perc = (ScValue)(it->second) * 100 / my_dots.size();
		c = cumperc + perc;
		if(cumperc <= 50 && c >= 50) {
			rleft = (ScValue)(it->first) * settings.binwidth_norm;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_norm / ( c - cumperc );
		}
		cumperc = c;

		#ifndef WIN_PYROSETTA
			if(settings.verbose) {
				char buf[128];
				snprintf(buf, sizeof(buf),
					"%.2f - %.2f\t%d\t%.1f\t%.1f",
					(ScValue)-it->first * settings.binwidth_norm - settings.binwidth_norm,
					(ScValue)-it->first * settings.binwidth_norm,
					it->second, perc, cumperc);
				tr << buf << std::endl;
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
		std::vector<DOT const*> const &their_dots)
{
	ScValue distmin = 999999.0, d;
	DOT const *neighbor = NULL;

	// Loop over the entire surface: find and flag neighbour of each point
	// that we're interested in and store nearest neighbour pointer

	for(std::vector<DOT const*>::const_iterator idot2 = their_dots.begin();
			idot2 < their_dots.end(); ++idot2) {
		DOT const &dot2 = **idot2;
		if(!dot2.buried)
			continue;
		d = dot2.coor.distance_squared(dot1.coor);
		if(d <= distmin) {
			distmin = d;
			neighbor = &dot2;
		}
	}
	return neighbor;
}

#ifdef USECUDA

#define cudaAssert(err) cudaThrowException(err, __FILE__, __LINE__)
#define UPPER_MULTIPLE(n,d) (((n)%(d)) ? (((n)/(d)+1)*(d)) : (n))

void ShapeComplementarityCalculator::cudaThrowException(int err, char const *fn, int line)
{
	if (cudaSuccess != err)
		throw ShapeComplementarityCalculatorException("CUDA Exception at %s:%d: %s", fn, line, cudaGetErrorString((cudaError_t)err));
}

core::Real inline ShapeComplementarityCalculator::GetTimerMs(clock_t &start)
{
	clock_t now = clock();
	core::Real d = (now - start)/(CLOCKS_PER_SEC/1000);
	return d;
}

void ShapeComplementarityCalculator::GPUInit()
{
	if(!settings.gpu)
		return;

	try {
		cudaDeviceProp deviceProp;

		// TODO: Can we cache this between instances? private static?
		// This is likely a per-physical host setting, not per thread of even
		// process.

		if(settings.gpu < 0) {
			// Detect GPU and available threads
			int dev, deviceCount;
			if (cudaGetDeviceCount(&deviceCount) == cudaSuccess) {
				for (dev = 0; dev < deviceCount; ++dev) {
					cudaGetDeviceProperties(&deviceProp, dev);
					if(deviceProp.maxThreadsPerBlock*deviceProp.multiProcessorCount >
						settings.gpu_threads*settings.gpu_proc) {
						settings.gpu_proc = deviceProp.multiProcessorCount;
						settings.gpu_threads = settings.gpu_threads ?
							MIN(settings.gpu_threads, deviceProp.maxThreadsPerBlock) :
							deviceProp.maxThreadsPerBlock;
						settings.gpu = dev + 1;
					}
				}
			}
		} else if(settings.gpu > 0) {
			// Detect GPU and available threads
			cudaAssert( cudaGetDeviceProperties(&deviceProp, settings.gpu -1) );
			settings.gpu_proc = deviceProp.multiProcessorCount;
			settings.gpu_threads = settings.gpu_threads ?
				MIN(settings.gpu_threads, deviceProp.maxThreadsPerBlock) :
				deviceProp.maxThreadsPerBlock;
		}

		if(!settings.gpu_threads) {
			if(settings.gpu < 0)
				tr << "No GPU available. GPU acceleration disabled." << std::endl;
			settings.gpu = 0;
		} else {
			VERBOSE("GPU support enabled: " << deviceProp.name <<
					" [" << (deviceProp.clockRate/1000) << " MHz, capability " <<
					deviceProp.major << "." << deviceProp.minor << "] with " <<
					deviceProp.multiProcessorCount << " multi processors, " <<
					settings.gpu_threads <<" threads." <<
					std::endl);
		}

	} catch(ShapeComplementarityCalculatorException e) {
		tr << e.error << std::endl;
		settings.gpu = 0;
	}
}

ShapeComplementarityCalculator::ScValue ShapeComplementarityCalculator::CudaTrimPeripheralBand(
		std::vector<DOT> const &dots,
		std::vector<DOT const*> &trimmed_dots)
{
	int n, nBur, nAcc;
	int threads;
	ScValue area = 0;
	clock_t timer;

	threads = MIN(1024, settings.gpu_threads);
	n = dots.size();
	timer = clock();

	// Host and device (GPU) memory pointers for dot coordinates and results
	float3 *hAccDotCoords, *phAccDotCoords;
	float3 *hBurDotCoords, *phBurDotCoords;
	float3 *dAccDotCoords;
	float3 *dBurDotCoords;
	char *hDotColl;
	char *dDotColl;

	// Allocate host memory
	// Use cudaHostAlloc for DMA zero-copy
	cudaAssert( cudaHostAlloc((void**)&hAccDotCoords, n * sizeof(*hAccDotCoords), cudaHostAllocDefault) );
	cudaAssert( cudaHostAlloc((void**)&hBurDotCoords, n * sizeof(*hBurDotCoords), cudaHostAllocDefault) );

	// Make CUDA copy of (x, y, z) buried and accessible coordinates
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
			*phBurDotCoords++ = *((float3*)&idot->coor.x());
		else
			*phAccDotCoords++ = *((float3*)&idot->coor.x());
#endif
	}
	nBur = phBurDotCoords - hBurDotCoords;
	nAcc = phAccDotCoords - hAccDotCoords;

	// Allocate host memory for results (detected collisions)
	cudaAssert( cudaHostAlloc((void**)&hDotColl, nBur * sizeof(*hDotColl), cudaHostAllocDefault) );

	// Allocate GPU memory
	cudaAssert( cudaMalloc((void **)&dBurDotCoords, UPPER_MULTIPLE(nBur, threads) * sizeof(*dBurDotCoords)) );
	cudaAssert( cudaMalloc((void **)&dAccDotCoords, nAcc * sizeof(*dAccDotCoords)) );
	cudaAssert( cudaMalloc((void **)&dDotColl, UPPER_MULTIPLE(nBur, threads) * sizeof(*dDotColl)) );

	// Copy data from host to GPU
	cudaAssert( cudaMemcpy(dBurDotCoords, hBurDotCoords, nBur * sizeof(*dBurDotCoords), cudaMemcpyHostToDevice) );
	cudaAssert( cudaMemcpy(dAccDotCoords, hAccDotCoords, nAcc * sizeof(*dAccDotCoords), cudaMemcpyHostToDevice) );

	// Run kernel in multi-threaded blocks on GPU
	::_cuda_sccalc_TrimPeripheralBand(UPPER_MULTIPLE(nBur, threads)/threads, threads,
				dAccDotCoords, nAcc,
				dBurDotCoords,
				dDotColl, pow(settings.band, 2));

	// Wait for threads to finish and copy back memory from GPU to host
	cudaAssert( cudaThreadSynchronize() );
	cudaAssert( cudaMemcpy(hDotColl, dDotColl,
			nBur * sizeof(*dDotColl), cudaMemcpyDeviceToHost) );

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

	// Free GPU and host memory
	cudaAssert( cudaFree(dBurDotCoords) );
	cudaAssert( cudaFree(dAccDotCoords) );
	cudaAssert( cudaFree(dDotColl) );

	cudaAssert( cudaFreeHost(hBurDotCoords) );
	cudaAssert( cudaFreeHost(hAccDotCoords) );
	cudaAssert( cudaFreeHost(hDotColl) );

	VERBOSE("Peripheral trimming GPU processing time: " << GetTimerMs(timer) << " ms" << std::endl);

	return area;
}

int ShapeComplementarityCalculator::CudaFindClosestNeighbors(
	std::vector<DOT const*> const &my_dots,
	std::vector<DOT const*> const &their_dots,
	std::vector<DOT const*> &neighbors)
{
	int nMyDots, nTheirDots, nNeighbors;
	int threads;
	clock_t timer;

	timer = clock();
	threads = MIN(512, settings.gpu_threads);

	// Memory pointers for my and their dot coordinate arrays, CPU and GPU
	float3 *hMyDotCoords, *phMyDotCoords;
	float3 *hTheirDotCoords, *phTheirDotCoords;
	float3 *dMyDotCoords, *dTheirDotCoords;

	// Dot point pointer map
	DOT const **hTheirDots, **phTheirDots;

	// Neighbor ID memory pointers
	::uint *hNeighbors;
	::uint *dNeighbors;

	nMyDots = my_dots.size();
	nTheirDots = their_dots.size();
	nNeighbors = nMyDots;

	// Allocate host memory
	cudaAssert( cudaHostAlloc((void**)&hMyDotCoords, nMyDots * sizeof(*hMyDotCoords), cudaHostAllocDefault) );
	cudaAssert( cudaHostAlloc((void**)&hTheirDotCoords, nTheirDots * sizeof(*hTheirDotCoords), cudaHostAllocDefault) );
	cudaAssert( cudaHostAlloc((void**)&hTheirDots, nTheirDots * sizeof(*hTheirDots), cudaHostAllocDefault) );
	cudaAssert( cudaHostAlloc((void**)&hNeighbors, nNeighbors * sizeof(*hNeighbors), cudaHostAllocDefault) );

	// Make CUDA copy of (x, y, z) dot coordinates for my dots
	phMyDotCoords = hMyDotCoords;
	for(std::vector<DOT const*>::const_iterator idot = my_dots.begin(); idot < my_dots.end(); ++idot) {
#ifdef SC_PRECISION_REAL
		phMyDotCoords->x = (*idot)->coor.x();
		phMyDotCoords->y = (*idot)->coor.y();
		phMyDotCoords->z = (*idot)->coor.z();
		++phMyDotCoords;
#else
		// Quick copy
		*phMyDotCoords++ = *((float3*)&(*idot)->coor);
#endif
	}
	nMyDots = phMyDotCoords - hMyDotCoords;

	// Make CUDA copy of (x, y, z) dot coordinates for their dots and keep a map
	phTheirDotCoords = hTheirDotCoords;
	phTheirDots = hTheirDots;
	for(std::vector<DOT const*>::const_iterator idot = their_dots.begin(); idot < their_dots.end(); ++idot) {
		if(!(*idot)->buried)
			continue;
#ifdef SC_PRECISION_REAL
		phTheirDotCoords->x = *(*idot)->coor.x();
		phTheirDotCoords->y = *(*idot)->coor.y();
		phTheirDotCoords->z = *(*idot)->coor.z();
		++phTheirDotCoords;
#else
		// Quick copy
		*phTheirDotCoords++ = *((float3*)&(*idot)->coor);
#endif
		*phTheirDots++ = *idot;
	}
	nTheirDots = phTheirDotCoords - hTheirDotCoords;

	// Allocate GPU memory and copy data there
	cudaAssert( cudaMalloc((void **)&dMyDotCoords, UPPER_MULTIPLE(nMyDots, threads) * sizeof(*dMyDotCoords)) );
	cudaAssert( cudaMalloc((void **)&dTheirDotCoords, nTheirDots * sizeof(*dTheirDotCoords)) );
	cudaAssert( cudaMalloc((void **)&dNeighbors, UPPER_MULTIPLE(nNeighbors, threads) * sizeof(*dNeighbors)) );
	cudaAssert( cudaMemcpy(dMyDotCoords, hMyDotCoords, nMyDots * sizeof(*dMyDotCoords), cudaMemcpyHostToDevice) );
	cudaAssert( cudaMemcpy(dTheirDotCoords, hTheirDotCoords, nTheirDots * sizeof(*dTheirDotCoords), cudaMemcpyHostToDevice) );

	// Run kernel in multi-threaded blocks on GPU
	::_cuda_sccalc_FindClosestNeighbor(UPPER_MULTIPLE(nMyDots, threads)/threads, threads,
				dMyDotCoords,
				dTheirDotCoords, nTheirDots,
				dNeighbors);

	// Wait for threads to finish and copy back memory from GPU to host
	cudaAssert( cudaThreadSynchronize() );
	cudaAssert( cudaMemcpy(hNeighbors, dNeighbors, nMyDots * sizeof(*hNeighbors), cudaMemcpyDeviceToHost) );

	for(int i = 0; i < nNeighbors; ++i)
		neighbors.push_back( hTheirDots[hNeighbors[i]] );

	// Free memory
	cudaAssert( cudaFree(dMyDotCoords) );
	cudaAssert( cudaFree(dTheirDotCoords) );
	cudaAssert( cudaFree(dNeighbors) );

	cudaAssert( cudaFreeHost(hMyDotCoords) );
	cudaAssert( cudaFreeHost(hTheirDotCoords) );
	cudaAssert( cudaFreeHost(hTheirDots) );
	cudaAssert( cudaFreeHost(hNeighbors) );

	VERBOSE("Find Neighbors GPU processing time: " << GetTimerMs(timer) << " ms" << std::endl);

	return 1;
}

#endif // GPU

////////////////////////////////////////////////////////////////////////////
// Private helpers
////////////////////////////////////////////////////////////////////////////

Atom::Atom() : numeric::xyzVector < ShapeComplementarityCalculator::ScValue > (0.0)
{
	natom = 0;
	nresidue = 0;
	molecule = 0;
	radius = 0;
	density = 0;
	atten = 0;
	access = 0;

	memset(atom, 0, sizeof(atom));
	memset(residue, 0, sizeof(residue));
}

Atom::~Atom() {}

// The End
////////////////////////////////////////////////////////////////////////////

} // namespace sc
} // namespace scoring
} // namespace core

#endif

// END //
