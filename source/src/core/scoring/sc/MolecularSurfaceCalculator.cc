// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file     core/scoring/sc/ShapeComplementarityCalculator.cc
/// @brief    Implementation of molecular surface calculation for shape complementarity.
/// @detailed Lawrence & Coleman shape complementarity calculator (based on CCP4's sc)
/// @author   Luki Goldschmidt <luki@mbi.ucla.edu>, refactored by Alex Ford (fordas@uw.edu)

/// This code was ported from the original Fortran code found in CCP4:
/// Sc (Version 2.0): A program for determining Shape Complementarity
/// Copyright Michael Lawrence, Biomolecular Research Institute
/// 343 Royal Parade Parkville Victoria Australia
///
#ifndef INCLUDED_core_scoring_sc_MolecularSurfaceCalculator_cc
#define INCLUDED_core_scoring_sc_MolecularSurfaceCalculator_cc

// Project Headers
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.fwd.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C headers
#include <stdio.h>

// C++ headers
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

static basic::Tracer TR("core.scoring.sc.MolecularSurfaceCalculator");

using namespace core;

namespace core {
namespace scoring {
namespace sc {

std::vector<ATOM_RADIUS> MolecularSurfaceCalculator::radii_;  // static const

////////////////////////////////////////////////////////////////////////////
// Public class functions
////////////////////////////////////////////////////////////////////////////

/// @begin MolecularSurfaceCalculator::MolecularSurfaceCalculator()
/// @brief
/// MolecularSurfaceCalculator constructor, initializes default settings

MolecularSurfaceCalculator::MolecularSurfaceCalculator()
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

	Reset();
}

/// @begin MolecularSurfaceCalculator::Init()
/// @brief
/// Initializes calculation and GPU (if used)
/// Init() is also called implicitly by the static CalcSc() function.

int MolecularSurfaceCalculator::Init()
{
	if(radii_.empty())
		ReadScRadii();
	if(radii_.empty())
		return 0;

	return 1;
}

MolecularSurfaceCalculator::~MolecularSurfaceCalculator()
{
	Reset();
}

/// @begin MolecularSurfaceCalculator::Reset()
/// @brief
/// Reset calculator for another calculation.
/// Must be used when the MolecularSurfaceCalculator instance is re-used.
/// @detailed
/// Atom, probe and surface dot vectors are reset here. We don't clear them
/// after the calculation is finished in case the caller would like to use those
/// data elsewhere.

void MolecularSurfaceCalculator::Reset()
{
	// Free data
	for(std::vector<Atom>::iterator a = run_.atoms.begin(); a < run_.atoms.end(); ++a) {
		a->neighbors.clear();
		a->buried.clear();
	}

	// Clear structures
	run_.atoms.clear();
	run_.probes.clear();
	run_.dots[0].clear();
 	run_.dots[1].clear();

	memset(&run_.results, 0, sizeof(run_.results));
	memset(&run_.prevp, 0, sizeof(run_.prevp));
	run_.prevburied = 0;
}

/// @begin MolecularSurfaceCalculator::Calc(core::pose::Pose const & pose, core::Size jump_id = 0)
/// @brief Generate molecular surfaces for the given pose.
///// @detailed
/// This function initializes the calculator, adds all residues in the given pose, and generates molecular surfaces.
///
/// The pose is partitioned into separate molecules across the given jump. If the given jump is 0, the entire pose is
/// loaded as molecule 1.
/// To control what residues make up either surface, use the AddResidue() or even AddAtom() function instead.
/// Returns true on success. Results are retrieved with GetResults().
///
/// Example:
/// core::scoring::sc::MolecularSurfaceCalculator calc;
/// if(calc.Calc( pose ))
///   ... = calc.GetResults();
int MolecularSurfaceCalculator::Calc(core::pose::Pose const & pose, core::Size jump_id)
{
	if(!Init())
		return 0;

	if( jump_id > pose.num_jump() )
	{
		TR.Error << "Jump ID out of bounds (pose has " << pose.num_jump() << " jumps)" << std::endl;
		return 0;
	}

	// Partition pose by jump_id
	ObjexxFCL::FArray1D_bool is_upstream ( pose.total_residue(), true );

	if( jump_id > 0 )
	{
		pose.fold_tree().partition_by_jump( jump_id, is_upstream );
	}

	for(Size i = 1; i <= pose.n_residue(); ++i) {
		core::conformation::Residue const & residue = pose.residue(i);
		if(residue.type().name() == "VRT")
			continue;
		if(!AddResidue(is_upstream(i) ? 0 : 1, residue))
			return 0;
	}

	return Calc();
}

/// @begin MolecularSurfaceCalculator::Calc()
/// @brief Generate molecular surfaces for loaded atoms.
///// @detailed
/// This function generates molecular surfaces for atoms added via AddAtom and AddResidue.
///
/// Init() must be called before this function.
/// Returns true on success.
int MolecularSurfaceCalculator::Calc()
{
	try
	{
		basic::gpu::Timer timer(TR.Debug);

		run_.results.valid = 0;

		if(run_.atoms.empty())
			throw ShapeComplementarityCalculatorException("No atoms defined");

		// Determine assign the attention numbers for each atom
		AssignAttentionNumbers(run_.atoms);

		GenerateMolecularSurfaces();
		return 1;
	}
	catch(ShapeComplementarityCalculatorException e)
	{
		TR.Error << "Failed: " << e.error << std::endl;
	}

	return 0;
}


/// @begin MolecularSurfaceCalculator::GenerateMolecularSurfaces
/// @brief Generate untrimmed surfaces for the defined molecules.
///// @detailed
/// This function should be called within a try/catch block for ShapeComplementarityCalculatorException.
/// Raises exception on error.
void MolecularSurfaceCalculator::GenerateMolecularSurfaces()
{
	if(run_.atoms.empty())
  {
		throw ShapeComplementarityCalculatorException("No atoms defined");
  }

	// Now compute the surface for the atoms in the interface and its neighbours
	TR.Debug << "Generating molecular surface, " << settings.density << " dots/A^2" << std::endl;
	CalcDotsForAllAtoms(run_.atoms);

	if(TR.Debug.visible())
  {
		TR.Debug << "      Buried atoms (1): " << run_.results.surface[0].nBuriedAtoms << std::endl;
		TR.Debug << "      Buried atoms (2): " << run_.results.surface[1].nBuriedAtoms << std::endl;
		TR.Debug << "     Blocked atoms (1): " << run_.results.surface[0].nBuriedAtoms << std::endl;
		TR.Debug << "     Blocked atoms (2): " << run_.results.surface[1].nBuriedAtoms << std::endl;
		TR.Debug << "           Convex dots: " << run_.results.dots.convex << std::endl;
		TR.Debug << "         Toroidal dots: " << run_.results.dots.toroidal << std::endl;
		TR.Debug << "          Concave dots: " << run_.results.dots.concave << std::endl;
		TR.Debug << "Total surface dots (1): " << run_.dots[0].size() << std::endl;
		TR.Debug << "Total surface dots (2): " << run_.dots[1].size() << std::endl;
		TR.Debug << "    Total surface dots: " << (run_.dots[0].size()+run_.dots[1].size()) << std::endl;
	}
}

/// @begin MolecularSurfaceCalculator::ReadScRadii()
/// @brief Read atom radius definitions from file
/// @defailed
/// This function is implicitly called, but can be overloaded or
/// called explicitly for custom handling of the atom radii library.
/// Returns true on success

int MolecularSurfaceCalculator::ReadScRadii()
{
	char const *fn = "scoring/score_functions/sc/sc_radii.lib";
	ATOM_RADIUS radius;
	utility::io::izstream in;

	if(!basic::database::open(in, fn)) {
		TR.Error << "Failed to read " << fn << std::endl;
		return 0;
	}

	radii_.clear();

	while( in.good() ) {
		memset(&radius, 0, sizeof(radius));
		in >> radius.residue >> radius.atom >> radius.radius;
		if(*radius.residue && *radius.atom && radius.radius > 0) {
			TR.Trace << "Atom Radius: " << radius.residue << ":" << radius.atom << " = " << radius.radius << std::endl;
			radii_.push_back(radius);
		}
	}

	TR.Trace << "Atom radii read: " << radii_.size() << std::endl;

	return !radii_.empty();
}

/// @begin MolecularSurfaceCalculator::AddResidue()
/// @brief Add a rosetta residue to a specific molecule
/// @detailed
/// Call this function when explicitly defining which residues belong to
/// which the molecular surface. If partitioning by jump_id is sufficient
/// for your application, you may use the Calc() or CalcSc() functions
/// instead.
/// Returns number of atoms added for the specified residue.
///
/// Example:
/// core::scoring::sc::MolecularSurfaceCalculator calc;
/// core::Real sc;
/// calc.Init();
/// calc.Reset(); // Only needed when re-using the calculator
/// for(core::Size i = 1; i <= pose.n_residue(); i++)
///   calc.AddResidue((i < 100), pose.residue(i));
/// if(calc.Calc())
///   sc = calc.GetResults().sc;

core::Size MolecularSurfaceCalculator::AddResidue(
	int molecule,
	core::conformation::Residue const & residue)
{
	std::vector<Atom> scatoms;

	if(!Init())
		return 0;

	// Pass 1: Assign atom radii and check if we can add all atoms for this residue
	// Only use heavy atoms for SC calculation
	for(Size i = 1; i <= residue.nheavyatoms(); ++i) {
		// Skip virtual atoms
		if(residue.is_virtual(i))
			continue;

		Atom scatom;
		numeric::xyzVector<Real> xyz = residue.xyz(i);
		scatom.x(xyz.x());
		scatom.y(xyz.y());
		scatom.z(xyz.z());
		scatom.nresidue = residue.seqpos();
		scatom.radius = 0;
		strncpy(scatom.residue, residue.name3().c_str(), sizeof(scatom.residue)-1);
		strncpy(scatom.atom, residue.atom_name(i).c_str()+1, sizeof(scatom.atom)-1);

		if(!AssignAtomRadius(scatom)) {
			TR.Error << "Failed to add residue " << residue.name3() << " to surface - cannot find radius for " << residue.atom_name(i) << std::endl;
			return 0;
		}
		scatoms.push_back(scatom);
	}

	// Pass 2: Add all atoms for the residue
	int n =0;
	for(std::vector<Atom>::iterator it = scatoms.begin(); it != scatoms.end(); ++it) {
#if defined(WIN32) && !defined(WIN_PYROSETTA) // windows WINAPI has an AddAtom function
		if(AddAtomWIN32(molecule, *it)) ++n;
#else
		if(AddAtom(molecule, *it)) ++n;
#endif
	}

	return n;
}

/// @begin MolecularSurfaceCalculator::AddAtom()
/// @brief Add an atom to a molecule for computation.
/// @detailed
/// Add an core::scoring::sc::Atom to the molecule.
/// Normally this is called by AddResidue(). Explicit addition
/// of atoms via this function is rarely needed.
/// This function also looks-up the atom radius and density.
/// Returns true on success.
#if defined(WIN32) && !defined(WIN_PYROSETTA)  // windows WINAPI has an AddAtom function
int MolecularSurfaceCalculator::AddAtomWIN32(int molecule, Atom &atom)
#else
int MolecularSurfaceCalculator::AddAtom(int molecule, Atom &atom)
#endif
{
	if(atom.radius <= 0)
		AssignAtomRadius(atom);

	if(atom.radius > 0) {
		molecule = (molecule == 1);
		atom.density = settings.density;
		atom.molecule = molecule;
		atom.natom = ++run_.results.nAtoms;
		atom.access = 0;

		run_.atoms.push_back(atom);
		++run_.results.surface[molecule].nAtoms;
		return 1;
	} else {
		TR.Error << "Failed to assign atom radius for residue "
			<< atom.residue << ":" << atom.atom << "!" << std::endl;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
// Protected class functions
////////////////////////////////////////////////////////////////////////////

// Look up the atom radius for an atom
int MolecularSurfaceCalculator::AssignAtomRadius(Atom &atom)
{
	std::vector<ATOM_RADIUS>::const_iterator radius;

	// Assign radius with wildcard matching
	for(radius = radii_.begin(); radius != radii_.end(); ++radius) {
		if(WildcardMatch(atom.residue, radius->residue, sizeof(atom.residue)) &&
			WildcardMatch(atom.atom, radius->atom, sizeof(atom.atom))) {
				atom.radius = radius->radius;
				if(TR.Trace.visible()) {
					char buf[256];
#ifdef WIN32
					_snprintf(buf, sizeof(buf),
#else
					snprintf(buf, sizeof(buf),
#endif
						"Assigned atom radius to %s:%s at (%8.4f, %8.4f, %8.4f) = %.3f",
						atom.residue,
						atom.atom,
						atom.x(),
						atom.y(),
						atom.z(),
						atom.radius
					);
					TR.Trace << buf << std::endl;
				}
				return 1;
		}
	}

	return 0;
}

// Inline residue and atom name matching function
int MolecularSurfaceCalculator::WildcardMatch(
		char const *query,
		char const *pattern,
		int l)
{
	while(--l > 0) {
		bool match =
			(*query == *pattern) ||
			(*query && *pattern == '*') ||
			(*query == ' ' && !*pattern);
		if(!match)
			return 0;

		// Allow anything following a * in pattern
		if(*pattern == '*' && !pattern[1])
			return 1;

		if(*query)
			++query;
		if(*pattern)
			++pattern;
	}
	return 1;
}

// Assign attention numbers for each atom. By default attention numbers are assigned
// to compute surface dots using all defined atoms.
int MolecularSurfaceCalculator::AssignAttentionNumbers(std::vector<Atom> & atoms)
{
	std::vector<Atom>::iterator pAtom;

	for(pAtom = atoms.begin(); pAtom < atoms.end(); ++pAtom)
	{
		pAtom->atten = ATTEN_BURIED_FLAGGED;
		++run_.results.surface[pAtom->molecule].nBuriedAtoms;
	}

	return 1;
}

////////////////////////////////////////////////////////////////////////////
// Molecular surface calculation
////////////////////////////////////////////////////////////////////////////
// M. L. Connolly J. Appl. Crystallogr., 16, p548 - p558 (1983)
// Compute the surface for the atoms in the interface and its neighbours
////////////////////////////////////////////////////////////////////////////

int MolecularSurfaceCalculator::CalcDotsForAllAtoms(std::vector<Atom> & )
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
		{
			continue;
		}

		// Find neighbor
		if(!FindNeighbordsAndBuriedAtoms(atom1))
		{
			continue;
		}

		if(!atom1.access)
		{
			continue;
		}

		if(atom1.atten <= ATTEN_BLOCKER)
		{
			continue;
		}

		if(atom1.atten == ATTEN_6 && atom1.buried.empty())
		{
			continue;
		}

		// Generate convex surface
		GenerateConvexSurface(atom1);
	}

	// Concave surface generation
	if(settings.rp > 0)
		GenerateConcaveSurface();

	return 1;
}

/// @brief A small struct to report which of two atoms is closer to a given atom:
struct
CloserToAtom {
	bool operator() ( Atom * a1, Atom * a2 ) {
		MolecularSurfaceCalculator::ScValue d1 = ref_atom_->distance( * a1 );
		MolecularSurfaceCalculator::ScValue d2 = ref_atom_->distance( * a2 );
		return d1 < d2;
	}
	CloserToAtom( Atom * reference_atom ) {
		ref_atom_ = reference_atom;
	}
private:
	Atom * ref_atom_;
};

// Private atom distance sort callback used below
// Atom *_atom_distance_ref = NULL;
//int MolecularSurfaceCalculator::_atom_distance_cb(void *a1, void *a2)
//{
//	ScValue d1 = _atom_distance_ref->distance(*((Atom*)a1));
//	ScValue d2 = _atom_distance_ref->distance(*((Atom*)a2));
//	return d1 < d2 ? -1: (d2 > d1 ? 1 : 0);
//}

// Calculate surface dots around a single atom (main loop in original code)
int MolecularSurfaceCalculator::FindNeighbordsAndBuriedAtoms(Atom &atom1)
{
	if(!FindNeighborsForAtom(atom1))
		return 0;

	// sort neighbors by distance from atom1
	//_atom_distance_ref = &atom1;


	CloserToAtom closer_to_atom1( &atom1 );
	std::sort(atom1.neighbors.begin(), atom1.neighbors.end(), closer_to_atom1 );
	//_atom_distance_ref = NULL;

	SecondLoop(atom1);

	return atom1.neighbors.size();
}

// Make a list of neighboring atoms from atom1
// (loop to label 100 in original code)
int MolecularSurfaceCalculator::FindNeighborsForAtom(Atom &atom1)
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
int MolecularSurfaceCalculator::SecondLoop(Atom &atom1)
{
	Vec3 uij, tij;
	ScValue erj, eri, rij, /*density,*/ dij, asymm, _far_, contain;
	int between;
	std::vector<Atom*> &neighbors = atom1.neighbors;

	eri = atom1.radius + settings.rp;

	for(std::vector<Atom*>::iterator iAtom2 = neighbors.begin(); iAtom2 < neighbors.end(); ++iAtom2) {
		Atom &atom2 = **iAtom2;

		if(atom2 <= atom1)
			continue;

		erj = atom2.radius + settings.rp;
		//density = (atom1.density + atom2.density) / 2;  // set but never used ~Labonte
		dij = atom1.distance(atom2);

		uij = (atom2 - atom1) / dij;
		asymm = (eri * eri - erj * erj) / dij;
		between = (ABS(asymm) < dij);

		tij = ((atom1 + atom2) * 0.5f) + (uij * (asymm * 0.5f));

		_far_ = (eri + erj) * (eri + erj) - dij * dij;
		if (_far_ <= 0.0)
			continue;

		_far_ = sqrt(_far_);

		contain = dij * dij - ((atom1.radius - atom2.radius) * (atom1.radius - atom2.radius));
		if (contain <= 0.0)
			continue;

		contain = sqrt(contain);
		rij = 0.5 * _far_ * contain / dij;

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
int MolecularSurfaceCalculator::ThirdLoop(
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
int MolecularSurfaceCalculator::CheckAtomCollision2(
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
int MolecularSurfaceCalculator::GenerateConvexSurface(Atom const &atom1)
{
	std::vector<Atom*> const &neighbors = atom1.neighbors;
	Atom const * neighbor;
	Vec3 north(0, 0, 1);
	Vec3 south(0, 0, -1);
	Vec3 eqvec(1, 0, 0);
	Vec3 vtemp, vql, uij, tij, pij, cen, pcen;
	ScValue dt, erj, dij, eri, _far_, contain;
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
		_far_ = pow(eri + erj, 2) - dij*dij;
		if(_far_ <= 0.0)
			throw ShapeComplementarityCalculatorException("Imaginary _far_ for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		_far_ = sqrt(_far_);

		contain = pow(dij, 2) - pow(ri - rj, 2);
		if(contain <= 0.0)
			throw ShapeComplementarityCalculatorException("Imaginary contain for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		contain = sqrt(contain);
		rij = 0.5 * _far_ * contain / dij;
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
int MolecularSurfaceCalculator::CheckPointCollision(
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
int MolecularSurfaceCalculator::GenerateToroidalSurface(
		Atom &atom1,
		Atom &atom2,
		Vec3 const uij,
		Vec3 const tij,
		ScValue rij,
		int between)
{
	std::vector<Atom*> &neighbors = atom1.neighbors;
	ScValue density, /*ri,*/ /*rj,*/ rb, rci, rcj, rs, e, edens, eri, erj, erl, dtq, pcusp, /*anglei,*/ /*anglej,*/ dt, ts, ps, area;
	Vec3 pi, pj, axis, dij, pqi, pqj, qij, qjk, qj;

	std::vector<Vec3> subs;

	// following Fortran original
	// will be optimized by compiler
	density = (atom1.density + atom2.density) / 2;
	//ri = atom1.radius;  // set but never used ~Labonte
	//rj = atom2.radius;  // set but never used ~Labonte
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
		//anglei = acosf(dt);  // set but never used ~Labonte

		dt = pqj.dot(pj);
		if(dt >= 1.0 || dt <= -1.0)
			return 0;
		//anglej = acosf(dt);  // set but never used ~Labonte

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
int MolecularSurfaceCalculator::GenerateConcaveSurface()
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
int MolecularSurfaceCalculator::CheckProbeCollision(
		Vec3 const &point,
		std::vector<PROBE const *> const nears,
		ScValue const r2)
{
	for(std::vector<const PROBE*>::const_iterator _near_ = nears.begin();
			_near_ < nears.end(); ++_near_) {
		if(point.distance_squared((*_near_)->point) < r2)
			// Collision
			return 1;
	}
	return 0;
}

// Add a molecular dot
void MolecularSurfaceCalculator::AddDot(
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
MolecularSurfaceCalculator::ScValue MolecularSurfaceCalculator::DistancePointToLine(
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
MolecularSurfaceCalculator::ScValue MolecularSurfaceCalculator::SubArc(
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
MolecularSurfaceCalculator::ScValue MolecularSurfaceCalculator::SubDiv(
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
MolecularSurfaceCalculator::ScValue MolecularSurfaceCalculator::SubCir(
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

///////////////////////////////////////////////////////////////////////////
// Private helpers
////////////////////////////////////////////////////////////////////////////

Atom::Atom() : numeric::xyzVector < MolecularSurfaceCalculator::ScValue > (0.0)
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

#endif // INCLUDED_core_scoring_sc_MolecularSurfaceCalculator_cc

// END //
