// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Definition of class FilterBySASA
/// @author Andrea Bazzoli

#include <devel/constel/FilterBySASA.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <map>


namespace devel {
namespace constel {

/// @brief A table listing, for each amino acid type, the atoms whose SASA value
/// 	is relevant to filtering.
std::map< char, utility::vector1<std::string> > FilterBySASA::aa_sasa_atoms;

/// @brief A table holding the sasa values of all atoms in the pose to which
/// 	constellations belong.
core::id::AtomID_Map<Real> FilterBySASA::atom_sasa;

/// @brief Maximum allowed sasa value for a constellation atom.
double FilterBySASA::MAX_ATOM_SASA;


///
/// @brief Initializes the data structures needed to filter out constellations
/// 	based on SASA.
///
/// @param[in] smax maximum allowed SASA value for an atom.
/// @param[in] ps pose to which constellations belong.
///
void FilterBySASA::init( Real const smax, Pose const& ps ) {

	MAX_ATOM_SASA = smax;

	//// initialize table of atoms relevant to filtering
	const std::string AA1 = "ACDEFGHIKLMNPQRSTVWY";

	utility::vector1<std::string> v;
	v.push_back("All");
	for(unsigned int i=0; i<AA1.length(); ++i)
		aa_sasa_atoms[ AA1[i] ] = v;
/*
	// Asp
	v.clear();
	v.push_back("OD1");
	v.push_back("OD2");
	aa_sasa_atoms['D'] = v;

	// Glu
	v.clear();
	v.push_back("OE1");
	v.push_back("OE2");
	aa_sasa_atoms['E'] = v;

	// Lys
	v.clear();
	v.push_back("NZ");
	aa_sasa_atoms['K'] = v;*/


	//// compute SASA for all atoms in the pose
	utility::vector1<Real> rsd_sasa( ps.n_residue(), 0.0 );
	core::scoring::calc_per_atom_sasa( ps, atom_sasa, rsd_sasa, 1.4);
}


///
/// @brief Tells whether a constellation has a sufficiently low per-atom SASA.
///
/// @param[in] ps pose to which all residues in the constellation belong.
/// @param[in] cnl indexes in the pose of the residues forming the
/// 	constellation.
///
/// @return true if the constellation has a sufficiently low per-atom SASA;
/// 	false otherwise.
///
/// @details In the current implementation, the function returns true if each
/// 	atom whose SASA is relevant to filtering has a SASA lower than or equal to
/// 	a given cutoff, specified in MAX_ATOM_SASA; the function returns false
/// 	otherwise.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
///   	the atoms that belong to the constellation. This is guaranteed if the
/// 		residues forming the constellation had their indexes previously passed
/// 		as arguments to function "SingResCnlCrea::zero_occ_for_deleted_atoms()".
/// 	2. Pose ps has the same per-atom solvent accessibility as the pose passed
///    	to function init().
///
bool FilterBySASA::has_low_per_atom_sasa(Pose const& ps,
	utility::vector1<Size> const& cnl) {

	for(Size i=1; i <= cnl.size(); ++i) {

		// get residue info
		Size const ri = cnl[i];
		core::conformation::Residue const & rsd( ps.residue(ri) );
		char const aat = core::chemical::oneletter_code_from_aa( ps.aa(ri) );
		utility::vector1<std::string> const atoms = aa_sasa_atoms[ aat ];

		// look for atoms with too high a sasa
		if(atoms[1] == "All") {
			for( Size j=1; j<=rsd.nheavyatoms(); ++j )
				if( ps.pdb_info()->occupancy( ri, j ) ) {
					core::id::AtomID aid(j, ri);
					if( atom_sasa[aid] > MAX_ATOM_SASA )
						return false;
				}
		}
		else
			for( Size j=1; j<=atoms.size(); ++j ) {
				Size ai = rsd.atom_index( atoms[j] );
				if( ps.pdb_info()->occupancy( ri, ai ) ) {
					core::id::AtomID aid(ai, ri);
					if( atom_sasa[aid] > MAX_ATOM_SASA )
						return false;
				}
			}
	}

	return true;
}

} // constel
} // devel
