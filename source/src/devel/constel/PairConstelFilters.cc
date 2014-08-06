// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Definition of filters for constellations formed by a pair of residues.
/// @author Andrea Bazzoli

#include <devel/constel/PairConstelFilters.hh>
#include <devel/constel/Primitives.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <string>

using core::Real;


namespace devel {
namespace constel {

///
/// @brief tells whether a constellation satisfies the requirements for being
/// 	 rescued by compounds containing an indole group and a carboxylic group.
///
/// @return true if the constellation satisfies such requirements; false
/// 	otherwise.
///
/// @param[in] ps pose to which all residues in the constellation belong.
/// @param[in] cnl indexes in the pose of the residues forming the
/// 	constellation.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
///			the atoms that belong to the constellation. This is guaranteed if the
///			residues forming the constellation had their indexes previously passed
/// 		as arguments to function "zero_occ_for_deleted_atoms()".
///
/// @details
/// 	1. For constellations formed by other than 2 residues, the function
/// 		currently returns false.
///
bool FilterByIndoleCOO::sat(Pose const& ps,
                            utility::vector1<Size> const& cnl) {

	using core::chemical::AA;
	using core::chemical::aa_trp;
	using core::chemical::aa_asp;
	using core::chemical::aa_glu;
	using core::conformation::Residue;
	using utility::vector1;

	if(cnl.size() != 2)
		return false;

	// check that one residue is Trp and the other is Asp or Glu
	vector1<AA> aa_typs(2);
	vector1<Size> aa_idxs(2);
	aa_typs[1] = aa_trp;
	aa_typs[2] = aa_asp;
	if(!PresenceCommon::are_aa_pres(ps, cnl, aa_typs, aa_idxs)) {
		aa_typs[2] = aa_glu;
		if(!PresenceCommon::are_aa_pres(ps, cnl, aa_typs, aa_idxs))
			return false;
	}
	Size const iw = aa_idxs[1];
	Size const io = aa_idxs[2];

	// check that Trp is not being mutated into Leu
	Residue const& wres = ps.residue(cnl[iw]);
	vector1<std::string> wpoi_nam(2);
	wpoi_nam[1] = "CB";
	wpoi_nam[2] = "CD1";
	vector1<Size> wpoi_idx(2);
	if(!PresenceCommon::are_atoms_pres(wres, wpoi_nam, wpoi_idx))
		return false;
	if(!ps.pdb_info()->occupancy( cnl[iw], wpoi_idx[2]))
		return false;

	// check that Trp's point of interest is further away from the other
	// residues's oxygens than from the other residue's branching carbon
	numeric::xyzVector<Real> wpoi_xyz = wres.xyz(wpoi_idx[1]);

	Residue const& ores = ps.residue(cnl[io]);
	vector1<std::string> o_nams(3);
	vector1<Size> o_idxs(3);
	if(aa_typs[2] == aa_asp) {
		o_nams[1] = "CG";
		o_nams[2] = "OD1";
		o_nams[3] = "OD2";
	}
	else {
		o_nams[1] = "CD";
		o_nams[2] = "OE1";
		o_nams[3] = "OE2";
	}
	if(!PresenceCommon::are_atoms_pres(ores, o_nams, o_idxs))
		return false;
	numeric::xyzVector<Real> obc_xyz = ores.xyz(o_idxs[1]);
	numeric::xyzVector<Real> oo1_xyz = ores.xyz(o_idxs[2]);
	numeric::xyzVector<Real> oo2_xyz = ores.xyz(o_idxs[3]);

	bool pts_ok = OrientCommon::is_closer_to_tgt(obc_xyz, oo1_xyz, wpoi_xyz) &&
	              OrientCommon::is_closer_to_tgt(obc_xyz, oo2_xyz, wpoi_xyz);

	// check that the other residue is hydrogen bonded
	vector1<std::string> hb_atoms;
	if(aa_typs[2] == aa_asp) {
		hb_atoms.push_back("OD1");
		hb_atoms.push_back("OD2");
	}
	else {
		hb_atoms.push_back("OE1");
		hb_atoms.push_back("OE2");
	}
	return pts_ok && HBondCommon::is_rmoi_hbonded(ps, cnl, io, false, hb_atoms);
}


///
/// @brief Tells whether a constellation satisfies the requirements for being
/// 	 rescued by tryptamine.
///
/// @return true if the constellation satisfies such requirements; false
/// 	otherwise.
///
/// @param[in] ps pose to which all residues in the constellation belong.
/// @param[cnl] indexes in the pose of the residues forming the constellation.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
/// 		the atoms that belong to the constellation. This is guaranteed if the
/// 		residues forming the constellation had their indexes previously passed
/// 		as arguments to function "zero_occ_for_deleted_atoms()".
///
/// @details
/// 	1. For constellations formed by other than 2 residues, the function
/// 	currently returns false.
///
bool FilterByTryptamine::sat(Pose const& ps,
	utility::vector1<Size> const& cnl) {

	using core::chemical::AA;
	using core::chemical::aa_trp;
	using core::chemical::aa_lys;
	using core::conformation::Residue;
	using utility::vector1;

	if( cnl.size() != 2 )
		return false;

	// check that one residue is Trp and the other is Lys
	vector1<AA> aa_typs(2);
	vector1<Size> aa_idxs(2);
	aa_typs[1] = aa_trp;
	aa_typs[2] = aa_lys;
	if(!PresenceCommon::are_aa_pres(ps, cnl, aa_typs, aa_idxs))
		return false;
	Size const iw = aa_idxs[1];
	Size const ik = aa_idxs[2];

	// check that Trp is not being mutated into Leu
	Residue const& wres = ps.residue(cnl[iw]);
	vector1<std::string> wpoi_nam(2);
	wpoi_nam[1] = "CG";
	wpoi_nam[2] = "CD1";
	vector1<Size> wpoi_idx(2);
	if(!PresenceCommon::are_atoms_pres(wres, wpoi_nam, wpoi_idx))
		return false;
	if(!ps.pdb_info()->occupancy( cnl[iw], wpoi_idx[2]))
		return false;

	// check that Trp's point of interest is farther away from Lys's NZ than from
	// Lys's CE
	numeric::xyzVector<Real> wpoi_xyz = wres.xyz(wpoi_idx[1]);

	Residue const& kres = ps.residue(cnl[ik]);
	vector1<std::string> k_nams(2);
	vector1<Size> k_idxs(2);
	k_nams[1] = "CE";
	k_nams[2] = "NZ";
	if(!PresenceCommon::are_atoms_pres(kres, k_nams, k_idxs))
		return false;
	numeric::xyzVector<Real> kce_xyz = kres.xyz(k_idxs[1]);
	numeric::xyzVector<Real> knz_xyz = kres.xyz(k_idxs[2]);

	if(!OrientCommon::is_closer_to_tgt(kce_xyz, knz_xyz, wpoi_xyz))
		return false;

	// check that Lys's NZ donates at least one hydrogen bond
	vector1<std::string> hb_atoms;
	hb_atoms.push_back("1HZ");
	hb_atoms.push_back("2HZ");
	hb_atoms.push_back("3HZ");
	return HBondCommon::is_rmoi_hbonded(ps, cnl, ik, true, hb_atoms);
}


///
/// @brief Tells whether a constellation satisfies the requirements for being
/// 	 rescued by amphetamine.
///
/// @return true if the constellation satisfies such requirements; false
/// 	otherwise.
///
/// @param[in] ps pose to which all residues in the constellation belong.
/// @param[cnl] indexes in the pose of the residues forming the constellation.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
/// 		the atoms that belong to the constellation. This is guaranteed if the
/// 		residues forming the constellation had their indexes previously passed
/// 		as arguments to function "zero_occ_for_deleted_atoms()".
///
/// @details
/// 	1. For constellations formed by other than 2 residues, the function
/// 	currently returns false.
///
bool FilterByAmphetamine::sat(Pose const& ps,
                              utility::vector1<Size> const& cnl) {

	using core::chemical::AA;
	using core::chemical::aa_phe;
	using core::chemical::aa_lys;
	using core::conformation::Residue;
	using utility::vector1;

	if( cnl.size() != 2 )
		return false;

	// check that one residue is Phe and the other is Lys
	vector1<AA> aa_typs(2);
	vector1<Size> aa_idxs(2);
	aa_typs[1] = aa_phe;
	aa_typs[2] = aa_lys;
	if(!PresenceCommon::are_aa_pres(ps, cnl, aa_typs, aa_idxs))
		return false;
	Size const ip = aa_idxs[1];
	Size const ik = aa_idxs[2];

	// check that Phe is not being mutated into Leu
	Residue const& pres = ps.residue(cnl[ip]);
	vector1<std::string> ppoi_nam(2);
	ppoi_nam[1] = "CG";
	ppoi_nam[2] = "CD1";
	vector1<Size> ppoi_idx(2);
	if(!PresenceCommon::are_atoms_pres(pres, ppoi_nam, ppoi_idx))
		return false;
	if(!ps.pdb_info()->occupancy( cnl[ip], ppoi_idx[2]))
		return false;

	// check that Phe's point of interest is farther away from Lys's NZ than from
	// Lys's CE
	numeric::xyzVector<Real> ppoi_xyz = pres.xyz(ppoi_idx[1]);

	Residue const& kres = ps.residue(cnl[ik]);
	vector1<std::string> k_nams(2);
	vector1<Size> k_idxs(2);
	k_nams[1] = "CE";
	k_nams[2] = "NZ";
	if(!PresenceCommon::are_atoms_pres(kres, k_nams, k_idxs))
		return false;
	numeric::xyzVector<Real> kce_xyz = kres.xyz(k_idxs[1]);
	numeric::xyzVector<Real> knz_xyz = kres.xyz(k_idxs[2]);

	if(!OrientCommon::is_closer_to_tgt(kce_xyz, knz_xyz, ppoi_xyz))
		return false;

	// check that Lys's NZ donates at least one hydrogen bond
	vector1<std::string> hb_atoms;
	hb_atoms.push_back("1HZ");
	hb_atoms.push_back("2HZ");
	hb_atoms.push_back("3HZ");
	return HBondCommon::is_rmoi_hbonded(ps, cnl, ik, true, hb_atoms);
}


///
/// @brief Tells whether a constellation satisfies the requirements for being
/// 	 rescued by histamine.
///
/// @return true if the constellation satisfies such requirements; false
/// 	otherwise.
///
/// @param[in] ps pose to which all residues in the constellation belong.
/// @param[cnl] indexes in the pose of the residues forming the constellation.
///
/// @remarks It is assumed that:
/// 	1. The residues forming the constellation have non-zero occupancy only for
/// 		the atoms that belong to the constellation. This is guaranteed if the
/// 		residues forming the constellation had their indexes previously passed
/// 		as arguments to function "zero_occ_for_deleted_atoms()".
///		2. All residues in the pose have all heavy atoms present (whether original
/// 		or reconstructed).
///
/// @details
/// 	1. For constellations formed by other than 2 residues, the function
/// 	currently returns false.
///
bool FilterByHistamine::sat(Pose const& ps, utility::vector1<Size> const& cnl) {

	using core::chemical::AA;
	using core::chemical::aa_his;
	using core::chemical::aa_lys;
	using core::conformation::Residue;
	using utility::vector1;

	if( cnl.size() != 2 )
		return false;

	// check that one residue is His and the other is Lys
	vector1<AA> aa_typs(2);
	vector1<Size> aa_idxs(2);
	aa_typs[1] = aa_his;
	aa_typs[2] = aa_lys;
	if(!PresenceCommon::are_aa_pres(ps, cnl, aa_typs, aa_idxs))
		return false;
	Size const ih = aa_idxs[1];
	Size const ik = aa_idxs[2];

	// check that His's point of interest is farther away from Lys's NZ than from
	// Lys's CE
	numeric::xyzVector<Real> hpoi_xyz = ps.residue(cnl[ih]).xyz("CD2");
	numeric::xyzVector<Real> kce_xyz = ps.residue(cnl[ik]).xyz("CE");
	numeric::xyzVector<Real> knz_xyz = ps.residue(cnl[ik]).xyz("NZ");

	if(!OrientCommon::is_closer_to_tgt(kce_xyz, knz_xyz, hpoi_xyz))
		return false;

	// check that Lys's NZ donates at least one hydrogen bond
	vector1<std::string> hb_atoms;
	hb_atoms.push_back("1HZ");
	hb_atoms.push_back("2HZ");
	hb_atoms.push_back("3HZ");
	return HBondCommon::is_rmoi_hbonded(ps, cnl, ik, true, hb_atoms);
}


} // constel
} // devel 
