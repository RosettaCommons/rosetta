// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Definition of the search functions declared in SearchOption.hh
/// @author jk
/// @author Andrea Bazzoli

#include <devel/constel/SingResCnlCrea.hh>
#include <devel/constel/MasterFilter.hh>
#include <devel/constel/NeighTeller.hh>
#include <devel/constel/Primitives.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <iostream>
#include <fstream>

namespace devel {
namespace constel {

static THREAD_LOCAL basic::Tracer TR( "src.devel.constel.SearchOptions" );


/// @brief Searches pair-constellations by target residue.
///
/// @details Extracts all the constellations formed by a target residue. Each
///  constellation involves a pair of residues: the target residue and one of
///  its neighbors. Different constellations are printed to different files.
///
/// @param[in] target_pdb_number residue number of the target residue in the
///  input PDB file.
/// @param[in] target_pdb_chain chain to which the target belongs in the PDB
///  file.
/// @param[in] pose_init pose to which the target belongs.
///
void pair_constel_set(int const target_pdb_number, char const target_pdb_chain,
	Pose& pose_init) {

	// set target_rosetta_resnum to Rosetta internal resid for the residue to be mutated
	core::Size target_rosetta_resnum = get_pose_resnum(target_pdb_number, target_pdb_chain, pose_init);

	// make a list of residues with vdw atr contacts (to the sidechain only, if possible?)
	utility::vector1<bool> interacting_residue( pose_init.total_residue(), false );
	mk_neigh_list(target_rosetta_resnum, interacting_residue, pose_init);

	// for each (target, neighbor) pair, loop over all allowed mutation combinations
	for ( Size j = 1; j <= pose_init.total_residue(); ++j ) {
		if ( interacting_residue.at(j) ) {
			pair_constel_set_idx2(target_rosetta_resnum, j, pose_init);
		}
	}
}


/// @brief Searches pair-constellations by pair of mutations.
///
/// @details Extracts from a pose all the constellations that correspond to a
///  given pair of mutations of amino acid types.
///
/// @param[in] tgtmuts an AB_CD string identifying the desired pair of
///  mutations, (A->B, C->D), where A, B, C, and D are amino acid types
///  expressed in one-letter code.
/// @param[in] pose_init pose from which constellations have to be extracted.
///
void pair_constel_set( std::string const& tgtmuts, Pose& pose_init ) {

	Size const TOTRES = pose_init.total_residue();

	NeighTeller nt(pose_init);

	char saa1 = tgtmuts[0];
	char saa2 = tgtmuts[3];

	for ( Size i=1; i<TOTRES; ++i ) {

		char aai = core::chemical::oneletter_code_from_aa( pose_init.aa( i ) );

		if ( ( aai == saa1 ) || ( aai == saa2 ) ) {

			core::conformation::Residue resi = pose_init.residue(i);
			char chi = pose_init.pdb_info()->chain(i);
			Size pdbi = pose_init.pdb_info()->number(i);

			for ( Size j=i+1; j<=TOTRES; ++j ) {

				char aaj = core::chemical::oneletter_code_from_aa( pose_init.aa( j ) );

				if ( ( ( aaj == saa2 ) && ( aai == saa1) ) ||
						( ( aaj == saa1 ) && ( aai == saa2) ) ) {

					core::conformation::Residue resj = pose_init.residue(j);
					char chj = pose_init.pdb_info()->chain(j);
					Size pdbj = pose_init.pdb_info()->number(j);

					if ( nt.isneigh( resi, resj, pose_init ) ) {

						char aaie = ( aai == saa1 ) ? tgtmuts[1] : tgtmuts[4];
						char aaje = ( aai == saa1 ) ? tgtmuts[4] : tgtmuts[1];

						Pose pose_mut = pose_init;
						SingResCnlCrea::zero_occ_for_deleted_atoms( pose_mut, i, aaie );
						SingResCnlCrea::zero_occ_for_deleted_atoms( pose_mut, j, aaje );
						utility::vector1<Size> cnl;
						cnl.push_back(i);
						cnl.push_back(j);
						if ( MasterFilter::is_constel_valid( pose_mut, cnl ) ) {

							ResMut muti( aai, aaie, chi, pdbi, i );
							ResMut mutj( aaj, aaje, chj, pdbj, j );

							out_pair_constel( muti, mutj, -1, pose_mut );
						}
					}
				}
			}
		}
	}
}


/// @brief Searches for the triple-constellations of a target residue.
///
/// @details Each constellation comprises a triple of spatially contiguous
///  residues, among which is the target. Different constellations are printed
///  to different files.
///
/// @param[in] target_pdb_number residue number of the target residue in the input
///  PDB file.
/// @param[in] target_pdb_chain chain to which the target belongs in the PDB file.
/// @param[in] pose_init: pose to which the target belongs.
///
void triple_constel_set(int const target_pdb_number,
	char const target_pdb_chain, Pose& pose_init) {

	using utility::vector1;
	using core::chemical::oneletter_code_from_aa;
	using core::Size;

	// set target_rosetta_resnum to Rosetta internal resid for target residue
	Size target_rosetta_resnum = get_pose_resnum(target_pdb_number,
		target_pdb_chain, pose_init);

	// make a list of residues with vdw atr contacts with target
	vector1<bool> interacting_residue( pose_init.total_residue(), false);
	mk_neigh_list(target_rosetta_resnum, interacting_residue, pose_init);

	// search for (target, neighbor, neighbor) triples
	Size UJ = pose_init.total_residue() - 1;
	Size UK = UJ+1;

	for ( Size j = 1; j <= UJ; ++j ) {
		if ( interacting_residue.at(j) ) {
			for ( Size k = j+1; k <= UK; ++k ) {
				if ( interacting_residue.at(k) ) {
					triple_constel_set_idx3(target_rosetta_resnum, j, k, pose_init );
				}
			}
		}
	}

	// search for (target, neighbor, other) triples
	NeighTeller ngbtel(pose_init);
	UJ = UK;
	for ( Size j = 1; j <= UJ; ++j ) {
		if ( interacting_residue.at(j) ) {
			core::conformation::Residue const& rj = pose_init.residue(j);
			for ( Size k = 1; k <= UK; ++k ) {
				if ( k!=target_rosetta_resnum ) {
					if ( !interacting_residue.at(k) ) {
						if ( ngbtel.isneigh(rj, pose_init.residue(k), pose_init) ) {
							if ( j<k ) {
								triple_constel_set_idx3(target_rosetta_resnum, j, k, pose_init );
							} else {
								triple_constel_set_idx3(target_rosetta_resnum, k, j, pose_init );
							}
						}
					}
				}
			}
		}
	}
}


/// @brief Searches for a target constellation
///
/// @param[in] tgtcnl_fil file specifying the target constellation.
/// @param[in] ps the pose
///
/// @details The format of the input file is as follows:
///
///  CID{1} RNU{1} ICO{1} AASTA{1} AAEND{1}
///  ...
///  CID{N} RNU{N} ICO{N} AASTA{N} AAEND{N}
///
///  where N is the number of residues forming the constellation and CID{i},
///  RNU{i}, ICO{i}, AASTA{i}, and AAEND{i} are the chain id, residue number,
///  insertion code, start amino acid type, and end amino acid type of the ith
///  residue contributing to the constellation (i=1,...,N). A blank chain
///  identifier in the PDB file is denoted by ',' (comma); a blank insertion
///  code in the PDB file is denoted by '_' (underscore).
///
/// @details The function exits if any input lines cannot identify any residue.
///
void target_constel(std::string &tgtcnl_fil, Pose & ps) {

	using core::chemical::oneletter_code_from_aa;

	// identify constellation according to input file
	std::ifstream ifs(tgtcnl_fil.c_str());
	if ( !ifs ) {
		TR << "can't open " << tgtcnl_fil << std::endl;
		return;
	}

	char cid;
	int rnum;
	char ico;
	char aasta;
	char aaend;

	Size const NRES = ps.total_residue();
	core::pose::PDBInfoCOP pdb_info = ps.pdb_info();

	// cnl[i] and resmut[i] contain info about the contributing residue specified
	// on input line i (i=1,...,N, where N is the number of lines; not tested).
	utility::vector1<Size> cnl;
	utility::vector1<ResMut> resmut;

	while ( ifs >> cid ) { // (proved for any number of lines)

		ifs >> rnum >> ico >> aasta >> aaend;

		if ( cid == ',' ) {
			cid = ' ';
		}

		if ( ico == '_' ) {
			ico = ' ';
		}

		bool found=false;

		// find residue that makes current contribution
		for ( Size i=1; i<=NRES; ++i ) {
			if (
					(pdb_info->number(i) == rnum) &&
					(pdb_info->chain(i) == cid) &&
					(pdb_info->icode(i) == ico) &&
					(oneletter_code_from_aa( ps.aa( i )) == aasta)
					) {
				SingResCnlCrea::zero_occ_for_deleted_atoms(ps, i, aaend);
				cnl.push_back(i);
				resmut.push_back(ResMut(aasta, aaend, cid, rnum, i));
				found = true;
				break;
			}
		}

		if ( !found ) {
			TR << "can't find residue contribution to constellation: " <<
				cid << ' ' << rnum << ' ' << ico << ' ' << aasta << std::endl;
			return;
		}
	}

	// if constellation is valid, print it to file
	Size const N = cnl.size();
	if ( (N != 2) && (N != 3) ) {
		TR << "can only handle constellations of 2 or 3 residues" <<
			std::endl;
		return;
	}

	NeighTeller nt(ps);
	if ( !nt.is_neigh_tree(cnl, ps) ) {
		return;
	}

	if ( MasterFilter::is_constel_valid(ps, cnl) ) {
		if ( N==2 ) {
			out_pair_constel(resmut[1], resmut[2], -1, ps);
		} else {
			out_triple_constel(resmut[1], resmut[2], resmut[3], -1, ps);
		}
	}
}

} // constel
} // devel
