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

/// @brief Definition of the search functions declared in SearchOption.hh
/// @author jk
/// @author Andrea Bazzoli

#include <devel/constel/MasterFilter.hh>
#include <devel/constel/Primitives.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AA.hh>
#include <core/types.hh>

namespace devel {
namespace constel {

///
/// @brief Search by target residue.
///
/// @details Extracts all the constellations formed by a target residue. Each
/// 	constellation involves a pair of residues: the target residue and one of
/// 	its neighbors. Different constellations are printed to different files.
///
/// @param[in] target_pdb_number residue number of the target residue in the
/// 	input PDB file.
/// @param[in] target_pdb_chain chain to which the target belongs in the PDB
/// 	file.
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
	core::Size constellation_number = 0;
	char aat = core::chemical::oneletter_code_from_aa(pose_init.aa(target_rosetta_resnum));
	utility::vector1<char> allowable_target_mutations =	list_allowable_mutations(aat);

	for ( Size j = 1; j <= pose_init.total_residue(); ++j )
		if ( interacting_residue.at(j) ) {

			char aaj = core::chemical::oneletter_code_from_aa(pose_init.aa(j));
			int j_pdb_number = pose_init.pdb_info()->number(j);
			char j_pdb_chain = pose_init.pdb_info()->chain(j);
			utility::vector1<char> allowable_secondary_mutations = list_allowable_mutations(aaj);

			for ( Size tmut=1; tmut <= allowable_target_mutations.size(); ++tmut ) {

				char aa_tmut = allowable_target_mutations.at(tmut);
				Pose target_mut_pose = pose_init;
				zero_occ_for_deleted_atoms( target_mut_pose, target_rosetta_resnum, aa_tmut);

				for ( Size jmut=1; jmut <= allowable_secondary_mutations.size(); ++jmut ) {

					char aa_jmut = allowable_secondary_mutations.at(jmut);
					Pose secondary_mut_pose = target_mut_pose;
					zero_occ_for_deleted_atoms( secondary_mut_pose, j, aa_jmut);

					utility::vector1<Size> cnl;
					cnl.push_back(target_rosetta_resnum);
					cnl.push_back(j);
					if( MasterFilter::is_constel_valid( secondary_mut_pose, cnl ) ) {
						++constellation_number;

						// print the atoms that would be removed by these mutations to a pdb file
						ResMut mut1(aat, aa_tmut, target_pdb_chain, target_pdb_number, target_rosetta_resnum);
						ResMut mut2(aaj, aa_jmut, j_pdb_chain, j_pdb_number, j);
						out_pair_constel(mut1, mut2, constellation_number, secondary_mut_pose);
					}
				}
			}
		}
}


///
/// @brief Search by pair of amino acid mutations.
///
/// @details Extracts from a pose all the constellations that correspond to a
/// 	given pair of mutations of amino acid types.
///
/// @param[in] tgtmuts an AB_CD string identifying the desired pair of
/// 	mutations, (A->B, C->D), where A, B, C, and D are amino acid types
/// 	expressed in one-letter code.
/// @param[in] pose_init pose from which constellations have to be extracted.
///
void pair_constel_set( std::string const& tgtmuts, Pose& pose_init ) {

	Size const TOTRES = pose_init.total_residue();

	NeighTeller nt(pose_init);

	char saa1 = tgtmuts[0];
	char saa2 = tgtmuts[3];

	for( Size i=1; i<TOTRES; ++i ) {

		char aai = core::chemical::oneletter_code_from_aa( pose_init.aa( i ) );

		if( ( aai == saa1 ) || ( aai == saa2 ) ) {

			core::conformation::Residue resi = pose_init.residue(i);
			char chi = pose_init.pdb_info()->chain(i);
			Size pdbi = pose_init.pdb_info()->number(i);

			for( Size j=i+1; j<=TOTRES; ++j ) {

				char aaj = core::chemical::oneletter_code_from_aa( pose_init.aa( j ) );

				if( ( ( aaj == saa2 ) && ( aai == saa1) ) ||
				    ( ( aaj == saa1 ) && ( aai == saa2) ) ) {

					core::conformation::Residue resj = pose_init.residue(j);
					char chj = pose_init.pdb_info()->chain(j);
					Size pdbj = pose_init.pdb_info()->number(j);

					if( nt.isneigh( resi, resj, pose_init ) ) {

						char aaie = ( aai == saa1 ) ? tgtmuts[1] : tgtmuts[4];
						char aaje = ( aai == saa1 ) ? tgtmuts[4] : tgtmuts[1];

						Pose pose_mut = pose_init;
						zero_occ_for_deleted_atoms( pose_mut, i, aaie );
						zero_occ_for_deleted_atoms( pose_mut, j, aaje );
						utility::vector1<Size> cnl;
						cnl.push_back(i);
						cnl.push_back(j);
						if( MasterFilter::is_constel_valid( pose_mut, cnl ) ) {

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

} // constel
} // devel
