// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Definition of the class to create single-residue constellations
/// @author jk
/// @author Andrea Bazzoli (bazzoli@ku.edu)
#include <devel/constel/SingResCnlCrea.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <sstream>
#include <stdlib.h>

static basic::Tracer TR("devel.constel.SingResCnlCrea");

namespace devel {
namespace constel {

using core::Size;

/// @brief True if constellations must be output to file in stripped form
bool SingResCnlCrea::stripped_;


/// @brief mutation-specific lists of constellation atoms whose occupancy must
/// 	be zeroed	when the constellation is output to file.
std::map<SingResCnlCrea::MutTyp, std::string> SingResCnlCrea::zeroed_occ_;


/// @brief initializes the static members of this class
/// @param[in] stripped true if one wants to output stripped constellations;
/// 	false if one wants to output full constellations.
///
void SingResCnlCrea::init(bool const stripped) {

	if(stripped) {

		stripped_ = true;

		zeroed_occ_[MutTyp('P','A')] = "CG CD";
		zeroed_occ_[MutTyp('C','A')] = "SG";
		zeroed_occ_[MutTyp('I','V')] = "CD1";
		zeroed_occ_[MutTyp('I','A')] = "CG1 CG2";
		zeroed_occ_[MutTyp('F','L')] = "CE1 CE2";
		zeroed_occ_[MutTyp('S','A')] = "OG";
		zeroed_occ_[MutTyp('T','S')] = "CG2";
		zeroed_occ_[MutTyp('T','A')] = "OG1 CG2";
		zeroed_occ_[MutTyp('Y','F')] = "OH";
		zeroed_occ_[MutTyp('Y','L')] = "CE1 CE2";
		zeroed_occ_[MutTyp('W','L')] = "CE2 CE3 NE1";
		zeroed_occ_[MutTyp('V','A')] = "CG1 CG2";
	}
	else
		stripped_ = false;
}


/// @brief Returns the list of amino acid types that a given amino acid
/// 	type can be reduced to.
///
/// @param[in] starting_aa the given amino acid type.
///
utility::vector1<char> SingResCnlCrea::list_allowable_mutations(
	char const starting_aa ) {

	utility::vector1<char> allowed_list;

	if ( starting_aa == 'G' ) return allowed_list;

	// general cases
	// anything other than Gly can become Gly, anything other than Gly/Ala can become Ala
	allowed_list.push_back('G');
	if ( starting_aa != 'A' ) {
		allowed_list.push_back('A');
	}

	// special cases
	// Thr can become Ser
	if ( starting_aa == 'T' ) {
		allowed_list.push_back('S');
	}
	// Ile can become Val
	if ( starting_aa == 'I' ) {
		allowed_list.push_back('V');
	}
	// Tyr can become Phe and Leu
	if ( starting_aa == 'Y' ) {
		allowed_list.push_back('F');
		allowed_list.push_back('L');
	}
	// Phe can become Leu
	if ( starting_aa == 'F' ) {
		allowed_list.push_back('L');
	}
	// Trp can become Leu
	if ( starting_aa ==	'W' ) {
		allowed_list.push_back('L');
	}

	return allowed_list;
}


/// @brief sets occupancy to 0 for a residue's backbone atoms and hydrogen
/// 	atoms. Sets occupancy to 1 for the residue's remaining atoms.
///
/// @param[in] ps pose to which the residue belongs.
/// @param[in] seqpos index of the residue in the pose.
///
void SingResCnlCrea::zero_occ_bb_h(Pose& ps, core::Size seqpos) {

	core::conformation::Residue const & rsd( ps.residue(seqpos) );

	// we never need to print backbone atoms or hydrogens
	for ( Size i=1; i <= rsd.natoms(); ++i ) {
		ps.pdb_info()->occupancy( seqpos, i, 1 ); // to overwrite (rare) zero occ.
		if ( rsd.atom_is_hydrogen(i) || rsd.atom_is_backbone(i) ) {
			ps.pdb_info()->occupancy( seqpos, i, 0 );
		}
	}
}


/// @brief Sets occupancy to zero for a residue's non-constellation atoms.
///
/// @details Given a residue to be reduced to a smaller amino acid type, sets
/// 	occupancy to zero for all atoms that are not to be printed (those
/// 	forming the new residue too) and sets occupancy to 1 for all atoms
/// 	that are to be printed (those forming the constellation).
///
/// @param[out] pose pose to which the residue belongs.
/// @param[in] seqpos index of the residue in the pose.
/// @param[in] target_aa amino acid type that the residue has to be mutated
///     into.
///
void SingResCnlCrea::zero_occ_for_deleted_atoms(Pose & pose, core::Size seqpos,
	char const target_aa) {

	char const starting_aa = core::chemical::oneletter_code_from_aa( pose.aa(seqpos) );

	core::conformation::Residue const & rsd( pose.residue(seqpos) );

	zero_occ_bb_h(pose, seqpos);

	// if it's a mutation to Gly, we're done
	if ( target_aa == 'G' ) return;

	// for anything other than a mutation to Gly, suppress the C-beta
	pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CB"), 0. );

	// if it's a mutation to Ala, we're done
	if ( target_aa == 'A' ) return;

	// I-->V: suppress everything other than CD1
	if ( ( starting_aa == 'I' ) && ( target_aa == 'V' ) ) {
		Size atom_inx_to_keep = rsd.atom_index("CD1");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// T-->S: suppress everything other than CG2
	if ( ( starting_aa == 'T' ) && ( target_aa == 'S' ) ) {
		Size atom_inx_to_keep = rsd.atom_index("CG2");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// Y-->F: suppress everything other than OH
	if ( ( starting_aa == 'Y' ) && ( target_aa == 'F' ) ) {
		Size atom_inx_to_keep = rsd.atom_index("OH");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
    		return;
  	}

	// Y-->L: suppress everything other than OH, CE1, CE2, CZ
	if ( ( starting_aa == 'Y' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("OH");
		Size inx2 = rsd.atom_index("CE1");
		Size inx3 = rsd.atom_index("CE2");
		Size inx4 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i !=inx4) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// F-->L: suppress everything other than CE1, CE2, CZ
	if ( ( starting_aa == 'F' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("CE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
    		return;
  	}

	// W-->L: suppress everything other than NE1, CE2, CE3, CZ2, CZ3, CH2
	if ( ( starting_aa == 'W' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("NE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CE3");
		Size inx4 = rsd.atom_index("CZ2");
		Size inx5 = rsd.atom_index("CZ3");
		Size inx6 = rsd.atom_index("CH2");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i != inx4) && (i != inx5) && (i != inx6) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	TR << "ERROR: cannot mutate " << starting_aa << " into " << target_aa
		<< std::endl;
	exit(1);
}


/// @brief sets occupancy to zero for constellation atoms that are not
/// 	to be printed on output.
///
/// @param[out] pose pose to which the residue belongs.
/// @param[in] seqpos index of the residue in the pose.
/// @param[in] target_aa amino acid type that the residue has to be mutated
///     into.
///
/// @details: it is assumed that function zero_occ_for_deleted_atoms(pose,
/// 	seqpos, target_aa) was previously called.
///
void SingResCnlCrea::strip_atoms(Pose & pose, core::Size seqpos,
	char const target_aa) {

	char const starting_aa = core::chemical::oneletter_code_from_aa( pose.aa(seqpos) );

	core::conformation::Residue const & rsd( pose.residue(seqpos) );

	// if target is Gly, always suppress CB.
	if(target_aa == 'G') {
		pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CB"), 0. );
		if(starting_aa == 'P')
			pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CD"), 0. );
		return;
	}

	// if target is Ala, suppress CB and CG for several amino acid types
	if(target_aa == 'A') {
		switch(starting_aa) {
			case 'R':
			case 'N':
			case 'D':
			case 'Q':
			case 'E':
			case 'H':
			case 'L':
			case 'K':
			case 'M':
			case 'F':
			case 'W':
			case 'Y':
				pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CG"), 0. );
				return;
		}
	}

	// if it's another kind of mutation, suppress its corresponding set of atoms
	MutTyp mt(starting_aa, target_aa);
	std::map<MutTyp, std::string>::const_iterator mti = zeroed_occ_.find(mt);
	if(mti != zeroed_occ_.end()) {
		std::stringstream zatoms(mti->second);
		std::string atom;
		while(zatoms >> atom) {
			pose.pdb_info()->occupancy(seqpos, rsd.atom_index(atom), 0. );
		}
	}
	else {
		TR << "ERROR: cannot mutate " << starting_aa << " into " << target_aa
			<< std::endl;
		exit(1);
	}
}

} // namespace constel
} // namespace devel
