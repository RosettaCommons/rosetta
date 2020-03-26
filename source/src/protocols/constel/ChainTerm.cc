// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/ChainTerm.cc
/// @brief implementation of class ChainTerm
/// @author Andrea Bazzoli

#include <protocols/constel/ChainTerm.hh>
#include <core/io/util.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <iomanip>

using std::setw;

namespace protocols {
namespace constel {

using core::Size;

/// @brief Prints this chain's info to tracer 't'.
///
/// @param[out] t output tracer.
///
void ChainTerm::print(basic::Tracer& t) const {

	t << setw(4) << cid_ << setw(8) << n_ps_ << setw(8) << c_ps_
		<< setw(8) << n_pdb_ << setw(8) << c_pdb_ << std::endl;
}


/// @brief identifies each chain and its N- and C-terminal residues in a pose.
///
/// @param[in] ps the pose.
/// @param[out] chains vector to be filled with ChainTerm items for the chains
///  in the pose.
///
/// @details chains[i] represents the ith chain in the pose built using only the
///  protein residues of pose ps (i=1,...,N, where N is the number of such chains)
///
void get_chain_terms(core::pose::Pose const &ps, utility::vector1<ChainTerm> &chains) {

	// create vector of protein residues (protein[i] <-- ith protein residue in the pose,
	// for i=1,...,TOTRES, where TOTRES is the number of protein residues in the pose)
	utility::vector1<core::Size> protein;
	for ( core::Size i=1; i<=ps.size(); ++i ) {
		if ( ps.residue(i).is_protein() ) {
			protein.push_back(i);
		}
	}
	if ( protein.size() < 1 ) {
		utility_exit_with_message( "get_chain_terms(): pose does not contain any protein residue" );
		exit(1);
	}

	// collect chain termini: chains[i] stores the termini of the ith chain in protein_pose
	// (i=1,...,N, where N is the number of such chains)
	core::pose::Pose protein_pose;
	core::io::pose_from_pose(protein_pose, ps, protein);
	for ( core::Size i=1; i<=protein_pose.num_chains(); ++i ) {

		core::Size chain_begin = protein_pose.chain_begin(i);
		char beg_cid = protein_pose.pdb_info()->chain(chain_begin);
		int beg_num = protein_pose.pdb_info()->number(chain_begin);
		char beg_ico = protein_pose.pdb_info()->icode(chain_begin);
		core::Size chain_begin_orig = ps.pdb_info()->pdb2pose(beg_cid, beg_num, beg_ico);

		core::Size chain_end = protein_pose.chain_end(i);
		int end_num = protein_pose.pdb_info()->number(chain_end);
		char end_ico = protein_pose.pdb_info()->icode(chain_end);
		core::Size chain_end_orig = ps.pdb_info()->pdb2pose(beg_cid, end_num, end_ico);

		chains.push_back(ChainTerm(beg_cid, chain_begin_orig, chain_end_orig, beg_num, end_num));
	}
}


/// @brief Prints the N- and C-terminal residues of all chains.
///
/// @param[in] chains vector of ChainTerm items, each representing a
///     different chain.
/// @param[out] t output tracer
///
/// @details The ith output line represents the ith chain in 'chains',
///     for i=1,...,N, where N is the size of 'chains'.
///
void print_chains(utility::vector1<ChainTerm> const &chains, basic::Tracer &t) {

	t << setw(4) << "cid_" << setw(8) << "n_ps_" << setw(8) << "c_ps_" <<
		setw(8) << "n_pdb_" << setw(8) << "c_pdb_" << std::endl;

	for ( core::Size i=1; i<=chains.size(); ++i ) {
		chains[i].print(t);
	}
}


} // constel
} // protocols
