// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief implementation of class ChainTerm
/// @author Andrea Bazzoli

#include <devel/constel/ChainTerm.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <iomanip>

using std::setw;


namespace devel {
namespace constel {


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
/// @details chains[i] represents the ith chain in the pose (i=1,...,N, where N
///  is the number of chains in the pose).
///
void get_chain_terms(core::pose::Pose const &ps, utility::vector1<ChainTerm> &chains) {

	Size const TOTRES = ps.size();

	Size nps = 1;
	int npdb = ps.pdb_info()->number(nps);
	char cid = ps.pdb_info()->chain(nps);
	for ( Size i=2; i<=TOTRES; ++i ) {
		char rcid = ps.pdb_info()->chain(i);
		if ( rcid != cid ) {
			chains.push_back(ChainTerm(cid, nps, i-1, // add chain
				npdb, ps.pdb_info()->number(i-1)));
			cid = rcid;
			nps = i;
			npdb = ps.pdb_info()->number(nps);
		}
	}
	chains.push_back(ChainTerm(cid, nps, TOTRES, // add last chain
		npdb, ps.pdb_info()->number(TOTRES)));
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

	for ( Size i=1; i<=chains.size(); ++i ) {
		chains[i].print(t);
	}
}


} // constel
} // devel
