// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief class identifying the N- and C-terminal residues of a chain
/// @author Andrea Bazzoli

#ifndef INCLUDED_ChainTerm_hh
#define INCLUDED_ChainTerm_hh

#include <core/pose/Pose.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>


namespace devel {
namespace constel {

using core::Size;

class ChainTerm {

	char cid_; // chain id

	Size n_ps_; // index of N-terminal residue in pose
	Size c_ps_; // index of C-terminal residue in pose

	int n_pdb_; // index of N-terminal residue in PDB file
	int c_pdb_; // index of C-terminal residue in PDB file

public:

	/// @brief constructor
	ChainTerm(char cid, Size nps, Size cps, int npdb, int cpdb) :
		cid_(cid), n_ps_(nps), c_ps_(cps), n_pdb_(npdb), c_pdb_(cpdb) {}

	/// @brief prints this chain's info to tracer t
	void print(basic::Tracer& t) const;

	/// @brief cid_ accessor
	char get_cid() const {return cid_;}

	/// @brief n_ps_ accessor
	Size get_nps() const {return n_ps_;}

	/// @brief c_ps_ accessor
	Size get_cps() const {return c_ps_;}
};


/// @brief identifies each chain and its N- and C-terminal residues in a pose.
void get_chain_terms(core::pose::Pose const &ps, utility::vector1<ChainTerm> &chains);

/// @brief prints the N- and C-terminal residues of all chains.
void print_chains(utility::vector1<ChainTerm> const &chains, basic::Tracer &t);


} // constel
} // devel

#endif
