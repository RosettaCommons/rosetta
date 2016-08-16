// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/SparseMatrixIndex.h
/// @brief  Sparse matrix index class to work with AminoAcidNeighborSparseMatrix class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_SparseMatrixIndex_hh
#define INCLUDED_core_pack_interaction_graph_SparseMatrixIndex_hh

namespace core {
namespace pack {
namespace interaction_graph {


class SparseMatrixIndex;

class SparseMatrixIndex
{
public:
	SparseMatrixIndex() : aa_type_( 0 ), state_ind_for_this_aa_type_( 0 ) {}

	void set_aa_type( int aa ) { aa_type_ = aa; }
	int get_aa_type() const    { return aa_type_;}

	void set_state_ind_for_this_aa_type(int state_ind ) { state_ind_for_this_aa_type_ = state_ind;}
	int  get_state_ind_for_this_aa_type() const         { return state_ind_for_this_aa_type_;}

protected:
	int aa_type_;
	int state_ind_for_this_aa_type_;
};

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif //INCLUDED_core_pack_interaction_graph_SparseMatrixIndex_HH
