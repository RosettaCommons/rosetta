// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/FlexbbSparseMatrixIndex.hh
/// @brief  Declaration of class for maintaining sparse matrix indices for flexbb states.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_FlexbbSparseMatrixIndex_hh
#define INCLUDED_protocols_flexpack_interaction_graph_FlexbbSparseMatrixIndex_hh

#include <core/pack/interaction_graph/SparseMatrixIndex.hh>

namespace protocols {
namespace flexpack {
namespace interaction_graph {

class FlexbbSparseMatrixIndex : public core::pack::interaction_graph::SparseMatrixIndex
{
public:
	typedef core::pack::interaction_graph::SparseMatrixIndex parent;

public:
	FlexbbSparseMatrixIndex() : parent(), bb_( 0 ) {}

	void set_bb( int bb) {
		bb_ = bb;
	}

	int get_bb() const {
		return bb_;
	}

private:
	int bb_;
};

}
}
}

#endif
