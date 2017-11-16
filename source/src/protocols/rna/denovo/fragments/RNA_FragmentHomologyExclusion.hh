// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author  watkins

#ifndef INCLUDED_protocols_farna_RNA_FragmentHomologyExclusion_HH
#define INCLUDED_protocols_farna_RNA_FragmentHomologyExclusion_HH

#include <protocols/rna/denovo/fragments/RNA_Fragments.fwd.hh>
#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.hh>
#include <core/types.hh>
#include <set>

namespace protocols {
namespace rna {
namespace denovo {
namespace fragments {

class RNA_FragmentHomologyExclusion {
public:
	RNA_FragmentHomologyExclusion( RNA_Fragments const & all_rna_fragments );
	RNA_FragmentHomologyExclusion( RNA_FragmentHomologyExclusion const & rhs ) : 
		fragment_lines_( rhs.fragment_lines_ )
	{}
	RNA_FragmentHomologyExclusion() {}

	std::set< core::Size > const & get_fragment_lines() const { return fragment_lines_; }

private:
	std::set< core::Size > fragment_lines_;

};


// For map key
inline bool operator<( RNA_FragmentHomologyExclusion const & lhs, RNA_FragmentHomologyExclusion const & rhs ) {
	// Bad but ok for now.
	return lhs.get_fragment_lines().size() < rhs.get_fragment_lines().size();
}

} //fragments
} //denovo
} //rna
} //protocols

#endif

