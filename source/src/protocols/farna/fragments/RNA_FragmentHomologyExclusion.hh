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

#include <protocols/farna/fragments/RNA_Fragments.fwd.hh>
#include <protocols/farna/fragments/FullAtomRNA_Fragments.hh>
#include <core/types.hh>
#include <set>

namespace protocols {
namespace farna {
namespace fragments {

class RNA_FragmentHomologyExclusion {
public:
	RNA_FragmentHomologyExclusion( RNA_Fragments const & all_rna_fragments );
	
	std::set< core::Size > const & get_fragment_lines() const { return fragment_lines_; }
					
private:
	std::set< core::Size > fragment_lines_;

};

}
}
}

#endif

