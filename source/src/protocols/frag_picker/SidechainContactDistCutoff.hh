// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/SidechainContactDistCutoff.hh
/// @brief  Defines sidechain contact distance cutoffs.
/// @author David E. Kim (dekim@u.washington.edu)


#ifndef INCLUDED_protocols_frag_picker_SidechainContactDistCutoff_hh
#define INCLUDED_protocols_frag_picker_SidechainContactDistCutoff_hh

// unit headers
#include <protocols/frag_picker/SidechainContactDistCutoff.fwd.hh>

// project headers
#include <core/types.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <map>

namespace protocols {
namespace frag_picker {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

/// @brief defines sidechain contact distance cutoffs.
/// @details provides amino acid pair specific distance cutoffs.
class SidechainContactDistCutoff: public utility::pointer::ReferenceCount {
public:

	SidechainContactDistCutoff();

	SidechainContactDistCutoff(Real scale_factor);

	virtual ~SidechainContactDistCutoff(); // auto-removing definition from header{};

	Real get_cutoff(char aa_i, char aa_j);

	Real get_cutoff_squared(char aa_i, char aa_j);

	Real scale_factor();

private:

	void initialize();

private:
	utility::vector1<char> aa_map_;
	std::map<char,Size> aa_to_index_map_;
	utility::vector1<utility::vector1<Real> >  cutoff_;
	utility::vector1<utility::vector1<Real> >  cutoff_squared_;
	Real scale_factor_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_SidechainContactDistCutoff_HH */
