// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_src_devel_blab_classic_frags_TorsionFragment_FWD_HH
#define INCLUDED_src_devel_blab_classic_frags_TorsionFragment_FWD_HH


// Rosetta Headers

// ObjexxFCL Headers

// Utility Headers
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ Headers

namespace protocols {
namespace frags {

class TorsionFragment;
typedef utility::pointer::shared_ptr< TorsionFragment > TorsionFragmentOP;

class SingleResidueTorsionFragmentLibrary;
typedef utility::pointer::shared_ptr< SingleResidueTorsionFragmentLibrary > SingleResidueTorsionFragmentLibraryOP;

class TorsionFragmentLibrary;
typedef utility::pointer::shared_ptr< TorsionFragmentLibrary > TorsionFragmentLibraryOP;
typedef utility::pointer::shared_ptr< TorsionFragmentLibrary const > TorsionFragmentLibraryCOP;

class FragLib;
typedef utility::pointer::shared_ptr< FragLib > FragLibOP;
typedef utility::pointer::shared_ptr< FragLib const > FragLibCOP;

class TorsionFragmentMover;
typedef utility::pointer::shared_ptr< TorsionFragmentMover > TorsionFragmentMoverOP;


} // ns frags
} // ns protocols
#endif
