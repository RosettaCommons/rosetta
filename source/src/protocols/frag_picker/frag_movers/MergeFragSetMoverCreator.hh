// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/frag_picker/fragment_movers/MergeFragSetMoverCreator
///
/// @brief      Merges two fragment sets
/// @details
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note

#ifndef INCLUDED_protocols_frag_picker_fragment_movers_MergeFragSetMoverCreator_hh
#define INCLUDED_protocols_frag_picker_fragment_movers_MergeFragSetMoverCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace frag_picker {
namespace frag_movers {

class MergeFragSetMoverCreator : public moves::MoverCreator
{
public:
	virtual moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static  std::string mover_name();
};

}//frag_movers
}//frag_picker
}//protocols

#endif
