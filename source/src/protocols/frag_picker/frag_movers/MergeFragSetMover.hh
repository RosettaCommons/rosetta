// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/frag_picker/fragment_movers/MergeFragSetMover
///
/// @brief      Merges two fragment sets
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note
#ifndef INCLUDED_protocols_frag_picker_frag_movers_MergeFragSetMover_hh
#define INCLUDED_protocols_frag_picker_frag_movers_MergeFragSetMover_hh

// Project Headers
//
#include <protocols/frag_picker/frag_movers/MergeFragSetMover.fwd.hh>

#include <basic/datacache/DataMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>

#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
// C++ Headers
#include <string>
#include <map>

namespace protocols {
namespace frag_picker {
namespace frag_movers {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
using namespace fragment;

class MergeFragSetMover : public protocols::moves::Mover {
public:
	MergeFragSetMover();
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new MergeFragSetMover( *this ) ); }
	//virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
};

}//frag_movers
}//frag_picker
}//protocols

#endif
