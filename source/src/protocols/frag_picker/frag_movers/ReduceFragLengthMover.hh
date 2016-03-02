// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/frag_picker/frag_movers/ReduceFragLengthMover
///
/// @brief      Converts a larger frag set such as 9mer to a smaller one such as 3mers
/// @details
///
/// @author     TJ Brunette (tjbrunette@gmail.com)
/// @note
#ifndef INCLUDED_protocols_frag_picker_frag_movers_ReduceFragLengthMover_hh
#define INCLUDED_protocols_frag_picker_frag_movers_ReduceFragLengthMover_hh

// Project Headers
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

using namespace core;
using namespace fragment;

class ReduceFragLengthMover : public protocols::moves::Mover {
public:
	ReduceFragLengthMover();
	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	moves::MoverOP clone() const { return moves::MoverOP( new ReduceFragLengthMover( *this ) ); }
	//virtual void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
    virtual void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
};

}//frag_movers
}//frag_picker
}//protocols

#endif
