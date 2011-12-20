// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/FoldTreeFromLoops.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_loops_FoldTreeFromLoops_hh
#define INCLUDED_protocols_loops_FoldTreeFromLoops_hh
#include <protocols/loops/FoldTreeFromLoopsWrapper.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/loops/Loops.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace protocols {
namespace loops {

class FoldTreeFromLoops : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	FoldTreeFromLoops();
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const { return protocols::moves::MoverOP( new FoldTreeFromLoops( *this ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new FoldTreeFromLoops ); }
		void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~FoldTreeFromLoops();
	void loop_str( std::string const str ){ loop_str_ = str; }
	std::string loop_str() const { return loop_str_; }
	void loops( LoopsOP const l ) { loops_ = l; }
	LoopsOP loops() const{ return loops_; }
private:
	std::string loop_str_; // loaded at parsetime but only realized at apply
	LoopsOP loops_; // a different interface into FoldTreeFromLoops, which takes precedence over loop_str_;
};


} // loops
} // protocols


#endif /*INCLUDED_protocols_loops_FoldTreeFromLoops_HH*/
