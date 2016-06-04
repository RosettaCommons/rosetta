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

#ifndef INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_HH
#define INCLUDED_protocols_loops_FoldTreeFromLoopsWrapper_HH

// Unit headers
#include <protocols/loops/FoldTreeFromLoopsWrapper.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/loops/Loops.fwd.hh>

namespace protocols {
namespace loops {

class FoldTreeFromLoops : public protocols::moves::Mover
{
public:
	FoldTreeFromLoops();
	virtual ~FoldTreeFromLoops();

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap &,
		Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const &
	);

	void loop_str( std::string const & str );

	std::string loop_str() const;

	void loops( LoopsOP const l );

	LoopsOP loops() const;

private:
	std::string loop_str_; // loaded at parsetime but only realized at apply
	LoopsOP loops_; // a different interface into FoldTreeFromLoops, which takes precedence over loop_str_;
};

} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_FoldTreeFromLoops_HH
