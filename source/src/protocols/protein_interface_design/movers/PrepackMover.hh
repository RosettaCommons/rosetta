// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PrepackMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PrepackMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PrepackMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class PrepackMover : public protocols::simple_moves::PackRotamersMover
{
public:
	PrepackMover();
	PrepackMover( core::scoring::ScoreFunctionCOP scorefxn, core::Size jump_num );
	virtual ~PrepackMover();

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	bool min_bb() const;
	void min_bb( bool const m );
	core::kinematics::MoveMapOP mm() const;
	void mm( core::kinematics::MoveMapOP mm );
private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Size jump_num_;
	bool min_bb_; //dflt false
	core::kinematics::MoveMapOP mm_; // only activated if min_bb is on.
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_PrepackMover_HH*/
