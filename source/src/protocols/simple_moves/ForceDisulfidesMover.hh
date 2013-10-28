// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   ForceDisulfidesMover.hh
/// @brief
/// @author Sarel Fleishman

#ifndef INCLUDED_protocols_simple_moves_ForceDisulfidesMover_hh
#define INCLUDED_protocols_simple_moves_ForceDisulfidesMover_hh

// Unit headers
#include <protocols/simple_moves/ForceDisulfidesMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {


/// @brief simple mover that fixes disulfides according to a defined list and then simultaneously repacks within 6A shells around each affected cystein residue.
class ForceDisulfidesMover : public moves::Mover {

public:
	ForceDisulfidesMover();
	~ForceDisulfidesMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;

	void disulfides( utility::vector1< std::pair< core::Size, core::Size > > );
	utility::vector1< std::pair< core::Size, core::Size > > disulfides() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP s );
private:
	utility::vector1< std::pair< core::Size, core::Size > > disulfides_;
	core::scoring::ScoreFunctionOP scorefxn_;// for repacking
};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_ForceDisulfidesMover_HH
