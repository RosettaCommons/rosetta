// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_rbsegment_relax_make_star_topology_hh
#define INCLUDED_protocols_rbsegment_relax_make_star_topology_hh

#include <protocols/rbsegment_relax/MakeStarTopologyCreator.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/electron_density/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>


namespace protocols {
namespace rbsegment_relax {

class MakeStarTopologyMover : public moves::Mover {
public:
	MakeStarTopologyMover() : Mover(), mode_(""), restore_(false) {}

	// XRW TEMP  std::string get_name() const override { return MakeStarTopologyMoverCreator::mover_name(); }
	moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new MakeStarTopologyMover( *this ) ) ); }

	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	std::string mode_;

	// restore a perturbed FT
	bool restore_;
	std::string tag_;
	core::kinematics::FoldTreeOP ft_restore_;
};


}
}

#endif
