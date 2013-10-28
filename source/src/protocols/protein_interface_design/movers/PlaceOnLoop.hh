// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceOnLoop.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlaceOnLoop_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlaceOnLoop_hh
#include <protocols/protein_interface_design/movers/PlaceOnLoop.fwd.hh>


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {


class PlaceOnLoop : public simple_moves::DesignRepackMover
{
public:
	PlaceOnLoop();
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void parse_my_tag( utility::tag::TagCOP const tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	void set_kinematic_defaults();
	bool minimize_toward_stub( core::pose::Pose & pose ) const;
	void add_bb_csts_to_loop( core::pose::Pose & pose ) const;
	bool position_stub( core::pose::Pose & pose ) const;
	void ala_pose_loop( core::pose::Pose & pose ) const;
	bool loop_length( core::pose::Pose & pose );
	virtual ~PlaceOnLoop();
private:
	core::Size loop_begin_, loop_end_, curr_loop_end_;
	core::scoring::ScoreFunctionOP hires_scorefxn_, lores_scorefxn_;
	utility::vector1< int > delta_length_;
	core::Size chain_closing_attempts_;
	core::Size host_chain_;
	protocols::hotspot_hashing::HotspotStubSetOP stub_set_;
	bool minimize_toward_stub_;
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinematic_mover_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_PlaceOnLoop_HH*/
