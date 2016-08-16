// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PlacementMinimizationMover.hh
/// @brief definition of a class for making the placement auction used by PlaceSimultaneouslyMover
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlacementMinimizationMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlacementMinimizationMover_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMover.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief a simple rb-minimization in a bb-stub constraint biased forcefield.
/// Note that this mover is dependent on a placement mover for setting its stubsets
class PlacementMinimizationMover : public simple_moves::DesignRepackMover
{
public:
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
public:
	PlacementMinimizationMover();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void refresh_bbstub_constraints( core::pose::Pose & pose );
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	virtual void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	//// mutators for Placement movers to copy their internals onto auctionMover
	void host_chain( core::Size const hc );
	//void max_cb_cb_dist( core::Real const mccd );
	void cb_force( core::Real const cf );
	void stub_sets( utility::vector1< StubSetStubPos > const & stub_sets );
	~PlacementMinimizationMover();
private:
	core::Size host_chain_;
	utility::vector1< StubSetStubPos > stub_sets_;
	core::Real cb_force_;
};

} //movers
} //protein_interface_design
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_PlacementMinimization_HH*/
