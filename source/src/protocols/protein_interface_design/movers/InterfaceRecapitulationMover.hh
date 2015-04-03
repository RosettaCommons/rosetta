// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/InterfaceRecapitulationMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_InterfaceRecapitulationMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_InterfaceRecapitulationMover_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief a pure virtual base class for movers which redesign and repack the interface
class InterfaceRecapitulationMover : public protocols::moves::Mover
{
public:
	InterfaceRecapitulationMover();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const;
	void set_reference_pose( core::pose::PoseOP );
	core::pose::PoseCOP get_reference_pose() const;
	virtual void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~InterfaceRecapitulationMover();
private:
	core::pose::PoseOP saved_pose_;
	simple_moves::DesignRepackMoverOP design_mover_;
	protocols::simple_moves::PackRotamersMoverOP design_mover2_;//ugly adaptation for the PackRotamers baseclass
	bool pssm_;
};

} // movers
} // protein_interface_design
} // protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_InterfaceRecapitulationMover_HH*/
