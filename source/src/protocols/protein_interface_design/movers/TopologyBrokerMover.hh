// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/TopologyBrokerMover.hh
/// @author Lei Shi (shilei@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_TopologyBrokerMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_TopologyBrokerMover_hh
#include <protocols/protein_interface_design/movers/TopologyBrokerMover.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>


//Auto Headers
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class TopologyBrokerMover : public calc_taskop_movers::DesignRepackMover
{
public:
	typedef core::pose::Pose Pose;
public:
	TopologyBrokerMover();
	void apply( Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return utility::pointer::make_shared< TopologyBrokerMover >(); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & ) override;
	~TopologyBrokerMover() override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	bool align_;
	// std::string target_filename_;
	// core::pose::PoseOP targetA_;
	// core::pose::PoseOP motif_;
	// std::string motif_filename_;
	core::Size start_;
	core::Size end_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_TopologyBrokerMover_HH*/
