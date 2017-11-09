// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/LoopFinder.hh
/// @brief Header for class to find loops (parseable options to control how loops are found)
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_LoopFinder_hh
#define INCLUDED_protocols_protein_interface_design_movers_LoopFinder_hh

#include <core/pose/Pose.fwd.hh>
//#include <protocols/moves/Mover.hh>
#include <core/types.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// the following must be included because LoopFinder is derived from DRMover
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class LoopFinder : public simple_moves::DesignRepackMover
{
public:
	LoopFinder();
	LoopFinder(
		bool const interface,
		bool const ch1,
		bool const ch2,
		core::Size const min_length,
		core::Size const max_length,
		core::Size const mingap,
		std::string const & resnum,
		core::Real const ca_ca_distance,
		core::Real const iface_cutoff,
		protocols::loops::LoopsOP loops
	);
	virtual ~LoopFinder();

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new LoopFinder ); }
	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	Size interface_;
	bool ch1_, ch2_;
	core::Size min_length_, max_length_, mingap_;
	std::string resnum_;
	core::Real ca_ca_distance_;
	core::Real iface_cutoff_;
	//basic::datacache::DataMapOP data_;
	protocols::loops::LoopsOP loops_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_LoopFinder_HH*/
