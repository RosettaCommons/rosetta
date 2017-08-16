// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/task_operations/StoreResidueSubsetMover.hh
/// @brief Headers for StoreResidueSubsetMover class
/// @author Tom Linsky ( tlinsky at uw dot edu )

#ifndef INCLUDED_protocols_residue_selectors_StoreResidueSubsetMover_hh
#define INCLUDED_protocols_residue_selectors_StoreResidueSubsetMover_hh

//unit headers
#include <protocols/residue_selectors/StoreResidueSubsetMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Core Headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

namespace protocols {
namespace residue_selectors {

/// @brief mover that can be used to save or restore a task at an arbitrary
/// point during a rosetta scripts protocol. other task operations, movers,
/// or filters can be set up to access tasks saved by this mover during their
/// apply calls.
class StoreResidueSubsetMover : public protocols::moves::Mover {

public:
	StoreResidueSubsetMover();
	StoreResidueSubsetMover(
		core::select::residue_selector::ResidueSelectorCOP selector,
		std::string subset_name, // move-constructed
		bool const overwrite );
	~StoreResidueSubsetMover() override;

	void apply( core::pose::Pose & pose  ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data_map,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	std::string subset_name_;
	bool overwrite_;
};


} // residue_selectors
} // protocols

#endif

