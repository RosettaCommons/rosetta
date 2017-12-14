// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/movers/SetupHotspotConstraintsMover.hh
/// @brief Derived classes from DockDesign for dock design
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_SetupHotspotConstraintsMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_SetupHotspotConstraintsMover_hh

// Project Headers
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class SetupHotspotConstraintsMover : public protocols::moves::Mover {

public:

	SetupHotspotConstraintsMover();
	SetupHotspotConstraintsMover(
		protocols::hotspot_hashing::HotspotStubSetCOP hotspot_stub_set,
		core::Size const chain_to_design,
		core::Real const & CB_force_constant,
		core::Real const & worst_allowed_stub_bonus,
		bool const apply_self_energies,
		core::Real const & bump_cutoff,
		bool const apply_ambiguous_constraints,
		bool const colonyE,
		std::string stub_energy_fxn
	);
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	SetupHotspotConstraintsMover( SetupHotspotConstraintsMover const & init );

	void apply( core::pose::Pose & pose ) override;

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	~SetupHotspotConstraintsMover();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_set_;
	core::Size chain_to_design_;
	core::Real CB_force_constant_;
	core::Real worst_allowed_stub_bonus_;
	bool apply_self_energies_;
	core::Real bump_cutoff_;
	bool apply_ambiguous_constraints_;
	bool colonyE_;
	std::string stub_energy_fxn_;
};

} // movers
} // protein_interface_design
} // devel

#endif /*INCLUDED_protocols_protein_interface_design_movers_SetupHotspotConstraintsMover_HH*/
