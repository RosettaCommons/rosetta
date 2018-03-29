// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/BuildAlaPose.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_BuildAlaPose_hh
#define INCLUDED_protocols_protein_interface_design_movers_BuildAlaPose_hh
#include <protocols/protein_interface_design/movers/BuildAlaPose.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class BuildAlaPose : public calc_taskop_movers::DesignRepackMover
{
public:
	typedef core::pose::Pose Pose;
public:
	BuildAlaPose();
	BuildAlaPose( bool const partner1, bool const partner2, core::Real interface_distance_cutoff=8.0 , std::string AA="ALA");
	void apply( Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new BuildAlaPose ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	virtual ~BuildAlaPose();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real interface_distance_cutoff_;
	std::string AA_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_BuildAlaPose_HH*/
