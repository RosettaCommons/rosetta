// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/BackrubDDMover.hh
/// @brief Derived class from mover for dock design with backrub moves
/// @author Sarel Fleishman (sarelf@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_BackrubDDMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_BackrubDDMover_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


// C++ headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

class BackrubDDMover : public simple_moves::DesignRepackMover
{
public:
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::pose::Pose Pose;
public:
	BackrubDDMover();
	BackrubDDMover(
		ScoreFunctionCOP scorefxn,
		bool const backrub_partner1,
		bool const backrub_partner2,
		core::Real const interface_distance_cutoff,
		core::Size const backrub_moves,
		core::Real const mc_kt,
		core::Real const sidechain_move_prob,
		utility::vector1< core::Size > const & residues );
	void apply( Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new BackrubDDMover ); }
	virtual ~BackrubDDMover();

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:
	/// @brief Adds the given ResidueSelector to the selection. (As in all residues in both selections.)
	void add_selector( core::select::residue_selector::ResidueSelectorOP const & sele  );

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	bool backrub_partner1_;
	bool backrub_partner2_;
	core::Real interface_distance_cutoff_;
	core::Size backrub_moves_;
	core::Real mc_kt_;
	core::Real sidechain_move_prob_;
	core::Real small_move_prob_;
	core::Real bbg_move_prob_;
	core::select::residue_selector::ResidueSelectorOP residues_;
};

} // movers
} // protein_interface_design
} // protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_BackrubDDMover_HH*/
