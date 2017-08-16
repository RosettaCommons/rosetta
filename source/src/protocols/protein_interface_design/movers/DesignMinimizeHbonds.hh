// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/DesignMinimizeHbonds.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_DesignMinimizeHbonds_hh
#define INCLUDED_protocols_protein_interface_design_movers_DesignMinimizeHbonds_hh
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief used to design a protein to hbond preferentially to a set of target residues on the partner.
/// Hbonds involving backbone or sidechain on the target can be counted, and whether to design donors or
/// acceptors can also be defined.
class DesignMinimizeHbonds : public simple_moves::DesignRepackMover
{
public:
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::pose::Pose Pose;
public:
	DesignMinimizeHbonds();
	DesignMinimizeHbonds( ScoreFunctionCOP scorefxn_repack, ScoreFunctionCOP scorefxn_minimize,
		utility::vector1< core::Size > const & target_residues, bool const donors, bool const acceptors,
		bool const bb_hbond, bool const sc_hbond, core::Real const hbond_energy_threshold,
		core::Real interface_distance_cutoff=8.0,
		bool const repack_partner1=true, bool const repack_partner2=false, bool const repack_non_ala = true );
	DesignMinimizeHbonds( ScoreFunctionOP scorefxn_repack, ScoreFunctionOP scorefxn_minimize,
		core::Size const target_residue, bool const donors, bool const acceptors,
		bool const bb_hbond, bool const sc_hbond, core::Real const hbond_energy_threshold,
		core::Real interface_distance_cutoff=8.0,
		bool const repack_partner1=true, bool const repack_partner2=false, bool const repack_non_ala=true );
	virtual ~DesignMinimizeHbonds();
	void apply( Pose & pose ) override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new DesignMinimizeHbonds ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	bool donors_, acceptors_;
	bool bb_hbond_, sc_hbond_;
	core::Real hbond_energy_threshold_, interface_distance_cutoff_;
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_DesignMinimizeHbonds_HH*/
