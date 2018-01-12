// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/PrepackMover.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PrepackMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PrepackMover_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/select/movemap/MoveMapFactory.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

class PrepackMover : public protocols::minimization_packing::PackRotamersMover
{
public:
	PrepackMover();
	PrepackMover( core::scoring::ScoreFunctionCOP scorefxn, core::Size jump_num );
	virtual ~PrepackMover();

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;
	void apply( core::pose::Pose & pose ) override;
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	bool min_bb() const;
	void min_bb( bool const m );
	core::kinematics::MoveMapOP mm( core::pose::Pose const & pose ) const;

	void mmf( core::select::movemap::MoveMapFactoryCOP mmf );
	void mm( core::kinematics::MoveMapOP mm );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::Size jump_num_;
	bool min_bb_; //dflt false

	core::select::movemap::MoveMapFactoryCOP mmf_;
	core::kinematics::MoveMapOP mm_; // only activated if min_bb is on.
};

} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_PrepackMover_HH*/
