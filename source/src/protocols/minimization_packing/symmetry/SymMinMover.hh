// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_hh

// Unit headers
#include <protocols/minimization_packing/symmetry/SymMinMover.fwd.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


// Package headers
//#include <protocols/moves/Mover.hh>

//#include <core/kinematics/MoveMap.fwd.hh>
//#include <core/optimization/MinimizerOptions.fwd.hh>
//#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>
//#include <core/types.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace minimization_packing {
namespace symmetry {
///////////////////////////////////////////////////////////////////////////////
class SymMinMover : public protocols::minimization_packing::MinMover
{
public:

	// default constructor
	SymMinMover();

	SymMinMover( std::string const & );

	~SymMinMover();

	// constructor with arguments
	SymMinMover(
		core::kinematics::MoveMapOP movemap_in,
		ScoreFunctionCOP scorefxn_in,
		std::string const & min_type_in,
		Real tolerance_in,
		bool use_nb_list_in,
		bool deriv_check_in = false,
		bool deriv_check_verbose_in = false
	);

	virtual void apply( core::pose::Pose & pose_ ) override;

	SymMinMover( MinMover const & src );

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		filters::Filters_map const &filters,
		moves::Movers_map const &movers,
		core::pose::Pose const & pose ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


};

} // symmetry
} // moves
} // rosetta
#endif
