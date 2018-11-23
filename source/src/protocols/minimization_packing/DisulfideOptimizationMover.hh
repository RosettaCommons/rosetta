// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/minimization_packing/DisulfideOptimizationMover.hh
/// @brief A Mover to jointly optimize the geometry of a pair of disulfide-bonded residues.
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_minimization_packing_DisulfideOptimizationMover_HH
#define INCLUDED_protocols_minimization_packing_DisulfideOptimizationMover_HH

// Unit headers
#include <protocols/minimization_packing/DisulfideOptimizationMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace minimization_packing {

///@brief A Mover to jointly optimize the geometry of a pair of disulfide-bonded residues.
class DisulfideOptimizationMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	DisulfideOptimizationMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	DisulfideOptimizationMover( DisulfideOptimizationMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~DisulfideOptimizationMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//DisulfideOptimizationMover & operator=( DisulfideOptimizationMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void set_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ) {
		selector_ = selector;
	}

	void set_score_function( core::scoring::ScoreFunctionCOP const & sfxn ) {
		sfxn_ = sfxn;
	}

	void set_final_optimization_n_iter( Size const final_optimization_n_iter ) {
		final_optimization_n_iter_ = final_optimization_n_iter;
	}

private: // methods

	void break_repack_reform( Pose & pose, utility::vector1< Size > const & cys_pos );

private: // data
	core::select::residue_selector::ResidueSelectorCOP selector_ = nullptr;
	core::scoring::ScoreFunctionCOP sfxn_ = nullptr;
	Size final_optimization_n_iter_ = 0;

};

std::ostream &
operator<<( std::ostream & os, DisulfideOptimizationMover const & mover );

} //protocols
} //minimization_packing

#endif //protocols_minimization_packing_DisulfideOptimizationMover_HH
