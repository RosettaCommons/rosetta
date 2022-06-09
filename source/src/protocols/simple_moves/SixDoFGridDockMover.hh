// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SixDoFGridDockMover.hh
/// @brief Class for enumerating docked orientations (3 translations and 3 rotations) between two chains
/// @author Odessa Goudy (oda@email.unc.edu)

#ifndef INCLUDED_protocols_simple_moves_SixDoFGridDockMover_HH
#define INCLUDED_protocols_simple_moves_SixDoFGridDockMover_HH

// Unit headers
#include <protocols/simple_moves/SixDoFGridDockMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace simple_moves {

///@brief Class for enumerating docked orientations (3 translations and 3 rotations) between two chains
class SixDoFGridDockMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	SixDoFGridDockMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	SixDoFGridDockMover( SixDoFGridDockMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~SixDoFGridDockMover() override;


public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


public:

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

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

	utility::vector1< core::Real >
	parse_range( utility::tag::TagCOP tag, std::string range_option_name );

	core::Real
	lex_position( utility::vector1< utility::vector1 < core::Real > > dof_values, core::Size index, core::Size dof);
	//utility::vector1< core::Real >
	//test_dof_1_range() const;
	// utility::vector1< core::Size >
	// res_selector() const;

	// void
	// test_dof_1_range( utility::vector1< core::Real > const & v );

	// utility::vector1< core::Size >
	// res_selector( core::pose::Pose pose, core::select::residue_selector::ResidueSelectorCOP selector);

	core::Size
	dof_sample_index() const;

	void
	dof_sample_index(core::Size setting);

private: // methods
	// void parse_range( utility::tag::TagCOP tag , std::string range_option_name ) {


private: // data

	core::select::residue_selector::ResidueSelectorCOP selector1_;
	core::select::residue_selector::ResidueSelectorCOP selector2a_;
	core::select::residue_selector::ResidueSelectorCOP selector2b_;
	int max_samples_;
	bool degree_check_;
	core::Size dof_sample_index_;
	utility::vector1< utility::vector1< core::Real > > dof_values_;
};

std::ostream &
operator<<( std::ostream & os, SixDoFGridDockMover const & mover );

} //simple_moves
} //protocols

#endif //protocols_simple_moves_SixDoFGridDockMover_HH
