// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/RandomizeBBByRamaPrePro.hh
/// @brief A simple mover to randomize a backbone, or a portion of a backbone, biased by the rama_prepro score of each residue.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_backbone_moves_RandomizeBBByRamaPrePro_HH
#define INCLUDED_protocols_backbone_moves_RandomizeBBByRamaPrePro_HH

// Unit headers
#include <protocols/backbone_moves/RandomizeBBByRamaPrePro.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace backbone_moves {

///@brief A simple mover to randomize a backbone, or a portion of a backbone, biased by the rama_prepro score of each residue.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class RandomizeBBByRamaPrePro : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	RandomizeBBByRamaPrePro();

	/// @brief Copy constructor.
	RandomizeBBByRamaPrePro( RandomizeBBByRamaPrePro const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~RandomizeBBByRamaPrePro() override;

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
		basic::datacache::DataMap & data
	) override;

	//RandomizeBBByRamaPrePro & operator=( RandomizeBBByRamaPrePro const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	/// @brief Get the mover name.
	std::string
	get_name() const override;

	/// @brief Get the mover name.
	static
	std::string
	mover_name();

	/// @brief Provide schema information about this mover.
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Set the residue selector that this mover will use.
	/// @details The selector is cloned.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in );

private: // methods

private: // data

	/// @brief A ResidueSelector to select a subset of residues.
	core::select::residue_selector::ResidueSelectorOP selector_;


};

std::ostream &
operator<<( std::ostream & os, RandomizeBBByRamaPrePro const & mover );

} //protocols
} //backbone_moves

#endif //protocols_backbone_moves_RandomizeBBByRamaPrePro_HH
