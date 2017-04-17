// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/residue_selectors/HBondSelector.hh
/// @brief  A ResidueSelector that selects residues forming hydrogen bonds with specified residues. By default it ignores backbone-backbone hydrogen bonds.
/// @author Sharon Guffy (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_residue_selectors_HBondSelector_HH
#define INCLUDED_protocols_residue_selectors_HBondSelector_HH


// Unit headers
// Unit headers
#include <protocols/residue_selectors/HBondSelector.fwd.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>

//Tracer headers

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace residue_selectors {

class HBondSelector : public core::select::residue_selector::ResidueSelector {

public:

	/// @brief Constructor.
	///
	HBondSelector();
	HBondSelector( HBondSelector const & src);

	/// @brief Destructor.
	///
	~HBondSelector() override;

	/// @brief Clone function.
	/// @details Copy this object and return owning pointer to the copy (created on the heap).
	core::select::residue_selector::ResidueSelectorOP clone() const override;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const override;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) override ;

	/// @brief Get the mover class name.
	///
	std::string get_name() const override;

	/// @brief Get the mover class name.
	///
	static std::string class_name();

	/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
	///
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


public:

	// ---------- SETTERS AND GETTERS ---------------

	/// @brief Set the threshold for considering something to be a
	/// hydrogen bond.
	void set_hbond_energy_cutoff( core::Real input_value );

	/// @brief Get the threshold for considering something to be a
	/// hydrogen bond.
	core::Real get_hbond_energy_cutoff() const;

	/// @brief Set whether backbone-backbone hydrogen bonds are included
	void set_include_bb_bb ( bool const input_setting ) ;

	/// @brief Get whether backbone-backbone hydrogen bonds are included
	bool get_include_bb_bb() const ;

	/// @brief Set the scorefunction.
	/// @details Clones the input.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in) ;

	/// @brief Get the scorefunction.
	///
	core::scoring::ScoreFunctionCOP get_scorefxn() const ;

	std::string
	get_input_set_str() const;

	void
	set_input_set_str( std::string );

	core::select::residue_selector::ResidueSelectorCOP
	get_input_set_selector() const;

	void
	set_input_set_selector( core::select::residue_selector::ResidueSelectorCOP );

	//We need these two getters for the copy constructor only (they shouldn't be set externally)
	bool get_input_set_defined() const;

	bool get_use_input_set_selector() const;

private: //methods

	void compute_input_set( core::pose::Pose const &, std::set< core::Size > &) const;

private: //data


	/// @brief The energy cutoff for considering something to be a hydrogen
	//bond}.
	/// @details Defaults to -0.5.
	core::Real hbond_energy_cutoff_;

	/// @brief Should backbone-backbone hydrogen bonds be included? Default false
	bool include_bb_bb_;

	/// @brief The scorefunction to use for hydrogen bond scoring.
	/// @details If no scorefunction is provided, then the default scorefunction is used.
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief Ways of specifying for which residues we'll be searching for hydrogen bond partners
	/// @details If none are provided, all residues in the pose forming hydrogen bonds will be selected.
	core::select::residue_selector::ResidueSelectorCOP input_set_selector_;
	std::string input_set_str_;

	//We'll see if these variables are really necessary
	bool input_set_defined_;
	bool use_input_set_selector_;



#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //protocols
} //residue_selectors

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_residue_selectors_HBondSelector )
#endif // SERIALIZATION

#endif //INCLUDEDcore/select/residue_selector_HBondSelector_hh
