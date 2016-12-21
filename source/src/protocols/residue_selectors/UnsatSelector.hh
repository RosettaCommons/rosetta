// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/UnsatSelector.hh
/// @brief  A ResidueSelector that selects hydrogen bond acceptors or selectors that are not satisfied with h-bond
/// @author Parisa Hosseinzadeh (parisah@uw.edu)

#ifndef INCLUDED_protocols_residue_selectors_UnsatSelector_HH
#define INCLUDED_protocols_residue_selectors_UnsatSelector_HH


// Unit headers
// Unit headers
#include <protocols/residue_selectors/UnsatSelector.fwd.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

//Tracer headers
#include <basic/Tracer.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <utility/assert.hh>


namespace protocols {
namespace residue_selectors {

/// @brief A ResidueSelector that selects alpha-amino acids that are either in the positive phi or negative phi region of Ramachandran space (depending on user preferences).
class UnsatSelector : public core::select::residue_selector::ResidueSelector {

public:

	/// @brief Constructor.
	///
	UnsatSelector();

	/// @brief Destructor.
	///
	~UnsatSelector() override;

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

	/// @brief Set the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor.
	void set_mode( bool const input_setting );

	/// @brief Get the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor.
	bool mode() const ;

	/// @brief Set the threshold for considering something to be a
	/// hydrogen bond.
	void set_hbond_energy_cutoff( core::Real const &input_value );

	/// @brief Get the threshold for considering something to be a
	/// hydrogen bond.
	core::Real const & hbond_energy_cutoff() const ;

	/// @brief Set whether we only consider mainchain hydrogen bond
	/// donors and acceptors.
	void set_consider_mainchain_only ( bool const input_setting ) ;

	/// @brief Get whether we only consider mainchain hydrogen bond
	/// donors and acceptors.
	bool consider_mainchain_only() const ;

	/// @brief Set the scorefunction.
	/// @details Clones the input.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in) ;

	/// @brief Get the scorefunction.
	///
	core::scoring::ScoreFunctionCOP scorefxn() const ;

private:

	// ---------- PRIVATE FUNCTIONS -----------------

	/// @brief The function that actually calculates the value that this filter returns, called by the apply(),
	/// report(), and report_sm() functions.
	/// @details Returns the number of atoms receiving more than the allowed number of hydrogen bonds.
	utility::vector1< utility::vector1 < core::Size > > compute( core::pose::Pose const &pose ) const;

	// ---------- PRIVATE MEMBER VARIABLES ----------

	/// @brief What is the maximum allowed number of instances of an oversaturated
	/// hydrogen bond acceptor?
	/// @details Defaults to 0.

	/// @brief The energy cutoff for considering something to be a hydrogen
	//bond}.
	/// @details Defaults to -0.1.
	core::Real hbond_energy_cutoff_;

	/// @brief Should we only consider mainchain hydrogen bond donors and acceptors?
	/// @details Defaults to true.
	bool consider_mainchain_only_;
	bool acceptors_;

	/// @brief The scorefunction to use for hydrogen bond scoring.
	/// @details If no scorefunction is provided, then the default scorefunction is used.
	core::scoring::ScoreFunctionOP scorefxn_;


};


} //protocols
} //residue_selectors


#endif //INCLUDEDcore/select/residue_selector_UnsatSelector_hh
