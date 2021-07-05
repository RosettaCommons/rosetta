// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/metal_interface/RemoveMetalConnectionsMover.hh
/// @brief A mover that removes the connections to metals that were added by the SetupMetalsMover
/// or by the -auto_setup_metals flag.
/// @details This mover:
///     - Removes the bonds between metals and metal-binding residues.
///     - Reverts metal-liganding residues back to their pre-bonded types.
///     - Reverts metals back to their pre-bonded types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMover_HH
#define INCLUDED_protocols_metal_interface_RemoveMetalConnectionsMover_HH

// Unit headers
#include <protocols/metal_interface/RemoveMetalConnectionsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace metal_interface {

/// @brief A mover that removes the connections to metals that were added by the SetupMetalsMover
/// or by the -auto_setup_metals flag.
/// @details This mover:
///     - Removes the bonds between metals and metal-binding residues.
///     - Reverts metal-liganding residues back to their pre-bonded types.
///     - Reverts metals back to their pre-bonded types.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class RemoveMetalConnectionsMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	RemoveMetalConnectionsMover();

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~RemoveMetalConnectionsMover() override;

public:

	/////////////////////
	/// Mover Methods ///
	/////////////////////


	/// @brief Apply the mover.
	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Summarize this mover's configuration.
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

	//RemoveMetalConnectionsMover & operator=( RemoveMetalConnectionsMover const & src );

	/// @brief required in the context of the parser/scripting scheme.
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme.
	protocols::moves::MoverOP
	clone() const override;

	/// @brief Returns "RemoveMetalConnectionsMover".
	std::string
	get_name() const override;

	/// @brief Returns "RemoveMetalConnectionsMover".  Useful in situations in which
	/// an instance is not available.
	static
	std::string
	mover_name();

	/// @brief Provide a user- or machine-readable description of this mover's XML interface.
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Sets the residue selector for selecting residues from which metal bonds
	/// will be removed.
	/// @details If a polymer residue AND a metal that it binds are selected, the corresponding
	/// metal bonds and variants are removed.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP const & selector_in );

	/// @brief Gets the residue selector for selecting residues from which metal bonds
	/// will be removed.
	/// @details If a polymer residue AND a metal that it binds are selected, the corresponding
	/// metal bonds and variants are removed.
	core::select::residue_selector::ResidueSelectorCOP residue_selector() const;

public: //Function overrides needed for the citation manager:

	/// @brief Provide citations to the passed CitationCollectionList.
	/// @details This mover is unpublished.  It returns Vikram K. Mulligan as its author.
	void provide_citation_info(basic::citation_manager::CitationCollectionList & citations) const override;

private: // methods

	/// @brief Remove one or two variant types from a pose residue.
	/// @details Based on name lookup, so hopefully should work for on-the-fly variants.
	void
	remove_variant_types_from_res(
		core::pose::Pose & pose,
		core::Size const resindex,
		std::string const & this_atom,
		bool const this_is_metal
	) const;

private: // data

	/// @brief Residue selector for residues to remove bonds from.
	/// @details If a polymer residue AND a metal that it binds are selected, the corresponding
	/// metal bonds and variants are removed.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

};

std::ostream &
operator<<( std::ostream & os, RemoveMetalConnectionsMover const & mover );

} //metal_interface
} //protocols

#endif //protocols_metal_interface_RemoveMetalConnectionsMover_HH
