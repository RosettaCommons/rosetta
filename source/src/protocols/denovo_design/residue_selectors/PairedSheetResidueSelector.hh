// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.hh
/// @brief  Selects residues that are involved in strand-strand pairings
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_residue_selectors_PairedSheetResidueSelector_HH
#define INCLUDED_protocols_denovo_design_residue_selectors_PairedSheetResidueSelector_HH

// Unit headers
#include <protocols/denovo_design/residue_selectors/PairedSheetResidueSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace residue_selectors {

/// @brief Selects residues that are involved in strand-strand pairings
class PairedSheetResidueSelector : public core::select::residue_selector::ResidueSelector {
public:
	typedef core::select::residue_selector::ResidueSelectorOP ResidueSelectorOP;
	typedef core::select::residue_selector::ResidueSubset ResidueSubset;

public:
	/// @brief Constructor.
	PairedSheetResidueSelector();

	/// @brief Destructor.
	virtual
	~PairedSheetResidueSelector();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual
	ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual
	ResidueSubset apply( core::pose::Pose const & pose ) const;

	/// @brief XML parse.
	/// @details Parse RosettaScripts tags and set up this mover.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Get the mover class name.
	virtual
	std::string
	get_name() const;

	/// @brief Get the mover class name.
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public:
	// mutators

	/// @brief Sets sheet topology to be filtered
	void
	set_sheet_topology( std::string const & sheet_topology );

	/// @brief Sets secondary structure to be used in calculations
	void
	set_secstruct( std::string const & secstruct );

private:
	/// @brief Gets secondary structure string to be used in computation
	/// @param[in] pose input pose
	/// @returns String representing secondary structure (e.g. "LEEEEELLEEEEEL")
	/// @details If secstruct_ is given, that is returned.
	///          If use_dssp_ is true, DSSP secstruct is returned
	///          If use_dssp_ is false and secstruct_ is not given, return pose's secondary structure
	std::string
	get_secstruct( core::pose::Pose const & pose ) const;

	/// @brief Gets sheet topology string to be used in computation
	/// @param[in] pose input pose
	/// @returns String representing sheet topology (e.g. "1-2.A.0;2-3.P.1")
	/// @details If sheet_topology_ is given, that is returned.
	///          If sheet_topology_ is not given, and the pose has a cached StructureData,
	///          the topology string is derived from that.
	///          If a sheet topology cannot be found by the above methods, an error is
	///          thrown
	std::string
	get_sheet_topology( core::pose::Pose const & pose ) const;

private:
	std::string secstruct_;
	std::string sheet_topology_;
	bool use_dssp_;

};


} //protocols
} //denovo_design
} //residue_selectors


#endif //INCLUDEDprotocols/denovo_design/residue_selectors_PairedSheetResidueSelector_hh
