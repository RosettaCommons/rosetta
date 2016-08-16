// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/antibody/residue_selector/CDRResidueSelector.hh
/// @brief  Select CDR residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_residue_selector_CDRResidueSelector_HH
#define INCLUDED_protocols_antibody_residue_selector_CDRResidueSelector_HH

// Unit headers
#include <protocols/antibody/residue_selector/CDRResidueSelector.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>

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
namespace antibody {
namespace residue_selector {

/// @brief Select CDR residues.
class CDRResidueSelector : public core::select::residue_selector::ResidueSelector {
public:

	/// @brief Constructor.
	CDRResidueSelector();

	/// @brief Constructor giving AntibodyInfo
	CDRResidueSelector(AntibodyInfoCOP ab_info);

	/// @brief Constructor Specifying CDRs
	CDRResidueSelector( AntibodyInfoCOP ab_info, utility::vector1< CDRNameEnum > cdrs );

	/// @brief Constructor Specifying CDRs
	CDRResidueSelector( AntibodyInfoCOP ab_info, utility::vector1< bool > cdrs );

	///@brief Copy Constructor
	CDRResidueSelector( CDRResidueSelector const & src);

public:

	void
	set_cdrs( utility::vector1< bool > cdrs );

	void
	set_cdrs( utility::vector1< CDRNameEnum > cdrs );

	void
	set_ab_info(AntibodyInfoCOP ab_info);

public:

	/// @brief Destructor.
	virtual
	~CDRResidueSelector();

	/// @brief Clone operator.
	/// @details Copy the current object (creating the copy on the heap) and return an owning pointer
	/// to the copy.  All ResidueSelectors must implement this.
	virtual
	core::select::residue_selector::ResidueSelectorOP clone() const;

	/// @brief "Apply" function.
	/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
	/// indicating whether each residue is selected ("true") or not ("false").
	virtual
	core::select::residue_selector::ResidueSubset apply( core::pose::Pose const & pose ) const;

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
	static
	std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	void
	set_defaults();


private:

	AntibodyInfoCOP ab_info_;
	utility::vector1< bool > cdrs_;
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;

};


} //protocols
} //antibody
} //residue_selector


#endif //INCLUDEDprotocols/antibody/residue_selector_CDRResidueSelector_hh
