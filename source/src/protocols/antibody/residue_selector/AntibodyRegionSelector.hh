// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody/residue_selector/AntibodyRegionSelector.hh
/// @brief  A simple selector to select residues of particular antibody regions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_residue_selector_AntibodyRegionSelector_HH
#define INCLUDED_protocols_antibody_residue_selector_AntibodyRegionSelector_HH

// Unit headers
#include <protocols/antibody/residue_selector/AntibodyRegionSelector.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>

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

/// @brief A simple selector to select residues of particular antibody regions.
class AntibodyRegionSelector : public core::select::residue_selector::ResidueSelector {
public:

	/// @brief Constructor.
	AntibodyRegionSelector();
	
	/// @brief Constructor Passing AntibodyInfo
	AntibodyRegionSelector( AntibodyInfoCOP ab_info );
	
	/// @brief Constructor giving the AntibodyRegion to select on.
	AntibodyRegionSelector( AntibodyInfoCOP ab_info, AntibodyRegionEnum region );
	
	/// @brief Copy Constructor
	AntibodyRegionSelector( AntibodyRegionSelector const & src);
	
public:

	void
	set_region( AntibodyRegionEnum region );
	
	void
	set_ab_info( AntibodyInfoCOP ab_info);
	
	
public:

	/// @brief Destructor.
	virtual
	~AntibodyRegionSelector();

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
	static std::string
	class_name();

	/// @brief Provide XSD information, enabling mechanical validation of input XML.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	
	void
	set_defaults();
	
private:

	AntibodyInfoCOP ab_info_;
	AntibodyRegionEnum region_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;
	
};


} //protocols
} //antibody
} //residue_selector


#endif //INCLUDEDprotocols/antibody/residue_selector_AntibodyRegionSelector_hh
