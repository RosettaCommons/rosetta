// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/DesignRestrictions
/// @brief Combine a set of ResidueSelectors and residue-level TaskOperations concisely, using boolean logic to join selectors concisely
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_task_operation_DesignRestrictions_hh
#define INCLUDED_core_pack_task_operation_DesignRestrictions_hh

// Unit headers
#include <core/pack/task/operation/DesignRestrictions.fwd.hh>

// Package headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace core {
namespace pack {
namespace task {
namespace operation {

///@brief Combine a set of ResidueSelectors and ResidueLevelTaskOperations concisely, using boolean logic to join selectors concisely
///
/// @details Basically, we want to be able to encode layer design as follows:
///   <RESIDUE_SELECTORS>
///    
///      <Layer name="surface" select_core="false" select_boundary="false" select_surface="true" use_sidechain_neighbors="true"/>
///      <Layer name="boundary" select_core="false" select_boundary="true" select_surface="false" use_sidechain_neighbors="true"/>
///      <Layer name="core" select_core="true" select_boundary="false" select_surface="false" use_sidechain_neighbors="true"/>
///      <SecondaryStructure name="helix" overlap="0" minH="3" include_terminal_loops="false" use_dssp="true" ss="H" />
///      <SecondaryStructure name="sheet" overlap="0" minE="3" include_terminal_loops="false" use_dssp="true" ss="E" />
///      <SecondaryStructure name="loop" overlap="0" minH="3" minE="3" include_terminal_loops="true" use_dssp="true" ss="L" />
///   </RESIDUE_SELECTORS>
///   <TASKOPERATIONS>
///    
///      <DesignRestrictions name="layer_design">
///         <Action selector_logic="surface AND ( helix OR sheet) "  aas="DEHKNQRST"/>
///         <Action selector_logic="surface AND loop"                aas="DEGHKNPQRST"/>
///         <Action selector_logic="boundary AND helix"              aas="ADEIKLMNQRSTVWY"/>
///         <Action selector_logic="boundary AND sheet"              aas="DEFIKLNQRSTVWY"/>
///         <Action selector_logic="boundary AND loop"               aas="ADEFGIKLMNPQRSTVWY"/>
///         <Action selector_logic="core AND helix"                  aas="AFILMVWY"/>
///         <Action selector_logic="core AND sheet"                  aas="FILVWY"/>
///         <Action selector_logic="core AND loop"                   aas="AFILMPVWY"/>
///      </DesignRestrictions>
///    
///   </TASKOPERATIONS>
/// So that it then becomes easy to adapt it to change the logic; e.g. if residues 19, 21, and 23 should keep their
/// native amino acid, and still be allowed to pack, then it will look like this
///
///   <RESIDUE_SELECTORS>
///     ... the surface, boundary, core, helix, sheet and loop selectors listed above ...
///     <ResidueIndexSelector name="special_res" resnums="19,21,23"/>
///   <RESIDUE_LEVEL_TASK_OPERATIONS>
///     <RestrictToRepackingRLT name="rtr"/>
///   <RESIDUE_LEVEL_TASK_OPERATIONS>
///   <TASKOPERATIONS>
///    
///      <DesignRestrictions name="custom_layer_design">
///         <Action residue_selector="special_res"                                    residue_level_operations="rtr"/>
///         <Action selector_logic="!special_res AND surface AND ( helix OR sheet) "  aas="DEHKNQRST"/>
///         <Action selector_logic="!special_res AND surface AND loop"                aas="DEGHKNPQRST"/>
///         <Action selector_logic="!special_res AND boundary AND helix"              aas="ADEIKLMNQRSTVWY"/>
///         <Action selector_logic="!special_res AND boundary AND sheet"              aas="DEFIKLNQRSTVWY"/>
///         <Action selector_logic="!special_res AND boundary AND loop"               aas="ADEFGIKLMNPQRSTVWY"/>
///         <Action selector_logic="!special_res AND core AND helix"                  aas="AFILMVWY"/>
///         <Action selector_logic="!special_res AND core AND sheet"                  aas="FILVWY"/>
///         <Action selector_logic="!special_res AND core AND loop"                   aas="AFILMPVWY"/>
///      </DesignRestrictions>
///    
///   </TASKOPERATIONS>
/// Where the second "Action" above is saying "apply this to residues besides the special residues that are on surface loops"
class DesignRestrictions: public core::pack::task::operation::TaskOperation {
public:

	DesignRestrictions();

	DesignRestrictions(DesignRestrictions const & src);

	~DesignRestrictions() override;

	core::pack::task::operation::TaskOperationOP
	clone() const override;

	void
	add_selector_rlto_pair(
		select::residue_selector::ResidueSelectorCOP selector,
		ResLvlTaskOperationCOP rlto );

	/// @brief Configure from a RosettaScripts XML tag.
	void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & ) override;

	//////////////////////

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;

	/// @brief Return the name used to construct this TaskOperation from an XML file
	static std::string keyname();

	/// @brief Describe the format of XML file used to initialize this TaskOperation
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

	typedef std::pair< select::residue_selector::ResidueSelectorCOP, std::list< ResLvlTaskOperationCOP > > action_pair;
	std::list< action_pair > actions_;
};

} //core
} //pack
} //task
} //operation

#endif //INCLUDED_core/pack/task/operation_DesignRestrictions_fwd_hh
