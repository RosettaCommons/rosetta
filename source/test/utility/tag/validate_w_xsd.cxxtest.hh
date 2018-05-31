// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/pose_inputters/JD2ResourceManager.cxxtest.hh
/// @brief  test suite for protocols::jd2::JD2ResourceManager and protocols::resource_manager::LazyResourceManager
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

#include <libxml/xmlschemas.h>
#include <libxml/xmlstring.h>
#include <libxml/tree.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdarg>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>

// TEMP
#include <protocols/task_operations/DsspDesignOperation.hh>
#include <protocols/task_operations/LinkResidues.hh>
#include <protocols/task_operations/RestrictToAlignedSegments.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.hh>

#include <core/select/residue_selector/ResidueSelectorFactory.hh>

#include <boost/bind.hpp>

#include <basic/Tracer.hh>

static basic::Tracer TR( "XMLSchemaValidationTests" );

using namespace utility::tag;

class XMLSchemaValidationTests : public CxxTest::TestSuite {
public:

	void setUp() {
	}

	std::string example_xsd() {
		return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:group name=\"residue_selector\">\n"
			" <xs:choice>\n"
			"  <xs:element name=\"And\" type=\"rs_AndType\"/>\n"
			"  <xs:element name=\"Chain\" type=\"rs_ChainType\"/>\n"
			"  <xs:element name=\"ClashBasedRepackShell\" type=\"rs_ClashBasedRepackShellType\"/>\n"
			"  <xs:element name=\"Index\" type=\"rs_IndexType\"/>\n"
			"  <xs:element name=\"InterfaceByVector\" type=\"rs_InterfaceByVectorType\"/>\n"
			"  <xs:element name=\"JumpDownstream\" type=\"rs_JumpDownstreamType\"/>\n"
			"  <xs:element name=\"JumpUpstream\" type=\"rs_JumpUpstreamType\"/>\n"
			"  <xs:element name=\"Layer\" type=\"rs_LayerType\"/>\n"
			"  <xs:element name=\"Neighborhood\" type=\"rs_NeighborhoodType\"/>\n"
			"  <xs:element name=\"Not\" type=\"rs_NotType\"/>\n"
			"  <xs:element name=\"NumNeighbors\" type=\"rs_NumNeighborsType\"/>\n"
			"  <xs:element name=\"Or\" type=\"rs_OrType\"/>\n"
			"  <xs:element name=\"ResiduePDBInfoHasLabel\" type=\"rs_ResiduePDBInfoHasLabelType\"/>\n"
			"  <xs:element name=\"SecondaryStructure\" type=\"rs_SecondaryStructureType\"/>\n"
			" </xs:choice>\n"
			"</xs:group>\n"
			"\n"
			"<xs:complexType name=\"rs_AndType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"2\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:simpleType name=\"chain_cslist\">\n"
			" <xs:restriction base=\"xs:string\">\n"
			"  <xs:pattern value=\"[A-Z](,[A-Z])*\"/>\n"
			" </xs:restriction>\n"
			"</xs:simpleType>\n"
			"\n"
			"<xs:complexType name=\"rs_ChainType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"chains\" type=\"chain_cslist\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_ClashBasedRepackShellType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"bump_overlap_factor\" type=\"xs:decimal\" default=\"0.5\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:simpleType name=\"int_cslist\">\n"
			" <xs:restriction base=\"xs:string\">\n"
			"  <xs:pattern value=\"[0-9]+(,[0-9]+)*\"/>\n"
			" </xs:restriction>\n"
			"</xs:simpleType>\n"
			"\n"
			"<xs:complexType name=\"rs_IndexType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"resnums\" type=\"int_cslist\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_InterfaceByVectorType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"2\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"cb_dist_cut\" type=\"xs:decimal\" default=\"11.0\"/>\n"
			" <xs:attribute name=\"nearby_atom_cut\" type=\"xs:decimal\" default=\"5.5\"/>\n"
			" <xs:attribute name=\"vector_angle_cut\" type=\"xs:decimal\" default=\"75.0\"/>\n"
			" <xs:attribute name=\"vector_dist_cut\" type=\"xs:decimal\" default=\"9.0\"/>\n"
			" <xs:attribute name=\"grp1_selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"grp2_selector\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_JumpDownstreamType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"jump\" type=\"xs:integer\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_JumpUpstreamType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"jump\" type=\"xs:integer\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_LayerType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"select_core\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"select_boundary\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"select_surface\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"use_sidechain_neighbors\" type=\"xs:boolean\" default=\"true\"/>\n"
			" <xs:attribute name=\"ball_radius\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_dist_midpoint\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_denominator\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_angle_shift_factor\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_angle_exponent\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_dist_exponent\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"core_cutoff\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"surface_cutoff\" type=\"xs:decimal\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NeighborhoodType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"resnums\" type=\"int_cslist\"/>\n"
			" <xs:attribute name=\"distance\" type=\"xs:decimal\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NotType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NumNeighborsType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"count_water\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"threshold\" type=\"xs:integer\" default=\"17\"/>\n"
			" <xs:attribute name=\"distance_cutoff\" type=\"xs:decimal\" default=\"10.0\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_OrType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"2\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_ResiduePDBInfoHasLabelType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"property\" type=\"xs:string\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_SecondaryStructureType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"overlap\" type=\"xs:integer\"/>\n"
			" <xs:attribute name=\"include_terminal_loops\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"pose_secstruct\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"ss\" type=\"xs:string\" use=\"required\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:element name=\"RESIDUE_SELECTORS\">\n"
			" <xs:complexType>\n"
			"  <xs:choice>\n"
			"   <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"unbounded\"/>\n"
			"  </xs:choice>\n"
			" </xs:complexType>\n"
			"</xs:element>\n"
			"\n"
			"<xs:element name=\"ROSETTASCRIPTS\">\n"
			" <xs:complexType>\n"
			"  <xs:sequence>\n"
			"   <xs:element ref=\"RESIDUE_SELECTORS\" minOccurs=\"0\"/>\n"
			"  </xs:sequence>\n"
			" </xs:complexType>\n"
			"</xs:element>\n"
			"\n"
			"</xs:schema>\n";
	}

	std::string example_bad_xsd() {
		return
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:group name=\"residue_selector\">\n"
			" <xs:choice>\n"
			"  <xs:element name=\"And\" type=\"rs_AndType\"/>\n"
			"  <xs:element name=\"Chain\" type=\"rs_ChainType\"/>\n"
			"  <xs:element name=\"ClashBasedRepackShell\" type=\"rs_ClashBasedRepackShellType\"/>\n"
			"  <xs:element name=\"Index\" type=\"rs_IndexType\"/>\n"
			"  <xs:element name=\"InterfaceByVector\" type=\"rs_InterfaceByVectorType\"/>\n"
			"  <xs:element name=\"JumpDownstream\" type=\"rs_JumpDownstreamType\"/>\n"
			"  <xs:element name=\"JumpUpstream\" type=\"rs_JumpUpstreamType\"/>\n"
			"  <xs:element name=\"Layer\" type=\"rs_LayerType\"/>\n"
			"  <xs:element name=\"Neighborhood\" type=\"rs_NeighborhoodType\"/>\n"
			"  <xs:element name=\"Not\" type=\"rs_NotType\"/>\n"
			"  <xs:element name=\"NumNeighbors\" type=\"rs_NumNeighborsType\"/>\n"
			"  <xs:element name=\"Or\" type=\"rs_OrType\"/>\n"
			"  <xs:element name=\"ResidueName\" type=\"rs_ResidueNameType\"/>\n"
			"  <xs:element name=\"ResiduePDBInfoHasLabel\" type=\"rs_ResiduePDBInfoHasLabelType\"/>\n"
			"  <xs:element name=\"SecondaryStructure\" type=\"rs_SecondaryStructureType\"/>\n"
			"  <xs:element name=\"Task\" type=\"rs_TaskType\"/>\n"
			" </xs:choice>\n"
			"</xs:group>\n"
			"\n"
			"<xs:complexType name=\"rs_AndType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"1\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"selectors\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:simpleType name=\"chain_cslist\">\n"
			" <xs:restriction base=\"xs:string\">\n"
			"  <xs:pattern value=\"[A-Z](,[A-Z])*\"/>\n"
			" </xs:restriction>\n"
			"</xs:simpleType>\n"
			"\n"
			"<xs:complexType name=\"rs_ChainType\" mixed=\"true\">\n"
			" <xs:attribute name=\"chains\" type=\"chain_cslist\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_ClashBasedRepackShellType\" mixed=\"true\">\n"
			" <xs:attribute name=\"bump_overlap_factor\" type=\"real\" default=\"0.5\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:simpleType name=\"int_cslist\">\n"
			" <xs:restriction base=\"xs:string\">\n"
			"  <xs:pattern value=\"[0-9]+(,[0-9]+)*\"/>\n"
			" </xs:restriction>\n"
			"</xs:simpleType>\n"
			"\n"
			"<xs:complexType name=\"rs_IndexType\" mixed=\"true\">\n"
			" <xs:attribute name=\"resnums\" type=\"int_cslist\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_InterfaceByVectorType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"cb_dist_cut\" type=\"xs:decimal\" default=\"11.0\"/>\n"
			" <xs:attribute name=\"nearby_atom_cut\" type=\"xs:decimal\" default=\"5.5\"/>\n"
			" <xs:attribute name=\"vector_angle_cut\" type=\"xs:decimal\" default=\"75.0\"/>\n"
			" <xs:attribute name=\"vector_dist_cut\" type=\"xs:decimal\" default=\"9.0\"/>\n"
			" <xs:attribute name=\"grp1_selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"grp2_selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_JumpDownstreamType\" mixed=\"true\">\n"
			" <xs:attribute name=\"jump\" type=\"xs:integer\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_JumpUpstreamType\" mixed=\"true\">\n"
			" <xs:attribute name=\"jump\" type=\"xs:integer\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_LayerType\" mixed=\"true\">\n"
			" <xs:attribute name=\"select_core\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"select_boundary\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"select_surface\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"cache_selection\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"use_sidechain_neighbors\" type=\"xs:boolean\" default=\"true\"/>\n"
			" <xs:attribute name=\"ball_radius\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_dist_midpoint\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_denominator\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_angle_shift_factor\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_angle_exponent\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"sc_neighbor_dist_exponent\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"core_cutoff\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"surface_cutoff\" type=\"xs:decimal\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NeighborhoodType\" mixed=\"true\">\n"
			" <xs:all>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"1\"/>\n"
			" </xs:all>\n"
			" <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"resnums\" type=\"int_cslist\"/>\n"
			" <xs:attribute name=\"distance\" type=\"xs:decimal\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NotType\" mixed=\"true\">\n"
			" <xs:all>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"1\"/>\n"
			" </xs:all>\n"
			" <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_NumNeighborsType\" mixed=\"true\">\n"
			" <xs:attribute name=\"count_water\" type=\"xs:boolean\" default=\"false\"/>\n"
			" <xs:attribute name=\"threshold\" type=\"xs:integer\" default=\"17\"/>\n"
			" <xs:attribute name=\"distance_cutoff\" type=\"xs:decimal\" default=\"10.0\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_OrType\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"1\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"selectors\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_ResidueNameType\" mixed=\"true\">\n"
			" <xs:attribute name=\"residue_name3\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"residue_names\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_ResiduePDBInfoHasLabelType\" mixed=\"true\">\n"
			" <xs:attribute name=\"property\" type=\"xs:string\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_SecondaryStructureType\" mixed=\"true\">\n"
			" <xs:attribute name=\"overlap\" type=\"xs:integer\"/>\n"
			" <xs:attribute name=\"include_terminal_loops\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"pose_secstruct\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"ss\" type=\"xs:string\" use=\"required\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"rs_TaskType\" mixed=\"true\">\n"
			" <xs:attribute name=\"task_operations\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"designable\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"packable\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"repackable\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"fixed\" type=\"xs:boolean\"/>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"RESIDUE_SELECTORS_Type\" mixed=\"true\">\n"
			" <xs:choice>\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"unbounded\"/>\n"
			" </xs:choice>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:element name=\"ROSETTASCRIPTS\">\n"
			" <xs:complexType mixed=\"true\">\n"
			"  <xs:sequence>\n"
			"   <xs:element name=\"RESIDUE_SELECTORS\" type=\"RESIDUE_SELECTORS_Type\" minOccurs=\"0\"/>\n"
			"  </xs:sequence>\n"
			" </xs:complexType>\n"
			"</xs:element>\n"
			"\n"
			"</xs:schema>\n";
	}


	std::string example_xml_doc1() {
		return
			"<ROSETTASCRIPTS>\n"
			" <RESIDUE_SELECTORS>\n"
			"  <And name=\"and_selector\">\n"
			"   <Chain chains=\"A,B\"/>\n"
			"   <Neighborhood resnums=\"1,2,3\" distance=\"10.0\"/>\n"
			"  </And>\n"
			"  <JumpDownstream name=\"ds_jump_2\" jump=\"2\"/>\n"
			" </RESIDUE_SELECTORS>\n"
			"</ROSETTASCRIPTS>\n";
	}

	// Here's an xml file that does not match the schema, and so is invalid
	std::string example_xml_doc2() {
		return
			"<ROSETTASCRIPTS>\n"
			" <RESIDUE_SELECTORS>\n"
			"  <And name=\"and_selector\">\n"
			"   <Chain chains=\"A,B\"/>\n"
			"   <Neighborhood resnums=\"1,2,3,\" distance=\"10.0\"/>\n"
			"  </And>\n"
			"  <JumpDownstream name=\"ds_jump_2\" jump=\"2\"/>\n"
			" </RESIDUE_SELECTORS>\n"
			"</ROSETTASCRIPTS>\n";
	}

	// here's an xml file that is missing quotes and so it is not properly
	// formatted -- the error message we get from having a bad XML file
	// should be different (but definitely present!) from the error message
	// of having an invalid xml file (i.e. one that does not match the schema
	std::string example_xml_doc3() {
		return
			"<ROSETTASCRIPTS>\n"
			" <RESIDUE_SELECTORS>\n"
			"  <And name=and_selector>\n"
			"   <Chain chains=A,B/>\n"
			"   <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"  </And>\n"
			"  <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" </RESIDUE_SELECTORS>\n"
			"</ROSETTASCRIPTS>\n";
	}

	void dont_test_residue_selector_schemas() {
		TR << "XSD for ResidueSelectors" << std::endl;
		XMLSchemaDefinition xsd;
		core::select::residue_selector::ResidueSelectorFactory::get_instance()
			->define_residue_selector_xml_schema( xsd );

		XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
		XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );

		XMLSchemaElementOP residue_selectors_element( new XMLSchemaElement );
		//XMLSchemaModelGroupOP residue_selectors_choice( new XMLSchemaModelGroup( xsmgt_choice ));
		XMLSchemaComplexTypeOP residue_selectors_type( new XMLSchemaComplexType );

		XMLSchemaModelGroupOP group_ref( new XMLSchemaModelGroup( xsmgt_group ) );
		group_ref->group_name( "residue_selector" );
		group_ref->min_occurs( 0 );
		group_ref->max_occurs( xsminmax_unbounded );
		residue_selectors_type->name( "RESIDUE_SELECTORS_Type" );
		residue_selectors_type->set_model_group( group_ref );

		xsd.add_top_level_element( *residue_selectors_type );

		residue_selectors_element->name( "RESIDUE_SELECTORS" );
		residue_selectors_element->type_name( "RESIDUE_SELECTORS_Type" );
		residue_selectors_element->min_occurs( 0 );

		rosetta_scripts_seq->append_particle( residue_selectors_element );
		rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

		rosetta_scripts_element->name( "ROSETTASCRIPTS" );
		rosetta_scripts_element->element_type_def( rosetta_scripts_type );

		xsd.add_top_level_element( *rosetta_scripts_element );

		//TR << xsd.full_definition() << std::endl;

		XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
		TS_ASSERT( ! output.valid() );
		TR << "output.warnings:\n" << output.warning_messages() << std::endl;
		TR << "output.errors:\n" << output.error_messages() << std::endl;
	}

	void dont_test_task_op_schemas() {
		TR << "XSD for DsspDesignOperation: " << std::endl;
		XMLSchemaDefinition xsd;
		protocols::task_operations::DsspDesignOperation::provide_xml_schema( xsd );
		TR << xsd.full_definition() << std::endl;

		TR << "XSD for LinkResidues" << std::endl;
		XMLSchemaDefinition xsd2;
		protocols::task_operations::LinkResidues::provide_xml_schema( xsd2 );
		TR << xsd2.full_definition() << std::endl;

		TR << "XSD for RestrictToAlignedSegments" << std::endl;
		XMLSchemaDefinition xsd3;
		protocols::task_operations::RestrictToAlignedSegmentsOperation::provide_xml_schema( xsd3 );
		TR << xsd3.full_definition() << std::endl;

		TR << "XSD for SeqprofConsensusOperation" << std::endl;
		XMLSchemaDefinition xsd4;
		protocols::task_operations::SeqprofConsensusOperation::provide_xml_schema( xsd4 );
		TR << xsd4.full_definition() << std::endl;

		TS_ASSERT( xsd.full_definition().size() > 0 );
	}

	void test_validate_good_xml_file_against_xsd_using_libxml2() {
		XMLValidationOutput output = validate_xml_against_xsd( example_xml_doc1(), example_xsd() );
		TS_ASSERT( output.valid() );
		TS_ASSERT( output.errors().empty() );
		TS_ASSERT( output.warnings().empty() );

		////std::string xsd_fname( argv[1] );
		////std::cerr << "loading XML Schema from " << xsd_fname << std::endl;
		////std::ifstream xsd_input( xsd_fname.c_str() );
		////xsd_string = slurp( xsd_input );
		//std::string xsd_string = example_xsd();
		//
		//xmlChar * xsd_xmlchar_string = xmlCharStrdup( xsd_string.c_str() );
		////xmlChar * xml_version = xmlCharStrdup( "1.0" );
		//xmlDoc * xsd_doc = xmlParseDoc( xsd_xmlchar_string );
		//
		//xmlSchemaParserCtxtPtr schema_parser_context = xmlSchemaNewDocParserCtxt( xsd_doc );
		//xmlSchemaPtr schema = xmlSchemaParse( schema_parser_context );
		////TR << "Parsed the schema" << std::endl;
		//xmlSchemaValidCtxtPtr schema_validator = xmlSchemaNewValidCtxt( schema );
		////TR << "Created schema validator" << std::endl;
		//
		//XMLErrorHandler handler;
		//xmlSchemaSetValidErrors( schema_validator, handle_xml_error, handle_xml_warning, &handler );
		//
		//// std::string xml_fname( argv[2] );
		//// std::cerr << "loading XML file from " << xml_fname << std::endl;
		//// std::ifstream xml_input_stream( xml_fname.c_str() );
		//std::string xml_input_string = example_xml_doc1();
		//xmlChar * xml_input_xmlchar = xmlCharStrdup( xml_input_string.c_str() );
		//xmlDoc * xml_doc = xmlParseDoc( xml_input_xmlchar );
		//
		////TR << "Validating XML document" << std::endl;
		//int validation_output = xmlSchemaValidateDoc( schema_validator, xml_doc );
		//
		////TR << "Validation output: " << validation_output << std::endl;
		//
		//TS_ASSERT_EQUALS( validation_output, 0 );
		//
		//// clean up
		//free( xsd_xmlchar_string );
		//free( xml_input_xmlchar );
		//xmlSchemaFree( schema );
		//xmlSchemaFreeParserCtxt( schema_parser_context );
		//xmlSchemaFreeValidCtxt( schema_validator );
		//xmlFreeDoc( xsd_doc );

	}

	void test_validate_invalid_xml_file_against_xsd_using_libxml2() {
		XMLValidationOutput output = validate_xml_against_xsd( example_xml_doc2(), example_xsd() );
		TS_ASSERT( ! output.valid() );
		TS_ASSERT_EQUALS( output.errors().size(), 2 );
		TS_ASSERT( output.warnings().empty() );

		std::string gold_errors =
			"From line 5:\n"
			"Error: Element 'Neighborhood', attribute 'resnums': [facet 'pattern'] The value '1,2,3,' is not accepted by the pattern '[0-9]+(,[0-9]+)*'.\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=\"and_selector\">\n"
			" 4:    <Chain chains=\"A,B\"/>\n"
			" 5:    <Neighborhood resnums=\"1,2,3,\" distance=\"10.0\"/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=\"ds_jump_2\" jump=\"2\"/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"From line 5:\n"
			"Error: Element 'Neighborhood', attribute 'resnums': '1,2,3,' is not a valid value of the atomic type 'int_cslist'.\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=\"and_selector\">\n"
			" 4:    <Chain chains=\"A,B\"/>\n"
			" 5:    <Neighborhood resnums=\"1,2,3,\" distance=\"10.0\"/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=\"ds_jump_2\" jump=\"2\"/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n";
		TS_ASSERT_EQUALS( output.error_messages(), gold_errors );

		//TR << "Errors:\n" << output.error_messages() << std::endl;

	}

	void test_validate_bad_xml_file_against_xsd_using_libxml2() {
		XMLValidationOutput output = validate_xml_against_xsd( example_xml_doc3(), example_xsd() );
		TS_ASSERT( ! output.valid() );
		TS_ASSERT_EQUALS( output.errors().size(), 15 );
		TS_ASSERT( output.warnings().empty() );

		std::string gold_errors =
			"Error: AttValue: \" or \' expected\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"Error: attributes construct error\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"Error: Couldn\'t find end of Start Tag And line 3\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"Error: AttValue: \" or \' expected\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"9: </ROSETTASCRIPTS>\n"
			"Error: attributes construct error\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"9: </ROSETTASCRIPTS>\n"
			"Error: Couldn\'t find end of Start Tag Chain line 4\n"
			"\n"
			"1: <ROSETTASCRIPTS>\n"
			"2:  <RESIDUE_SELECTORS>\n"
			"3:   <And name=and_selector>\n"
			"4:    <Chain chains=A,B/>\n"
			"5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			"6:   </And>\n"
			"7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			"8:  </RESIDUE_SELECTORS>\n"
			"9: </ROSETTASCRIPTS>\n"
			"Error: AttValue: \" or \' expected\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: attributes construct error\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: Couldn\'t find end of Start Tag Neighborhood line 5\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: Opening and ending tag mismatch: RESIDUE_SELECTORS line 2 and And\n"
			"\n"
			" 1: <ROSETTASCRIPTS>\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: AttValue: \" or \' expected\n"
			"\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: attributes construct error\n"
			"\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: Couldn\'t find end of Start Tag JumpDownstream line 7\n"
			"\n"
			" 2:  <RESIDUE_SELECTORS>\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: Opening and ending tag mismatch: ROSETTASCRIPTS line 1 and RESIDUE_SELECTORS\n"
			"\n"
			" 3:   <And name=and_selector>\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n"
			"Error: Extra content at the end of the document\n"
			"\n"
			" 4:    <Chain chains=A,B/>\n"
			" 5:    <Neighborhood resnums=1,2,3 distance=10.0/>\n"
			" 6:   </And>\n"
			" 7:   <JumpDownstream name=ds_jump_2 jump=2/>\n"
			" 8:  </RESIDUE_SELECTORS>\n"
			" 9: </ROSETTASCRIPTS>\n"
			"10: \n";

		TS_ASSERT_EQUALS( output.error_messages(), gold_errors );

		//TR << "Errors:\n" << output.error_messages() << std::endl;

	}

	void test_validate_xml_file_against_bad_xsd_using_libxml2() {
		XMLValidationOutput output = validate_xml_against_xsd( example_xml_doc1(), example_bad_xsd() );
		TS_ASSERT( ! output.valid() );
		TS_ASSERT_EQUALS( output.errors().size(), 3 );
		TS_ASSERT( output.warnings().empty() );

		std::string gold_errors =
			"Your XML Schema failed to self-validate.  This is essentially a compile-time error, but it cannot be"
			" exposed until runtime.  You've edited some class’s provide_XML_schema function recently in such a way that it is"
			" now providing an invalid or incomplete schema.  The most likely problem is that you put <, >, or & in an attribute's"
			" description.  Email the devel list if the error message below does not help you locate and fix the problem.\n"
			"From line 102:\n"
			"Error: Element \'{http://www.w3.org/2001/XMLSchema}all\': The content is not valid. Expected is (annotation?, (annotation?, element*).\n"
			"\n"
			" 97:  <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" 98: </xs:complexType>\n"
			" 99: \n"
			"100: <xs:complexType name=\"rs_NeighborhoodType\" mixed=\"true\">\n"
			"101:  <xs:all>\n"
			"102:   <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"1\"/>\n"
			"103:  </xs:all>\n"
			"104:  <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			"105:  <xs:attribute name=\"resnums\" type=\"int_cslist\"/>\n"
			"106:  <xs:attribute name=\"distance\" type=\"xs:decimal\" use=\"required\"/>\n"
			"107:  <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"From line 112:\n"
			"Error: Element \'{http://www.w3.org/2001/XMLSchema}all\': The content is not valid. Expected is (annotation?, (annotation?, element*).\n"
			"\n"
			"107:  <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"108: </xs:complexType>\n"
			"109: \n"
			"110: <xs:complexType name=\"rs_NotType\" mixed=\"true\">\n"
			"111:  <xs:all>\n"
			"112:   <xs:group ref=\"residue_selector\" minOccurs=\"0\" maxOccurs=\"1\"/>\n"
			"113:  </xs:all>\n"
			"114:  <xs:attribute name=\"selector\" type=\"xs:string\"/>\n"
			"115:  <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			"116: </xs:complexType>\n"
			"117: \n";

		TS_ASSERT_EQUALS( output.error_messages(), gold_errors );

		//TR << "Errors:\n" << output.error_messages() << std::endl;

	}


	void dont_test_validate_bad_xml_file_from_disk() {
		//std::string xsd_fname( argv[1] );
		//std::cerr << "loading XML Schema from " << xsd_fname << std::endl;
		//std::ifstream xsd_input( xsd_fname.c_str() );
		//xsd_string = slurp( xsd_input );
		std::string xsd_string = example_xsd();

		xmlChar * xsd_xmlchar_string = xmlCharStrdup( xsd_string.c_str() );
		//xmlChar * xml_version = xmlCharStrdup( "1.0" );
		xmlDoc * xsd_doc = xmlParseDoc( xsd_xmlchar_string );

		xmlSchemaParserCtxtPtr schema_parser_context = xmlSchemaNewDocParserCtxt( xsd_doc );
		xmlSchemaPtr schema = xmlSchemaParse( schema_parser_context );
		//TR << "Parsed the schema" << std::endl;
		xmlSchemaValidCtxtPtr schema_validator = xmlSchemaNewValidCtxt( schema );
		//TR << "Created schema validator" << std::endl;

		XMLErrorHandler handler;
		xmlSetStructuredErrorFunc( & handler, handle_structured_xml_error );

		//xmlSchemaSetValidErrors( schema_validator, handle_xml_error, handle_xml_warning, &handler );

		// std::string xml_fname( argv[2] );
		// std::cerr << "loading XML file from " << xml_fname << std::endl;
		// std::ifstream xml_input_stream( xml_fname.c_str() );
		//std::string xml_input_string = example_xml_doc2();
		//xmlChar * xml_input_xmlchar = xmlCharStrdup( xml_input_string.c_str() );
		//xmlDoc * xml_doc = xmlParseDoc( xml_input_xmlchar );

		//TR << "Validating XML document" << std::endl;
		//int validation_output = xmlSchemaValidateDoc( schema_validator, xml_doc );
		int validation_output = xmlSchemaValidateFile( schema_validator, "protocols/jd3/example_bad.xml", 0 );

		//TR << "Validation output: " << validation_output << std::endl;

		TS_ASSERT_EQUALS( validation_output, 1824 );
		TS_ASSERT_EQUALS( handler.errors().size(), 2 );
		if ( handler.errors().size() != 2 ) return;

		utility::vector1< std::string > verrors;
		for ( XMLErrorHandler::Errors::const_iterator iter = handler.errors().begin(), iter_end = handler.errors().end();
				iter != iter_end; ++iter ) {
			verrors.push_back( *iter );
		}
		TS_ASSERT_EQUALS( verrors[1], std::string("Element 'Neighborhood', attribute 'resnums': [facet 'pattern'] The value '1,2,3,' is not accepted by the pattern '[0-9]+(,[0-9]+)*'.\n"));
		TS_ASSERT_EQUALS( verrors[2], std::string("Element 'Neighborhood', attribute 'resnums': '1,2,3,' is not a valid value of the atomic type 'int_cslist'.\n" ));

		// clean up
		free( xsd_xmlchar_string );
		//free( xml_input_xmlchar );
		xmlSchemaFree( schema );
		xmlSchemaFreeParserCtxt( schema_parser_context );
		xmlSchemaFreeValidCtxt( schema_validator );
		//xmlFreeDoc( xsd_doc );


	}


	void test_XMLValidator_validate_one_file() {
		using namespace utility::tag;
		XMLValidator validator;
		{
			XMLValidationOutput output = validator.set_schema( example_xsd() );
			TS_ASSERT( output.valid() );
		}

		{
			XMLValidationOutput output = validator.validate_xml_against_schema( example_xml_doc1() );
			TS_ASSERT( output.valid() );
		}
	}

	void test_XMLValidator_validate_several_files() {
		using namespace utility::tag;
		XMLValidator validator;
		{
			XMLValidationOutput output = validator.set_schema( example_xsd() );
			TS_ASSERT( output.valid() );
		}

		{
			XMLValidationOutput output = validator.validate_xml_against_schema( example_xml_doc1() );
			TS_ASSERT( output.valid() );
		}

		{
			XMLValidationOutput output = validator.validate_xml_against_schema( example_xml_doc2() );
			TS_ASSERT( ! output.valid() );
		}

		{
			XMLValidationOutput output = validator.validate_xml_against_schema( example_xml_doc3() );
			TS_ASSERT( ! output.valid() );
		}

		{
			XMLValidationOutput output = validator.validate_xml_against_schema( example_xml_doc1() );
			TS_ASSERT( output.valid() );
		}
	}

};
