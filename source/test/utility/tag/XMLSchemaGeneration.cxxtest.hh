// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/utility/tag/XMLSchemaGeneration.cxxtest.hh
/// @brief  test suite for the XMLSchemaGeneration classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <utility/tag/XMLSchemaGeneration.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/bind.hpp>

// C++ headers
#include <sstream>

using namespace utility::tag;

/// @brief Little utility function for identifying the first position at which two strings
/// differ.  Returns 1 greater than the shorter of the two strings if there is no difference
/// until the first string ends -- compare what is returned against the max of the two strings'
/// sizes() to tell if both strings are equal.
platform::Size
first_difference( std::string const & s1, std::string const & s2 ) {
	for ( platform::Size ii = 0; ii < std::min( s1.size(), s2.size() ); ++ii ) {
		if ( s1[ii] != s2[ii] ) { return ii; }
	}
	return std::min( s1.size(), s2.size() );
}

/// @brief Dummy class used to test the XMLSchemaDefinition class
class ResSelectCT : public XMLSchemaTopLevelElement {
public:
	ResSelectCT() : name_( "ResidueSelectorType" ) {}
	virtual std::string const & element_name() const { return name_; }
	virtual void write_definition( int indentation, std::ostream & os ) const {
		for ( int ii = 1; ii <= indentation; ++ii ) os << " ";
		os <<
			"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
			" <xs:group ref=\"ResidueSelector\"/>\n"
			"</xs:complexType>\n";
	}
	virtual void prepare_for_output( XMLSchemaDefinition & ) const {}

private:
	std::string name_;
};

/// @brief Dummy class for testing the XMLSchemaDefinition class
class ResSelectGrp : public XMLSchemaTopLevelElement {
public:
	ResSelectGrp() : name_( "ResidueSelector" ) {}
	virtual std::string const & element_name() const { return name_; }
	virtual void write_definition( int indentation, std::ostream & os ) const {
		for ( int ii = 1; ii <= indentation; ++ii ) os << " ";
		os << "<xs:group name=\"ResidueSelector\">\n"
			" <xs:choice>\n"
			"  <xs:element name=\"AndResidueSelector\" type=\"AndResidueSelectorType\"/>\n"
			"  <xs:element name=\"OrResidueSelector\" type=\"OrResidueSelectorType\"/>\n"
			"  <xs:element name=\"NotResidueSelector\" type=\"NotResidueSelectorType\"/>\n"
			"  <xs:element name=\"ChainResidueSelector\" type=\"ChainResidueSelectorType\"/>\n"
			"  <xs:element name=\"NeighborhoodResidueSelector\" type=\"NeighborhoodResidueSelectorType\"/>\n"
			" </xs:choice>\n"
			"</xs:group>\n";
	}
	virtual void prepare_for_output( XMLSchemaDefinition & ) const {}

private:
	std::string name_;
};

class XMLSchemaGenerationTests : public CxxTest::TestSuite
{
public:

	void compare_line_by_line( std::string const & gold, std::string const & actual )
	{
		std::vector< std::string > gold_lines = utility::split_by_newlines( gold );
		std::vector< std::string > xsd_lines  = utility::split_by_newlines( actual );
		TS_ASSERT_EQUALS( gold_lines.size(), xsd_lines.size() );
		for ( platform::Size ii = 0; ii < std::min( gold_lines.size(), xsd_lines.size()); ++ii ) {
			TS_ASSERT_EQUALS( gold_lines[ ii ], xsd_lines[ ii ] );
		}
	}

	void test_xmlschema_data_type() {
		XMLSchemaType xsstr; xsstr.type( xs_string );
		TS_ASSERT_EQUALS( std::string( "xs:string" ), xsstr.type_name() );

		XMLSchemaType xsdec; xsdec.type( xs_decimal );
		TS_ASSERT_EQUALS( std::string( "xs:decimal" ), xsdec.type_name() );

		XMLSchemaType xsint; xsint.type( xs_integer );
		TS_ASSERT_EQUALS( std::string( "xs:integer" ), xsint.type_name() );

		XMLSchemaType xsbool; xsbool.type( xs_boolean );
		TS_ASSERT_EQUALS( std::string( "xs:boolean" ), xsbool.type_name() );

		XMLSchemaType xsdate; xsdate.type( xs_date );
		TS_ASSERT_EQUALS( std::string( "xs:date" ), xsdate.type_name() );

		XMLSchemaType xstime; xstime.type( xs_time );
		TS_ASSERT_EQUALS( std::string( "xs:time" ), xstime.type_name() );

		XMLSchemaType cs_intlist;
		cs_intlist.custom_type_name( "int_cslist" );
		TS_ASSERT_EQUALS( std::string( "int_cslist" ), cs_intlist.type_name() );
	}

	void test_xs_attribute_custom() {
		XMLSchemaAttribute att;
		att.name( "residues" );
		XMLSchemaType cs_intlist;
		cs_intlist.custom_type_name( "int_cslist" );
		att.type( cs_intlist );
		std::ostringstream oss;
		att.write_definition( 0, oss );
		TS_ASSERT_EQUALS( std::string( "<xs:attribute name=\"residues\" type=\"int_cslist\"/>\n" ), oss.str() );
	}

	void test_xs_attribute_string() {
		XMLSchemaAttribute att;
		att.name( "chain" );
		XMLSchemaType xsstr( xs_string );
		att.type( xsstr );
		std::ostringstream oss;
		att.write_definition( 0, oss );
		TS_ASSERT_EQUALS( std::string( "<xs:attribute name=\"chain\" type=\"xs:string\"/>\n" ), oss.str() );
	}

	void test_xs_attribute_string_required() {
		XMLSchemaAttribute att;
		att.name( "chain" );
		XMLSchemaType xsstr( xs_string );
		att.type( xsstr );
		att.is_required( true );
		std::ostringstream oss;
		att.write_definition( 0, oss );
		TS_ASSERT_EQUALS( std::string( "<xs:attribute name=\"chain\" type=\"xs:string\" use=\"required\"/>\n" ), oss.str() );
	}

	void test_xs_attribute_int_default() {
		XMLSchemaAttribute att;
		att.name( "jump" );
		att.type( XMLSchemaType( xs_integer ) );
		att.default_value( "1" );
		std::ostringstream oss;
		att.write_definition( 0, oss );
		TS_ASSERT_EQUALS( std::string( "<xs:attribute name=\"jump\" type=\"xs:integer\" default=\"1\"/>\n" ), oss.str() );
	}

	void test_xmlschema_restriction_int_range() {
		XMLSchemaRestriction int_1_to_100;
		int_1_to_100.name( "grade" );
		int_1_to_100.base_type( XMLSchemaType( xs_integer ) );
		int_1_to_100.add_restriction( xsr_minInclusive, "1" );
		int_1_to_100.add_restriction( xsr_maxInclusive, "100" );
		std::string definition = "<xs:simpleType name=\"grade\">\n <xs:restriction base=\"xs:integer\">\n  <xs:minInclusive value=\"1\"/>\n"
			"  <xs:maxInclusive value=\"100\"/>\n </xs:restriction>\n</xs:simpleType>\n";
		std::ostringstream oss;
		int_1_to_100.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );

		std::string definition_indent2 = "  <xs:simpleType name=\"grade\">\n   <xs:restriction base=\"xs:integer\">\n    <xs:minInclusive value=\"1\"/>\n"
			"    <xs:maxInclusive value=\"100\"/>\n   </xs:restriction>\n  </xs:simpleType>\n";
		std::ostringstream oss2;
		int_1_to_100.write_definition( 2, oss2 );
		TS_ASSERT_EQUALS( definition_indent2, oss2.str() );

	}

	void test_xmlschema_restriction_enumeration() {
		XMLSchemaRestriction enumeration;
		enumeration.name( "color_channel" );
		enumeration.base_type( XMLSchemaType( xs_string ) );
		enumeration.add_restriction( xsr_enumeration, "red" );
		enumeration.add_restriction( xsr_enumeration, "green" );
		enumeration.add_restriction( xsr_enumeration, "blue" );
		std::string definition = "<xs:simpleType name=\"color_channel\">\n <xs:restriction base=\"xs:string\">\n  <xs:enumeration value=\"red\"/>\n"
			"  <xs:enumeration value=\"green\"/>\n  <xs:enumeration value=\"blue\"/>\n </xs:restriction>\n</xs:simpleType>\n";
		std::ostringstream oss;
		enumeration.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xmlschema_restriction_pattern() {
		XMLSchemaRestriction comma_separated_int_list;
		comma_separated_int_list.name( "int_cslist" );
		comma_separated_int_list.base_type( XMLSchemaType( xs_string ) );
		comma_separated_int_list.add_restriction( xsr_pattern, "[0-9]+(,[0-9]+)*" );
		std::string definition = "<xs:simpleType name=\"int_cslist\">\n <xs:restriction base=\"xs:string\">\n  <xs:pattern value=\"[0-9]+(,[0-9]+)*\"/>\n"
			" </xs:restriction>\n</xs:simpleType>\n";
		std::ostringstream oss;
		comma_separated_int_list.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xmlschema_restriction_type_names() {
		TS_ASSERT_EQUALS( restriction_type_name( xsr_enumeration ), "xs:enumeration" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_fractionDigits ), "xs:fractionDigits" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_length ), "xs:length" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_maxExclusive ), "xs:maxExclusive" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_maxInclusive ), "xs:maxInclusive" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_maxLength ), "xs:maxLength" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_minExclusive ), "xs:minExclusive" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_minInclusive ), "xs:minInclusive" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_minLength ), "xs:minLength" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_pattern ), "xs:pattern" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_totalDigits ), "xs:totalDigits" );
		TS_ASSERT_EQUALS( restriction_type_name( xsr_whitespace ), "xs:whitespace" );
	}

	void test_xml_complex_type_no_children_no_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "to_RestrictToRepackingType" );
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( oss.str(), "<xs:complexType name=\"to_RestrictToRepackingType\" mixed=\"true\">\n</xs:complexType>\n" );
	}

	void test_xml_complex_type_sequence1() {
		XMLSchemaComplexType ct;
		ct.name( "seq_of_one" );
		XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		seq->append_particle( elem1 );
		ct.set_model_group( seq );
		std::string definition = "<xs:complexType name=\"seq_of_one\" mixed=\"true\">\n <xs:sequence>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n </xs:sequence>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_sequence2() {
		XMLSchemaComplexType ct;
		ct.name( "seq_of_two" );
		XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		seq->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( seq );
		std::string definition = "<xs:complexType name=\"seq_of_two\" mixed=\"true\">\n <xs:sequence>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:sequence>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_sequence_w_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "seq_of_two" );
		XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		seq->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( seq );
		ct.add_attribute( XMLSchemaAttribute::attribute_w_default( "attr1", xs_string, "",  "1"  )); // default value constructor
		ct.add_attribute( XMLSchemaAttribute::required_attribute( "attr2", xs_boolean , "" )); // is_required constructor

		std::string definition = "<xs:complexType name=\"seq_of_two\" mixed=\"true\">\n <xs:sequence>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:sequence>\n"
			" <xs:attribute name=\"attr1\" type=\"xs:string\" default=\"1\"/>\n"
			" <xs:attribute name=\"attr2\" type=\"xs:boolean\" use=\"required\"/>\n"
			"</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_empty_w_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "empty_w_two_attrs" );
		ct.add_attribute( XMLSchemaAttribute::attribute_w_default( "attr1", xs_string, "",  "1"  )); // default value constructor
		ct.add_attribute( XMLSchemaAttribute::required_attribute( "attr2", xs_boolean , "" )); // is_required constructor

		std::string definition = "<xs:complexType name=\"empty_w_two_attrs\" mixed=\"true\">\n"
			" <xs:attribute name=\"attr1\" type=\"xs:string\" default=\"1\"/>\n"
			" <xs:attribute name=\"attr2\" type=\"xs:boolean\" use=\"required\"/>\n"
			"</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	/// @brief make sure that "xs:sequence" doesn't get written out if there are no sub-elements
	void test_xml_complex_type_empty_but_seq_w_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "empty_w_two_attrs" );
		XMLSchemaModelGroupOP seq( new XMLSchemaModelGroup( xsmgt_sequence ));
		ct.set_model_group( seq );
		ct.add_attribute( XMLSchemaAttribute::attribute_w_default( "attr1", xs_string, "",  "1"  )); // default value constructor
		ct.add_attribute( XMLSchemaAttribute::required_attribute( "attr2", xs_boolean , "" )); // is_required constructor

		std::string definition = "<xs:complexType name=\"empty_w_two_attrs\" mixed=\"true\">\n"
			" <xs:attribute name=\"attr1\" type=\"xs:string\" default=\"1\"/>\n"
			" <xs:attribute name=\"attr2\" type=\"xs:boolean\" use=\"required\"/>\n"
			"</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	// choice
	void test_xml_complex_type_choice1() {
		XMLSchemaComplexType ct;
		ct.name( "choice_of_one" );
		XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		choice->append_particle( elem1 );
		ct.set_model_group( choice );
		std::string definition = "<xs:complexType name=\"choice_of_one\" mixed=\"true\">\n <xs:choice>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n </xs:choice>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_choice2() {
		XMLSchemaComplexType ct;
		ct.name( "choice_of_two" );
		XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		choice->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( choice );
		std::string definition = "<xs:complexType name=\"choice_of_two\" mixed=\"true\">\n <xs:choice>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:choice>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_choice_w_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "choice_of_two" );
		XMLSchemaModelGroupOP choice( new XMLSchemaModelGroup( xsmgt_choice ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		choice->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( choice );
		ct.add_attribute( XMLSchemaAttribute::attribute_w_default( "attr1", xs_string, "",  "1"  )); // default value constructor
		ct.add_attribute( XMLSchemaAttribute::required_attribute( "attr2", xs_boolean , "" )); // is_required constructor

		std::string definition = "<xs:complexType name=\"choice_of_two\" mixed=\"true\">\n <xs:choice>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"attr1\" type=\"xs:string\" default=\"1\"/>\n"
			" <xs:attribute name=\"attr2\" type=\"xs:boolean\" use=\"required\"/>\n"
			"</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	// all
	void test_xml_complex_type_all1() {
		XMLSchemaComplexType ct;
		ct.name( "all_of_one" );
		XMLSchemaModelGroupOP all( new XMLSchemaModelGroup( xsmgt_all ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		all->append_particle( elem1 );
		ct.set_model_group( all );

		std::string definition = "<xs:complexType name=\"all_of_one\" mixed=\"true\">\n <xs:all>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n </xs:all>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_all2() {
		XMLSchemaComplexType ct;
		ct.name( "seq_of_two" );
		XMLSchemaModelGroupOP all( new XMLSchemaModelGroup( xsmgt_all ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		all->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( all );

		std::string definition = "<xs:complexType name=\"seq_of_two\" mixed=\"true\">\n <xs:all>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:all>\n</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_complex_type_all_w_attributes() {
		XMLSchemaComplexType ct;
		ct.name( "seq_of_two" );
		XMLSchemaModelGroupOP all( new XMLSchemaModelGroup( xsmgt_all ));
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		all->append_particle( elem1 ).append_particle( elem2 );
		ct.set_model_group( all );
		ct.add_attribute( XMLSchemaAttribute::attribute_w_default( "attr1", xs_string, "",  "1"  )); // default value constructor
		ct.add_attribute( XMLSchemaAttribute::required_attribute( "attr2", xs_boolean , "" )); // is_required constructor

		std::string definition = "<xs:complexType name=\"seq_of_two\" mixed=\"true\">\n <xs:all>\n"
			"  <xs:element name=\"elem1\" type=\"xs:string\"/>\n  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n"
			" </xs:all>\n"
			" <xs:attribute name=\"attr1\" type=\"xs:string\" default=\"1\"/>\n"
			" <xs:attribute name=\"attr2\" type=\"xs:boolean\" use=\"required\"/>\n"
			"</xs:complexType>\n";
		std::ostringstream oss;
		ct.write_definition( 0, oss );
		TS_ASSERT_EQUALS( definition, oss.str() );
	}


	void test_complex_type_group_sequence() {
		XMLSchemaModelGroupOP seqtype( new XMLSchemaModelGroup( xsmgt_sequence ) );
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		seqtype->append_particle( elem1 ).append_particle( elem2 );

		XMLSchemaModelGroup grouptype( xsmgt_group );
		grouptype.group_name( "group_seq" );
		grouptype.append_particle( seqtype );

		std::ostringstream oss;
		grouptype.write_definition( 0, oss );

		std::string definition = "<xs:group name=\"group_seq\">\n <xs:sequence>\n  <xs:element name=\"elem1\" type=\"xs:string\"/>\n"
			"  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n </xs:sequence>\n</xs:group>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );

	}

	void test_complex_type_group_choice() {
		XMLSchemaModelGroupOP seqtype( new XMLSchemaModelGroup( xsmgt_choice ) );
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		seqtype->append_particle( elem1 ).append_particle( elem2 );

		XMLSchemaModelGroup grouptype( xsmgt_group );
		grouptype.group_name( "group_seq" );
		grouptype.append_particle( seqtype );

		std::ostringstream oss;
		grouptype.write_definition( 0, oss );

		std::string definition = "<xs:group name=\"group_seq\">\n <xs:choice>\n  <xs:element name=\"elem1\" type=\"xs:string\"/>\n"
			"  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n </xs:choice>\n</xs:group>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );

	}

	void test_complex_type_group_all() {
		XMLSchemaModelGroupOP seqtype( new XMLSchemaModelGroup( xsmgt_all ) );
		XMLSchemaElementOP elem1( new XMLSchemaElement ); elem1->name( "elem1" ); elem1->type_name( xs_string );
		XMLSchemaElementOP elem2( new XMLSchemaElement ); elem2->name( "elem2" ); elem2->type_name( xs_integer );
		seqtype->append_particle( elem1 ).append_particle( elem2 );

		XMLSchemaModelGroup grouptype( xsmgt_group );
		grouptype.group_name( "group_seq" );
		grouptype.append_particle( seqtype );

		std::ostringstream oss;
		grouptype.write_definition( 0, oss );

		std::string definition = "<xs:group name=\"group_seq\">\n <xs:all>\n  <xs:element name=\"elem1\" type=\"xs:string\"/>\n"
			"  <xs:element name=\"elem2\" type=\"xs:integer\"/>\n </xs:all>\n</xs:group>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );

	}

	void test_xml_element_type_reference() {
		XMLSchemaElement reslist;
		reslist.name( "residues" );
		reslist.type_name( "int_cslist" );
		std::ostringstream oss;
		reslist.write_definition( 2, oss );
		std::string definition = "  <xs:element name=\"residues\" type=\"int_cslist\"/>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_element_type_reference_w_occurs() {
		XMLSchemaElement reslist;
		reslist.name( "residues" );
		reslist.type_name( "int_cslist" );
		reslist.min_occurs( 3 );
		reslist.max_occurs( 10 );
		std::ostringstream oss;
		reslist.write_definition( 2, oss );
		std::string definition = "  <xs:element name=\"residues\" type=\"int_cslist\" minOccurs=\"3\" maxOccurs=\"10\"/>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_element_group_reference() {
		XMLSchemaModelGroup selectors;
		selectors.group_name( "residue_selector" );
		std::ostringstream oss;
		selectors.write_definition( 2, oss );
		std::string definition = "  <xs:group ref=\"residue_selector\"/>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );
	}


	void test_xml_element_group_reference_w_occurs() {
		XMLSchemaModelGroup selectors;
		selectors.group_name( "residue_selector" );
		selectors.min_occurs( 2 );
		selectors.max_occurs( xsminmax_unbounded );
		std::ostringstream oss;
		selectors.write_definition( 2, oss );
		std::string definition = "  <xs:group ref=\"residue_selector\" minOccurs=\"2\" maxOccurs=\"unbounded\"/>\n";
		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	void test_xml_element_w_complex_type() {
		XMLSchemaElement and_selector;
		and_selector.name( "AndSelector" );

		XMLSchemaModelGroupOP selectors( new XMLSchemaModelGroup );
		selectors->group_name( "residue_selector" );
		selectors->min_occurs( 2 );
		selectors->max_occurs( xsminmax_unbounded );

		XMLSchemaComplexTypeOP and_type( new XMLSchemaComplexType );
		and_type->set_model_group( selectors );
		and_type->add_attribute( XMLSchemaAttribute( "some_attribute", xs_string , "" ));

		and_selector.element_type_def( and_type );

		std::ostringstream oss;
		and_selector.write_definition( 0, oss );
		std::string definition = "<xs:element name=\"AndSelector\">\n <xs:complexType mixed=\"true\">\n"
			"  <xs:group ref=\"residue_selector\" minOccurs=\"2\" maxOccurs=\"unbounded\"/>\n"
			"  <xs:attribute name=\"some_attribute\" type=\"xs:string\"/>\n"
			" </xs:complexType>\n</xs:element>\n";

		TS_ASSERT_EQUALS( definition, oss.str() );
	}

	//void dont_test_xml_schema_element_w_restriction_type() {
	//	XMLSchemaElement element;
	//	element.name( "residues" );
	//	XMLSchemaRestrictionOP cs_intlist( new XMLSchemaRestriction );
	//	cs_intlist->base_type( xs_string );
	//	cs_intlist->add_restriction( xsr_pattern, "[0-9]+(,[0-9]+)*" );
	//	element.restriction_type_def( cs_intlist );
	//
	//	std::ostringstream oss;
	//	element.write_definition( 0, oss );
	//
	//	std::string definition = "<xs:element name=\"residues\">\n <xs:simpleType>\n  <xs:restriction base=\"xs:string\">\n"
	//		"   <xs:pattern value=\"[0-9]+(,[0-9]+)*\"/>\n"
	//		"  </xs:restriction>\n </xs:simpleType>\n</xs:element>\n";
	//
	//	TS_ASSERT_EQUALS( definition, oss.str() );
	//}

	void test_xml_schema_definition_detect_name_collision() {
		XMLSchemaDefinition xsd;
		XMLSchemaElement rs_str; rs_str.name( "RosettaScripts" ); rs_str.type_name( "xs:string" );
		XMLSchemaElement rs_int; rs_int.name( "RosettaScripts" ); rs_int.type_name( "xs:integer" );
		xsd.add_top_level_element( rs_str );
		//xsd.add_top_level_element( "RosettaScripts", "<xs:element name=\"RosettaScripts\" type=\"xs:string\"/>\n" );
		try {
			xsd.add_top_level_element( rs_int );
			//xsd.add_top_level_element( "RosettaScripts", "<xs:element name=\"RosettaScripts\" type=\"xs:integer\"/>\n" );
			TS_ASSERT( false ); // we should not have reached this point
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			TS_ASSERT_EQUALS( e.msg(), "Name collision in creation of XML Schema definition: top level element with name \"RosettaScripts\" already has been defined.\n"
				"Old definition:\n<xs:element name=\"RosettaScripts\" type=\"xs:string\"/>\nNew definition:\n"
				"<xs:element name=\"RosettaScripts\" type=\"xs:integer\"/>\n");
		}
	}

	void test_xml_schema_definition_detect_name_mismatch() {
		XMLSchemaDefinition xsd;
		class imposter : public XMLSchemaTopLevelElement {
		public:
			imposter() : name_( "imposter" ) {}
			virtual std::string const & element_name() const { return name_; }
			virtual void write_definition( int indentation, std::ostream & os ) const {
				for ( int ii = 1; ii <= indentation; ++ii ) os << " ";
				os << "<xs:element name=\"RosettaScripts\" type=\"xs:integer\"/>\n" << std::endl;
			}
			virtual void prepare_for_output( XMLSchemaDefinition & ) const {}

		private:
			std::string name_;
		} imp;

		try {
			xsd.add_top_level_element( imp );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			TS_ASSERT_EQUALS( e.msg(), "top level tag with presumed name \"imposter\" has an actual name of \"RosettaScripts\"" );
		}
	}

	void test_xml_schema_definition_detect_missing_name() {
		XMLSchemaDefinition xsd;
		class imposter : public XMLSchemaTopLevelElement {
		public:
			imposter() : name_( "RosettaScripts" ) {}
			virtual std::string const & element_name() const { return name_; }
			virtual void write_definition( int indentation, std::ostream & os ) const {
				for ( int ii = 1; ii <= indentation; ++ii ) os << " ";
				os << "<xs:element type=\"xs:integer\"/>\n" << std::endl;
			}
			virtual void prepare_for_output( XMLSchemaDefinition & ) const {}

		private:
			std::string name_;
		} imp;
		try {
			xsd.add_top_level_element( imp );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			TS_ASSERT_EQUALS( e.msg(), "top level tag \"RosettaScripts\" does not have a name attribute (a.k.a. option)." );
		}
	}

	void test_xml_schema_definition_write_xsd() {
		XMLSchemaDefinition xsd;

		ResSelectCT res_selector_ct;
		xsd.add_top_level_element( res_selector_ct );

		ResSelectGrp res_selector_grp;
		xsd.add_top_level_element( res_selector_grp );

		std::string xsd_expected_output =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
			" <xs:group ref=\"ResidueSelector\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:group name=\"ResidueSelector\">\n"
			" <xs:choice>\n"
			"  <xs:element name=\"AndResidueSelector\" type=\"AndResidueSelectorType\"/>\n"
			"  <xs:element name=\"OrResidueSelector\" type=\"OrResidueSelectorType\"/>\n"
			"  <xs:element name=\"NotResidueSelector\" type=\"NotResidueSelectorType\"/>\n"
			"  <xs:element name=\"ChainResidueSelector\" type=\"ChainResidueSelectorType\"/>\n"
			"  <xs:element name=\"NeighborhoodResidueSelector\" type=\"NeighborhoodResidueSelectorType\"/>\n"
			" </xs:choice>\n"
			"</xs:group>\n"
			"\n"
			"</xs:schema>\n";

		//std::cout << xsd_expected_output << "\nand\n" << xsd.full_definition() << std::endl;


		TS_ASSERT_EQUALS( xsd.full_definition(), xsd_expected_output );

	}

	void test_xml_schema_definition_write_xsd_w_identical_duplicate() {
		XMLSchemaDefinition xsd;

		ResSelectCT res_selector_ct;
		xsd.add_top_level_element( res_selector_ct );

		ResSelectGrp res_selector_grp;
		xsd.add_top_level_element( res_selector_grp );

		// now for the duplicated element.  This addition should not throw
		// an exception because the definition should match perfectly.
		ResSelectCT res_selector_ct2;
		xsd.add_top_level_element( res_selector_ct2 );

		std::string xsd_expected_output =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
			" <xs:group ref=\"ResidueSelector\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:group name=\"ResidueSelector\">\n"
			" <xs:choice>\n"
			"  <xs:element name=\"AndResidueSelector\" type=\"AndResidueSelectorType\"/>\n"
			"  <xs:element name=\"OrResidueSelector\" type=\"OrResidueSelectorType\"/>\n"
			"  <xs:element name=\"NotResidueSelector\" type=\"NotResidueSelectorType\"/>\n"
			"  <xs:element name=\"ChainResidueSelector\" type=\"ChainResidueSelectorType\"/>\n"
			"  <xs:element name=\"NeighborhoodResidueSelector\" type=\"NeighborhoodResidueSelectorType\"/>\n"
			" </xs:choice>\n"
			"</xs:group>\n"
			"\n"
			"</xs:schema>\n";

		//std::cout << xsd_expected_output << "\nand\n" << xsd.full_definition() << std::endl;


		TS_ASSERT_EQUALS( xsd.full_definition(), xsd_expected_output );

	}

	void test_xml_schema_definition_write_xsd_w_nonidentical_duplicate() {
		XMLSchemaDefinition xsd;

		ResSelectCT res_selector_ct;
		xsd.add_top_level_element( res_selector_ct );

		ResSelectGrp res_selector_grp;
		xsd.add_top_level_element( res_selector_grp );



		try {
			class ResSelectCT2 : public XMLSchemaTopLevelElement {
			public:
				ResSelectCT2() : name_( "ResidueSelectorType" ) {}
				virtual std::string const & element_name() const { return name_; }
				virtual void write_definition( int indentation, std::ostream & os ) const {
					for ( int ii = 1; ii <= indentation; ++ii ) os << " ";
					os <<
						"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
						" <xs:group ref=\"TResidueSelector\"/>\n" // typo! TResidueSelector instead of ResidueSelector
						"</xs:complexType>\n";
				}
				virtual void prepare_for_output( XMLSchemaDefinition & ) const {}
			private:
				std::string name_;
			} res_select_ct2;

			// now for the duplicated element
			xsd.add_top_level_element( res_select_ct2 );
			TS_ASSERT( false ); // this should not be reached
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::string expected_message = "Name collision in creation of XML Schema definition: top level element with "
				"name \"ResidueSelectorType\" already has been defined.\n"
				"Old definition:\n"
				"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
				" <xs:group ref=\"ResidueSelector\"/>\n"
				"</xs:complexType>\n"
				"New definition:\n"
				"<xs:complexType name=\"ResidueSelectorType\" mixed=\"true\">\n"
				" <xs:group ref=\"TResidueSelector\"/>\n"
				"</xs:complexType>\n";
			TS_ASSERT_EQUALS( e.msg(), expected_message );
			// std::cout << std::endl << std::endl;
			// std::cout << e.msg() << std::endl << std::endl;
			// std::cout << expected_message << std::endl << std::endl;
			// platform::Size first_diff = first_difference( e.msg(), expected_message );
			// if ( first_diff < e.msg().size() && first_diff < expected_message.size() ) {
			// 	std::cout << "first difference between messages: pos " << first_diff+1 << " '" << e.msg()[ first_diff ] << "' vs '"
			// 						<< expected_message[ first_diff ] << "'" << std::endl;
			// } else {
			// 	std::cout << "first diff pos " << first_diff << " vs sizes " << e.msg().size() << " and " << expected_message.size() << std::endl;
			// }
		}
	}

	static std::string dummy_ct_naming_func( std::string const & name ) { return "test_" + name + "Type"; }
	static std::string dummy_ct_naming_func2( std::string const & name ) { return "test2_" + name + "Type"; }
	static std::string dummy_group_name_func() { return "test_group"; }

	void test_xml_complex_type_schema_generator_simple_ct() {
		XMLSchemaComplexTypeGenerator ctgen;
		ctgen.element_name( "Simple" );
		ctgen.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );
		AttributeList attributes;
		attributes.push_back( XMLSchemaAttribute( "name", xs_string , "" ));
		attributes.push_back( XMLSchemaAttribute( "chains", xs_string , "" ));
		ctgen.add_attributes( attributes );
		ctgen.description( "" );

		XMLSchemaDefinition xsd;
		ctgen.write_complex_type_to_schema( xsd );

		//std::cout << "test_xml_complex_type_schema_generator_simple_ct:\n" << xsd.full_definition() << std::endl;
		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"test_SimpleType\" mixed=\"true\">\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"chains\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";
		TS_ASSERT_EQUALS( gold, xsd.full_definition() );
	}

	void test_xml_complex_type_schema_generator_subelement_ctref() {
		XMLSchemaComplexTypeGenerator ctgen;
		ctgen.element_name( "Testing123" );
		ctgen.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );

		XMLSchemaSimpleSubelementList subelements;
		subelements.add_already_defined_subelement( "Helix", boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );
		subelements.add_already_defined_subelement( "Sheet", boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );
		ctgen.set_subelements_single_appearance_required( subelements );
		ctgen.description( "" );

		XMLSchemaDefinition xsd;
		ctgen.write_complex_type_to_schema( xsd );

		//std::cout << "test_xml_complex_type_schema_generator_simple_ct:\n" << xsd.full_definition() << std::endl;
		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"test_Testing123Type\" mixed=\"true\">\n"
			" <xs:all>\n"
			"  <xs:element name=\"Helix\" type=\"test_HelixType\"/>\n"
			"  <xs:element name=\"Sheet\" type=\"test_SheetType\"/>\n"
			" </xs:all>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";
		TS_ASSERT_EQUALS( gold, xsd.full_definition() );
	}


	void test_xml_complex_type_schema_generator_group_ct() {
		XMLSchemaComplexTypeGenerator ctgen;
		ctgen.element_name( "Group" );
		ctgen.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );

		AttributeList attributes;
		attributes.push_back( XMLSchemaAttribute( "name", xs_string , "" ));
		attributes.push_back( XMLSchemaAttribute( "chains", xs_string , "" ));
		ctgen.add_attributes( attributes );
		ctgen.description( "" );

		AttributeList subelement_attributes;
		subelement_attributes.push_back( XMLSchemaAttribute( "aa", xs_string , "" ));
		subelement_attributes.push_back( XMLSchemaAttribute( "appends", xs_string , "" ));
		subelement_attributes.push_back( XMLSchemaAttribute( "removes", xs_string , "" ));

		XMLSchemaSimpleSubelementList subelements;
		subelements.add_simple_subelement( "Helix", subelement_attributes , "" );
		subelements.add_simple_subelement( "Strand", subelement_attributes , "" );
		subelements.add_simple_subelement( "HelixNTerm", subelement_attributes , "" );
		subelements.add_simple_subelement( "HelixCTerm", subelement_attributes , "" );
		subelements.add_simple_subelement( "Loop", subelement_attributes , "" );
		subelements.add_simple_subelement( "all", subelement_attributes , "" );
		subelements.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func2, _1 ));

		ctgen.set_subelements_repeatable( subelements );

		XMLSchemaDefinition xsd;
		ctgen.write_complex_type_to_schema( xsd );

		//std::cout << "test_xml_complex_type_schema_generator_group_ct:\n" << xsd.full_definition() << std::endl;
		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"test2_HelixType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test2_StrandType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test2_HelixNTermType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test2_HelixCTermType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test2_LoopType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test2_allType\" mixed=\"true\">\n"
			" <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"test_GroupType\" mixed=\"true\">\n"
			" <xs:choice minOccurs=\"0\" maxOccurs=\"unbounded\">\n"
			"  <xs:element name=\"Helix\" type=\"test2_HelixType\"/>\n"
			"  <xs:element name=\"Strand\" type=\"test2_StrandType\"/>\n"
			"  <xs:element name=\"HelixNTerm\" type=\"test2_HelixNTermType\"/>\n"
			"  <xs:element name=\"HelixCTerm\" type=\"test2_HelixCTermType\"/>\n"
			"  <xs:element name=\"Loop\" type=\"test2_LoopType\"/>\n"
			"  <xs:element name=\"all\" type=\"test2_allType\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"chains\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";
		TS_ASSERT_EQUALS( gold, xsd.full_definition() );

	}

	void test_xml_complex_type_schema_generator_group_ct_no_inner_naming_func() {
		XMLSchemaComplexTypeGenerator ctgen;
		ctgen.element_name( "Group" );
		ctgen.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func, _1 ) );

		AttributeList attributes;
		attributes.push_back( XMLSchemaAttribute( "name", xs_string , "" ));
		attributes.push_back( XMLSchemaAttribute( "chains", xs_string , "" ));
		ctgen.add_attributes( attributes );
		ctgen.description( "testing 123" );

		AttributeList subelement_attributes;
		subelement_attributes.push_back( XMLSchemaAttribute( "aa", xs_string , "" ));
		subelement_attributes.push_back( XMLSchemaAttribute( "appends", xs_string , "" ));
		subelement_attributes.push_back( XMLSchemaAttribute( "removes", xs_string , "" ));

		XMLSchemaSimpleSubelementList subelements;
		subelements.add_simple_subelement( "Helix", subelement_attributes , "" );
		subelements.add_simple_subelement( "Strand", subelement_attributes , "" );
		subelements.add_simple_subelement( "HelixNTerm", subelement_attributes , "" );
		subelements.add_simple_subelement( "HelixCTerm", subelement_attributes , "" );
		subelements.add_simple_subelement( "Loop", subelement_attributes , "" );
		subelements.add_simple_subelement( "all", subelement_attributes , "" );
		// do not set a naming function for the subelements
		// subelements.complex_type_naming_func( boost::bind( XMLSchemaGenerationTests::dummy_ct_naming_func2, _1 ));

		ctgen.set_subelements_repeatable( subelements );

		XMLSchemaDefinition xsd;
		ctgen.write_complex_type_to_schema( xsd );

		//std::cout << "test_xml_complex_type_schema_generator_group_ct:\n" << xsd.full_definition() << std::endl;

		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"test_GroupType\" mixed=\"true\">\n"
			" <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"  testing 123\n"
			" </xs:documentation></xs:annotation>\n"
			" <xs:choice minOccurs=\"0\" maxOccurs=\"unbounded\">\n"
			"  <xs:element name=\"Helix\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"Strand\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"HelixNTerm\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"HelixCTerm\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"Loop\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"all\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"aa\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"appends\" type=\"xs:string\"/>\n"
			"    <xs:attribute name=\"removes\" type=\"xs:string\"/>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"name\" type=\"xs:string\"/>\n"
			" <xs:attribute name=\"chains\" type=\"xs:string\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";

		TS_ASSERT_EQUALS( gold, xsd.full_definition() );
		if ( gold != xsd.full_definition() ) {
			compare_line_by_line( gold, xsd.full_definition() );
		}

	}

	void test_complex_type_generator_seq_of_sets1()
	{
		XMLSchemaComplexTypeGenerator ctgen;

		XMLSchemaSimpleSubelementList input_element;
		input_element.add_already_defined_subelement( "Input", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_required( input_element );

		XMLSchemaSimpleSubelementList output_element;
		output_element.add_already_defined_subelement( "Output", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_required( output_element );

		XMLSchemaSimpleSubelementList options_element;
		options_element.add_already_defined_subelement( "Options", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_optional( options_element );

		XMLSchemaSimpleSubelementList task_ops_and_sfxns;
		task_ops_and_sfxns.add_already_defined_subelement( "TaskOperations", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		task_ops_and_sfxns.add_already_defined_subelement( "ScoreFunctions", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_repeatable( task_ops_and_sfxns );

		XMLSchemaSimpleSubelementList packer_types;
		packer_types.add_already_defined_subelement( "PackRotamersMover", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		packer_types.add_already_defined_subelement( "MinPackMover", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_pick_one( packer_types  );

		XMLSchemaSimpleSubelementList min_types;
		min_types.add_already_defined_subelement( "MinMover", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		min_types.add_already_defined_subelement( "WildMinMover", & XMLSchemaGenerationTests::dummy_ct_naming_func );
		ctgen.add_ordered_subelement_set_as_pick_one_or_none( min_types  );

		TS_ASSERT( ctgen.subelement_behavior() == se_ordered_sets );
		XMLSchemaDefinition xsd;
		ctgen
			.element_name( "Job" )
			.complex_type_naming_func( & XMLSchemaGenerationTests::dummy_ct_naming_func )
			.add_attribute( XMLSchemaAttribute( "nstruct", xs_integer , "" ) )
			.description( "testing 123" )
			.write_complex_type_to_schema( xsd );

		//std::cout << "test_complex_type_generator_seq_of_sets1:\n" << xsd.full_definition() << std::endl;
		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"test_JobType\" mixed=\"true\">\n"
			" <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"  testing 123\n"
			" </xs:documentation></xs:annotation>\n"
			" <xs:sequence>\n"
			"  <xs:element name=\"Input\" type=\"test_InputType\" minOccurs=\"1\" maxOccurs=\"1\"/>\n"
			"  <xs:element name=\"Output\" type=\"test_OutputType\" minOccurs=\"1\" maxOccurs=\"1\"/>\n"
			"  <xs:element name=\"Options\" type=\"test_OptionsType\" minOccurs=\"0\" maxOccurs=\"1\"/>\n"
			"  <xs:choice minOccurs=\"0\" maxOccurs=\"unbounded\">\n"
			"   <xs:element name=\"TaskOperations\" type=\"test_TaskOperationsType\"/>\n"
			"   <xs:element name=\"ScoreFunctions\" type=\"test_ScoreFunctionsType\"/>\n"
			"  </xs:choice>\n"
			"  <xs:choice minOccurs=\"1\" maxOccurs=\"1\">\n"
			"   <xs:element name=\"PackRotamersMover\" type=\"test_PackRotamersMoverType\"/>\n"
			"   <xs:element name=\"MinPackMover\" type=\"test_MinPackMoverType\"/>\n"
			"  </xs:choice>\n"
			"  <xs:choice minOccurs=\"0\" maxOccurs=\"1\">\n"
			"   <xs:element name=\"MinMover\" type=\"test_MinMoverType\"/>\n"
			"   <xs:element name=\"WildMinMover\" type=\"test_WildMinMoverType\"/>\n"
			"  </xs:choice>\n"
			" </xs:sequence>\n"
			" <xs:attribute name=\"nstruct\" type=\"xs:integer\"/>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";
		TS_ASSERT_EQUALS( gold, xsd.full_definition() );
		if ( gold != xsd.full_definition() ) {
			compare_line_by_line( gold, xsd.full_definition() );
		}

	}

	static std::string root_mangler( std::string const & element ) { return "root_mangler_" + element + "_type"; }
	static std::string mangler1( std::string const & element ) { return "mangler1_" + element + "_type"; }
	static std::string mangler2( std::string const & element ) { return "mangler2_" + element + "_type"; }
	//static std::string mangler3( std::string const & element ) { return "mangler3_" + element + "_type"; }
	//static std::string mangler4( std::string const & element ) { return "mangler4_" + element + "_type"; }
	//static std::string mangler5( std::string const & element ) { return "mangler5_" + element + "_type"; }

	void test_repeatable_complex_type_node()
	{
		// <Element1 att1="foo" att2="blah">
		//   <Element2>
		//     <Element3/>
		//     <Element4/>
		//   </Element2>
		//   <Element5/>
		// <Element1>
		//
		// faux parse_my_tag( tag ):
		//   for ( subtag : tag->getTags() )
		//      if subtag->getName() == "Element2" :
		//         for ( subsubtag : subtag->getTags() )
		//            if subsubtag->getName() == "Element3" :
		//               ... do stuff
		//            else if subsubtag->getName() == "Element4" :
		//               ... do other stuff
		//      else if subtag->getName() == "Element5" :
		//         ... do completely different stuff
		typedef XMLSchemaAttribute Attr;
		AttributeList attlist1, attlist2, attlist3, attlist4, attlist5;
		attlist1 + Attr( "att1", xs_string, "testing" ) + Attr( "att2", xs_string, "testing" );
		attlist2 + Attr( "att3", xs_string, "testing" ) + Attr( "att4", xs_string, "testing" );
		attlist3 + Attr( "att5", xs_string, "testing" ) + Attr( "att6", xs_string, "testing" );
		attlist4 + Attr( "att7", xs_string, "testing" ) + Attr( "att8", xs_string, "testing" );
		attlist5 + Attr( "att9", xs_string, "testing" ) + Attr( "att10", xs_string, "testing" );

		XMLSchemaRepeatableCTNodeOP node1( new XMLSchemaRepeatableCTNode );
		XMLSchemaRepeatableCTNodeOP node2( new XMLSchemaRepeatableCTNode );
		XMLSchemaRepeatableCTNodeOP node3( new XMLSchemaRepeatableCTNode );
		XMLSchemaRepeatableCTNodeOP node4( new XMLSchemaRepeatableCTNode );
		XMLSchemaRepeatableCTNodeOP node5( new XMLSchemaRepeatableCTNode );

		node1->set_element_w_attributes( "Element1", attlist1, "testing" ).set_kids_naming_func( & mangler1 ).set_root_node_naming_func( & root_mangler );
		node2->set_element_w_attributes( "Element2", attlist2, "testing" ).set_kids_naming_func( & mangler2 );
		node3->set_element_w_attributes( "Element3", attlist3, "testing" );
		node4->set_element_w_attributes( "Element4", attlist4, "testing" );
		node5->set_element_w_attributes( "Element5", attlist5, "testing" );

		node1->add_child( node2 ).add_child( node5 );
		node2->add_child( node3 ).add_child( node4 );

		XMLSchemaDefinition xsd;
		node1->recursively_write_ct_to_schema( xsd );

		//std::cout << "XMLSchemaRepeatableCTNode schema:\n" << xsd.full_definition() << std::endl;

		std::string gold =
			"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
			"<xs:schema xmlns:xs=\"http://www.w3.org/2001/XMLSchema\">\n"
			"\n"
			"<xs:complexType name=\"mangler1_Element2_type\" mixed=\"true\">\n"
			" <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"  testing\n"
			" </xs:documentation></xs:annotation>\n"
			" <xs:choice minOccurs=\"0\" maxOccurs=\"unbounded\">\n"
			"  <xs:element name=\"Element3\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"att5\" type=\"xs:string\">\n"
			"     <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"      testing\n"
			"      </xs:documentation></xs:annotation>\n"
			"    </xs:attribute>\n"
			"    <xs:attribute name=\"att6\" type=\"xs:string\">\n"
			"     <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"      testing\n"
			"      </xs:documentation></xs:annotation>\n"
			"    </xs:attribute>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			"  <xs:element name=\"Element4\">\n"
			"   <xs:complexType mixed=\"true\">\n"
			"    <xs:attribute name=\"att7\" type=\"xs:string\">\n"
			"     <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"      testing\n"
			"      </xs:documentation></xs:annotation>\n"
			"    </xs:attribute>\n"
			"    <xs:attribute name=\"att8\" type=\"xs:string\">\n"
			"     <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"      testing\n"
			"      </xs:documentation></xs:annotation>\n"
			"    </xs:attribute>\n"
			"   </xs:complexType>\n"
			"  </xs:element>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"att3\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			" <xs:attribute name=\"att4\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"mangler1_Element5_type\" mixed=\"true\">\n"
			" <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"  testing\n"
			" </xs:documentation></xs:annotation>\n"
			" <xs:attribute name=\"att9\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			" <xs:attribute name=\"att10\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			"</xs:complexType>\n"
			"\n"
			"<xs:complexType name=\"root_mangler_Element1_type\" mixed=\"true\">\n"
			" <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"  testing\n"
			" </xs:documentation></xs:annotation>\n"
			" <xs:choice minOccurs=\"0\" maxOccurs=\"unbounded\">\n"
			"  <xs:element name=\"Element2\" type=\"mangler1_Element2_type\"/>\n"
			"  <xs:element name=\"Element5\" type=\"mangler1_Element5_type\"/>\n"
			" </xs:choice>\n"
			" <xs:attribute name=\"att1\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			" <xs:attribute name=\"att2\" type=\"xs:string\">\n"
			"  <xs:annotation><xs:documentation xml:lang=\"en\">\n"
			"   testing\n"
			"   </xs:documentation></xs:annotation>\n"
			" </xs:attribute>\n"
			"</xs:complexType>\n"
			"\n"
			"</xs:schema>\n";

		TS_ASSERT_EQUALS( gold, xsd.full_definition() );
		if ( gold != xsd.full_definition() ) {
			compare_line_by_line( gold, xsd.full_definition() );
		}

	}

};
