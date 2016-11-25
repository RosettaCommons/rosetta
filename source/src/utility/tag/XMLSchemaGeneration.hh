// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/tag/XMLSchemaGeneration.hh
/// @brief  Declaration of the classes used to define an XML Schema
///
/// @details A Brief Overview of XML Schema and how it will be used in Rosetta:
///
/// XML Schema describes the structure of an XML file.  The objects that
/// appear in the angle-brackets are called "Elements", e.g.
///
/// \code{.xml}
/// <MinMover name="minmover" scorefxn="talaris2014"/>
/// \endcode
///
/// is an element.
/// (aka "Tags" in Rosetta speak). The structure of the elements can be
/// described by types and these types can be user defined, e.g. the two elements:
///
/// \code{.xml}
/// <MinMover name="minmover1" scorefxn="talaris2014"/>, and
/// <MinMover name="minmover2" scorefxn="talaris2014_cst"/>
/// \endcode
///
/// could both be described as two instances that have the same type.
/// At least the way that they are used in rosetta, element types
/// almost universally fall into the XML Schema category of "complex types."
/// Complex types can contain "attributes" (e.g. the "name" and "scorefxn"
/// attributes, ) and they can also contain sub-elements organized
/// into blocks of "model groups" (sequence, choice, all & group).  Attributes
/// (aka "options" in Rosetta speak) are usually described by simple types,
/// and XML Schema defines many of them by default, (e.g. "string" and "integer"
/// types); they can also be described by "restrictions" to existing types, e.g.
/// you could define a type that is a string composed of integers separated
/// by commas (e.g. "1,34,42").
///
/// In a picture:
///
/// \verbatim
/// XMLSchemaElement --> XMLSchemaComplexType --> XMLSchemaAttribute --> XMLSchemaType --> Primative Type
///                 \--> XMLSchemaType       \--> XMLSchemaModelGroup                  \--> Restriction --> Primative Type
///                                                         \----> XMLSchemaElement
///                                                          \---> XMLSchemaModelGroup
/// \endverbatim
///
/// The chunk of XMLSchema that would describe the MinMover element above would
/// look like this: (If you thought XML was verbose, you will be overwhemed by
/// XML Schema!)
///
/// \verbatim
/// <xs:complexType name="MinMoverType" mixed="True">
///  <xs:attribute name="name" type="xs:string"/>
///  <xs:attribute name="scorefxn" type="xs:string"/>
/// </xs:complexType>
///
/// <xs:element name="MinMover" type="MinMoverType"/>
/// \endverbatim
///
/// In a XMLSchema, you can define restrictions and complex types as free standing
/// blocks, and these blocks can be given in any order.  Somewhere in the XMLSchema,
/// there needs to be a free-standing XML schema element, and this element will be
/// interpretted as the root of a tree so that an "order" can be discerned from what
/// is otherwise a jumble of types.
///
/// The way I envision schemas being defined in Rosetta is that an "XMLSchemaDefinition"
/// object is passed around between static methods of classes that are seeking to define
/// the grammar of their definitions, and other non-member utility functions (e.g. one that
/// defines the comma-separated-integer-list XMLRestriction), and that perhaps these utility
/// functions are called multiple times (with the XMLSchemaDefinition ignoring duplicate
/// definitions). After every class has described the structure of its organization, the
/// XMLSchemaDefinition can be written to a file.  Ideally, the construction of an XML schema
/// will not require the instantiation of any of the classes described in the XSD; e.g.
/// to describe the schema of the MinMover, the MoverFactory will ask the MinMoverCreator
/// for its XSD, which will defer to a static method of the MinMover. There would be no
/// need to actually instantiate a MinMover in this process.
///
/// One note about the way we use XML in rosetta:
/// When most people use XML, they embed information between two tags, e.g.
/// \code{.xml}
/// <FirstName> Andrew </FirstName>
/// \endcode
/// and we instead embed all our information as attributes (aka options) within
/// a single tag, e.g.
/// \code{.xml}
/// <FirstName value="Andrew"/>
/// \endcode
/// So our XML is a little funny. Instead of storing data between two tags, we ignore
/// everything outside of a tag.  This is allowed in XML, and so we can still used
/// XML Schema to describe our XML files.  In order to do so, we have to use the
/// feature of describing complex types where we describe the tags as "mixed". See
///
/// http://www.w3schools.com/xml/schema_complex_mixed.asp
///
/// The XMLSchemaComplexType implemented here describes all complex types as mixed;
/// this may not be appropriate for non-Rosetta XML.
///
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_tag_XMLSchemaGeneration_HH
#define INCLUDED_utility_tag_XMLSchemaGeneration_HH

// Unit headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// Boost headers
#include <boost/function.hpp>

// C++ headers
#include <list>
#include <vector>
#include <map>
#include <iosfwd>
#include <sstream>

namespace utility {
namespace tag {

/// @brief An enum listing lots of commonly-used simple types so they don't have to be added manually
/// in lots of places.  If you create an XMLSchemaType using this enum, then the corresponding
/// restriction will be added to the XMLSchemaDefintion automatically
enum XMLSchemaCommonType {
	xsct_none,
	xsct_int_cslist,
	xsct_int_wsslist,
	xsct_nnegative_int_cslist,
	xsct_nnegative_int_wsslist,
	xsct_real,
	xsct_real_cslist,
	xsct_real_cslist_w_ws,
	xsct_real_wsslist,
	xsct_bool_cslist,
	xsct_bool_wsslist,
	xsct_positive_integer_cslist,
	xsct_positive_integer_wsslist,
	xsct_non_negative_integer,
	xsct_positive_integer,
	xsct_rosetta_bool,
	xsct_residue_number,
	xsct_residue_number_cslist,
	xsct_refpose_enabled_residue_number,
	xsct_refpose_enabled_residue_number_cslist,
	xsct_char,
	xsct_minimizer_type,
	xsct_size_cs_pair,
	xsct_chain_cslist,
	xsct_dssp_string
};
std::string residue_number_string();
std::string real_regex_pattern();
std::string refpose_enabled_residue_number_string();
std::string name_for_common_type( XMLSchemaCommonType common_type );
std::ostream & operator << ( std::ostream & os, XMLSchemaCommonType common_type );

std::string chr_chains_nonrepeated();

/// @brief class %XMLSchemaType represents the name of a defined type that can
/// be used to describe either an XMLElement or an XMLAttribute.  It may refer
/// to either a complex type or to a primative type or to a simple type.
class XMLSchemaType : public utility::pointer::ReferenceCount {
public:
	XMLSchemaType();
	XMLSchemaType( XMLSchemaDataType setting );
	XMLSchemaType( XMLSchemaCommonType setting );
	XMLSchemaType( std::string const & custom_type );
	XMLSchemaType( char const * custom_type );
	void type( XMLSchemaDataType setting );
	void common_type( XMLSchemaCommonType setting );
	void custom_type_name( std::string const & setting );

	std::string type_name() const;
	XMLSchemaDataType type() const;
	XMLSchemaCommonType common_type() const;


private:
	XMLSchemaDataType type_;
	XMLSchemaCommonType common_type_;
	std::string custom_type_name_;
};

/// @brief The %XMLSchemaTopLevelElement class signifies a class that could be written out
/// as an XML element, with possible sub-elemets, in an XML schema. When generating a
/// schema, the developer will hand an instance of a class derived from an %XMLSchemaTopLevelElement
/// to an XMLSchemaDefinition object.
class XMLSchemaTopLevelElement : public utility::pointer::ReferenceCount
{
public:
	virtual std::string const & element_name() const = 0;
	virtual void write_definition( int indentation, std::ostream & os ) const = 0;
	virtual void prepare_for_output( XMLSchemaDefinition & xsd ) const = 0;
};

/// @brief class %XMLSchemaAttribute represents what we refer to in Rosetta as an
/// option for a tag.  An attribute would reside inside of a tag, such as "scorefxn"
/// in this tag (this XML Element): <MinMover name="min" scorefxn="talaris2014"/>
///
/// @details The "name" of an attribute is what will live on the left-hand side of the
/// equals sign, e.g. "scorefxn" in the MinMover example above.  The type of the attribute
/// refers to the structure of the data on the right-hand side of the equals sign.
/// Attributes may also define a default value, which is a useful way to communicate default
/// behaviors to users.  Finally, an attribute can be listed as being required or not.
/// By default, attributes are not required, so not specifying whether the attribute is
/// required means that it is not required.  An example of a required attribute would be
/// the "name" attribute for the MinMover.
///
/// Once all of the details of an XMLAttribute have been set, it will write its definition
/// to an output stream in its write_definition method.
class XMLSchemaAttribute : public XMLSchemaTopLevelElement {
public:
	XMLSchemaAttribute();
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, std::string const & desc );
	//XMLSchemaAttribute( std::string const & name, XMLSchemaType type, std::string const & default_value ); // temp -- find attribtues w/ default values
	//static XMLSchemaAttribute required_attribute( std::string const & name, XMLSchemaType type );
	static XMLSchemaAttribute required_attribute( std::string const & name, XMLSchemaType type, std::string const & desc );
	//static XMLSchemaAttribute attribute_w_default( std::string const & name, XMLSchemaType type, std::string const & default_value );
	static XMLSchemaAttribute attribute_w_default( std::string const & name, XMLSchemaType type, std::string const & desc, std::string const & default_value );

	XMLSchemaAttribute & name( std::string const & setting );
	XMLSchemaAttribute & type( XMLSchemaType setting );
	XMLSchemaAttribute & default_value( std::string const & setting );
	XMLSchemaAttribute & description( std::string const & setting );
	XMLSchemaAttribute & is_required( bool setting );

	XMLSchemaType const & type() const;
	std::string const & element_name() const override;
	void write_definition( int indentation, std::ostream & os ) const override;
	void prepare_for_output( XMLSchemaDefinition & xsd ) const override;

private:
	std::string name_;
	XMLSchemaType type_;
	std::string default_value_;
	bool is_required_;
	std::string description_;
};

AttributeList &
operator + ( AttributeList & attributes, XMLSchemaAttribute const & attribute_to_append );

/// @brief class %XMLSchemaRestriction describes a refinement on the behavior of existing
/// types.  For example, one could define a restriction representing a list of residue indexes
/// separating commas: "15,44,102" and then describe an attribute of a complex type as having
/// to conform to that restriction.  An xml-schema validator would be able to say that an input
/// file with "fifteen,fortyfour,onehundredandtwo" did not meet the schema.
///
/// @details See the description of XML schema restrictions here:
/// http://www.w3schools.com/xml/schema_facets.asp
/// Restrictions are given as pairs of restriction types and then values for those restriction
/// types.  The restriction types are defined in an enumeration in XMLSchemaGeneration.fwd.hh.
///
/// An example of a useful XMLSchemaRestriction is oen to define a comma-separated list of integers.
/// The elements in the XMLSchema that would define such a restriction look like this:
///
/// \verbatim
/// <xs:simpleType name="int_cslist">
///  <xs:restriction base="xs:string">
///   <xs:pattern value="[0-9]+(,[0-9]+)*"/>
///  </xs:restriction>
/// </xs:simpleType>
/// \endverbatim

class XMLSchemaRestriction : public XMLSchemaTopLevelElement {
public:
	XMLSchemaRestriction();
	void name( std::string const & setting );
	void base_type( XMLSchemaType setting );
	void add_restriction( XMLSchemaRestrictionType type, std::string const & value );

	std::string const & element_name() const override;
	void write_definition( int indentation, std::ostream & os ) const override;
	void prepare_for_output( XMLSchemaDefinition & xsd ) const override;

private:
	std::string name_;
	XMLSchemaType base_type_;
	std::list< std::pair< XMLSchemaRestrictionType, std::string > > restrictions_;
};

/// @brief This abstract class is meant to represent either an XMLSChemaElement
/// or an XMLSchemaModelGroup so that the interchangable set of these objects
/// in a ModelGroup can be represented.  I may be misusing the term "particle"
/// in the way that it is meant within XML Schema -- so the mapping of this term
/// to the term used in XML Schema is probably imperfect.
class XMLSchemaParticle : public XMLSchemaTopLevelElement {
public:
	XMLSchemaParticle();
	~XMLSchemaParticle() override;

	XMLSchemaParticle & min_occurs( int setting );
	XMLSchemaParticle & max_occurs( int setting );

	int min_occurs() const;
	int max_occurs() const;

private:
	int min_occurs_;
	int max_occurs_;

};

/// @brief The ModelGroup covers four XML Schema model groups: xs:sequence,
/// xs:choice, xs:all, and xs:group.  A ModelGroup may contain any number of
/// XMLSchemaParticles, interchanging between XMLElements and XMLModelGroups,
/// BUT there are some fairly heavy restrictions on the xs:all model group,
/// marginally enforced by this class and by XMLSchemaComplexType.
/// This class is not exactly a top-level element, in that only xs:group is
/// allowed to appear at the top level, but it is definitely worthwhile for
/// this class to implement the write_definition funciton.
class XMLSchemaModelGroup : public  XMLSchemaParticle {
public:
	XMLSchemaModelGroup();
	~XMLSchemaModelGroup() override;
	XMLSchemaModelGroup( XMLSchemaModelGroupType type );
	XMLSchemaModelGroup( XMLSchemaModelGroupType type, std::list< XMLSchemaParticleCOP > const & particles );
	XMLSchemaModelGroup( std::string const & group_name );

	XMLSchemaModelGroup & type( XMLSchemaModelGroupType type );
	XMLSchemaModelGroup & append_particle( XMLSchemaParticleCOP particle );
	XMLSchemaModelGroup & append_particles( std::list< XMLSchemaParticleCOP > const & particles );
	XMLSchemaModelGroup & group_name( std::string const & name );

	std::string const & element_name() const override;
	void write_definition( int indentation, std::ostream & os ) const override;
	void prepare_for_output( XMLSchemaDefinition & xsd ) const override;

private:
	void validate_content() const;
	bool is_group_holding_all() const;

private:
	XMLSchemaModelGroupType type_;
	std::list< XMLSchemaParticleCOP > particles_;
	std::string group_name_; // if this is an xs:group, what is the name for it? Any group with no particles is treated as a group reference
};

/// @brief class %XMLSchemaComplexType represents the definition of the type for an element -- that is,
/// the structure of a set of elements with the same name. If an XMLSchemaElement is analogous to a utility::tag::Tag,
/// an %XMLSchemaComplexType is analogous to the wiki page describing the valid format for an instance of that Tag.
///
/// @details %XMLSchemaComplexTypes can be named or unnamed: if they are named, then the are inserted into the
/// XMLSchemaDefinition on their own, and can be referred to by name in other XMLSchemaElements; if they are
/// unnamed, then they are given as an in-line definition for an XMLSchemaElement.  So a name is not required
/// for an %XMLSchemaComplexType.  That said: the recommended structure for defining an XML Schema is that
/// complex types should be named and to live on their own -- elements should refer to these complex types by name.
///
/// (There are some who recommend that instead elements be defined at "global scope" -- i.e. directly beneath the
/// xs:schema tag -- and that complex types be unnamed.  The principle disadvantage of that idea is that you cannot
/// have two elements in different contexts with the same name but different structures. This problem of "name
/// collisions" is addressed in XML Schema through the use of "namespacing," which is horrendous, so we won't use
/// it. By giving complex types names, we can avoid name collision by just giving the complex types different
/// names, which has no effect on the apperance of the XML file itself.)
///
/// Complex types must be one of the five options in the XMLSchemaComplexTypeType.  They are either:
/// * empty (xsctt_empty) meaning that they do not contain any sub-elements, though they may contain attributes.
///   Most of the complex types in Rosetta will be declared as empty.
/// * sequence (xsctt_sequence) meaning that they may contain attributes, and that they expect every
///   one of the listed sub-elements to appear inside of them in the order defined.
/// * all (xsctt_all) meaning that they may contain attributes, and that they expect every one of the
///   listed sub-elements to appear inside of them, but that these sub-elements may appear in any order
/// * choice (xsctt_choice) meaning that they may contain attributes, and that they expect exactly one
///   of the listed sub-elements to appear, or
/// * group (xsctt_group) in which case, attributes are not allowed, and a list of subelements are given.
///   The group category is perhaps the most confusing: a complex type with a group type will contain a
///   sub-complex-type with either a sequence, all, or choice designation.  A group can be written out
///   by name, and then referred to within one or more complex types.  Interestingly, the way to define
///   flexible, recursive data types is to use a group.  For example, if the ParsedProtocolMover type can
///   contain any mover, including a ParsedProtocolMover, then to represent this in XML schema requires
///   defining a group for all movers -- and then saying that instances of this group can appear within
///   the ParsedProtocolMoverType.  Similarly, residue selectors can contain other residue selectors.
///   Look at ResidueSelectorFactory::define_residue_selector_xml_schema for an example of how
///   a group is defined.
///
/// A complex type may refer to a reference type, and then extend that type, e.g. by adding extra attributes.
///
/// In the case of a "group" complex type, the complex type needs to contain a subtype.  E.g. the group
/// complex-type defined by ResidueSelectorFactory holds a choice complex-type.
///
/// Complex types may contain subelements; these subelements can be added through a model group
/// (e.g. "xs:all") that itself contains elements and possibly other model groups.
///
/// Complex types may contain attributes; these attributes can be added through the add_attribute method.
///
/// Complex types may be repeated, e.g., in the sequence or group cases.  You may specify the minimum
/// and maximum number of times the complex type is to appaer.
class XMLSchemaComplexType : public XMLSchemaTopLevelElement {
public:
	XMLSchemaComplexType();
	XMLSchemaComplexType & name( std::string const & setting );
	XMLSchemaComplexType & description( std::string const & setting );
	XMLSchemaComplexType & set_model_group( XMLSchemaModelGroupCOP model_group );
	XMLSchemaComplexType & add_attribute( XMLSchemaAttribute attribute );
	XMLSchemaComplexType & add_attributes( AttributeList const & attributes );

	std::string const & element_name() const override;
	void write_definition( int indentation, std::ostream & os ) const override;
	void prepare_for_output( XMLSchemaDefinition & xsd ) const override;

private:
	std::string name_;
	std::string desc_;
	XMLSchemaModelGroupCOP model_group_; // 0 if this is an empty complex type
	std::list< XMLSchemaAttribute > attributes_;
};


/// @brief An XMLSchema element, e.g. the FixbbMover "tag" <FixbbMover name="fixbb" task_operations="ex1,ex2"/> or
/// the And element which contains a sub-element: <And><Chain id="A"/></And>.  An element can be defined with an
/// unnamed complex type in-line, or it can say that its type is that of a complex type defined elsewhere.
/// (In XML naming a tag is a single block beginning with "<" and ending with ">", so both "<And>" and "</And>"
/// are tags. The block between the opening and closing tags is called an element.  The utility::tag::Tag class
/// represents a complete Element, and not simply the opening or closing tag of an element.  This is certainly
/// confusing.
///
/// @details An element may either have its type definition given inline (i.e. as a sub-tag / sub-element of the
/// element declaration) as either a complex type or a restriction -- or it may have its type definition given by
/// reference to a (named) type defined elsewhere in the XML Schema.  In Rosetta's style of XML, the type of an
/// element is universally a complex type; if the element's type were given as either a primitive type, or a
/// restriction, then the data would have to appear between two tags, e.g.
///
/// <FirstName> Andrew </FirstName>
///
/// and in Rosetta's XML, we ignore everything that appears outside of a tag, and store all of the data in
/// attributes (aka options) inside of tags.  The elements in Rosetta either are "empty" or they contain
/// sub-elements: both of these are classified as complex types in XML Schema.
///
/// Elements may appear multiple times, and the number of times they must/may appear are controlled
/// by the min_occurs and max_occurs functions.  If you are setting no limit on max_occurs, use the
/// "xsminmax_unbounded" value of the XMLSchemaMinOccursMaxOccurs enumeration.
class XMLSchemaElement : public XMLSchemaParticle {
public:
	XMLSchemaElement();
	XMLSchemaElement & name( std::string const & setting );
	XMLSchemaElement & set_abstract();
	XMLSchemaElement & substitution_group( std::string const & setting );
	XMLSchemaElement & type_name( XMLSchemaType setting );
	XMLSchemaElement & reference_name( std::string const & setting );
	XMLSchemaElement & element_type_def( XMLSchemaComplexTypeOP setting );
	//XMLSchemaElement & restriction_type_def( XMLSchemaRestrictionOP setting );

	std::string const & element_name() const override;
	void write_definition( int indentation, std::ostream & os ) const override;
	void prepare_for_output( XMLSchemaDefinition & xsd ) const override;

private:

	XMLSchemaElementCategory category_;
	std::string name_;
	XMLSchemaType type_name_;                     // used if category_ is xs_element_is_type_reference
	XMLSchemaComplexTypeOP complex_type_;         // used if category_ is xs_element_is_complex_type_w_definition
	std::string substitution_group_;
};

/// @brief The %XMLSchemaDefinition class's purpose is to collect all of the elements that go into an XML Schema,
/// to filter out the elements that are repeated (e.g. a restriction such as the "int_cslist" given in
/// the description for XMLSchemaRestriction above may be reported twice to the %XMLSchemaDefinition by several
/// attributes that rely upon it), to detect non-identical duplicates that have the same name, and to
/// write out the elements that it has been handed into a single string.  The %XMLSchemaDefinition is intended
/// to be passed between static functions / non-class-member functions as a container for the XML Schema
/// representations that these functions define.  Such functions will always take an %XMLSchemaDefinition reference
/// as one of their input parameters. It is perfectly legitimate / recommended for one XML-schema-defining
/// function that relies on a complexType or restriction that it does not itself define to pass its input
/// %XMLSchemaDefinition to the function that does define that complexType or restriction.
class XMLSchemaDefinition : public utility::pointer::ReferenceCount
{
public:
	XMLSchemaDefinition();
	~XMLSchemaDefinition() override;
	void add_top_level_element( XMLSchemaTopLevelElement const & element );
	bool has_top_level_element( std::string const & element_name ) const;
	std::string full_definition() const;

private:
	void validate_new_top_level_element( std::string const & element_name, std::string const & definition );

	std::list< std::string > elements_in_order_;
	std::map< std::string, std::string > top_level_elements_;
};


/// @brief The XMLSchemaSimpleSubelementList class defines an interface that can be used by those wishing
/// to define a complexType that contains sub-elements. The structure of the XML Schema for the sub-elements
/// will determined by how this list is given to the XMLSchemaComplexTypeGenerator.
/// "simple" subelements are those which themselves contain no subelements (but may contain attributes).
/// Also allowed are subelements that refer to previously-defined complexTypes or those that refer to
/// previously defined xs:groups.
class XMLSchemaSimpleSubelementList : public utility::pointer::ReferenceCount
{
public:
	// A function that returns the name for an object given the name of another related object
	// e.g. the complex_type_name_for_task_op function returns the name for the xs:complexType
	// for a task operation by modifying the name of that task operation.
	typedef boost::function< std::string ( std::string const & ) > DerivedNameFunction;

	// A function that returns the name for a particular object, e.g.
	// ResidueSelectorFactory::residue_selector_xml_schema_group_name gives the name
	// for the xs:group listing all of the ResidueSelectors accessible through
	// the factory (and thus all the ResidueSelectors accessible through RosettaScripts).
	typedef boost::function< std::string () >                      NameFunction;

	struct ElementSummary {
		ElementSummary();

		enum { ct_simple, ct_ref, ct_group } element_type;
		std::string element_name;
		std::string ct_name;
		std::string description;
		AttributeList attributes;
		bool min_or_max_occurs_set;
		int min_occurs;
		int max_occurs;
	};

public:
	XMLSchemaSimpleSubelementList();
	~XMLSchemaSimpleSubelementList() override;
	XMLSchemaSimpleSubelementList( XMLSchemaSimpleSubelementList const & src );
	XMLSchemaSimpleSubelementList & operator = ( XMLSchemaSimpleSubelementList const & rhs );

public:

	// Functions that the users should interact with

	/// @brief the naming function is required if this subelement list is going to be repeatable
	XMLSchemaSimpleSubelementList & complex_type_naming_func( DerivedNameFunction const & naming_function );

	/// @brief Add a new subelement to the growing list of subelements, defined by its name and a (possibly empty) list of attributes
	/// This subelement may not itself contain other subelements, however.  The name for the complexType of this element will be
	/// derived from the function given in complex_type_naming_func, if provided, and otherwise, the complexType for this element
	/// will be listed within the element declaration and will be unnamed.
	XMLSchemaSimpleSubelementList & add_simple_subelement( std::string const & name, AttributeList const &, std::string const & description  );

	/// @brief Add a new subelement to the growing list of subelements, but where the minimum and maximum number of
	/// occurrences for this subelement has been set; note: usually, the min and max are set by the
	/// XMLSchemaComplexTypeGenerator, and so tension could arise between the generator and your settings.
	/// Consider this an advanced function.
	XMLSchemaSimpleSubelementList & add_simple_subelement(
		std::string const & name,
		AttributeList const &,
		std::string const & description,
		int min_occurs,
		int max_occurs = xsminmax_unspecified
	);

	/// @brief Add a new subelement to the growing list of subelements, but one whose complexType definition has been provided
	/// elsewhere and where the name of the complexType for this subelement is derived from the subelement's name using the
	/// provided function.  (Why pass a function to define the name for the complex type instead of just passing the
	/// result of the function? Because such an interface would also let you bypass a function entirely and let you pass
	/// in a raw string, and that would lead to brittle code. If we should want to change the complexType naming rule
	/// for a particular class of elements in the schema, we want to be able to change only a single function and
	/// have that change ripple outwards through all the parts of the schema that interacted with those elements.)
	XMLSchemaSimpleSubelementList & add_already_defined_subelement(
		std::string const & name,
		DerivedNameFunction const & ct_naming_function
	);

	/// @brief Add a new subelement to the growing list of subelements, but where the minimum and maximum number of
	/// occurrences for this subelement has been set; note: usually, the min and max are set by the
	/// XMLSchemaComplexTypeGenerator, and so tension could arise between the generator and your settings.
	/// Consider this an advanced function.
	XMLSchemaSimpleSubelementList & add_already_defined_subelement(
		std::string const & name,
		DerivedNameFunction const & ct_naming_function,
		int min_occurs,
		int max_occurs = xsminmax_unspecified
	);

	/// @brief Add a subelement with a different name than the complex type was created with so that to get the
	/// correct name for the complex type, the original name must also be provided. E.g. The PROTOCOLS subelement
	/// takes its structure from the ParsedProtocol mover. Instead of duplicating the complex type for that mover,
	/// instead, create an element whose type will be given by the protocols::moves::complex_type_name_for_mover
	/// function when given the "ParsedProtocol" element name.
	XMLSchemaSimpleSubelementList & add_already_defined_subelement_w_alt_element_name(
		std::string const & alt_element_name, // the new element that you want to add
		std::string const & reference_element_name, // the original name that when combined with the ct_naming_function gives the correct complex type name
		DerivedNameFunction const & ct_naming_function
	);

	/// @brief Add a subelement with a different name than the complex type was created with so that to get the
	/// correct name for the complex type, the original name must also be provided. E.g. The PROTOCOLS subelement
	/// takes its structure from the ParsedProtocol mover. Instead of duplicating the complex type for that mover,
	/// instead, create an element whose type will be given by the protocols::moves::complex_type_name_for_mover
	/// function when given the "ParsedProtocol" element name.
	XMLSchemaSimpleSubelementList & add_already_defined_subelement_w_alt_element_name(
		std::string const & alt_element_name, // the new element that you want to add
		std::string const & reference_element_name, // the original name that when combined with the ct_naming_function gives the correct complex type name
		DerivedNameFunction const & ct_naming_function,
		int min_occurs,
		int max_occurs = xsminmax_unspecified
	);

	/// @brief Add a new subelement to the growing list of subelements which refers to an already existing group
	/// whose name is given by the input function.
	XMLSchemaSimpleSubelementList & add_group_subelement(
		NameFunction const & group_name_function
	);

	/// @brief Add a new subelement to the growing list of subelements, but where the minimum and maximum number of
	/// occurrences for this subelement has been set; note: usually, the min and max are set by the
	/// XMLSchemaComplexTypeGenerator, and so tension could arise between the generator and your settings.
	/// Consider this an advanced function.
	XMLSchemaSimpleSubelementList & add_group_subelement(
		NameFunction const & group_name_function,
		int min_occurs,
		int max_occurs = xsminmax_unspecified
	);

	/// @brief Really only intended to be used by the XMLSchemaRepeatableCTNode
	XMLSchemaSimpleSubelementList &
	add_subelement(
		ElementSummary const & summary
	);

public:

	// Functions interacted with by the XMLSchemaComplexTypeGenerator only

	bool simple_element_naming_func_has_been_set() const;
	std::string complex_typename_for_element( std::string const & element_name ) const;
	DerivedNameFunction naming_func() const;
	std::list< ElementSummary > const & element_list() const;

public:
	// Static functions for creating ElementSummary objects

	static
	ElementSummary
	element_summary_as_simple_subelement(
		std::string const & name,
		AttributeList const & attributes,
		std::string const & description
	);

	static
	ElementSummary
	element_summary_as_simple_subelement(
		std::string const & name,
		AttributeList const & attributes,
		std::string const & description,
		int min_occurs,
		int max_occurs
	);

	static
	ElementSummary
	element_summary_as_already_defined_subelement(
		std::string const & name,
		DerivedNameFunction const & ct_naming_function
	);

	static
	ElementSummary
	element_summary_as_already_defined_subelement(
		std::string const & name,
		DerivedNameFunction const & ct_naming_function,
		int min_occurs,
		int max_occurs
	);

	static
	ElementSummary
	element_summary_as_already_defined_subelement_w_alt_element_name(
		std::string const & alt_element_name, // the new element that you want to add
		std::string const & reference_element_name, // the original name that when combined with the ct_naming_function gives the correct complex type name
		DerivedNameFunction const & ct_naming_function
	);

	static
	ElementSummary
	element_summary_as_already_defined_subelement_w_alt_element_name(
		std::string const & alt_element_name, // the new element that you want to add
		std::string const & reference_element_name, // the original name that when combined with the ct_naming_function gives the correct complex type name
		DerivedNameFunction const & ct_naming_function,
		int min_occurs,
		int max_occurs /*= xsminmax_unspecified*/
	);

	static
	ElementSummary
	element_summary_as_group_subelement(
		NameFunction const & group_name_function
	);

	static
	ElementSummary
	element_summary_as_group_subelement(
		NameFunction const & group_name_function,
		int min_occurs,
		int max_occurs
	);

private:
	DerivedNameFunction ct_naming_func_for_simple_subelements_;
	std::list< ElementSummary > elements_;

};

/// @brief The %XMLComplexTypeSchemaGenerator is used to define the schema for a complex type
/// as they typically occurr in Rosetta.
///
/// @details There are seven main categories of complexTypes that regularly appear in Rosetta's
/// XML Schemas, listed in order of most common to least common:
///
/// 1) complexTypes that contain only attributes, but no subelements
///
/// 2) complexTypes that contain subelements that can repeat and that need not appear at all
///    e.g.
///    <DsspDesignOperation blueprint="somefile">
///      <Strand aa="ACDEFGH"/>
///      <Helix aa="ACDEFG"/>
///      <all append="IKLMNP"/>
///      <Helix exclude="P"/>
///    </DsspDesignOperation>
///
///    (Perhaps a repreat presence of "Helix" could be avoided, but the way DsspDesignOperation
///    works currently, repeat appearances are perfectly legal and sensical).
///
/// 3) complexTypes that contain subelements that must appear once, where the subelements are
///    allowed to appear in any order,
///
/// 4) complexTypes that contain subelements that can can appear zero or one times, where
///    the subelements are allowed to appear in any order,
///
/// 5) complexTypes where one of a set of subelements must be chosen,
///
/// 6) complexTypes where zero or one of a set of subelements can be chosen,
///
/// 7) complexTypes where elements must appear in a particular order (which is required if you have
///    two "xs:group" model groups you want both of, as is the case for the OperateOnResidueSubset
///    task operation which can hold both a "residue_selector" group and a "res_lvl_task_op" group),
///    and
///
/// 8) complexTypes that have ordered sets of subelements where each set is either:
///        a) repeatable,
///        b) optional,
///        c) required,
///        d) pick-one-from-a-list, or
///        e) pick-one-or-none-from-a-list.
///    In these complexTypes, no elements may appear in multiple sets.
///
///    For example, it's ok if you want to say "the <Option> subtag is first and it's either there or
///    not, and then afterwards you can have as many <ScoreFunction> and <TaskOperation> tags as you
///    want in any order, and finally, you can put either a single <PackRotamersMover> tag or a
///    <RotamerTrialsMover> tag." This would be like saying that the first set has an <Option>
///    element, and that it is optional, that the second set has both a <ScoreFunction> element
///    and a <TaskOperation> element and the elements of this set are repeatable, and that finally,
///    a third set contains the <PackRotamersMover> element and the <RotamerTrialsMover>, and that
///    neither has to appear but only one of the two may appear, and only once.
///
///    However, you cannot say
///    "The <Option> subtag is first, and its either there or not, and then afterwards, you
///    can have as many <ScoreFunction> and <TaskOperation> tags as you want, and then afterwards
///    you may have another <Option> subtag," because then the <Option> tag would appear in
///    two sets (the first and the third).
///    (The reason for this is restriction is to avoid violating the unique particle attribution
///    rule of XML Schema: https://en.wikipedia.org/wiki/Unique_Particle_Attribution )
///
/// Clearly XML Schema can support more kinds of elements than the eight that are supported by
/// this class, but these are the ones that seem to appear the most frequently. All
/// functionality provided by this class can also be accomplished using the XMLSchemaComplexType
/// and XMLSchemaElement classes, but would require a deeper understanding of XML Schema. The
/// niche that this class fills is to make it easy for developers to provide XSDs for these
/// most-common case. Adding additional functionality to this class comes at the expense of
/// making it less clear how the class should be used.  Write me if you feel that there
/// categories that should be added.
///
/// A description must be provided for every complex type even if that description
/// is an empty string; do not forget to provide a description (your code will compile
/// but will not run).
class XMLSchemaComplexTypeGenerator : public utility::pointer::ReferenceCount
{
public:
	// A function that returns the name for an object given the name of another related object
	// e.g. the complex_type_name_for_task_op function returns the name for the xs:complexType
	// for a task operation by modifying the name of that task operation.
	typedef boost::function< std::string ( std::string const & ) > DerivedNameFunction;

	// A function that returns the name for a particular object, e.g.
	// ResidueSelectorFactory::residue_selector_xml_schema_group_name gives the name
	// for the xs:group listing all of the ResidueSelectors accessible through
	// the factory (and thus all the ResidueSelectors accessible through RosettaScripts).
	typedef boost::function< std::string () >                      NameFunction;


public:

	XMLSchemaComplexTypeGenerator();
	~XMLSchemaComplexTypeGenerator() override;
	XMLSchemaComplexTypeGenerator( XMLSchemaComplexTypeGenerator const & src ) = delete;
	XMLSchemaComplexTypeGenerator & operator = (  XMLSchemaComplexTypeGenerator const & rhs ) = delete;

	XMLSchemaComplexTypeGenerator & element_name( std::string const & );

	/// @brief Provide the description for this complex type -- this is a function that must be called, even
	/// if you are passing in an empty string for the description.
	XMLSchemaComplexTypeGenerator & description( std::string const & );
	XMLSchemaComplexTypeGenerator & complex_type_naming_func( DerivedNameFunction const & naming_function );
	XMLSchemaComplexTypeGenerator & add_attribute( XMLSchemaAttribute const & attribute );
	// @brief Add a required "name" attribute of type "xs_string." Note that most name attributes are optional
	XMLSchemaComplexTypeGenerator & add_required_name_attribute( std::string const & desc = "" );
	// @brief Add an optional "name" attribute of type "xs_string." Note that most name attributes are optional
	XMLSchemaComplexTypeGenerator & add_optional_name_attribute( std::string const & desc = "" );
	XMLSchemaComplexTypeGenerator & add_attributes( AttributeList const & attributes );

	/// @brief set subelements as repeatable (and optional); setting the sub-elements replaces any sub-elements that were
	/// previously set. These repeatable subelements are allowed to appear in any order from the point of view of the XML
	/// Schema, which is not to say that the order in which they appear cannot matter to the code reading these subelements.
	/// This function represents case 2 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_repeatable(
		XMLSchemaSimpleSubelementList const & subelements,
		int min_occurs = 0,
		int max_occurs = xsminmax_unbounded
	);

	/// @brief set subelements as single-appearance (and required); setting the sub-elements replaces any sub-elements that were previously set.
	/// These single-appearance subelements are allowed to appear in any order from the point of view of the XML Schema, which is
	/// not to say that the order in which they appear cannot matter to the code reading these subelements.
	/// Due to limitations of XML Schema, if you have a "group" subelement, then it must appear alone; it may not appear with
	/// any other group subelements, and it may not appear with other "regular" subelements.
	/// This function represents case 3 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_single_appearance_required( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set subelements as single-appearance (and optional); setting the sub-elements replaces any sub-elements that were previously set.
	/// These single-appearance subelements are allowed to appear in any order from the point of view of the XML Schema, which is
	/// not to say that the order in which they appear cannot matter to the code reading these subelements.
	/// This function represents case 4 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_single_appearance_optional( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set the set of subelements with the stipulation that exactly one must be chosen.
	/// This function represents case 5 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_pick_one( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set the set of subelements with the stipulation that one or none should be chosen.
	/// This function represents case 6 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_pick_one_or_none( XMLSchemaSimpleSubelementList const & subelements );

	/// @brief set subelements as single-appearance (and required) with a specified order;
	/// setting the sub-elements replaces any sub-elements that were previously set.
	/// This function represents case 7 in the list of behaviors above.
	XMLSchemaComplexTypeGenerator & set_subelements_single_appearance_required_and_ordered(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being repeatable.
	/// This function corresponds to case 8a in the list of behaviors above and calling it
	/// overrides any of the functions for cases 1-7.
	XMLSchemaComplexTypeGenerator & add_ordered_subelement_set_as_repeatable(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being optional.  There should be only a single element in the
	/// input subelement list.
	/// This function corresponds to case 8b in the list of behaviors above and calling it
	/// overrides any of the functions for cases 1-7.
	XMLSchemaComplexTypeGenerator & add_ordered_subelement_set_as_optional(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as being requried.  There should be only a single element in the
	/// input subelement list.
	/// This function corresponds to case 8c in the list of behaviors above and calling it
	/// overrides any of the functions for cases 1-7.
	XMLSchemaComplexTypeGenerator & add_ordered_subelement_set_as_required(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as "pick exactly one".  There should be more than one element in the
	/// input subelement list.
	/// This function corresponds to case 8d in the list of behaviors above and calling it
	/// overrides any of the functions for cases 1-7.
	XMLSchemaComplexTypeGenerator & add_ordered_subelement_set_as_pick_one(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Add a subelement list to a growing set of ordered subelements, where elements in
	/// this set are labeled as "pick one or none".  There should be more than one element in the
	/// input subelement list.
	/// This function corresponds to case 8e in the list of behaviors above and calling it
	/// overrides any of the functions for cases 1-7.
	XMLSchemaComplexTypeGenerator & add_ordered_subelement_set_as_pick_one_or_none(
		XMLSchemaSimpleSubelementList const & subelements
	);

	/// @brief Return what the internal settings of the instance have been set to;
	/// this lets one class give another class access to a generator, and then inspect it
	/// on a crude level when it gets it back.
	CTGenSubelementBehavior subelement_behavior() const;

	void write_complex_type_to_schema( XMLSchemaDefinition & xsd );


private:
	class XMLSchemaComplexTypeGeneratorImpl * pimpl_;

};

class XMLSchemaRepeatableCTNode : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< XMLSchemaRepeatableCTNode >
{
public:
	typedef XMLSchemaSimpleSubelementList::DerivedNameFunction DerivedNameFunction;
	typedef XMLSchemaSimpleSubelementList::NameFunction NameFunction;

public:
	XMLSchemaRepeatableCTNode();
	virtual ~XMLSchemaRepeatableCTNode();

	/// self pointers
	inline XMLSchemaRepeatableCTNodeCOP get_self_ptr() const      { return shared_from_this(); }
	inline XMLSchemaRepeatableCTNodeOP  get_self_ptr()            { return shared_from_this(); }
	inline XMLSchemaRepeatableCTNodeCAP get_self_weak_ptr() const { return XMLSchemaRepeatableCTNodeCAP( shared_from_this() ); }
	inline XMLSchemaRepeatableCTNodeAP  get_self_weak_ptr()       { return XMLSchemaRepeatableCTNodeAP( shared_from_this() ); }

	XMLSchemaRepeatableCTNode &
	set_element_w_attributes(
		std::string const & name,
		AttributeList const & atts,
		std::string const & description
	);

	//XMLSchemaRepeatableCTNode &
	//set_element_w_attributes(
	// std::string const & name,
	// AttributeList const & atts,
	// std::string const & description,
	// int min_occurs,
	// int max_occurs = xsminmax_unspecified
	//);

	XMLSchemaRepeatableCTNode &
	set_already_defined_element(
		std::string const & name,
		DerivedNameFunction const & naming_func
	);

	//XMLSchemaRepeatableCTNode &
	//set_already_defined_element(
	// std::string const & name,
	// DerivedNameFunction const & naming_func,
	// int min_occurs,
	// int max_occurs = xsminmax_unspecified
	//);

	XMLSchemaRepeatableCTNode &
	set_already_defined_element_w_alt_name(
		std::string const & name,
		std::string const & reference_element_name,
		DerivedNameFunction const & naming_func
	);

	//XMLSchemaRepeatableCTNode &
	//set_already_defined_element_w_alt_name(
	// std::string const & name,
	// DerivedNameFunction const & naming_func,
	// int min_occurs,
	// int max_occurs = xsminmax_unspecified
	//);

	XMLSchemaRepeatableCTNode &
	set_group_subelement(
		NameFunction const & group_name_function
	);

	//XMLSchemaRepeatableCTNode &
	//set_group_subelement(
	// NameFunction const & group_name_function,
	// int min_occurs,
	// int max_occurs = xsminmax_unspecified
	//);

	XMLSchemaRepeatableCTNode &
	set_kids_naming_func( DerivedNameFunction const & naming_func );

	XMLSchemaRepeatableCTNode &
	set_root_node_naming_func( DerivedNameFunction const & naming_func );

	XMLSchemaRepeatableCTNode &
	add_child( XMLSchemaRepeatableCTNodeOP child_element );

	void
	recursively_write_ct_to_schema( XMLSchemaDefinition & xsd );

private:
	std::list< XMLSchemaRepeatableCTNodeOP > children_;
	XMLSchemaRepeatableCTNodeAP parent_;
	DerivedNameFunction kids_ct_naming_func_;
	DerivedNameFunction my_naming_func_; // for the root node only!
	XMLSchemaSimpleSubelementList::ElementSummary element_;
};


void indent_w_spaces( int indentation, std::ostream & os );
std::string restriction_type_name( XMLSchemaRestrictionType type );
std::ostream & operator << ( std::ostream & os, XMLSchemaRestrictionType type );
std::string xs_model_group_name( XMLSchemaModelGroupType xsmgt );
std::ostream & operator << ( std::ostream & os, XMLSchemaModelGroupType type );


/// @brief Convenience function for defining an inclusive range restriction
XMLSchemaRestriction
integer_range_restriction( std::string const & name, int lower_inclusive, int upper_inclusive );

/// @brief This function creates an attribute named "name" of type xs_string that is optional;
/// naming is not always required -- it is not even mostly required. This is probably the
/// function you need; very few classes actually need to have a name -- the name is often
/// only used by the function reading in the name, but there are reasonable times when a Mover,
/// e.g. could go nameless -- for instance, in the fixbb_jd3 application, a PackRotamersMover
/// can be given as a subtag of the <Job> tag. In this case, the PackRotamersMover doesn't need
/// to be given a name.
XMLSchemaAttribute
optional_name_attribute( std::string const & description = "" );

/// @brief This function creates an attribute named "name" of type xs_string that is required;
/// naming is not always required -- it is not even mostly required. This function is probably
/// not what you need. Use with care. See comments for the "optional_name_attribute" function
/// above.
XMLSchemaAttribute
required_name_attribute( std::string const & description = "" );

/// @brief append an attribute of "name" with type "xs:string" along with the other attributes
/// in the input attribute list to the input complex type.
void
append_name_and_attributes_to_complex_type(
	AttributeList const & attributes,
	XMLSchemaComplexType & type_definition
);

/// @brief append a required attribute of "name" with type "xs:string" along with the
/// other attributes in the input attribute list to the input complex type.
void
append_required_name_and_attributes_to_complex_type(
	AttributeList const & attributes,
	XMLSchemaComplexType & type_definition
);

/// @brief check whether an attribute with name attname already exists in AttributeList attlist to avoid collisions
bool
attribute_w_name_in_attribute_list(std::string const& attname, AttributeList const& attlist);

}  // namespace tag
}  // namespace utility

#endif
