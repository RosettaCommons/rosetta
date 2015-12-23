// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/tag/XMLSchemaGeneration.fwd.hh
/// @brief  forward declaration of the classes used to define an XML Schema
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
/// attributes, ) and they can also contain sub-elements.  Attributes
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
///                 \--> XMLSchemaType       \--> XMLSchemaElement                    \--> Restriction --> Primative Type
/// \endverbatic
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

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <vector>
#include <map>
#include <iosfwd>
#include <sstream>

namespace utility {
namespace tag {


/// @brief class %XMLSchemaType represents the name of a defined type that can
/// be used to describe either an XMLElement or an XMLAttribute.  It may refer
/// to either a complex type or to a primative type or to a simple type.
class XMLSchemaType : public utility::pointer::ReferenceCount {
public:
	XMLSchemaType();
	XMLSchemaType( XMLSchemaDataType setting );
	XMLSchemaType( std::string const & custom_type );
	XMLSchemaType( char const * custom_type );
	void type( XMLSchemaDataType setting );
	void custom_type_name( std::string const & setting );

	std::string type_name() const;

private:
	XMLSchemaDataType type_;
	std::string custom_type_name_;
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
class XMLSchemaAttribute : public utility::pointer::ReferenceCount {
public:
	XMLSchemaAttribute();
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, bool is_required=false );
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, std::string const & default_value );
	XMLSchemaAttribute( std::string const & name, XMLSchemaType type, char const * default_value );
	void name( std::string const & setting );
	void type( XMLSchemaType setting );
	void default_value( std::string const & setting );
	void is_required( bool setting );

	XMLSchemaType const & type() const;
	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaType type_;
	std::string default_value_;
	bool is_required_;
};

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

class XMLSchemaRestriction : public utility::pointer::ReferenceCount {
public:
	XMLSchemaRestriction();
	void name( std::string const & setting );
	void base_type( XMLSchemaType setting );
	void add_restriction( XMLSchemaRestrictionType type, std::string const & value );

	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaType base_type_;
	std::list< std::pair< XMLSchemaRestrictionType, std::string > > restrictions_;
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
/// Complex types may contain subelements; these subelements can be added through the add_subelement method.
///
/// Complex types may contain attributes; these attributes can be added through the add_attribute method.
///
/// Complex types may be repeated, e.g., in the sequence or group cases.  You may specify the minimum
/// and maximum number of times the complex type is to appaer.
class XMLSchemaComplexType : public utility::pointer::ReferenceCount {
public:
	XMLSchemaComplexType();
	void name( std::string const & setting );
	void type( XMLSchemaComplexTypeType setting );
	void reference_type( std::string const & setting ); // used if this is simply a reference to an existing type
	void subtype( XMLSchemaComplexTypeOP setting ); // used iff type_ is xsctt_group and reference_type_ is the empty string.
	void add_subelement( XMLSchemaElementOP  subelement );
	void add_attribute( XMLSchemaAttribute attribute );
	void min_occurs( int setting );
	void max_occurs( int setting );

	void write_definition( int indentation, std::ostream & os ) const;

private:
	std::string name_;
	XMLSchemaComplexTypeType type_;
	std::string reference_type_; // used if this is simply a reference to an existing type
	XMLSchemaComplexTypeOP subtype_; // used iff type_ is xsctt_group and reference_type_ is the empty string.
	std::list< XMLSchemaElementOP > subelements_;
	std::list< XMLSchemaAttribute > attributes_;
	int min_occurs_;
	int max_occurs_;
	bool suppress_complex_type_output_;
};


/// @brief An XMLSchema element, e.g. the FixbbMover "tag" <FixbbMover name="fixbb" task_operations="ex1,ex2"/> or
/// the And element which contains a sub-element: <And><Chain id="A"/></And>.  An element can be defined with an
/// unnamed complex type in-line, or it can say that its type is that of a complex type defined elsewhere.
/// (In XML naming a tag is a single block beginning with "<" and ending with ">", so both "<And>" and "</And>"
/// are tags.  The block between the opening and closing tags is called an element.  The utility::tag::Tag class
/// represents a complete Element, and not simply the opening or closing tag of an element.  This could easily
/// be confusing.
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
/// Elements which refer to existing complex types take two forms:
/// references where the type is listed by name right in the element tag
/// (e.g. <xs:element name="MinMover" type="MinMoverType/>
/// and elements whose type is actually a group and that are embedded in a complexType.
/// This second form is how recursive type definitions can be implemented in XML Schema.
/// For example (see the code snipet below), the "And" residue selector can take a list of
/// of ResidueSelectors that it will use as sub-selectors, including additional And selectors.
/// The way this is implemented in XML Schema is to define a group that lists every one of the
/// residue selector types, and then to say that the AndType is a "choice" complex type,
/// and the choice of elements is really one element that it may choose two or more times.
/// The critical thing to note here is that the sub-element beneath AndType's choice block isn't
/// an xs:element, but rather, an xs:group.
///
/// \verbatim
/// <xs:group name="residue_selector">
///  <xs:choice>
///   <xs:element name="And" type="AndType"/>
///   <xs:element name="Chain" type="ChainType"/>
///   <xs:element name="ClashBasedRepackShell" type="ClashBasedRepackShellType"/>
///   <xs:element name="DummyResidueSelector" type="DummyResidueSelectorType"/>
///   <xs:element name="Index" type="IndexType"/>
///   <xs:element name="InterfaceByVector" type="InterfaceByVectorType"/>
///   <xs:element name="JumpDownstream" type="JumpDownstreamType"/>
///   <xs:element name="JumpUpstream" type="JumpUpstreamType"/>
///   <xs:element name="Neighborhood" type="NeighborhoodType"/>
///   <xs:element name="Not" type="NotType"/>
///   <xs:element name="NumNeighbors" type="NumNeighborsType"/>
///   <xs:element name="Or" type="OrType"/>
///   <xs:element name="ResidueName" type="ResidueNameType"/>
///  </xs:choice>
/// </xs:group>
///
/// <xs:complexType name="AndType" mixed="true">
///  <xs:choice>
///   <xs:group ref="residue_selector" minOccurs="2" maxOccurs="unbounded"/>
///  </xs:choice>
///  <xs:attribute name="name" type="xs:string"/>
/// </xs:complexType>
/// \endverbatim
///
/// So the way this is implemented is that first an %XMLSchemaElement is created, and the "group_name"
/// method is set to "residue_selector".  This element's min_occurs is set to 2 and its max_occurs
/// is set to xsminmax_unbounded.  Then the XMLSchemaComplexType for "AndType" is set to be an
/// xsctt_choice, and its add_element method is called with the previosly created %XMLSchemaElement.
/// Again, the weirdness here is that the XMLSchemaElement does not report itself in the .xsd as an
/// "xs:element" but rather as an "xs:group".
///
/// Elements may appear multiple times, and the number of times they must/may appear are controlled
/// by the min_occurs and max_occurs functions.  If you are setting no limit on max_occurs, use the
/// "xsminmax_unbounded" value of the XMLSchemaMinOccursMaxOccurs enumeration.
class XMLSchemaElement : public utility::pointer::ReferenceCount {
public:
	XMLSchemaElement();
	void name( std::string const & setting );
	void type_name( XMLSchemaType setting );
	void group_name( std::string const & setting );
	void element_type_def( XMLSchemaComplexTypeOP setting );
	void restriction_type_def( XMLSchemaRestrictionOP setting );
	void min_occurs( int setting );
	void max_occurs( int setting );

	void write_definition( int indentation, std::ostream & os ) const;

private:

	std::string name_;
	XMLSchemaElementCategory category_;
	XMLSchemaType type_name_;                     // used if category_ is xs_element_is_type_reference
	std::string ref_name_;                        // used if category_ is xs_element_is_group_reference
	XMLSchemaComplexTypeOP complex_type_;         // used if category_ is xs_element_is_complex_type_w_definition
	XMLSchemaRestrictionOP restriction_type_def_; // used if category_ is xs_element_is_restriction_w_definition
	int min_occurs_;
	int max_occurs_;
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
	virtual ~XMLSchemaDefinition();
	void add_top_level_element( std::string const & element_name, std::string const & definition );
	bool has_top_level_element( std::string const & element_name ) const;
	std::string full_definition() const;

private:
	void validate_new_top_level_element( std::string const & element_name, std::string const & definition );

	std::list< std::string > elements_in_order_;
	std::map< std::string, std::string > top_level_elements_;
};


void indent_w_spaces( int indentation, std::ostream & os );
std::string restriction_type_name( XMLSchemaRestrictionType type );
std::ostream & operator << ( std::ostream & os, XMLSchemaRestrictionType type );
std::string xs_complex_type_type_name( XMLSchemaComplexTypeType xsctt );
std::ostream & operator << ( std::ostream & os, XMLSchemaComplexTypeType type );


}
}

#endif
