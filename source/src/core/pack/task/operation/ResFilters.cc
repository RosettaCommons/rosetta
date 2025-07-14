// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/ResFilters.hh
/// @brief  core-level (very general) classes that take a pose and a residue index, and return true or false
/// @author ashworth

// Unit Headers
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResFilterCreators.hh>

// Package headers
#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/task_op_schemas.hh>


// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <set>

#include <utility/tag/XMLSchemaGeneration.hh> // AUTO IWYU For operator+, XMLSchemaComplexTypeGenerator, XMLSchema...


namespace core {
namespace pack {
namespace task {
namespace operation {

static basic::Tracer TR( "core.pack.task.operation.ResFilters" );

ResFilterComposition::ResFilterComposition() :
	parent(),
	sub_filters_()
{}

ResFilterComposition::ResFilterComposition(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(),
	sub_filters_(sub_filters)
{}

void ResFilterComposition::parse_tag(TagCOP tag)
{
	parse_sub_filters_tag(tag);

	if ( sub_filters_.size() == 0 ) {
		TR.Debug << "ResFilterComposition without sub-filters defined: " << *tag << std::endl;
	}
}

void ResFilterComposition::parse_sub_filters_tag(TagCOP tag)
{
	utility::vector0< TagCOP > const & subtags( tag->getTags() );

	for ( auto const & subtag : subtags ) {
		std::string const & type( subtag->getName() );

		ResFilterFactory * res_filter_factory = ResFilterFactory::get_instance();
		if ( res_filter_factory && res_filter_factory->has_type( type ) ) {
			ResFilterOP filter = res_filter_factory->newResFilter( type );
			filter->parse_tag( subtag );
			sub_filters_.push_back(filter);

			continue;
		}
	}
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
ResFilterComposition::define_composition_schema( utility::tag::XMLSchemaDefinition & )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & ResFilterFactory::res_filter_xml_schema_group_name );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen
		->complex_type_naming_func( & complex_type_name_for_res_filter )
		.set_subelements_repeatable( subelements );

	return ct_gen;
}

AnyResFilter::AnyResFilter() :
	parent()
{}

AnyResFilter::AnyResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP AnyResFilter::clone() const { return utility::pointer::make_shared< AnyResFilter >( *this ); }

bool AnyResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( auto const & sub_filter : sub_filters_ ) {
		if ( (*sub_filter)(pose, index) ) {
			return true;
		}
	}

	return false;
}

std::string AnyResFilter::keyname() { return "AnyResFilter"; }

void AnyResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() )
		.description( "Compound filter. Any number of subfilters may be declared." )
		.write_complex_type_to_schema( xsd );
}

ResFilterOP AnyResFilterCreator::create_res_filter() const
{
	return utility::pointer::make_shared< AnyResFilter >();
}

std::string AnyResFilterCreator::keyname() const { return AnyResFilter::keyname(); }

void AnyResFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AnyResFilter::provide_xml_schema( xsd );
}

AllResFilter::AllResFilter() :
	parent()
{}

AllResFilter::AllResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP AllResFilter::clone() const { return utility::pointer::make_shared< AllResFilter >( *this ); }

bool AllResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( auto const & sub_filter : sub_filters_ ) {
		if ( !(*sub_filter)(pose, index) ) {
			return false;
		}
	}

	return true;
}

std::string AllResFilter::keyname() { return "AllResFilter"; }

void AllResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {using namespace utility::tag;
	utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() )
		.description( "Compound filter. Select all residues." )
		.write_complex_type_to_schema( xsd );
}

ResFilterOP AllResFilterCreator::create_res_filter() const
{
	return utility::pointer::make_shared< AllResFilter >();
}

std::string AllResFilterCreator::keyname() const { return AllResFilter::keyname(); }

void AllResFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AllResFilter::provide_xml_schema( xsd );
}

NoResFilter::NoResFilter() :
	parent()
{}

NoResFilter::NoResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP NoResFilter::clone() const { return utility::pointer::make_shared< NoResFilter >( *this ); }

bool NoResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( auto const & sub_filter : sub_filters_ ) {
		if ( (*sub_filter)(pose, index) ) {
			return false;
		}
	}

	return true;
}

std::string NoResFilter::keyname() { return "NoResFilter"; }

void NoResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() )
		.description( "Compound filter. Excludes residues specified in sub-tags from a task." )
		.write_complex_type_to_schema( xsd );
}

ResFilterOP NoResFilterCreator::create_res_filter() const
{
	return utility::pointer::make_shared< NoResFilter >();
}

std::string NoResFilterCreator::keyname() const { return NoResFilter::keyname(); }

void NoResFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NoResFilter::provide_xml_schema( xsd );
}

ResidueTypeFilter::ResidueTypeFilter() :
	parent(),
	polar_(false), apolar_(false), aromatic_(false), charged_(false)
{}

ResidueTypeFilter::ResidueTypeFilter(bool polar, bool apolar, bool aromatic, bool charged) :
	parent(),
	polar_(polar), apolar_(apolar), aromatic_(aromatic), charged_(charged)
{}

ResFilterOP ResidueTypeFilter::clone() const { return utility::pointer::make_shared< ResidueTypeFilter >( *this ); }

bool ResidueTypeFilter::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	core::conformation::Residue const & residue = pose.residue(index);

	return  (polar_ && residue.is_polar()) ||
		(apolar_ && residue.is_apolar()) ||
		(aromatic_ && residue.is_aromatic()) ||
		(charged_ && residue.is_charged());
}

void ResidueTypeFilter::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("polar") ) polar_ = tag->getOption<bool>("polar");
	if ( tag->hasOption("apolar") ) apolar_ = tag->getOption<bool>("apolar");
	if ( tag->hasOption("aromatic") ) aromatic_ = tag->getOption<bool>("aromatic");
	if ( tag->hasOption("charged") ) charged_ = tag->getOption<bool>("charged");
}

std::string ResidueTypeFilter::keyname() { return "ResidueType"; }

void ResidueTypeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "polar",    xsct_rosetta_bool , "Select polar residues." )
		+ XMLSchemaAttribute::required_attribute( "apolar",   xsct_rosetta_bool , "Select apolar residues." )
		+ XMLSchemaAttribute::required_attribute( "aromatic", xsct_rosetta_bool , "Select aromatic residues." )
		+ XMLSchemaAttribute::required_attribute( "charged",  xsct_rosetta_bool , "Select chared residues." );

	res_filter_schema_w_attributes( xsd, keyname(), attributes );
}

ResFilterOP ResidueTypeFilterCreator::create_res_filter() const
{
	return utility::pointer::make_shared< ResidueTypeFilter >();
}

std::string ResidueTypeFilterCreator::keyname() const { return ResidueTypeFilter::keyname(); }

void ResidueTypeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueTypeFilter::provide_xml_schema( xsd );
}



////////////////////////////////////////////////////////////////////////////////////////////////////
// begin ResidueHasProperty
ResidueHasProperty::ResidueHasProperty()
: parent()
{}

ResidueHasProperty::ResidueHasProperty( std::string const & str )
: parent(),
	property_( str )
{}

bool ResidueHasProperty::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	return pose.residue_type(index).has_property( property_ );
}

ResFilterOP ResidueHasProperty::clone() const { return utility::pointer::make_shared< ResidueHasProperty >( *this ); }

void ResidueHasProperty::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("property") ) property_ = tag->getOption<std::string>("property");
}

std::string ResidueHasProperty::keyname() { return "ResidueHasProperty"; }

void ResidueHasProperty::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes() );
}

utility::tag::AttributeList
ResidueHasProperty::get_xml_schema_attributes()
{
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute( "property", utility::tag::xs_string , "Select residues based on the give properties (DNA, PROTEIN, POLAR, CHARGED))" );
	return attributes;
}


ResFilterOP
ResidueHasPropertyCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueHasProperty >();
}

std::string ResidueHasPropertyCreator::keyname() const { return ResidueHasProperty::keyname(); }

void ResidueHasPropertyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueHasProperty::provide_xml_schema( xsd );
}

// begin ResidueLacksProperty
ResidueLacksProperty::ResidueLacksProperty()
: parent()
{}

ResidueLacksProperty::ResidueLacksProperty( std::string const & str )
: parent( str )
{}

bool ResidueLacksProperty::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResidueLacksProperty::clone() const { return utility::pointer::make_shared< ResidueLacksProperty >( *this ); }

std::string ResidueLacksProperty::keyname() { return "ResidueLacksProperty"; }

void ResidueLacksProperty::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResidueLacksPropertyCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueLacksProperty >();
}

std::string ResidueLacksPropertyCreator::keyname() const { return ResidueLacksProperty::keyname(); }

void ResidueLacksPropertyCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueLacksProperty::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//begin ResiduePDBInfoHasLabel
ResiduePDBInfoHasLabel::ResiduePDBInfoHasLabel()
: parent()
{}

ResiduePDBInfoHasLabel::ResiduePDBInfoHasLabel( std::string const & str )
: parent(),
	property_( str )
{}

bool ResiduePDBInfoHasLabel::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	return pose.pdb_info()->res_haslabel( index, property_ );
}

ResFilterOP ResiduePDBInfoHasLabel::clone() const { return utility::pointer::make_shared< ResiduePDBInfoHasLabel >( *this ); }

void ResiduePDBInfoHasLabel::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("property") ) property_ = tag->getOption<std::string>("property");
}

std::string ResiduePDBInfoHasLabel::keyname() { return "ResiduePDBInfoHasLabel"; }

void ResiduePDBInfoHasLabel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes() );
}

utility::tag::AttributeList
ResiduePDBInfoHasLabel::get_xml_schema_attributes()
{
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute( "property", utility::tag::xs_string , "Select residues based on ther pose::ResidueRecord.label" );
	return attributes;
}


ResFilterOP
ResiduePDBInfoHasLabelCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResiduePDBInfoHasLabel >();
}

std::string ResiduePDBInfoHasLabelCreator::keyname() const { return ResiduePDBInfoHasLabel::keyname(); }

void ResiduePDBInfoHasLabelCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResiduePDBInfoHasLabel::provide_xml_schema( xsd );
}

//end ResiduePDBInfoHasLabel

//begin ResiduePDBInfoLacksLabel
ResiduePDBInfoLacksLabel::ResiduePDBInfoLacksLabel()
: parent()
{}

ResiduePDBInfoLacksLabel::ResiduePDBInfoLacksLabel( std::string const & str )
: parent( str )
{}

bool ResiduePDBInfoLacksLabel::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResiduePDBInfoLacksLabel::clone() const { return utility::pointer::make_shared< ResiduePDBInfoLacksLabel >( *this ); }

std::string ResiduePDBInfoLacksLabel::keyname() { return "ResiduePDBInfoLacksLabel"; }

void ResiduePDBInfoLacksLabel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResiduePDBInfoLacksLabelCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResiduePDBInfoLacksLabel >();
}

std::string ResiduePDBInfoLacksLabelCreator::keyname() const { return ResiduePDBInfoLacksLabel::keyname(); }

void ResiduePDBInfoLacksLabelCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResiduePDBInfoLacksLabel::provide_xml_schema( xsd );
}

//end ResiduePDBInfoLacksLabel

////////////////////////////////////////////////////////////////////////////////////////////////////
// begin ResidueName3Is
ResidueName3Is::ResidueName3Is()
: parent(),
	name3_set()
{}

ResidueName3Is::ResidueName3Is( std::string const & str )
: parent(),
	name3_set()
{
	name3_set.insert(str);
}

ResidueName3Is::ResidueName3Is( std::set<std::string> const & strs )
: parent(),
	name3_set(strs)
{}

bool ResidueName3Is::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	return name3_set.count(pose.residue_type(index).name3()) != 0;
}
ResFilterOP ResidueName3Is::clone() const { return utility::pointer::make_shared< ResidueName3Is >( *this ); }

void ResidueName3Is::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("name3") ) {
		utility::vector1<std::string> names = utility::string_split(tag->getOption<std::string>("name3"), ',');
		name3_set.insert(names.begin(), names.end());
	}
}


std::string ResidueName3Is::keyname() { return "ResidueName3Is"; }

void ResidueName3Is::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes() );
}

utility::tag::AttributeList
ResidueName3Is::get_xml_schema_attributes()
{
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute( "name3", utility::tag::xs_string , "Select resides by three letter code, e.g ARG,LYS" );
	return attributes;
}


ResFilterOP
ResidueName3IsCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueName3Is >();
}

std::string ResidueName3IsCreator::keyname() const { return ResidueName3Is::keyname(); }

void ResidueName3IsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueName3Is::provide_xml_schema( xsd );
}

// begin ResidueName3Isnt
ResidueName3Isnt::ResidueName3Isnt()
: parent()
{}

ResidueName3Isnt::ResidueName3Isnt( std::string const & str )
: parent( str )
{}

ResidueName3Isnt::ResidueName3Isnt( std::set<std::string> const & strs )
: parent( strs )
{}

bool ResidueName3Isnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResidueName3Isnt::clone() const { return utility::pointer::make_shared< ResidueName3Isnt >( *this ); }

std::string ResidueName3Isnt::keyname() { return "ResidueName3Isnt"; }

void ResidueName3Isnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResidueName3IsntCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueName3Isnt >();
}

std::string ResidueName3IsntCreator::keyname() const { return ResidueName3Isnt::keyname(); }

void ResidueName3IsntCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueName3Isnt::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// begin ResidueIndexIs
ResidueIndexIs::ResidueIndexIs()
: parent()
{}

ResidueIndexIs::ResidueIndexIs( Size index )
: parent()
{
	indices_.push_back( index );
}

ResidueIndexIs::ResidueIndexIs( utility::vector1< Size > const & indices )
: parent(),
	indices_( indices )
{}

utility::vector1< Size > const &
ResidueIndexIs::indices() const { return indices_; }

bool ResidueIndexIs::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	return std::find( indices_.begin(), indices_.end(), index ) != indices_.end();
}

ResFilterOP ResidueIndexIs::clone() const { return utility::pointer::make_shared< ResidueIndexIs >( *this ); }

void ResidueIndexIs::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("indices") ) {
		indices_.clear();
		std::string field( tag->getOption<std::string>("indices") );
		utility::vector1< std::string > values( utility::string_split( field, ',' ) );
		for ( std::string const & value : values ) {
			std::istringstream ss( value );
			Size index;
			ss >> index;
			indices_.push_back( index );
		}
	}

	TR << "ResidueIndex with indices:";
	for ( Size const index : indices_ ) {
		TR << " " << index;
	}
	TR << std::endl;

}

std::string ResidueIndexIs::keyname() { return "ResidueIndexIs"; }

void ResidueIndexIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes() );
}

utility::tag::AttributeList
ResidueIndexIs::get_xml_schema_attributes()
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "indices", xsct_int_cslist , "Comma-separated list of residues to be selected." );
	return attributes;
}

// static utility::tag::AttributeList get_xml_schema_attributes();



ResFilterOP
ResidueIndexIsCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueIndexIs >();
}

std::string ResidueIndexIsCreator::keyname() const { return ResidueIndexIs::keyname(); }

void ResidueIndexIsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueIndexIs::provide_xml_schema( xsd );
}


// begin ResidueIndexIsnt
ResidueIndexIsnt::ResidueIndexIsnt()
: parent()
{}

ResidueIndexIsnt::ResidueIndexIsnt( Size index )
: parent( index )
{}

ResidueIndexIsnt::ResidueIndexIsnt( utility::vector1< Size > const & indices )
: parent( indices )
{}

bool ResidueIndexIsnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResidueIndexIsnt::clone() const { return utility::pointer::make_shared< ResidueIndexIsnt >( *this ); }

std::string ResidueIndexIsnt::keyname() { return "ResidueIndexIsnt"; }

void ResidueIndexIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResidueIndexIsntCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResidueIndexIsnt >();
}

std::string ResidueIndexIsntCreator::keyname() const { return ResidueIndexIsnt::keyname(); }

void ResidueIndexIsntCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResidueIndexIsnt::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// ResiduePDBIndexIs
ResiduePDBIndexIs::ResiduePDBIndexIs()
: parent()
{}

ResiduePDBIndexIs::ResiduePDBIndexIs( std::string const & chain, int pos )
: parent()
{
	indices_.push_back( ChainPos(chain,pos) );
}

ResiduePDBIndexIs::ResiduePDBIndexIs( utility::vector1< ChainPos > const & indices )
: parent(),
	indices_( indices )
{}

utility::vector1< ResiduePDBIndexIs::ChainPos > const &
ResiduePDBIndexIs::indices() const { return indices_; }

bool ResiduePDBIndexIs::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	if ( ! pose.pdb_info() ) {
		utility_exit_with_message("can't apply ResiduePDBIndexIs filter on pose without pdb info");
	}
	ChainPos cp( pose.pdb_info()->chain(index), pose.pdb_info()->number(index) );
	return std::find( indices_.begin(), indices_.end(), cp ) != indices_.end();
}


ResFilterOP ResiduePDBIndexIs::clone() const { return utility::pointer::make_shared< ResiduePDBIndexIs >( *this ); }

/// @brief the expected format for the 'indices' option is: indices=A.2,B.3,Z.-20
void ResiduePDBIndexIs::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("indices") ) {
		indices_.clear();
		std::string field( tag->getOption<std::string>("indices") );
		// split option value by comma
		utility::vector1< std::string > values( utility::string_split( field, ',' ) );
		for ( std::string const & value : values ) {
			// split pdb chain.pos token by '.'
			utility::vector1< std::string > chainposstr( utility::string_split( value, '.' ) );
			if ( chainposstr.size() != 2 ) utility_exit_with_message("can't parse pdb index " + value);
			std::string chain( chainposstr.front() );
			std::istringstream ss( chainposstr.back() );
			int pdbnum;
			ss >> pdbnum;
			indices_.push_back( ChainPos(chain,pdbnum) );
		}
	}

	TR << "ResiduePDBIndex with indices:";
	for ( auto const & chainpos : indices_ ) {
		TR << " " << chainpos.chain_ << '.' << chainpos.pos_;
	}
	TR << std::endl;

}

std::string ResiduePDBIndexIs::keyname() { return "ResiduePDBIndexIs"; }

void ResiduePDBIndexIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes( xsd ) );
}

utility::tag::AttributeList
ResiduePDBIndexIs::get_xml_schema_attributes( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction pdb_index_cslist;
	pdb_index_cslist.name( "period_separated_pdb_index_cslist" );
	pdb_index_cslist.base_type( xs_string );
	pdb_index_cslist.add_restriction( xsr_pattern, "[a-zA-Z]\\.-?[0-9]+(,[a-zA-Z]\\.-?[0-9]+)*" );

	xsd.add_top_level_element( pdb_index_cslist );
	std::cout << xsd.full_definition() << std::endl;

	AttributeList attributes;
	attributes + XMLSchemaAttribute( "indices", "period_separated_pdb_index_cslist" , "Select residues as \" chainID \".\" residue pdb number \", "
		"e.g. A.2 means chain A, residue at position 2." );
	return attributes;
}


ResFilterOP
ResiduePDBIndexIsCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResiduePDBIndexIs >();
}

std::string ResiduePDBIndexIsCreator::keyname() const { return ResiduePDBIndexIs::keyname(); }

void ResiduePDBIndexIsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResiduePDBIndexIs::provide_xml_schema( xsd );
}

// ResiduePDBIndexIsnt
ResiduePDBIndexIsnt::ResiduePDBIndexIsnt()
: parent()
{}

ResiduePDBIndexIsnt::ResiduePDBIndexIsnt( std::string const & chain, int pos )
: parent( chain, pos )
{}

ResiduePDBIndexIsnt::ResiduePDBIndexIsnt( utility::vector1< ChainPos > const & indices )
: parent( indices )
{}

bool ResiduePDBIndexIsnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResiduePDBIndexIsnt::clone() const { return utility::pointer::make_shared< ResiduePDBIndexIs >( *this ); }

std::string ResiduePDBIndexIsnt::keyname() { return "ResiduePDBIndexIsnt"; }

void ResiduePDBIndexIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes( xsd ) );
}

ResFilterOP
ResiduePDBIndexIsntCreator::create_res_filter() const {
	return utility::pointer::make_shared< ResiduePDBIndexIsnt >();
}

std::string ResiduePDBIndexIsntCreator::keyname() const { return ResiduePDBIndexIsnt::keyname(); }

void ResiduePDBIndexIsntCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResiduePDBIndexIsnt::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// ChainIs
ChainIs::ChainIs()
: parent()
{}

ChainIs::ChainIs( std::string const & chain )
: parent(),
	chain_( chain )
{}

bool ChainIs::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.size() );
	if ( ! pose.pdb_info() ) {
		utility_exit_with_message("Can't apply ChainIs filter on pose without pdb info");
	}
	return( (pose.pdb_info()->chain(index) == chain_) ? true : false );
}

ResFilterOP ChainIs::clone() const { return utility::pointer::make_shared< ChainIs >( *this ); }

void ChainIs::parse_tag( TagCOP tag )
{
	chain_ = tag->getOption<std::string>("chain", "A");
}

std::string ChainIs::keyname() { return "ChainIs"; }

void ChainIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes( xsd ) );
}

utility::tag::AttributeList
ChainIs::get_xml_schema_attributes( utility::tag::XMLSchemaDefinition & )
{
	using namespace utility::tag;

	AttributeList attributes;
	attributes + XMLSchemaAttribute::required_attribute( "chain", xs_string , "Select chain by chain ID." );
	return attributes;
}


ResFilterOP
ChainIsCreator::create_res_filter() const {
	return utility::pointer::make_shared< ChainIs >();
}

std::string ChainIsCreator::keyname() const { return ChainIs::keyname(); }

void ChainIsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ChainIs::provide_xml_schema( xsd );
}

// ChainIsnt
ChainIsnt::ChainIsnt()
: parent()
{}

ChainIsnt::ChainIsnt( std::string const & chain )
: parent( chain )
{}

bool ChainIsnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ChainIsnt::clone() const { return utility::pointer::make_shared< ChainIsnt >( *this ); }

std::string ChainIsnt::keyname() { return "ChainIsnt"; }

void ChainIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes( xsd ) );
}

ResFilterOP
ChainIsntCreator::create_res_filter() const {
	return utility::pointer::make_shared< ChainIsnt >();
}

std::string ChainIsntCreator::keyname() const { return ChainIsnt::keyname(); }

void ChainIsntCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ChainIsnt::provide_xml_schema( xsd );
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
