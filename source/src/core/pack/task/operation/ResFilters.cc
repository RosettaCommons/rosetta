// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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


namespace core {
namespace pack {
namespace task {
namespace operation {

static THREAD_LOCAL basic::Tracer TR( "core.pack.task.operation.ResFilters" );

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

	for ( utility::vector0< TagCOP >::const_iterator subtag( subtags.begin() ), end( subtags.end() );
			subtag != end;
			++subtag ) {
		std::string const type( (*subtag)->getName() );

		ResFilterFactory * res_filter_factory = ResFilterFactory::get_instance();
		if ( res_filter_factory && res_filter_factory->has_type( type ) ) {
			ResFilterOP filter = res_filter_factory->newResFilter( type );
			filter->parse_tag( *subtag );
			sub_filters_.push_back(filter);

			continue;
		}
	}
}

utility::tag::XMLComplexTypeSchemaGeneratorOP
ResFilterComposition::define_composition_schema( utility::tag::XMLSchemaDefinition & )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	subelements.add_group_subelement( & ResFilterFactory::res_filter_xml_schema_group_name );

	XMLComplexTypeSchemaGeneratorOP ct_gen( new XMLComplexTypeSchemaGenerator );
	ct_gen
		->complex_type_naming_func( & complex_type_name_for_res_filter )
		.set_subelements_repeatable( subelements, & ResFilterFactory::res_filter_xml_schema_group_name );

	return ct_gen;
}

AnyResFilter::AnyResFilter() :
	parent()
{}

AnyResFilter::AnyResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP AnyResFilter::clone() const { return ResFilterOP( new AnyResFilter( *this ) ); }

bool AnyResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( core::Size i = 1; i <= sub_filters_.size(); ++i ) {
		if ( (*(sub_filters_[i]))(pose, index) ) {
			return true;
		}
	}

	return false;
}

std::string AnyResFilter::keyname() { return "AnyResFilter"; }

void AnyResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::XMLComplexTypeSchemaGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() );
	ct_gen->write_complex_type_to_schema( xsd );
}

ResFilterOP AnyResFilterCreator::create_res_filter() const
{
	return ResFilterOP( new AnyResFilter() );
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

ResFilterOP AllResFilter::clone() const { return ResFilterOP( new AllResFilter( *this ) ); }

bool AllResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( core::Size i = 1; i <= sub_filters_.size(); ++i ) {
		if ( !(*(sub_filters_[i]))(pose, index) ) {
			return false;
		}
	}

	return true;
}

std::string AllResFilter::keyname() { return "AllResFilter"; }

void AllResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {using namespace utility::tag;
	utility::tag::XMLComplexTypeSchemaGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() );
	ct_gen->write_complex_type_to_schema( xsd );
}

ResFilterOP AllResFilterCreator::create_res_filter() const
{
	return ResFilterOP( new AllResFilter() );
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

ResFilterOP NoResFilter::clone() const { return ResFilterOP( new NoResFilter( *this ) ); }

bool NoResFilter::operator() ( Pose const & pose, Size index ) const
{
	for ( core::Size i = 1; i <= sub_filters_.size(); ++i ) {
		if ( (*(sub_filters_[i]))(pose, index) ) {
			return false;
		}
	}

	return true;
}

std::string NoResFilter::keyname() { return "NoResFilter"; }

void NoResFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::XMLComplexTypeSchemaGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( keyname() );
	ct_gen->write_complex_type_to_schema( xsd );
}

ResFilterOP NoResFilterCreator::create_res_filter() const
{
	return ResFilterOP( new NoResFilter() );
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

ResFilterOP ResidueTypeFilter::clone() const { return ResFilterOP( new ResidueTypeFilter( *this ) ); }

bool ResidueTypeFilter::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.total_residue() );
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
	activate_common_simple_type( xsd, "zero_or_one" );

	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "polar", "zero_or_one" ));
	attributes.push_back( XMLSchemaAttribute::required_attribute( "apolar", "zero_or_one" ));
	attributes.push_back( XMLSchemaAttribute::required_attribute( "aromatic", "zero_or_one" ));
	attributes.push_back( XMLSchemaAttribute::required_attribute( "charged", "zero_or_one" ));

	res_filter_schema_w_attributes( xsd, keyname(), attributes );
}

ResFilterOP ResidueTypeFilterCreator::create_res_filter() const
{
	return ResFilterOP( new ResidueTypeFilter() );
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
	runtime_assert( index > 0 && index <= pose.total_residue() );
	return pose.residue_type(index).has_property( property_ );
}

ResFilterOP ResidueHasProperty::clone() const { return ResFilterOP( new ResidueHasProperty( *this ) ); }

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
	attributes.push_back( utility::tag::XMLSchemaAttribute( "property", utility::tag::xs_string ));
	return attributes;
}


ResFilterOP
ResidueHasPropertyCreator::create_res_filter() const {
	return ResFilterOP( new ResidueHasProperty );
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

ResFilterOP ResidueLacksProperty::clone() const { return ResFilterOP( new ResidueLacksProperty( *this ) ); }

std::string ResidueLacksProperty::keyname() { return "ResidueLacksProperty"; }

void ResidueLacksProperty::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResidueLacksPropertyCreator::create_res_filter() const {
	return ResFilterOP( new ResidueLacksProperty );
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
	runtime_assert( index > 0 && index <= pose.total_residue() );
	return pose.pdb_info()->res_haslabel( index, property_ );
}

ResFilterOP ResiduePDBInfoHasLabel::clone() const { return ResFilterOP( new ResiduePDBInfoHasLabel( *this ) ); }

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
	attributes.push_back( utility::tag::XMLSchemaAttribute( "property", utility::tag::xs_string ));
	return attributes;
}


ResFilterOP
ResiduePDBInfoHasLabelCreator::create_res_filter() const {
	return ResFilterOP( new ResiduePDBInfoHasLabel );
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

ResFilterOP ResiduePDBInfoLacksLabel::clone() const { return ResFilterOP( new ResiduePDBInfoLacksLabel( *this ) ); }

std::string ResiduePDBInfoLacksLabel::keyname() { return "ResiduePDBInfoLacksLabel"; }

void ResiduePDBInfoLacksLabel::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResiduePDBInfoLacksLabelCreator::create_res_filter() const {
	return ResFilterOP( new ResiduePDBInfoLacksLabel );
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
	runtime_assert( index > 0 && index <= pose.total_residue() );
	return name3_set.count(pose.residue_type(index).name3()) != 0;
}
ResFilterOP ResidueName3Is::clone() const { return ResFilterOP( new ResidueName3Is( *this ) ); }

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
	attributes.push_back( utility::tag::XMLSchemaAttribute( "name3", utility::tag::xs_string ));
	return attributes;
}


ResFilterOP
ResidueName3IsCreator::create_res_filter() const {
	return ResFilterOP( new ResidueName3Is );
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

ResFilterOP ResidueName3Isnt::clone() const { return ResFilterOP( new ResidueName3Isnt( *this ) ); }

std::string ResidueName3Isnt::keyname() { return "ResidueName3Isnt"; }

void ResidueName3Isnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes() );
}

ResFilterOP
ResidueName3IsntCreator::create_res_filter() const {
	return ResFilterOP( new ResidueName3Isnt );
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
	runtime_assert( index > 0 && index <= pose.total_residue() );
	return std::find( indices_.begin(), indices_.end(), index ) != indices_.end();
}

ResFilterOP ResidueIndexIs::clone() const { return ResFilterOP( new ResidueIndexIs( *this ) ); }

void ResidueIndexIs::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("indices") ) {
		indices_.clear();
		std::string field( tag->getOption<std::string>("indices") );
		utility::vector1< std::string > values( utility::string_split( field, ',' ) );
		for ( utility::vector1< std::string >::const_iterator it( values.begin() ),
				end( values.end() ); it != end; ++it ) {
			std::istringstream ss( *it );
			Size index;
			ss >> index;
			indices_.push_back( index );
		}
	}

	TR << "ResidueIndex with indices:";
	for ( utility::vector1< Size >::const_iterator it( indices_.begin() ), end( indices_.end() );
			it != end; ++it ) {
		TR << " " << *it;
	}
	TR << std::endl;

}

std::string ResidueIndexIs::keyname() { return "ResidueIndexIs"; }

void ResidueIndexIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes( xsd ) );
}

utility::tag::AttributeList
ResidueIndexIs::get_xml_schema_attributes( utility::tag::XMLSchemaDefinition & xsd )
{
	utility::tag::activate_common_simple_type( xsd, "int_cslist" );

	utility::tag::AttributeList attributes;
	attributes.push_back( utility::tag::XMLSchemaAttribute( "indices", "int_cslist" ));
	return attributes;
}

// static utility::tag::AttributeList get_xml_schema_attributes();



ResFilterOP
ResidueIndexIsCreator::create_res_filter() const {
	return ResFilterOP( new ResidueIndexIs );
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

ResFilterOP ResidueIndexIsnt::clone() const { return ResFilterOP( new ResidueIndexIsnt( *this ) ); }

std::string ResidueIndexIsnt::keyname() { return "ResidueIndexIsnt"; }

void ResidueIndexIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes( xsd ) );
}

ResFilterOP
ResidueIndexIsntCreator::create_res_filter() const {
	return ResFilterOP( new ResidueIndexIsnt );
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

ResiduePDBIndexIs::ResiduePDBIndexIs( char chain, int pos )
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
	runtime_assert( index > 0 && index <= pose.total_residue() );
	if ( ! pose.pdb_info() ) {
		utility_exit_with_message("can't apply ResiduePDBIndexIs filter on pose without pdb info");
	}
	ChainPos cp( pose.pdb_info()->chain(index), pose.pdb_info()->number(index) );
	return std::find( indices_.begin(), indices_.end(), cp ) != indices_.end();
}


ResFilterOP ResiduePDBIndexIs::clone() const { return ResFilterOP( new ResiduePDBIndexIs( *this ) ); }

/// @brief the expected format for the 'indices' option is: indices=A.2,B.3,Z.-20
void ResiduePDBIndexIs::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("indices") ) {
		indices_.clear();
		std::string field( tag->getOption<std::string>("indices") );
		// split option value by comma
		utility::vector1< std::string > values( utility::string_split( field, ',' ) );
		for ( utility::vector1< std::string >::const_iterator it( values.begin() ),
				end( values.end() ); it != end; ++it ) {
			// split pdb chain.pos token by '.'
			utility::vector1< std::string > chainposstr( utility::string_split( *it, '.' ) );
			if ( chainposstr.size() != 2 ) utility_exit_with_message("can't parse pdb index " + *it);
			char chain( *chainposstr.front().begin() );
			std::istringstream ss( chainposstr.back() );
			int pdbnum;
			ss >> pdbnum;
			indices_.push_back( ChainPos(chain,pdbnum) );
		}
	}

	TR << "ResiduePDBIndex with indices:";
	for ( utility::vector1< ChainPos >::const_iterator it( indices_.begin() ), end( indices_.end() );
			it != end; ++it ) {
		TR << " " << it->chain_ << '.' << it->pos_;
	}
	TR << std::endl;

}

std::string ResiduePDBIndexIs::keyname() { return "ResiduePDBIndexIs"; }

void ResiduePDBIndexIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(),get_xml_schema_attributes( xsd ) );
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

	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "indices", "period_separated_pdb_index_cslist" ));
	return attributes;
}


ResFilterOP
ResiduePDBIndexIsCreator::create_res_filter() const {
	return ResFilterOP( new ResiduePDBIndexIs );
}

std::string ResiduePDBIndexIsCreator::keyname() const { return ResiduePDBIndexIs::keyname(); }

void ResiduePDBIndexIsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ResiduePDBIndexIs::provide_xml_schema( xsd );
}

// ResiduePDBIndexIsnt
ResiduePDBIndexIsnt::ResiduePDBIndexIsnt()
: parent()
{}

ResiduePDBIndexIsnt::ResiduePDBIndexIsnt( char chain, int pos )
: parent( chain, pos )
{}

ResiduePDBIndexIsnt::ResiduePDBIndexIsnt( utility::vector1< ChainPos > const & indices )
: parent( indices )
{}

bool ResiduePDBIndexIsnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ResiduePDBIndexIsnt::clone() const { return ResFilterOP( new ResiduePDBIndexIs( *this ) ); }

std::string ResiduePDBIndexIsnt::keyname() { return "ResiduePDBIndexIsnt"; }

void ResiduePDBIndexIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes( xsd ) );
}

ResFilterOP
ResiduePDBIndexIsntCreator::create_res_filter() const {
	return ResFilterOP( new ResiduePDBIndexIsnt );
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

ChainIs::ChainIs( char const & chain )
: parent(),
	chain_( chain )
{}

bool ChainIs::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.total_residue() );
	if ( ! pose.pdb_info() ) {
		utility_exit_with_message("Can't apply ChainIs filter on pose without pdb info");
	}
	return( (pose.pdb_info()->chain(index) == chain_) ? true : false );
}

ResFilterOP ChainIs::clone() const { return ResFilterOP( new ChainIs( *this ) ); }

void ChainIs::parse_tag( TagCOP tag )
{
	chain_ = tag->getOption<char>("chain", 'A');
}

std::string ChainIs::keyname() { return "ChainIs"; }

void ChainIs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), get_xml_schema_attributes( xsd ) );
}

utility::tag::AttributeList
ChainIs::get_xml_schema_attributes( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction chain_character;
	chain_character.name( "chain_character" );
	chain_character.base_type( xs_string );
	chain_character.add_restriction( xsr_pattern, "[a-zA-Z]" );

	xsd.add_top_level_element( chain_character );

	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute::required_attribute( "property", "chain_character" ));
	return attributes;
}


ResFilterOP
ChainIsCreator::create_res_filter() const {
	return ResFilterOP( new ChainIs );
}

std::string ChainIsCreator::keyname() const { return ChainIs::keyname(); }

void ChainIsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ChainIs::provide_xml_schema( xsd );
}

// ChainIsnt
ChainIsnt::ChainIsnt()
: parent()
{}

ChainIsnt::ChainIsnt( char const & chain )
: parent( chain )
{}

bool ChainIsnt::operator() ( Pose const & pose, Size index ) const
{
	return ( ! parent::operator()( pose, index ) );
}

ResFilterOP ChainIsnt::clone() const { return ResFilterOP( new ChainIsnt( *this ) ); }

std::string ChainIsnt::keyname() { return "ChainIsnt"; }

void ChainIsnt::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	res_filter_schema_w_attributes( xsd, keyname(), parent::get_xml_schema_attributes( xsd ) );
}

ResFilterOP
ChainIsntCreator::create_res_filter() const {
	return ResFilterOP( new ChainIsnt );
}

std::string ChainIsntCreator::keyname() const { return ChainIsnt::keyname(); }

void ChainIsntCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ChainIsnt::provide_xml_schema( xsd );
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
