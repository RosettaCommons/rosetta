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
#include <core/pack/task/operation/ResFilterFactory.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
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

static basic::Tracer TR("core.pack.task.operation.ResFilters");

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

	if (sub_filters_.size() == 0)
	{
		TR.Debug << "ResFilterComposition without sub-filters defined: " << *tag << std::endl;
	}
}

void ResFilterComposition::parse_sub_filters_tag(TagCOP tag)
{
  utility::vector0< TagCOP > const & subtags( tag->getTags() );

  for (utility::vector0< TagCOP >::const_iterator subtag( subtags.begin() ), end( subtags.end() );
			subtag != end;
			++subtag )
	{
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

AnyResFilter::AnyResFilter() :
	parent()
{}

AnyResFilter::AnyResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP AnyResFilter::clone() const { return new AnyResFilter( *this ); }

bool AnyResFilter::operator() ( Pose const & pose, Size index ) const
{
	for (core::Size i = 1; i <= sub_filters_.size(); ++i)
	{
		if((*(sub_filters_[i]))(pose, index))
		{
			return true;
		}
	}

	return false;
}

ResFilterOP AnyResFilterCreator::create_res_filter() const
{
	return new AnyResFilter();
}

AllResFilter::AllResFilter() :
	parent()
{}

AllResFilter::AllResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP AllResFilter::clone() const { return new AllResFilter( *this ); }

bool AllResFilter::operator() ( Pose const & pose, Size index ) const
{
	for (core::Size i = 1; i <= sub_filters_.size(); ++i)
	{
		if(!(*(sub_filters_[i]))(pose, index))
		{
			return false;
		}
	}

	return true;
}

ResFilterOP AllResFilterCreator::create_res_filter() const
{
	return new AllResFilter();
}

NoResFilter::NoResFilter() :
	parent()
{}

NoResFilter::NoResFilter(utility::vector1<ResFilterCOP> const & sub_filters) :
	parent(sub_filters)
{}

ResFilterOP NoResFilter::clone() const { return new NoResFilter( *this ); }

bool NoResFilter::operator() ( Pose const & pose, Size index ) const
{
	for (core::Size i = 1; i <= sub_filters_.size(); ++i)
	{
		if((*(sub_filters_[i]))(pose, index))
		{
			return false;
		}
	}

	return true;
}

ResFilterOP NoResFilterCreator::create_res_filter() const
	{
	return new NoResFilter();
}

ResidueTypeFilter::ResidueTypeFilter() :
	parent(),
	polar_(false), apolar_(false), aromatic_(false), charged_(false)
{}

ResidueTypeFilter::ResidueTypeFilter(bool polar, bool apolar, bool aromatic, bool charged) :
	parent(),
	polar_(polar), apolar_(apolar), aromatic_(aromatic), charged_(charged)
{}

ResFilterOP ResidueTypeFilter::clone() const { return new ResidueTypeFilter( *this ); }

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

ResFilterOP ResidueTypeFilterCreator::create_res_filter() const
	{
	return new ResidueTypeFilter();
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

ResFilterOP
ResidueHasPropertyCreator::create_res_filter() const {
	return new ResidueHasProperty;
}

ResFilterOP ResidueHasProperty::clone() const { return new ResidueHasProperty( *this ); }

void ResidueHasProperty::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("property") ) property_ = tag->getOption<std::string>("property");
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

ResFilterOP
ResidueLacksPropertyCreator::create_res_filter() const {
  return new ResidueLacksProperty;
}

ResFilterOP ResidueLacksProperty::clone() const { return new ResidueLacksProperty( *this ); }

////////////////////////////////////////////////////////////////////////////////////////////////////
//begin ResiduePDB_InfoHasLabel
ResiduePDB_InfoHasLabel::ResiduePDB_InfoHasLabel()
  : parent()
{}

ResiduePDB_InfoHasLabel::ResiduePDB_InfoHasLabel( std::string const & str )
  : parent(),
    property_( str )
{}

bool ResiduePDB_InfoHasLabel::operator() ( Pose const & pose, Size index ) const
{
  runtime_assert( index > 0 && index <= pose.total_residue() );
  return pose.pdb_info()->res_haslabel( index, property_ );
}

ResFilterOP
ResiduePDB_InfoHasLabelCreator::create_res_filter() const {
  return new ResiduePDB_InfoHasLabel;
}

ResFilterOP ResiduePDB_InfoHasLabel::clone() const { return new ResiduePDB_InfoHasLabel( *this ); }

void ResiduePDB_InfoHasLabel::parse_tag( TagCOP tag )
{
  if ( tag->hasOption("property") ) property_ = tag->getOption<std::string>("property");
}
//end ResiduePDB_InfoHasLabel

//begin ResiduePDB_InfoLacksLabel
ResiduePDB_InfoLacksLabel::ResiduePDB_InfoLacksLabel()
  : parent()
{}

ResiduePDB_InfoLacksLabel::ResiduePDB_InfoLacksLabel( std::string const & str )
  : parent( str )
{}

bool ResiduePDB_InfoLacksLabel::operator() ( Pose const & pose, Size index ) const
{
  return ( ! parent::operator()( pose, index ) );
}

ResFilterOP
ResiduePDB_InfoLacksLabelCreator::create_res_filter() const {
  return new ResiduePDB_InfoLacksLabel;
}

ResFilterOP ResiduePDB_InfoLacksLabel::clone() const { return new ResiduePDB_InfoLacksLabel( *this ); }
//end ResiduePDB_InfoLacksLabel

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

ResFilterOP
ResidueName3IsCreator::create_res_filter() const {
	return new ResidueName3Is;
}

ResFilterOP ResidueName3Is::clone() const { return new ResidueName3Is( *this ); }

void ResidueName3Is::parse_tag( TagCOP tag )
{
	if ( tag->hasOption("name3") )
	{
		utility::vector1<std::string> names = utility::string_split(tag->getOption<std::string>("name3"), ',');
		name3_set.insert(names.begin(), names.end());
	}
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

ResFilterOP
ResidueName3IsntCreator::create_res_filter() const {
	return new ResidueName3Isnt;
}

ResFilterOP ResidueName3Isnt::clone() const { return new ResidueName3Isnt( *this ); }

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

ResFilterOP
ResidueIndexIsCreator::create_res_filter() const {
	return new ResidueIndexIs;
}

ResFilterOP ResidueIndexIs::clone() const { return new ResidueIndexIs( *this ); }

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

ResFilterOP
ResidueIndexIsntCreator::create_res_filter() const {
	return new ResidueIndexIsnt;
}

ResFilterOP ResidueIndexIsnt::clone() const { return new ResidueIndexIsnt( *this ); }

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

ResFilterOP
ResiduePDBIndexIsCreator::create_res_filter() const {
	return new ResiduePDBIndexIs;
}

ResFilterOP ResiduePDBIndexIs::clone() const { return new ResiduePDBIndexIs( *this ); }

///@brief the expected format for the 'indices' option is: indices=A.2,B.3,Z.-20
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

ResFilterOP
ResiduePDBIndexIsntCreator::create_res_filter() const {
	return new ResiduePDBIndexIsnt;
}

ResFilterOP ResiduePDBIndexIsnt::clone() const { return new ResiduePDBIndexIs( *this ); }

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

ResFilterOP
ChainIsCreator::create_res_filter() const {
	return new ChainIs;
}

ResFilterOP ChainIs::clone() const { return new ChainIs( *this ); }

void ChainIs::parse_tag( TagCOP tag )
{
	chain_ = tag->getOption<char>("chain", 'A');
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

ResFilterOP
ChainIsntCreator::create_res_filter() const {
	return new ChainIsnt;
}

ResFilterOP ChainIsnt::clone() const { return new ChainIsnt( *this ); }

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
