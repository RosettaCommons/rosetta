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

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

static basic::Tracer TR("core.pack.task.operation.ResFilters");

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

void ResidueHasProperty::parse_tag( TagPtr tag )
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
// begin ResidueName3Is
ResidueName3Is::ResidueName3Is()
	: parent()
{}

ResidueName3Is::ResidueName3Is( std::string const & str )
	: parent(),
		name3_( str )
{}

bool ResidueName3Is::operator() ( Pose const & pose, Size index ) const
{
	runtime_assert( index > 0 && index <= pose.total_residue() );
	return pose.residue_type(index).name3() == name3_;
}

ResFilterOP
ResidueName3IsCreator::create_res_filter() const {
	return new ResidueName3Is;
}

ResFilterOP ResidueName3Is::clone() const { return new ResidueName3Is( *this ); }

void ResidueName3Is::parse_tag( TagPtr tag )
{
	if ( tag->hasOption("name3") ) name3_ = tag->getOption<std::string>("name3");
}

// begin ResidueName3Isnt
ResidueName3Isnt::ResidueName3Isnt()
	: parent()
{}

ResidueName3Isnt::ResidueName3Isnt( std::string const & str )
	: parent( str )
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

void ResidueIndexIs::parse_tag( TagPtr tag )
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
void ResiduePDBIndexIs::parse_tag( TagPtr tag )
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

void ChainIs::parse_tag( TagPtr tag )
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
