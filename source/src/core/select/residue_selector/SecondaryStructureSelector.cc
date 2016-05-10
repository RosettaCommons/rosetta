// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/SecondaryStructureSelector.cc
/// @brief  The SecondaryStructureSelector selects residues based on secondary structure
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit headers
#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.SecondaryStructureSelector" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace residue_selector {

SecondaryStructureSelector::SecondaryStructureSelector() :
	ResidueSelector(),
	pose_secstruct_( "" ),
	overlap_( 0 ),
	include_terminal_loops_( false ),
	always_use_dssp_( true ),
	selected_ss_()
{
	selected_ss_.clear();
}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP SecondaryStructureSelector::clone() const { return ResidueSelectorOP( new SecondaryStructureSelector(*this) ); }

SecondaryStructureSelector::SecondaryStructureSelector( std::string const & selected ) :
	ResidueSelector(),
	pose_secstruct_( "" ),
	overlap_( 0 ),
	include_terminal_loops_( false ),
	always_use_dssp_( true ),
	selected_ss_()
{
	if ( !check_ss( selected ) ) {
		throw utility::excn::EXCN_BadInput( "SecondaryStructureSelector: invalid ss to select: " + selected + " -- only H, E and L are allowed.\n" );
	}
	set_selected_ss( selected );
}

SecondaryStructureSelector::~SecondaryStructureSelector()
{
}

void
SecondaryStructureSelector::set_pose_secstruct( std::string const & ss )
{
	pose_secstruct_ = ss;
}

IntervalVec
subset_to_intervals( ResidueSubset const & subset )
{
	IntervalVec intervals;
	bool inrange = false;
	Interval current = std::make_pair( Size( 0 ), Size( 0 ) );
	for ( core::Size i=1, endi=subset.size(); i<=endi; ++i ) {
		if ( !inrange && subset[ i ] ) {
			current.first = i;
			inrange = true;
		}
		if ( inrange && !subset[ i ] ) {
			current.second = i - 1;
			inrange = false;
			intervals.push_back( current );
			current = std::make_pair( Size( 0 ), Size( 0 ) );
		}
	}
	return intervals;
}

bool is_poly_l( std::string const & ss )
{
	for ( std::string::const_iterator c = ss.begin(), endc = ss.end(); c != endc; ++c ) {
		if ( *c != 'L' ) {
			return false;
		}
	}
	return true;
}

ResidueSubset
SecondaryStructureSelector::apply( core::pose::Pose const & pose ) const
{
	// selected secondary structure can't be empty
	if ( selected_ss_.empty() ) {
		std::stringstream err;
		err << "SecondaryStructureSelector: no secondary structure types are selected -- you must specify one or more secondary structure types for selection." << std::endl;
		throw utility::excn::EXCN_BadInput( err.str() );
	}

	// first check pose secstruct, otherwise use dssp
	std::string ss = "";
	if ( pose_secstruct_.empty() ) {
		if ( (!always_use_dssp_) && check_ss( pose.secstruct() ) && !is_poly_l( pose.secstruct() ) ) {
			ss = pose.secstruct();
		} else {
			core::scoring::dssp::Dssp dssp( pose );
			ss = dssp.get_dssp_secstruct();
		}
	} else {
		ss = pose_secstruct_;
	}

	if ( !check_ss( ss ) ) {
		std::stringstream err;
		err << "SecondaryStructureSelector: the secondary structure string " << ss << " is invalid." << std::endl;
		throw utility::excn::EXCN_BadInput( err.str() );
	}

	TR << "Using secondary structure " << ss << std::endl;
	ResidueSubset matching_ss( pose.total_residue(), false );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( selected_ss_.find( ss[ i - 1 ] ) != selected_ss_.end() ) {
			TR.Debug << "Found ss match at position " << i << std::endl;
			matching_ss[ i ] = true;
		}
	}

	add_overlap( matching_ss, pose, ss );

	return matching_ss;
}

void SecondaryStructureSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & )
{
	if ( tag->hasOption( "overlap" ) ) {
		set_overlap( tag->getOption< core::Size >( "overlap" ) );
	}
	if ( tag->hasOption( "include_terminal_loops" ) ) {
		set_include_terminal_loops( tag->getOption< bool >( "include_terminal_loops" ) );
	}

	if ( tag->hasOption( "use_dssp" ) ) {
		set_always_use_dssp( tag->getOption< bool >( "use_dssp" ) );
	}

	if ( tag->hasOption( "pose_secstruct" ) ) {
		set_pose_secstruct( tag->getOption< std::string >( "pose_secstruct" ) );
	}

	if ( !tag->hasOption( "ss" ) ) {
		std::stringstream err;
		err << "SecondaryStructureSelector: the ss option is required to specify which secondary structure elements are being selected." << std::endl;
		err << "Usage: ss=\"([HEL])*\" (e.g. ss=\"L\" for loops only, ss=\"HE\" for helices and strands)" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( err.str() );
	}
	std::string const selected = tag->getOption< std::string >( "ss" );

	if ( !check_ss( selected ) ) {
		std::stringstream err;
		err << "SecondaryStructureSelector: an invalid character was specified in the ss option -- only H, E and L are allowed: " << selected << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( err.str() );
	}

	set_selected_ss( selected );
}

std::string
SecondaryStructureSelector::get_name() const
{
	return SecondaryStructureSelector::class_name();
}

std::string
SecondaryStructureSelector::class_name()
{
	return "SecondaryStructure";
}

void
SecondaryStructureSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "overlap",                xs_integer )
		+ XMLSchemaAttribute( "include_terminal_loops", xs_boolean )
		+ XMLSchemaAttribute( "pose_secstruct",         xs_string )
		+ XMLSchemaAttribute::required_attribute( "ss", xs_string );
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}

bool
SecondaryStructureSelector::check_ss( std::string const & ss ) const
{
	for ( std::string::const_iterator c = ss.begin(); c != ss.end(); ++c ) {
		if ( *c == 'L' ) {
			continue;
		}
		if ( *c == 'E' ) {
			continue;
		}
		if ( *c == 'H' ) {
			continue;
		}
		return false;
	}
	return true;
}

void
SecondaryStructureSelector::add_overlap(
	ResidueSubset & matching_ss,
	pose::Pose const & pose,
	std::string const & ss ) const
{
	IntervalVec intervals = subset_to_intervals( matching_ss );
	for ( IntervalVec::iterator i = intervals.begin(); i != intervals.end(); ++i ) {
		Size count = Size( 0 );
		while ( (count < overlap_) && !pose::pose_residue_is_terminal( pose, i->first ) ) {
			i->first -= 1;
			matching_ss[ i->first ] = true;
			count += 1;
		}
		count = Size( 0 );
		while ( (count < overlap_) && !pose::pose_residue_is_terminal( pose, i->second ) ) {
			i->second += 1;
			matching_ss[ i->second ] = true;
			count += 1;
		}
		TR.Debug << "Interval: " << i->first << " " << i->second << " " << ss[ i->first - 1 ] << std::endl;

		// get rid of one-residue terminal "loops"
		if ( include_terminal_loops_ ) {
			continue;
		}

		if ( ( i->first == i->second ) && ( ss[ i->first - 1 ] == 'L' )
				&& ( pose::pose_residue_is_terminal( pose, i->first ) ) ) {
			matching_ss[ i->first ] = false;
		}
		if ( ( i->first + 1 == i->second ) && ( ss[ i->first - 1 ] == 'L' )
				&& ( pose::pose_residue_is_terminal( pose, i->first ) ) ) {
			matching_ss[ i->first ] = false;
		}
		if ( ( i->first + 1 == i->second ) && ( ss[ i->second - 1 ] == 'L' )
				&& ( pose::pose_residue_is_terminal( pose, i->second ) ) ) {
			matching_ss[ i->second ] = false;
		}
	}
}

void
SecondaryStructureSelector::set_always_use_dssp( bool const always_use_dssp )
{
	always_use_dssp_ = always_use_dssp;
}

void
SecondaryStructureSelector::set_overlap( core::Size const overlapval )
{
	overlap_ = overlapval;
}

void
SecondaryStructureSelector::set_selected_ss( std::string const & selected )
{
	selected_ss_.clear();
	selected_ss_ = std::set< char >( selected.begin(), selected.end() );
}

/// @brief if true, one-residue terminal "loops" will be included (default=false)
void
SecondaryStructureSelector::set_include_terminal_loops( bool const inc_term )
{
	include_terminal_loops_ = inc_term;
}

ResidueSelectorOP
SecondaryStructureSelectorCreator::create_residue_selector() const
{
	return ResidueSelectorOP( new SecondaryStructureSelector );
}

std::string
SecondaryStructureSelectorCreator::keyname() const
{
	return SecondaryStructureSelector::class_name();
}

void
SecondaryStructureSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SecondaryStructureSelector::provide_xml_schema( xsd );
}

} //namespace residue_selector
} //namespace select
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::SecondaryStructureSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( pose_secstruct_ ) ); // std::string
	arc( CEREAL_NVP( overlap_ ) ); // core::Size
	arc( CEREAL_NVP( include_terminal_loops_ ) ); // _Bool
	arc( CEREAL_NVP( always_use_dssp_ ) ); // _Bool
	arc( CEREAL_NVP( selected_ss_ ) ); // std::set<char>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::SecondaryStructureSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( pose_secstruct_ ); // std::string
	arc( overlap_ ); // core::Size
	arc( include_terminal_loops_ ); // _Bool
	arc( always_use_dssp_ ); // _Bool
	arc( selected_ss_ ); // std::set<char>
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::SecondaryStructureSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::SecondaryStructureSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_SecondaryStructureSelector )
#endif // SERIALIZATION
