// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/HelixKinkFilter.cc
/// @brief
/// @details filter structures out, which have kink helices
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HelixKinkFilter.hh>
#include <protocols/fldsgn/filters/HelixKinkFilterCreator.hh>

// Package Headers
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>


#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#ifdef WIN32
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#endif


//// C++ headers
static basic::Tracer TR( "protocols.fldsgn.filters.HelixKinkFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
HelixKinkFilter::HelixKinkFilter():
	Filter( "HelixKink" ),
	bend_angle_( 20.0 ),
	secstruct_( "" ),
	string_resnums_( "" ),
	select_resnums_( false ),
	select_range_( false ),
	helix_start_( 1 ),
	helix_end_( 1 )
{}


/// @brief Set the helix bend angle.
void
HelixKinkFilter::set_bend_angle( Real const bend_angle )
{
	bend_angle_ = bend_angle;
}

/// @brief Set the secondary structure string.
void
HelixKinkFilter::set_secstruct( String const &secstruct )
{
	secstruct_ = secstruct;
}

/// @brief Set the resnums (comma-separated string to be evaluated)
/// @details This function sets string_resnums_ and select_resnums_.  select_resnums_ will be set to arguments'
/// value.  If an empty string is passed (default), select_resnums_ will be set to false.  Otherwise,
/// select_resnums_ will be set to true.
void
HelixKinkFilter::set_resnums( String const &resnums )
{
	if ( resnums.empty() ) {
		select_resnums_ = false;
	} else {
		TR << "Only filter helix that contains residues specified in resnums" << std::endl;
		select_resnums_ = true;
		string_resnums_ = resnums;
	}
}

/// @brief Set helix_start_, helix_end_, and select_range_
/// @details This function sets helix_start_, helix_end_, and sets select_range_ to true
void
HelixKinkFilter::set_helix_range( core::Size const helix_start, core::Size const helix_end )
{
	if ( helix_start >= helix_end )  {
		utility_exit_with_message("helix_start_ is greater than or equal to helix_end_");
	}
	TR << "Only filter helix that contains range specified by helix_start to helix_end" << std::endl;
	select_range_ = true;
	helix_start_ = helix_start;
	helix_end_ = helix_end;
}

/// @brief Set select_range_
void
HelixKinkFilter::set_select_range( bool const select_range )
{
	select_range_ = select_range;
}

/// @brief Get the helix bend angle.
core::Real
HelixKinkFilter::get_bend_angle() const
{
	return bend_angle_;
}

/// @brief Get the secondary structure string.
std::string
HelixKinkFilter::get_secstruct() const
{
	return secstruct_;
}

/// @brief Get the resnum string.
std::string
HelixKinkFilter::get_string_resnums() const
{
	return string_resnums_;
}

/// @brief Get select_resnums_
bool
HelixKinkFilter::get_select_resnums() const
{
	return select_resnums_;
}

/// @brief Get helix_start_
core::Size
HelixKinkFilter::get_helix_start() const
{
	return helix_start_;
}

/// @brief Get helix_end_
core::Size
HelixKinkFilter::get_helix_end() const
{
	return helix_end_;
}

/// @brief Get select_range_
bool
HelixKinkFilter::get_select_range() const
{
	return select_range_;
}


/// @brief
bool
HelixKinkFilter::apply( Pose const & pose ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::check_kink_helix;
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	// set SS_Info
	String secstruct( pose.secstruct() );
	if ( secstruct_ != "" ) {
		secstruct = secstruct_;
	}
	SS_Info2_OP  ss_info( new SS_Info2( pose, secstruct ) );
	Helices const & helices( ss_info->helices() );

	if ( helices.size() < 1 ) {
		TR << "There is no helix definition in pose. " << std::endl;
		return true;
	}

	if ( ! pose.energies().data().has( HBOND_SET ) ) {
		TR << " Pose does not have HBOND_SET. Checking hbonds will be skipped. " << std::endl;
	}


	//vector to help quickly identify if a helix contains residues of interest
	utility::vector1<bool> residues_to_check;
	if ( select_resnums_ ) {
		utility::vector1< core::Size > const res_set_vec (core::pose::get_resnum_list_ordered( string_resnums_, pose ));
		residues_to_check.resize(pose.size(),false);
		TR << "filter residues contain: ";
		for ( core::Size i_res_vec = 1; i_res_vec <= res_set_vec.size(); ++i_res_vec ) {
			residues_to_check[res_set_vec[ i_res_vec ]]=true;
			TR << res_set_vec[ i_res_vec ] << " ";
		}
		TR << std::endl;
	} else {
		residues_to_check.resize(pose.size(),true);
	}


	//This checks if there is broken helix in the range already
	if ( select_range_ ) {
		for ( core::Size ii=1; ii<=helices.size(); ++ii ) {
			for ( core::Size it=helices[ ii ]->begin(), ite=helices[ ii ]->end(); it <= ite; ++it ) {
				TR.Debug << "Marking helix that begins at residue " << helices[ ii ]->begin() << " and ends at " << helices[ ii ]->end() << std::endl;
				residues_to_check[it] = false;
			}
		}

		for ( core::Size i=helix_start_; i<=helix_end_; ++i ) {
			if ( residues_to_check[i] == true ) {
				TR << "In range " << helix_start_ <<"-"<<helix_end_ << " contains broken helix at " << i << " already! Skip Kink" << std::endl;
				return false;
			}
		}
	}

	// check kink
	for ( core::Size ii=1; ii<=helices.size(); ++ii ) {
		bool check = false;
		if ( select_resnums_ ) {
			for ( core::Size it=helices[ ii ]->begin(), ite=helices[ ii ]->end(); it != ite; ++it ) {
				if ( residues_to_check[it] ) {
					check = true;
					//TR << "Helix " << ii << ", " << helices[ ii ]->begin() << "-" << helices[ ii ]->end() << ", is considered" << std::endl;
					break;
				}
				//TR << "Helix " << ii << ", " << helices[ ii ]->begin() << "-" << helices[ ii ]->end() << ", NOT considered" << std::endl;
			}
		} else if ( select_range_ ) {
			//This will check if the select_range_ is not broken into multiple helix
			if ( helices[ ii ]->begin() <= helix_start_ && helices[ ii ]->end() >= helix_end_ ) {
				TR << "Helix " << ii << ", res " << helices[ ii ]->begin() << "-" << helices[ ii ]->end() << "contains specified range: " << helix_start_ <<"-" << helix_end_ << std::endl;
				check = true;
			}
		} else {
			check = true; // TL: we always want to check if no residue numbers are specified
		}

		if ( !check ) continue;

		TR << "Helix " << ii << ", res " << helices[ ii ]->begin() << "-" << helices[ ii ]->end() << ", ";
		// check helix bend
		if ( helices[ ii ]->bend() > bend_angle_ ) {
			TR << "is bended angle=" << helices[ ii ]->bend() << std::endl;
			return false;
		} else {
			TR << "is bended angle=" << helices[ ii ]->bend() << std::endl;
		}

		// check broken hydrogen within helix
		/// @brief check kink of helix, return number of loosen hydrogen
		if ( pose.energies().data().has( HBOND_SET ) ) {
			core::Size hbond_start = helices[ ii ]->begin();
			if ( (hbond_start > 1) &&
					(!pose.residue(hbond_start).is_lower_terminus()) &&
					(pose.chain(hbond_start) == pose.chain(hbond_start-1)) ) {
				--hbond_start;
			}
			core::Size hbond_end = helices[ ii ]->end();
			if ( hbond_end >= hbond_start + 5 ) {
				hbond_end = hbond_end - 5;
			} else {
				hbond_end = hbond_start;
			}
			core::Size broken_hbonds( check_kink_helix( pose, hbond_start, hbond_end ) );
			if ( broken_hbonds > 0 ) {
				TR << "is kinked, hbonds are broken. " << std::endl;
				return false;
			}
		}

		TR << "is OK." << std::endl;
	}

	// check broken hydrogen bond
	TR << " Filter success ! " << std::endl;

	return true;

}


/// @brief parse xml
void
HelixKinkFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &
)
{
	bend_angle_  = tag->getOption<Real>( "bend",  20.0 );
	// secondary strucuture info
	// blueprint, or secstruct can be provided, but not both
	if ( tag->hasOption( "blueprint" ) && tag->hasOption( "secstruct" ) ) {
		std::stringstream msg;
		msg << "HelixKinkFilter::parse_my_tag(): Both blueprint and secstruct cannot be specified at the same time." << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( !blueprint.empty() ) {
		protocols::parser::BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
	}
	String const secstruct = tag->getOption<String>( "secstruct", "" );
	if ( !secstruct.empty() ) secstruct_ = secstruct;

	// residues that need to contained in helix
	if ( tag->hasOption( "resnums" ) ) {
		set_resnums(tag->getOption< std::string >( "resnums" ));
	} else {
		select_resnums_=false;
	}

	// residues that need to contained in a single helix
	if ( tag->hasOption( "helix_start" ) && tag->hasOption( "helix_end" ) ) {
		set_helix_range( tag->getOption< core::Size >( "helix_start" ), tag->getOption< core::Size >( "helix_end" ) );
	} else {
		select_range_=false;
	}
}



std::string HelixKinkFilter::name() const {
	return class_name();
}

std::string HelixKinkFilter::class_name() {
	return "HelixKink";
}

void HelixKinkFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "bend", xsct_real, "XRW TO DO", "20.0" )
		//Blueprint and secstruct are mutually exclusive
		+ XMLSchemaAttribute( "blueprint", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "secstruct", xsct_dssp_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "resnums", xsct_refpose_enabled_residue_number_cslist, "XRW TO DO" )
		+ XMLSchemaAttribute( "helix_start", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "helix_end", xsct_non_negative_integer, "XRW TO DO" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string HelixKinkFilterCreator::keyname() const {
	return HelixKinkFilter::class_name();
}

protocols::filters::FilterOP
HelixKinkFilterCreator::create_filter() const {
	return utility::pointer::make_shared< HelixKinkFilter >();
}

void HelixKinkFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HelixKinkFilter::provide_xml_schema( xsd );
}


} // filters
} // fldsgn
} // protocols
