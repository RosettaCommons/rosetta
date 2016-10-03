// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/HSSTripletFilter.cc
/// @brief filter structures by hsstriplets angle and distance
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HSSTripletFilter.hh>
#include <protocols/fldsgn/filters/HSSTripletFilterCreator.hh>

// Package Headers
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

//// C++ headers
static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.filters.HSSTripletFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

/// @brief default constructor
HSSTripletFilter::HSSTripletFilter():
	Filter( "HSSTriplet" ),
	hss3set_(),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -12.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	ignore_helix_direction_( false ),
	use_dssp_( true )
{}

/// @brief constructor with arguments
HSSTripletFilter::HSSTripletFilter( HSSTriplets const & hss3s ):
	Filter( "HSSTriplet" ),
	hss3set_( HSSTripletSetOP( new HSSTripletSet( hss3s ) ) ),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -12.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	ignore_helix_direction_( false ),
	use_dssp_( true )
{}

/// @brief constructor with arguments
HSSTripletFilter::HSSTripletFilter( String const & hss3s ):
	Filter( "HSSTriplet" ),
	hss3set_( HSSTripletSetOP( new HSSTripletSet( hss3s ) ) ),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -12.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	ignore_helix_direction_( false ),
	use_dssp_( true )
{}

/// @brief copy constructor -- required because we clone the hss3set_ pointer
HSSTripletFilter::HSSTripletFilter( HSSTripletFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	hss3set_(),
	secstruct_( rval.secstruct_ ),
	filter_min_dist_( rval.filter_min_dist_ ),
	filter_max_dist_( rval.filter_max_dist_ ),
	filter_min_angle_( rval.filter_min_angle_ ),
	filter_max_angle_( rval.filter_max_angle_ ),
	output_id_( rval.output_id_ ),
	output_type_( rval.output_type_ ),
	ignore_helix_direction_( rval.ignore_helix_direction_ ),
	use_dssp_( rval.use_dssp_ )
{
	if ( rval.hss3set_ ) {
		hss3set_ = HSSTripletSetOP( new HSSTripletSet( *rval.hss3set_ ) );
	}
}


/// @brief set filtered HSSTriplets from string
void HSSTripletFilter::add_hsstriplets( String const & hss3s )
{
	HSSTripletSet const hss3set( hss3s );
	add_hsstriplets( hss3set.hss_triplets() );
}

// @brief set filtered HSSTriplets
void HSSTripletFilter::add_hsstriplets( HSSTriplets const & hss3s )
{
	if ( !hss3set_ ) hss3set_ = HSSTripletSetOP( new HSSTripletSet );
	hss3set_->add_hsstriplets( hss3s );
}

// @brief set secondary structure
void HSSTripletFilter::secstruct( String const & ss )
{
	secstruct_ = ss;
}

// @brief minimum distance for filtering
void
HSSTripletFilter::filter_min_dist( Real const r )
{
	filter_min_dist_ = r;
}

/// @brief maximum distance for filtering
void
HSSTripletFilter::filter_max_dist( Real const r )
{
	filter_max_dist_ = r;
}

/// @brief miniimum angle for filtering
void
HSSTripletFilter::filter_min_angle( Real const r )
{
	filter_min_angle_ = r; } /// @brief maximum angle for filtering
void
HSSTripletFilter::filter_max_angle( Real const r )
{
	filter_max_angle_ = r;
}

/// @brief maximum angle for filtering
void
HSSTripletFilter::output_id( Size const i )
{
	output_id_ = i;
}

/// @brief
void
HSSTripletFilter::output_type( String const & s )
{
	output_type_ = s;
}

/// @brief
HSSTripletFilter::Real
HSSTripletFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}

// @brief returns true if the given pose passes the filter, false otherwise.
bool
HSSTripletFilter::apply( Pose const & pose ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::SS_Info2_COP;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::HSSTripletOP;
	using protocols::fldsgn::topology::HSSTriplets;

	// get secondary structure
	std::string const secstruct = get_secstruct( pose );
	runtime_assert( secstruct.length() == pose.size() );

	// get HSSTriplets to filter
	HSSTriplets const hss3s = get_hss3s( pose );

	// set SS_Info
	SS_Info2_COP const ss_info( SS_Info2_OP( new SS_Info2( pose, secstruct ) ) );

	// check conformation hsstriplets
	bool filter( true );
	Size current_id( 0 );
	for ( HSSTriplets::const_iterator it=hss3s.begin(); it!=hss3s.end(); ++it ) {
		current_id ++;
		HSSTriplet hss = **it;

		if ( !check_elements( hss, *ss_info ) ) {
			return false;
		}

		TR << hss << " ";
		hss.calc_geometry( ss_info );
		core::Real const dist_val = compute_dist( hss );
		core::Real const angle_val = compute_angle( hss.hs_angle() );

		TR << "hsheet_dist=" << dist_val << ", hs_angle=" << angle_val
			<< ", hs_dist1=" << hss.hs1_dist() << ", hs_dist2=" << hss.hs2_dist() << std::endl;

		if ( dist_val < filter_min_dist_ || dist_val > filter_max_dist_ ) {
			filter = false;
		}

		if ( angle_val < filter_min_angle_ || angle_val > filter_max_angle_ ) {
			filter = false;
		}

		if ( hss.hs1_dist() < filter_min_dist_ || hss.hs1_dist() > filter_max_dist_ ) {
			filter = false;
		}

		if ( hss.hs2_dist() < filter_min_dist_ || hss.hs2_dist() > filter_max_dist_ ) {
			filter = false;
		}
	}

	if ( filter ) {
		TR << " Filter success ! " << std::endl;
	} else {
		TR << " Filter failed ! " << std::endl;
	}

	return filter;

} // apply_filter

// @brief returns computed value
HSSTripletFilter::Real
HSSTripletFilter::compute( Pose const & pose ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::SS_Info2_COP;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::HSSTripletOP;
	using protocols::fldsgn::topology::HSSTriplets;

	// get secondary structure
	std::string const secstruct = get_secstruct( pose );
	runtime_assert( secstruct.length() == pose.size() );

	// get HSSTriplets to filter
	HSSTriplets const hss3s = get_hss3s( pose );
	if ( output_id_ > hss3s.size() ) {
		TR.Error << "[ERROR] The value of output_id is more than the number of input hsstriplets " << std::endl;
		runtime_assert( false );
	}

	// set SS_Info
	SS_Info2_COP const ss_info( SS_Info2_OP( new SS_Info2( pose, secstruct ) ) );

	// check conformation hsstriplets
	Size current_id( 0 );
	for ( HSSTriplets::const_iterator it=hss3s.begin(); it!=hss3s.end(); ++it ) {
		++current_id;
		if ( current_id != output_id_ ) continue;
		HSSTriplet hss( **it );

		if ( !check_elements( hss, *ss_info ) ) {
			utility_exit();
		}

		hss.calc_geometry( ss_info );

		if ( output_type_ == "dist" ) {
			return hss.hsheet_dist();
		} else if ( output_type_ == "angle" ) {
			return compute_angle( hss.hs_angle() );
		} else {
			std::stringstream msg;
			msg << "HSSTripletFilter::compute(): Invalid output_type " << output_type_ << std::endl;
			utility_exit_with_message( msg.str() );
		}
	}

	return 0.0;
}

/// @brief parse xml
void
HSSTripletFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::jd2::parser::BluePrint;

	// set filtered helix_pairings
	String const hss3s = tag->getOption<String>( "hsstriplets", "" );
	if ( hss3s != ""  ) {
		add_hsstriplets( hss3s );
	} else { }

	// secondary strucuture info
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
		if ( !hss3set_ ) {
			add_hsstriplets( blue.hss_triplets() );
		} else {
			TR.Info << "hsstriplets are already spectified to HSSTripletFilter, so the definitions in the blueprint file will be ignored." << std::endl;
		}
	}

	filter_min_dist_  = tag->getOption<Real>( "min_dist", filter_min_dist_ );
	filter_max_dist_  = tag->getOption<Real>( "max_dist", filter_max_dist_ );
	filter_min_angle_ = tag->getOption<Real>( "min_angle", filter_min_angle_ );
	filter_max_angle_ = tag->getOption<Real>( "max_angle", filter_max_angle_ );

	ignore_helix_direction_ = tag->getOption<bool>( "ignore_helix_direction", ignore_helix_direction_ );
	use_dssp_ = tag->getOption<bool>( "use_dssp", use_dssp_ );
	output_id_ = tag->getOption<Size>( "output_id", output_id_ );
	output_type_ = tag->getOption<String>( "output_type", output_type_ );

	if ( output_type_ != "dist" && output_type_ != "angle" ) {

		TR << "Invalid type of output_type, choose either dist or angle. " << std::endl;

	} else {

		if ( hss3set_ ) {
			if ( output_id_ > hss3set_->hss_triplets().size() ) {
				TR << "[ERROR] The value of output_id is more than the number of input hsstriplets " << std::endl;
				runtime_assert( false );
			}

			TR << "HSSTriplet " << *hss3set_->hss_triplets()[output_id_] << ", "
				<< output_type_ << " is used for output value. " << std::endl;
		} else {
			TR << "HSSTriplet filter will determing triplets and check output_id at run time" << std::endl;
		}
	}
}

/// @brief    checks secondary structure elements in the triplet, returns false if invalid
/// @details  If the pose doesn't contain the helix or strands, returns false
///           If the length of the helix is < 5, returns false
///           If the length of either strand is < 2, returns false
bool
HSSTripletFilter::check_elements(
	topology::HSSTriplet const & hss,
	topology::SS_Info2 const & ss_info ) const
{
	topology::Helices const & helices = ss_info.helices();
	topology::Strands const & strands = ss_info.strands();

	if ( !( helices.size() >= hss.helix() ) || helices[ hss.helix() ]->length() < 5 ) {
		TR << "Helix " << hss.helix() << " dones not exist, or is too short. " << std::endl;
		return false;
	}
	if ( !( strands.size() >= hss.strand1() ) || strands[ hss.strand1() ]->length() < 2 ) {
		TR << "Strand1 " << hss.strand1() << " dones not exist, or is too short. " << std::endl;
		return false;
	}
	if ( !( strands.size() >= hss.strand2() ) || strands[ hss.strand2() ]->length() < 2 ) {
		TR << "Strand2 " << hss.strand2() << " dones not exist, or is too short. " << std::endl;
		return false;
	}
	return true;
}

/// @brief given an HSS triplet, compute the distance from helix to sheet
HSSTripletFilter::Real
HSSTripletFilter::compute_dist( HSSTriplet const & hss ) const
{
	//Real const hsheet_dist = hss.hsheet_dist();
	Real const hsheet_dist = ( std::abs( hss.hs1_dist() ) + std::abs( hss.hs2_dist() ) ) / 2.0;
	// THe HSS triplet code can return a negative distance if the helix
	// is not oriented the expected way.  This is a patch for that, where if
	// we don't care about the direction of the helix, the distance is always positive.
	if ( hsheet_dist < 0.0 ) return -hsheet_dist;
	return hsheet_dist;
}

/// @brief given an hs-angle, return a valid angle accounting for ignore_helix_direction_
///        if ignore_helix_direction_ is true, this basically makes angle periodic from
///        -90 to 90
HSSTripletFilter::Real
HSSTripletFilter::compute_angle( core::Real const angle ) const
{
	// ignore_helix_direction essentially makes the sheet-helix angle periodic from -90 to 90
	Real angle_val = angle;
	if ( ignore_helix_direction_ ) {
		while ( angle_val < -90.0 ) {
			angle_val += 180;
		}
		while ( angle_val > 90.0 ) {
			angle_val -= 180;
		}
	}
	return angle_val;
}

/// @brief    computes and returns secondary structure string to use for this filter.
/// @details  if secstruct_ is non-empty, returns that
///           if use_dssp_ is true, use DSSP to compute secstruct
///           otherwise, use secstruct stored in the pose
HSSTripletFilter::String
HSSTripletFilter::get_secstruct( Pose const & pose ) const
{
	if ( !secstruct_.empty() ) return secstruct_;
	if ( !use_dssp_ ) return pose.secstruct();
	core::scoring::dssp::Dssp dssp( pose );
	return dssp.get_dssp_secstruct();
}

/// @brief    returns HSSTriplets object to use for this filter
/// @details  if hss triplets are given prior to apply time, returns those
///           otherwise, look for HSS info in the pose's StructureData and return that
HSSTripletFilter::HSSTriplets
HSSTripletFilter::get_hss3s( Pose const & pose ) const
{
	if ( hss3set_ ) return hss3set_->hss_triplets();

	// try to get triplets from the pose
	using protocols::denovo_design::components::SegmentPairing;
	using protocols::denovo_design::components::StructureData;
	using protocols::denovo_design::components::StructureDataFactory;

	if ( StructureDataFactory::get_instance()->has_cached_data( pose ) ) {
		StructureData const & sd = StructureDataFactory::get_instance()->get_from_const_pose( pose );
		std::string const hss3_str = SegmentPairing::get_hss_triplets( sd );
		TR << "Determining hss triplet string from StructureData: " << hss3_str << std::endl;
		HSSTripletSet const hss3( hss3_str );

		return hss3.hss_triplets();
	}

	std::stringstream msg;
	msg << "HSSTripletFilter::get_hss3s(): No user-specified HSS triplets were found, "
		<< "and no StructureData was found in the pose cache. Therefore, the desired HSS triplets coult not be determined."
		<< std::endl;
	utility_exit_with_message( msg.str() );
	return HSSTriplets();
}

protocols::filters::FilterOP
HSSTripletFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HSSTripletFilter ); }

std::string
HSSTripletFilterCreator::keyname() const { return "HSSTriplet"; }

} // filters
} // fldsgn
} // protocols
