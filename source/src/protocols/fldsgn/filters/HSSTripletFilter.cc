// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/HSSTripletFilter.cc
/// @brief filter structures by hsstriplets angle and distance
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HSSTripletFilter.hh>
#include <protocols/fldsgn/filters/HSSTripletFilterCreator.hh>

// Package Headers
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/BB_Pos.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//// C++ headers
static basic::Tracer TR("protocols.fldsgn.filters.HSSTripletFilter");

namespace protocols {
namespace fldsgn {
namespace filters {

// @Brief default constructor
HSSTripletFilter::HSSTripletFilter():
	Filter( "HSSTriplet" ),
	hss3set_( new HSSTripletSet ),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -12.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	output_value_( -990.0 )
{}

// @brief constructor with arguments
HSSTripletFilter::HSSTripletFilter( HSSTriplets const & hss3s ):
	Filter( "HSSTriplet" ),
	hss3set_( new HSSTripletSet( hss3s ) ),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -12.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	output_value_( -990.0 )
{}

// @brief constructor with arguments
HSSTripletFilter::HSSTripletFilter( String const & hss3s ):
	Filter( "HSSTriplet" ),
	hss3set_( new HSSTripletSet( hss3s ) ),
	secstruct_( "" ),
	filter_min_dist_(  7.5 ),
	filter_max_dist_( 13.0 ),
	filter_min_angle_( -2.5 ),
	filter_max_angle_( 90.0 ),
	output_id_( 1 ),
	output_type_( "dist" ),
	output_value_( -990.0 )
{}

// @brief copy constructor
HSSTripletFilter::HSSTripletFilter( HSSTripletFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	hss3set_( rval.hss3set_ ),
	secstruct_( rval.secstruct_ ),
	filter_min_dist_( rval.filter_min_dist_ ),
	filter_max_dist_( rval.filter_max_dist_ ),
	filter_min_angle_( rval.filter_min_angle_ ),
	filter_max_angle_( rval.filter_max_angle_ ),
	output_id_( rval.output_id_ ),
	output_type_( rval.output_type_ ),
	output_value_( rval.output_value_ )
{}

// @brief set filtered HSSTriplets
void HSSTripletFilter::add_hsstriplets( HSSTriplets const & hss3s )
{
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
	filter_min_angle_ = r;
}

/// @brief maximum angle for filtering
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
	using core::Vector;
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::BB_Pos;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::Strand;
	using protocols::fldsgn::topology::Strands;
	using protocols::fldsgn::topology::HSSTripletOP;

	// set secondary structure
	if( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	runtime_assert( secstruct_.length() == pose.total_residue() );

	// set SS_Info
	SS_Info2_OP  ss_info = new SS_Info2( pose, secstruct_ );
	Helices const & helices( ss_info->helices() );
	Strands const & strands( ss_info->strands() );

	// check conformation hsstriplets
	bool filter( true );
	Size current_id( 0 );
	for( HSSTriplets::const_iterator
				 it=hss3set_->hss_triplets().begin(), ite=hss3set_->hss_triplets().end(); it != ite; ++it ) {

		current_id ++;
		HSSTripletOP const & hssop( *it );

		if( ! helices.size() >= hssop->helix() || helices[ hssop->helix() ]->length() < 5 ) {
			TR << "Helix " << hssop->helix() << " dones not exist, or is too short. " << std::endl;
			return false;
		}
		if( ! strands.size() >= hssop->strand1() || strands[ hssop->strand1() ]->length() < 2 ) {
			TR << "Strand1 " << hssop->strand1() << " dones not exist, or is too short. " << std::endl;
			return false;
		}
		if( ! strands.size() >= hssop->strand2() || strands[ hssop->strand2() ]->length() < 2 ) {
			TR << "Strand2 " << hssop->strand2() << " dones not exist, or is too short. " << std::endl;
			return false;
		}
		
		TR << *hssop << " ";
		hssop->calc_geometry( ss_info );

		TR << "hsheet_dist=" << hssop->hsheet_dist() << ", hs_angle=" << hssop->hs_angle()
		   << ", hs_dist1=" << hssop->hs1_dist() << ", hs_dist2=" << hssop->hs2_dist() << std::endl;

		if( hssop->hsheet_dist() < filter_min_dist_ || hssop->hsheet_dist() > filter_max_dist_ ) {
			filter = false;
		}
		if( hssop->hs_angle() < filter_min_angle_ || hssop->hs_angle() > filter_max_angle_ ) {
			filter = false;
		}

		if( output_id_ == current_id ) {
			if ( output_type_ == "dist" ) {
				output_value_ = hssop->hsheet_dist();
			} else if ( output_type_ == "angle" ) {
				output_value_ = hssop->hs_angle();
			}
		}

	}

	if( filter ) {
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
	apply( pose );
	return output_value_;
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
	bool triplets_specified( false );
  String const hss3s = tag->getOption<String>( "hsstriplets", "" );
	if( hss3s != ""  ) {
		HSSTripletSet hss3set( hss3s );
		add_hsstriplets( hss3set.hss_triplets() );
		triplets_specified = true;
	} else {
	}

	// secondary strucuture info
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if( blueprint != "" ) {
		BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
		if ( ! triplets_specified ) {
			HSSTripletSet hss3set( blue.hss_triplets() );
			add_hsstriplets( hss3set.hss_triplets() );
			triplets_specified = true;
		} else {
			TR.Info << "hsstriplets are spectified in the XML, so the definitions in the blueprint file will be ignored." << std::endl;
		}
	}

	// check to make sure we have added hss triplets and throw an error if not
	if ( ! triplets_specified ) {
		TR << "[ERROR] no input of hsstriplets in xml " << std::endl;
		runtime_assert( false );
	}

  filter_min_dist_  = tag->getOption<Real>( "min_dist", 7.5 );
  filter_max_dist_  = tag->getOption<Real>( "max_dist", 13.0);
  filter_min_angle_ = tag->getOption<Real>( "min_angle", -12.5 );
  filter_max_angle_ = tag->getOption<Real>( "max_angle", 90.0 );

	output_id_ = tag->getOption<Size>( "output_id", 1 );
	output_type_ = tag->getOption<String>( "output_type", "dist" );

	if( output_type_ != "dist" && output_type_ != "angle" ) {

		TR << "Invalid type of output_type, choose either dist or angle. " << std::endl;

	} else {

		if( output_id_ > hss3set_->hss_triplets().size() ) {
			TR << "[ERROR] The value of output_id is more than the number of input hsstriplets " << std::endl;
			runtime_assert( false );
		}

		TR << "HSSTriplet " << hss3set_->hss_triplet( output_id_ ) << ", "
			 << output_type_ << " is used for output value. " << std::endl;
	}

}

protocols::filters::FilterOP
HSSTripletFilterCreator::create_filter() const { return new HSSTripletFilter; }

std::string
HSSTripletFilterCreator::keyname() const { return "HSSTriplet"; }

} // filters
} // fldsgn
} // protocols
