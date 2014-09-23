// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/HelixPairingFilter.cc
/// @brief filter structures by sheet topology
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HelixPairingFilter.hh>
#include <protocols/fldsgn/filters/HelixPairingFilterCreator.hh>

// Package Headers
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
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

#ifdef WIN32
	#include <protocols/fldsgn/topology/HSSTriplet.hh>
#endif


//// C++ headers
static thread_local basic::Tracer TR( "protocols.fldsgn.filters.HelixPairingFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
HelixPairingFilter::HelixPairingFilter():
	Filter( "HelixPairing" ),
	secstruct_( "" ),
	dist_cutoff_( 15.0 ),
	bend_angle_( 20.0 ),
	cross_angle_( 45.0 ),
	align_angle_( 25.0 )
{}

// @Brief constructor with arguments
HelixPairingFilter::HelixPairingFilter( HelixPairings const & hpairs ):
	Filter( "HelixPairing" ),
	secstruct_( "" ),
	dist_cutoff_( 15.0 ),
	bend_angle_( 20.0 ),
	cross_angle_( 45.0 ),
	align_angle_( 25.0 ),
	hpairset_( HelixPairingSetOP( new HelixPairingSet( hpairs ) ) )
{}

// @brief constructor with arguments
HelixPairingFilter::HelixPairingFilter( String const & hpairs ):
	Filter( "HelixPairing" ),
	secstruct_( "" ),
	dist_cutoff_( 15.0 ),
	bend_angle_( 20.0 ),
	cross_angle_( 45.0 ),
	align_angle_( 25.0 )
{
	hpairset_ = HelixPairingSetOP( new HelixPairingSet( hpairs ) );
}

// @brief copy constructor
HelixPairingFilter::HelixPairingFilter( HelixPairingFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	secstruct_( rval.secstruct_ ),
	dist_cutoff_( rval.dist_cutoff_ ),
	bend_angle_( rval.bend_angle_ ),
	cross_angle_( rval.cross_angle_ ),
	align_angle_( rval.align_angle_ ),
	hpairset_( rval.hpairset_ )
{}

// @brief set filtered sheet_topology by HelixPairings
void HelixPairingFilter::helix_pairings( HelixPairings const & hpairs )
{
	hpairset_ = HelixPairingSetOP( new HelixPairingSet( hpairs ) );
}

// @brief set filtered sheet_topology by SrandPairingSetOP
void HelixPairingFilter::helix_pairings( String const & hpairs )
{
	hpairset_ = HelixPairingSetOP( new HelixPairingSet( hpairs ) );
}

/// @brief
void HelixPairingFilter::secstruct( String const & ss )
{
	secstruct_ = ss;
}

/// @brief
HelixPairingFilter::Real
HelixPairingFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}

/// @brief
bool
HelixPairingFilter::apply( Pose const & pose ) const
{
	using core::Vector;
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;


	HelixPairingSet hpairset( *hpairset_ );
	HelixPairings helix_pairings = hpairset.helix_pairings();

	//
	runtime_assert( ! helix_pairings.empty() );

	// set secondary structure
	if( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	runtime_assert( secstruct_.length() == pose.total_residue() );

	// set SS_Info
	SS_Info2_OP  ss_info( new SS_Info2( pose, secstruct_ ) );
	Helices const & helices( ss_info->helices() );

	// calc geometry of helixpairing
	hpairset.calc_geometry( ss_info );

	// check conformation of helical pairings
	bool filter( true );
	for( HelixPairings::const_iterator it=helix_pairings.begin(), ite=helix_pairings.end(); it != ite; ++it ) {
		HelixPairing const & hpair( **it );

		if( helices.size() < hpair.h1() ) {
			TR << "Helix " << hpair.h1() << " does not exist ! " << std::endl;
			return false;
		}

		if( helices.size() < hpair.h2() ) {
			TR << "Helix " << hpair.h2() << " does not exist ! " << std::endl;
			return false;
		}

		// bend check
		if( bend_angle_ >= 0.0 ) {
			if ( helices[ hpair.h1() ]->bend() > bend_angle_ ) {
				TR << "Helix " << hpair.h1() << "is bend, angle=" << helices[ hpair.h1() ]->bend() << std::endl;
				filter=false;
			}
			if ( helices[ hpair.h2() ]->bend() > bend_angle_ ) {
				TR << "Helix bend " << hpair.h2() << " " << helices[ hpair.h2() ]->bend() << std::endl;
				filter=false;
			}
		}

		//
		if( hpair.dist() > dist_cutoff_ && dist_cutoff_ >= 0 ) filter = false;
		if( hpair.cross_angle() > cross_angle_ && cross_angle_ >= 0 ) filter = false;

		if( hpair.align_angle() < 0 ) continue;
		if( hpair.align_angle() > align_angle_ && align_angle_ >= 0 ) filter = false;

		if( filter == false ) break;
	}

	TR << " Filter condition: " << std::endl;
	TR << " bend ( intra helix ) <= " << bend_angle_ << std::endl;
	TR << " dist <= " << dist_cutoff_  << std::endl;
	TR << " cross <= " << cross_angle_  << std::endl;
	TR << " align <= " << align_angle_	<< std::endl;

	TR << hpairset;

	if( filter ) {
		TR << " Filter success ! " << std::endl;
	} else {
		TR << " Filter failed ! " << std::endl;
	}

	return filter;

}

// @brief returns true if the given pose passes the filter, false otherwise.
HelixPairingFilter::Real
HelixPairingFilter::compute( Pose const & pose ) const
{
	using core::Vector;
	using core::scoring::dssp::Dssp;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;

	HelixPairingSet hpairset( *hpairset_ );
	HelixPairings helix_pairings = hpairset.helix_pairings();

	//
	runtime_assert( ! helix_pairings.empty() );

	// set secondary structure
	if( secstruct_ == "" ) {
		Dssp dssp( pose );
		secstruct_ = dssp.get_dssp_secstruct();
	}
	runtime_assert( secstruct_.length() == pose.total_residue() );

	// set SS_Info
	SS_Info2_OP  ss_info( new SS_Info2( pose, secstruct_ ) );
	// Helices const & helices( ss_info->helices() ); // Unused variable causes warning.

	// calc geometry of helixpairing
	hpairset.calc_geometry( ss_info );

	Real value( 0.0 );
	if ( output_type_ == "dist" ) {
		value = hpairset.helix_pairing( output_id_ )->dist();
	} else if ( output_type_ == "cross" ) {
		value = hpairset.helix_pairing( output_id_ )->cross_angle();
	} else if ( output_type_ == "align" ) {
		value = hpairset.helix_pairing( output_id_ )->align_angle();
	} else {
		TR << "Invalid type for output_type, choose either dist or cross or align. " << std::endl;
		runtime_assert( false );
	}

	return value;

}

/// @brief parse xml
void
HelixPairingFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::jd2::parser::BluePrint;

	// set filtered helix_pairings
	String const hpairs = tag->getOption<String>( "helix_pairings", "" );
	if( hpairs != ""  ) helix_pairings( hpairs );

	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// HHPAIR line is read for the topology of helix pairings
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if( blueprint != "" ) {
		BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();

		if( ! blue.helix_pairings().empty() ) {
			if( hpairs == "" ) {
				helix_pairings( blue.helix_pairings() );
			} else {
				TR << " HHPAIR in blueprint will be igonared " << std::endl;
			}
		}
	}

	if( hpairset_->helix_pairings().empty() ){
		TR.Error << "Error!, Option of helix_pairings is empty." << std::endl;
		runtime_assert( false );
	}

	dist_cutoff_ = tag->getOption<Real>( "dist",  15 );
	bend_angle_  = tag->getOption<Real>( "bend",  20 );
	cross_angle_ = tag->getOption<Real>( "cross", 45 );
	align_angle_ = tag->getOption<Real>( "align", 25 );

	output_id_ = tag->getOption<Size>( "output_id", 1 );
	output_type_ = tag->getOption<String>( "output_type", "dist" );

	if( output_type_ != "dist" && output_type_ != "cross" && output_type_ != "align" ) {
		TR << "Invalid type for output_type, choose among dist or cross or align. " << std::endl;
	} else {
		TR << "HelixPairing " << *hpairset_->helix_pairing( output_id_ ) << ", "
			 << output_type_ << " is used for output value. " << std::endl;
	}

}

protocols::filters::FilterOP
HelixPairingFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HelixPairingFilter ); }

std::string
HelixPairingFilterCreator::keyname() const { return "HelixPairing"; }

} // filters
} // fldsgn
} // protocols
