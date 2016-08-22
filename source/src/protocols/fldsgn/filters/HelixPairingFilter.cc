// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/HelixPairingFilter.cc
/// @brief filter structures by sheet topology
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HelixPairingFilter.hh>
#include <protocols/fldsgn/filters/HelixPairingFilterCreator.hh>

// Package Headers
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
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

// C++ headers
#include <set>

static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.filters.HelixPairingFilter" );

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
	align_angle_( 25.0 ),
	hpairset_(),
	output_id_( 1 ),
	output_type_( "dist" ),
	use_dssp_( true )
{}

// @Brief constructor with arguments
HelixPairingFilter::HelixPairingFilter( HelixPairings const & hpairs ):
	Filter( "HelixPairing" ),
	secstruct_( "" ),
	dist_cutoff_( 15.0 ),
	bend_angle_( 20.0 ),
	cross_angle_( 45.0 ),
	align_angle_( 25.0 ),
	hpairset_( HelixPairingSetOP( new HelixPairingSet( hpairs ) ) ),
	output_id_( 1 ),
	output_type_( "dist" ),
	use_dssp_( true )
{}

// @brief constructor with arguments
HelixPairingFilter::HelixPairingFilter( String const & hpairs ):
	Filter( "HelixPairing" ),
	secstruct_( "" ),
	dist_cutoff_( 15.0 ),
	bend_angle_( 20.0 ),
	cross_angle_( 45.0 ),
	align_angle_( 25.0 ),
	hpairset_( new HelixPairingSet( hpairs ) ),
	output_id_( 1 ),
	output_type_( "dist" ),
	use_dssp_( true )
{}

// @brief copy constructor
HelixPairingFilter::HelixPairingFilter( HelixPairingFilter const & rval ):
	Super( rval ),
	secstruct_( rval.secstruct_ ),
	dist_cutoff_( rval.dist_cutoff_ ),
	bend_angle_( rval.bend_angle_ ),
	cross_angle_( rval.cross_angle_ ),
	align_angle_( rval.align_angle_ ),
	hpairset_(),
	output_id_( rval.output_id_ ),
	output_type_( rval.output_type_ ),
	use_dssp_( rval.use_dssp_ )
{
	// clone hpairset_ so there are no problems with sharing non-const data
	if ( rval.hpairset_ ) {
		hpairset_ = HelixPairingSetOP( new HelixPairingSet( *rval.hpairset_ ) );
	}
}


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

/// @brief sets distance cutoff
void HelixPairingFilter::dist( Real const dist_val )
{
	dist_cutoff_ = dist_val;
}

/// @brief sets max helix bend
void HelixPairingFilter::bend_angle( Real const bend_angle_val )
{
	bend_angle_ = bend_angle_val;
}

/// @brief sets max cross angle
void HelixPairingFilter::cross_angle( Real const cross_angle_val )
{
	cross_angle_ = cross_angle_val;
}

/// @brief sets max alignment angle
void HelixPairingFilter::align_angle( Real const align_angle_val )
{
	align_angle_ = align_angle_val;
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

	// get secondary structure information
	std::string const pose_secstruct = secstruct( pose );
	runtime_assert( pose_secstruct.length() == pose.total_residue() );
	SS_Info2_OP ss_info( new SS_Info2( pose, pose_secstruct ) );
	Helices const helices( ss_info->helices() );

	// compute desired helix pairings
	HelixPairingSet hpairset = compute_helix_pairing_set( pose );
	HelixPairings const helix_pairings = hpairset.helix_pairings();
	if ( helix_pairings.empty() ) {
		TR << "No pairings found, returning true" << std::endl;
		return true;
	}

	// find missing helices
	std::set< core::Size > const missing_helices = find_missing_helices( *ss_info, helix_pairings );
	if ( !missing_helices.empty() ) {
		TR << "The helix pairing set contains helices " << missing_helices
			<< ", but according to the secondary structure, there are only " << helices.size()
			<< " helices in the pose. Returning false." << std::endl;
		return false;
	}

	// calc geometry of helixpairing
	hpairset.calc_geometry( ss_info );

	// check conformation of helical pairings
	bool filter( true );
	for ( HelixPairings::const_iterator it=helix_pairings.begin(), ite=helix_pairings.end(); it != ite; ++it ) {
		HelixPairing const & hpair( **it );

		// bend check
		if ( bend_angle_ >= 0.0 ) {
			if ( helices[ hpair.h1() ]->bend() > bend_angle_ ) {
				TR << "Helix " << hpair.h1() << "is bent, angle=" << helices[ hpair.h1() ]->bend() << std::endl;
				filter=false;
			}
			if ( helices[ hpair.h2() ]->bend() > bend_angle_ ) {
				TR << "Helix bend " << hpair.h2() << " " << helices[ hpair.h2() ]->bend() << std::endl;
				filter=false;
			}
		}


		if ( hpair.dist() > dist_cutoff_ && dist_cutoff_ >= 0 ) filter = false;
		if ( hpair.cross_angle() > cross_angle_ && cross_angle_ >= 0 ) filter = false;

		if ( hpair.align_angle() < 0 ) continue;
		if ( hpair.align_angle() > align_angle_ && align_angle_ >= 0 ) filter = false;

		if ( filter == false ) break;
	}

	TR << " Filter condition: " << std::endl;
	TR << " bend ( intra helix ) <= " << bend_angle_ << std::endl;
	TR << " dist <= " << dist_cutoff_  << std::endl;
	TR << " cross <= " << cross_angle_  << std::endl;
	TR << " align <= " << align_angle_ << std::endl;

	TR << hpairset;

	if ( filter ) {
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
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;

	// get secondary structure information
	std::string const pose_secstruct = secstruct( pose );
	runtime_assert( pose_secstruct.length() == pose.total_residue() );
	SS_Info2_OP ss_info( new SS_Info2( pose, pose_secstruct ) );

	// compute desired helix pairings
	HelixPairingSet hpairset = compute_helix_pairing_set( pose );
	HelixPairings const helix_pairings = hpairset.helix_pairings();
	if ( helix_pairings.empty() ) {
		TR << "No pairings found, returning true" << std::endl;
		return 0.0;
	}

	// find missing helices
	std::set< core::Size > const missing_helices = find_missing_helices( *ss_info, helix_pairings );
	if ( !missing_helices.empty() ) {
		TR << "The helix pairing set contains helices " << missing_helices
			<< ", but according to the secondary structure, there are only " << ss_info->helices().size()
			<< " helices in the pose. Returning false." << std::endl;
		return 0.0;
	}

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
	if ( hpairs != ""  ) helix_pairings( hpairs );

	// Blueprint is for giving secondary structure information, otherwise dssp will run for ss definition
	// HHPAIR line is read for the topology of helix pairings
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if ( blueprint != "" ) {
		BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();

		if ( ! blue.helix_pairings().empty() ) {
			if ( hpairs == "" ) {
				helix_pairings( blue.helix_pairings() );
			} else {
				TR << " HHPAIR in blueprint will be igonared " << std::endl;
			}
		}
	}

	secstruct_ = tag->getOption< std::string >( "secstruct", secstruct_ );
	use_dssp_ = tag->getOption< bool >( "use_dssp", use_dssp_ );
	/*
	if ( hpairset_->helix_pairings().empty() ) {
	TR.Error << "Error!, Option of helix_pairings is empty." << std::endl;
	runtime_assert( false );
	}
	*/

	dist_cutoff_ = tag->getOption<Real>( "dist",  dist_cutoff_ );
	bend_angle_  = tag->getOption<Real>( "bend",  bend_angle_ );
	cross_angle_ = tag->getOption<Real>( "cross", cross_angle_ );
	align_angle_ = tag->getOption<Real>( "align", align_angle_ );

	output_id_ = tag->getOption<Size>( "output_id", output_id_ );
	output_type_ = tag->getOption<String>( "output_type", output_type_ );

	if ( output_type_ != "dist" && output_type_ != "cross" && output_type_ != "align" ) {
		TR << "Invalid type for output_type, choose among dist or cross or align. " << std::endl;
	} else {
		TR << output_type_ << " is used for output value. " << std::endl;
	}

}

/// @brief returns secondary structure to be used for finding helices
/// @details If secstruct_ is set, returns that.
///          If use_dssp_ is true, returns secstruct from DSSP
///          Otherwise, returns the pose.secstruct()
std::string
HelixPairingFilter::secstruct( core::pose::Pose const & pose ) const
{
	using core::scoring::dssp::Dssp;
	if ( !secstruct_.empty() ) return secstruct_;

	if ( use_dssp_ ) {
		Dssp dssp( pose );
		return dssp.get_dssp_secstruct();
	}

	return pose.secstruct();
}

/// @brief Returns the Helix pairing set to filter based on.
/// @details If hpairset_ is set, returns that.
///          Otherwise, look up the pairing set using the StructureData cached in the pose
topology::HelixPairingSet
HelixPairingFilter::compute_helix_pairing_set( core::pose::Pose const & pose ) const
{
	using namespace protocols::denovo_design::components;

	if ( hpairset_ ) {
		return *hpairset_;
	} else {
		StructureData const sd = StructureDataFactory::get_instance()->create_from_pose( pose );
		std::stringstream helix_str;
		for ( SegmentPairingCOPs::const_iterator p=sd.pairings_begin(); p!=sd.pairings_end(); ++p ) {
			if ( p != sd.pairings_begin() ) helix_str << ';';
			helix_str << (*p)->pairing_string( sd );
		}
		TR << "Found helix pairings: " << helix_str.str() << std::endl;
		return HelixPairingSet( helix_str.str() );
	}
}

std::set< core::Size >
HelixPairingFilter::find_missing_helices(
	topology::SS_Info2 const & ss_info,
	topology::HelixPairings const & hpairs ) const
{
	core::Size const n_helices = ss_info.helices().size();

	std::set< core::Size > missing;
	for ( topology::HelixPairings::const_iterator hp=hpairs.begin(); hp!=hpairs.end(); ++hp ) {
		if ( n_helices < (*hp)->h1() ) missing.insert( (*hp)->h1() );
		if ( n_helices < (*hp)->h2() ) missing.insert( (*hp)->h2() );
	}
	return missing;
}

protocols::filters::FilterOP
HelixPairingFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HelixPairingFilter ); }

std::string
HelixPairingFilterCreator::keyname() const { return "HelixPairing"; }

} // filters
} // fldsgn
} // protocols
