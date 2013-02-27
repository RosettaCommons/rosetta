// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/SSShapeComplementarityFilter.cc
/// @brief Tom's Denovo design protocol
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/filters/SSShapeComplementarityFilter.hh>
#include <protocols/denovo_design/filters/SSShapeComplementarityFilterCreator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>
#include <protocols/toolbox/SelectResiduesByLayer.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

//C++ Headers


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif


static basic::Tracer TR("protocols.denovo_design.SSShapeComplementarityFilter");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {

std::string
SSShapeComplementarityFilterCreator::keyname() const
{
	return SSShapeComplementarityFilterCreator::filter_name();
}

protocols::filters::FilterOP
SSShapeComplementarityFilterCreator::create_filter() const {
	return new SSShapeComplementarityFilter();
}

std::string
SSShapeComplementarityFilterCreator::filter_name()
{
	return "SSShapeComplementarity";
}

///  ---------------------------------------------------------------------------------
///  SSShapeComplementarityFilter main code:
///  ---------------------------------------------------------------------------------
SSShapeComplementarityFilter::SSShapeComplementarityFilter() :
	Filter( "SSShapeComplementarityFilter" ),
	sc_( new simple_filters::ShapeComplementarityFilter() ),
	blueprint_( NULL )
{
}

SSShapeComplementarityFilter::SSShapeComplementarityFilter( SSShapeComplementarityFilter const & rval ) :
	Filter( rval ),
	sc_( new simple_filters::ShapeComplementarityFilter( *(rval.sc_) ) ),
	blueprint_( rval.blueprint_ )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
SSShapeComplementarityFilter::~SSShapeComplementarityFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
SSShapeComplementarityFilter::clone() const
{
	return new SSShapeComplementarityFilter(*this);
}

protocols::filters::FilterOP
SSShapeComplementarityFilter::fresh_instance() const
{
	return new SSShapeComplementarityFilter();
}

void
SSShapeComplementarityFilter::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	std::string const bp_filename( tag->getOption< std::string >( "blueprint", "" ) );
	if ( bp_filename != "" ) {
		blueprint_ = new jd2::parser::BluePrint( bp_filename );
	}
	if ( ! blueprint_ ) {
		utility_exit_with_message( "BluePrint must be specified to SSShapeComplementarityFilter" );
	}
	sc_->parse_my_tag( tag, data, filters, movers, pose );
}

std::string
SSShapeComplementarityFilter::get_name() const
{
	return "SSShapeComplementarity";
}

void
SSShapeComplementarityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	sc_->report( out, pose );
}

core::Real
SSShapeComplementarityFilter::report_sm( core::pose::Pose const & pose ) const
{
	fldsgn::topology::SS_Info2 ss_info( blueprint_->secstruct() );

	// we will average out the shape complementarity from HSS triplets and Helix-Helix pairings
	fldsgn::topology::HSSTriplets hss_triplets( fldsgn::topology::HSSTripletSet( blueprint_->hss_triplets() ).hss_triplets() );
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=hss_triplets.size(); ++i ) {
		setup_sc_hss( pose, ss_info, hss_triplets[i] );
		sum += sc_->report_sm( pose );
		TR << "SUM=" << sum << std::endl;
	}

	fldsgn::topology::HelixPairings helix_pairs( fldsgn::topology::HelixPairingSet( blueprint_->helix_pairings() ).helix_pairings() );
	for ( core::Size i=1; i<=helix_pairs.size(); ++i ) {
		setup_sc_hh( pose, ss_info, helix_pairs[i] );
		sum += sc_->report_sm( pose );
		TR << "SUM=" << sum << std::endl;
	}
	TR << "Returning " << sum /( hss_triplets.size() + helix_pairs.size() ) << std::endl;
	return sum /( hss_triplets.size() + helix_pairs.size() );
}

/// @brief Does the SSShapeComplementarity Filtering
bool
SSShapeComplementarityFilter::apply( core::pose::Pose const & pose ) const
{
	return report_sm( pose );
}

/// @brief sets up the underlying shapecomplementarity filter to work based on secondary structure elements
void
SSShapeComplementarityFilter::setup_sc_hss( core::pose::Pose const & pose,
																						fldsgn::topology::SS_Info2 const & ss_info,
																						fldsgn::topology::HSSTripletCOP hss_triplet ) const
{
	runtime_assert( blueprint_ );
	runtime_assert( sc_ );
	utility::vector1< core::Size > set1, set2;
	// set 1 is the helix
	fldsgn::topology::HelixCOP const helix( ss_info.helix( hss_triplet->helix() ) );
	for ( core::Size i=helix->begin(); i<=helix->end(); ++i ) {
		set1.push_back( i );
	}
	sc_->residues1( set1 );
	fldsgn::topology::StrandCOP const strand( ss_info.strand( hss_triplet->strand1() ) );
	for ( core::Size i=strand->begin(); i<=strand->end(); ++i ) {
		set2.push_back( i );
	}
	fldsgn::topology::StrandCOP const strand2( ss_info.strand( hss_triplet->strand2() ) );
	for ( core::Size i=strand2->begin(); i<=strand2->end(); ++i ) {
		set2.push_back( i );
	}
	sc_->residues2( set2 );
	TR.Debug << "Set residues1 to ";
	for ( core::Size i=1; i<=set1.size(); ++i )
		TR.Debug << set1[i] << " ";
	TR.Debug << std::endl << "Set residues2 to ";
	for ( core::Size i=1; i<=set2.size(); ++i )
		TR.Debug << set2[i] << " ";
	TR.Debug << std::endl;
}

/// @brief sets up the underlying shapecomplementarity filter to work based on secondary structure elements
void
SSShapeComplementarityFilter::setup_sc_hh( core::pose::Pose const & pose,
																					 fldsgn::topology::SS_Info2 const & ss_info,
																					 fldsgn::topology::HelixPairingCOP helix_pair ) const
{
	runtime_assert( blueprint_ );
	runtime_assert( sc_ );
	utility::vector1< core::Size > set1, set2;
	// set 1 is helix1, set2 is helix2
	fldsgn::topology::HelixCOP const helix( ss_info.helix( helix_pair->h1() ) );
	for ( core::Size i=helix->begin(); i<=helix->end(); ++i ) {
		set1.push_back( i );
	}
	sc_->residues1( set1 );

	fldsgn::topology::HelixCOP const helix2( ss_info.helix( helix_pair->h2() ) );
	for ( core::Size i=helix2->begin(); i<=helix2->end(); ++i ) {
		set2.push_back( i );
	}
	sc_->residues2( set2 );
	TR.Debug << "Set residues1 to ";
	for ( core::Size i=1; i<=set1.size(); ++i )
		TR.Debug << set1[i] << " ";
	TR.Debug << std::endl << "Set residues2 to ";
	for ( core::Size i=1; i<=set2.size(); ++i )
		TR.Debug << set2[i] << " ";
	TR.Debug << std::endl;
}

} // namespace filters
} // namespace denovo_design
} // namespace protocols


