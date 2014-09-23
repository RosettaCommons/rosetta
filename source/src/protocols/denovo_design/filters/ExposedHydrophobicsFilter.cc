// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/ExposedHydrophobicsFilter.cc
/// @brief Tom's Denovo design protocol
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilter.hh>
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilterCreator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/DsspMover.hh>
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


static thread_local basic::Tracer TR( "protocols.denovo_design.ExposedHydrophobicsFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {

std::string const ExposedHydrophobicsFilter::hydrophobic_residues_ = "MLIVFAWY";

std::string
ExposedHydrophobicsFilterCreator::keyname() const
{
	return ExposedHydrophobicsFilterCreator::filter_name();
}

protocols::filters::FilterOP
ExposedHydrophobicsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter() );
}

std::string
ExposedHydrophobicsFilterCreator::filter_name()
{
	return "ExposedHydrophobics";
}

///  ---------------------------------------------------------------------------------
///  ExposedHydrophobicsFilter main code:
///  ---------------------------------------------------------------------------------
ExposedHydrophobicsFilter::ExposedHydrophobicsFilter() :
	Filter( "ExposedHydrophobicsFilter" ),
	threshold_( -1.0 ),
	sasa_cutoff_( 20.0 )
{}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ExposedHydrophobicsFilter::~ExposedHydrophobicsFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
ExposedHydrophobicsFilter::clone() const
{
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter(*this) );
}

protocols::filters::FilterOP
ExposedHydrophobicsFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter() );
}

void
ExposedHydrophobicsFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	sasa_cutoff_ = tag->getOption< core::Real >( "sasa_cutoff", sasa_cutoff_ );
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );
}

std::string
ExposedHydrophobicsFilter::get_name() const
{
	return "ExposedHydrophobics";
}

void
ExposedHydrophobicsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "ExposedHydrophobicsFilter returning " << compute( pose ) << std::endl;
}

core::Real
ExposedHydrophobicsFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
ExposedHydrophobicsFilter::compute( core::pose::Pose const & pose ) const
{
	toolbox::SelectResiduesByLayer srbl( true, true, true );
	core::Real const surface_sasa( 40 );
	srbl.initialize( sasa_cutoff_, surface_sasa );
	// do dssp on copy of pose to get secondary structure
	core::pose::Pose posecopy( pose );
	moves::DsspMover dssp;
	dssp.apply( posecopy );
	srbl.compute( pose, pose.secstruct() );
	// need to brainstorm the best way of calculating this. There is a cutoff where we don't care. We want to weight bad ones highly. Therefore, maybe a temperature-weighted boltzmann sum is appropriate.
	core::Size residue_count( 0 );
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		//check to see if this residue is a protein residue, and if it's found in the hydrophobic residue list
		if ( pose.residue( i ).is_protein() &&
				 hydrophobic_residues_.find( pose.residue( i ).name1() ) != std::string::npos ) {
			++residue_count;
			// a residue only counts as exposed if it has SASA >= sasa_cutoff_
			if ( srbl.rsd_sasa( i ) >= sasa_cutoff_ ) {
				sum += srbl.rsd_sasa( i ) - sasa_cutoff_;
			}
			TR.Debug << "Resi " << pose.residue( i ).name3() << i << " = " << srbl.rsd_sasa( i ) << ", total sum=" << sum << ", count=" << residue_count << std::endl;
		}
	}
	TR << "ExposedHydrophobics value=" << sum << std::endl;
	return sum;
}

/// @brief Does the ExposedHydrophobics Filtering
bool
ExposedHydrophobicsFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const value( compute( pose ) );
	if ( threshold_ < 0 || value <= threshold_ ) {
		return true;
	} else {
		return false;
	}
}


} // namespace filters
} // namespace denovo_design
} // namespace protocols


