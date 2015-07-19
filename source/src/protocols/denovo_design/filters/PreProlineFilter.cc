// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/filters/PreProlineFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/filters/PreProlineFilter.hh>
#include <protocols/denovo_design/filters/PreProlineFilterCreator.hh>

// Project Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/sequence/ABEGOManager.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

//C++ Headers

static thread_local basic::Tracer TR( "protocols.denovo_design.PreProlineFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {

std::string
PreProlineFilterCreator::keyname() const
{
	return PreProlineFilterCreator::filter_name();
}

protocols::filters::FilterOP
PreProlineFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new PreProlineFilter() );
}

std::string
PreProlineFilterCreator::filter_name()
{
	return "PreProline";
}

///  ---------------------------------------------------------------------------------
///  PreProlineFilter main code:
///  ---------------------------------------------------------------------------------
PreProlineFilter::PreProlineFilter() :
	Filter( "PreProlineFilter" ),
	threshold_( 0.0 ),
	selector_()
{}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
PreProlineFilter::~PreProlineFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
PreProlineFilter::clone() const
{
	return protocols::filters::FilterOP( new PreProlineFilter(*this) );
}

protocols::filters::FilterOP
PreProlineFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new PreProlineFilter() );
}

void
PreProlineFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

std::string
PreProlineFilter::get_name() const
{
	return "PreProline";
}

void
PreProlineFilter::set_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector )
{
	selector_ = selector;
}

void
PreProlineFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "PreProlineFilter returning " << compute( pose ) << std::endl;
}

core::Real
PreProlineFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
PreProlineFilter::compute( core::pose::Pose const & pose ) const
{
	// find selected residues
	utility::vector1< bool > selection( pose.total_residue(), true );
	if ( selector_ ) {
		selection = selector_->apply( pose );
	}

	// as a simple first trial, only "B" torsion spaces should be allowed
	// TODO: Use statistical potentials to evaluate specific phi/psi combinations
	core::Size bad_count = core::Size( 0.0 );
	core::Size pro_count = core::Size( 0 );
	std::string const & sequence = pose.sequence();
	std::string const abegos = abego_str( core::sequence::get_abego( pose, 1 ) );
	core::Size resi = 1;
	for ( std::string::const_iterator a = abegos.begin(), s = sequence.begin() + 1;
			( a != abegos.end() ) && ( s != sequence.end() );
			++a, ++s, ++resi ) {
		// ignore if proline is not selected
		if ( ! selection[ resi + 1 ] )
			continue;

		// ignore if this is the end of the chain
		if ( core::pose::is_upper_terminus( pose, resi ) )
			continue;

		if ( *s == 'P' ) {
			++pro_count;
			TR.Debug << "Res " << resi << " " << pose.residue( resi ).name() << " " << pose.residue( resi + 1 ).name() << " " << *a << std::endl;
			if ( *a != 'B' ) {
				++bad_count;
			}
		}
	}
	TR << "Prolines in pose: " << pro_count << " Bad pre-proline torsions: " << bad_count << std::endl;
	return static_cast< core::Real >( bad_count );
}

/// @brief Does the PreProline Filtering
bool
PreProlineFilter::apply( core::pose::Pose const & pose ) const
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


