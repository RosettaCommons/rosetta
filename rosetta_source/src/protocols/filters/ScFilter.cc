// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/ScFilter.cc
/// @brief  filter structures by shape complementarity
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

// Unit Headers
#include <protocols/filters/ScFilter.hh>
#include <protocols/filters/ScFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

//Auto Headers


//// C++ headers
static basic::Tracer tr("protocols.filters.ScFilter");

namespace protocols {
namespace filters {

// @brief default constructor
ScFilter::ScFilter():
	Filter( "Sc" ),
	filtered_sc_( 0.50 ),
	jump_id_( 1 ),
	quick_( false ),
	verbose_( false )
{}


// @brief constructor with arguments
ScFilter::ScFilter( Real const & filtered_sc, Size const & jump_id, Size const & quick, Size const & verbose ):
	Filter( "Sc" ),
	filtered_sc_( filtered_sc ),
	jump_id_( jump_id ),
	quick_( quick ),
	verbose_( verbose )
{}

// @brief copy constructor
ScFilter::ScFilter( ScFilter const & rval ):
	Super( rval ),
	filtered_sc_( rval.filtered_sc_ ),
	jump_id_( rval.jump_id_ ),
	quick_( rval.quick_ ),
	verbose_( rval.verbose_ )
{}

void ScFilter::filtered_sc( Real const & filtered_sc ) { filtered_sc_ = filtered_sc; }
void ScFilter::jump_id( Size const & jump_id ) { jump_id_ = jump_id; }
void ScFilter::quick( Size const & quick ) { quick_ = quick; }
void ScFilter::verbose( Size const & verbose ) { verbose_ = verbose; }

/// @brief
ScFilter::Real
ScFilter::compute( Pose const & pose ) const
{
	// We could save the ShapeComplementarityCalculator instance between
	// calculations, but initialization is quick between runs...
	core::scoring::sc::ShapeComplementarityCalculator sc;
	if(quick_)
		sc.settings.density = 5.0;
	if(!sc.Init())
		return -1;
	if(!sc.Calc( pose, jump_id_ ))
		return -1;

	core::scoring::sc::RESULTS const &r = sc.GetResults();
	if(verbose_) {

		// Verbose view
		tr << "==================================================" << std::endl;
		tr << std::endl;
		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << ":" << std::endl;
			else
				tr << "Total/Average for both molecules:" << std::endl;

			tr << "          Total Atoms: " << r.surface[i].nAtoms << std::endl;
			tr << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << std::endl;
			tr << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << std::endl;
			tr << "           Total Dots: " << r.surface[i].nAllDots << std::endl;
			tr << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << std::endl;
			tr << "         Trimmed Area: " << r.surface[i].trimmedArea << " (avg) " << std::endl;
			tr << std::endl;
    }
		tr << std::endl;

		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << std::endl;
			else
				tr << "Average for both molecules:" << std::endl;
			tr << "      Mean Separation: " << r.surface[i].d_mean << std::endl;
			tr << "    Median Separation: " << r.surface[i].d_median << std::endl;
			tr << "    Mean Shape Compl.: " << r.surface[i].s_mean << std::endl;
			tr << "  Median Shape Compl.: " << r.surface[i].s_median << std::endl;
			tr << std::endl;
		}

	}

	tr << "Shape complementarity: " << r.sc << std::endl;

	return r.sc;
}

/// @brief
ScFilter::Real
ScFilter::report_sm( Pose const & pose ) const
{
	return compute( pose );
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose has high enough shape
// complementarity.
bool ScFilter::apply( Pose const & pose ) const
{
	Real sc = compute( pose );
	if( sc >= filtered_sc_ ){
		tr << "Successfully filtered: " << sc << std::endl;
		return true;
	}else{
		tr << "Filter failed current/threshold=" << sc << "/" << filtered_sc_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
ScFilter::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	filtered_sc_ = tag->getOption<Real>( "threshold", 0.50 );
	verbose_ = tag->getOption<Size>( "verbose", false );
	quick_ = tag->getOption<Size>( "quick", false );
	jump_id_ = tag->getOption<Size>( "jumpid", 1 );
	tr << "Structures with shape complementarity <= " << filtered_sc_ << " will be filtered." << std::endl;
	if(quick_)
		tr << "Calculating shape complementarity in quick mode with less accuracy." << std::endl;
	if(jump_id_ != 1)
		tr << "Using Jump ID " << jump_id_ << " to define surfaces." << std::endl;
}

FilterOP
ScFilterCreator::create_filter() const { return new ScFilter; }

std::string
ScFilterCreator::keyname() const { return "Sc"; }


} // filters
} // protocols
