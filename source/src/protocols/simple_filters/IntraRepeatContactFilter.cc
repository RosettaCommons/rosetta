// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/IntraRepeatContactFilter
/// @brief filter structures by IntraRepeatContacts
/// @detailed
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/IntraRepeatContactFilter.hh>
#include <protocols/simple_filters/IntraRepeatContactFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//// C++ headers
static basic::Tracer tr("protocols.filters.IntraRepeatContactFilter");

namespace protocols {
namespace simple_filters {

// @brief default constructor
IntraRepeatContactFilter::IntraRepeatContactFilter():
    Filter( "IntraRepeatContactsPerResidue"),
    filtered_value_( -1.0 ),
    numbRepeats_(4),
    sequenceSep_(6),
    distThresh_(10)
{}

// @brief copy constructor
IntraRepeatContactFilter::IntraRepeatContactFilter( IntraRepeatContactFilter const & rval ):
    Super( rval ),
    filtered_value_( rval.filtered_value_ ),
    numbRepeats_(rval.numbRepeats_),
    sequenceSep_(rval.sequenceSep_),
    distThresh_(rval.distThresh_)
{}

// @brief set filtered value
void IntraRepeatContactFilter::filtered_value( Real const & value )
{
    filtered_value_ = value;
}

/// @brief
IntraRepeatContactFilter::Real
IntraRepeatContactFilter::report_sm( Pose const & pose ) const
{
    return  compute( pose );
}

/// @brief
void
IntraRepeatContactFilter::report( std::ostream & out, Pose const & pose ) const
{
    out << "IntraRepeatContacts: " <<  compute( pose ) << std::endl;
}

/// @brief
IntraRepeatContactFilter::Real
IntraRepeatContactFilter::compute( Pose const & pose ) const
{
    using numeric::xyzVector;
    Size repeatLength = pose.total_residue()/numbRepeats_;
    Size resStart = 1;
    Size resEnd = repeatLength;
    Size contactCt = 0;
    for(Size ii=resStart; ii<=resEnd; ++ii){
        for(Size jj=ii+sequenceSep_; jj<=resEnd; ++jj){
            xyzVector<double> a = pose.xyz(core::id::NamedAtomID("CA", ii));
            xyzVector<double> b = pose.xyz(core::id::NamedAtomID("CA", jj));
            if(a.distance(b)<distThresh_){
                contactCt+=1;
            }
        }
    }
  return((Real(contactCt))/Real(repeatLength));
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool IntraRepeatContactFilter::apply( Pose const & pose ) const
{
    Real value = compute( pose );
    if( value > filtered_value_ ){
        tr << "Successfully filtered: " << value << std::endl;
        return true;
    }else{
        tr << "Filter failed current/threshold=" << value << "/" << filtered_value_ << std::endl;
        return false;
    }
} // apply_filter

/// @brief parse xml
void
IntraRepeatContactFilter::parse_my_tag(
    TagCOP const tag,
    basic::datacache::DataMap &,
    filters::Filters_map const &,
    Movers_map const &,
    Pose const & )
{
    // set threshold
    filtered_value_ = tag->getOption<Real>( "threshold", 0.0 );
    numbRepeats_ = tag->getOption<Size>("numb_repeats",4);
    sequenceSep_ = tag->getOption<Size>("sequenceSeperation",6);
    distThresh_ = tag->getOption<Size>("distanceThreshold",10.0);
    tr << "Structures which have IntraRepeatContacts less than " << filtered_value_ << " will be filtered." << std::endl;
}

filters::FilterOP
IntraRepeatContactFilterCreator::create_filter() const { return protocols::filters::FilterOP(new IntraRepeatContactFilter);}

std::string
IntraRepeatContactFilterCreator::keyname() const { return "IntraRepeatContactsPerResidue"; }
} // filters
} // protocols
