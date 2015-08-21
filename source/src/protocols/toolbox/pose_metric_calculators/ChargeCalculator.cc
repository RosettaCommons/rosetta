// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculator/ChargeCalculator.cc
/// @brief  calculator to compute nonlocal/tertiary contacts in a given pose
/// @author Florian Richter

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/ChargeCalculator.hh>

//#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>

// Utility headers
//#include <core/util/Tracer.hh>
#include <utility/exit.hh>
#include <basic/MetricValue.hh>


#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;


namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

ChargeCalculator::ChargeCalculator()
: total_charge_(0.0), total_pos_charges_(0),
	total_neg_charges_(0), SR_total_charge_(0.0),
	SR_total_pos_charges_(0), SR_total_neg_charges_(0)
{}

ChargeCalculator::ChargeCalculator(
	std::set< core::Size > const & special_region
) : total_charge_(0.0), total_pos_charges_(0),
	total_neg_charges_(0), SR_total_charge_(0.0),
	SR_total_pos_charges_(0), SR_total_neg_charges_(0),
	special_region_(special_region)
{}


ChargeCalculator::~ChargeCalculator(){}


void
ChargeCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{

	if ( key == "total_charge" ) {
		basic::check_cast( valptr, &total_charge_, "total_charge expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_charge_ );

	} else if ( key == "total_pos_charges" ) {
		basic::check_cast( valptr, &total_pos_charges_, "total_pos_charges expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( total_pos_charges_ );

	} else if ( key == "total_neg_charges" ) {
		basic::check_cast( valptr, &total_neg_charges_, "total_neg_charges expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( total_neg_charges_ );

	} else if ( key == "SR_total_charge" ) {
		basic::check_cast( valptr, &SR_total_charge_, "SR_total_charge expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( SR_total_charge_ );

	} else if ( key == "SR_total_pos_charges" ) {
		basic::check_cast( valptr, &SR_total_pos_charges_, "SR_total_pos_charges expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( SR_total_pos_charges_ );

	} else if ( key == "SR_total_neg_charges" ) {
		basic::check_cast( valptr, &SR_total_neg_charges_, "SR_total_neg_charges expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( SR_total_neg_charges_ );
	} else {
		basic::Error() << "ChargeCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


std::string
ChargeCalculator::print( std::string const & key ) const
{


	basic::Error() << "ChargeCalculator print function not written yet, developer too lazy, cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

} //print


/// @brief simple: go through sequence and figure out how many charged residue
/// there are. Note: code in here doesn't check for (de)protonated residue types
/// at the moment.
void
ChargeCalculator::recompute( Pose const & this_pose )
{
	using namespace core::chemical;

	total_charge_ = 0.0; total_pos_charges_ = 0; total_neg_charges_ = 0;
	SR_total_charge_ = 0.0; SR_total_pos_charges_ = 0; SR_total_neg_charges_ = 0;


	for ( core::Size i = 1; i <= this_pose.total_residue(); ++i ) {
		if ( !this_pose.residue_type(i).is_protein() ) continue;
		AA i_aa( this_pose.residue_type(i).aa() );

		if ( (i_aa == aa_glu) || (i_aa == aa_asp) ) {
			total_charge_ -= 1.0;
			total_neg_charges_++;
			if ( special_region_.find( i ) != special_region_.end() ) {
				SR_total_charge_ -= 1.0;
				SR_total_neg_charges_++;
			}
		} else if ( (i_aa == aa_arg) || (i_aa == aa_lys) ) {
			total_charge_ += 1.0;
			total_pos_charges_++;
			if ( special_region_.find( i ) != special_region_.end() ) {
				SR_total_charge_ += 1.0;
				SR_total_pos_charges_++;
			}
		}

	} //loop over residues
} //recompute


} //namespace PoseMetricCalculators
} //namespace toolbox
} //namespace protocols
