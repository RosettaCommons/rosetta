// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/EnergiesByTaskCalculator.cc
/// @brief  EnergiesByTaskCalculator class
/// @author Colin A. Smith

// Unit headers
#include <protocols/mean_field/EnergiesByTaskCalculator.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/scoring/Energies.hh>

// Utility headers
#include <utility>
#include <utility/exit.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <cassert>

#include <utility/vector1.hh>
#include <set>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace utility;
using namespace ObjexxFCL::format;

namespace protocols {
namespace mean_field {


EnergiesByTaskCalculator::EnergiesByTaskCalculator( core::pack::task::PackerTaskCOP task ) :
	task_(std::move( task )),
	total_score_( 0.0 )
{}

EnergiesByTaskCalculator::EnergiesByTaskCalculator(
	EnergiesByTaskCalculator const & calculator
) :
	EnergyDependentCalculator(),
	task_(calculator.task()),
	total_score_(calculator.total())
{}

core::pose::metrics::PoseMetricCalculatorOP
EnergiesByTaskCalculator::clone() const {
	return core::pose::metrics::PoseMetricCalculatorOP( new EnergiesByTaskCalculator(*this) );
}

void EnergiesByTaskCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "total" ) {
		basic::check_cast( valptr, &total_score_, "total expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_score_ );

	} else if ( key == "summary" ) {
		std::ostringstream sstream;
		show(sstream);
		std::string const & summary(sstream.str());
		basic::check_cast( valptr, &summary, "summary expects to return a std::string" );
		(static_cast<basic::MetricValue<std::string> *>(valptr))->set( summary );

	} else {
		basic::Error() << "EnergiesByTaskCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}
}


std::string EnergiesByTaskCalculator::print( std::string const & key ) const {

	if ( key == "total" ) {
		return utility::to_string( total_score_ );
	} else if ( key == "summary" ) {
		std::ostringstream sstream;
		show(sstream);
		return sstream.str();
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}

void
EnergiesByTaskCalculator::recompute(
	Pose const & this_pose
)
{

	runtime_assert(this_pose.energies().energies_updated());

	core::Real total = 0.0;
	core::scoring::EnergyMap weights = this_pose.energies().weights();

	using ScoreTypeVec = utility::vector1<core::scoring::ScoreType>;
	ScoreTypeVec score_types;
	for ( int i = 1; i <= core::scoring::n_score_types; ++i ) {
		auto ii = core::scoring::ScoreType(i);
		if ( weights[ii] != 0 ) score_types.push_back(ii);
	}

	for ( core::Size i = 1, end_i = this_pose.size(); i <= end_i; ++i ) {

		if ( task_->being_packed( i ) || task_->being_designed( i ) ) {
			//core::Real rsd_total = 0.0;
			foreach ( core::scoring::ScoreType score_type, score_types ) {
				core::Real score = (weights[score_type] * this_pose.energies().residue_total_energies(i)[ score_type ]);
				total += score;
			}

			//total += this_pose.energies().residue_total_energy(i);
		}
	}

	total_score_ = total;
}

void
EnergiesByTaskCalculator::show(
	std::ostream & out
) const
{

	out << "Total number of residues being designed or packed " << task_->num_to_be_packed() << std::endl;
	out << "Score " << F(9, 3, total() ) << std::endl;
}

} // mean_field
} // protocols
