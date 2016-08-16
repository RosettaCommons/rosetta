// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit headers
#include <devel/matdes/SymmetrizerSampler.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


static THREAD_LOCAL basic::Tracer TR( "devel.matdes.SymmetrizerSampler" );

namespace devel {
namespace matdes {

#if defined MULTI_THREADED && defined CXX11
std::atomic< SymmetrizerSampler * > SymmetrizerSampler::instance_( 0 );
#else
SymmetrizerSampler * SymmetrizerSampler::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex SymmetrizerSampler::singleton_mutex_;

std::mutex & SymmetrizerSampler::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
SymmetrizerSampler & SymmetrizerSampler::get_instance()
{
	boost::function< SymmetrizerSampler * () > creator = boost::bind( &SymmetrizerSampler::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return *instance_;
}

SymmetrizerSampler *
SymmetrizerSampler::create_singleton_instance()
{
	return new SymmetrizerSampler;
}

SymmetrizerSampler::SymmetrizerSampler():
	angle_min_(0),
	angle_max_(0),
	angle_step_(0),
	radial_disp_min_(0),
	radial_disp_max_(0),
	radial_disp_step_(0),
	current_angle_(0),
	current_radial_disp_(0)
{ }

void
SymmetrizerSampler::set_angle_range(Real angle_min, Real angle_max, Real angle_step) {
	runtime_assert_msg(angle_step < angle_max - angle_min, "angle_step has to be smaller than (angle_max - angle_min)");
	angle_min_ = angle_min;
	angle_max_ = angle_max;
	angle_step_ = angle_step;
	current_angle_ = angle_min;
	TR << "angle range set to [" << angle_min << ", " << angle_max << "]" << std::endl;
}

void
SymmetrizerSampler::set_radial_disp_range(Real radial_disp_min, Real radial_disp_max, Real radial_disp_step ) {
	runtime_assert_msg(radial_disp_step < radial_disp_max - radial_disp_min, "radial_disp_step has to be smaller than (radial_disp_max - radial_disp_min");
	radial_disp_min_ = radial_disp_min;
	radial_disp_max_ = radial_disp_max;
	radial_disp_step_ = -1*radial_disp_step;
	current_radial_disp_ = radial_disp_min;
	TR << "radial_disp range set to [" << radial_disp_min << ", " << radial_disp_max << "]" << std::endl;
}

// loops between radial_disp_min_ and radial_disp_max_ and
// between angle_min and angle_max for each current_radial_disp_
void
SymmetrizerSampler::step() {
	TR << "Performing step:" << std::endl;
	TR << "\tbefore: radial_disp = " << current_radial_disp_ << " angle = " << current_angle_ << std::endl;
	Real new_angle = current_angle_ + angle_step_;
	Real new_radial_disp = current_radial_disp_ + radial_disp_step_;

	if ( new_angle < angle_max_ ) {
		current_angle_ = new_angle;
	} else {
		current_angle_ = angle_min_;
		if ( new_radial_disp < radial_disp_max_ ) {
			current_radial_disp_ = new_radial_disp;
		} else {
			current_radial_disp_ = radial_disp_min_;
		}
	}
	TR << "\tafter : radial_disp = " << current_radial_disp_ << " angle = " << current_angle_ << std::endl;
}


} // matdes
} // devel

