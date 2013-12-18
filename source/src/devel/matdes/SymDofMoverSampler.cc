// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @brief
/// @author Jacob Bale ( balej@uw.edu )

// Unit headers
#include <devel/matdes/SymDofMoverSampler.hh>

// Project headers
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Basic headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <string>

static basic::Tracer TR("devel.matdes.SymDofMoverSampler");

namespace devel {
namespace matdes {

SymDofMoverSampler * SymDofMoverSampler::instance_( 0 );

SymDofMoverSampler::SymDofMoverSampler():
	sym_dof_names_(),
	angles_(),
	radial_disps_(),
	angles_range_min_(),
	angles_range_max_(),
	angle_steps_(),
	radial_disps_range_min_(),
	radial_disps_range_max_(),
	radial_disp_steps_(),
	current_angles_(),
	current_radial_disps_()
	{ }


#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex SymDofMoverSampler::singleton_mutex_;

std::mutex & SymDofMoverSampler::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
SymDofMoverSampler & SymDofMoverSampler::get_instance()
{
	boost::function< SymDofMoverSampler * () > creator = boost::bind( &SymDofMoverSampler::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return *instance_;
}

SymDofMoverSampler *
SymDofMoverSampler::create_singleton_instance()
{
	return new SymDofMoverSampler;
}

void
SymDofMoverSampler::set_sym_dof_names( utility::vector1<std::string> sym_dof_names) {
	sym_dof_names_ = sym_dof_names;
}

void
SymDofMoverSampler::set_angles( utility::vector1<Real> angles) {
	angles_ = angles;
}

void
SymDofMoverSampler::set_radial_disps( utility::vector1<Real> radial_disps) {
	radial_disps_ = radial_disps;
}

void
SymDofMoverSampler::set_angle_ranges(utility::vector1<Real> angles_range_min, utility::vector1<Real> angles_range_max, utility::vector1<Real> angle_steps) {
	if ( angles_range_min.size() != angles_range_max.size() || angles_range_min.size() != angle_steps.size() ) {
		utility_exit_with_message("SymDofMoverSampler angle lists are not all the same length");
	}
	for (Size i = 1; i <= angles_range_min.size(); i++) {
		if( angle_steps[i] > angles_range_max[i] - angles_range_min[i]) {
			utility_exit_with_message("angle_step has to be less than or equal to (angle_max - angle_min)");
		}
		angles_range_min_.push_back(angles_range_min[i]);
		current_angles_.push_back(angles_range_min[i]);
		angles_range_max_.push_back(angles_range_max[i]);
		angle_steps_.push_back(angle_steps[i]);
		TR << "angle range set to [" << angles_range_min[i] << ", " << angles_range_max[i] << "] with a step of " << angle_steps[i] << std::endl;
	}
}

void
SymDofMoverSampler::set_radial_disp_ranges(utility::vector1<Real> radial_disps_range_min, utility::vector1<Real> radial_disps_range_max, utility::vector1<Real> radial_disp_steps ) {
	if ( radial_disps_range_min.size() != radial_disps_range_max.size() || radial_disps_range_min.size() != radial_disp_steps.size() ) {
		utility_exit_with_message("SymDofMoverSampler radial disps lists are not all the same length");
	}
	for (Size i = 1; i <= radial_disps_range_min.size(); i++) {
		if( radial_disp_steps[i] > radial_disps_range_max[i] - radial_disps_range_min[i]) {
			utility_exit_with_message("radial_disp_step has to be less than or equal to (radial_disp_max - radial_disp_min)");
		}
		radial_disps_range_min_.push_back(radial_disps_range_min[i]);
		current_radial_disps_.push_back(radial_disps_range_min[i]);
		radial_disps_range_max_.push_back(radial_disps_range_max[i]);
		radial_disp_steps_.push_back(radial_disp_steps[i]);
		TR << "disp range set to [" << radial_disps_range_min[i] << ", " << radial_disps_range_max[i] << "] with a step of " << radial_disp_steps[i] << std::endl;
	}
}

// loops between radial_disps_range_min_ and radial_disps_range_max_ and
// between angles_range_min and anglse_range_max for each of the current_radial_disps_
void
SymDofMoverSampler::step() {
	TR << "Performing step:" << std::endl;
	utility::vector1<Real> new_angles;
	utility::vector1<Real> new_radial_disps;
	std::string current_values = "\t" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb:";
	for (Size i = 1; i <= radial_disp_steps_.size(); i++) {
		current_values.append(" [radial_disp " + ObjexxFCL::string_of(i) + " = " + ObjexxFCL::string_of(current_radial_disps_[i]) + "] [angle " + ObjexxFCL::string_of(i) + " = " + ObjexxFCL::string_of(current_angles_[i]) + "]");
		new_angles.push_back(current_angles_[i] + angle_steps_[i]);
		new_radial_disps.push_back(current_radial_disps_[i] + radial_disp_steps_[i]);
	}

	TR << current_values << std::endl;

	if(new_angles[1] <= angles_range_max_[1] && !(angles_range_min_[1] == angles_range_max_[1]) ) {
		current_angles_[1] = new_angles[1];
  }	else {
		current_angles_[1] = angles_range_min_[1];
		if( new_radial_disps[1] <= radial_disps_range_max_[1] && !(radial_disps_range_min_[1] == radial_disps_range_max_[1]) ) {
			current_radial_disps_[1] = new_radial_disps[1];
		} else {
			current_radial_disps_[1] = radial_disps_range_min_[1];
			if( angle_steps_.size() == 2) {
				if( new_angles[2] <= angles_range_max_[2] && !(angles_range_min_[2] == angles_range_max_[2]) ) {
					current_angles_[2] = new_angles[2];
				} else {
					current_angles_[2] = angles_range_min_[2];
					if( new_radial_disps[2] <= radial_disps_range_max_[2] && !(radial_disps_range_min_[2] == radial_disps_range_max_[2]) ) {
						current_radial_disps_[2] = new_radial_disps[2];
					} else {
						current_radial_disps_[2] = radial_disps_range_min_[2];
					}
				}
			}
		}
	}

	std::string new_values = "\tNext:";
	for (Size i = 1; i <= radial_disp_steps_.size(); i++) {
		new_values.append(" [radial_disp " + ObjexxFCL::string_of(i) + " = " + ObjexxFCL::string_of(current_radial_disps_[i]) + "] [angle " + ObjexxFCL::string_of(i) + " = " + ObjexxFCL::string_of(current_angles_[i]) + "]");
	}
	TR << new_values << std::endl;
}


} // matdes
} // devel

