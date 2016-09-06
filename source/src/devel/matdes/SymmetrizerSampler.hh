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

#ifndef INCLUDED_devel_matdes_SymmetrizerSampler_HH
#define INCLUDED_devel_matdes_SymmetrizerSampler_HH

#include <core/types.hh>

#ifdef MULTI_THREADED
#include <atomic>
#include <mutex>
#endif

namespace devel {
namespace matdes {

class SymmetrizerSampler
{
	typedef core::Real Real;

public:
	static SymmetrizerSampler& get_instance();
	void set_angle_range(Real min_angle, Real max_angle, Real angle_step);
	void set_radial_disp_range(Real min_radial_disp, Real max_radial_disp, Real radial_disp_step );
	Real get_angle() { return current_angle_; }
	Real get_radial_disp() { return current_radial_disp_; }
	void step();

#ifdef MULTI_THREADED
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif

private:
	// Don't implement the methods belowed, this class is a singleton.
	SymmetrizerSampler();
	SymmetrizerSampler(SymmetrizerSampler const&);
	void operator=(SymmetrizerSampler const&);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SymmetrizerSampler * create_singleton_instance();

private:
#if defined MULTI_THREADED
	static std::atomic< SymmetrizerSampler * > instance_;
#else
	static SymmetrizerSampler * instance_;
#endif


	Real angle_min_;
	Real angle_max_;
	Real angle_step_;
	Real radial_disp_min_;
	Real radial_disp_max_;
	Real radial_disp_step_;
	Real current_angle_;
	Real current_radial_disp_;
};

} //matdes
} // devel

#endif
