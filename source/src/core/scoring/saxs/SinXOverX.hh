// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/scoring/saxs/SAXSEnergyCreatorFA.hh
/// @brief  Declaration for the class that connects SAXSEnergyCreator with the ScoringManager
/// @author Dominik Gront dgront@chem.uw.edu.pl

#ifndef INCLUDED_core_scoring_saxs_SinXOverX_hh
#define INCLUDED_core_scoring_saxs_SinXOverX_hh

#include <core/types.hh>
#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace core {
namespace scoring {
namespace saxs {

class SinXOverX {
public:

	static SinXOverX* get_instance();

	core::Real evaluate(core::Real x) const {
		core::Size tmp_i = ((Size) (x * 100) + 1);
		return sin_x_over_x_[tmp_i];
	}

#ifdef MULTI_THREADED
#ifdef CXX11
public:

	/// @brief This public method is meant to be used only by the
	/// utility::thread::safely_create_singleton function and not meant
	/// for any other purpose.  Do not use.
	static std::mutex & singleton_mutex();

private:
	static std::mutex singleton_mutex_;
#endif
#endif

private:

	SinXOverX();
	void fill_sin_x_over_x_table();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static SinXOverX * create_singleton_instance();

private:
	static SinXOverX* instance_;
	static utility::vector1<Real> sin_x_over_x_;

};

}
}
}

#endif
