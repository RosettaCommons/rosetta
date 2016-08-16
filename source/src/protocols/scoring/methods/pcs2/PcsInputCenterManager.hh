// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsInputCenterManager.hh
///
/// @brief Singleton that hold everything about the input PCS
/// This avoid multiple reading of the input files.
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsInputCenterManager_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsInputCenterManager_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsInputCenter.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <map>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <atomic>
#include <mutex>
#endif
#endif

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsInputCenterManager {

private:
	PcsInputCenterManager(); //Construct

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PcsInputCenterManager * create_singleton_instance();

private:
#if defined MULTI_THREADED && defined CXX11
	static std::atomic< PcsInputCenterManager * > instance_;
#else
	static PcsInputCenterManager * instance_;
#endif

	std::map<std::string, PcsInputCenter> PcsInputCenter_all_;

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

public:
	static PcsInputCenterManager *
	get_instance();

	/// @ Re init the singleton to default value
	void
	re_init();

	/// @brief Give me the PcsInputCenter given the vector of filename and vector of weight
	PcsInputCenter
	get_PcsInputCenter_for(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & vec_weight);

	/// @brief Output myself on the stream
	friend std::ostream &
	operator<<(std::ostream& out, const PcsInputCenterManager &me);

};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
