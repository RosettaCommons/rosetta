// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh
 ///
 /// @brief
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified February 2010
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsEnergyParameterManager_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsEnergyParameterManager_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameter.hh>

#include <utility/vector1.hh>

#ifdef MULTI_THREADED
#ifdef CXX11
// C++11 Headers
#include <thread>
#endif
#endif

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsEnergyParameterManager {

public:
	static PcsEnergyParameterManager *
	get_instance();

	friend std::ostream &
	operator<<(std::ostream& out, const PcsEnergyParameterManager &me);


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

	PcsEnergyParameterManager();

	~PcsEnergyParameterManager();

	PcsEnergyParameterManager(PcsEnergyParameterManager const & other);

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PcsEnergyParameterManager * create_singleton_instance();

	PcsEnergyParameterManager&
	operator=( PcsEnergyParameterManager const & other );

	static PcsEnergyParameterManager * instance_;
	utility::vector1<PcsEnergyParameter> pcs_e_p_all_;
	utility::vector1<std::string> vec_filename_all_;
	utility::vector1<core::Real> vec_individual_weight_all_;

public:

	/// @brief Re init the singleton
	void
	re_init();

	/// @brief Give me the number of paramagnetic center
	core::Size
	get_n_multi_data() const;

	/// @brief Add a new paramagnetic center
	void
	incremente_n_multi_data();

	/// @brief Give me the PcsEnergyParameter number i_multi_data
	PcsEnergyParameter &
	get_PcsEnergyParameter_for(core::Size i_multi_data);

};

} //PCS
} //methods
} //scoring
} //core

#endif
