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
 /// @file protocols/scoring/methods/pcs2/PcsGridSearchParameter.cc
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


// Unit Headers
#include <protocols/scoring/methods/pcs2/PcsGridSearchParameterManager.hh>

// Package Headers
//#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>
//#include <protocols/scoring/methods/pcs2/PcsInputCenterManager.hh>

// Project Headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// C++ headers
#include <iostream>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static thread_local basic::Tracer TR_PcsGridSearchParameterManager( "protocols.scoring.methods.pcs.PcsGridSearchParameterManager" );

// Singleton instance and mutex static data members
namespace utility {

using protocols::scoring::methods::pcs2::PcsGridSearchParameterManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< PcsGridSearchParameterManager > ::singleton_mutex_;
template <> std::atomic< PcsGridSearchParameterManager * > utility::SingletonBase< PcsGridSearchParameterManager >::instance_( 0 );
#else
template <> PcsGridSearchParameterManager * utility::SingletonBase< PcsGridSearchParameterManager >::instance_( 0 );
#endif

}

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

PcsGridSearchParameterManager *
PcsGridSearchParameterManager::create_singleton_instance()
{
	return new PcsGridSearchParameterManager;
}


PcsGridSearchParameterManager::PcsGridSearchParameterManager() {}

/// @brief WARNING WARNING WARNING! THREAD UNSAFE!
void
PcsGridSearchParameterManager::re_init(){

	core::Size n(grid_s_p_all_.size());
	core::Size i(1);


	for(i = 1; i <= n; ++i){
		grid_s_p_all_.pop_back();
	}
	std::cerr <<"CHECKING zz 0 = " << grid_s_p_all_.size() << std::endl ;

}

core::Size
PcsGridSearchParameterManager::get_n_multi_data() const{
	return(grid_s_p_all_.size());
}

/// @brief WARNING WARNING WARNING! THREAD UNSAFE!
void
PcsGridSearchParameterManager::incremente_n_multi_data(){
	PcsGridSearchParameter g;
	grid_s_p_all_.push_back(g);
}

/// @brief WARNING WARNING WARNING! THREAD UNSAFE!
/// THIS CODE COULD INSTEAD RETURN A CONST REFERENCE
PcsGridSearchParameter &
PcsGridSearchParameterManager::get_grid_search_parameters(core::Size i_multi_data){

	core::Size n_multi_data(get_n_multi_data());
	if ( (i_multi_data > n_multi_data) ){
		TR_PcsGridSearchParameterManager << "Problem in get_grid_search_parameters: i_multi_data = "<<i_multi_data<<" and n_multi_data = "<<n_multi_data << std::endl;
		utility_exit_with_message("There is a coding problem");
	}

	return(grid_s_p_all_[i_multi_data]);
}




}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

