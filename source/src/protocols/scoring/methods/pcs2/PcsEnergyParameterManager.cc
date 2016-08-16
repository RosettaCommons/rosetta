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
/// @file protocols/scoring/methods/pcs2/PcsEnergyParameterManager.cc
///
/// @brief
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


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.hh>

// Package headers
//#include <protocols/scoring/methods/pcs2/PcsGridSearchParameterManager.hh>
//#include <protocols/scoring/methods/pcs2/PcsInputCenterManager.hh>

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <iostream>

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::scoring::methods::pcs2::PcsEnergyParameterManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< PcsEnergyParameterManager >::singleton_mutex_{};
template <> std::atomic< PcsEnergyParameterManager * > utility::SingletonBase< PcsEnergyParameterManager >::instance_( 0 );
#else
template <> PcsEnergyParameterManager * utility::SingletonBase< PcsEnergyParameterManager >::instance_( 0 );
#endif

}

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static THREAD_LOCAL basic::Tracer TR_PcsEnergyParameterManager( "protocols.scoring.methods.pcs.PcsEnergyParameterManager" );

PcsEnergyParameterManager *
PcsEnergyParameterManager::create_singleton_instance()
{
	return new PcsEnergyParameterManager;
}

PcsEnergyParameterManager::PcsEnergyParameterManager(){

}

PcsEnergyParameterManager::~PcsEnergyParameterManager(){
}

void
PcsEnergyParameterManager::re_init(){
	core::Size i(1);
	core::Size n(pcs_e_p_all_.size());

	vec_filename_all_.clear();
	std::cerr <<"CHECKING a 0 = " << vec_filename_all_.size() << std::endl;
	vec_individual_weight_all_.clear();
	std::cerr <<"CHECKING b 0 = " << vec_individual_weight_all_.size() << std::endl;

	for ( i = 1; i <= n; ++i ) {
		pcs_e_p_all_.pop_back();
	}

	std::cerr <<"CHECKING c 0 = " << pcs_e_p_all_.size() << std::endl;
}

// singletons shouldn't have copy constructors PcsEnergyParameterManager::PcsEnergyParameterManager(PcsEnergyParameterManager const & other)
// singletons shouldn't have copy constructors {
// singletons shouldn't have copy constructors  pcs_e_p_all_ = other.pcs_e_p_all_;
// singletons shouldn't have copy constructors  vec_filename_all_ = other.vec_filename_all_;
// singletons shouldn't have copy constructors  vec_individual_weight_all_ = other.vec_individual_weight_all_;
// singletons shouldn't have copy constructors }

// singletons shouldn't have assignment operators PcsEnergyParameterManager&
// singletons shouldn't have assignment operators PcsEnergyParameterManager::operator=( PcsEnergyParameterManager const & other ){
// singletons shouldn't have assignment operators  if ( this != &other ) {
// singletons shouldn't have assignment operators   instance_ = other.instance_;
// singletons shouldn't have assignment operators   pcs_e_p_all_ = other.pcs_e_p_all_;
// singletons shouldn't have assignment operators   vec_filename_all_ = other.vec_filename_all_;
// singletons shouldn't have assignment operators   vec_individual_weight_all_ = other.vec_individual_weight_all_;
// singletons shouldn't have assignment operators  }
// singletons shouldn't have assignment operators  return *this;
// singletons shouldn't have assignment operators }

std::ostream &
operator << ( std::ostream& out, const PcsEnergyParameterManager &me ){
	core::Size i;

	out << "************************ Printing PcsEnergyParameterManager" << std::endl;

	out << "vec_filename_all_: ";
	for ( i = 1; i <= me.vec_filename_all_.size(); ++i ) {
		out << me.vec_filename_all_[i]<< " ";
	}
	out << std::endl;

	out << "vec_individual_weight_all_: ";
	for ( i = 1; i <= me.vec_individual_weight_all_.size(); ++i ) {
		out << me.vec_individual_weight_all_[i]<< " ";
	}
	out << std::endl;

	out << "pcs_e_p_all_: ";
	for ( i = 1; i <= me.pcs_e_p_all_.size(); ++i ) {
		out << me.pcs_e_p_all_[i] << std::endl;
	}
	out << std::endl;

	return(out);

}

core::Size
PcsEnergyParameterManager::get_n_multi_data() const{
	return(pcs_e_p_all_.size());
}

void
PcsEnergyParameterManager::incremente_n_multi_data(){
	PcsEnergyParameter pcs_e_p;
	pcs_e_p_all_.push_back(pcs_e_p);
}

PcsEnergyParameter &
PcsEnergyParameterManager::get_PcsEnergyParameter_for(core::Size i_multi_data){

	if ( (i_multi_data > get_n_multi_data()) ) {
		TR_PcsEnergyParameterManager << "Problem in get_PcsEnergyParameter_for: i_multi_data = "<<i_multi_data<<" and n_multi_data = "<< get_n_multi_data() << std::endl;
		utility_exit_with_message("There is a coding problem in get_PcsEnergyParameter_for");
	}
	return(pcs_e_p_all_[i_multi_data]);
}

} // PCS
} // methods
} // scoring
} // core
