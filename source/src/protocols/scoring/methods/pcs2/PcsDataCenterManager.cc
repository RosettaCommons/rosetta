// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsDataCenterManager.cc
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
#include <protocols/scoring/methods/pcs2/PcsDataCenterManager.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers

#include <utility/vector1.hh>


// Numeric headers

// Objexx headers

// C++ headers

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static thread_local basic::Tracer TR_PcsDataCenterManager( "protocols.scoring.methods.pcs.PcsDataCenterManager" );

PcsDataCenterManager::PcsDataCenterManager(){

	// TR_PcsDataCenterManager << "Empty constructor called" << std::endl;
	// utility_exit_with_message( "You shouldn't call the empty constructor for PcsDataCenterManager class" );
}

PcsDataCenterManager::~PcsDataCenterManager(){
}

PcsDataCenterManager &
PcsDataCenterManager::operator=( PcsDataCenterManager const &other )
{

	// TR_PcsDataCenterManager << " = called" << std::endl;

	if ( this != &other ) {
		PCS_data_all_ = other.PCS_data_all_;
	}
	return *this;
}

PcsDataCenterManager::PcsDataCenterManager(PcsDataCenterManager const &other):
	CacheableData()
{
	// TR_PcsDataCenterManager << " () called" << std::endl;
	PCS_data_all_ = other.PCS_data_all_;
}

utility::vector1<PcsDataCenter> &
PcsDataCenterManager::get_PCS_data_all() {
	return (PCS_data_all_);
}

basic::datacache::CacheableDataOP
PcsDataCenterManager::clone() const {
	// TR_PcsDataCenterManager << "clone called" << std::endl;
	return basic::datacache::CacheableDataOP( new PcsDataCenterManager( *this ) );
}

std::ostream &
operator<<(std::ostream& out, const PcsDataCenterManager & m){
	core::Size i;

	out << "n paramagnetic center: " << m.get_n_multi_data() << std::endl;
	for ( i = 1 ; i <= m.get_n_multi_data(); ++i ) {
		out << m.PCS_data_all_[i] << std::endl;
	}
	return out;
}

core::Size
PcsDataCenterManager::get_n_multi_data() const{
	return (PCS_data_all_.size());
}

}//namespcacs PCS
}//namespace methods
}//namespace scoring
}//namespace protocols
