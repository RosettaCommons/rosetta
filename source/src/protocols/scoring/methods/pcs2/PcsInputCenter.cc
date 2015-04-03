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
 /// @file protocols/scoring/methods/pcs2/PcsInputCenter.cc
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
#include <protocols/scoring/methods/pcs2/PcsInputCenter.hh>
#include <protocols/scoring/methods/pcs2/PcsInputFile.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// Numeric headers

// Objexx headers

// C++ headers

#include <utility/vector1.hh>


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

static thread_local basic::Tracer TR_PcsInputCenter( "protocols.scoring.methods.pcs.PcsInputCenter" );

PcsInputCenter::PcsInputCenter(){
	utility_exit_with_message( "You shouldn't call the empty constructor for PcsInputCenter class" );
}

PcsInputCenter::~PcsInputCenter(){
}

PcsInputCenter::PcsInputCenter(PcsInputCenter const & other):
ReferenceCount()
{
	//	TR_PcsInputCenter << " () called" << std::endl;
	PcsInputFile_all_ = other.PcsInputFile_all_;
}

PcsInputCenter &
PcsInputCenter::operator=( PcsInputCenter const & other ){
	//	TR_PcsInputCenter << " = called" << std::endl;
	if ( this != &other ) {
		PcsInputFile_all_ = other.PcsInputFile_all_;
	}
	return *this;
}

std::map< std::string, PcsInputFile > &
PcsInputCenter::get_PcsInputFile_all(){
	return  PcsInputFile_all_;
}

	PcsInputCenter::PcsInputCenter(utility::vector1<std::string> const & filenames, utility::vector1<core::Real> const & weight){
	//	TR_PcsInputCenter << " constructor called" << std::endl;
	core::Real weight_sum;
	core::Size i;

	weight_sum = 0;
	for (i = 1; i <= filenames.size(); i++){
		weight_sum += weight[i];
	}

	for (i = 1; i <= filenames.size(); i++){
		//core::Real my_weight(weight[i]/weight_sum);
		//TODO correct the weighting scheme. For the moment it is one automatically
		core::Real my_weight(weight[i]);
		PcsInputFile pcs_i_f_temp(filenames[i], my_weight);
		PcsInputFile_all_.insert ( std::pair< std::string, PcsInputFile >(filenames[i], pcs_i_f_temp) );
	}
}

std::ostream &
 operator<<(std::ostream & out,  const PcsInputCenter &me ){

		std::map< std::string, PcsInputFile >::iterator it;
		std::map< std::string, PcsInputFile > mymap;
		mymap = me.PcsInputFile_all_;

		for ( it = mymap.begin(); it != mymap.end(); ++it ) {
			out << "For the file '" << it->first << std::endl;
			out << it->second;
			out << "The relative weight is " << (it->second).get_weight() << std::endl;
		}
		return out;
}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
