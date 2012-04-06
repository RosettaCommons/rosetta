// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineBetaBarrel.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/RefineBetaBarrel.hh>

#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.antibody2.RefineBetaBarrel");

using namespace core;

namespace protocols {
namespace antibody2 {
        
        
        
        
// default constructor
RefineBetaBarrel::RefineBetaBarrel() : Mover() {
    user_defined_ = false;
}

RefineBetaBarrel::~RefineBetaBarrel() {}
    
    
void RefineBetaBarrel::set_default(){
}
    
    
std::string RefineBetaBarrel::get_name() const {
    return "RefineBetaBarrel";
}
    
    
    
void RefineBetaBarrel::apply( pose::Pose & pose ) {

}



}// namespace antibody2
}// namespace protocols





