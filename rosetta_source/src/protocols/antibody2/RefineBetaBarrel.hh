// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineBetaBarrel.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_RefineBetaBarrel_hh
#define INCLUDED_protocols_antibody2_RefineBetaBarrel_hh

#include <protocols/antibody2/RefineBetaBarrel.fwd.hh>
#include <protocols/moves/Mover.hh>



namespace protocols {
namespace antibody2 {
        
        
        
class RefineBetaBarrel: public moves::Mover {
            
            
    public:
        /// @brief default constructor
        RefineBetaBarrel();
        
        /// @brief default destructor
        ~RefineBetaBarrel();
        
        void set_default();
        virtual void apply( core::pose::Pose & pose_in );
        virtual std::string get_name() const;
            
    private:
        bool user_defined_;
        
        
};
    
    
}// namespace antibody2
}// namespace protocols




#endif


