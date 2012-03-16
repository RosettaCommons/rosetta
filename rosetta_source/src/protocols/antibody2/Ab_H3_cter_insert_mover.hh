// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_H3_cter_insert_mover.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )




#ifndef INCLUDED_protocols_antibody2_Ab_H3_cter_insert_mover_hh
#define INCLUDED_protocols_antibody2_Ab_H3_cter_insert_mover_hh


#include <protocols/antibody2/Ab_H3_cter_insert_mover.fwd.hh>
#include <core/fragment/FragData.hh>
#include <protocols/antibody2/Ab_Info.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>




using namespace core;
namespace protocols {
namespace antibody2 {

    
    
    
    
//////////////////////////////////////////////////////////////////////////
/// @brief H3 CDR, Fragment Insertion and CCD
/// @details
class Ab_H3_cter_insert_mover : public protocols::moves::Mover {
    
public:
    /// @brief default constructor
	Ab_H3_cter_insert_mover(antibody2::Ab_Info & ab_info);
    
	/// @brief constructor with arguments
	Ab_H3_cter_insert_mover(antibody2::Ab_Info & ab_info, bool camelid );
	
    
	/// @brief default destructor
	~Ab_H3_cter_insert_mover();
    
    void set_default();
    
	virtual void apply(pose::Pose & pose );
    virtual std::string get_name() const;

    
    
    // read CDR H3 C-terminal fragments (size: 4)
    void read_H3_cter_fragment(
                               antibody2::Ab_Info & ab_info,
                               bool is_camelid
                               );
    
    
private:
    
    // CDR H3 C-terminal fragments
	utility::vector1< core::fragment::FragData > H3_base_library_;
   
    Ab_InfoOP ab_info_;

    /// @brief insert C-terminal fragments
    void antibody_modeling_insert_ter( core::pose::Pose & pose) ;
    
    bool user_defined_;
        
    /// @brief benchmark flag
	bool benchmark_;
    
    /// @brief is camelid antibody without light chain
	bool is_camelid_;
    
        
    void init(antibody2::Ab_Info & ab_info, bool camelid, bool benchmark);
//    void setup_objects();
//    void finalize_setup( core::pose::Pose & pose );

    
    std::string H3_ter_library_filename_;
    
};
    
    
    
    
    

}//antibody2
}//protocols

#endif

