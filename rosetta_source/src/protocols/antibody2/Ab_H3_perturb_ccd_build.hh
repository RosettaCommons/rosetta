// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_H3_perturb_ccd_build.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_Ab_H3_perturb_ccd_build_hh
#define INCLUDED_protocols_antibody2_Ab_H3_perturb_ccd_build_hh






#include <protocols/moves/Mover.hh>

#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/antibody2/Ab_Info.fwd.hh>
#include <protocols/antibody2/Ab_H3_perturb_ccd_build.fwd.hh>
#include <protocols/antibody2/Ab_H3_cter_insert_mover.fwd.hh>





namespace protocols {
namespace antibody2 {
        
    
    
class Ab_H3_perturb_ccd_build: public moves::Mover {
            
            
public:
    
    /// @brief default constructor
	Ab_H3_perturb_ccd_build();
    
	/// @brief constructor with arguments
    Ab_H3_perturb_ccd_build(
                            bool current_loop_is_H3, 
                            bool H3_filter, 
                            bool is_camelid, 
                            Ab_InfoOP & antibody_in);
    
        
    virtual protocols::moves::MoverOP clone() const;
    
	/// @brief default destructor
	~Ab_H3_perturb_ccd_build();
    


    virtual void apply( core::pose::Pose & pose );
    
    virtual std::string get_name() const;
    void read_and_store_fragments( core::pose::Pose & pose );
    
    
private:

    Ab_InfoOP ab_info_;
    

    
    bool user_defined_;
    bool is_camelid_;
    
    
    void set_default();
    void init( bool current_loop_is_H3, bool H3_filter,bool is_camelid, Ab_InfoOP & antibody_in);
    void setup_objects();
    void finalize_setup( core::pose::Pose & pose );

    
    /// @brief Build centroid mode CDR H3 loop
	void build_centroid_loop( core::pose::Pose & pose );
    
    void scored_frag_close(
                           core::pose::Pose & pose_in,
                           loops::Loop const trimmed_cdr_h3 );

    
    core::Size max_cycle_;
    
    /// @brief Number of ADDITIONAL residues modeled from H3_CTERM
	///        These residues range from H:n-2,n-1,n,n+1 of H3
	core::Size c_ter_stem_;
    
    
    /// @brief size of loop above which 9mer frags are used
	core::Size cutoff_9_; // default 16
    
	/// @brief size of loop above which 3mer frags are used
	core::Size cutoff_3_; // default 6
    
    core::scoring::ScoreFunctionOP lowres_scorefxn_;
    
    Ab_H3_cter_insert_moverOP ab_h3_cter_insert_mover_;
    
    core::Real cen_cst_;
    core::pose::Pose hfr_pose_;

    
    
    /// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;
    /// @brief actually enables H3 filter for H3 operations
	bool H3_filter_;
    
    utility::vector1< core::fragment::FragSetOP > cdr_h3_frags_;
    
};
    
    
    
} // namespace antibody2
} // namespace protocols

#endif




