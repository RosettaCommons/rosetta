// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_Relax_a_CDR_FullAtom.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_Ab_Relax_a_CDR_FullAtom_hh
#define INCLUDED_protocols_antibody2_Ab_Relax_a_CDR_FullAtom_hh






#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/Ab_Info.hh>
#include <core/fragment/FragSet.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


#include <protocols/antibody2/Ab_Relax_a_CDR_FullAtom.fwd.hh>




namespace protocols {
namespace antibody2 {
        
class Ab_Relax_a_CDR_FullAtom: public moves::Mover {
            
            
public:
    
    /// @brief default constructor
	Ab_Relax_a_CDR_FullAtom();
    
	/// @brief constructor with arguments
    Ab_Relax_a_CDR_FullAtom(bool current_loop_is_H3, bool H3_filter);
    
    /// @brief constructor with arguments
    Ab_Relax_a_CDR_FullAtom(bool current_loop_is_H3, bool H3_filter, bool is_camelid);
    


        
    virtual protocols::moves::MoverOP clone() const;
    
	/// @brief default destructor
	~Ab_Relax_a_CDR_FullAtom();
    
    void pass_start_pose(core::pose::Pose & start_pose);

    virtual void apply( core::pose::Pose & pose );
    
    virtual std::string get_name() const;
    
    
private:

    Ab_Info ab_info_;
    // Ab_InfoOP
    
    bool user_defined_;
    bool benchmark_;
    bool is_camelid_;
    loops::Loops all_loops_; 
    
    
    void set_default();
    void init();
    void setup_objects();
    
    

    
    /// @brief Build fullatom mode CDR H3 loop
	void build_fullatom_loop( core::pose::Pose & pose );
    
	void loop_fa_relax(
                       core::pose::Pose & pose_in,
                       core::Size const loop_begin,
                       core::Size const loop_end
                       );

	utility::vector1< core::fragment::FragSetOP > cdr_h3_frags_;
       
    // score functions
	core::scoring::ScoreFunctionOP highres_scorefxn_;
    
    /// @brief Fullatom mode loop building
    bool apply_fullatom_mode_;
    
	/// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;
    
	/// @brief actually enables H3 filter for H3 operations
	bool H3_filter_;
    
	/// @brief build H3 only
	bool antibody_build_;
    
	/// @brief refine H3 only
	bool antibody_refine_;
    
	/// @brief lower amplitude during base relaxation
	bool min_base_relax_;
    
	/// @brief use random cutpoints for h3 modeling
	bool h3_random_cut_;
    
	/// @brief cutpoint whose separation is computed in scorefile
	Size decoy_loop_cutpoint_;

 	/// @brief number of flanking residues:default 5
	core::Size h3_flank_;

	/// @brief relax flanking regions of h3
	bool flank_relax_;
    
	/// @brief freeze h3 during all cdr relax and local refine
	bool freeze_h3_;

    /// @brief enable docking local refine of LH chains & simultaneous H3 min
	bool snug_fit_;
    
    
    core::pose::Pose start_pose_;
    
    antibody2::Ab_Info antibody_in_;
    
	//packer task
	core::pack::task::TaskFactoryOP tf_;

    /// @brief just refine input loop
	bool refine_input_loop_;

	core::Size max_cycle_;
    
    core::Size base_;

    

};
    
    
    
    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif








