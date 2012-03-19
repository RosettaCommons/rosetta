// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/Ab_Info.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef INCLUDED_protocols_antibody2_Ab_Info_hh
#define INCLUDED_protocols_antibody2_Ab_Info_hh

#include <protocols/antibody2/Ab_Info.fwd.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
// C++ Headers

// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace antibody2 {

/// antibody2 definition
class Ab_Info : public utility::pointer::ReferenceCount {

public:
	typedef std::map < std::string, loops::LoopOP > LoopMap;
    typedef std::map <std::string, core::Size> CDR_Numbering_Begin_Map;
    typedef std::map <std::string, core::Size> CDR_Numbering_End_Map;

    
	/// default constructor
	Ab_Info();

	/// constructor with arguments
	Ab_Info( core::pose::Pose & pose );
	Ab_Info( core::pose::Pose & pose, bool camelid );
	Ab_Info( core::pose::Pose & pose, std::string cdr_name );

	void setup_CDR_loops( core::pose::Pose & pose, bool camelid );

	void all_cdr_fold_tree( core::pose::Pose & pose );
//	void cdr_h3_fold_tree( core::pose::Pose & pose );

	/// @brief return the loop of a certain loop type
	loops::LoopOP get_CDR_loop( std::string loop );

	// return kinked/extended
	bool is_kinked()   { return kinked_H3_;   }
	bool is_extended() { return extended_H3_; }
    bool is_camelid()  { return is_camelid_;     }
    utility::vector1<char> get_Fv_sequence() { return Fv_sequence_;}

	/// align current Fv to native.Fv
	void align_to_native( core::pose::Pose & pose, 
                          antibody2::Ab_Info & native, 
                          core::pose::Pose & native_pose );
    
    void load_CDR_query_info_to_check();

    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const Ab_Info & ab_info );

    
    
    // Start coordinates of active loop
	core::Size current_start;
	// End coordinates of active loop
	core::Size current_end;

    
	loops::Loops all_cdr_loops_;
    
    void detect_and_set_CDR_H3_stem_type( core::pose::Pose & pose );
	void detect_and_set_camelid_CDR_H3_stem_type();
	void detect_and_set_regular_CDR_H3_stem_type( core::pose::Pose & pose );
    void get_CDRs_numbering();
	void set_default( bool camelid );
    
    //bool is_my_pose_antibody(core::pose::Pose & pose);
    //bool is_my_antibody_camelid(core::pose::Pose & pose);


    
private:
	// cdr loops
	LoopMap loops_;
	loops::LoopOP L1_, L2_, L3_;
    loops::LoopOP H1_, H2_, H3_;
    std::string L1_seq_, L2_seq_, L3_seq_;
    std::string H1_seq_,H2_seq_,H3_seq_;

    
	core::Size hfr_[7][3]; // array of framework residues for alignment

    CDR_Numbering_Begin_Map CDR_numbering_begin_;
    CDR_Numbering_End_Map CDR_numbering_end_;

    
	bool is_camelid_;
	bool kinked_H3_;
	bool extended_H3_;
    
    utility::vector1< char > Fv_sequence_;
 
    
    
    

};


} //namespace antibody2
} //namespace protocols


#endif //INCLUDED_protocols_loops_Ab_Info_HH




