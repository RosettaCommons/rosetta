// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/GraftCDRLoopsProtocol.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody2_GraftCDRLoopsProtocol_hh
#define INCLUDED_protocols_antibody2_GraftCDRLoopsProtocol_hh


#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/Ab_TemplateInfo.fwd.hh>
#include <protocols/antibody2/GraftCDRLoopsProtocol.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>

#include <iostream>




namespace protocols {
namespace antibody2 {

class GraftCDRLoopsProtocol: public moves::Mover {
public:
    typedef std::map < std::string, bool > GraftMap;

	// default constructor
	GraftCDRLoopsProtocol();

	// default destructor
	~GraftCDRLoopsProtocol();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	virtual void apply( core::pose::Pose & pose );

	// simple inline setters
	void set_graft_l1( bool graft_l1 ) { graft_l1_ = graft_l1; }
	void set_graft_l2( bool graft_l2 ) { graft_l2_ = graft_l2; }
	void set_graft_l3( bool graft_l3 ) { graft_l3_ = graft_l3; }
	void set_graft_h1( bool graft_h1 ) { graft_h1_ = graft_h1; }
	void set_graft_h2( bool graft_h2 ) { graft_h2_ = graft_h2; }
	void set_graft_h3( bool graft_h3 ) { graft_h3_ = graft_h3; }
	void set_camelid( bool camelid ) { camelid_ = camelid; }
	void set_camelid_constraints( bool camelid_constraints ) { camelid_constraints_ = camelid_constraints; }
	void set_benchmark( bool benchmark ) { benchmark_ = benchmark; }
    void set_cst_weight(core::Real cst_weight){ cst_weight_ = cst_weight; }

	virtual std::string get_name() const;

    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    
	void display_constraint_residues( core::pose::Pose & pose );
    
    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const GraftCDRLoopsProtocol & ab_m_2 );

    
    
private:

	// Modeling H3 options
	bool graft_l1_;
	bool graft_l2_;
	bool graft_l3_;
	bool graft_h1_;
	bool graft_h2_;
	bool graft_h3_;
	bool camelid_;
	bool camelid_constraints_;

    
    GraftMap grafts_ ;
    

    
	// Benchmark mode for shorter_cycles
	bool benchmark_;

	bool user_defined_; // for constructor options passed to init

	// flag for one time fragment initialization
	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions
	core::scoring::ScoreFunctionOP scorefxn_;

	// external objects
	AntibodyInfoOP ab_info_;
    Ab_TemplateInfoOP ab_t_info_ ;


	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();
    void set_packer_default( core::pose::Pose & pose, bool include_current ) ;


    
    // movers
    protocols::moves::SequenceMoverOP graft_sequence_ ;
    protocols::simple_moves::PackRotamersMoverOP packer_ ;
    protocols::moves::PyMolMoverOP pymol_ ;


	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///		or via the command line.
	void init();

	void setup_objects();

}; // class GraftCDRLoopsProtocol

    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif





