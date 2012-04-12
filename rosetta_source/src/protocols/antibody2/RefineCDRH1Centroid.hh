// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineCDRH1Centroid.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_RefineCDRH1Centroid_hh
#define INCLUDED_protocols_antibody2_RefineCDRH1Centroid_hh


#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/loops/Loops.hh>

#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>
#include <protocols/moves/ChangeFoldTreeMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/RefineCDRH1Centroid.fwd.hh>

#ifdef PYROSETTA
	#include <protocols/moves/PyMolMover.hh>
#endif

using namespace core;
namespace protocols {
namespace antibody2 {

class RefineCDRH1Centroid: public moves::Mover {


public:

    /// @brief default constructor
	RefineCDRH1Centroid();

    /// @brief constructor with arguments
    RefineCDRH1Centroid(AntibodyInfoOP antibody_info, std::string loop_name);

    /// @brief constructor with arguments
    RefineCDRH1Centroid( loops::LoopOP a_cdr_loop);

    /// @brief default destructor
	~RefineCDRH1Centroid();
    
	void set_default();
	virtual void apply( pose::Pose & pose_in );
    virtual std::string get_name() const;



private:
    void init(loops::LoopOP a_cdr_loop);

    void finalize_setup( core::pose::Pose & pose );
	void loop_centroid_relax(
                             pose::Pose & pose_in,
                             Size const loop_begin,
                             Size const loop_end );
    
    loops::LoopOP the_cdr_loop_;
    
    bool benchmark_;
    scoring::ScoreFunctionOP lowres_scorefxn_;
	/// @brief refine H3 only
	bool antibody_refine_;
    
    /// @brief enable docking local refine of LH chains & simultaneous H3 min
	bool snug_fit_;
    
    /// @brief just refine input loop
	bool refine_input_loop_;
    
    bool use_pymol_diy_;

};







} // namespace antibody2
} // namespace protocols

#endif








