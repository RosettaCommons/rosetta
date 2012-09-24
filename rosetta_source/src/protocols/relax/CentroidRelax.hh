// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/relax/CentroidRelax.hh
/// @brief Centroid-based relax, using Frank Dimaio's updated centroid stats
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_relax_CENTROIDRELAX_HH
#define	INCLUDED_protocols_relax_CENTROIDRELAX_HH


//Core
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/kinematics/MoveMap.hh>
//Protocols

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/moves/Mover.hh>

//Basic
#include <utility/vector1.hh>

//Unoptimized, not completely tested class.  In development.
namespace protocols{
    namespace relax{
        using namespace core::scoring;
        using namespace protocols::simple_moves;
        using namespace core::kinematics;
        using namespace core::pack::task;
        using core::Size;
        using std::string;
        using core::Real;
        using core::pose::Pose;
        using utility::vector1;
           
        
class CentroidRelax : public RelaxProtocolBase {
public:
    CentroidRelax();
    CentroidRelax(MoveMapOP mm);
    CentroidRelax(MoveMapOP mm, ScoreFunctionOP cen_scorefxn_in);
        
        //Specific functions for Mover...
    virtual ~CentroidRelax();
    
    virtual string get_name() const;
        //virtual bool reinitialize_for_new_input() const;
        //virtual bool reinitialize_for_each_job() const;
    virtual protocols::moves::MoverOP clone() const;
    
    //Yes, I like options...Though their interactions are getting confusing...
    
    void set_rounds(Size rounds);
    ///@brief sets to use Rama2b instead of Rama
    void set_use_rama2b(bool use);
    void set_ramp_rama(bool use);
    void set_ramp_vdw(bool use);
    void set_score_function(ScoreFunctionOP cen_score);
    void set_fa_score_function(ScoreFunctionOP fa_score);
    //void set_packer_task(PackerTaskOP task);
    
    void set_min_type(string min);
    void set_cartesian(bool cart);
    void set_final_min_fullatom(bool min);
    void set_final_min_sidechains(bool min);
    void set_final_repack_sidechains(bool repack);
    
    ///@brief Load the default parameters from the default file.
    void set_default_parameters();
    ///@brief Nessessary if you pass a centroid rep and centroid only mode is false.
    void set_recover_sidechains_mover(ReturnSidechainMoverOP recover_sc);
    
    
    ///@brief Switch to full atom mode during final step, doing any repacking/minimization. 
    void set_return_full_atom_from_centroid(bool fa_ret = false);
    ///@brief Sets the Mover to only do centroid-based operations, returning a centroid representation
    void set_centroid_only_mode(bool cen_only=false);
    
    ///@brief controls whether to repack side chains intermittently during protocol(Using that the fa_scorefxn_ score). If centroid only mode is false
    void set_repack_intermittent(bool use);
    
    void set_defaults();
    
    ///@brief Makes DNA Rigid exactly like FastRelax.
    void makeDnaRigid(Pose & pose);
    ///@brief Applies the protocol, See notes
    virtual void apply( Pose & pose );
    //Notes:
        //Pass centroid with default settings: Will give error
        //Pass centroid with set_final_min_sidechains and set_final_repack_sidechains false: Will run and return Centroid
        //Pass centroid + set_centroid_only mode: Pose returned is centroid
        //Pass centroid + set recover_sidechains_mover: Pose returned is centroid
        //Pass centroid + set recover_sidechains_mover + set fa_return: Pose returned is FA
        
        //Pass full_atom: will return full_atom structure.
        //Pass full_atom + set_centroid_only mode: Pose returned is centroid
    
        //Setting ramp_rama and ramp_vdw to false switches to the BASIC protocol which is rounds of centroid minmover
    
    
    
private:
    void initialize_objects();
    //@brief used internally to setup extra stuff for movemap.
    void setup_class_movemap(Pose & pose);
    void movemap_to_taskfactory(Pose & pose, PackerTaskOP task);
    
    bool ramp_rama_;
    bool ramp_vdw_;
    Size rounds_;
    bool repack_intermittent_; //Do an intermittent repack of the SC true in movemap
    bool cartesian_mode_;
    bool final_min_fa_; //Minimize full atom structure @end, before repack and sc_min
    bool final_min_sc_; //Minimize Sidechains after repack @end.
    bool final_repack_sc_;  //Rpack at end?
    bool centroid_only_mode_; //Never switch to FA for anything
    bool fa_return_; //Return FA if recover_sc is set
    bool constrain_to_original_coords_;
    
    ScoreFunctionOP scorefxn_;
    ScoreFunctionOP fa_scorefxn_;
    SwitchResidueTypeSetMoverOP to_fa_;
    SwitchResidueTypeSetMoverOP to_cen_;
    ReturnSidechainMoverOP recover_sc_;
    bool sc_recoverable_;
    MoveMapOP movemap_;
    ScoreType rama_type_;
    
    struct parameters{
        vector1< Real > vdw_params;
        vector1< Real > rama_params;
        vector1< Real > min_params;
        vector1< Real > cst_params;
        };
        
    parameters def_parameters;    
        
};
    }
}

#endif	//#ifndef INCLUDED_protocols/relax_CENTROIDRELAX_HH
