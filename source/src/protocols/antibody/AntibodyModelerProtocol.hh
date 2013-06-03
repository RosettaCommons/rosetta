// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody/AntibodyModelerProtocol.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody_AntibodyModelerProtocol_hh
#define INCLUDED_protocols_antibody_AntibodyModelerProtocol_hh

#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyModelerProtocol.fwd.hh>
#include <protocols/simple_moves/ConstraintSetMover.fwd.hh>



using namespace core;
namespace protocols {
namespace antibody {

class AntibodyModelerProtocol: public moves::Mover {
public:

	// default constructor
	AntibodyModelerProtocol();

	// default destructor
	~AntibodyModelerProtocol();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	virtual void apply( pose::Pose & pose );

	virtual std::string get_name() const;
    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    
	// simple inline setters
    void set_BenchMark(bool setting) {
        benchmark_ = setting;
    }
    void set_ModelH3(bool setting) {
        model_h3_ = setting; 
    }
	void set_SnugFit(bool setting) {
        snugfit_ = setting;  
    }
    void set_refine_h3(bool setting){
        refine_h3_ = setting;
    }
    void set_H3Filter(bool setting) {
        h3_filter_ = setting;
    }
    void set_CterInsert (bool setting) {
        cter_insert_ = setting;
    }
    void set_sc_min(bool setting) {
        sc_min_ = setting ;
    }
    void set_rt_min(bool setting) {
        rt_min_ = setting ;
    }
    void set_flank_residue_min (bool setting) {
        flank_residue_min_ = setting;
    }
    void set_packonly_after_graft (bool setting){
        packonly_after_graft_ = setting;
    }
    void set_perturb_type(std::string remodel) {
        h3_perturb_type_ = remodel;
    }
    void set_refine_type (std::string refine)  {
        h3_refine_type_ = refine;
    }
    void set_H3Filter_Tolerance(core::Size const number){
        h3_filter_tolerance_ = number;
    }
    void set_cst_weight ( core::Real const cst_weight){
        cst_weight_ = cst_weight;
    }
    void set_use_constraints(bool const use_csts){
        use_csts_ = use_csts;
    }
	void set_constrain_cter(bool const setting){
		constrain_cter_ = setting;
	}
	void set_constrain_vlvh_qq(bool const setting){
		constrain_vlvh_qq_ = setting;
	}
    void set_flank_residue_size(core::Real const flank_residue_size) {
        flank_residue_size_ = flank_residue_size;
    }
    void set_middle_pack_min( bool middle_pack_min) {
        middle_pack_min_ = middle_pack_min;
    }
    void set_bad_nter(bool setting){
        bad_nter_ = setting;
    }


	void display_constraint_residues( pose::Pose & pose );

	void show( std::ostream & out=std::cout ) const;
	friend std::ostream & operator<<(std::ostream& out, const AntibodyModelerProtocol & ab_m );


private:
    bool model_h3_;
    bool snugfit_;
    bool refine_h3_;
    bool h3_filter_;
    bool cter_insert_;
    bool LH_repulsive_ramp_;
    bool sc_min_;
    bool rt_min_;
    bool camelid_;
    bool camelid_constraints_;
    bool flank_residue_min_;
    bool flank_residue_size_;
    bool middle_pack_min_;
    bool packonly_after_graft_;
    std::string h3_perturb_type_;
    std::string h3_refine_type_;
    core::Real cen_cst_, high_cst_;
    bool bad_nter_;
    bool use_csts_;
	bool constrain_vlvh_qq_, constrain_cter_;
    
	// Benchmark mode for shorter_cycles
	bool benchmark_;

	bool user_defined_; // for constructor options passed to init

	bool flags_and_objects_are_in_sync_;
    core::Size h3_filter_tolerance_;


	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions    
    core::scoring::ScoreFunctionOP loop_scorefxn_highres_;
    core::scoring::ScoreFunctionOP loop_scorefxn_centroid_;
    core::scoring::ScoreFunctionOP dock_scorefxn_highres_;
    core::scoring::ScoreFunctionOP pack_scorefxn_;
    
    // constraint set mover
	protocols::simple_moves::ConstraintSetMoverOP cdr_constraint_;
    
	// external objects
	AntibodyInfoOP ab_info_;


	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///		or via the command line.
	void init();
    
	void setup_objects();
    
    /// @brief Output of various metrics of final model
    void echo_metrics_to_jd2( core::pose::Pose & pose, protocols::jd2::JobOP job );

}; // class
    
    
    
    
} // namespace antibody
} // namespace protocols

#endif

