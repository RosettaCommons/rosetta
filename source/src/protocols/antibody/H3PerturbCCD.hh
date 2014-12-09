// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/H3PerturbCCD.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody_H3PerturbCCD_hh
#define INCLUDED_protocols_antibody_H3PerturbCCD_hh


#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/fragment/FragSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/H3PerturbCCD.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>


using namespace core;
namespace protocols {
namespace antibody {



class H3PerturbCCD: public moves::Mover {


public:

	/// @brief default constructor
	H3PerturbCCD();

	/// @brief constructor with arguments
	H3PerturbCCD(AntibodyInfoOP antibody_in);

	/// @brief constructor with arguments
	H3PerturbCCD(AntibodyInfoOP antibody_in,
	             core::scoring::ScoreFunctionCOP highres_scorefxn);

	virtual protocols::moves::MoverOP clone() const;

	/// @brief default destructor
	~H3PerturbCCD();



	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const {
		return "H3PerturbCCD";
	}

	void read_and_store_fragments();


	void pass_the_loop(loops::Loop & input_loop) {
		input_loop_ = input_loop;
	}

	void turn_off_H3_filter() {
		H3_filter_ = false;
	}

	/// @brief set scorefunction for low resolution of CDR H3 modeling
	void set_lowres_score_func(scoring::ScoreFunctionCOP lowres_scorefxn ) {
		lowres_scorefxn_ = lowres_scorefxn->clone();
	}

private:

	AntibodyInfoOP ab_info_;

	loops::Loop input_loop_;

	bool user_defined_;
	bool is_camelid_;



	void set_default();
	void init();

	void finalize_setup( core::pose::Pose & pose );


	/// @brief Build centroid mode CDR H3 loop
	void build_centroid_loop( core::pose::Pose & pose );

	void scored_frag_close(
	    core::pose::Pose & pose_in,
	    loops::Loop const trimmed_cdr_h3 );



	/// @brief size of loop above which 9mer frags are used
	core::Size cutoff_9_; // default 16

	/// @brief size of loop above which 3mer frags are used
	core::Size cutoff_3_; // default 6

	core::scoring::ScoreFunctionOP lowres_scorefxn_;

	core::Real cen_cst_;

	/// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;

	/// @brief actually enables H3 filter for H3 operations
	bool H3_filter_;

	core::Size num_cycles1_;
	core::Size max_ccd_cycles_;

	utility::vector1< core::fragment::FragSetOP > cdr_h3_frags_;

	core::Real h3_fraction_;
	core::Real ccd_threshold_;
	core::Real Temperature_;

	protocols::moves::MonteCarloOP mc_, outer_mc_;

};



} // namespace antibody
} // namespace protocols

#endif




