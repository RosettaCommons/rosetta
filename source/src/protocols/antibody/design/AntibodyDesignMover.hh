// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.hh
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh

#include <protocols/antibody/design/AntibodyDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>
#include <protocols/antibody/design/AntibodyCDRDesigner.hh>
#include <protocols/antibody/AntibodyInfo.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>
namespace protocols {
namespace antibody {
namespace design {

using namespace protocols::antibody;	
using std::string;
using utility::vector1;


///@brief Main AntibodyDesignMover, tieing together individual movers and classes.  Main mover for application.
///
class AntibodyDesignMover : public protocols::moves::Mover {
	
public:

	AntibodyDesignMover();

	virtual ~AntibodyDesignMover();
	
	virtual string 
	get_name() const;
	
	virtual protocols::moves::MoverOP clone() const;

	
	////////////////////////////////////////////////////////////////////////////
	// Optional Custom Settings
	//
	//
	
	void
	set_graft_designer(AntibodyGraftDesignerOP graft_designer);
	
	void
	set_sequence_designer(AntibodyCDRDesignerOP seq_designer);
	
	void
	set_modeler(AntibodyDesignModelerOP modeler);
	
	void
	set_scorefxn(ScoreFunctionOP scorefxn);
	
	void
	set_ab_info(AntibodyInfoOP ab_info);

	
	////////////////////////////////////////////////////////////////////////////
	// Algorithm Settings
	//
	//
	
	///@brief Use low-resolution graft designer for structural sampling.  Default true.
	void
	set_use_graft_designer(bool setting);
	
	///@brief Use high-resolution sequence designer for sequence sampling. Default true.
	void
	set_use_sequence_designer(bool setting);
	
	///@brief Run dock/min modeling step after the graft design step if run.
	void
	set_do_post_graft_design_modeling(bool setting);
	
	///@brief Run dock/min modeling step after sequence design if run.
	void
	set_do_post_design_modeling(bool setting);
	
	////////////////////////////////////////////////////////////////////////////
	// Main Method.  All options read from cmd-line or set through individual classes.
	//
	//
	
	virtual void 
	apply( core::pose::Pose & pose );
	
private:
	
	///@brief Post-graft step modeling.  If no graft step, no need to do post-graft modeling.  Default false.
	void model_post_graft(core::pose::Pose & pose);
	
	///@brief Post-design step modeling.  Less aggressive, more high resolution.  Default false.
	void model_post_design(core::pose::Pose & pose);
	
	void read_options();
	
	void init_on_new_input( core::pose::Pose const & pose );

	AntibodyDatabaseManagerOP cdr_db_parser_;
	
	AntibodyGraftDesignerOP graft_designer_;
	AntibodyCDRDesignerOP seq_designer_;
	AntibodyDesignModelerOP modeler_;
	
	AntibodyInfoOP ab_info_;
	
	ScoreFunctionOP scorefxn_;
	
	bool run_graft_designer_;
	bool run_sequence_designer_;
	bool run_post_graft_modeling_;
	bool run_post_design_modeling_;
	
};
} //design
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
