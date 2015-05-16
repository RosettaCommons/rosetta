// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignProtocol.hh
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_hh

#include <protocols/antibody/design/AntibodyDesignProtocol.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyDesignMover.hh>

#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>

#include <protocols/antibody/constraints/CDRDihedralConstraintMover.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace antibody {
namespace design {

///@brief Main AntibodyDesignProtocol, tieing together individual movers and classes.  Main mover for application.
///
class AntibodyDesignProtocol : public protocols::moves::Mover {
	
public:

	AntibodyDesignProtocol();

	virtual ~AntibodyDesignProtocol();
	
	virtual std::string 
	get_name() const;
	
	//virtual protocols::moves::MoverOP clone() const;

	///@brief Parse my tag for RosettaScripts.  Main RS interface is in AntibodyDesignMover.
	/// This is just a small implementation, controlled mainly through cmd-line flags.
	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	);
	
	////////////////////////////////////////////////////////////////////////////
	// Optional Custom Settings
	//
	//
	
	///@brief Set the global scorefunction.  Used for snugdock/relax/total energy/etc.
	/// Not Used for Docking within the protocol.  
	void
	set_scorefxn(core::scoring::ScoreFunctionOP scorefxn);
	
	///@brief Set the min scorefunction used by the AntibodyDesignMover during the Minimization step.
	/// Should include dihedral_constraint weights if using default CDR constraints added during the protocol.
	void
	set_scorefxn_min(core::scoring::ScoreFunctionOP scorefxn);
	
	///@brief Set the instruction file path instead of reading it from the cmd-line options.
	void
	set_instruction_file_path(std::string instruction_file);
	
	
	////////////////////////////////////////////////////////////////////////////
	// Algorithm Settings
	//
	//
	

	///@brief Run SnugDock after main design runs
	void
	set_run_snugdock(bool setting);
	
	///@brief Run Dualspace Relax after main design runs
	void
	set_run_relax(bool setting);
	
	protocols::moves::MoverOP
	fresh_instance() const;

	////////////////////////////////////////////////////////////////////////////
	// Main Method.  All options read from cmd-line or set through individual classes.
	//
	//
	
	virtual void 
	apply( core::pose::Pose & pose );
	
private:
	
	void
	read_cmd_line_options();
	
	void
	setup_design_mover();
	
	///@brief Set constraint and chainbreak score on scorefunction if not already set.
	void
	setup_scorefxns();
	
	
	////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////
	
	
	
	///@brief Post-design step modeling.  Less aggressive, more high resolution.  Default false.
	void
	model_post_design(core::pose::Pose & pose);
	
	
	void
	init_on_new_input( core::pose::Pose const & pose );
	
	///@brief Used to output ongoing current ensembles during the protocol.  Specify a range in the vector to output
	void
	output_ensemble( utility::vector1< core::pose::PoseOP > ensemble, core::Size range_start, core::Size range_end, std::string prefix);
	
	
        void
	reorder_poses(utility::vector1<core::pose::PoseOP> & poses);
	
private:
	
	AntibodyDesignMoverOP graft_designer_;
	
	constraints::CDRDihedralConstraintMoverOP cdr_dihedral_cst_mover_;
	
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_min_;
	
	AntibodyInfoOP ab_info_;

	//AntibodyDatabaseManagerOP cdr_db_parser_;

	bool run_graft_designer_;
	bool run_snugdock_;
	bool run_relax_;
	
	bool remove_antigen_;
	
	std::string instruction_file_;
	utility::vector1<CDRNameEnum> design_override_;
	
};
} //design
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_design_AntibodyDesignProtocol_hh
