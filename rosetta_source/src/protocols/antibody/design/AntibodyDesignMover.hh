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


#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/design/AntibodyDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyGraftDesigner.hh>

#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace antibody {
namespace design {

using namespace protocols::antibody;	
using std::string;
using utility::vector1;


class AntibodyDesignMover : public protocols::moves::Mover {
	
public:

	AntibodyDesignMover();
	//AntibodyDesignMover( AntibodyDesignMover const & rhs );
	//AntibodyDesignMover & operator=( AntibodyDesignMover const & rhs );
        
	//Specific functions for Mover...
	virtual ~AntibodyDesignMover();
	virtual void apply( core::pose::Pose & pose );
	virtual string get_name() const;
	//virtual bool reinitialize_for_new_input() const;
	//virtual bool reinitialize_for_each_job() const;
	virtual protocols::moves::MoverOP clone() const;
	//virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief Reads options from options system

	//////////////////////// Optional Setters ////////////////////////////////////////////////////////////////
	void
	set_graft_designer(AntibodyGraftDesignerOP & graft_designer);
	
	void
	set_scorefxn(ScoreFunctionOP & scorefxn);
	
	void
	set_ab_info(AntibodyInfoOP & ab_info);
	
private:
	void read_options();
	void model_starting_pose(core::pose::Pose & pose);
	void run_graft_designer(core::pose::Pose & pose);
	void init_on_new_input( core::pose::Pose const & pose );
	///@brief CDR database parser for design and grafting.
	AntibodyDatabaseManagerOP cdr_db_parser_;
	AntibodyGraftDesignerOP graft_designer_;
	AntibodyInfoOP ab_info_;
	AntibodyDesignModelerOP modeler_;
	ScoreFunctionOP scorefxn_;
	
};
} //design
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
