// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/AntibodySeqDesignMover.hh
/// @brief  Mover class that designs cdrs using a number of different strategies.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodySeqDesignMover_hh
#define INCLUDED_protocols_antibody_design_AntibodySeqDesignMover_hh

#include <protocols/antibody/design/AntibodySeqDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.fwd.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/util.hh>

#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.fwd.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/Mover.hh>

namespace protocols {
namespace antibody {
namespace design {

///@brief Designs antibody CDR structures using a variety of different methods.  Main methods involve using cluster-based sequence probability statistics and conservative design strategies 
/// to limit unknown structural effects caused by aggressive design.  Uses North_AHO numbering and cluster-based design.  Part of Rosetta Antibody Designer.
/// ->  Will be replaced by SnugDesign Mover in the application, but still useful in it's own right.
///
class AntibodySeqDesignMover : public protocols::moves::Mover {

public:
	
	//Default constructor for RS?
	//AntibodySeqDesignMover();
	
	AntibodySeqDesignMover(AntibodyInfoOP ab_info);
	
	AntibodySeqDesignMover(AntibodyInfoOP ab_info, std::string instruction_path);
	
	virtual ~AntibodySeqDesignMover();
	
	///@brief Sets defaults, Reads defaultInstruction file, read user-passed instruction file.
	void
	set_defaults();
	
	virtual
	void
	apply(core::pose::Pose & pose);
	
	virtual
	std::string
	get_name() const;
	
	
public:
	
	////////////////////////////////////////////////////////////////////////////
	// General Settings
	//
	//
	
	///@brief Set the CDR-specific settings for SeqDesign
	void
	set_seq_designer_options(AntibodyCDRSeqDesignOptions options);
	
	
	///@brief Set any paratope CDRs.  If not set, will use all CDRs as the paratope where needed.
	void
	set_paratope_cdrs(utility::vector1<bool> const & cdrs);
	
	///@brief Set any epitope residues in PDB numbering.  If not set, they will be detected automatically.
	void
	set_epitope_residues(utility::vector1<PDBNumbering > epitope_residues);
	
	
	///@brief Repack neighbors of CDR's being designed within this distance.
	void
	neighbor_detection_dis(core::Real const neighbor_distance);
	
	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);
	
	
	void
	set_rounds(core::Size const rounds);
	
	
	////////////////////////////////////////////////////////////////////////////
	// Design Type Settings
	// -These could each be a (configurable) step instead-
	//
	
	///@brief Type of design method used.  Acceptable options - relaxed_design, fixbb, flxbb
	void
	set_design_method(DesignTypeEnum const design_method);
	
	///@brief Sets the designer to not use our fancy taskops and statistics.  Only use basic design with all residues turned on. 
	/// Use for benchmarking.
	void
	set_basic_design(bool const setting);
	
	///@brief Should we use the cdr cluster CircularHarmonic constraints?  Default true - this is mainly used for benchmarking.
	void
	set_use_cluster_constraints(bool const setting);
	
	
private:
	
	///@brief Set cluster-based harmonic dihedral constraints if known cluster or set coordinate constraints.
	void
	setup_cdr_constraints(core::pose::Pose & pose);
	
	void
	read_command_line_options();
	
	///@brief Call DesignInstructionParser to determine settings for the design run.
	void
	setup_options_classes();
	
	void
	setup_paratope_epitope_constraints(core::pose::Pose & pose);
	
	void
	setup_scorefxn();
	
private:
	AntibodyCDRSeqDesignOptions cdr_seq_design_options_;
	
	AntibodyInfoOP ab_info_;
	DesignTypeEnum design_method_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Real neighbor_dis_;
	
	core::Size rounds_;
	
	bool basic_design_; //Do not use fancy taskops and data.  Simply design the CDRs using chosen method.
	bool full_cdr_min_;
	bool use_cluster_constraints_;
	bool use_epitope_constraints_;
	
	//std::map<CDRNameEnum, CDRDesignInstructions> instructions_;
	std::string instruction_path_;
	std::string default_instruction_path_;
	
	AntibodyDesignModelerOP modeler_;
	protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverOP paratope_epitope_cst_mover_;
	protocols::antibody::constraints::ParatopeSiteConstraintMoverOP paratope_cst_mover_;
	
	vector1<bool> paratope_cdrs_;
	utility::vector1<PDBNumbering > epitope_residues_; //Vector of resnum, chain pairs.
	
	
	//vector1 <bool> cdrs_;
	
	//vector1<CDRNameEnum> cdrs_wo_prob_data_; //CDRs for which no probability data could be assigned.
};

	
} //design
} //antibody
} //protocols

#endif	//INCLUDED_protocols_antibody_design_AntibodySeqDesignMover.hh

