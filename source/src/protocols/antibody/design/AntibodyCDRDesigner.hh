// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyCDRDesigner.hh
/// @brief  Mover class that designs cdrs using a number of different strategies.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodyCDRDesigner_hh
#define INCLUDED_protocols_antibody_design_AntibodyCDRDesigner_hh

#include <protocols/antibody/design/AntibodyCDRDesigner.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>

#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace antibody {
namespace design {
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using namespace core::scoring;
	
	struct CDRDesignInstructions{
		bool design;
		bool conservative_design;
	};

///@brief Designs antibody CDR structures using a variety of different methods.  Main methods involve using cluster-based sequence probability statistics and conservative design strategies 
/// to limit unknown structural effects caused by aggressive design.  Uses North_AHO numbering and cluster-based design.  Part of Rosetta Antibody Designer.
///
class AntibodyCDRDesigner : public protocols::moves::Mover {

public:
	
	//Default constructor for RS?
	//AntibodyCDRDesigner();
	
	AntibodyCDRDesigner(AntibodyInfoOP ab_info);
	
	AntibodyCDRDesigner(AntibodyInfoOP ab_info, std::string instruction_path);
	
	virtual ~AntibodyCDRDesigner();
	
	virtual
	std::string
	get_name() const;
	
	
public:
	
	////////////////////////////////////////////////////////////////////////////
	// General Settings
	//
	//
	
	///@brief Set to design CDRs.  Default is all of them true.
	void
	set_cdr(CDRNameEnum const cdr, const bool setting);
	
	///@brief Set to model only one cdr, or all others but one.
	void
	set_cdr_only(CDRNameEnum const cdr, bool const setting);
	
	///@brief Set a range of CDRs.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool const setting);
	
	///@brief Repack neighbors of CDR's being designed within this distance.
	void
	set_neighbor_detection_dis(core::Real const neighbor_distance);
	
	void
	set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn);

	////////////////////////////////////////////////////////////////////////////
	// Accessors
	//
	//
	
	//core::pack::task::TaskFactoryOP
	//get_task_factory();
	
	////////////////////////////////////////////////////////////////////////////
	// Tweaking Settings
	//
	//
	
	
	void
	set_rounds(core::Size const rounds);
	
	///@brief Set a conservative method for design, including only using general conserved aminoacids, and specific amino acids for identified turn types.
	/// Default true for H3.  Any probabilities not found for a particular cluster will result in conservative design.
	///
	void
	set_use_conservative_design(CDRNameEnum const cdr, bool const setting);
	
	void
	set_use_conservative_design_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool const setting);
	
	
	///@brief Try to conserve turn structure using known turn-based conservative mutations during conservative design.
	void
	set_use_turn_conservation(bool const setting);
	
	
	///@brief Use conservative mutations (or alternative method) instead of using cluster sequence probabilities for design
	/// if the number of sequences in the particular CDR's cluster probability data is lower than this cutoff. Default is 10.  This is why we try and stay in type 1 clusters during graft.
	void
	set_probability_data_cutoff(core::Size const cutoff);
	
	///@brief Use these weights during probabilistic design for data that is normally zero.
	void
	set_zero_prob_weight_at(core::Real const weight);
	
	///@brief keep proline fixed.  Maybe should be not an option - and just do it.
	void
	no_design_proline(bool const setting);
	
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
	
public:
	
	virtual
	void
	apply(core::pose::Pose & pose);

private:
	
	///@brief Reads from database, gets prob_set for ResidueProbDesignOperation + sets cdrs_wo_prob_data_.
	std::map< core::Size, std::map< core::chemical::AA, core::Real > >
	setup_probability_data(core::pose::Pose & pose);
	
	///@brief Setup the Designer's taskfactory.
	core::pack::task::TaskFactoryOP
	setup_task_factory(core::pose::Pose & pose);
	
	///@brief Get a list of residues where conservative design will be used. 
	vector1<core::Size>
	get_conservative_design_residues(core::pose::Pose & pose);
	
	///@brief Turns off CDRs for design that are set to off.
	void
	disable_design_cdrs(core::pack::task::TaskFactoryOP tf, core::pose::Pose & pose);
	
	void
	disable_design_cdr(CDRNameEnum cdr, core::pack::task::TaskFactoryOP tf, core::pose::Pose & pose);
	
	//void
	//disable_packing_cdr(CDRNameEnum cdr, core::pack::task::TaskFactoryOP tf);
	
private:
	
	///@brief Set cluster-based harmonic dihedral constraints if known cluster or set coordinate constraints.
	void
	setup_constraints(core::pose::Pose & pose);
	
	///@Does The CDR have constraints of this type.  A bit hacky.  Used for going form graftdesigner to this, without regenerating the 'wrong' constraints.
	///@details probably only works 'correctly' for dihedral/coordinate constraints.  THe function assumes that constraints were set for the whole CDR.
	bool
	cdr_has_constraints(core::pose::Pose const & pose, CDRNameEnum const cdr, std::string const constraint_type);
	
	void
	read_command_line_options();
	
	///@brief Call DesignInstructionParser to determine settings for the design run.
	void
	read_instructions(std::string instruction_path);
	
	///@brief Removes  residues from prob_set from instruction settings.  Used so that we speed task generation instead of overwriting these residues.
	void
	remove_conservative_design_residues_from_prob_set(vector1<core::Size> const & positions, std::map< core::Size, std::map< core::chemical::AA, core::Real > > & prob_set);
	
	AntibodyInfoOP ab_info_;
	DesignTypeEnum design_method_;
	ScoreFunctionOP scorefxn_;
	core::Real zero_prob_weight_;
	core::Real neighbor_dis_;
	
	core::Size prob_cutoff_;
	core::Size rounds_;
	
	bool turn_conservation_; //Setting to use turn-based conservation during conservative design
	bool no_design_proline_;
	bool basic_design_; //Do not use fancy taskops and data.  Simply design the CDRs using chosen method.
	bool use_cluster_constraints_;
	
	std::map<CDRNameEnum, CDRDesignInstructions> instructions_;
	std::string instruction_path_;
	
	//vector1 <bool> cdrs_;
	
	//vector1<CDRNameEnum> cdrs_wo_prob_data_; //CDRs for which no probability data could be assigned.
};

	
} //design
} //antibody
} //protocols

#endif	//INCLUDED_protocols_antibody_design_AntibodyCDRDesigner.hh

