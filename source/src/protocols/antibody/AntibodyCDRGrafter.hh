// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyCDRGrafter.hh
/// @brief Class to graft CDR loops from an antibody to a new antibody or from a CDR pose into a different antibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyCDRGrafter_hh
#define INCLUDED_protocols_antibody_AntibodyCDRGrafter_hh

#include <protocols/antibody/AntibodyCDRGrafter.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace antibody {



///@brief Class to graft CDR loops from an antibody to a new antibody or from a CDR pose into a different antibody.
///  Independant of antibody and CDR modeling.
///  Results in 100 percent loop closure if using both graft algorithms.
///
///@details Recommended to use cluster-based or general dihedral constraints on the CDR with a min or relax to optimize graft.
/// Use optimiize_cdrs option to optimize the CDRs that were grafted and any neighbor CDRs using these dihedral constraints during relax.
///
/// By default, uses two residues on N and C terminal of insert CDR and scaffold to close loop.
/// Use getters of primary and secondary graft mover to change this if needed.
/// By default, stops after the graft is closed!
///
class AntibodyCDRGrafter : public protocols::moves::Mover {

public:

	AntibodyCDRGrafter();

	///@brief Default constructor
	AntibodyCDRGrafter( AntibodyInfoOP ab_info);

	///@brief Constructor with most needed arguments.
	///
	///@details Donor structure is what we will be grafting from.  Either a CDR or a whole structure or chain.
	///  cdrs_to_graft: Boolean vector of CDRNameEnums (1-8) (includes CDR4/DE loop)
	///  overhangs: Set the N and Cterminal overhang used for superposition before the graft.
	///
	AntibodyCDRGrafter( AntibodyInfoOP ab_info,
		core::pose::Pose const & donor_structure,
		utility::vector1<bool> const & cdrs_to_graft,
		core::Size nter_overhang = 3,
		core::Size cter_overhang = 3 );



	AntibodyCDRGrafter( AntibodyCDRGrafter const & src );

	virtual ~AntibodyCDRGrafter();

	virtual void
	apply( core::pose::Pose & pose );



public:

	///@brief Set the CDRs to graft.  Boolean vector of CDRNameEnums (1-8) (includes CDR4/DE loop)
	void
	set_cdrs(utility::vector1<bool> const & cdrs);

	///@brief Set a single CDR we are going to graft in.
	void
	set_cdr_only( CDRNameEnum cdr );

	///@brief Set a donor structure.  Culd be a single CDR, multiple CDRs or a whole other antibody.
	void
	set_donor_structure( core::pose::Pose const & pose );


	///@brief Set the N and Cterminal overhang used for superposition before the graft.
	void
	set_overhang(core::Size nter_overhang, core::Size cter_overhang);

	///@brief Set a boolean for whether to stop after closure or not.  Default TRUE.
	void
	set_stop_after_closure( bool stop_after_closure );


	/// @brief Set a boolean for whether to use the secondary graft mover is the graft is not closed using the first.
	///  Default FALSE.
	///
	/// @details
	///  Using this results in 100% graft closure, however, in order to fix the structure, a dihedral constrained relax is
	///  done on the CDR while repacking neighbor residues.
	///  85% Graft closure without (slight chainbreaks  - not as bad as what used to be used for antibody modeling.)
	///
	void
	set_use_secondary_graft_mover_if_needed( bool use_secondary_graft_mover );


	///@brief Set to optimize any grafted and neighbor CDRs using dihdedral constrained relax on them. Default FALSE.
	///  Recommended if using secondary graft mover.
	void
	set_optimize_cdrs( bool optimize_cdrs );

	///@brief Set to include the DE loop, or CDR4 if it is a neighbor to a grafted CDR and optimization is set to on.
	/// Recommended if optimizing.  Default TRUE.
	void
	set_include_cdr4_in_optimization( bool include_cdr4 );


	///@brief Set the dihedral constraint weight used if optimization is done and the scorefxn dihedral_constraint
	/// weight is zero.  Default = 2.0
	void
	set_dihedral_constraint_weight ( core::Real dih_cst_wt );

	///@brief Idealize the insert.  Default TRUE.
	void
	set_idealize_insert(bool idealize_insert);


	///@brief Set the (all-atom) scorefunction used by the graft movers for packing
	/// Otherwise use cmd-line default
	void
	set_scorefxn_pack( core::scoring::ScoreFunctionCOP scorefxn );

	///@brief Set the (low-res) scorefunction used by the graft movers to close loop.
	///  Default is to use Frank Dimaio's smoothed low-res terms.
	void
	set_scorefxn_low( core::scoring::ScoreFunctionCOP scorefxn );

public:

	///@brief Get a reference of the main graft mover to tweak settings.
	protocols::grafting::CCDEndsGraftMover &
	get_primary_graft_mover();

	///@brief Get a reference of the secondary graft mover to tweak settings.
	protocols::grafting::AnchoredGraftMover &
	get_secondary_graft_mover();


public:

	///@brief Apply graft mover to pose, inserting the cdr_region.
	///@details Return success or failure of the graft and the position of the failure.
	/// Public in case one wants to use it by itself...
	std::pair< bool, core::Size >
	apply_to_cdr(
		core::pose::Pose & pose,
		core::pose::Pose & cdr_region_with_overhang,
		CDRNameEnum const cdr,
		protocols::grafting::AnchoredGraftMoverOP grafter );


public:

	std::string
	get_name() const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	protocols::moves::MoverOP
	clone() const;

	//AntibodyCDRGrafter & operator=( AntibodyCDRGrafter const & src );

	virtual moves::MoverOP fresh_instance() const;


private:

	void
	setup_classes();

	void
	set_defaults();





private:

	AntibodyInfoOP ab_info_;
	core::pose::PoseOP donor_structure_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_low_;

	utility::vector1<bool> cdrs_to_graft_;

	protocols::grafting::CCDEndsGraftMoverOP graft_mover_;
	protocols::grafting::AnchoredGraftMoverOP anchored_graft_mover_;



	bool use_secondary_graft_mover_;
	bool stop_after_closure_;

	core::Size nter_overhang_;
	core::Size cter_overhang_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;
	CDRDefinitionEnum cdr_definition_;

	bool optimize_cdrs_;
	bool include_cdr4_;
	utility::vector1< utility::vector1< bool > > neighbor_cdrs_;  //Neighbor cdrs not including self for optimization.

	core::Real dihedral_cst_weight_;
};


}//protocols
}//antibody


#endif //protocols/antibody_AntibodyCDRGrafter_hh







