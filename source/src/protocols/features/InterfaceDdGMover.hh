// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/InterfaceDdGMover.hh
/// @brief Calculates ddG of binding
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_InterfaceDdGMover_hh
#define INCLUDED_protocols_features_InterfaceDdGMover_hh

// Unit headers
#include <protocols/features/InterfaceDdGMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/SavePoseMover.fwd.hh>
#include <protocols/features/ReportToDB.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace features {

///@brief Calculates ddG of binding
class InterfaceDdGMover : public protocols::moves::Mover {

public:
	InterfaceDdGMover();

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~InterfaceDdGMover();

	void
	apply( core::pose::Pose & pose );

	void
	const_pose_apply( core::pose::Pose const & pose );

public:
	virtual std::string
	get_name() const;

	static
	std::string
	mover_name();

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	/// @brief "unbinds" a pose based on degrees of freedom previously defined
	void unbind(core::pose::Pose & pose) const;

	// Getters/setters
	const utility::vector1<core::Size> & get_chain_ids() const;

	///@brief appends chain_id
	void add_chain_id( core::Size chain_id, core::pose::Pose const & pose );

	///@brief converts chain_name to a chain_id, and then appends
	void add_chain_name( std::string chain_name, core::pose::Pose const & pose );

	///@brief gets chain_id from jump id, and then appends
	void add_jump_id( core::Size jump, core::pose::Pose const & pose );

	/// @brief Sets the scorefunction.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Stores the wildtype pose
	void set_wt_pose( core::pose::Pose const & pose );
	/// @brief Stores the mutant pose
	void set_mut_pose( core::pose::Pose const & pose );
	/// @brief Stores the wildtype pose
	void set_wt_pose( core::pose::PoseOP pose );
	/// @brief Stores the mutant pose
	void set_mut_pose( core::pose::PoseOP pose );

	/// @brief Clones db_reporter into all member variable database reporters, and sets batch names accordingly
	void set_db_reporter( protocols::features::ReportToDBOP db_reporter );

	core::pose::PoseCOP get_unbound_wt_pose() const;
	core::pose::PoseCOP get_bound_wt_pose() const;
	core::pose::PoseCOP get_unbound_mut_pose() const;
	core::pose::PoseCOP get_bound_mut_pose() const;

	utility::vector1< std::pair<std::string, core::Real> > get_all_scores() const;

	///@brief Determines which residues, if any, are different between bound_wt_pose_ and bound_mut_pose_
	utility::vector1< std::tuple<std::string, core::Size, std::string> > mutation_list() const;

private:

	///@brief Scores the cached poses, caches the scores
	core::Real compute();

	utility::vector1<core::Size> chain_ids_;

	core::scoring::ScoreFunctionOP scorefxn_;

	///@brief distance in angstroms to separate moledules
	core::Real translate_by_;

	///@brief unbound wildtype pose state; stored for later use
	core::pose::PoseOP unbound_wt_pose_;
	///@brief bound wildtype pose state; stored for later use
	core::pose::PoseOP bound_wt_pose_;

	///@brief unbound mutant pose state; stored for later use
	core::pose::PoseOP unbound_mut_pose_;
	///@brief bound mutant pose state; stored for later use
	core::pose::PoseOP bound_mut_pose_;

	// Cached scores
	core::Real unbound_wt_score_;
	core::Real bound_wt_score_;
	core::Real unbound_mut_score_;
	core::Real bound_mut_score_;
	core::Real ddG_score_;

	///@brief Keeps track of whether the pose passed to the apply function will be considered to be the wildtype of mutant structure
	bool apply_pose_is_mutant_;

	///@brief Used in Rosetta scripts context to save a reference pose of the wildtype pose state
	protocols::rosetta_scripts::SavePoseMoverOP wt_save_pose_mover_;
	///@brief Used in Rosetta scripts context to save a reference pose of the mutant pose state
	protocols::rosetta_scripts::SavePoseMoverOP mut_save_pose_mover_;

	// The following database reporters are used for features database reporting, and are cloned from "db_reporter" in parse_my_tag or set_db_reporter()
	// Their batch_names will reflect their variable names
	protocols::features::ReportToDBOP bound_wt_db_reporter_;
	protocols::features::ReportToDBOP unbound_wt_db_reporter_;
	protocols::features::ReportToDBOP bound_mut_db_reporter_;
	protocols::features::ReportToDBOP unbound_mut_db_reporter_;

	bool report_to_db_;

};

} //protocols
} //features

#endif //INCLUDED_protocols_features_InterfaceDdGMover_hh
