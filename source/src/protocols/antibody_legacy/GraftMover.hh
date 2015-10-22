// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/antibody_legacy/GraftMover.hh
/// @brief
/// @author Aroop Sircar (aroopsircar@yahoo.com)


#ifndef INCLUDED_protocols_antibody_legacy_GraftMover_hh
#define INCLUDED_protocols_antibody_legacy_GraftMover_hh


// Rosetta headers
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <utility/vector1.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace antibody_legacy {

class GraftMover;
typedef utility::pointer::shared_ptr< GraftMover > GraftMoverOP;
typedef utility::pointer::shared_ptr< const GraftMover > GraftMoverCOP;

//////////////////////////////////////////////////////////////////////////
/// @brief Grafts a series of CDR onto a framework
/// @details
class GraftMover : public moves::Mover {
public:
	// default constructor
	GraftMover();

	// default destructor
	~GraftMover();

	inline void enable_graft_l1( bool setting ) {
		graft_l1_ = setting;
	}
	inline void enable_graft_l2( bool setting ) {
		graft_l2_ = setting;
	}
	inline void enable_graft_l3( bool setting ) {
		graft_l3_ = setting;
	}
	inline void enable_graft_h1( bool setting ) {
		graft_h1_ = setting;
	}
	inline void enable_graft_h2( bool setting ) {
		graft_h2_ = setting;
	}
	inline void enable_graft_h3( bool setting ) {
		graft_h3_ = setting;
	}
	inline void set_camelid( bool setting ) {
		camelid_ = setting;
	}

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

	/// @brief relax optimized CDR grafted regions
	void relax_optimized_CDR_grafts( core::pose::Pose & pose_in );

	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

private:
	// Grafting options
	bool graft_l1_;
	bool graft_l2_;
	bool graft_l3_;
	bool graft_h1_;
	bool graft_h2_;
	bool graft_h3_;
	/// @brief benchmark flag
	bool benchmark_;
	bool camelid_;

	// Packer
	protocols::simple_moves::PackRotamersMoverOP packer_;

	void set_packer_default(
		core::pose::Pose & pose_in,
		core::scoring::ScoreFunctionOP scorefxn,
		bool include_current
	);
}; // class GraftMover

class GraftOneMover;
typedef utility::pointer::shared_ptr< GraftOneMover > GraftOneMoverOP;
typedef utility::pointer::shared_ptr<const GraftOneMover> GraftOneMoverCOP;

//////////////////////////////////////////////////////////////////////////
/// @brief Grafts only one CDR onto a framework
/// @details
class GraftOneMover : public moves::Mover {
public:
	// default constructor
	//GraftOneMover();
	// constructor with arguments
	GraftOneMover(
		core::Size query_start,
		core::Size query_end,
		std::string template_name
	);

	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

private:
	// Limits of query loop
	core::Size query_start_;
	core::Size query_end_;
	// Limits of template loop
	core::Size template_start_;
	core::Size template_end_;
	std::string template_name_;
	core::pose::Pose template_pose_;
	/// @brief benchmark flag
	bool benchmark_;

}; // class GraftOneMover

class CloseOneMover;
typedef utility::pointer::shared_ptr< CloseOneMover > CloseOneMoverOP;
typedef utility::pointer::shared_ptr<const CloseOneMover> CloseOneMoverCOP;

/// @brief Closes only one CDR onto a framework
class CloseOneMover : public moves::Mover {
public:
	// default constructor
	// CloseOneMover();
	// constructor with arguments
	CloseOneMover(
		core::Size query_start,
		core::Size query_end  );

	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	void close_one_loop_stem (
		core::pose::Pose & pose_in,
		core::Size cutpoint_in,
		bool nter );

	void close_one_loop_stem (
		core::pose::Pose & pose_in,
		core::Size loop_begin,
		core::Size loop_end,
		core::Size cutpoint );

	void close_one_loop_stem_helper(
		core::Size loop_begin,
		core::Size loop_end,
		core::Size cutpoint,
		core::pose::Pose & pose_in,
		core::kinematics::MoveMapOP loop_map
	);

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

private:
	// Limits of query loop
	core::Size loop_start_;
	core::Size loop_end_;
	core::Real allowed_separation_;
	core::Size flanking_residues_;
	/// @brief benchmark flag
	bool benchmark_;

}; // class CloseOneMover

class LoopRlxMover;
typedef utility::pointer::shared_ptr< LoopRlxMover > LoopRlxMoverOP;
typedef utility::pointer::shared_ptr<const LoopRlxMover> LoopRlxMoverCOP;

/// @brief Closes only one CDR onto a framework
class LoopRlxMover : public moves::Mover {
public:
	// default constructor
	// LoopRlxMover();
	// constructor with arguments
	LoopRlxMover(
		core::Size query_start,
		core::Size query_end  );

	void set_default();

	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	void setup_packer_task( core::pose::Pose & pose_in );

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

private:
	// Limits of query loop
	core::Size loop_start_;
	core::Size loop_end_;
	/// @brief benchmark flag
	bool benchmark_;

	// score functions
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

}; // class LoopRlxMover

} // antibody
} // protocols


#endif
