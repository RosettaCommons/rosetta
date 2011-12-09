// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file antibody2/moves/GraftMover2.hh
/// @brief
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_moves_GraftMover2_hh
#define INCLUDED_protocols_antibody2_moves_GraftMover2_hh

#include <protocols/antibody2/GraftMover2.fwd.hh>

// Rosetta headers
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/antibody2/AntibodyInfo.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace antibody2 {

	//////////////////////////////////////////////////////////////////////////
	/// @brief Grafts a series of CDR onto a framework
	/// @details
	class GraftMover2 : public protocols::moves::Mover {
	public:
		typedef std::map < std::string, bool > GraftMap;
		// default constructor
		GraftMover2();

		/// @brief constructor with arguments
		GraftMover2( bool l1, bool l2, bool l3, bool h1, bool h2, bool h3, bool camelid, bool benchmark );

		// default destructor
		~GraftMover2();

		void init( bool l1, bool l2, bool l3, bool h1, bool h2, bool h3, bool camelid, bool benchmark );

		inline void enable_graft_l1( bool setting ) { graft_l1_ = setting; }
		inline void enable_graft_l2( bool setting ) { graft_l2_ = setting; }
		inline void enable_graft_l3( bool setting ) { graft_l3_ = setting; }
		inline void enable_graft_h1( bool setting ) { graft_h1_ = setting; }
		inline void enable_graft_h2( bool setting ) { graft_h2_ = setting; }
		inline void enable_graft_h3( bool setting ) { graft_h3_ = setting; }
		inline void set_camelid( bool setting ) { camelid_ = setting; }

		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

		/// @brief relax optimized CDR grafted regions
		void relax_optimized_CDR_grafts( core::pose::Pose & pose );

		void set_default();
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;

	private:
		// Grafting options
		bool graft_l1_;
		bool graft_l2_;
		bool graft_l3_;
		bool graft_h1_;
		bool graft_h2_;
		bool graft_h3_;

		GraftMap grafts_;

		/// @brief benchmark flag
		bool benchmark_;
		bool camelid_;

		bool user_defined_;
		bool first_apply_with_current_setup_;

		// movers
		protocols::moves::SequenceMoverOP graft_sequence_, relax_sequence_;
		protocols::simple_moves::PackRotamersMoverOP packer_;
		protocols::moves::PyMolMoverOP pymol_;

		core::scoring::ScoreFunctionOP scorefxn_;

		void finalize_setup( AntibodyInfo & ab_info );
		void set_packer_default(
			core::pose::Pose & pose,
			bool include_current );
	}; // class GraftMover2

	class GraftOneMover;
	typedef utility::pointer::owning_ptr< GraftOneMover > GraftOneMoverOP;
	typedef utility::pointer::owning_ptr<const GraftOneMover> GraftOneMoverCOP;



















	//////////////////////////////////////////////////////////////////////////
	/// @brief Grafts only one CDR onto a framework
	/// @details
	class GraftOneMover : public protocols::moves::Mover {
	public:
		// default constructor
		//GraftOneMover();
		// constructor with arguments
		GraftOneMover(
			core::Size query_start,
			core::Size query_end,
			std::string template_name,
			core::scoring::ScoreFunctionOP scorefxn );

		void set_default( std::string template_name );
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

		std::string template_name_;

		// Limits of template loop
		core::Size template_start_;
		core::Size template_end_;
		core::pose::Pose template_pose_;

		/// @brief benchmark flag
		bool benchmark_;

		core::scoring::ScoreFunctionOP scorefxn_;

	}; // class GraftOneMover















	class CloseOneMover;
	typedef utility::pointer::owning_ptr< CloseOneMover > CloseOneMoverOP;
	typedef utility::pointer::owning_ptr<const CloseOneMover> CloseOneMoverCOP;

	/// @brief Closes only one CDR onto a framework
	class CloseOneMover : public protocols::moves::Mover {
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

		void set_pymol( protocols::moves::PyMolMoverOP pymol ) { pymol_ = pymol; }

		/// @brief enable benchmark mode
		inline void enable_benchmark_mode( bool setting ) {
			benchmark_ = setting;
		}

	private:
		// Limits of query loop
		core::Size loop_start_, cdr_loop_start_;
		core::Size loop_end_, cdr_loop_end_;
		core::Real allowed_separation_;
		core::Size flanking_residues_;

		/// @brief benchmark flag
		bool benchmark_;

		core::kinematics::MoveMapOP movemap_;
		protocols::moves::PyMolMoverOP pymol_;

	}; // class CloseOneMover














class LoopRlxMover;
typedef utility::pointer::owning_ptr< LoopRlxMover > LoopRlxMoverOP;
typedef utility::pointer::owning_ptr<const LoopRlxMover> LoopRlxMoverCOP;














/// @brief Closes only one CDR onto a framework
class LoopRlxMover : public protocols::moves::Mover {
public:
	// default constructor
	// LoopRlxMover();
	// constructor with arguments
	LoopRlxMover(
		core::Size query_start,
		core::Size query_end  );

	void set_default();
	void setup_objects( core::pose::Pose & pose );

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
	core::Size inner_cycles_, outer_cycles_;
	core::Real temperature_, gamma_;
	/// @brief benchmark flag
	bool benchmark_;

	protocols::moves::MonteCarloOP mc_;
	protocols::simple_moves::MinMoverOP loop_min_mover_;
	protocols::moves::SequenceMoverOP wiggle_loop_;

	protocols::loops::LoopOP one_loop_;

	// score functions
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	core::kinematics::MoveMapOP movemap_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;

}; // class LoopRlxMover

} // antibody2
} // protocols


#endif
