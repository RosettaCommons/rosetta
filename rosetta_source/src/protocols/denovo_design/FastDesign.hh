// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/FastDesign.hh
/// @brief The FastDesign Protocol
/// @detailed
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_FastDesign_hh
#define INCLUDED_protocols_denovo_design_FastDesign_hh

// Unit headers
#include <protocols/denovo_design/FastDesign.fwd.hh>
#include <protocols/relax/RelaxProtocolBase.hh>

// Protocol headers
//#include <protocols/denovo_design/RestrictWorstRegion.fwd.hh>
//#include <protocols/flxbb/filters/HighestEnergyRegion.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/filters/Filter.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>



namespace protocols {
namespace denovo_design {

struct FilterParams {
	FilterParams() {
		filter = 0;
		use_tolerance = false;
		use_threshold = false;
		threshold = 0.0;
		tolerance = 0.0;
		sample_low = true;
	}
	filters::FilterOP filter;
	bool use_tolerance;
	bool use_threshold;
	core::Real threshold;
	core::Real tolerance;
	bool sample_low;
};

struct RelaxScriptCommand {
	RelaxScriptCommand(){
		command="";
		param1=0;
		param2=0;
		param3=0;
		param4=0;
		nparams=0;
	}
	std::string command;
	core::Real  param1;
	core::Real  param2;
	core::Real  param3;
	core::Real  param4;
	core::Size  nparams;
};

class FastDesign : public relax::RelaxProtocolBase {
public:

	/// @brief default constructor
	FastDesign();

	/// @brief copy constructor
	FastDesign( FastDesign const & rval );

  /// @brief virtual constructor to allow derivation
	virtual ~FastDesign();

  /// @brief Parses the FastDesignTags
	void parse_my_tag(
	  utility::tag::TagPtr const tag,
	  protocols::moves::DataMap & data,
	  protocols::filters::Filters_map const &,
	  protocols::moves::Movers_map const &,
	  core::pose::Pose const &
	);

  virtual void parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & score_fxns,
		utility::lua::LuaObject const & tasks,
		protocols::moves::MoverCacheSP cache );

  /// @brief Return the name of this mover.
  virtual std::string get_name() const;

  /// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

  /// @brief Apply the FastDesign. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	/// @brief sets the movemap to not allow DNA to move during relax.
	void makeDnaRigid( core::pose::Pose & pose, core::kinematics::MoveMapOP mm );

	/// @brief sets the number of repeats
	void set_repeats( core::Size const repeats );

	/// @brief initializes cache of allowed amino acids
	void initialize_aa_cache( core::pose::Pose const & pose );

protected:

  void cmd_accept_to_best(
    const core::scoring::ScoreFunctionOP local_scorefxn,
    core::pose::Pose &pose,
    core::pose::Pose &best_pose,
    const core::pose::Pose &start_pose,
    core::Real       &best_score,
    core::Size       &accept_count
  );

  void do_minimize(
    core::pose::Pose &pose,
    core::Real tolerance,
    core::kinematics::MoveMapOP local_movemap,
    core::scoring::ScoreFunctionOP local_scorefxn
  );

	/// @brief outputs a pdb file with sequential numbering, with the specified prefix
	void dump_pdb( core::pose::Pose const & pose, std::string const & prefix = "dump" );

private:

	void read_script_file( const std::string &script_file, core::Size standard_repeats = 5  );

	/// @brief initializes and creates the packer task
	core::pack::task::PackerTaskOP
	create_packer_task( core::pose::Pose const & pose,
											core::kinematics::MoveMapCOP local_movemap,
											utility::vector1< core::Size > const & allow_repacking ) const;

	/// @brief checks each mutation at each designable position, and disallows those mutations that hurt the filter score
	utility::vector1< utility::vector1< bool > >
	check_and_disallow_mutations_by_filter( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, filters::FilterOP filter ) const;


	/// @brief increments the number of times a particular residue (and aa) have been designed on
	void increment_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos );

	/// @brief checks how many times the residue at seqpos has been redesigned. If the value is larger than than the maximum allowed number of redesigns, return false, otherwise return true
	bool check_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos ) const;

	/// @brief Creates and returns a pointer to a taskoperation that designs the worst regions of a protein only. Returns NULL if we don't need the taskop
	toolbox::task_operations::DesignAroundOperationOP
	create_worst_region_operation( core::pose::Pose const & pose );

private:   // options

	/// @brief Number of repeats if using default script
  core::Size default_repeats_;

  /// @brief Apply Ramady(tm) fix to residues with bad rama scores ?
  bool ramady_;

  /// @brief Number of bad ramas to fix per call to 'fix_worst_bad_ramas'
	core::Size ramady_num_rebuild_;

  /// @brief Reject ramady changes perturbing structure more than this amount
	core::Real ramady_rms_limit_;

  /// @brief Cutoff for calling a rama 'bad'
	core::Real ramady_cutoff_;

  /// @brief Force ramady to be run (normally skip rate of 10%)
	bool ramady_force_;

  /// @brief Allow Chi angles to move ?
	bool repack_;

	/// @brief Do only a few test_cycles ?
	bool test_cycles_;

	/// @brief Dump pdb after repack, min, or ramp_repack_min?
	bool dumpall_;

  /// @brief Quit after this many accepts ?// limits themaximum number of accepts, default is 1000000
	core::Size script_max_accept_;

  /// @brief Do a symmetric RMSD calculation rather then a monomeric one.
  bool symmetric_rmsd_;

	/// @brief Allow design as well as repacking? default=false
	bool allow_design_;

		/// @brief Should we design only the worst-scoring 8A shell around a residue? default=false
	bool design_worst_;

	/// @brief Should we design only the shell around the worst psipred residue? default=false
	bool design_by_psipred_;

	/// @brief should we design only the shell around the worst 9-mer fragment? default=false
	bool design_by_frag_qual_;

	/// @brief should we only design residues with type that has changed since the last time this mover was called?
	bool only_design_changes_;

	/// @brief Blueprint file used for psipred predictions
	std::string blueprint_file_;

	/// @brief Should we ramp down design constraints as the run is proceeding? default=false
	bool ramp_design_constraints_;

	/// @brief Should we clear designable residues (sets them to ALA so that they can be rebuilt)
	bool clear_designable_residues_;

	/// @brief scorefunction used to score and rank poses (NOT for packing)
	core::scoring::ScoreFunctionOP rank_scorefxn_;

	/// @brief should we prevent the native amino acid, if it's one of the worst in the pose?
	bool disallow_bad_natives_;

	/// @brief the maximum number of times we should redesign a particular position and amino acid before disallowing that amino acid
	core::Size max_redesigns_;

	/// @brief resfile
	std::string resfile_;

private:   // other data

	protocols::checkpoint::CheckPointer checkpoints_;

	std::vector <RelaxScriptCommand> script_;
	utility::tag::TagPtr movemap_tag_; // this cannot be parsed before apply b/c the fold tree is likely to change during a run
	core::Size dump_counter_; // counter for number of pdb dumps
	utility::vector1< FilterParams > filters_; // filters to be applied at each step
	utility::vector1< utility::vector1 < bool > > allowed_aas_; // list of allowed amino acids at each position
	utility::vector1< utility::vector1< core::Size > > num_redesigns_; // list of the number of times each position has been redesigned
	//flxbb::filters::DesignByFragmentQualityOperationOP frag_qual_op_;
	core::Size regions_to_design_;
	core::Size run_count_;
	std::string cached_sequence_;
	//RestrictWorstRegionOP worst_region_mover_;
};




}
} // protocols

#endif
