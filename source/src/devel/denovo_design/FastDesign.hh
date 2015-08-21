// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/FastDesign.hh
/// @brief The FastDesign Protocol
/// @details
/// @author Tom Linsky


#ifndef INCLUDED_devel_denovo_design_FastDesign_hh
#define INCLUDED_devel_denovo_design_FastDesign_hh

// Unit headers
#include <devel/denovo_design/FastDesign.fwd.hh>
#include <protocols/relax/FastRelax.hh>

// Protocol headers
//#include <protocols/denovo_design/RestrictWorstRegion.fwd.hh>
//#include <protocols/flxbb/filters/HighestEnergyRegion.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/filters/Filter.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace denovo_design {

struct FilterParams {
	FilterParams() {
		filter.reset();
		use_tolerance = false;
		use_threshold = false;
		threshold = 0.0;
		tolerance = 0.0;
		sample_low = true;
	}
	protocols::filters::FilterOP filter;
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

class FastDesign : public protocols::relax::FastRelax {
public:

	/// @brief default constructor
	FastDesign();

	/// @brief virtual constructor to allow derivation
	virtual ~FastDesign();

	/// @brief Parses the FastDesignTags
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);

	/// @brief Return the name of this mover.
	virtual std::string get_name() const;

	/// @brief return a fresh instance of this class in an owning pointer
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the FastDesign. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	/// @brief initializes cache of allowed amino acids
	void initialize_aa_cache( core::pose::Pose const & pose );

protected:

	/// @brief sets constraint weights -- used with constraint ramping
	virtual void
	set_constraint_weight(
		core::scoring::ScoreFunctionOP local_scorefxn,
		core::scoring::EnergyMap const & full_weights,
		core::Real const weight,
		core::pose::Pose & pose ) const;

	/// @brief outputs a pdb file with sequential numbering, with the specified prefix
	void dump_pdb( core::pose::Pose const & pose, std::string const & prefix = "dump" );

private:

	/// @brief checks each mutation at each designable position, and disallows those mutations that hurt the filter score
	utility::vector1< utility::vector1< bool > >
	check_and_disallow_mutations_by_filter( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, protocols::filters::FilterOP filter ) const;


	/// @brief increments the number of times a particular residue (and aa) have been designed on
	void increment_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos );

	/// @brief checks how many times the residue at seqpos has been redesigned. If the value is larger than than the maximum allowed number of redesigns, return false, otherwise return true
	bool check_num_redesigns( core::pose::Pose const & pose, core::Size const seqpos ) const;

	/// @brief Creates and returns a pointer to a taskoperation that designs the worst regions of a protein only. Returns NULL if we don't need the taskop
	protocols::toolbox::task_operations::DesignAroundOperationOP
	create_worst_region_operation( core::pose::Pose const & pose );

private:   // options

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

	/// @brief Should we clear designable residues (sets them to ALA so that they can be rebuilt)
	bool clear_designable_residues_;

	/// @brief scorefunction used to score and rank poses (NOT for packing)
	core::scoring::ScoreFunctionOP rank_scorefxn_;

	/// @brief should we prevent the native amino acid, if it's one of the worst in the pose?
	bool disallow_bad_natives_;

	/// @brief the maximum number of times we should redesign a particular position and amino acid before disallowing that amino acid
	core::Size max_redesigns_;

private:   // other data

	core::Size dump_counter_; // counter for number of pdb dumps
	utility::vector1< FilterParams > filters_; // filters to be applied at each step
	utility::vector1< utility::vector1 < bool > > allowed_aas_; // list of allowed amino acids at each position
	utility::vector1< utility::vector1< core::Size > > num_redesigns_; // list of the number of times each position has been redesigned
	//flxbb::filters::DesignByFragmentQualityOperationOP frag_qual_op_;
	core::Size regions_to_design_;
	core::Size run_count_;
	std::string cached_sequence_;
	utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > rcgs_;
	//RestrictWorstRegionOP worst_region_mover_;
};


} // denovo_design
} // devel

#endif
