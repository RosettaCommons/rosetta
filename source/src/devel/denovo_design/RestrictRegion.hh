// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/RestrictRegion.hh
/// @brief Tom's Denovo Protocol. This is freely mutable and used for playing around with stuff
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_RestrictRegion_hh
#define INCLUDED_devel_denovo_design_RestrictRegion_hh

// Unit headers
#include <devel/denovo_design/RestrictRegion.fwd.hh>
#include <devel/denovo_design/task_operations/DesignBySecondaryStructure.hh>
// Project headers
#include <protocols/moves/Mover.hh>

#include <core/chemical/AA.hh>

#include <core/fragment/Frame.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

//// C++ headers
#include <string>

#include <core/io/silent/silent.fwd.hh>
#include <utility/vector1.hh>



namespace devel {
namespace denovo_design {

class RestrictRegion : public protocols::moves::Mover {
public:
  /// @brief Initialize RestrictRegion
  RestrictRegion();

	/// @brief copy constructor
	RestrictRegion( RestrictRegion const & rval );

  /// @brief virtual constructor to allow derivation
	virtual ~RestrictRegion();

  /// @brief Parses the RestrictRegion tags
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

  /// @brief Apply the RestrictRegion. Overloaded apply function from mover base class.
	virtual void apply( core::pose::Pose & pose );

	/// @brief static method that tells the last residue affected by any RestrictRegion mover. This method also makes sure that the pose length hasn't changed. If this is called before any instance of the mover::apply() is called, returns 0.
	static utility::vector1< core::Size > const &
	last_residues_restricted();

	/// @brief static method that tells the last pose processed by any RestrictRegion mover. If this is called before any instance of the mover::apply() is called, returns 0.
	static core::pose::PoseCOP
	previous_pose();

	/// @brief Permanently restricts design at the specified position such that the existing amino acid is not allowed.
	static void
	permanently_restrict_aa( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, core::Size const seqpos );

private: // private methods
	/// @brief initialize the resfile -- takes an input resfile, but converts it to a resfile in the output directory
	void initialize_resfile( std::string const & orig_resfile );

	/// @brief writes a resfile from the task
	void write_resfile( core::pose::Pose const & pose, core::pack::task::PackerTaskCOP task ) const;

	/// @brief restricts design at the specified position such that the existing amino acid is not allowed, and writes this action to a resfile so that other movers can use it
	bool restrict_aa( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, core::Size const seqpos );

	/// @brief tells whether a given fragment is compatible with the task
	/// basically, this function makes sure every amino acid in the sequence is allowed at the position given
  core::Size compatible_with_task( core::pack::task::PackerTaskOP task,
														 core::Size const frag_id,
														 core::fragment::FrameOP frame ) const;


	/// @brief restricts design at the specified position such that only the existing amino acid is allowed, and writes this to a resfile
	void
	preserve_residue( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task,  core::Size const seqpos ) const;

	/// @brief function that checks to see whether an amino acid is restricted at a certain position
	bool is_restricted( char const aa, core::Size const seqpos ) const;

	/// @brief function that tells the one letter codes of residues allowed at a certain position
	std::string
	residues_allowed( core::pack::task::PackerTaskCOP task, core::Size const seqpos ) const;

private:   // options
	/// @brief type of calculation to do default="score"
	std::string type_;

	/// @brief should we permanently restrict residues when a successful trial occurs?
	bool permanent_restriction_;

	/// @brief resfile to be modified/used
	std::string resfile_;

	/// @brief blueprint file to be used for psipred
	std::string blueprint_file_;

	/// @brief command to run for psipred predictions (required for frag_qual and psipred)
	std::string psipred_cmd_;

	/// @brief maximum number of tryptophans allowed -- if this number is reached, TRP will be disallowed
	bool enable_max_trp_;
	core::Size max_trp_;

private:   // other data
	/// @brief task factory
	core::pack::task::TaskFactoryOP task_factory_;

	/// @brief Stores lists of restricted residues for each sequence position
	utility::vector1< std::string > restricted_residues_;

	/// @brief Scorefunction to use for score-based worst region design
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief number of regions to mutate
	core::Size regions_to_mutate_;

	/// @brief previous type of worst region restriction applied -- useful for bookeeping to see which metrics are successful
	std::string last_type_;

	/// @brief stats for each possible metric type -- the pair is the number of accepts (first) total number of trials (second)
	std::map< std::string, std::pair< core::Size, core::Size > > metric_stats_;

	/// @brief Vector of residues that have been restricted at each position
	static utility::vector1< std::string > permanently_restricted_residues_;

	/// @brief This is the pose from the last apply() call -- if it is different from the current pose, the changes will be made permanent, otherwise they will be ignored.
	static core::pose::PoseOP previous_pose_;

	/// @brief Stores the last residue restricted by any RestrictRegion mover
	static utility::vector1< core::Size > last_residues_restricted_;

/// @brief strores a vector or the possible highestEnergyRegionOperations
	utility::vector1<task_operations::HighestEnergyRegionOperationOP> highestEnergyRegionOperation_ops_;
};

/// @brief function that tells whether a given amino acid is allowed at a certain position
bool
residue_is_allowed( core::pack::task::PackerTaskCOP task, core::Size const seqpos,  core::chemical::AA const aa );

/// @brief alter the resfile line such that the amino acid given is listed as allowed and return the modified version
std::string
allow_in_resfile_line( std::string const & cmd_orig, char const aa );

}
} // devel

#endif
