// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/denovo_design/task_operations/HighestEnergyRegion.hh
/// @brief Design regions with low per-residue energy
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_devel_denovo_design_task_operations_highestenergyregionoperation_hh
#define INCLUDED_devel_denovo_design_task_operations_highestenergyregionoperation_hh

// unit headers
#include <devel/denovo_design/task_operations/HighestEnergyRegion.fwd.hh>
// protocol headers

// project headers
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

namespace devel {
namespace denovo_design {
namespace task_operations {

class HighestEnergyRegionOperation : public core::pack::task::operation::TaskOperation {
public:
  /// @brief default constructor
  HighestEnergyRegionOperation();

  /// @brief copy constructor
  HighestEnergyRegionOperation( HighestEnergyRegionOperation const & rval );

  /// @brief destructor
  virtual ~HighestEnergyRegionOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

  /// @brief apply
  virtual void apply( Pose const & pose, core::pack::task::PackerTask & task ) const;

	/// @brief Runs the calculation and caches residues to design
	virtual void cache_result( core::pose::Pose const & pose );

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "Score"; }

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

public:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map );

	/// @brief tells this task operation whether it should use the cache when it is applied
	void set_use_cache( bool const use_cache );

	/// @brief sets the score function
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );

	/// @brief Gets a list of residues to design, but uses cached result if present
	utility::vector1< core::Size >
	residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief accessor function to access the cached pose, if it exists
	core::pose::PoseCOP
	cached_pose() const;

	/// @brief Gets list of allowed amino acids at position resi
	utility::vector1< bool > const &
	allowed_aas( core::Size const resi ) const;

	/// @brief Sets amino acid aa as disallowed at position resi
	void disallow_aa( core::Size const resi, char const aa );

	/// @brief Sets the number of regions to select for design
	void set_regions_to_design( core::Size const num_regions );

	/// @brief Returns the number of regions selected for design
	core::Size regions_to_design() const;

	/// @brief set the shell size for the regions to be examined and designed
	void set_region_shell( core::Real const region_shell );

	/// @brief accessor function for region_shell member variable
	core::Real region_shell() const;

	/// @brief mutator for repack_non_selected -- sets behavior for whether we should repack or fix areas outside the design region
	void repack_non_selected( bool const repack_non_selected );

	/// @brief accessor for repack_non_selected member var
	bool repack_non_selected() const;

private:
	/// @brief initializes cache of allowed amino acids
	void initialize_aa_cache( core::pose::Pose const & pose );


private:
	/// @brief Consider all residues within this many angstroms of target residue for energy calculations
	core::Real region_shell_;
	/// @brief tells how many regions should be set to designable. If set to n, the worst n residues by psipred prediction will be set to designable. default=1
	core::Size regions_to_design_;
	/// @brief tells whether residues not selected should be repaced or not (default = false, which means repack everything)
	bool repack_non_selected_;
	/// @brief Tells whether calls to the apply function should recompute energy for regions, or whether the result of the last calculations should be blindly applied. This is useful when a mover which changes the amino acids in the pose (e.g. polyx) is applied within another mover (e.g. flxbbdesign)
	bool use_cache_;
	/// @brief score function to use (default=score12)
	core::scoring::ScoreFunctionOP scorefxn_;
	/// @brief Result of the last call to apply()
	utility::vector1< core::Size > residues_to_design_;
	/// @brief pose by which cached results were generated
	core::pose::PoseCOP cached_pose_;
	/// @brief cached set of amino acids that have been present at the worst agreeing psipred positions
	/// if the "prevent_native_aa" option is used, all amino acids in this list at the position in question will be disallowed
	/// When any position is disallowed, it is added to this data structure
	utility::vector1< utility::vector1< bool > > allowed_aas_;

};

/// @brief class for finding the poorest fragment quality region and redesigning it
class DesignByPackStatOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignByPackStatOperation();

	/// @brief copy constructor
	DesignByPackStatOperation( DesignByPackStatOperation const & rval );

	/// @brief virtual destructor
	virtual ~DesignByPackStatOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "PackStat"; }

private:
};

/// @brief Simply chooses residues randomly for design
class DesignRandomRegionOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignRandomRegionOperation();

	/// @brief virtual destructor
	virtual ~DesignRandomRegionOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "RandomRegion"; }

private:
};

/// @brief class for finding residues with highest closeness centrality and designing them
class DesignByResidueCentralityOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignByResidueCentralityOperation();

	/// @brief virtual destructor
	virtual ~DesignByResidueCentralityOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "ResidueCentrality"; }

private:
};

/// @brief class for finding catalytic residues and designing around them
class DesignCatalyticResiduesOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignCatalyticResiduesOperation();

	/// @brief virtual destructor
	virtual ~DesignCatalyticResiduesOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "CatalyticResidues"; }

private:
};

/// @brief class for finding residues close to buried cavities and designing them
class DesignByCavityProximityOperation : public HighestEnergyRegionOperation {
public:
	/// @brief default constructor
	DesignByCavityProximityOperation();

	/// @brief virtual destructor
	virtual ~DesignByCavityProximityOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Gets a list of residues for design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "CavityProximity"; }

private:

	/// @brief given a cavity and a residue, tells how far Cb of the residue is from the edge of the cavity. Normalizes distance by cavity volume, such that residues around larger cavities should be preferred
	core::Real
	proximity_to_cavity( core::conformation::Residue const & res,
											 core::scoring::packstat::CavityBallCluster const & cluster ) const;
};


} // task_operations
} // denovo_design
} // devel

#endif
