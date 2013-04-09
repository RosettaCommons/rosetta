// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/filters/DesignBySecondaryStructure.hh
/// @brief Design residues with secondary structures that don't match the desired secondary structure
/// specified in the blueprint.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_flxbb_filters_designbysecondarystructure_hh
#define INCLUDED_protocols_flxbb_filters_designbysecondarystructure_hh

// unit headers
#include <devel/denovo_design/DesignBySecondaryStructure.fwd.hh>

// protocol headers
#include <protocols/denovo_design/filters/PsiPredInterface.fwd.hh>
#include <devel/denovo_design/HighestEnergyRegion.hh>


// project headers
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

namespace protocols {
namespace flxbb {
namespace filters {

class DesignBySecondaryStructureOperation : public flxbb::filters::HighestEnergyRegionOperation {
public:
  /// @brief default constructor
  DesignBySecondaryStructureOperation();

  /// @brief value constructor
  DesignBySecondaryStructureOperation( std::string const bp_file, std::string const cmd, core::Real const design_shell, core::Real const repack_shell, core::Size const nres_to_design, bool const prevent_native, bool const prevent_bad_point_mutations );

  /// @brief copy constructor
  DesignBySecondaryStructureOperation( DesignBySecondaryStructureOperation const & rval );

  /// @brief destructor
  virtual ~DesignBySecondaryStructureOperation();

  /// @brief make clone
  virtual core::pack::task::operation::TaskOperationOP clone() const;

  /// @brief apply
  virtual void apply( Pose const & pose, core::pack::task::PackerTask & task ) const;

	/// @brief Runs the calculation and caches residues to design
	virtual utility::vector1< core::Size >
	get_residues_to_design( core::pose::Pose const & pose ) const;

	/// @brief Returns the name of the class
	virtual std::string get_name() const { return "DesignBySecondaryStructureOperation"; }

public:
  void parse_tag( utility::tag::TagPtr tag );

	/// @brief computes and calculates the predicted secondary structure string
	std::string compute_and_cache_ss( core::pose::Pose const & pose );

	/// @brief sets the number of residues to select
	void set_nres_to_design( core::Size const nres_to_design );

	/// @brief sets the shell (in angstroms) around the target residue(s) to design
	void set_design_shell( core::Real const shell );

	/// @brief sets the shell (in angstroms) around the target residue(s) to repack
	void set_repack_shell( core::Real const shell );

	/// @brief opens the passed blueprint file and determines the desired secondary structure
	void initialize_blueprint_ss( std::string const blueprint_file );

private:
	std::string blueprint_ss_; // the blueprint secondary structure
	std::string pred_ss_; // cache of the predicted secondary structure
	utility::vector1< core::Real > psipred_prob_; // cache of the predicted probability of secondary structure

	/// @brief design all residues within this many angstroms of residues that don't match SS
	core::Real design_shell_;
	/// @brief repack all residues within this many angstroms of residues that don't match SS
	core::Real repack_shell_;
	/// @brief tells how many residues should be set to designable. If set to n, the worst n residues by psipred prediction will be set to designable. If set to 0, any residues that have non-matching secondary structure by psipred will be designed.
	core::Size nres_to_design_;
	/// @brief If set, the amino acid residue already in the pose will be disallowed default=false
	bool prevent_native_aa_;
	/// @brief If set, all mutations at all positions will be scanned in one pass, and those that cause worse psipred secondary structure agreement will be disallowed (default=false)
	bool prevent_bad_point_mutants_;
	/// @brief should we use the cached value if it exists, or call psipred again?
	bool use_cache_;
	/// @brief the object which directly communicates with psipred and parses psipred output
	denovo_design::filters::PsiPredInterfaceOP psipred_interface_;
};


} // filters
} // flxbb
} // protocols

#endif
