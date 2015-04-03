// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/pmut_scan/AlterSpecDisruptionDriver.hh
/// @brief A protocol that tries to find interface-disrupting mutations as phase 1 of an altered-specificity protocol
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_pmut_scan_AlterSpecDisruptionDriver_HH
#define INCLUDED_protocols_pmut_scan_AlterSpecDisruptionDriver_HH

// Project Headers
#include <protocols/pmut_scan/Mutant.fwd.hh>
#include <protocols/pmut_scan/PointMutScanDriver.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>

// Utility headers

// ObjexxFCL header

// C++
#include <set>

namespace protocols {
namespace pmut_scan {

/// @details this subclass of Ron's PointMutScanDriver exists to tweak one aspect of his code: instead of looking for total-energy-stabilizing mutations, it looks for binding-energy-DEstabilizing point & pair mutations.
class AlterSpecDisruptionDriver: public PointMutScanDriver {

public:
	AlterSpecDisruptionDriver( utility::vector1< std::string > & pdb_file_names, bool double_mutant_scan, std::string list_file, bool output_mutant_structures );
	~AlterSpecDisruptionDriver();

  /// @brief return a score that is a ddG of binding, rather than a ddG of the interface.  It returns a reversed value because this class wants to find DEstabilizing mutations.
  virtual core::Energy score(core::pose::Pose & pose);

	/// @brief offers a chance for child classes to inject mutant selection logic
	virtual bool reject_mutant( Mutant const & mutant, core::pose::Pose const & pose );

private:

	/// @brief reject Mutant based on chain IDs of constituent mutations
	bool reject_on_chains( Mutant const & mutant );

	/// @brief reject Mutant based on interface-ness of constituent mutations
	bool reject_on_interface( Mutant const & mutant, core::pose::Pose const & pose );

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	std::set<core::Size> interface_set_;

}; // class AlterSpecDisruptionDriver

} // namespace pmut_scan
} // namespace protocols

#endif //INCLUDED_protocols_pmut_scan_AlterSpecDisruptionDriver_HH
