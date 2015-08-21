// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/metal_interface/AddZincSiteConstraints.hh
/// @brief   Adds metal binding site geometry constraints to pose.
/// @details I typically work with zinc binding sites having tetrahedral coordination geometry, in which case there are 4 distance constraints, 4 angle constraints, 6 tetrahedral constraints, and 4 dihedral constraints.  The code is flexibile enough to accommodate 2, 3, or 4-residue metal binding sites.
/// @author Bryan Der

#ifndef INCLUDED_protocols_metal_interface_AddZincSiteConstraints_HH
#define INCLUDED_protocols_metal_interface_AddZincSiteConstraints_HH

#include <protocols/metal_interface/AddZincSiteConstraints.fwd.hh>
#include <protocols/metal_interface/MetalSiteResidue.fwd.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace metal_interface {

/// @brief Add metalsite geometry constraints to pose
class AddZincSiteConstraints : public utility::pointer::ReferenceCount {

public:

	AddZincSiteConstraints( utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr );
	virtual ~AddZincSiteConstraints();

	virtual void add_constraints ( pose::Pose & pose );
	virtual void evaluate_constraints( pose::Pose const & pose );
	virtual void view_constraints_in_pymol( pose::Pose const & pose );
	virtual void output_constraints_file( pose::Pose const & pose );

private:
	std::string pdbname_;
	utility::vector1< protocols::metal_interface::MetalSiteResidueOP > msr_;

	utility::vector1< core::scoring::constraints::AtomPairConstraintCOP > distance_constraints_;
	utility::vector1< core::scoring::constraints::AngleConstraintCOP >    angle_constraints_;
	utility::vector1< core::scoring::constraints::DihedralConstraintCOP > dihedral_constraints_;
	utility::vector1< core::scoring::constraints::AngleConstraintCOP >    tetrahedral_constraints_;


};//end AddZincSiteConstraints


}//namespace metal_interface
}//namespace protocols

#endif // INCLUDED_protocols_metal_interface_AddZincSiteConstraints_HH

