// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceStubMover.hh
/// @brief definition of classes for grafting hotspots into a pose
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlaceUtils_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlaceUtils_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>


// C++ headers

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

bool
test_res_res_aln( core::conformation::Residue const & res1, core::conformation::Residue const & res2, core::Real & C_N_angle, core::Real & CB_CA_angle  );

/*core::scoring::constraints::ConstraintCOPs
add_coordinate_constraints( core::pose::Pose & pose, core::Size const host_chain, core::Size const resnum, core::Real const coord_sdev, core::scoring::constraints::HarmonicFuncOP & coord_cst_func );
*/

core::scoring::constraints::ConstraintCOPs
add_coordinate_constraints( core::pose::Pose & pose, core::conformation::Residue const source, core::Size const host_chain, core::Size const resnum, core::Real const coord_sdev, core::scoring::func::HarmonicFuncOP & coord_cst_func );

void
generate_taskfactory_and_add_task_awareness( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers, basic::datacache::DataMap & data, core::pack::task::TaskFactoryOP & task_factory );

std::string nearest_atom_for_constraint( core::conformation::Residue const residue );

/// @brief find the nearest residue to a coordinate
core::Size
find_nearest_residue_to_coord( core::pose::Pose const & pose, numeric::xyzVector< core::Real > const coord, core::Size const host_chain );

/// @brief a utility function for parsing stubset information from a tag
utility::vector1< std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > >
parse_stub_sets( utility::tag::TagCOP tag, core::pose::Pose const & pose, core::Size const host_chain, basic::datacache::DataMap data );

core::scoring::ScoreFunctionOP make_stub_scorefxn();
} //movers
} //protein_interface_design
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_PlaceUtils_HH*/
