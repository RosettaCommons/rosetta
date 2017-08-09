// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/Hydrate.cc
/// @brief The Hydrate Protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_hydrate_functions_HH
#define INCLUDED_protocols_hydrate_hydrate_functions_HH

// Protocols
#include <protocols/hydrate/Hydrate.hh>
#include <protocols/moves/Mover.hh>

// Core
#include <core/chemical/AA.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh> //yumeng 03/20/13
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

//wym
#include <core/pack/rotamer_set/RotamerSets.hh>

namespace protocols {
namespace hydrate {

// Read hyfile;
// Hyfile should be a list with two columns:
// 1st with one of the keywords://      "HYD" to hydrate that residue (only residues)
//      "ENF" to enforce that water molecule, to ensure it will be present at the end of the simulation (only water molecules)
// 2nd with the residue or water molecule number
void
read_hyfile(
	std::string const & hyfile,
	utility::vector1< bool > & enforced_V,
	utility::vector1< bool > & hydrate_V
);

// Add water molecules to the system, in the residues specified in the hyfile
void
hydrate_hyfile(
	core::pose::Pose & pose,
	utility::vector1< bool > const & hydrate_V,
	std::string const & resfile
);

// Move water molecule at anchor position, important for neighbor calculation
void
place_de_novo_wat_at_anchor(
	core::pose::Pose & pose
);

// When using explicit water molecules, the pose needs WATER_PACKING_INFO to handle them,
// also add de novo water molecules
void
set_water_info_and_add_de_novo_water(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & scorefxn
);

// This function defines an atom as hydratable when it has neighboring positions that are:
// large enough for a water molecule and they are inside.
bool
atom_is_hydratable(
	core::pose::Pose const & pose,
	core::Size const residue,
	std::string const & atom
);

// This function defines an atom as hydratable when it has neighboring positions that are:
// large enough for a water molecule and they are inside.
bool
atom_is_hydratable(
	core::pose::Pose const & pose,
	core::Size const residue,
	core::Size const atom
);

// This function adds de novo water molecules to all cavities with a potential anchor atom (polar, not engaed in hb)
void
hydrate_cavities(
	core::pose::Pose & pose
);

// Determines if the location of the core::Vector is "inside" the pose or not
bool
is_inside(
	core::pose::Pose const & pose,
	core::Vector const & xyz
);

// This function moves a fraction of the water molecules away from the protein, not to include them in the
// first step, with rotamers oprimizing two hb (double edge water, dew). They stay in a "buffer" for the second step
// where thay will build rotamers optimizing just one hb (single edge, sew)
void
set_dew_waters_not_to_be_included(
	core::pose::Pose & pose,
	core::Real const partial_hydrate_dew
);

//calculate whether the residue is near water
//considering all heavy atoms
bool  // yumeng
residue_near_water(
	core::pose::Pose const & pose,
	core::Size const ii
);

// This function sets the task to be used with the packer, and the movemap to be used with the minimizer
void
set_task_and_movemap(
	core::pose::Pose const & pose,
	std::string const & protein_flexibility,
	core::pack::task::PackerTaskOP & task,
	core::kinematics::MoveMap & mm,
	bool const minimize_bb_where_packing
);

// Removes water molecules with high energy. It accounts for water specific energy corrections
void
remove_high_energy_water_molecules(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & scorefxn
);

// Calculates the water overcoordinated or bifurcated hydrogen bonds correction
void
calculate_water_overcoordinated_hb_correction(
	core::pose::Pose const & pose,
	utility::vector1< core::Real > & water_hb_correction
);

// This function has all water molecules forced to stay near the protein (active)
void
enforce_all_waters(
	core::pose::Pose & pose
);

// Set the task for packing de novo water molecules that are away from the protein at this point by using just
// one optimized hydrogen bond (single edge, sew)
void
get_ready_for_sew_packing(
	core::pose::Pose & pose,
	core::pack::task::PackerTaskOP & task
);

// All water molecules will not have an anchor atom associated to them and they will not be enforced to stay
// active (near protein)
void
remove_all_anchors_and_ENF(
	core::pose::Pose & pose
);

// Set a bb movemap for the minimizer according to a mini_backbone_file
void // yumeng
set_bb_movemap(
	core::pose::Pose const & pose,
	std::string const & mini_backbone_file_name,
	core::kinematics::MoveMap & mm
);

// Remove all water molecules not considered buried
void
remove_non_buried_wat(
	core::pose::Pose & pose
);

// Add wat_overcoor_hb into the "scores" set and add it to the total_score. It will be printed in the pdb and sc file.
void
add_water_overcoordinated_hb_score(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn
);

// Display all the hydrogen bonds involving water molecules
void
show_water_hb_network(
	core::pose::Pose const & pose
);

// Output no_water hbond and water-specific hbond energies in the end of pdb file   wym
void
water_specific_hbond_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn
);

// Set task considering a resfile and packing also all de novo water molecules
// deprecated
void
set_task_with_de_novo_water_using_resfile(
	core::pose::Pose & pose,
	std::string resfile,
	core::pack::task::PackerTaskOP & task
);

// Function to print out which residues to zero out fa_sol (based on if a residue is near water)
void
print_residues_near_water(
	core::pose::Pose const & pose
);

// Function to add in a single water far away from structure
void
append_single_far_away_water(
	core::pose::Pose & pose
);



} // close hydrate
} // close hydrate

#endif
