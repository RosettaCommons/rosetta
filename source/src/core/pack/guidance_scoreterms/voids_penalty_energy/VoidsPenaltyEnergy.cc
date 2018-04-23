// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergy.cc
/// @brief An EnergyMethod intended for packing, which penalizes solutions in which the total volume to fill differs greatly
/// from the total volume of the current set of rotamers.
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

// Unit headers
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergy.hh>
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyEnergyCreator.hh>
#include <core/pack/guidance_scoreterms/voids_penalty_energy/VoidsPenaltyVoxelGrid.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <utility/numbers.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/symmetry/util.hh>

// Options system
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// File I/O
#include <basic/database/open.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>

// Other Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace voids_penalty_energy {

#define ENERGY_MULTIPLIER 0.01 //Reduce the strength of this energy a bit.

static basic::Tracer TR("core.pack.voids_penalty_energy.VoidsPenaltyEnergy");

/// @brief This must return a fresh instance of the VoidsPenaltyEnergy class, never an instance already in use.
///
core::scoring::methods::EnergyMethodOP
VoidsPenaltyEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &options ) const
{
	return core::scoring::methods::EnergyMethodOP( new VoidsPenaltyEnergy( options ) );
}

/// @brief Defines the score types that this energy method calculates.
///
core::scoring::ScoreTypes
VoidsPenaltyEnergyCreator::score_types_for_method() const
{
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::voids_penalty );
	return sts;
}

/// @brief Options constructor.
///
VoidsPenaltyEnergy::VoidsPenaltyEnergy ( core::scoring::methods::EnergyMethodOptions const &options ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new VoidsPenaltyEnergyCreator ) ),
	parent2( ),
	cone_dotproduct_cutoff_( options.voids_penalty_energy_cone_dotproduct_cutoff() ),
	cone_distance_cutoff_( options.voids_penalty_energy_cone_distance_cutoff() ),
	containing_cones_cutoff_( options.voids_penalty_energy_containing_cones_cutoff() ),
	voxel_size_( options.voids_penalty_energy_voxel_size() ),
	voxel_grid_padding_( options.voids_penalty_energy_voxel_grid_padding() ),
	disabled_except_during_packing_( options.voids_penalty_energy_disabled_except_during_packing() ),
	rotamer_volumes_(),
	volume_to_fill_(0.0),
	minimizing_(false)
{
	if ( TR.Debug.visible() ) report();
}

/// @brief Copy constructor.
///
VoidsPenaltyEnergy::VoidsPenaltyEnergy( VoidsPenaltyEnergy const &src ) :
	parent1( core::scoring::methods::EnergyMethodCreatorOP( new VoidsPenaltyEnergyCreator ) ),
	parent2( src ),
	cone_dotproduct_cutoff_( src.cone_dotproduct_cutoff_ ),
	cone_distance_cutoff_( src.cone_distance_cutoff_ ),
	containing_cones_cutoff_(src.containing_cones_cutoff_),
	voxel_size_( src.voxel_size_ ),
	voxel_grid_padding_( src.voxel_grid_padding_ ),
	disabled_except_during_packing_( src.disabled_except_during_packing_ ),
	rotamer_volumes_(src.rotamer_volumes_),
	volume_to_fill_(src.volume_to_fill_),
	minimizing_(src.minimizing_)
{}

/// @brief Default destructor.
///
VoidsPenaltyEnergy::~VoidsPenaltyEnergy() {}

/// @brief Clone: create a copy of this object, and return an owning pointer
/// to the copy.
core::scoring::methods::EnergyMethodOP VoidsPenaltyEnergy::clone() const {
	return core::scoring::methods::EnergyMethodOP( new VoidsPenaltyEnergy(*this) );
}

/// @brief VoidsPenaltyEnergy is context-independent and thus indicates that no context graphs need to be maintained by
/// class Energies.
void VoidsPenaltyEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
	//Do nothing.
	return;
}

/// @brief VoidsPenaltyEnergy is version 1.0 right now.
///
core::Size VoidsPenaltyEnergy::version() const
{
	return 1; // Initial versioning
}

/// @brief Actually calculate the total energy
/// @details Called by the scoring machinery.
/// @note VoidsPenaltyEnergy::finalize_total_energy() can return a slightly different energy than was computed
/// during packing.  This is because reachable volume cannot be computed (since we don't have a rotamer set), so
/// total buried volume is used in the calculation.  This will also be a whole-pose calculation, and won't just
/// focus on the designable region (since there's no "designable region" when this function is called).
void VoidsPenaltyEnergy::finalize_total_energy(
	core::pose::Pose &pose,
	core::scoring::ScoreFunction const &/*sfxn*/,
	core::scoring::EnergyMap &totals
) const {
	if ( disabled_except_during_packing_ ) return; //Do nothing if we're scoring and if evaluation except during packing is disabled.
	if ( minimizing_ ) return; //Do nothing during minimization trajectory.

	TR << "Scoring with VoidsPenaltyEnergy." << std::endl;

	VoidsPenaltyVoxelGrid voxel_grid; //Note: this goes out of scope and is destroyed at the end of this function.  It is recalculated with every call!
	configure_voxel_grid( voxel_grid );
	voxel_grid.set_up_voxel_grid_and_compute_burial( pose );
	totals[ core::scoring::voids_penalty ] += pow( voxel_grid.total_buried_volume() - voxel_grid.compute_total_volume_of_current_residues( pose ) , 2 )*ENERGY_MULTIPLIER;
}

/// @brief Calculate the total energy given a vector of const owning pointers to residues.
/// @details Called directly by the ResidueArrayAnnealingEvaluator during packer runs.  Requires
/// that setup_residuearrayannealablenergy_for_packing() be called first.
core::Real
VoidsPenaltyEnergy::calculate_energy(
	utility::vector1< core::conformation::ResidueCOP > const &resvect,
	core::Size const //substitution_position
) const {
	debug_assert(resvect.size() == rotamer_volumes_.size());

	core::Real total_rotamer_volume(0.0);
	for ( core::Size i(1), imax(resvect.size()); i<=imax; ++i ) {
		if ( !rotamer_volumes_[i].count( resvect[i] ) ) continue;
		total_rotamer_volume += rotamer_volumes_[i].at( resvect[i] );
	}

	total_rotamer_volume -= volume_to_fill_;

	return pow(total_rotamer_volume, 2)*ENERGY_MULTIPLIER;
}

/// @brief Get a summary of all loaded data.
///
void VoidsPenaltyEnergy::report() const {
	if ( !TR.Debug.visible() ) return; //Do nothing if I don't have a tracer.

	TR.Debug << std::endl << "Summary of data loaded by VoidsPenaltyEnergy object:" << std::endl;

	//TODO

	TR.Debug << std::endl;

	TR.Debug.flush();

	return;
}

/// @brief Cache data from the pose in this EnergyMethod in anticipation of scoring.
///
void
VoidsPenaltyEnergy::set_up_residuearrayannealableenergy_for_packing (
	core::pose::Pose &pose,
	core::pack::rotamer_set::RotamerSets const &rotamersets,
	core::scoring::ScoreFunction const &/*sfxn*/
) {
	if ( TR.visible() ) TR << "Setting up the VoidsPenaltyEnergy for packing." << std::endl;

	minimizing_ = false;

	// Setup for symmetry:
	core::conformation::symmetry::SymmetryInfoCOP symminfo;
	core::conformation::symmetry::SymmetricConformationCOP symmconf;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		symmconf = utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() );
		runtime_assert( symmconf != nullptr );
		symminfo = symmconf->Symmetry_Info();
		runtime_assert( symminfo != nullptr );
	} else {
		symminfo = nullptr;
		symmconf = nullptr;
	}

	VoidsPenaltyVoxelGrid voxel_grid; //Note: this goes out of scope and is destroyed at the end of packer setup.  Data needed for packing must be cached by the end of this function.
	configure_voxel_grid( voxel_grid );
	voxel_grid.set_up_voxel_grid_and_compute_burial(pose);
	voxel_grid.prune_voxels_for_fixed_residues( pose, rotamersets, symminfo, symmconf );
	voxel_grid.compute_volumes_of_buried_rotamers( pose, rotamersets, rotamer_volumes_, symminfo, symmconf );

	volume_to_fill_ = voxel_grid.reachable_buried_volume();

	return;
}

/// @brief Called at the beginning of atom tree minimization, this method
/// allows the derived class the opportunity to initialize pertinent data
/// that will be used during minimization.  During minimzation, the chemical
/// structure of the pose is constant, so assumptions on the number of atoms
/// per residue and their identities are safe so long as the pose's Energies
/// object's "use_nblist()" method returns true.
/// @details This just disables this score term during minimization, in the case of the VoidsPenaltyEnergy.
void
VoidsPenaltyEnergy::setup_for_minimizing(
	pose::Pose &,//pose,
	core::scoring::ScoreFunction const &,// sfxn,
	kinematics::MinimizerMapBase const &//minmap
) const {
	if ( TR.visible() && !minimizing_ ) TR << "Disabling VoidsPenaltyEnergy during minimization." << std::endl;
	minimizing_ = true;
}

/// @brief Called after minimization.
/// @details Re-enables the score term after minimization.
void
VoidsPenaltyEnergy::finalize_after_minimizing(
	pose::Pose &// pose
) const {
	if ( TR.visible() && minimizing_ ) TR << "Re-enabling VoidsPenaltyEnergy after minimization." << std::endl;
	minimizing_ = false;
}

//////////////////PRIVATE FUNCTIONS////////////////////////////////////

/// @brief given a voxel grid object, set its parameters using stored values (originally
/// from user settings or the options system).
void
VoidsPenaltyEnergy::configure_voxel_grid(
	VoidsPenaltyVoxelGrid &voxel_grid
) const {
	voxel_grid.set_cone_parameters( cone_dotproduct_cutoff_, cone_distance_cutoff_, containing_cones_cutoff_ );
	voxel_grid.set_voxel_size_and_padding( voxel_size_, voxel_grid_padding_ );
}

} // voids_penalty_energy
} // guidance_scoreterms
} // pack
} // core
