// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/TargetClashEnergy.cc
/// @brief  for clashing check
/// @author Longxing (longxing@uw.edu)


// Unit headers
#include <protocols/scoring/methods/TargetClashEnergy.hh>
#include <protocols/scoring/methods/TargetClashEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>

//symmetry

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

#include <core/id/AtomID_Map.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <numeric/geometry/hashing/MinimalClashHash.hh>
#include <core/pose/xyzStripeHashPose.hh>

#include <basic/Tracer.hh>


using namespace core;
using namespace core::scoring;
using namespace core::scoring::methods;

namespace protocols {
namespace scoring {
namespace methods {

static basic::Tracer TR("protocols.scoring.methods.TargetClashEnergy");

/// @details This must return a fresh instance of the TargetClashEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
TargetClashEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< TargetClashEnergy >(options);
}

ScoreTypes
TargetClashEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( target_clash );
	return sts;
}


/// c-tor
TargetClashEnergy::TargetClashEnergy( core::scoring::methods::EnergyMethodOptions const & options ) :
	parent( utility::pointer::make_shared< TargetClashEnergyCreator >() ),
	target_clash_pdb_( options.target_clash_pdb()),
	voxel_initialized_(false),
	clash_atom_scale_(0.5),
	clash_check_resn_("ALA"),
	context_clash_hash_(nullptr)
{
	initiate_voxel();
}


/// clone
EnergyMethodOP
TargetClashEnergy::clone() const
{
	return utility::pointer::make_shared< TargetClashEnergy >(*this);
}


void
TargetClashEnergy::initiate_voxel(){
	if ( target_clash_pdb_ == "" ) return;

	if ( ! voxel_initialized_ ) {

		pose::Pose context_pose;
		core::import_pose::pose_from_file(context_pose, target_clash_pdb_, false);

		// TODO: poly VAL ALA GLY
		if ( "NAT" != clash_check_resn_ ) {
			core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
			core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map(clash_check_resn_) );
			for ( core::Size ir = 1; ir <= context_pose.size(); ++ir ) {
				if (   context_pose.residue(ir).is_protein()   ) continue;
				if (   context_pose.residue(ir).name3()=="GLY" ) continue;
				if (   context_pose.residue(ir).name3()=="PRO" ) continue;
				if (   context_pose.residue(ir).name3()=="CYD" ) continue;
				context_pose.replace_residue( ir, *ala, true );
			}
		}

		core::id::AtomID_Map<core::Real> context_heavy_atoms = core::pose::make_atom_map( context_pose, core::pose::PoseCoordPickMode_HVY );
		utility::vector1<numeric::geometry::hashing::Ball> context_balls;
		core::pose::xyzStripeHashPose::extract_pose_balls( context_pose, context_balls, context_heavy_atoms );

		// Ok, so here's the deal with the atomic radii we choose here. 2.0 is carbon, 1.6 is oxygen.
		// But Daniel's clashing function uses half the atomic radii. So 1.6 is too small, but the
		// "correct" number might be something like 0.9
		context_clash_hash_ = utility::pointer::make_shared<numeric::geometry::hashing::MinimalClashHash>( 0.25f, clash_atom_scale_*1.8f, context_balls );



		voxel_initialized_ = true;
	}
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the RMS difference between native_pose_ (provided by
/// the option -in::file::native and the given Pose. The actual energy calculation
/// is the difference between the RMSD and the target RMSD. Target RMSD is specified
/// the option -score::rms_target.


void
TargetClashEnergy::setup_for_scoring( pose::Pose &, ScoreFunction const & ) const
{

}

void
TargetClashEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	if ( ! voxel_initialized_ ) {
		if ( target_clash_pdb_ == "" ) {
			utility_exit_with_message( "protocols.scoring.methods.TargetClashEnergy: Error! You must set: target_clash_pdb!");
		} else {
			utility_exit_with_message("protocols.scoring.methods.TargetClashEnergy: WTF, what happened?? Some bug??");
		}
	}
	core::id::AtomID_Map<core::Real> ca_atoms = core::pose::make_atom_map( pose, core::pose::PoseCoordPickMode_CA ) ;
	utility::vector1<numeric::geometry::hashing::Ball> balls;
	core::pose::xyzStripeHashPose::extract_pose_balls( pose, balls, ca_atoms );
	core::Size clash_count = context_clash_hash_->clash_check_balls( balls );
	emap[ target_clash ] = clash_count;
}

core::Size
TargetClashEnergy::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
