// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/bb_sampler/SugarBBSampler.cc
/// @brief Sample dihdrals using sugar_bb data.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>
#include <core/scoring/carbohydrates/util.hh>
#include <core/pose/carbohydrates/util.hh>

#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/angle.functions.hh>

#include <utility/vector1.functions.hh>
#include <utility/string_util.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.bb_sampler.SugarBBSampler" );

namespace protocols {
namespace simple_moves {
namespace bb_sampler {

using namespace core::id;
using namespace core::scoring::carbohydrates;
using namespace core::chemical::carbohydrates;
using namespace core::pose;
using namespace core::pose::carbohydrates;
using core::Angle;
using core::Probability;
using core::conformation::Residue;

SugarBBSampler::SugarBBSampler():
	BBDihedralSampler(core::id::phi_dihedral, probability),
	sampling_step_size_(0.1)
{

}

SugarBBSampler::SugarBBSampler( core::id::MainchainTorsionType torsion_type, BBSampleType sampling_type, core::Real sampling_step_size):
	BBDihedralSampler(torsion_type, sampling_type),
	sampling_step_size_(sampling_step_size)
{

}


SugarBBSampler::~SugarBBSampler(){}

SugarBBSampler::SugarBBSampler( SugarBBSampler const & src ):
	BBDihedralSampler( src ),
	sampling_step_size_(src.sampling_step_size_) {

}


SugarBBSamplerOP
SugarBBSampler::clone() const {
	return SugarBBSamplerOP( new SugarBBSampler( *this ) );
}


core::Real
SugarBBSampler::get_torsion(Pose const & pose, Size resnum ) const
{
	using core::scoring::ScoringManager;

	if ( ! pose.residue( resnum ).is_carbohydrate() ) {
		utility_exit_with_message("Resnum not a carbohydrate residue! "+ utility::to_string( resnum ));
	}

	if ( torsion_type_ == core::id::omega_dihedral ) {
		utility_exit_with_message("Cannot currently sample using SugarBB data for omega torsion.");
	}

	//This test for exocyclic takes too long - store in LinkageInfo
	if ( torsion_type_ == core::id::psi_dihedral && has_exocyclic_glycosidic_linkage( pose, resnum  ) ) {
		std::string msg = "No data for linkage.  Either this is psi and previous residue is not carbohydrate or we do not have a pyranose ring in the previous residue.";
		//Throw here.  We expect an angle, better to use exception then throw a bogus angle.
		//Catch this to keep going.
		TR << msg << std::endl;
		throw utility::excn::EXCN_Msg_Exception( msg );
	}

	CHIEnergyFunction const & sugar_bb =
		ScoringManager::get_instance()->get_CHIEnergyFunction( true /*setup for scoring*/, sampling_step_size_ );

	// Get linkage type.
	CHIEnergyFunctionLinkageType linkage_type =
		get_CHI_energy_function_linkage_type_for_residue_in_pose( torsion_type_, pose, resnum );

	if ( linkage_type == LINKAGE_NA ) {
		std::string msg = "No data for linkage.  Either this is psi and previous residue is not carbohydrate or we do not have a pyranose ring in the previous residue.";
		//Throw here.  We expect an angle, better to use exception then throw a bogus angle.
		//Catch this to keep going.
		TR << msg << std::endl;
		throw utility::excn::EXCN_Msg_Exception( msg );
	}
	CHIDihedralSamplingData const & sampling_data = sugar_bb.get_chi_sampling_data( linkage_type );

	//TR << "Optimizing resnum " << resnum << " dihedral " << Size(torsion_type_ ) << std::endl;

	//Sample Dihedral
	Angle new_dihedral = 0.0;
	if ( sampling_type_ == probability ) {

		numeric::random::WeightedSampler sampler;
		sampler.weights( sampling_data.probabilities );
		new_dihedral = sampling_data.angles[ sampler.random_sample( numeric::random::rg() ) ];
	} else if ( sampling_type_ == minima ) {
		new_dihedral = sampling_data.angles[ utility::arg_max( sampling_data.probabilities ) ];
	}


	//Fix Phi or Psi for L or D
	core::conformation::Residue const & rsd = pose.residue( resnum );
	if ( torsion_type_ == phi_dihedral && rsd.carbohydrate_info()->is_L_sugar() ) {
		// L-Sugars use the mirror image of the score functions.

		new_dihedral = -new_dihedral;

	} else if ( torsion_type_ == psi_dihedral ) {
		core::conformation::Residue const & prev_rsd( pose.residue( find_seqpos_of_saccharides_parent_residue( rsd ) ) );

		// Next, convert the psi to between 0 and 360 (because that's what the function expects).
		new_dihedral = numeric::nonnegative_principal_angle_degrees( new_dihedral );
		if ( prev_rsd.carbohydrate_info()->is_L_sugar() ) {
			// L-Sugars use the mirror image of the score functions.
			new_dihedral = 360 - new_dihedral;
		}
	}

	// Make sure we are between -180 and 180.
	// TODO:We really need to change this to auto-convert in whatever holds torsions...
	return basic::periodic_range( new_dihedral, 360);
}

void
SugarBBSampler::set_torsion_to_pose(Pose & pose, core::Size resnum) const{
	Angle torsion_angle = get_torsion(pose, resnum);
	set_glycosidic_torsion( torsion_type_, pose, resnum, torsion_angle );

}


}  // bb_sampler
}  // simple_moves
}  // protocols
