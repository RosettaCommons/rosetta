// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/kinematic_closure/KinematicPerturber.cc
/// @brief  implementations for KinematicPerturbers used by the kinematic mover
/// @author Florian Richter, floric@u.washington.edu, march 2009
/// @author Rhiju Das, rhiju@stanford.edu, 2011 -- options of cis/trans prolines, and turn off ca bond geometry variation.

//Unit headers
#include <protocols/moves/kinematic_closure/KinematicPerturber.hh>

// Project headers
#include <protocols/moves/KinematicMover.hh>

// Rosetta Headers
#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/AtomTree.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

// option key includes
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/kinematics/MoveMap.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {
namespace kinematic_closure {

static numeric::random::RandomGenerator RG(43134);
static basic::Tracer TR("protocols.moves.kinematic_closure.KinematicPerturber");

KinematicPerturber::KinematicPerturber()
	: max_sample_iterations_( basic::options::option[ basic::options::OptionKeys::loops::max_kic_perturber_samples ]() )
{}

KinematicPerturber::~KinematicPerturber(){}

void KinematicPerturber::set_movemap( core::kinematics::MoveMapCOP mm ) { movemap_ = mm; }

core::kinematics::MoveMapCOP KinematicPerturber::get_movemap() const { return movemap_; }

void
KinematicPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const &, //bond_ang,
	utility::vector1<core::Real> const &, //bond_len,
	bool //closure_successful
) const
{

	core::Size start( kinmover_->start_res() );

	for( core::Size res = 0; res < kinmover_->segment_length(); res++ ){

		pose.set_phi( start + res, torsions[ (3*(res+1)) + 1 ] );
		pose.set_psi( start + res, torsions[ (3*(res+1)) + 2 ] );
		if ( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ) pose.set_omega( start + res, torsions[ (3*(res+1)) + 3 ] );

	}
} //set_pose_after_closure

///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////TorsionSamplingKinematicPerturber////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

TorsionSamplingKinematicPerturber::TorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
	: KinematicPerturber(),
		vary_ca_bond_angles_(false),
		sample_vicinity_(false),
		degree_vicinity_(1.0 ),
		sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
		rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
{ set_kinmover( kinmover_in ); }

TorsionSamplingKinematicPerturber::~TorsionSamplingKinematicPerturber(){}

///@details randomly varies the torsions (and possibly the bond angles) for the loop.  Respects a MoveMap, if present, for torsions.  Does NOT respect the movemap for angles; does NOT cause any sort of interactions between the MoveMap and the KinematicMover's pivot residues.
void
TorsionSamplingKinematicPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
 	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & //bond_len
) const
{
	core::kinematics::MoveMapCOP mm(get_movemap());

	if( vary_ca_bond_angles_ ){

		core::Size pvatom3( (3* (kinmover_->segment_length() + 1)) - 1 );

		core::Real bangle_min( kinmover_->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_->BANGLE_SD() );

		//what is this iterating over?
		for( Size i = 5; i <= pvatom3; i+=3 ) {
			bond_ang[ i ] = bangle_min + RG.uniform() * bangle_sd;
		}
	}



	core::Size tor_end = torsions.size() - 3;

	for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){
		//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;

		if(!mm || mm->get_bb(cur_res)){ //if movemap is null, or (if not null) if movemap says mobile

			core::Real rama_phi, rama_psi;

			if( sample_vicinity_ ){
				rama_phi = pose.phi( cur_res ) + degree_vicinity_ * RG.gaussian();
				rama_psi = pose.psi( cur_res ) + degree_vicinity_ * RG.gaussian();
			}
			else rama_.random_phipsi_from_rama(pose.aa(cur_res), rama_phi, rama_psi);

			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi

			i++; // leave omega alone

		} else {
			i += 3; //ensure i indexing increases
		}

	}

	if (  sample_omega_for_pre_prolines_ && !sample_vicinity_ ) {
		// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.

		static const core::Real OMEGA_MEAN( 179.8 );

		for( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ){

			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {

				i++; //phi
				i++; //psi
				torsions[i++] = ( static_cast<int>( RG.uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)

			} else {
				i += 3;
			}
		}

	}


} //perturb_chain


void
TorsionSamplingKinematicPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
 	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful // what is this used for?
) const
{

	//	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );

	if( vary_ca_bond_angles_ ){

		core::Real offset( 0.0 );
		for (Size res=kinmover_->start_res(), atom=5; res<= kinmover_->end_res(); res++, atom+=3) {

			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
				numeric::conversions::radians(180 - bond_ang[atom]));

		}
	}
} //TorsionSamplingKinematicPerturber::set_pose_after_closure(


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////TorsionSweepkingKinematicPerturber//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

TorsionSweepingKinematicPerturber::TorsionSweepingKinematicPerturber()
	: KinematicPerturber()
{}

TorsionSweepingKinematicPerturber::~TorsionSweepingKinematicPerturber(){}

void TorsionSweepingKinematicPerturber::set_nonpivot_res_to_sweep( utility::vector1< Size > const & resids )
{
	nonpivot_res_to_sweep_ = resids;
}

void TorsionSweepingKinematicPerturber::set_nonpivot_bb_torsion_id( utility::vector1< Size > const & bbtorids )
{
	assert( nonpivot_res_to_sweep_.size() == bbtorids.size() );
	sweep_torsion_ids_ = bbtorids;
}

void TorsionSweepingKinematicPerturber::set_sweep_start_angle( utility::vector1< core::Real > const & angles_in_degrees )
{
	assert( nonpivot_res_to_sweep_.size() == angles_in_degrees.size() );
	sweep_nonpivot_torsion_starts_ = angles_in_degrees;
}

void TorsionSweepingKinematicPerturber::set_sweep_step_size( utility::vector1< core::Real > const & angle_steps_in_degrees )
{
	assert( nonpivot_res_to_sweep_.size() == angle_steps_in_degrees.size() );
	sweep_step_sizes_ = angle_steps_in_degrees;
}

/// @details Initializes the LexicographicalIterator
void TorsionSweepingKinematicPerturber::set_sweep_nsteps( utility::vector1< core::Size > const & nsteps )
{
	assert( nonpivot_res_to_sweep_.size() == nsteps.size() );
	sweep_iterator_.set_dimension_sizes( nsteps );
}


void
TorsionSweepingKinematicPerturber::perturb_chain(
	core::pose::Pose const &, //pose,
	utility::vector1<core::Real> & torsions,
		utility::vector1<core::Real> &,// bond_ang,
		utility::vector1<core::Real> & //bond_len
) const
{

	if ( sweep_iterator_.at_end() ) {
		utility_exit_with_message("TorsionSweepingKinematicPerturber asked to perturb chain even though sweep iterator is at end."); }

	core::Size start( kinmover_->start_res() );

	for ( Size ii = 1; ii <= nonpivot_res_to_sweep_.size(); ++ii ) {
		Size torsion_ind = 3 * ( nonpivot_res_to_sweep_[ ii ] - start + 1 ) + sweep_torsion_ids_[ ii ];
		torsions[ torsion_ind ] = sweep_nonpivot_torsion_starts_[ ii ] + sweep_step_sizes_[ ii ] * ( sweep_iterator_[ ii ] - 1 );
		//std::cout << " " <<  nonpivot_res_to_sweep_[ ii ] << " " << sweep_torsion_ids_[ ii ] << " " << torsion_ind << " = " << dt_ang[ torsion_ind ];
	}
	//std::cout << std::endl;
	++sweep_iterator_;
} //TorsionSweepingKinematicPerturber::perturb_chain(

} // end namespace kinematic_closure
} // end namespace moves
} // end namespace protocols
