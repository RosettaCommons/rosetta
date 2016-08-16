// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.cc
/// @brief  implementations for KinematicPerturbers used by the kinematic mover
/// @author Florian Richter, floric@u.washington.edu, march 2009
/// @author Rhiju Das, rhiju@stanford.edu, 2011 -- options of cis/trans prolines, and turn off ca bond geometry variation.
/// @author Amelie Stein, amelie.stein@ucsf.edu, Oct 2012 -- vicinity sampling refactoring & new perturbers

// Unit headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

// Project headers
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>

// Rosetta Headers
#include <core/chemical/AA.hh>
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/conversions.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace kinematic_closure {

using namespace std;
using core::Size;
using core::Real;

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_closure.kinematic_closure.KinematicPerturber" );

KinematicPerturber::KinematicPerturber() :
	max_sample_iterations_( basic::options::option[ basic::options::OptionKeys::loops::max_kic_perturber_samples ]() )
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

	KinematicMoverCOP kinmover_op( kinmover() );

	core::Size start( kinmover_op->start_res() );
	core::Size end( kinmover_op->end_res() );
	core::Size first_atom_index ( (pose.residue(start).is_lower_terminus() ?
		kinmover_op->count_bb_atoms_in_residue(pose, start) + 1 :
		kinmover_op->count_bb_atoms_in_residue(pose, start-1) + 1) );

	// AS: if any of the (potentially) omega-sampling flags is set, omega needs to be copied from the torsions vector
	bool also_copy_omega = ( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() || basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() ) || ( basic::options::option[ basic::options::OptionKeys::loops::restrict_kic_sampling_to_torsion_string ].user() || basic::options::option[ basic::options::OptionKeys::loops::derive_torsion_string_from_native_pose ]() ) ;

	for ( core::Size res = start, ia = first_atom_index; res <= end; res++ ) {
		if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) { //If this is a beta-amino acid, set phi, theta, psi, and possibly omega
			using namespace core;
			using namespace core::id;

			for ( core::Size i=1; i<=( res<end ? 3 : 2 /*don't set psi if it's the last residue*/); i++ ) pose.set_torsion( TorsionID( res, id::BB, i ), torsions[ ia++ ] ); //phi, theta, psi
			if ( also_copy_omega && res<end ) pose.set_torsion( TorsionID ( res, id::BB, 4 ), torsions[ia] ); //omega
			ia++; //Increment ia so that it now points to the first atom of the next residue.
		} else { //Default case -- alpha-amino acid
			pose.set_phi( res, torsions[ ia++ ] );
			pose.set_psi( res, torsions[ ia++ ] );
			if ( also_copy_omega && res<end ) pose.set_omega( res, torsions[ ia ] );
			ia++; //Increment atom counter in any case; now at first atom of next residue
		}
	} //looping through residues of segment

	pose.update_residue_neighbors();

} //set_pose_after_closure

KinematicMoverCAP
KinematicPerturber::kinmover() const {
	return kinmover_;
}


///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////TorsionSamplingKinematicPerturber////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

TorsionSamplingKinematicPerturber::TorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in ) :
	KinematicPerturber(),
	vary_ca_bond_angles_(false),
	sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
{
	set_kinmover( kinmover_in );
}

TorsionSamplingKinematicPerturber::~TorsionSamplingKinematicPerturber(){}

/// @details Helper function for TorsionSamplingKinematicPerturber::perturb_beta_residue.  This initializes the list of minima in the Ramachandran cube for beta-3-amino acids.
void TorsionSamplingKinematicPerturber::initialize_betaresidue_minima (
	utility::vector1 < core::Real > &philist, //outputs -- will be cleared by this function
	utility::vector1 < core::Real > &thetalist, //outputs -- will be cleared by this function
	utility::vector1 < core::Real > &psilist, //outputs -- will be cleared by this function
	const core::Size mode //mode 1 initializes for beta-3-amino acids.
) const {
	if ( mode==1 ) {
		philist.clear(); thetalist.clear(); psilist.clear();
		//The following minima are taken from Gunther and Hofmann.  (2001.)  Searching for Periodic Structures in Beta-Peptides.  J. Phys. Chem. B 105:5559-67.
		philist.push_back(-143.7); thetalist.push_back( -63.8); psilist.push_back(-144.6); //B1 from Table 2 in the above reference
		philist.push_back( -71.9); thetalist.push_back( 144.4); psilist.push_back( -80.9); //B2
		philist.push_back(  60.4); thetalist.push_back(-124.1); psilist.push_back(  86.7); //B2'
		philist.push_back( -63.7);  thetalist.push_back( -45.0); psilist.push_back( 111.2); //B3
		philist.push_back(  55.3); thetalist.push_back(  51.1); psilist.push_back(-116.0); //B3'
		philist.push_back(  63.5); thetalist.push_back(  59.3); psilist.push_back(-159.0); //B4
		philist.push_back(-160.6); thetalist.push_back(  56.0); psilist.push_back( 102.6); //B5
		philist.push_back(-111.9); thetalist.push_back(  60.2); psilist.push_back(  25.6); //B6
		philist.push_back( -75.6); thetalist.push_back( 168.6); psilist.push_back( 178.6); //B7
		philist.push_back(  61.7); thetalist.push_back( 163.5); psilist.push_back( 150.2); //B8
		philist.push_back( -91.4); thetalist.push_back(  49.8); psilist.push_back(  90.4); //B9
		philist.push_back(  78.8); thetalist.push_back( -47.0); psilist.push_back( -95.3); //B9'
		philist.push_back(-155.0); thetalist.push_back(  63.4); psilist.push_back(-129.3); //B10
		philist.push_back(-109.4); thetalist.push_back(  74.1); psilist.push_back( -90.8); //B11
		philist.push_back(-159.9); thetalist.push_back( -68.0); psilist.push_back(  29.5); //B12
		//Additional minimia from the same publication corresponding to certain periodic secondary structures (Table 3):
		//Note that these duplicate some of the minima above, from single-residue conformational sampling, but I don't mind biasing the result towards forming secondary structure, a bit.
		philist.push_back(-156.0); thetalist.push_back(  58.0); psilist.push_back(-119.0); //H14
		philist.push_back(-144.0); thetalist.push_back(  87.0); psilist.push_back( -66.0); //H12
		philist.push_back(  64.0); thetalist.push_back(  59.0); psilist.push_back(  75.0); //H10
		philist.push_back( -72.0); thetalist.push_back( 165.0); psilist.push_back( -87.0); //H8I
		philist.push_back(-123.0); thetalist.push_back(  56.0); psilist.push_back(  21.0); //H8II
		philist.push_back(  62.0); thetalist.push_back(  61.0); psilist.push_back( -87.0); //H8III
		philist.push_back(-159.0); thetalist.push_back(  60.0); psilist.push_back(  75.0); //H6II
		philist.push_back( -74.0); thetalist.push_back( 160.0); psilist.push_back( 161.0); //H_E
		//TODO -- add additional conformations consistent with various secondary structures.
	}
	return;
}

/// @details This randomly picks a minimum from the enumerated minimia of the beta-amino acid Ramachandran cube, then chooses phi/theta/psi values randomly in a Gaussian centered in that minimum.  The behaviour depends on beta_residue_type (0=pick everything randomly, 1=pick using beta-3-alanine Ramachandran cube).
void TorsionSamplingKinematicPerturber::perturb_beta_residue (
	core::Real &phi, //output value
	core::Real &theta, //output value
	core::Real &psi, //output value
	const core::Size beta_residue_type //0 = perturb all angles totally randomly; 1 = perturb using Ramachandran cube of beta-3-alanine; other behaviours may be added.
) const {

	const core::Real wellwidth = 12.0; //Standard deviation of the width of the well around each minimum, in degrees.  Arbitrarily chosen.
	utility::vector1 < core::Real > phi_minima;
	utility::vector1 < core::Real > theta_minima;
	utility::vector1 < core::Real > psi_minima;
	initialize_betaresidue_minima(phi_minima, theta_minima, psi_minima, beta_residue_type); //Calling this every time is a bit inefficient...  Still, it's not a lot of data that gets initialized, here.
	core::Size randindex = 0;

	switch (beta_residue_type) {
	case 1 : //Case 1 -- perturb using the minima of the Ramachandran cube of beta-3-alanine.
		randindex = numeric::random::rg().random_range(1, phi_minima.size()); //Pick a random minimum.  Note: all minima are treated as being equally probable.
		phi = phi_minima[randindex] + wellwidth*numeric::random::rg().gaussian();
		theta = theta_minima[randindex] + wellwidth*numeric::random::rg().gaussian();
		psi = psi_minima[randindex] + wellwidth*numeric::random::rg().gaussian();
		break;
	default : //Case zero -- perturb all angles totally randomly
		phi = numeric::random::rg().uniform() * 360.0;
		theta = numeric::random::rg().uniform() * 360.0;
		psi = numeric::random::rg().uniform() * 360.0;
		break;
	}

	if ( phi>360.0 ) phi-=360.0; else if ( phi<360.0 ) phi+=360.0;
	if ( theta>360.0 ) theta-=360.0; else if ( theta<360.0 ) theta+=360.0;
	if ( psi>360.0 ) psi-=360.0; else if ( psi<360.0 ) psi+=360.0;

	return;
}

/// @details randomly varies the torsions (and possibly the bond angles) for the loop.  Respects a MoveMap, if present, for torsions.  Does NOT respect the movemap for angles; does NOT cause any sort of interactions between the MoveMap and the KinematicMover's pivot residues.
void
TorsionSamplingKinematicPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & //bond_len
)
{
	core::kinematics::MoveMapCOP mm(get_movemap());
	KinematicMoverCOP kinmover_op( kinmover() );

	//Get the start and end residues:
	core::Size startres = kinmover_op->start_res();
	core::Size endres = kinmover_op->end_res();

	// Get the number of backbone atoms stored for the padding residues:
	// Number of backbone atoms for the first residue in the segment
	// (start_res_ - 1 or the prepended start_res_ if start_res_ is a terminus)
	core::Size start_minus_one_bb_atom_count = pose.residue(startres).is_lower_terminus() ?
		kinmover_op->count_bb_atoms_in_residue(pose, startres) :
		kinmover_op->count_bb_atoms_in_residue(pose, startres - 1);


	if ( vary_ca_bond_angles_ ) { //For now, ONLY CA bond angles of alpha-amino acids will be varied.

		core::Size pvatom1 = start_minus_one_bb_atom_count + 2; // Second backbone atom of start_res_ (CA if alpha or beta-amino acid).

		core::Real bangle_min( kinmover_op->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_op->BANGLE_SD() );

		//Looping over all CA angles:
		for ( core::Size ir = startres, curatom = pvatom1; ir<=endres; ++ir ) {
			if ( !kinmover_op->is_beta_aminoacid(pose.residue(ir)) /*Add checks here as other backbones are added*/ ) bond_ang[curatom] = bangle_min + numeric::random::rg().uniform() * bangle_sd; //Shouldn't this be bangle_avg() + RG.gaussian() * bangle_sd?
			curatom += kinmover_op->count_bb_atoms_in_residue(pose, ir);
		}
	} //if( vary_ca_bond_angles_ )

	core::Size tor_end = torsions.size() - 3;

	for ( core::Size i=start_minus_one_bb_atom_count + 1, cur_res = startres; i<= tor_end; cur_res++ ) {
		//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;

		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile

			if ( kinmover_op->is_beta_aminoacid(pose.residue(cur_res)) ) { //If this is a beta-amino acid, select backbone dihedral values:
				perturb_beta_residue(torsions[i], torsions[i+1], torsions[i+2], 1);
				i+=4; //Leave omega (CM-C-N-CA) alone and increment i to point at the first torsion of the next residue.
			} else { //Default case -- if this is an alpha-amino acid:
				core::Real rama_phi, rama_psi;
				rama_.random_phipsi_from_rama(pose.aa(cur_res), rama_phi, rama_psi);
				if ( pose.residue(cur_res).has_property( "D_AA" ) ) {
					rama_phi *= -1.0;
					rama_psi *= -1.0;
				}
				torsions[i++]=rama_phi; // phi
				torsions[i++]=rama_psi; // psi
				i++; // leave omega alone
			}

		} else {
			i += kinmover_op->count_bb_atoms_in_residue(pose, cur_res); //ensure i indexing increases
		}

	}

	// sample omegas?
	// all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.
	if (  sample_omega_for_pre_prolines_ ) {
		// old KIC standard: coin flip between cis and trans
		static const core::Real OMEGA_MEAN( 179.8 );

		for ( core::Size i=start_minus_one_bb_atom_count + 1, cur_res = startres; i<= tor_end; cur_res++ ) {

			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro || pose.aa( cur_res+1 ) == core::chemical::aa_dpr ) { //L-proline or D-proline.
				i++; //phi
				if ( kinmover_op->is_beta_aminoacid(pose.residue(cur_res)) ) i++; //theta
				i++; //psi
				torsions[i++] = ( static_cast<int>( numeric::random::rg().uniform()*2 ) ? (core::chemical::is_canonical_D_aa(pose.residue(cur_res).aa()) ? -1.0 : 1.0 ) * OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis).  If it's a D-amino acid, it's multiplied by -1.0 if it's a D-Pro.

			} else {
				i += kinmover_op->count_bb_atoms_in_residue(pose, cur_res); //Increment i by the number of torsion angles in this residue.
			}
		}
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() ) {
		// actual omega sampling, values from Berkholz et al., PNAS 2012
		static const core::Real OMEGA_MEAN( 179.1 );
		static const core::Real OMEGA_STDDEV( 6.3 );

		// cis or trans? -- cis is much less common than trans, according to data from Matt / top8000, so the current coinflip is probably overestimating cis
		// as a proxy, using 1/1000 for now
		static const core::Real cis_prob_threshold = 0.001; //VKM, 26 Aug 2013: this was 0.0001, but the comment above says 1/1000.  I'm switching this to 0.001, since I think that was the intent...

		for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {
			i++; //phi
			if ( kinmover_op->is_beta_aminoacid(pose.residue(cur_res)) ) i++; //theta
			i++; //psi
			core::Real trans_prob = 1;
			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro || pose.aa( cur_res+1 ) == core::chemical::aa_dpr ) { //L- or D-proline
				trans_prob = numeric::random::rg().uniform();
			}
			if ( trans_prob < cis_prob_threshold ) {
				torsions[i++] = 0; // there's very little variation -- currently not captured here at all
			} else { // trans
				torsions[i++] = (core::chemical::is_canonical_D_aa(pose.residue(cur_res).aa()) ? -1.0 : 1.0 ) * (OMEGA_MEAN + numeric::random::rg().gaussian() * OMEGA_STDDEV); //Multiply by -1 if the current residue is D
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
	bool closure_successful
) const
{

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful ); //This correctly handles beta-amino acids, now.
	KinematicMoverCOP kinmover_op( kinmover() );

	core::Size startres = kinmover_op->start_res(); //The first residue in the segment
	core::Size endres = kinmover_op->end_res(); //The last residue in the segment
	//The index of the first backbone atom of the first residue in the segment (i.e. index in the torsions, bond_ang, and bon_len lists):
	core::Size starting_atom = (pose.residue(startres).is_lower_terminus() ? kinmover_op->count_bb_atoms_in_residue(pose, startres) + 1 : kinmover_op->count_bb_atoms_in_residue(pose, startres-1) + 1 ) ;

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );

		// C-N-CA
		for ( Size res=startres, atom=starting_atom; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res);
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res-1);

				core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_C, atomid_N, atomid_CA, offset ); // DOFs canoot be set across jumps (segfault)
				if ( pose.has_dof(dof_of_interest) ) {
					pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
				}
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=startres, atom=starting_atom+1; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res);
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res);
				pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
					numeric::conversions::radians(180 - bond_ang[atom]));
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

		// CA-C-N
		for ( Size res=startres, atom=starting_atom+2; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res+1);
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res);
				core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_CA, atomid_C, atomid_N, offset );
				if ( pose.has_dof(dof_of_interest) ) { //In case this is a terminus
					pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
				}
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=startres, atom=starting_atom; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res);
				const core::id::AtomID atomid_CA(2, res);

				pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ), bond_len[atom]);
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

		// CA-C
		for ( Size res=startres, atom=starting_atom+1; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res);

				pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ), bond_len[atom]);
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

		// C-N
		for ( Size res=startres, atom=starting_atom+2; res<=endres; res++ ) {
			if ( kinmover_op->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_C (3, res);
				const core::id::AtomID atomid_N (1, res+1);

				core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
				if ( pose.has_dof(dof_of_interest) ) {
					pose.set_dof(dof_of_interest, bond_len[atom]);
				}
			}
			atom += kinmover_op->count_bb_atoms_in_residue(pose, res);
		}

	}

	pose.update_residue_neighbors();

} //TorsionSamplingKinematicPerturber::set_pose_after_closure(


///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////VicinitySamplingKinematicPerturber////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

VicinitySamplingKinematicPerturber::VicinitySamplingKinematicPerturber( KinematicMoverCAP kinmover_in )
: KinematicPerturber(),
	vary_ca_bond_angles_(false),
	degree_vicinity_( basic::options::option[ basic::options::OptionKeys::loops::vicinity_degree ]() ),
	sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ) // is this respected at all?
{ set_kinmover( kinmover_in ); }

VicinitySamplingKinematicPerturber::~VicinitySamplingKinematicPerturber(){}

/// @details small variation around the starting phi/psi angles -- order of magnitude is determined by degree_vicinity_
void
VicinitySamplingKinematicPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & //bond_len
)
{
	core::kinematics::MoveMapCOP mm(get_movemap());
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( vary_ca_bond_angles_ ) {

		core::Size pvatom3( (3* (kinmover_op->segment_length() + 1)) - 1 );

		core::Real bangle_min( kinmover_op->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_op->BANGLE_SD() );

		//what is this iterating over?
		for ( Size i = 5; i <= pvatom3; i+=3 ) {
			bond_ang[ i ] = bangle_min + numeric::random::rg().uniform() * bangle_sd;
		}
	}


	core::Size tor_end = torsions.size() - 3;

	for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {
		//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;

		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile

			core::Real rama_phi, rama_psi;

			rama_phi = pose.phi( cur_res ) + degree_vicinity_ * numeric::random::rg().gaussian();
			rama_psi = pose.psi( cur_res ) + degree_vicinity_ * numeric::random::rg().gaussian();

			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi

			i++; // leave omega alone

		} else {
			i += 3; //ensure i indexing increases
		}

	}

	/* [currently] no pre-pro-omega sampling in vicinity mode */

} //perturb_chain


void
VicinitySamplingKinematicPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful // what is this used for?
) const
{

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );


		// C-N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res-1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_C, atomid_N, atomid_CA, offset ); // DOFs canoot be set across jumps (segfault)
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}


		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
				numeric::conversions::radians(180 - bond_ang[atom]));
		}

		// CA-C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res+1);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_CA, atomid_C, atomid_N, offset );
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// CA-C
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_C (3, res);
			const core::id::AtomID atomid_N (1, res+1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_len[atom]));
			}
		}

	}

} //VicinitySamplingKinematicPerturber::set_pose_after_closure(


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
)
{

	if ( sweep_iterator_.at_end() ) {
		utility_exit_with_message("TorsionSweepingKinematicPerturber asked to perturb chain even though sweep iterator is at end.");
	}

	KinematicMoverCOP kinmover_op( kinmover() );
	core::Size start( kinmover_op->start_res() );

	for ( Size ii = 1; ii <= nonpivot_res_to_sweep_.size(); ++ii ) {
		Size torsion_ind = 3 * ( nonpivot_res_to_sweep_[ ii ] - start + 1 ) + sweep_torsion_ids_[ ii ];
		torsions[ torsion_ind ] = sweep_nonpivot_torsion_starts_[ ii ] + sweep_step_sizes_[ ii ] * ( sweep_iterator_[ ii ] - 1 );
		//std::cout << " " <<  nonpivot_res_to_sweep_[ ii ] << " " << sweep_torsion_ids_[ ii ] << " " << torsion_ind << " = " << dt_ang[ torsion_ind ];
	}
	//std::cout << std::endl;
	++sweep_iterator_;
} //TorsionSweepingKinematicPerturber::perturb_chain(


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////// NeighborDependentTorsionSamplingKinematicPerturber //////////////////
///////////////////////////////////////////////////////////////////////////////////////////

NeighborDependentTorsionSamplingKinematicPerturber::NeighborDependentTorsionSamplingKinematicPerturber( KinematicMoverCAP kinmover_in ) :
	KinematicPerturber(),
	vary_ca_bond_angles_(false),
	sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B() )
{
	set_kinmover( kinmover_in );
}

NeighborDependentTorsionSamplingKinematicPerturber::~NeighborDependentTorsionSamplingKinematicPerturber(){}

/// @details randomly varies the torsions (and possibly the bond angles) for the loop, using phi/psi combinations based on rama2b
void
NeighborDependentTorsionSamplingKinematicPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & //bond_len
)
{
	core::kinematics::MoveMapCOP mm(get_movemap());
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( vary_ca_bond_angles_ ) {

		core::Size pvatom3( (3* (kinmover_op->segment_length() + 1)) - 1 );

		core::Real bangle_min( kinmover_op->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_op->BANGLE_SD() );

		for ( Size i = 5; i <= pvatom3; i+=3 ) {
			bond_ang[ i ] = bangle_min + numeric::random::rg().uniform() * bangle_sd;
			//TR << "replacing CA bond angle at " << (kinmover()->start_res()+int((i-4)/3)) << std::endl;
		}
	}


	core::Size tor_end = torsions.size() - 3;

	for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {

		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile

			core::Real rama_phi, rama_psi;

			// warning: not safe for terminal positions -- however, KIC shouldn't actually include any termini anyway...

			// currently we don't really have data for both neighbors together
			// -- for now, do a coin flip on which one to use, though later we should implement this such that each side has an individual perturber, and we call both with equal likelihood (should have fewer ifs)
			static_cast<int>( numeric::random::rg().uniform()*2 ) ? rama_.random_phipsi_from_rama_left(pose.aa(cur_res-1), pose.aa(cur_res),rama_phi, rama_psi) : rama_.random_phipsi_from_rama_right(pose.aa(cur_res), pose.aa(cur_res+1), rama_phi, rama_psi);

			if ( pose.residue(cur_res).has_property( "D_AA" ) ) {
				rama_phi *= -1.0;
				rama_psi *= -1.0;
			}

			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi

			i++; // leave omega alone

		} else {
			i += 3; //ensure i indexing increases
		}

	}

	if (  sample_omega_for_pre_prolines_ ) {
		// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.

		static const core::Real OMEGA_MEAN( 179.8 );

		for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {

			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {

				core::Real rand_omega = ( static_cast<int>( numeric::random::rg().uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)

				i++; //phi
				i++; //psi
				torsions[i++] = rand_omega;

			} else {
				i += 3;
			}
		}

	}


} //perturb_chain


void
NeighborDependentTorsionSamplingKinematicPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful // what is this used for?
) const {

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );


		// C-N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res-1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_C, atomid_N, atomid_CA, offset ); // DOFs canoot be set across jumps (segfault)
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}


		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
				numeric::conversions::radians(180 - bond_ang[atom]));
		}

		// CA-C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res+1);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_CA, atomid_C, atomid_N, offset );
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// CA-C
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_C (3, res);
			const core::id::AtomID atomid_N (1, res+1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_len[atom]));
			}
		}

	}

} //NeighborDependentTorsionSamplingKinematicPerturber::set_pose_after_closure(


///////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////TorsionRestrictedKinematicPerturber////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

TorsionRestrictedKinematicPerturber::TorsionRestrictedKinematicPerturber(
	KinematicMoverCAP kinmover_in,
	core::conformation::torsion_bin_string const & torsion_bins
) :
	KinematicPerturber(),
	vary_ca_bond_angles_(false),
	// if the torsion string is lowercase we use cis, otherwise trans -- note that this flag is used by
	// the parent::set_pose_after_closure function and thus needs to be set to true, otherwise omega
	// will be ignored -- this should be fixed in the future though
	// sample_omega_for_pre_prolines_( true ),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
{
	set_kinmover( kinmover_in );
	predefined_torsions_ = torsion_bins; // store torsion string provided by the user or derived from the native structure
	//std::cerr << " predefined torsion string: " << predefined_torsions_ << std::endl;
}

TorsionRestrictedKinematicPerturber::~TorsionRestrictedKinematicPerturber(){}

/// @details randomly varies the torsions within the given torsion bin
void
TorsionRestrictedKinematicPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & //bond_len
)
{
	core::kinematics::MoveMapCOP mm(get_movemap());
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( vary_ca_bond_angles_ ) { // is this appropriate here?

		core::Size pvatom3( (3* (kinmover_op->segment_length() + 1)) - 1 );

		core::Real bangle_min( kinmover_op->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_op->BANGLE_SD() );

		for ( Size i = 5; i <= pvatom3; i+=3 ) {
			bond_ang[ i ] = bangle_min + numeric::random::rg().uniform() * bangle_sd;
			//TR << "replacing CA bond angle at " << (kinmover()->start_res()+int((i-4)/3)) << std::endl;
		}
	}


	core::Size tor_end = torsions.size() - 3;

	core::Size torsion_string_offset = kinmover_op->start_res() - kinmover_op->loop_begin(); //offset to fetch the torsion bin for the subsegment we're currently sampling -- the string is always for the entire loop -- maybe we should make sure that this doesn't generate a segmentation fault?

	for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {
		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile

			core::Real rama_phi, rama_psi;

			//here we only care about the phi/psi-based torsion bin, omega is set below
			core::conformation::ppo_torsion_bin remapped_torbin = remap_cis_omega_torsion_bins_to_trans( predefined_torsions_[(i-4)/3 + torsion_string_offset] );
			rama_.random_phipsi_from_rama_by_torsion_bin(pose.aa(cur_res), rama_phi, rama_psi, remapped_torbin );
			if ( pose.residue(cur_res).has_property( "D_AA" ) ) {
				rama_phi *= -1.0;
				rama_psi *= -1.0;
			}

			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi

			i++; // leave omega alone

		} else {
			i += 3; //ensure i indexing increases
		}

	}

	// the TorsionRestricted mover will automatically derive cis/trans from the torsion bin string -- sampling for trans omega can be activated via the -loops:kic_omega_sampling flag (but won't affect the cis/trans choice as it does in other perturbers)
	static const core::Real OLD_OMEGA_MEAN( 179.8 );
	// values from Berkholz et al., PNAS 2012
	static const core::Real OMEGA_MEAN( 179.1 );
	static const core::Real OMEGA_STDDEV( 6.3 );
	core::Real rand_omega;

	for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) { // note that this ignores any movemaps for now -- usally they just refer to phi/psi, but it'd be better to check

		if ( basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() || pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
			rand_omega = OLD_OMEGA_MEAN; // mainly for legacy purposes -- at least the means should be adjusted to match
			core::conformation::ppo_torsion_bin i_torbin = predefined_torsions_[(i-4)/3 + torsion_string_offset];
			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro &&
					i_torbin >= core::conformation::ppo_torbin_a &&
					i_torbin != core::conformation::ppo_torbin_X ) {
				rand_omega = 0; // the lowercase torbins (ppo_torbin_{a,b,e,g}) represent cis omega values.

				//std::cerr << "cis proline at " << cur_res << std::endl; // debug

			} else { // trans
				if ( basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() ) {
					rand_omega = OMEGA_MEAN + numeric::random::rg().gaussian() * OMEGA_STDDEV;
				}
			}
			i++; //phi
			i++; //psi
			torsions[i++] = rand_omega;
		} else { // not touching omega for this residue
			i += 3;
		}
	}

} //perturb_chain


void
TorsionRestrictedKinematicPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful // what is this used for?
) const
{

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );


		// C-N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res-1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_C, atomid_N, atomid_CA, offset ); // DOFs canoot be set across jumps (segfault)
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}


		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
				numeric::conversions::radians(180 - bond_ang[atom]));
		}

		// CA-C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res+1);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_CA, atomid_C, atomid_N, offset );
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// CA-C
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_C (3, res);
			const core::id::AtomID atomid_N (1, res+1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_len[atom]));
			}
		}

	}

} //TorsionRestrictedKinematicPerturber::set_pose_after_closure(


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// BaseTabooPerturber /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

// constructor in case we don't know the sequence (e.g. when remodeling multiple loops)
// note that you MUST initialize the sequence before calling apply (else segfaults guaranteed)
BaseTabooPerturber::BaseTabooPerturber( KinematicMoverCAP kinmover_in ) :
	KinematicPerturber(),
	vary_ca_bond_angles_(false),
	sample_omega_for_pre_prolines_( basic::options::option[ basic::options::OptionKeys::loops::sample_omega_at_pre_prolines ]() ),
	// large numbers could be inefficient here if there are lots of design steps (which invalidate
	// the current torsion strings & taboo map) -- however, generation of the random strings seems
	// to be very fast
	num_strings_(1000)
{
	set_kinmover( kinmover_in );
}

BaseTabooPerturber::~BaseTabooPerturber(){}

/// @details Repopulate the torsion-bin-string vector that's used to select the regions of the
/// Ramachandran map from which to draw the phi/psi angles.  This queries the kinmover_
/// to make sure that none of the candidate torsion bin strings have been tried in the past.
/// - for each position in the loop, fetch the corresponding residues, and then ask the Rama
///   lookup table for the relative populations of each torsion bin, and generate a string of
///   length n with the frequencies of each letter corresponding to the respective torsion
///   bin's frequency for this residue
/// - maybe the function to generate a corresponding string can even be part of the Rama function?
/// - fill up with X if the strings aren't long enough -- they must all have the same length
/// - for later applications we could think about adjusting the torsion bin frequencies by those
///   that have already been sampled --> we'd need another structure that holds the strings that
///   were already sampled (might be moved there directly the moment they're provided by this
///   function)
/// - randomly shuffle each string
/// - generate n random torsion strings for the loop by appending positions 1,2,3... of each of the
///   now randomized strings [for the positions] -- this reflects the respective torsion bin
///   frequencies (because the initial strings do, see above)

void
BaseTabooPerturber::refill_torsion_string_vector()
{
	using namespace core::conformation;

	// data structure to hold the vector of torsion bins, with frequencies dependent on the amino
	// acid type, for each position in the loop
	KinematicMoverCOP kinmover_op( kinmover() );
	utility::vector1< torsion_bin_string > torsion_bins_per_position;
	utility::vector1< core::chemical::AA > loop_sequence = kinmover_op->get_loop_sequence();
	torsion_bins_per_position.resize(loop_sequence.size()); // check if we need +1 here somewhere...
	core::Real ideal_freq, current_freq;
	for ( core::Size i = 1; i <= loop_sequence.size(); i++ ) {
		std::map< ppo_torsion_bin, core::Size > entries_per_torsion_bin = get_entries_per_torsion_bin( loop_sequence, i );
		torsion_bin_string torsion_bins_for_pos;
		torsion_bins_for_pos.reserve( num_strings_ );
		// then iterate over the map and fill the torsion_bins_for_pos vector accordingly
		Size entries_for_X = entries_per_torsion_bin[ ppo_torbin_X ];
		for ( std::map< ppo_torsion_bin, core::Size >::const_iterator mcs_i = entries_per_torsion_bin.begin(),
				end = entries_per_torsion_bin.end(); mcs_i != end; ++mcs_i ) {
			// we only keep X in here for the totals
			if ( mcs_i->first == ppo_torbin_X ) continue;

			ideal_freq = core::Real(mcs_i->second) / entries_for_X;
			current_freq = kinmover_op->frequency_in_taboo_map( i-1, mcs_i->first ); // calculated on a string, i.e., base is 0

			//TR << "for residue " << loop_sequence[i] << ", torsion bin " << mcs_i->first << " has " << mcs_i->second << " entries --> " << ideal_freq << " while current freq is " << current_freq << std::endl;

			if ( current_freq > 0 ) {
				// adjust frequency for what has been sampled already -- note that this means that the size
				// of the resulting string isn't necessarily num_strings_ any more
				ideal_freq /= current_freq;
			}
			core::Size sample_rate = static_cast< core::Size > ( ideal_freq * num_strings_ );
			for ( core::Size ii = 1; ii < sample_rate; ++ii ) { torsion_bins_for_pos.push_back( mcs_i->first ); }
		}

		// make sure they all have the same size -- with adjustments the size may vary (slightly)
		// APL: what happens if torsion_bins_for_pos has a size greater than num_strings_?
		while ( torsion_bins_for_pos.size() < num_strings_ ) {
			torsion_bins_for_pos.push_back( ppo_torbin_X );
		}

		numeric::random::random_permutation(torsion_bins_for_pos.begin(), torsion_bins_for_pos.end(), numeric::random::rg());
		torsion_bins_per_position[i] = torsion_bins_for_pos;
	}

	for ( core::Size j = 0; j < num_strings_; j++ ) {

		//std::string new_torsion_string; // = new std::string(sequence_.size()); // initialize to correct size
		torsion_bin_string new_torsion_string;
		for ( core::Size i = 1; i <= loop_sequence.size(); i++ ) {
			// make sure this is correct... though it will happily segfault if not
			new_torsion_string.push_back( torsion_bins_per_position[i][j] );
		}
		// check if we have already sampled this string -- if not, add it to the list
		if ( !( kinmover_op->is_in_taboo_map(new_torsion_string) ) ) {
			random_torsion_strings_.push_back(new_torsion_string);
		}

	}
}

core::conformation::torsion_bin_string
BaseTabooPerturber::next_torsion_string()
{

	/*
	// check if the queue still contains a string -- if so, return it;
	// if not, generate new torsion strings that should cover the entire space (ha!)
	*/

	if ( random_torsion_strings_.size() == 0 ) {
		random_torsion_strings_.clear();
		refill_torsion_string_vector();

		// when the taboo map gets fuller, it is possible that all 1000 random strings have already been tried,
		// and thus the map is still empty even after "refilling."
		// APL: This has all the hallmarks of an infinite loop.
		while ( random_torsion_strings_.size() == 0 ) {
			refill_torsion_string_vector();
		}
	}


	runtime_assert(random_torsion_strings_.size() > 0); // shouldn't happen any more...

	KinematicMoverCOP kinmover_op( kinmover() );

	// make sure the torsion string we're returning hasn't been tested yet
	while ( kinmover_op->is_in_taboo_map(random_torsion_strings_[random_torsion_strings_.size()]) ) {
		//TR << random_torsion_strings_[random_torsion_strings_.size()] << " has already been tested, next... " << std::endl;
		random_torsion_strings_.pop_back();
		while ( random_torsion_strings_.size() == 0 ) { // refill if necessary
			refill_torsion_string_vector();
		}
	}

	core::conformation::torsion_bin_string tb = random_torsion_strings_[ random_torsion_strings_.size() ];
	//TR << tb << " " << random_torsion_strings_.size() <<  std::endl; // debug

	random_torsion_strings_.pop_back();

	return tb;
}

void
BaseTabooPerturber::perturb_chain(
	core::pose::Pose const & pose,
	utility::vector1< core::Real > & torsions,
	utility::vector1< core::Real > & bond_ang,
	utility::vector1< core::Real > & //bond_len
)
{
	core::kinematics::MoveMapCOP mm(get_movemap());
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( vary_ca_bond_angles_ ) { // is this appropriate here?

		core::Size pvatom3( (3* (kinmover_op->segment_length() + 1)) - 1 );

		core::Real bangle_min( kinmover_op->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_op->BANGLE_SD() );

		//what is this iterating over?
		for ( Size i = 5; i <= pvatom3; i+=3 ) {
			bond_ang[ i ] = bangle_min + numeric::random::rg().uniform() * bangle_sd;
			//TR << "replacing CA bond angle at " << (kinmover()->start_res()+int((i-4)/3)) << std::endl;
		}
	}


	core::Size tor_end = torsions.size() - 3;

	// function to provide the torsion string to be sampled -- to be written
	core::conformation::torsion_bin_string torsion_string = next_torsion_string();

	// offset to fetch the torsion bin for the subsegment we're currently sampling -- the string
	// is always for the entire loop -- maybe we should make sure that this doesn't generate a
	// segmentation fault?
	core::Size torsion_string_offset = kinmover_op->start_res() - kinmover_op->loop_begin();

	//std::cerr << "Torsion string is " << torsion_string << ", offset is " << torsion_string_offset << std::endl;

	for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {
		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile

			core::Real rama_phi, rama_psi;

			//rama_.random_phipsi_from_rama_by_torsion_bin(pose.aa(cur_res), rama_phi, rama_psi, toupper(torsion_string[(i-4)/3 + torsion_string_offset])); //here we only care about the phi/psi-based torsion bin, omega is set below
			get_random_phi_psi_for_residue( pose, cur_res, torsion_string[ (i-4)/3 + torsion_string_offset ], rama_phi, rama_psi );
			// This function is given the AA. Therefore, it needs no correction.

			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi

			i++; // leave omega alone

		} else {
			i += 3; //ensure i indexing increases
		}

	}

	if (  sample_omega_for_pre_prolines_ ) {
		// sample omegas. all omegas before prolines are assumed to be fair game. Note that we're not using move-map, which is phi/psi-specific.

		static const core::Real OMEGA_MEAN( 179.8 );

		for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {

			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {

				core::Real rand_omega = ( static_cast<int>( numeric::random::rg().uniform()*2 ) ? OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis)

				/*
				if ( islower(torsion_string[(i-4)/3 + torsion_string_offset]) ) { // actually at the moment this is impossible... use random value instead?
				rand_omega = 0; // lowercase is for cis
				} else {
				rand_omega = OMEGA_MEAN;
				}
				//}
				*/
				i++; //phi
				i++; //psi
				torsions[i++] = rand_omega;

			} else {
				i += 3;
			}
		}

	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() ) {
		// actual omega sampling, values from Berkholz et al., PNAS 2012
		static const core::Real OMEGA_MEAN( 179.1 );
		static const core::Real OMEGA_STDDEV( 6.3 );

		// cis or trans? -- cis is much less common than trans, according to data from Matt / top8000, so the current coinflip is probably overestimating cis
		// as a proxy, using 1/1000 for now
		static const core::Real cis_prob_threshold = 0.0001;

		for ( core::Size i=4, cur_res = kinmover_op->start_res(); i<= tor_end; cur_res++ ) {
			i++; //phi
			i++; //psi
			core::Real trans_prob = 1;
			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro ) {
				trans_prob = numeric::random::rg().uniform();
			}
			if ( trans_prob < cis_prob_threshold ) {
				torsions[i++] = 0; // there's very little variation -- currently not captured here at all
			} else { // trans
				torsions[i++] = OMEGA_MEAN + numeric::random::rg().gaussian() * OMEGA_STDDEV;
			}
		}
	}


} //perturb_chain

void
BaseTabooPerturber::set_pose_after_closure(
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful // what is this used for?
) const
{

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful );
	KinematicMoverCOP kinmover_op( kinmover() );

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );


		// C-N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res-1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_C, atomid_N, atomid_CA, offset ); // DOFs canoot be set across jumps (segfault)
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}


		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
				numeric::conversions::radians(180 - bond_ang[atom]));
		}

		// CA-C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res+1);
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);
			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_angle_dof_id(atomid_CA, atomid_C, atomid_N, offset );
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_ang[atom]));
			}
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=kinmover_op->start_res(), atom=4; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_N (1, res);
			const core::id::AtomID atomid_CA(2, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// CA-C
		for ( Size res=kinmover_op->start_res(), atom=5; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_CA(2, res);
			const core::id::AtomID atomid_C (3, res);

			pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ),
				numeric::conversions::radians(180 - bond_len[atom]));
		}

		// C-N
		for ( Size res=kinmover_op->start_res(), atom=6; res<= kinmover_op->end_res(); res++, atom+=3 ) {
			const core::id::AtomID atomid_C (3, res);
			const core::id::AtomID atomid_N (1, res+1);

			core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
			if ( pose.has_dof(dof_of_interest) ) {
				pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_len[atom]));
			}
		}

	}

} //TabooSamplingKinematicPerturber::set_pose_after_closure(


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// TabooSamplingKinematicPerturber ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

TabooSamplingKinematicPerturber::TabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in ) :
	BaseTabooPerturber( kinmover_in ),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() )
{}

TabooSamplingKinematicPerturber::~TabooSamplingKinematicPerturber(){}


void
TabooSamplingKinematicPerturber::get_random_phi_psi_for_residue(
	core::pose::Pose const & pose,
	core::Size resid,
	core::conformation::ppo_torsion_bin torbin,
	core::Real & phi,
	core::Real & psi
) const
{
	core::conformation::ppo_torsion_bin remapped_torbin = core::conformation::remap_cis_omega_torsion_bins_to_trans( torbin );
	rama_.random_phipsi_from_rama_by_torsion_bin( pose.aa(resid), phi, psi, remapped_torbin );
	if ( pose.residue( resid ).has_property( "D_AA" ) ) {
		phi *= -1.0;
		psi *= -1.0;
	}

}


std::map< core::conformation::ppo_torsion_bin, core::Size >
TabooSamplingKinematicPerturber::get_entries_per_torsion_bin(
	utility::vector1< core::chemical::AA > loop_seq,
	core::Size resid
) const
{
	std::map< core::conformation::ppo_torsion_bin, core::Size > entries_per_torsion_bin;
	rama_.get_entries_per_torsion_bin( loop_seq[ resid ], entries_per_torsion_bin );
	return entries_per_torsion_bin;
}

///////////////////////////////////////////////////////////////////////////////////////////
//////////////////// NeighborDependentTabooSamplingKinematicPerturber /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

NeighborDependentTabooSamplingKinematicPerturber::NeighborDependentTabooSamplingKinematicPerturber( KinematicMoverCAP kinmover_in ) :
	BaseTabooPerturber( kinmover_in ),
	rama_( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B() )
{}

NeighborDependentTabooSamplingKinematicPerturber::~NeighborDependentTabooSamplingKinematicPerturber(){}

void
NeighborDependentTabooSamplingKinematicPerturber::get_random_phi_psi_for_residue(
	core::pose::Pose const & pose,
	core::Size resid,
	core::conformation::ppo_torsion_bin torbin,
	core::Real & phi,
	core::Real & psi
) const
{
	core::conformation::ppo_torsion_bin remapped_torbin = core::conformation::remap_cis_omega_torsion_bins_to_trans( torbin );
	KinematicMoverCOP kinmover_op( kinmover() );
	if ( resid == kinmover_op->loop_begin() ) {
		rama_.random_phipsi_from_rama_by_torsion_bin_right( pose.aa(resid), pose.aa(resid+1), phi, psi, remapped_torbin );
	} else if ( resid == kinmover_op->loop_end() ) {
		rama_.random_phipsi_from_rama_by_torsion_bin_left( pose.aa(resid-1), pose.aa(resid), phi, psi, remapped_torbin );
	} else if ( static_cast<int>( numeric::random::rg().uniform()*2 ) ) {
		rama_.random_phipsi_from_rama_by_torsion_bin_right( pose.aa(resid), pose.aa(resid+1), phi, psi, remapped_torbin );
	} else {
		rama_.random_phipsi_from_rama_by_torsion_bin_left( pose.aa(resid-1), pose.aa(resid), phi, psi, remapped_torbin );
	}
	if ( pose.residue( resid ).has_property( "D_AA" ) ) {
		phi *= -1.0;
		psi *= -1.0;
	}

}


std::map< core::conformation::ppo_torsion_bin, core::Size >
NeighborDependentTabooSamplingKinematicPerturber::get_entries_per_torsion_bin(
	utility::vector1< core::chemical::AA > loop_seq,
	core::Size resid
) const
{
	std::map< core::conformation::ppo_torsion_bin, core::Size > entries_per_torsion_bin, entries_left, entries_right;
	if ( resid == 1 ) {
		rama_.get_entries_per_torsion_bin_right( loop_seq[ resid ], loop_seq[ resid+1 ], entries_per_torsion_bin );
	} else if ( resid == loop_seq.size() ) {
		rama_.get_entries_per_torsion_bin_left( loop_seq[ resid-1 ], loop_seq[ resid ], entries_per_torsion_bin );
	} else {
		// take the minimum count from both sides -- especially if one of them is 0
		// we cannot ask for a random phi/psi combination as it won't exist
		rama_.get_entries_per_torsion_bin_left(  loop_seq[ resid-1 ], loop_seq[ resid ],   entries_left );
		rama_.get_entries_per_torsion_bin_right( loop_seq[ resid ],   loop_seq[ resid+1 ], entries_right );
		for ( std::map< core::conformation::ppo_torsion_bin, core::Size >::const_iterator
				iter = entries_left.begin();
				iter != entries_left.end();
				++iter ) {
			entries_per_torsion_bin[iter->first] = std::min( iter->second, entries_right[ iter->first ] );
		}
	}
	return entries_per_torsion_bin;
}

} // namespace kinematic_closure
} // namespace loop_closure
} // namespace loops
} // namespace protocols
