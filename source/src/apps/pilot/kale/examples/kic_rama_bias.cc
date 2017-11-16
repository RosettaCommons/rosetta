// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicPerturber.hh>

// Utility headers
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace core;
using namespace basic::options;
using core::pose::Pose;
using core::pose::make_pose_from_sequence;
using core::scoring::ScoreFunction;
using core::scoring::ScoringManager;
using protocols::moves::MonteCarlo;
using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
using protocols::loops::loop_closure::kinematic_closure::KinematicPerturber;
using numeric::random::uniform;

OPT_2GRP_KEY(File, kale, example, prefix)
OPT_2GRP_KEY(Integer, kale, example, iterations)

class UniformPerturber : public KinematicPerturber { // {{{1

public:
	typedef KinematicPerturber parent;

	UniformPerturber();

	string perturber_type() const { return "UniformPerturber"; }

	void perturb_chain(
		core::pose::Pose const & pose,
		utility::vector1< core::Real> & torsions,
		utility::vector1< core::Real> & bond_ang,
		utility::vector1< core::Real> & bond_len);

	void set_pose_after_closure(
		core::pose::Pose & pose,
		utility::vector1< core::Real> const & torsions,
		utility::vector1< core::Real> const & bond_ang,
		utility::vector1< core::Real> const & bond_len,
		bool closure_successful) const;

	void set_vary_ca_bond_angles( bool vary_ca_bond_angles ) {
		vary_ca_bond_angles_ = vary_ca_bond_angles; }

private:
	bool vary_ca_bond_angles_;
	bool sample_omega_for_pre_prolines_;
	core::scoring::Ramachandran const & rama_;

};

static basic::Tracer TR(
	"protocols.loops.loop_closure.kinematic_closure.KinematicPerturber");

UniformPerturber::UniformPerturber() // {{{1
: KinematicPerturber(),
	vary_ca_bond_angles_(false),
	sample_omega_for_pre_prolines_(false),
	rama_(ScoringManager::get_instance()->get_Ramachandran()) {}

void UniformPerturber::perturb_chain( // {{{1
	core::pose::Pose const & pose,
	utility::vector1<core::Real> & torsions,
	utility::vector1<core::Real> & bond_ang,
	utility::vector1<core::Real> & bond_len) {

	core::kinematics::MoveMapCOP mm(get_movemap());

	//Get the start and end residues:
	core::Size startres = kinmover_->start_res();
	core::Size endres = kinmover_->end_res();

	//Get the number of backbone atoms stored for the padding residues:
	core::Size start_minus_one_bb_atom_count = pose.residue(startres).is_lower_terminus() ? kinmover_->count_bb_atoms_in_residue(pose, startres) : kinmover_->count_bb_atoms_in_residue(pose, startres - 1); //Number of backbone atoms for the first residue in the segment (start_res_ - 1 or the prepended start_res_ if start_res_ is a terminus)
	core::Size end_plus_one_bb_atom_count = pose.residue(endres).is_upper_terminus() ? kinmover_->count_bb_atoms_in_residue(pose, endres) : kinmover_->count_bb_atoms_in_residue(pose, endres); //Number of backbone atoms for the last residue in the segment (end_res_ + 1 or the appended end_res_ if end_res_ is a terminus)


	if ( vary_ca_bond_angles_ ) { //For now, ONLY CA bond angles of alpha-amino acids will be varied.

		core::Size pvatom1 = start_minus_one_bb_atom_count + 2; // Second backbone atom of start_res_ (CA if alpha or beta-amino acid).
		core::Size pvatom3 = bond_ang.size() - end_plus_one_bb_atom_count - kinmover_->count_bb_atoms_in_residue(pose, endres) + 2; // Second backbone atom of end_res_ (CA if alpha or beta-amino acid).

		core::Real bangle_min( kinmover_->BANGLE_MIN() );
		core::Real bangle_sd( kinmover_->BANGLE_SD() );

		//Looping over all CA angles:
		for ( core::Size ir = startres, curatom = pvatom1; ir<=endres; ++ir ) {
			if ( !kinmover_->is_beta_aminoacid(pose.residue(ir)) /*Add checks here as other backbones are added*/ ) bond_ang[curatom] = bangle_min + numeric::random::rg().uniform() * bangle_sd; //Shouldn't this be bangle_avg() + RG.gaussian() * bangle_sd?
			curatom += kinmover_->count_bb_atoms_in_residue(pose, ir);
		}
	} //if( vary_ca_bond_angles_ )

	core::Size tor_end = torsions.size() - 3;

	for ( core::Size i=start_minus_one_bb_atom_count + 1, cur_res = startres; i<= tor_end; cur_res++ ) {
		//if(mm) TR << "current residue " << cur_res << "mm reports " << mm->get_bb(cur_res) << std::endl;

		if ( !mm || mm->get_bb(cur_res) ) { //if movemap is null, or (if not null) if movemap says mobile
			core::Real rama_phi, rama_psi;
			//rama_.uniform_phipsi_from_allowed_rama(pose.aa(cur_res), rama_phi,
			//rama_psi);
			rama_phi = 360 * numeric::random::uniform();
			rama_psi = 360 * numeric::random::uniform();
			torsions[i++]=rama_phi; // phi
			torsions[i++]=rama_psi; // psi
			i++; // leave omega alone
		} else {
			i += kinmover_->count_bb_atoms_in_residue(pose, cur_res); //ensure i indexing increases
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
				if ( kinmover_->is_beta_aminoacid(pose.residue(cur_res)) ) i++; //theta
				i++; //psi
				torsions[i++] = ( static_cast<int>( numeric::random::rg().uniform()*2 ) ? (core::chemical::is_D_aa(pose.residue(cur_res).aa()) ? -1.0 : 1.0 ) * OMEGA_MEAN : 0.0 );  // flip a coin -- either 179.8 (trans) or 0.0 (cis).  If it's a D-amino acid, it's multiplied by -1.0 if it's a D-Pro.

			} else {
				i += kinmover_->count_bb_atoms_in_residue(pose, cur_res); //Increment i by the number of torsion angles in this residue.
			}
		}
	} else if ( basic::options::option[ basic::options::OptionKeys::loops::kic_omega_sampling ]() ) {
		// actual omega sampling, values from Berkholz et al., PNAS 2012
		static const core::Real OMEGA_MEAN( 179.1 );
		static const core::Real OMEGA_STDDEV( 6.3 );

		// cis or trans? -- cis is much less common than trans, according to data from Matt / top8000, so the current coinflip is probably overestimating cis
		// as a proxy, using 1/1000 for now
		static const core::Real cis_prob_threshold = 0.001; //VKM, 26 Aug 2013: this was 0.0001, but the comment above says 1/1000.  I'm switching this to 0.001, since I think that was the intent...

		for ( core::Size i=4, cur_res = kinmover_->start_res(); i<= tor_end; cur_res++ ) {
			i++; //phi
			if ( kinmover_->is_beta_aminoacid(pose.residue(cur_res)) ) i++; //theta
			i++; //psi
			core::Real trans_prob = 1;
			if ( pose.aa( cur_res+1 ) == core::chemical::aa_pro || pose.aa( cur_res+1 ) == core::chemical::aa_dpr ) { //L- or D-proline
				trans_prob = numeric::random::rg().uniform();
			}
			if ( trans_prob < cis_prob_threshold ) {
				torsions[i++] = 0; // there's very little variation -- currently not captured here at all
			} else { // trans
				torsions[i++] = (core::chemical::is_D_aa(pose.residue(cur_res).aa()) ? -1.0 : 1.0 ) * (OMEGA_MEAN + numeric::random::rg().gaussian() * OMEGA_STDDEV); //Multiply by -1 if the current residue is D
			}
		}
	}
}

void UniformPerturber::set_pose_after_closure( // {{{1
	core::pose::Pose & pose,
	utility::vector1<core::Real> const & torsions,
	utility::vector1<core::Real> const & bond_ang,
	utility::vector1<core::Real> const & bond_len,
	bool closure_successful ) const {

	// parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful, sample_omega_for_pre_prolines_ );
	parent::set_pose_after_closure( pose, torsions, bond_ang, bond_len, closure_successful ); //This correctly handles beta-amino acids, now.

	core::Size startres = kinmover_->start_res(); //The first residue in the segment
	core::Size endres = kinmover_->end_res(); //The last residue in the segment
	//The index of the first backbone atom of the first residue in the segment (i.e. index in the torsions, bond_ang, and bon_len lists):
	core::Size starting_atom = (pose.residue(startres).is_lower_terminus() ? kinmover_->count_bb_atoms_in_residue(pose, startres) + 1 : kinmover_->count_bb_atoms_in_residue(pose, startres-1) + 1 ) ;

	if ( !closure_successful || vary_ca_bond_angles_ ) { // if the closure wasn't successful, we may need to overwrite previously idealized angles
		core::Real offset( 0.0 );

		// C-N-CA
		for ( Size res=startres, atom=starting_atom; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
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
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

		// N-CA-C -- these are all within the same residue, so jumps are not an issue
		for ( Size res=startres, atom=starting_atom+1; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res);
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res);
				pose.set_dof(pose.atom_tree().bond_angle_dof_id(atomid_N, atomid_CA, atomid_C, offset ),
					numeric::conversions::radians(180 - bond_ang[atom]));
			}
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

		// CA-C-N
		for ( Size res=startres, atom=starting_atom+2; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
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
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

	}

	// overwrite bond lengths, at least if the closure was not successful
	if ( !closure_successful ) { // if sampling of bond lengths is added, activate this section

		// N-CA
		for ( Size res=startres, atom=starting_atom; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_N (1, res);
				const core::id::AtomID atomid_CA(2, res);

				pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_N, atomid_CA ),
					numeric::conversions::radians(180 - bond_len[atom]));
			}
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

		// CA-C
		for ( Size res=startres, atom=starting_atom+1; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_CA(2, res);
				const core::id::AtomID atomid_C (3, res);

				pose.set_dof(pose.atom_tree().bond_length_dof_id(atomid_CA, atomid_C ),
					numeric::conversions::radians(180 - bond_len[atom]));
			}
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

		// C-N
		for ( Size res=startres, atom=starting_atom+2; res<=endres; res++ ) {
			if ( kinmover_->is_beta_aminoacid(pose.residue(res)) ) {
				//For now, do nothing for beta-amino acids.
			} else { //Default case -- alpha-amino acid:
				const core::id::AtomID atomid_C (3, res);
				const core::id::AtomID atomid_N (1, res+1);

				core::id::DOF_ID dof_of_interest = pose.atom_tree().bond_length_dof_id(atomid_C, atomid_N);
				if ( pose.has_dof(dof_of_interest) ) {
					pose.set_dof(dof_of_interest, numeric::conversions::radians(180 - bond_len[atom]));
				}
			}
			atom += kinmover_->count_bb_atoms_in_residue(pose, res);
		}

	}
	pose.update_residue_neighbors();
}

// }}}1

int main(int argc, char** argv) {
	NEW_OPT(kale::example::prefix, "Output prefix", "");
	NEW_OPT(kale::example::iterations, "Iterations", 720);

	devel::init(argc, argv);

	Size iterations = option[OptionKeys::kale::example::iterations]();

	string rama_path =
		option[OptionKeys::kale::example::prefix]().name() + "rama.dat";
	string uniform_path =
		option[OptionKeys::kale::example::prefix]().name() + "uniform.dat";

	ofstream rama_log; rama_log.open(rama_path.c_str());
	ofstream uniform_log; uniform_log.open(uniform_path.c_str());

	Pose rama_pose, uniform_pose;
	make_pose_from_sequence(rama_pose, "GGGGGG", core::chemical::FA_STANDARD, false);
	make_pose_from_sequence(uniform_pose, "GGGGGG", core::chemical::FA_STANDARD, false);

	ScoreFunction score_function;
	score_function.set_weight(core::scoring::rama, 1);

	MonteCarlo rama_monte_carlo(rama_pose, score_function, 1);
	MonteCarlo uniform_monte_carlo(uniform_pose, score_function, 1);

	KinematicMover rama_mover;
	rama_mover.set_pivots(2, 3, 5);

	KinematicMover uniform_mover;
	uniform_mover.set_pivots(2, 3, 5);
	uniform_mover.set_perturber(new UniformPerturber);

	for ( Size i = 1; i <= iterations; i++ ) {
		cerr << "\r[" << i << "/" << iterations << "]";
		Real phi = 360 * uniform();

		rama_pose.set_phi(3, phi);
		rama_mover.apply(rama_pose);
		rama_monte_carlo.boltzmann(rama_pose);
		rama_log
			<< rama_pose.phi(4) << " "
			<< rama_pose.psi(4) << endl;

		uniform_pose.set_phi(3, phi);
		uniform_mover.apply(uniform_pose);
		uniform_monte_carlo.boltzmann(uniform_pose);
		uniform_log
			<< uniform_pose.phi(4) << " "
			<< uniform_pose.psi(4) << endl;
	}
	cerr << endl;
}
