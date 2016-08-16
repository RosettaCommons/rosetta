// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file cyclic_sampler.cc
/// @brief run Monte Carlo for sampling protein conformation
/// @author Yuan Liu
/// @details
/// Modified backbone Gaussian sampler for cyclic peptide
/// 0. currently only work for full-atom
/// 1. randomly select N(>=4) hinge residues
/// 2. restrain the CYS-CYS distance (tune factor A/B)
/// 3. repack CYS CYS
/// 4. simulated annealing or cartmin

// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/Minimizer.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>

#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>

//jd2
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/InnerJob.hh>

void *my_main( void *);

static numeric::random::RandomGenerator RG(62331911);
basic::Tracer TR("pilot.wendao.cyclic");

using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY(Boolean, save_trajectory)
OPT_KEY(Integer, traj_interval)
OPT_KEY(Boolean, min_after_mc)
OPT_KEY(Integer, mc_steps)
OPT_KEY(Real, mc_temperature)
OPT_KEY(Real, mc_temp_annealing)
OPT_KEY(Integer, mc_temp_step)

int main( int argc, char *argv [] )
{
	NEW_OPT(save_trajectory, "save trajectory", false);
	NEW_OPT(traj_interval, "trajectory interval", 10000);
	NEW_OPT(min_after_mc, "minimize", false);
	NEW_OPT(mc_steps, "steps", 100000);
	NEW_OPT(mc_temperature, "temperature", 1.0);
	NEW_OPT(mc_temp_annealing, "annealing T", 1.0);
	NEW_OPT(mc_temp_step, "annealing step", 0);

	try {
		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );
	}
	catch ( utility::excn::EXCN_Base const &e) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

////////////////////////////////////////////////////
// BBGMover for cyclic peptide backbone
using namespace protocols::moves;
using namespace protocols::simple_moves;

class BBG_Cyclic_Mover : public BBGaussianMover
{
public:
	BBG_Cyclic_Mover():
		BBGaussianMover(1,28,4),
		dphi(utility::vector1<Real>(n_dof_angle_))
	{
		protocols::moves::Mover::type("BBG_Cyclic_Mover");

		factorA_ = option[bbg::factorA];
		factorB_ = option[bbg::factorB];

		// core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
		// main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
		// pack_mover_ = new protocols::simple_moves::PackRotamersMover;
		// pack_mover_->task_factory( main_task_factory );
		// pack_mover_->score_function( core::scoring::getScoreFunction() );
	}

	void apply( Pose &pose )
	{
		//setup locker
		if (locker_.size()==0) {
			setup_locker(pose);
			assert(locker_.size()==2);
			TR << "CYD: " << locker_[1] << "<-->" << locker_[2] << std::endl;
		}

		//select pivot residues
		// Real rA, rB, rC, rD;
		// 	rA = 2;
		// 	rB = 5;
		// 	rC = 8;
		// 	rD = 11;
		// std::cout << rA << "," << rB << "," << rC << "," << rD << std::endl;

		//bbg move
		//get_VdRdPhi(pose, rA, rB, rC, rD);
		get_VdRdPhi(pose);
		get_G();
		get_A();
		//Real W_old = get_L_move(pose, rA, rB, rC, rD);
		Real W_old = get_L_move(pose);

		//get_VdRdPhi(pose, rA, rB, rC, rD);
		get_VdRdPhi(pose);
		get_G();
		get_A();
		Real W_new = get_L_prime();

		// -- repack CYS and the end?

		//random tail move

		//pack_mover_->apply(pose);
		last_proposal_density_ratio_ = W_new / W_old;
	}

	Real last_proposal_density_ratio() {
		return last_proposal_density_ratio_;
	}

	void get_VdRdPhi(Pose const &pose)
	{
		//take all residues between CYD
		Size first = locker_[1];
		Size last = locker_[2];
		//the last CYD
		conformation::Residue const & lock( pose.residue(last) );
		Vector end_xyz = lock.atom("SG").xyz();

		Size ndof = 0;
		for (Size i=first; i<=last; i++) {
			conformation::Residue const & rsd( pose.residue( i ) );
			if (i>first) {
				//N-CA
				ndof++;
				matrix_dRdPhi[1][ndof] = get_dRdPhi(rsd.atom("N").xyz(), rsd.atom("CA").xyz(), end_xyz);
			}
			if (i<last) {
				//CA-C
				ndof++;
				matrix_dRdPhi[1][ndof] = get_dRdPhi(rsd.atom("CA").xyz(), rsd.atom("C").xyz(), end_xyz);
			}
		}
		//TR << "DOF: " << ndof << std::endl;
	}

	void get_G()
	{
		for (Size i=1; i<=n_dof_angle_; i++) {
			for (Size j=i; j<=n_dof_angle_; j++) {
				matrix_G[i][j] = matrix_dRdPhi[1][i].dot(matrix_dRdPhi[1][j]);
				if (i<j) matrix_G[j][i] = matrix_G[i][j];
			}
		}
	}

	void get_A()
	{
		for (Size i=1; i<=n_dof_angle_; i++) {
			for (Size j=i; j<=n_dof_angle_; j++) {
				matrix_A[i][j] = factorB_ * matrix_G[i][j];
				if (i==j) matrix_A[i][j] += 1.0;
				matrix_A[i][j] *= factorA_ / 2.0;
				if (i<j) matrix_A[j][i] = matrix_A[i][j];
			}
		}
	}

	Real get_L_move(Pose &pose)
	{
		//gerate a Gaussian dx vector
		utility::vector1<Real> delta(n_dof_angle_);
		for (Size i=1; i<=n_dof_angle_; i++) delta[i]=RG.gaussian();
		//calculate d^2 = delta^2
		Real d2=0.0;
		for (Size i=1; i<=n_dof_angle_; i++) d2+=delta[i]*delta[i];
		//cholesky, get L^t, L^-1
		Real detL = cholesky_fw(matrix_A, n_dof_angle_, delta, dphi);
		
		//W_old *= exp(-d^2)
		Real W_old = detL*exp(-d2/2.0);
		//set the new phi,psi (above all called phi, actually 4 phi, 4 psi)
		Size first = locker_[1];
		Size last = locker_[2];
		Size ndof = 0;
		for (Size i=first; i<=last; i++) {
			conformation::Residue const & rsd( pose.residue( i ) );
			if (i>first) {
				//N-CA
				ndof++;
				pose.set_phi(i, basic::periodic_range( pose.phi(i)+dphi[ndof], 360.0 ) );
			}
			if (i<last) {
				//CA-C
				ndof++;
				pose.set_psi(i, basic::periodic_range( pose.psi(i)+dphi[ndof], 360.0 ) );
			}
		}
		return W_old;
	}

	Real get_L_prime()
	{
		//gerate a Gaussian dx vector
		utility::vector1<Real> delta(n_dof_angle_);
		//get L
		Real detL = cholesky_bw(matrix_A, n_dof_angle_, dphi, delta);
		//delta = L^t * dphi
		//calculate d^2 = delta^2
		Real d2=0.0;
		for (Size i=1; i<=n_dof_angle_; i++) d2 += delta[i]*delta[i];
		Real W_new = detL * exp(-d2/2.0);
		return W_new;
	}

	void setup_locker( Size begin, Size end )
	{
		locker_.erase(locker_.begin(),locker_.end());
		locker_.push_back(begin);
		locker_.push_back(end);
	}

	void setup_locker( Pose &pose )
	{
		bool fix_fail = false;
		while (true) {
			locker_.erase(locker_.begin(),locker_.end());

			//go through all residues, lock the first and the last CYD
			for (Size i=1; i<=pose.total_residue(); i++) {
				//TR << pose.residue(i).name().substr(0,3) << std::endl;
				if (pose.residue(i).name3() == "CYS" && pose.residue(i).is_disulfide_bonded() ) {
					locker_.push_back(i);
					TR << "Add CYD: " << i << std::endl;
				}
			}

			if (locker_.size()==2) {
				break;
			}
			
			if (fix_fail) utility_exit_with_message("Can not identify CYD-CYD correctly!");
			fix_disulf(pose);
			fix_fail = true; //just try once
		}
	}

	void fix_disulf( Pose &pose )
	{
		if (option[ in::fix_disulf ].user()) {
			core::pose::initialize_disulfide_bonds(pose);
		}
		else {
			utility_exit_with_message("Can not identify CYD-CYD correctly!");
		}
	}

	virtual std::string get_name() const {
		return "BBG_Cyclic_Mover";
	}

private:
	utility::vector1< Size > locker_;
	utility::vector1< Real > dphi;
	utility::vector1< Size > pivots_;
	Real factorA_;
	Real factorB_;
	//protocols::simple_moves::PackRotamersMoverOP pack_mover_;
};

class MyProtocol : public Mover
{
public:
	MyProtocol()
	{
		protocols::moves::Mover::type("MyProtocol");
	}

	virtual std::string get_name() const {
		return "MyProtocol";
	}

	void apply( Pose &pose )
	{
		using namespace core::pack::task;
		using namespace protocols::moves;
		using namespace protocols::simple_moves;

		core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

		//task
		
		TaskFactoryOP main_task_factory = new TaskFactory;
		main_task_factory->push_back( new operation::InitializeFromCommandline );
		operation::RestrictToRepackingOP rtrop = new operation::RestrictToRepacking;
		main_task_factory->push_back( rtrop );

		//sc
		sidechain_moves::SidechainMoverOP sidechainmover = new sidechain_moves::SidechainMover();
		sidechainmover->set_task_factory(main_task_factory);
		//mark: if don't set, seg fault, why?
		sidechainmover->set_prob_uniform(0.1);
		sidechainmover->set_prob_withinrot(0.6);
		sidechainmover->set_prob_random_pert_current(0);
		sidechainmover->set_preserve_detailed_balance(true);

		//Random
		RandomMoverOP randmover = new RandomMover();
		//add bb
		randmover->add_mover(new BBG_Cyclic_Mover(), 1.0);
		randmover->add_mover(sidechainmover, 4.0);

		//mc 2->0.2
		Size nsteps = option[mc_steps]();
		Real temp_start = option[mc_temperature]();
		Real beta = 1.0/temp_start;
		Real temp_end = temp_start;
		Real deltaB = 0;
		Size sa_step = option[mc_temp_step]();
		Size sa_intv = nsteps / (sa_step+1);
		if (option[mc_temp_annealing].user() && sa_step>0) {
			temp_end = option[mc_temp_annealing]();
			deltaB = (1.0/temp_end - 1.0/temp_start) / sa_step;
		}
		TR << "Set temperature to " << 1.0/beta << std::endl;
		MonteCarloOP mc = new MonteCarlo(pose, *score_fxn, temp_start);

		//traj
		std::string output_tag(protocols::jd2::current_output_name());
		bool traj_flag = option[save_trajectory]();
		Size traj_intv = option[traj_interval]();
		protocols::canonical_sampling::PDBTrajectoryRecorder trajectory;
		if (traj_flag) {
			trajectory.file_name(output_tag + "_traj.pdb");
			trajectory.stride(traj_intv);
			trajectory.reset(*mc);
		}
		
		for (Size i=1; i<=nsteps; i++) {

			randmover->apply(pose);
			Real proposal_density_ratio=randmover->last_proposal_density_ratio();

			//TR << "Pratio: " << randmover->type() << " " << proposal_density_ratio << std::endl;

			mc->boltzmann(pose, randmover->type(), proposal_density_ratio);
			if (traj_flag) trajectory.update_after_boltzmann(*mc);
			if ( sa_step>0 && i%sa_intv==0 ) {
				mc->show_counters();
				mc->reset_counters();

				if (i<nsteps) {
					beta += deltaB;
					TR << "Set temperature to " << 1.0/beta << std::endl;
					mc->set_temperature(1.0/beta);
				}
			}
		}

		//report if haven't
		if (sa_step==0) {
			mc->show_counters();
		}

		//Trail mover
		//TrialMoverOP trial;
		//trial = new TrialMover( randmover, mc );
		//RepeatMoverOP do_cyclic = new RepeatMover( trial, 100000 );
		//do_cyclic->apply(pose);
	}
};

void *my_main( void * )
{	
	protocols::jd2::JobDistributor::get_instance()->go( new MyProtocol() );
	return 0;
}

