// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/rblinker.cc
/// @brief samples rigid bodies connected by linker


#include <basic/options/after_opts.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/DOF_ID.hh>
#include <devel/init.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/MPIFileBufJobDistributor.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.fwd.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.hh>
#include <protocols/moves/mc_convergence_checks/MPIBPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
//#include <devel/Gaussian/BBGaussianMover.hh>

#ifdef USEMPI
#include <mpi.h>
#endif


using core::Real;
using core::Size;
using core::id::AtomID;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR( "cycgly_canonical" );


OPT_1GRP_KEY(Real,probabilities,localbb)
OPT_1GRP_KEY(Real, probabilities,sc)
OPT_1GRP_KEY(Real, probabilities, sc_prob_uniform)
OPT_1GRP_KEY(Real, probabilities, sc_prob_withinrot)
OPT_1GRP_KEY(Real, probabilities, sc_prob_perturbcurrent)
OPT_1GRP_KEY(Boolean, probabilities, MPI_sync_pools)
OPT_1GRP_KEY(Boolean, probabilities, MPI_bcast)
OPT_1GRP_KEY(Boolean, probabilities, fast_sc_moves)
OPT_1GRP_KEY(Real, probabilities, fast_sc_moves_ntrials)
OPT_1GRP_KEY(Boolean,probabilities,no_jd2_output)
OPT_1GRP_KEY(Boolean,probabilities,use_hierarchical_clustering)
OPT_1GRP_KEY(Integer, probabilities, hierarchical_max_cache_size)


struct AbsFunc : public core::scoring::constraints::Func {
	AbsFunc( Real const x0_in, Real const sd_in ): x0_( x0_in ), sd_( sd_in ){}
	core::scoring::constraints::FuncOP
	clone() const { return new AbsFunc( *this ); }
	Real func( Real const x ) const {
		Real const z = ( x-x0_ )/sd_;
		if(z < 0) return -z;
		else return z;
	}
	Real dfunc( Real const x ) const {
		if(x-x0_ < 0) return -1.0/sd_;
		else return 1.0/sd_;
	}
	void read_data( std::istream & in ){ in >> x0_ >> sd_;  }
	void show_definition( std::ostream &out ) const { out << "ABS " << x0_ << " " << sd_ << std::endl; }
	Real x0() const { return x0_; }
	Real sd() const { return sd_; }
	void x0( Real x ) { x0_ = x; }
	void sd( Real sd ) { sd_ = sd; }
private:
	Real x0_;
	Real sd_;
};


Real mod360(Real x) {
	while(x >  180.0) x -= 360.0;
	while(x < -180.0) x += 360.0;
	return x;
}

class CycBBMover : public protocols::moves::Mover {
	Size nres_;
	Size copyres_;
	Real mag_;
public:
	CycBBMover(Pose const & pose, Real mag) : nres_(pose.n_residue()-2),copyres_(pose.n_residue()-1),mag_(mag) {}
	void apply(core::pose::Pose & pose) {
		Size i = std::ceil(numeric::random::uniform()*nres_);
		if(     numeric::random::uniform()<0.5) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		else                                    pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// if(     numeric::random::uniform()<0.45) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		// else if(numeric::random::uniform()<0.90) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// else                                     pose.set_omega(i,pose.omega(i)+numeric::random::gaussian()*mag_/10.0);
		// if     (numeric::random::uniform()<0.495) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		// else if(numeric::random::uniform()<0.990) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// else                                      pose.set_omega(i,mod360(pose.psi(i)+180.0));
		// make sure end res is identical
		if( 1 == i ) {
			pose.set_phi(copyres_,pose.phi(1));
			pose.set_psi(copyres_,pose.psi(1));
		}
		if( 2 == i ) {
			pose.set_phi(copyres_+1,pose.phi(1));
			pose.set_psi(copyres_+1,pose.psi(1));
		}
	}
	std::string get_name() const { return "CycBBMover"; }
};


void bb_sample(Pose & pose, ScoreFunctionOP sf, Size niter) {
	protocols::moves::MoverOP bbmove = new CycBBMover(pose,10.0);
	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sf, 2.0 );
	mc->set_autotemp( true, 2.0 );
	mc->set_temperature( 2.0 );
	protocols::moves::RepeatMover( new protocols::moves::TrialMover(bbmove,mc), niter ).apply( pose );
}

void minimize(Pose & pose, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true);
	movemap->set_chi(true);
	movemap->set_jump(true);
	protocols::simple_moves::MinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

Pose cyclic_perm(Pose const & orig, Size start) {
	Pose pose;
	pose.append_residue_by_jump(orig.residue(start),1);
	for(Size i = 1; i <= orig.n_residue()-1; ++i) {
		// std::cout << "appending res " << (i+start-1)%orig.n_residue()+1 << std::endl;
		pose.append_residue_by_bond(orig.residue((start+i-1)%orig.n_residue()+1));
	}
	return pose;
}

Real cyclic_all_atom_rms(Pose const & pose, Pose const & other) {
	Real mr = 9e9;
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		Real r = core::scoring::all_atom_rmsd( cyclic_perm(pose,i), other );
		if( r < mr ) mr = r;
	}
	return mr;
}

void cyclic_superimpose(Pose & move, Pose const & ref) {
	Real mr = 9e9;
	Size am = 0;
	for(Size i = 1; i <= move.n_residue(); ++i) {
		Real r = core::scoring::CA_rmsd( cyclic_perm(move,i), ref );
		if( r < mr ) {
			mr = r;
			am = i;
		}
	}
	move = cyclic_perm(move,am);
	core::scoring::calpha_superimpose_pose(move,ref);
}

int main( int argc, char * argv [] ) {

	try {


	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;

	NEW_OPT(probabilities::sc, "probability of making a side chain move", 0.25);
	NEW_OPT(probabilities::localbb, "probability of making a small move", 0.75);
	NEW_OPT(probabilities::sc_prob_uniform, "probability of uniformly sampling chi angles", 0.0);
	NEW_OPT(probabilities::sc_prob_withinrot, "probability of sampling within the current rotamer", 0.0);
	NEW_OPT(probabilities::sc_prob_perturbcurrent, "probability of perturbing the current rotamer", 0.9);
	NEW_OPT(probabilities::MPI_sync_pools, "use MPI to synchronize pools and communicate between nodes", false );
	NEW_OPT(probabilities::MPI_bcast, "use broadcasting in syncing", false );
	NEW_OPT(probabilities::fast_sc_moves, "use the fast SidechainMCMover", false);
	NEW_OPT(probabilities::fast_sc_moves_ntrials, "specify the number of ntrials for each call of scmover apply", 1000);
	NEW_OPT(probabilities::no_jd2_output, "do not write to silent-file specified by -out:file:silent", false );
	NEW_OPT(probabilities::use_hierarchical_clustering, "use the HierarchicalLevel class",false);
	NEW_OPT(probabilities::hierarchical_max_cache_size, "set the max-cache size of the hierarchy", 100);
	protocols::canonical_sampling::CanonicalSamplingMover::register_options();
	BBG8T3AMover::register_options();
	devel::init(argc,argv);

	std::string seq = "G"; while((int)seq.size() < option[cyclic::nres]()) seq += "G";

	// score functions
	ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	                sf->set_weight(core::scoring::rama,1.0);
	                sf->set_weight(core::scoring::omega,1.0);
	ScoreFunctionOP sfc = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	                sfc->set_weight(core::scoring::rama,1.0);
	                sfc->set_weight(core::scoring::atom_pair_constraint,10.0);
	                sfc->set_weight(core::scoring::omega,1.0);
	ScoreFunctionOP sffa = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	                sffa->set_weight(core::scoring::atom_pair_constraint,10.0);
	                sffa->set_weight(core::scoring::omega,1.0);
	ScoreFunctionOP sffastd = core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	                sffastd->set_weight(core::scoring::omega,1.0);

	core::io::silent::SilentFileData sfd;
	Pose ref;

	// liz stuff
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
	//protocols::canonical_sampling::CanonicalSamplingMoverOP csm(new protocols::canonical_sampling::CanonicalSamplingMover(sfxn,pool_ptr,1000));
	protocols::canonical_sampling::CanonicalSamplingMoverOP csm(new CanonicalSamplingMover);
	csm->set_scorefunction(sfxn);
	csm->use_hierarchical_clustering(true);

	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb( true );
	movemap->set_chi( true );
	BBG8T3AMoverOP bbg8t3mover = new BBG8T3AMover();
	bbg8t3mover->movemap(movemap);
	csm->add_mover(bbg8t3mover,option[probabilities::localbb]);

	mc_convergence_checks::Pool_RMSD_OP pool_ptr = new mc_convergence_checks::Pool_RMSD();
   csm->set_poolrmsd(pool_ptr);

	// end liz stuff

	TR << "setup pose & constraints" << std::endl;
	Pose pose;
	core::pose::make_pose_from_sequence(pose,seq,core::chemical::CENTROID,false);
	core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",       1        );
	core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",pose.n_residue());
	pose.conformation().declare_chemical_bond( 1, "N", pose.n_residue(), "C" );
	// Size const nbb( lower_rsd.mainchain_atoms().size() );
	// total_dev +=
	// 	( upper_rsd.atom( upper_rsd.mainchain_atoms()[  1] ).xyz().distance_squared( lower_rsd.atom( "OVL1" ).xyz() ) +
	// 	  upper_rsd.atom( upper_rsd.mainchain_atoms()[  2] ).xyz().distance_squared( lower_rsd.atom( "OVL2" ).xyz() ) +
	// 	  lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz().distance_squared( upper_rsd.atom( "OVU1" ).xyz() ) );
	TR << " " << std::endl;
	{
		AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.n_residue()).atom_index("OVL1"), pose.n_residue() );
		AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.n_residue()).atom_index("OVL2"), pose.n_residue() );
		AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.n_residue()).atom_index(   "C"), pose.n_residue() );
		pose.add_constraint(new AtomPairConstraint(a1,a2,new AbsFunc(0.0,0.01)));
		pose.add_constraint(new AtomPairConstraint(b1,b2,new AbsFunc(0.0,0.01)));
		pose.add_constraint(new AtomPairConstraint(c1,c2,new AbsFunc(0.0,0.01)));
	}
	// gen structure
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(numeric::random::uniform() < 0.0) pose.set_omega(i,  0.0);
		else                                 pose.set_omega(i,180.0);
	}

	TR << "bb sample for reasonable closed starting point" << std::endl;
	bb_sample(pose,sf ,100);
	bb_sample(pose,sfc,1000);


	TR << "switch to fa" << std::endl;
	core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_UPPER",1);
	core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_LOWER",pose.n_residue());
	protocols::toolbox::switch_to_residue_type_set(pose,"fa_standard");
	core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",       1        );
	core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",pose.n_residue());
	{
		AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.n_residue()).atom_index("OVL1"), pose.n_residue() );
		AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.n_residue()).atom_index("OVL2"), pose.n_residue() );
		AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.n_residue()).atom_index(   "C"), pose.n_residue() );
		pose.remove_constraints();
		pose.add_constraint(new AtomPairConstraint(a1,a2,new AbsFunc(0.0,0.01)));
		pose.add_constraint(new AtomPairConstraint(b1,b2,new AbsFunc(0.0,0.01)));
		pose.add_constraint(new AtomPairConstraint(c1,c2,new AbsFunc(0.0,0.01)));
	}
	TR << "minimize" << std::endl;
	minimize(pose,sffa);

	TR << "enter canonical sampling mover" << std::endl;
	csm->apply(pose);


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

