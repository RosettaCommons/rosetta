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

#include <basic/basic.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>
#include <basic/options/keys/bbg.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/cyclic.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
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
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>
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
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.fwd.hh>
#include <protocols/canonical_sampling/CanonicalSamplingMover.hh>
#include <protocols/moves/mc_convergence_checks/MPIBPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIHPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/MPIPool_ConvergenceCheck.hh>
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
// #include <protocols/moves/ReplicaExchangeMC.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <bitset>
//#include <devel/Gaussian/BBGaussianMover.hh>

#include <time.h>

#ifdef USEMPI
#include <mpi.h>
#endif


using core::Real;
using core::Size;
using core::id::AtomID;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR( "cycgly_bbg8ta3" );


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
	CycBBMover(Pose const & pose, Real mag) : nres_(pose.size()-2),copyres_(pose.size()-1),mag_(mag) {}
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

void BBG8T3A_sample(Pose & pose, ScoreFunctionOP sf, Size niter, Real temp = 2.0) {
	protocols::simple_moves::BBG8T3AMoverOP bbmove = new protocols::simple_moves::BBG8T3AMover();
	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sf, temp );
	mc->set_autotemp( true, temp );
	mc->set_temperature(temp);
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

void linmin(Pose & pose, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_bb(true);
	movemap->set_chi(true);
	movemap->set_jump(true);
	protocols::simple_moves::MinMover m( movemap, sf, "linmin", 1e-3, true, false, false );
	m.apply(pose);
}

Pose cyclic_perm(Pose const & orig, Size start, bool mirror=false) {
	Pose pose;
	pose.append_residue_by_jump(orig.residue(start),1);
	for(Size i = 1; i <= orig.size()-1; ++i) {
		// std::cout << "appending res " << (i+start-1)%orig.size()+1 << std::endl;
		pose.append_residue_by_bond(orig.residue((start+i-1)%orig.size()+1));
	}
	if(mirror) {
		for(Size i = 1; i <= pose.size(); ++i) {
			for(Size j = 1; j <= pose.residue_type(i).natoms(); ++j) {
				numeric::xyzVector<Real> xyz = pose.xyz(AtomID(j,i));
				xyz.z() = - xyz.z();
				pose.set_xyz(AtomID(j,i),xyz);
			}
		}
	}
	return pose;
}

// TODO: the following functions could be sped up significantly....

Real cyclic_all_atom_rmsd(Pose const & pose, Pose const & other) {
	Real mr = 9e9;
	for(Size m = 0; m <= 1; ++m) { // true false
		for(Size i = 1; i <= pose.size(); ++i) {
			Real r = core::scoring::all_atom_rmsd( cyclic_perm(pose,i,m==1), other );
			if( r < mr ) mr = r;
		}
	}
	return mr;
}

Real cyclic_ca_rmsd(Pose const & pose, Pose const & other) {
	Real mr = 9e9;
	for(Size m = 0; m <= 1; ++m) { // true false
		for(Size i = 1; i <= pose.size(); ++i) {
			Real r = core::scoring::CA_rmsd( cyclic_perm(pose,i,m==1), other );
			if( r < mr ) mr = r;
		}
	}
	return mr;
}

void cyclic_superimpose(Pose & move, Pose const & ref) {
	Real mr = 9e9;
	Size am = 0;
	Size amm = 0;
	for(Size m = 0; m <= 1; ++m) { // true false
		for(Size i = 1; i <= move.size(); ++i) {
			// Pose tmp = cyclic_perm(move,i,m==1);
			// tmp.dump_pdb("cyc_sup_test_"+ObjexxFCL::string_of(m)+"_"+ObjexxFCL::string_of(i)+".pdb");
			// Real r = core::scoring::all_atom_rmsd( tmp, ref );
			Real r = core::scoring::all_atom_rmsd( cyclic_perm(move,i,m==1), ref );
			// TR << "RMS " << m << " " << i << " " << r << std::endl;
			if( r < mr ) {
				mr = r;
				am = i;
				amm = m;
				// if(m==1) TR << "CYC SUP MIRROR!" << std::endl;
			}
		}
	}
	move = cyclic_perm(move,am,amm==1);
	core::scoring::calpha_superimpose_pose(move,ref);
	// utility_exit_with_message("testing cyclic_superimpose");
}

typedef unsigned long BINTYPE;
BINTYPE const BINBITS = 2;
BINTYPE const BINSIZE = 4; // == 2 ** BINBITS
Size    const MAXRES  = 16;

std::string printbits(BINTYPE nn) {
	std::string s = "";
	for(Size i = 0; i < 8*sizeof(nn); ++i) {
		if(i%8==0) s += " ";
		if(nn & BINTYPE(1u) << (8*sizeof(nn)-1-i)) s += "1";
		else                                       s += "0";
	}
	return s;
}

std::string bin2string(BINTYPE bin, Size nres) {
	nres = min(nres,MAXRES);
	std::string s = "";
	for(Size i = 0; i < nres; ++i) {
		for(Size j = 0; j < 2; ++j) {
			int tmp = ( bin >> 2*BINBITS*i+BINBITS*j ) % BINSIZE;
			s += ObjexxFCL::lead_zero_string_of(tmp,1);
		}
		if(i+1 < nres) s += "-";
	}
	return s;
}


BINTYPE pose2bin(core::pose::Pose const & pose) {
	using namespace ObjexxFCL::format;
	int nres = min(pose.size(),MAXRES);
	BINTYPE bin = 0;
	for(int i = 0; i < nres; ++i) {
		// Real phid = pose.phi(i+1);
		// Real psid = pose.psi(i+1);
		numeric::xyzVector<Real> c0 = pose.residue((i-1+pose.size())%pose.size()+1).xyz("C" );
		numeric::xyzVector<Real> n  = pose.residue((i  +pose.size())%pose.size()+1).xyz("N" );
		numeric::xyzVector<Real> ca = pose.residue((i  +pose.size())%pose.size()+1).xyz("CA");
		numeric::xyzVector<Real> c  = pose.residue((i  +pose.size())%pose.size()+1).xyz("C" );
		numeric::xyzVector<Real> n2 = pose.residue((i+1+pose.size())%pose.size()+1).xyz("N" );
		Real phid = basic::unsigned_periodic_range(numeric::dihedral_degrees(c0,n,ca,c),360.0);
		Real psid = basic::unsigned_periodic_range(numeric::dihedral_degrees(n,ca,c,n2),360.0);
		BINTYPE phi = phid * BINSIZE / 360.0;
		BINTYPE psi = psid * BINSIZE / 360.0;
		// std::cout << "PHIPSI " << phid << " " << psid << " binned: " << phi << " " << psi << std::endl;
		phi = phi << 2*BINBITS*i;
		psi = psi << 2*BINBITS*i+BINBITS;
		bin += phi;
		bin += psi;
	}
	// TR << bin2string(bin,pose.size()) << std::endl;
	// std::exit(-1);
	return bin;
}


Size compute_num_bins(Size nres) {
	nres = min(nres,MAXRES);
	// assumes all BINSIZE possibilities are really possible (i.e. 4 for bits, not 3)
	float ntot = std::pow((float)BINSIZE,(int)(2*nres));
	TR << "computing num bins... " << ntot << std::endl;
	std::set<BINTYPE> uniq;
	for(BINTYPE bin = 0; bin < ntot; ++bin) {
		if( bin % 1000000 == 0 ) TR << "compute_num_bins " << uniq.size() << " of " << bin << std::endl;
		// uniq.insert( cyclic_unique_bin(bin,nres) );
	}
	TR << "TOTAL (not cyclicly unique) bins: " << uniq.size() << " of total non-unique " << ntot << std::endl;
	return uniq.size();
}


void fixH(core::pose::Pose & pose) {
	for(Size i = 1; i <= pose.size(); ++i) {
		if(!pose.residue(i).has("H")) continue;
		numeric::xyzVector<Real> n  = pose.residue(i).xyz("N");
		numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
		Size in = i-1;
		if(in == 0) in = pose.size();
		numeric::xyzVector<Real> c  = pose.residue(in).xyz("C");
		numeric::xyzVector<Real> h  = n + (n-(ca+c)/2.0).normalized()*1.01;
		pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
	}
}

void cyclize_pose(core::pose::Pose & pose) {
	Size N = pose.size();
	for(Size i = 1; i <= N; ++i) {
		if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
		if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
		if(pose.residue(i).has_variant_type("CUTPOINT_UPPER")) core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_UPPER",i);
		if(pose.residue(i).has_variant_type("CUTPOINT_LOWER")) core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_LOWER",i);
	}
	if(!pose.residue(1).has_variant_type("CUTPOINT_UPPER")) core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",1);
	if(!pose.residue(N).has_variant_type("CUTPOINT_LOWER")) core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",N);
	pose.conformation().declare_chemical_bond( 1, "N", N, "C" );
	fixH(pose);
	using namespace core::scoring::constraints;
	AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.size()).atom_index("OVL1"), pose.size() );
	AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.size()).atom_index("OVL2"), pose.size() );
	AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.size()).atom_index(   "C"), pose.size() );
	pose.remove_constraints();
	pose.add_constraint(new AtomPairConstraint(a1,a2,new HarmonicFunc(0.0,0.1)));
	pose.add_constraint(new AtomPairConstraint(b1,b2,new HarmonicFunc(0.0,0.1)));
	pose.add_constraint(new AtomPairConstraint(c1,c2,new HarmonicFunc(0.0,0.1)));
}


class MoveThenFixH : public protocols::moves::Mover {
	protocols::moves::MoverOP mover_;
	ScoreFunctionOP sf_;
public:
	MoveThenFixH(protocols::moves::MoverOP mover, ScoreFunctionOP sf) : mover_(mover),sf_(sf) {}
	void apply(core::pose::Pose & pose) {
		mover_->apply(pose);
		fixH(pose);
		if(numeric::random::uniform() < 0.01) minimize(pose,sf_);
		else                                  linmin  (pose,sf_);

	}
	std::string get_name() const { return "MoveThenFixH"; }
};


void* doit(void*) {
	using namespace core::pose;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	// score functions
	ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	                sf->set_weight(core::scoring::rama,1.0);
	                sf->set_weight(core::scoring::omega,1.0);
	ScoreFunctionOP sfc = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	                sfc->set_weight(core::scoring::rama,1.0);
	                sfc->set_weight(core::scoring::atom_pair_constraint,2.0);
	                sfc->set_weight(core::scoring::omega,1.0);
	ScoreFunctionOP sffa = core::scoring::get_score_function();
	                sffa->set_weight(core::scoring::atom_pair_constraint,5.0);
	                sffa->set_weight(core::scoring::omega,1.0);
	                // sffa->set_weight(core::scoring::hbond_sr_bb,2.0); // up from 1.17
	                // sffa->set_weight(core::scoring::hbond_lr_bb,2.0); // up from 1.17
	                sffa->set_weight(core::scoring::fa_intra_rep,0.4); // up from tiny
	ScoreFunctionOP sffastd = core::scoring::get_score_function();
	                sffastd->set_weight(core::scoring::omega,1.0);
	                // sffastd->set_weight(core::scoring::hbond_sr_bb,2.0); // up from 1.17
	                // sffastd->set_weight(core::scoring::hbond_lr_bb,2.0); // up from 1.17
	                sffastd->set_weight(core::scoring::fa_intra_rep,0.4); // up from tiny


	// compute_num_bins(N);
	// utility_exit_with_message("only computing num bins! remove these lines");

	core::io::silent::SilentFileData sfd;
	Pose ref;

	TR << "setup pose & constraints" << std::endl;
	Pose pose;
	// core::pose::make_pose_from_sequence(pose,seq,core::chemical::CENTROID,false);
	core::import_pose::pose_from_file(pose,option[in::file::s]()[1], core::import_pose::PDB_file);
	Size N = pose.size();
	if( N > 4 * (Size)sizeof(BINTYPE) ) {
		utility_exit_with_message("my stupid way of binning chis uses an unsigned short/long, and can't handle that many residues!\nthere is probably a way around this to get maybe 3 more residues easily, or by making the key a pair of longs o rthe bins smaller....");
	}

	if (option[parser::view]()) protocols::viewer::add_conformation_viewer(pose.conformation(),"cycsamp",1000,1000);

	for(Size i = 1; i <= pose.size(); ++i) {
		if(pose.residue(i).name3()=="GLY" || pose.residue(i).name3()=="PRO" || pose.residue(i).name3()=="DPR" ) continue;
		core::pose::replace_pose_residue_copying_existing_coordinates(pose,i,pose.residue(i).residue_type_set().name_map("ALA"));
	}

	core::pose::replace_pose_residue_copying_existing_coordinates(pose, 4,pose.residue( 4).residue_type_set().name_map("GLY"));
	core::pose::replace_pose_residue_copying_existing_coordinates(pose,11,pose.residue(11).residue_type_set().name_map("GLY"));
	core::pose::replace_pose_residue_copying_existing_coordinates(pose, 4,pose.residue( 4).residue_type_set().name_map("DALA"));
	core::pose::replace_pose_residue_copying_existing_coordinates(pose,11,pose.residue(11).residue_type_set().name_map("DALA"));
	AtomID d1( pose.residue(4).atom_index("CB"),4), d2( pose.residue(11).atom_index("CB"),11);
	pose.add_constraint(new AtomPairConstraint(d1,d2,new HarmonicFunc(4.5,0.3)));


	(*sffa)(pose);
	cyclize_pose(pose);
	(*sffa)(pose);
	TR << "minimize" << std::endl;
	pose.dump_pdb("init.pdb");
	minimize(pose,sffa);
	pose.dump_pdb("min.pdb");

	protocols::moves::MoverOP bbmove;
	Real temp = option[cyclic::temperature]();
	TR << "using bb gaussian moves, lowering chainbreak weight!!!" << std::endl;
	bbmove = new MoveThenFixH(new protocols::simple_moves::BBG8T3AMover(),sffa);

	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sffa, temp );
	utility::vector1<Real> temps; temps.push_back(2.0); temps.push_back(4.0); temps.push_back(8.0); temps.push_back(16.0); temps.push_back(32.0);
	// utility::vector1<Real> temps; temps.push_back(2.0); temps.push_back(4.0); temps.push_back(8.0); temps.push_back(16.0);
	// utility::vector1<Real> temps; temps.push_back(2.0);
	// protocols::moves::MonteCarloOP mc = new protocols::moves::ReplicaExchangeMC( pose, *sffa, temps, 100 );
	protocols::moves::TrialMoverOP trial = new protocols::moves::TrialMover(bbmove,mc);
	std::map<BINTYPE,core::pose::PoseOP> posebins;
	std::map<BINTYPE,uint> bincount;
	BINTYPE lastbin = 0;
	int nrecorded = 0, ntransitions = 0;
	// core::pose::Pose refpose;
	core::pose::Pose last_in_bin = pose;
	core::pose::Pose start_pose = pose;
	time_t prevt = clock();
	Size itemp = temps.size();
	for(int ITER=1; ITER <= option[out::nstruct](); ITER++) {

		trial->apply( pose );

		// for testing....
		// if(numeric::random::uniform() < 0.001) pose.dump_pdb("test_"+ObjexxFCL::string_of(ITER)+".pdb");

		// temperature changes:
		if(ITER%1000 == 0) {
			itemp--;
			if(itemp==0) itemp = temps.size();
			mc->set_temperature(temps[itemp]);
		}

		if( trial->num_accepts() > nrecorded ) {// last move was an accept
			BINTYPE bin = pose2bin(pose);
			bincount[bin]++;
			if(pose.energies().total_energy() <= option[cyclic::energy_cut]()) {
				if( posebins.find(bin)==posebins.end() || (*sffa)(pose) < (*sffa)(*posebins[bin]) ) {
					if( basic::options::option[basic::options::OptionKeys::out::file::o].user() ) {
						core::pose::PoseOP tmp = new core::pose::Pose(pose);
						//if( nrecorded == 0 ) refpose = *tmp;
						//else cyclic_superimpose(*tmp,refpose);
						posebins[bin] = tmp;
					}
				}
			}

			nrecorded++;
			if(bin != lastbin) {
				ntransitions++;
				lastbin = bin;
			}

			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			ss_out->fill_struct(pose,ObjexxFCL::string_of(ITER));
			ss_out->add_energy( "rms", core::scoring::CA_rmsd(pose,start_pose) );
			sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() + ".sc" );

		}
		if( ITER % option[cyclic::log_interval]() == 0) {
			TR << ObjexxFCL::format::I(10, ITER                ) << " ITERs "
			   << ObjexxFCL::format::I(10, trial->num_accepts()) << " accpets "
			   << ObjexxFCL::format::I(10, ntransitions        ) << " transitions "
			   << ObjexxFCL::format::I(10, bincount.size()     ) << " bins sampled "
			   << ObjexxFCL::format::I(10, posebins.size()     ) << " bins filled "
			   << ObjexxFCL::format::F(7,3, sffastd->score(pose) ) << " score "
			   << ObjexxFCL::format::F(7,3, core::scoring::CA_rmsd(pose,start_pose) ) << " rms "
			   << ObjexxFCL::format::F(7,3,   sqrt(pose.energies().total_energies()[core::scoring::atom_pair_constraint]/3.0)/5.0) << " mean apc violation "
				<< Real(clock()-prevt) / 100.0 << " time "
			   << std::endl;
			prevt = clock();
			{
				utility::io::ozstream out(option[basic::options::OptionKeys::out::file::o]() + "/coverage.dat");
				for(std::map<BINTYPE,uint>::const_iterator i = bincount.begin(); i != bincount.end(); ++i) {
					out << i->first << "          " << i->second << endl;
				}
			}
			bool is_out_iter = (option[cyclic::output_interval]() != 0 && ITER % option[cyclic::output_interval]() == 0);
			if( is_out_iter || ITER == option[out::nstruct]() ) {
				core::pose::PoseOP refpose = NULL;
				TR << "dumping structs..." << std::endl;
				for(std::map<BINTYPE,core::pose::PoseOP>::const_iterator i = posebins.begin(); i != posebins.end(); ++i) {
					(*sffa)(*(i->second));
					if(!refpose) refpose = i->second;
					Pose tmp = *(i->second);
					cyclic_superimpose(tmp,*refpose);
					cyclize_pose(tmp);
					// std::string tag = ObjexxFCL::string_of(pose.size()) +"-"+ bin2string(i->first,pose.size());
					std::string tag = ObjexxFCL::string_of(tmp.size()) +"-"+ ObjexxFCL::lead_zero_string_of(i->first,10);
					tmp.dump_scored_pdb( option[basic::options::OptionKeys::out::file::o]() + "/" + tag +".pdb", *sffa );
					core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
					ss_out->fill_struct( *(i->second) ,tag);
					ss_out->add_energy( "rms", core::scoring::CA_rmsd( *(i->second), start_pose ) );
					sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() + "_ITER_"+ObjexxFCL::string_of(ITER)+".sc" );
				}

			}
		}
	}
	TR << "FINAL " << trial->num_accepts() << " accpets " << ntransitions << " transitions " << posebins.size() << " bins filled " << std::endl;
	return NULL;
}


int main( int argc, char * argv [] ) {

	try {

	using namespace core::pose;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::BBG8T3AMover::register_options();
	devel::init(argc,argv);

	void* (*func)(void*) = &doit;

	if (option[parser::view]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


