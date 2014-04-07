// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/rblinker.cc
/// @brief samples rigid bodies connected by linker

#include <basic/basic.hh>
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
#include <core/import_pose/import_pose.hh>
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
#include <ObjexxFCL/FArray2D.hh>
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

basic::Tracer TR("cycgly_superimpose");

//
// OPT_1GRP_KEY(Real,probabilities,localbb)
// OPT_1GRP_KEY(Real, probabilities,sc)
// OPT_1GRP_KEY(Real, probabilities, sc_prob_uniform)
// OPT_1GRP_KEY(Real, probabilities, sc_prob_withinrot)
// OPT_1GRP_KEY(Real, probabilities, sc_prob_perturbcurrent)
// OPT_1GRP_KEY(Boolean, probabilities, MPI_sync_pools)
// OPT_1GRP_KEY(Boolean, probabilities, MPI_bcast)
// OPT_1GRP_KEY(Boolean, probabilities, fast_sc_moves)
// OPT_1GRP_KEY(Real, probabilities, fast_sc_moves_ntrials)
// OPT_1GRP_KEY(Boolean,probabilities,no_jd2_output)
// OPT_1GRP_KEY(Boolean,probabilities,use_hierarchical_clustering)
// OPT_1GRP_KEY(Integer, probabilities, hierarchical_max_cache_size)
//
//
//
// struct AbsFunc : public core::scoring::constraints::Func {
// 	AbsFunc( Real const x0_in, Real const sd_in ): x0_( x0_in ), sd_( sd_in ){}
// 	core::scoring::constraints::FuncOP
// 	clone() const { return new AbsFunc( *this ); }
// 	Real func( Real const x ) const {
// 		Real const z = ( x-x0_ )/sd_;
// 		if(z < 0) return -z;
// 		else return z;
// 	}
// 	Real dfunc( Real const x ) const {
// 		if(x-x0_ < 0) return -1.0/sd_;
// 		else return 1.0/sd_;
// 	}
// 	void read_data( std::istream & in ){ in >> x0_ >> sd_;  }
// 	void show_definition( std::ostream &out ) const { out << "ABS " << x0_ << " " << sd_ << std::endl; }
// 	Real x0() const { return x0_; }
// 	Real sd() const { return sd_; }
// 	void x0( Real x ) { x0_ = x; }
// 	void sd( Real sd ) { sd_ = sd; }
// private:
// 	Real x0_;
// 	Real sd_;
// };
//
//
// Real mod360(Real x) {
// 	while(x >  180.0) x -= 360.0;
// 	while(x < -180.0) x += 360.0;
// 	return x;
// }
//
// class CycBBMover : public protocols::moves::Mover {
// 	Size nres_;
// 	Size copyres_;
// 	Real mag_;
// public:
// 	CycBBMover(Pose const & pose, Real mag) : nres_(pose.n_residue()-2),copyres_(pose.n_residue()-1),mag_(mag) {}
// 	void apply(core::pose::Pose & pose) {
// 		Size i = std::ceil(numeric::random::uniform()*nres_);
// 		if(     numeric::random::uniform()<0.5) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
// 		else                                    pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
// 		// if(     numeric::random::uniform()<0.45) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
// 		// else if(numeric::random::uniform()<0.90) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
// 		// else                                     pose.set_omega(i,pose.omega(i)+numeric::random::gaussian()*mag_/10.0);
// 		// if     (numeric::random::uniform()<0.495) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
// 		// else if(numeric::random::uniform()<0.990) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
// 		// else                                      pose.set_omega(i,mod360(pose.psi(i)+180.0));
// 		// make sure end res is identical
// 		if( 1 == i ) {
// 			pose.set_phi(copyres_,pose.phi(1));
// 			pose.set_psi(copyres_,pose.psi(1));
// 		}
// 		if( 2 == i ) {
// 			pose.set_phi(copyres_+1,pose.phi(1));
// 			pose.set_psi(copyres_+1,pose.psi(1));
// 		}
// 	}
// 	std::string get_name() const { return "CycBBMover"; }
// };
//
//
// void bb_sample(Pose & pose, ScoreFunctionOP sf, Size niter) {
// 	protocols::moves::MoverOP bbmove = new CycBBMover(pose,10.0);
// 	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sf, 2.0 );
// 	mc->set_autotemp( true, 2.0 );
// 	mc->set_temperature( 2.0 );
// 	protocols::moves::RepeatMover( new protocols::moves::TrialMover(bbmove,mc), niter ).apply( pose );
// }
//
// void BBG8T3A_sample(Pose & pose, ScoreFunctionOP sf, Size niter, Real temp = 2.0) {
// 	protocols::simple_moves::BBG8T3AMoverOP bbmove = new protocols::simple_moves::BBG8T3AMover();
// 	protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sf, temp );
// 	mc->set_autotemp( true, temp );
// 	mc->set_temperature(temp);
// 	protocols::moves::RepeatMover( new protocols::moves::TrialMover(bbmove,mc), niter ).apply( pose );
// }
//
// void minimize(Pose & pose, ScoreFunctionOP sf) {
// 	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
// 	movemap->set_bb(true);
// 	movemap->set_chi(true);
// 	movemap->set_jump(true);
// 	protocols::simple_moves::MinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
// 	m.apply(pose);
// }

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

Real cyclic_ca_rms(Pose const & pose, Pose const & other) {
	Real mr = 9e9;
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		Real r = core::scoring::CA_rmsd( cyclic_perm(pose,i), other );
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

// typedef unsigned short BINTYPE;
// BINTYPE const BINBITS = 1;
// BINTYPE const BINSIZE = 2; // == 2 ** BINBITS
//
// std::string printbits(BINTYPE nn) {
// 	std::string s = "";
// 	for(Size i = 0; i < 8*sizeof(nn); ++i) {
// 		if(i%8==0) s += " ";
// 		if(nn & BINTYPE(1u) << (8*sizeof(nn)-1-i)) s += "1";
// 		else                                       s += "0";
// 	}
// 	return s;
// }
//
// std::string bin2string(BINTYPE bin, Size nres) {
// 	std::string s = "";
// 	for(Size i = 0; i < nres; ++i) {
// 		for(Size j = 0; j < 2; ++j) {
// 			int tmp = ( bin >> 2*BINBITS*i+BINBITS*j ) % BINSIZE;
// 			s += ObjexxFCL::lead_zero_string_of(tmp,1);
// 		}
// 		if(i+1 < nres) s += "-";
// 	}
// 	return s;
// }
//
// BINTYPE makemask(Size n) {
// 	BINTYPE ONE = 1u;
// 	BINTYPE mask = 0;
// 	for(Size i = 0; i < n; ++i) {
// 		mask = mask | ONE << i;
// 	}
// 	return mask;
// }
//
// BINTYPE cyclic_unique_bin(BINTYPE bin, Size nres) {
// 	// for debuggung....
// 	// TR << "testmask " << printbits( makemask(8) ) << std::endl;
// 	// TR << "testmask " << printbits( makemask(13) ) << std::endl;
// 	// TR << "testmask " << printbits( MASK ) << std::endl;
// 	// bin = 1UL<<0*9 | 1UL<<1*9 | 1UL<<2*9 | 1UL<<3*9 | 1UL<<4*9 | 1UL<<5*9;// | 1UL<<6*9 ;
// 	// TR << "bin " << printbits(bin) << std::endl;
// 	BINTYPE MASK = makemask(2*BINBITS*nres);
// 	BINTYPE mn = bin;
// 	for(Size i = 0; i < nres; ++i) {
// 	 	BINTYPE tmp = (bin << 2*BINBITS*i) | (bin >> (2*BINBITS*( nres - i) )); //
// 		tmp = tmp & MASK;
// 		// TR << i << " " << printbits(tmp) << std::endl;
// 		if(tmp < mn) mn = tmp;
// 	}
// 	// TR << "m " << printbits(mn) << std::endl;
// 	// utility_exit_with_message("debug pose2bin");
// 	return mn;
// }
//
// BINTYPE pose2bin(core::pose::Pose const & pose) {
// 	using namespace ObjexxFCL::format;
// 	BINTYPE bin = 0;
// 	for(int i = 0; i < (int)pose.n_residue(); ++i) {
// 		// Real phid = pose.phi(i+1);
// 		// Real psid = pose.psi(i+1);
// 		numeric::xyzVector<Real> c0 = pose.residue((i-1+pose.n_residue())%pose.n_residue()+1).xyz("C" );
// 		numeric::xyzVector<Real> n  = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("N" );
// 		numeric::xyzVector<Real> ca = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("CA");
// 		numeric::xyzVector<Real> c  = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("C" );
// 		numeric::xyzVector<Real> n2 = pose.residue((i+1+pose.n_residue())%pose.n_residue()+1).xyz("N" );
// 		Real phid = basic::unsigned_periodic_range(numeric::dihedral_degrees(c0,n,ca,c),360.0);
// 		Real psid = basic::unsigned_periodic_range(numeric::dihedral_degrees(n,ca,c,n2),360.0);
// 		// std::cout << "PHIPSI " << phid << " " << psid << std::endl;
// 		BINTYPE phi = phid / 180.0;
// 		BINTYPE psi = psid / 180.0;
// 		phi = phi << 2*BINBITS*i;
// 		psi = psi << 2*BINBITS*i+BINBITS;
// 		bin += phi;
// 		bin += psi;
// 	}
// 	// TR << bin2string(bin,pose.n_residue()) << std::endl;
// 	// std::exit(-1);
// 	return cyclic_unique_bin(bin,pose.n_residue());
// }

// std::string printbits(BINTYPE nn) {
// 	std::string s = "";
// 	for(Size i = 0; i < 8*sizeof(nn); ++i) {
// 		if(i%8==0) s += " ";
// 		if(nn & BINTYPE(1) << (8*sizeof(nn)-1-i)) s += "1";
// 		else                               s += "0";
// 	}
// 	return s;
// }
//
// std::string bin2string(BINTYPE bin, Size nres) {
// 	std::string s = "";
// 	for(Size i = 0; i < nres; ++i) {
// 		for(Size j = 0; j < 2; ++j) {
// 			int tmp = ( bin >> 2*BINBITS*i+BINBITS*j ) % BINSIZE;
// 			s += ObjexxFCL::lead_zero_string_of(tmp,1);
// 		}
// 		if(i+1 < nres) s += "-";
// 	}
// 	return s;
// }
//
// BINTYPE pose2bin(core::pose::Pose const & pose) {
// 	using namespace ObjexxFCL::format;
// 	BINTYPE bin = 0;
// 	for(int i = 0; i < (int)pose.n_residue(); ++i) {
// 		// Real phid = pose.phi(i+1);
// 		// Real psid = pose.psi(i+1);
// 		numeric::xyzVector<Real> c0 = pose.residue((i-1+pose.n_residue())%pose.n_residue()+1).xyz("C" );
// 		numeric::xyzVector<Real> n  = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("N" );
// 		numeric::xyzVector<Real> ca = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("CA");
// 		numeric::xyzVector<Real> c  = pose.residue((i  +pose.n_residue())%pose.n_residue()+1).xyz("C" );
// 		numeric::xyzVector<Real> n2 = pose.residue((i+1+pose.n_residue())%pose.n_residue()+1).xyz("N" );
// 		Real phid = numeric::dihedral_degrees(c0,n,ca,c);
// 		Real psid = numeric::dihedral_degrees(n,ca,c,n2);
// 		// TR << phid << " " << psid << std::endl;
// 		BINTYPE phi = basic::unsigned_periodic_range(phid,360.0) / 180.0;
// 		BINTYPE psi = basic::unsigned_periodic_range(psid,360.0) / 180.0;
// 		phi = phi << 2*BINBITS*i;
// 		psi = psi << 2*BINBITS*i+BINBITS;
// 		bin += phi;
// 		bin += psi;
// 	}
// 	// TR << bin2string(bin,pose.n_residue()) << std::endl;
// 	// std::exit(-1);
// 	// for debuggung....
// 	// bin = 1UL<<0*9 | 1UL<<1*9 | 1UL<<2*9 | 1UL<<3*9 | 1UL<<4*9 | 1UL<<5*9 | 1UL<<6*9 ;
// 	BINTYPE mn = bin;
// 	for(Size i = 0; i < pose.n_residue(); ++i) {
// 	 	BINTYPE tmp = (bin << 2*BINBITS*i) | (bin >> (2*BINBITS*pose.n_residue() - 2*BINBITS*i));
// 		tmp = tmp & ~(BINTYPE(BINSIZE*BINSIZE-1)<<((sizeof(BINTYPE)-1)*2*BINBITS));
// 		// TR << i << " " << printbits(tmp) << std::endl;
// 		if(tmp < mn) mn = tmp;
// 	}
// 	// TR << "m " << printbits(mn) << std::endl;
// 	// utility_exit_with_message("debug pose2bin");
// 	return mn;
// }

// Size compute_num_bins(Size nres) {
// 	// assumes all BINSIZE possibilities are really possible (i.e. 4 for bits, not 3)
// 	float ntot = std::pow((float)BINSIZE,(int)(2*nres));
// 	std::set<BINTYPE> uniq;
// 	for(BINTYPE bin = 0; bin < ntot; ++bin) {
// 		uniq.insert( cyclic_unique_bin(bin,nres) );
// 	}
// 	TR << "TOTAL cyclicly unique bins: " << uniq.size() << " of total non-unique " << ntot << std::endl;
// 	return uniq.size();
// }

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

	utility::vector1<Pose> targets_in;
	targets_in = core::import_pose::poses_from_pdbs( option[in::file::s]() );
	for(Size i = 1; i <= targets_in.size(); i++) {
		core::pose::remove_lower_terminus_type_from_pose_residue(targets_in[i],           1             );
		core::pose::remove_upper_terminus_type_from_pose_residue(targets_in[i],targets_in[i].n_residue());
	}
	ObjexxFCL::FArray2D_float rmsmat(targets_in.size(),targets_in.size(),0.0);
	TR << "computing central structure out of " << targets_in.size() << std::endl;
	Real mn = 9e9;
	Size argmin = 1;
	for(Size i = 1; i <= targets_in.size(); i++) {
		Real x = float(i)/float(targets_in.size());
		TR << "computing rms for " << i << ", " << (2*x-x*x)*100 << " percent complete" << std::endl;
		for(Size j = i+1; j <= targets_in.size(); j++) {
			rmsmat(i,j) = cyclic_all_atom_rms( targets_in[i], targets_in[j] );
			rmsmat(j,i) = rmsmat(i,j);
		}
	}

	for(Size i = 1; i <= targets_in.size(); i++) {
		Real meanrms = 0;
		for(Size j = 1; j <= targets_in.size(); j++) {
			meanrms += rmsmat(i,j);
		}
		meanrms /= float(targets_in.size());
		if( meanrms < mn ) {
			mn = meanrms;
			argmin = i;
		}
	}

	for(Size i = 1; i <= targets_in.size(); ++i) {
		cyclic_superimpose(targets_in[i],targets_in[argmin]);
		targets_in[i].dump_pdb("aligned_"+ObjexxFCL::lead_zero_string_of(i,6)+".pdb");
	}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}






















