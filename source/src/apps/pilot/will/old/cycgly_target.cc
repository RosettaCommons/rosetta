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
#include <core/pack/task/PackerTask.hh>
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
#include <protocols/simple_moves/PackRotamersMover.hh>
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

static basic::Tracer TR( "cycgly_target" );


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
		if ( z < 0 ) return -z;
		else return z;
	}
	Real dfunc( Real const x ) const {
		if ( x-x0_ < 0 ) return -1.0/sd_;
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
	while ( x >  180.0 ) x -= 360.0;
	while ( x < -180.0 ) x += 360.0;
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
		if (     numeric::random::uniform()<0.5 ) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		else                                    pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// if(     numeric::random::uniform()<0.45) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		// else if(numeric::random::uniform()<0.90) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// else                                     pose.set_omega(i,pose.omega(i)+numeric::random::gaussian()*mag_/10.0);
		// if     (numeric::random::uniform()<0.495) pose.set_phi(i,pose.phi(i)+numeric::random::gaussian()*mag_);
		// else if(numeric::random::uniform()<0.990) pose.set_psi(i,pose.psi(i)+numeric::random::gaussian()*mag_);
		// else                                      pose.set_omega(i,mod360(pose.psi(i)+180.0));
		// make sure end res is identical
		if ( 1 == i ) {
			pose.set_phi(copyres_,pose.phi(1));
			pose.set_psi(copyres_,pose.psi(1));
		}
		if ( 2 == i ) {
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

Pose cyclic_perm(Pose const & orig, Size start, bool mirror=false) {
	Pose pose;
	pose.append_residue_by_jump(orig.residue(start),1);
	for ( Size i = 1; i <= orig.size()-1; ++i ) {
		// std::cout << "appending res " << (i+start-1)%orig.size()+1 << std::endl;
		pose.append_residue_by_bond(orig.residue((start+i-1)%orig.size()+1));
	}
	if ( mirror ) {
		for ( Size i = 1; i <= pose.size(); ++i ) {
			for ( Size j = 1; j <= pose.residue_type(i).natoms(); ++j ) {
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
	for ( Size m = 0; m <= 1; ++m ) { // true false
		for ( Size i = 1; i <= pose.size(); ++i ) {
			Real r = core::scoring::all_atom_rmsd( cyclic_perm(pose,i,m==1), other );
			if ( r < mr ) mr = r;
		}
	}
	return mr;
}

Real cyclic_ca_rmsd(Pose const & pose, Pose const & other) {
	Real mr = 9e9;
	for ( Size m = 0; m <= 1; ++m ) { // true false
		for ( Size i = 1; i <= pose.size(); ++i ) {
			Real r = core::scoring::CA_rmsd( cyclic_perm(pose,i,m==1), other );
			if ( r < mr ) mr = r;
		}
	}
	return mr;
}

void cyclic_superimpose(Pose & move, Pose const & ref) {
	Real mr = 9e9;
	Size am = 0;
	Size amm = 0;
	for ( Size m = 0; m <= 1; ++m ) { // true false
		for ( Size i = 1; i <= move.size(); ++i ) {
			// Pose tmp = cyclic_perm(move,i,m==1);
			// tmp.dump_pdb("cyc_sup_test_"+ObjexxFCL::string_of(m)+"_"+ObjexxFCL::string_of(i)+".pdb");
			// Real r = core::scoring::all_atom_rmsd( tmp, ref );
			Real r = core::scoring::all_atom_rmsd( cyclic_perm(move,i,m==1), ref );
			// TR << "RMS " << m << " " << i << " " << r << std::endl;
			if ( r < mr ) {
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

std::string printbits(BINTYPE nn) {
	std::string s = "";
	for ( Size i = 0; i < 8*sizeof(nn); ++i ) {
		if ( i%8==0 ) s += " ";
		if ( nn & BINTYPE(1u) << (8*sizeof(nn)-1-i) ) s += "1";
		else                                       s += "0";
	}
	return s;
}

std::string bin2string(BINTYPE bin, Size nres) {
	std::string s = "";
	for ( Size i = 0; i < nres; ++i ) {
		for ( Size j = 0; j < 2; ++j ) {
			int tmp = ( bin >> 2*BINBITS*i+BINBITS*j ) % BINSIZE;
			s += ObjexxFCL::lead_zero_string_of(tmp,1);
		}
		if ( i+1 < nres ) s += "-";
	}
	return s;
}

BINTYPE makemask(Size n) {
	BINTYPE ONE = 1u;
	BINTYPE mask = 0;
	for ( Size i = 0; i < n; ++i ) {
		mask = mask | ONE << i;
	}
	return mask;
}

BINTYPE mirror_bin(BINTYPE bin, Size nres) {
	BINTYPE two = 2u;
	BINTYPE m = 0;
	for ( Size i = 0; i < nres; ++i ) {
		for ( Size j = 0; j < 2; ++j ) {
			BINTYPE tmp = (( bin >> 2*BINBITS*i+BINBITS*j ) ) % BINSIZE;
			// std::cerr << "tmp        " << tmp << std::endl;
			tmp = BINSIZE - tmp - 1;
			// std::cerr << "tmp mirror " << tmp << std::endl;
			m += (tmp << 2*BINBITS*i+BINBITS*j);
		}
	}
	// TR << "MIRROR TEST " << bin << std::endl;
	// TR << "MIRROR TEST " << m   << std::endl;
	return m;
}

BINTYPE cyclic_unique_bin(BINTYPE bin, Size nres) {
	// for debuggung....
	// TR << "testmask " << printbits( makemask(8) ) << std::endl;
	// TR << "testmask " << printbits( makemask(13) ) << std::endl;
	// TR << "testmask " << printbits( MASK ) << std::endl;
	// bin = 1UL<<0*9 | 1UL<<1*9 | 1UL<<2*9 | 1UL<<3*9 | 1UL<<4*9 | 1UL<<5*9;// | 1UL<<6*9 ;
	// TR << "bin " << printbits(bin) << std::endl;
	BINTYPE MASK = makemask(2*BINBITS*nres);
	BINTYPE mn = bin;
	for ( Size m = 0; m <= 1; ++m ) { // true false
		for ( Size i = 0; i < nres; ++i ) {
			BINTYPE tmp = (bin << 2*BINBITS*i) | (bin >> (2*BINBITS*( nres - i) )); //
			tmp = tmp & MASK;
			if ( 1==m ) tmp = mirror_bin(tmp,nres);
			// TR << i << " " << printbits(tmp) << std::endl;
			if ( tmp < mn ) mn = tmp;
		}
	}
	// TR << "m " << printbits(mn) << std::endl;
	// utility_exit_with_message("debug pose2bin");
	return mn;
}

BINTYPE pose2bin(core::pose::Pose const & pose) {
	using namespace ObjexxFCL::format;
	BINTYPE bin = 0;
	for ( int i = 0; i < (int)pose.size(); ++i ) {
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
	return cyclic_unique_bin(bin,pose.size());
}


Size compute_num_bins(Size nres) {
	// assumes all BINSIZE possibilities are really possible (i.e. 4 for bits, not 3)
	float ntot = std::pow((float)BINSIZE,(int)(2*nres));
	TR << "computing num bins... " << ntot << std::endl;
	std::set<BINTYPE> uniq;
	for ( BINTYPE bin = 0; bin < ntot; ++bin ) {
		if ( bin % 1000000 == 0 ) TR << "compute_num_bins " << uniq.size() << " of " << bin << std::endl;
		uniq.insert( cyclic_unique_bin(bin,nres) );
	}
	TR << "TOTAL cyclicly unique bins: " << uniq.size() << " of total non-unique " << ntot << std::endl;
	return uniq.size();
}


void fixH(core::pose::Pose & pose) {
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( !pose.residue(i).has("H") ) continue;
		numeric::xyzVector<Real> n  = pose.residue(i).xyz("N");
		numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
		Size in = i-1;
		if ( in == 0 ) in = pose.size();
		numeric::xyzVector<Real> c  = pose.residue(in).xyz("C");
		numeric::xyzVector<Real> h  = n + (n-(ca+c)/2.0).normalized()*1.01;
		pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
	}
}

void cyclize_pose(core::pose::Pose & pose) {
	Size N = pose.size();
	for ( Size i = 1; i <= N; ++i ) {
		if ( pose.residue(i).is_lower_terminus() ) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
		if ( pose.residue(i).is_upper_terminus() ) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
		if ( pose.residue(i).has_variant_type("CUTPOINT_UPPER") ) core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_UPPER",i);
		if ( pose.residue(i).has_variant_type("CUTPOINT_LOWER") ) core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_LOWER",i);
	}
	if ( !pose.residue(1).has_variant_type("CUTPOINT_UPPER") ) core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",1);
	if ( !pose.residue(N).has_variant_type("CUTPOINT_LOWER") ) core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",N);
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
public:
	MoveThenFixH(protocols::moves::MoverOP mover) : mover_(mover) {}
	void apply(core::pose::Pose & pose) {
		mover_->apply(pose);
		fixH(pose);
	}
	std::string get_name() const { return "MoveThenFixH"; }
};

class MoveThenFixHMin : public protocols::moves::Mover {
	protocols::moves::MoverOP mover_,minmover_,packmover_;
public:
	MoveThenFixHMin(core::pose::Pose const & pose, ScoreFunctionOP sf, protocols::moves::MoverOP mover) : mover_(mover) {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(true);
		minmover_ = new protocols::simple_moves::MinMover( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );

		using namespace core::pack::task;
		PackerTaskOP task = TaskFactory::create_packer_task(pose);
		task->initialize_extra_rotamer_flags_from_command_line();
		for ( Size i = 1; i <= pose.size(); ++i ) {
			task->nonconst_residue_task(i).restrict_to_repacking();
		}
		// TR << *task << std::endl;
		packmover_ = new protocols::simple_moves::PackRotamersMover( sf, task );
	}
	void apply(core::pose::Pose & pose) {
		mover_->apply(pose);
		fixH(pose);
		if ( numeric::random::uniform() < 0.1 ) {
			packmover_->apply(pose);
			minmover_->apply(pose);
		}
	}
	std::string get_name() const { return "MoveThenFixHMin"; }
};

// void target_constraints_lcc(core::pose::Pose & pose, core::pose::Pose const & refpose, utility::vector1<utility::vector1<Size> > refres) {
// LocalCoordinateConstraint(
//   id::AtomID const & a1,
//  id::StubID const & fixed_stub_in,
//  Vector const & xyz_target_in,
//   FuncOP func,
// using namespace core::scoring::constraints;
// for(utility::vector1<utility::vector1<Size> >::iterator rr = refres.begin(); rr != refres.end(); ++rr) {
//  for(Size i = 1; i <= rr->size(); ++i) {
//   for(Size j = i+1; j <= rr->size(); ++j) {
//    Size ir = (*rr)[i];
//    Size jr = (*rr)[j];
//    for(Size ia = 1; ia <= 3; ia++) {
//     for(Size ja = 1; ja <= 3; ja++) {
//      AtomID a1(ia,ir),a2(ja,jr);
//      Real d = refpose.xyz(a1).distance(refpose.xyz(a2));
//      pose.add_constraint(new AtomPairConstraint(a1,a2, new HarmonicFunc(d,1.0) ));
//      TR << "add cst " << ir << "," << ia << " " << jr << "," << ja << " " << d << std::endl;
//     }
//    }
//   }
//  }
// }
// }

void target_constraints(core::pose::Pose & pose, core::pose::Pose const & refpose, utility::vector1<utility::vector1<Size> > poseres, utility::vector1<utility::vector1<Size> > refres) {
	using namespace core::scoring::constraints;
	runtime_assert(poseres.size()==refres.size());
	for ( Size ii = 1; ii <= poseres.size(); ++ii ) {
		runtime_assert(poseres[ii].size()==refres[ii].size());
		for ( Size i = 1; i <= poseres[ii].size(); ++i ) {
			for ( Size j = i+1; j <= poseres[ii].size(); ++j ) {
				Size pir = poseres[ii][i];
				Size pjr = poseres[ii][j];
				Size rir = refres [ii][i];
				Size rjr = refres [ii][j];
				for ( Size ia = 1; ia <= 4; ia++ ) {
					for ( Size ja = 1; ja <= 4; ja++ ) {
						AtomID pa1(ia,pir),pa2(ja,pjr),ra1(ia,rir),ra2(ja,rjr);
						Real d = refpose.xyz(ra1).distance(refpose.xyz(ra2));
						pose.add_constraint(new AtomPairConstraint(pa1,pa2, new HarmonicFunc(d,4.0) ));
						TR << "add cst " << pir << "," << ia << " " << pjr << "," << ja << " " << d << std::endl;
					}
				}
			}
		}
	}
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
		Size N = option[cyclic::nres]();

		// score functions
		ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		sf->set_weight(core::scoring::rama,1.0);
		sf->set_weight(core::scoring::omega,1.0);
		ScoreFunctionOP sfc = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		sfc->set_weight(core::scoring::rama,1.0);
		sfc->set_weight(core::scoring::atom_pair_constraint,10.0);
		sfc->set_weight(core::scoring::omega,1.0);
		ScoreFunctionOP sffa = core::scoring::get_score_function();
		sffa->set_weight(core::scoring::atom_pair_constraint,1.0);
		sffa->set_weight(core::scoring::omega,1.0);
		sffa->set_weight(core::scoring::hbond_sr_bb,2.0); // up from 1.17
		ScoreFunctionOP sffastd = core::scoring::get_score_function();
		sffastd->set_weight(core::scoring::omega,1.0);
		sffastd->set_weight(core::scoring::hbond_sr_bb,2.0); // up from 1.17


		utility::vector1<Pose> targets_in;
		std::set<BINTYPE> targetbins;
		std::map<BINTYPE,Pose> targets;
		std::map<BINTYPE,Pose> closest,furthest;
		if ( option[in::file::s].user() ) {
			targets_in = core::import_pose::poses_from_files( option[in::file::s]() , core::import_pose::PDB_file);
			for ( Size i = 1; i <= targets_in.size(); i++ ) {
				if ( targets_in[i].size() != N ) {
					utility_exit_with_message("wrong number of residues (not "+ObjexxFCL::string_of(N)+") in target " + std::string(option[in::file::s]()[i]) );
				}
				cyclize_pose( targets_in[i] );
				BINTYPE bin = pose2bin(targets_in[i]);
				targetbins.insert( bin );
				targets[bin] = targets_in[i];
				// TR << "target: " << option[in::file::s]()[1] << " " << bin << " " << bin2string(bin,targets[bin].size()) << std::endl;
				std::string fn = basic::options::option[basic::options::OptionKeys::out::file::o]() + "/target_"+ObjexxFCL::lead_zero_string_of(bin,10)+".pdb";
				targets[bin].dump_pdb(fn);

				// TR << "BIN orig " << printbits(bin) << std::endl;
				// bin = mirror_bin( bin, N );
				// TR << "BIN mirr " << printbits(bin) << std::endl;
				// bin = mirror_bin( bin, N );
				// TR << "BIN orig " << printbits(bin) << std::endl;
				//
				// TR << "ORIG:   " << std::endl;
				// BINTYPE test1 = pose2bin(targets_in[i]);
				// TR << "BIN test1 " << printbits(test1) << std::endl;
				// Pose mirror = cyclic_perm( targets_in[i], 1, true );
				// mirror.dump_pdb(fn+"_mirror1.pdb");
				//
				// TR << "MIRROR: " << std::endl;
				// BINTYPE test2 = pose2bin(mirror);
				// TR << "BIN test2 " << printbits(test2) << std::endl;
			}
			// utility_exit_with_message("will is testing a new feature! oh boy!");
		}


		if ( N > 4 * (Size)sizeof(BINTYPE) ) {
			utility_exit_with_message("my stupid way of binning chis uses an unsigned short/long, and can't handle that many residues!\nthere is probably a way around this to get maybe 3 more residues easily, or by making the key a pair of longs o rthe bins smaller....");
		}

		// compute_num_bins(N);
		// utility_exit_with_message("only computing num bins! remove these lines");

		std::string seq = "G"; while ( (Size)seq.size() < N ) seq += "G";

		core::io::silent::SilentFileData sfd;
		Pose ref;

		TR << "setup pose & constraints" << std::endl;
		Pose pose = targets_in[1];
		// core::pose::make_pose_from_sequence(pose,seq,core::chemical::CENTROID,false);
		core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",       1        );
		core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",pose.size());
		pose.conformation().declare_chemical_bond( 1, "N", pose.size(), "C" );
		// Size const nbb( lower_rsd.mainchain_atoms().size() );
		// total_dev +=
		//  ( upper_rsd.atom( upper_rsd.mainchain_atoms()[  1] ).xyz().distance_squared( lower_rsd.atom( "OVL1" ).xyz() ) +
		//    upper_rsd.atom( upper_rsd.mainchain_atoms()[  2] ).xyz().distance_squared( lower_rsd.atom( "OVL2" ).xyz() ) +
		//    lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz().distance_squared( upper_rsd.atom( "OVU1" ).xyz() ) );
		TR << " " << std::endl;
		{
			AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.size()).atom_index("OVL1"), pose.size() );
			AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.size()).atom_index("OVL2"), pose.size() );
			AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.size()).atom_index(   "C"), pose.size() );
			pose.add_constraint(new AtomPairConstraint(a1,a2,new AbsFunc(0.0,0.02)));
			pose.add_constraint(new AtomPairConstraint(b1,b2,new AbsFunc(0.0,0.02)));
			pose.add_constraint(new AtomPairConstraint(c1,c2,new AbsFunc(0.0,0.02)));
		}

		// // gen structure
		// for(Size i = 1; i <= pose.size(); ++i) {
		//  pose.set_phi(i,0.0);
		//  pose.set_psi(i,0.0);
		//  if(numeric::random::uniform() < 0.0) pose.set_omega(i,  0.0);
		//  else                                 pose.set_omega(i,180.0);
		// }
		//
		// TR << "bb sample for reasonable closed starting point" << std::endl;
		// bb_sample(pose,sf ,100);
		// bb_sample(pose,sfc,1000);
		//
		//
		// TR << "switch to fa" << std::endl;
		// core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_UPPER",1);
		// core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_LOWER",pose.size());
		// protocols::toolbox::switch_to_residue_type_set(pose,"fa_standard");
		// cyclize_pose(pose);
		// TR << "minimize" << std::endl;
		// minimize(pose,sffa);
		// // BBG8T3A_sample(pose,sffa,100000,50.0);
		// // minimize(pose,sffa);
		// // pose.dump_pdb("min.pdb");

		protocols::moves::MoverOP bbmove;
		Real temp = option[cyclic::temperature]();
		// if( basic::options::option[ basic::options::OptionKeys::backrub::backrub ].user() ) {
		//  TR << "using backrub moves" << std::endl;
		//  core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_UPPER",1);
		//  core::pose::remove_variant_type_from_pose_residue(pose,"CUTPOINT_LOWER",pose.size());
		//  core::pose::add_lower_terminus_type_to_pose_residue(pose,1);
		//  core::pose::add_upper_terminus_type_to_pose_residue(pose,pose.size());
		//  protocols::backrub::BackrubMoverOP brmove = new protocols::backrub::BackrubMover();
		//  brmove->branchopt().read_database();
		//  brmove->clear_segments();
		//  core::pose::PoseOP pose_ptr( new core::pose::Pose() );
		//  *pose_ptr = pose;
		//  brmove->set_input_pose( pose_ptr );
		//  utility::vector1< Size > pep_nbrs;
		//  //backrub all pep nbrs plus one up/downstream
		//  for( Size ii = 1; ii <= pose.size(); ++ii ) pep_nbrs.push_back( ii );
		//  utility::vector1< std::string > backrub_atomnames( pep_nbrs.size(), "CA" );
		//  brmove->add_mainchain_segments( pep_nbrs, backrub_atomnames, 3, 5 );
		//  brmove->optimize_branch_angles( pose );
		//  core::pose::remove_lower_terminus_type_from_pose_residue(pose,1);
		//  core::pose::remove_upper_terminus_type_from_pose_residue(pose,pose.size());
		//  core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_UPPER",       1        );
		//  core::pose::add_variant_type_to_pose_residue(pose,"CUTPOINT_LOWER",pose.size());
		//  TR << "done setting up BackrubMover" << std::endl;
		//  bbmove = brmove;
		// } else {
		TR << "using bb gaussian moves, lowering chainbreak weight!!!" << std::endl;
		bbmove = new MoveThenFixHMin(pose,sffa,new protocols::simple_moves::BBG8T3AMover());
		// AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.size()).atom_index("OVL1"), pose.size() );
		// AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.size()).atom_index("OVL2"), pose.size() );
		// AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.size()).atom_index(   "C"), pose.size() );
		// pose.remove_constraints();
		// pose.add_constraint(new AtomPairConstraint(a1,a2,new AbsFunc(0.0,0.01)));
		// pose.add_constraint(new AtomPairConstraint(b1,b2,new AbsFunc(0.0,0.01)));
		// pose.add_constraint(new AtomPairConstraint(c1,c2,new AbsFunc(0.0,0.01)));
		// }

		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CHEAT!
		pose = targets_in[1];
		utility::vector1<utility::vector1<Size> > prr(2);
		for ( Size i =  4; i <=  9; ++i ) prr[1].push_back(i);
		for ( Size i =  1; i <=  2; ++i ) prr[2].push_back(i);
		for ( Size i = 11; i <= 14; ++i ) prr[2].push_back(i);
		// utility::vector1<utility::vector1<Size> > rr(1);
		// for(Size i =  1; i <=  14; ++i) rr[1].push_back(i);
		target_constraints(pose,targets_in[1],prr,prr);
		sffa->show(pose);
		// pose.constraint_set()->show(std::cout);

		protocols::moves::MonteCarloOP mc = new protocols::moves::MonteCarlo( pose, *sffa, temp );
		utility::vector1<Real> temps; temps.push_back(1.0); temps.push_back(2.0); temps.push_back(4.0); temps.push_back(8.0); temps.push_back(16.0);
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


		for ( int ITER=1; ITER <= option[out::nstruct](); ITER++ ) {

			trial->apply( pose );

			// for testing....
			// if(numeric::random::uniform() < 0.001) pose.dump_pdb("test_"+ObjexxFCL::string_of(ITER)+".pdb");

			// temperature changes:
			if ( ITER%100 == 0 ) {
				itemp--;
				if ( itemp==0 ) itemp = temps.size();
				mc->set_temperature(temps[itemp]);
			}

			// do cyclic perms at some interval
			// if(ITER%10000 == 0) {
			//  Size news = (Size)std::ceil(numeric::random::uniform()*N);
			//  pose = cyclic_perm( pose, news, numeric::random::uniform() < 0.5 );
			//  cyclize_pose(pose);
			// }

			if ( trial->num_accepts() > nrecorded ) { // last move was an accept
				BINTYPE bin = pose2bin(pose);
				bincount[bin]++;
				if ( pose.energies().total_energy() <= option[cyclic::energy_cut]() ) {
					// if( posebins.find(bin)==posebins.end() || (*sffa)(pose) < (*sffa)(*posebins[bin]) ) {
					if ( posebins.find(bin)==posebins.end() || core::scoring::CA_rmsd(pose,targets_in[1]) > core::scoring::CA_rmsd(*posebins[bin],targets_in[1]) ) {
						if ( basic::options::option[basic::options::OptionKeys::out::file::o].user() ) {
							core::pose::PoseOP tmp = new core::pose::Pose(pose);
							//if( nrecorded == 0 ) refpose = *tmp;
							//else cyclic_superimpose(*tmp,refpose);
							posebins[bin] = tmp;
						}
					}
				}

				// if( option[cyclic::single_bin]() ) {
				//  if( bin == lastbin ) {
				//   last_in_bin = pose;
				//   pose.dump_pdb("SINGLE_BIN_"+ObjexxFCL::lead_zero_string_of(bin,10)+"_"+ObjexxFCL::lead_zero_string_of(ITER,10)+".pdb");
				//  } else {
				//                 pose = last_in_bin;
				//  }
				// }
				// if( option[cyclic::output_all_bin].user() ) {
				//  if( bin == BINTYPE(option[cyclic::output_all_bin]()) ) {
				//   pose.dump_pdb(basic::options::option[basic::options::OptionKeys::out::file::o]() + "/ALL_IN_BIN_"+ObjexxFCL::lead_zero_string_of(bin,10)+"_"+ObjexxFCL::lead_zero_string_of(ITER,10)+".pdb");
				//  }
				// }
				if ( ITER%10 == 0 ) {
					// write silent struct
					sffastd->score(pose);
					core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
					ss_out_all->fill_struct(pose,ObjexxFCL::string_of(ITER));
					ss_out_all->add_energy( "tgt_rms", core::scoring::CA_rmsd(pose,targets_in[1]) );
					sfd.write_silent_struct( *ss_out_all, option[ out::file::silent ]() );
				}

				nrecorded++;
				if ( bin != lastbin ) {
					ntransitions++;
					lastbin = bin;
				}
				// only use if with targets... for testing
				// if( targetbins.count(bin) > 0 ) {
				//  std::string fn = basic::options::option[basic::options::OptionKeys::out::file::o]() + "/target_hit_"+ObjexxFCL::lead_zero_string_of(bin,10)+"_"+ObjexxFCL::lead_zero_string_of(ITER,10)+".pdb";
				//  Real rms = cyclic_all_atom_rmsd(pose,targets[bin]);
				//  if( closest.find(bin)==closest.end() || rms < cyclic_all_atom_rmsd(closest[bin],targets[bin])) {
				//   TR << "CLOSEST TARGET HIT! recording " << bin << " " << rms << std::endl;
				//   Pose tmp = pose;
				//   // cyclic_superimpose(tmp,targets[bin]);
				//   closest[bin] = tmp;
				//  }
				//  if( furthest.find(bin)==furthest.end() || rms > cyclic_all_atom_rmsd(furthest[bin],targets[bin])) {
				//   TR << "FURTHEST TARGET HIT! recording " << bin << " " << rms << std::endl;
				//   Pose tmp = pose;
				//   // cyclic_superimpose(tmp,targets[bin]);
				//   furthest[bin] = tmp;
				//  }
				// }
			}
			if ( ITER % option[cyclic::log_interval]() == 0 ) {
				TR << ObjexxFCL::format::I(10, ITER                ) << " ITERs "
					<< ObjexxFCL::format::I(10, trial->num_accepts()) << " accpets "
					<< ObjexxFCL::format::I(10, ntransitions        ) << " transitions "
					<< ObjexxFCL::format::I(10, bincount.size()     ) << " bins sampled "
					<< ObjexxFCL::format::I(10, posebins.size()     ) << " bins filled "
					// << ObjexxFCL::format::F(7,3, cyclic_all_atom_rmsd(pose,start_pose) ) << " cyc_rms "
					<< ObjexxFCL::format::F(7,3, core::scoring::CA_rmsd(pose,targets_in[1]) ) << " rms "
					// << ObjexxFCL::format::F(7,3, core::scoring::all_atom_rmsd(pose,start_pose) ) << " rms "
					<< ObjexxFCL::format::F(7,3,  pose.energies().total_energies()[core::scoring::atom_pair_constraint]) << " cst "
					<< Real(clock()-prevt) / 100.0 << " time "
					<< std::endl;
				prevt = clock();
				{
					utility::io::ozstream out(option[basic::options::OptionKeys::out::file::o]() + "/coverage.dat");
					for ( std::map<BINTYPE,uint>::const_iterator i = bincount.begin(); i != bincount.end(); ++i ) {
						out << i->first << "          " << i->second << endl;
					}
				}
				bool is_out_iter = (option[cyclic::output_interval]() != 0 && ITER % option[cyclic::output_interval]() == 0);
				if ( is_out_iter || ITER == option[out::nstruct]() ) {
					core::pose::PoseOP refpose = NULL;
					TR << "dumping structs..." << std::endl;
					for ( std::map<BINTYPE,core::pose::PoseOP>::const_iterator i = posebins.begin(); i != posebins.end(); ++i ) {
						(*sffa)(*(i->second));
						if ( !refpose ) refpose = i->second;
						Pose tmp = *(i->second);
						cyclic_superimpose(tmp,*refpose);
						cyclize_pose(tmp);
						// std::string tag = ObjexxFCL::string_of(pose.size()) +"-"+ bin2string(i->first,pose.size());
						std::string tag = ObjexxFCL::string_of(tmp.size()) +"-"+ ObjexxFCL::lead_zero_string_of(i->first,10);
						tmp.dump_pdb( option[basic::options::OptionKeys::out::file::o]() + "/" + tag +".pdb");
						core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
						ss_out->fill_struct( *(i->second) ,tag);
						sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() + "_ITER_"+ObjexxFCL::string_of(ITER)+".sc" );
					}
					for ( std::map<BINTYPE,core::pose::Pose>::const_iterator i = closest.begin(); i != closest.end(); ++i ) {
						std::string tag = ObjexxFCL::string_of(pose.size()) +"-"+ ObjexxFCL::lead_zero_string_of(i->first,10);
						Pose tmp = i->second;
						cyclic_superimpose(tmp,targets[i->first]);
						cyclize_pose(tmp);
						tmp.dump_pdb( basic::options::option[basic::options::OptionKeys::out::file::o]() + "/closest_" + tag +".pdb");
					}
					for ( std::map<BINTYPE,core::pose::Pose>::const_iterator i = furthest.begin(); i != furthest.end(); ++i ) {
						std::string tag = ObjexxFCL::string_of(pose.size()) +"-"+ ObjexxFCL::lead_zero_string_of(i->first,10);
						Pose tmp = i->second;
						cyclic_superimpose(tmp,targets[i->first]);
						cyclize_pose(tmp);
						tmp.dump_pdb( basic::options::option[basic::options::OptionKeys::out::file::o]() + "/furthest_" + tag +".pdb");
					}

				}
			}
		}
		TR << "FINAL " << trial->num_accepts() << " accpets " << ntransitions << " transitions " << posebins.size() << " bins filled " << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
