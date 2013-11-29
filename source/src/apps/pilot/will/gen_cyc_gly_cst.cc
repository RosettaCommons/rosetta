// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/rblinker.cc
/// @brief samples rigid bodies connected by linker

#include <basic/options/option.hh>
#include <basic/options/keys/cyclic.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <devel/init.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/NumericTraits.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

using core::Real;
using core::Size;
using core::id::AtomID;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;

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


class BBMover : public protocols::moves::Mover {
	Size start_,stop_;
	Real mag_;
public:
	BBMover(Size start, Size stop, Real mag) : start_(start),stop_(stop),mag_(mag) {}
	void apply(core::pose::Pose & pose) {
		using namespace numeric::random;
		Size i = start_-1 + std::ceil(uniform()*(stop_-start_+1));
		if(uniform()<0.5) pose.set_phi(i,pose.phi(i)+gaussian()*mag_);
		else              pose.set_psi(i,pose.psi(i)+gaussian()*mag_);
	}
	std::string get_name() const { return "BBMover"; }
};


void bb_sample(Pose & pose, ScoreFunctionOP sf, Size niter) {
  protocols::moves::MoverOP bbmove = new BBMover(1,pose.n_residue(),10.0);
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
  protocols::simple_moves::MinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
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

  const Real PI = numeric::NumericTraits<Real>::pi();

  devel::init(argc,argv);

  std::string seq = "G"; while((int)seq.size() < option[cyclic::nres]()) seq += "G";

  // score functions
  ScoreFunctionOP sf = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
  sf->set_weight(core::scoring::rama,1.0);
  ScoreFunctionOP sfc = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
  sfc->set_weight(core::scoring::rama,1.0);
  sfc->set_weight(core::scoring::atom_pair_constraint,1.0);
  sfc->set_weight(core::scoring::angle_constraint    ,1.0);
  sfc->set_weight(core::scoring::dihedral_constraint ,1.0);
  ScoreFunctionOP sffa = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
  sffa->set_weight(core::scoring::atom_pair_constraint,1.0);
  sfc->set_weight(core::scoring::angle_constraint     ,1.0);
  sfc->set_weight(core::scoring::dihedral_constraint  ,1.0);
  sffa->set_weight(core::scoring::omega,2.0);
  ScoreFunctionOP sffastd = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
  sffa->set_weight(core::scoring::omega,2.0);

  for(Size ITER = 1; ITER <= (Size)option[out::nstruct](); ++ITER) {
    while(true) {

      // setup pose
      Pose pose;
      core::pose::make_pose_from_sequence(pose,seq,"centroid",false);
			core::pose::add_variant_type_to_pose_residue(pose,"VIRTUAL_GLY",1);
			Size N = pose.n_residue();
			pose.add_constraint(new AtomPairConstraint(            AtomID(1,1),AtomID(3,N)            ,new core::scoring::constraints::HarmonicFunc(1.3       ,0.001)));
			pose.add_constraint(new AngleConstraint   (AtomID(5,1),AtomID(1,1),AtomID(3,N)            ,new core::scoring::constraints::HarmonicFunc(2.08      ,0.0001)));
			pose.add_constraint(new AngleConstraint   (            AtomID(1,1),AtomID(3,N),AtomID(4,N),new core::scoring::constraints::HarmonicFunc(2.1467    ,0.0001)));
			pose.add_constraint(new DihedralConstraint(AtomID(5,1),AtomID(1,1),AtomID(3,N),AtomID(4,N),new core::scoring::constraints::HarmonicFunc(PI   ,0.0001)));

      // gen structure
      for(Size i = 1; i <= pose.n_residue(); ++i) pose.set_omega(i,180.0);
      bb_sample(pose,sf,100);
      bb_sample(pose,sfc,1000);
      protocols::toolbox::switch_to_residue_type_set(pose,"fa_standard");
      minimize(pose,sffa);

      if( (*sffastd)(pose)/(pose.n_residue()) > -1.0 ) {
				std::cout << "retry!" << (*sffastd)(pose) << std::endl;
				continue;
      }

      //pose.delete_polymer_residue(pose.n_residue());
      //pose.delete_polymer_residue(pose.n_residue());
      pose.dump_pdb("cyc_gly_"+ObjexxFCL::string_of(pose.n_residue())+"_"+ObjexxFCL::string_of(ITER)+".pdb");

      std::cout << "finish " << ITER << std::endl;
      break;
    }
  }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

