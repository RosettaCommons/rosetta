// -*- mode:c++;tab-width:1;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
// AUTO-REMOVED #include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AngleConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/MultiConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
// AUTO-REMOVED #include <core/scoring/func/XYZ_Func.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <protocols/simple_moves/FragmentMover.hh>
// AUTO-REMOVED #include <numeric/kinematic_closure/bridgeObjects.hh>
// AUTO-REMOVED #include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
// AUTO-REMOVED #include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SymMinMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/TrialMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

#include <protocols/moves/MoverStatistics.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::Real;
using core::Size;
using core::pose::Pose;
using core::kinematics::Stub;
using core::conformation::Residue;
using core::conformation::ResidueOP;
using protocols::scoring::ImplicitFastClashCheck;
using std::string;
using utility::vector1;
using numeric::min;
using core::import_pose::pose_from_pdb;
using basic::options::option;
using numeric::min;
using numeric::max;
using utility::io::ozstream;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;

static basic::Tracer TR("IKFGDB");

inline Real const sqr(Real const r) { return r*r; }

inline Vec xyz(Pose const & p, Size const & ia, Size const & ir) {
return p.xyz(AtomID(ia,ir));
}

Pose & alapose(Pose & pose) {
for(Size i=1; i<=pose.n_residue(); ++i) {
core::pose::replace_pose_residue_copying_existing_coordinates(pose,i,pose.residue(i).residue_type_set().name_map("ALA"));
}
return pose;
}

vector1<Size> allifnone(vector1<Size> v, Size n) {
if( v.size()==0 ) {
v.resize(n);
for(Size i = 1; i <=n; ++i) v[i] = i;
}
return(v);
}

class RotSel : public utility::pointer::ReferenceCount {
protected:
  Pose const pose_;
  vector1<Size> const pos_;
  ImplicitFastClashCheck const ifc_;
  utility::vector1<Stub> stb_;
  utility::vector1<Vec> cb_;
  core::chemical::ResidueTypeSetCAP frs_;
public:
  RotSel( Pose p_in, vector1<Size> & pos )
    : pose_(alapose(p_in)), pos_(allifnone(pos,pose_.n_residue())), ifc_(pose_,2.7)
  {
    frs_ = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
    for(vector1<Size>::const_iterator i=pos_.begin(); i!=pos_.end(); ++i) {
      stb_.push_back(Stub(xyz(pose_,5,*i),xyz(pose_,2,*i),xyz(pose_,1,*i)));
      cb_.push_back(xyz(pose_,5,*i));
    }
  }
  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const = 0;
};

class BruteRotSel : public RotSel {
public:
  BruteRotSel( Pose p_in, vector1<Size> pos=vector1<Size>() ) : RotSel(p_in,pos) {}
  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const {
    vector1<Size>::const_iterator ip = pos_.begin();
    vector1<Vec> ::const_iterator ic = cb_ .begin();
    Pose asn;
    make_pose_from_sequence(asn,"N","fa_standard",false);
    if(asn.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(asn,1);
    if(asn.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(asn,1);

    for(vector1<Stub>::const_iterator is = stb_.begin(); is!=stb_.end(); ++is,++ip,++ic) {
      Real cbd = cen0.distance_squared(*ic);
      if(cbd > 36.0) continue;
						ResidueOP bestrsd;
						Real bestsc = 9e9;
      for(Real ch1 = 0; ch1 <= 360; ch1+=1.0) {
        core::conformation::ResidueOP rsd = core::conformation::ResidueFactory::create_residue(frs_->name_map("ASN"),pose_.residue(*ip),pose_.conformation());
        rsd->set_chi(1,ch1);
        Vec cb = rsd->xyz("CB");
        Vec cg = rsd->xyz("CG");
        Vec cen = cb+2.5*(cg-cb);
								Real const dsq = cen.distance_squared(cen0);
        if( dsq > disth*disth) continue;
								Real const dt = (cg-cb).normalized().dot(axs0);
        if( dt < dotth ) continue;
        if( dsq-10*dt < bestsc ) {
          bestrsd = rsd;
										bestsc = dsq-10*dt;
        }

      }
      if(bestsc < 9e8) hits.push_back(bestrsd);
    }
  }
};

class KinRotSel : public RotSel {
public:
  KinRotSel( Pose p_in, vector1<Size> pos=vector1<Size>() ) : RotSel(p_in,pos) {}

  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const {
    Real const & pi  (numeric::constants::d::pi);
    Real const & pi_2(numeric::constants::d::pi_2);
    Real const dis2ub = sqr(        3.75894266112+disth );
    Real const dis2lb = sqr(max(0.0,3.75894266112-disth));
    Real const cendotang( asin(disth/3.75894266112) );
    Real const cendotlb(cos(1.985117-cendotang));
    Real const cendotub(cos(1.985117+cendotang));

    vector1<Size>::const_iterator ip = pos_.begin();
    vector1<Vec> ::const_iterator ic = cb_ .begin();
    for(vector1<Stub>::const_iterator is = stb_.begin(); is!=stb_.end(); ++is,++ip,++ic) {
						if(*ip != 43) continue;
      Real cbd = cen0.distance_squared(*ic);
      if( dis2lb > cbd || cbd > dis2ub ) {TR<<"cb dis fail"<<std::endl;continue;}
      //TR << *ip << " " << cbd << " " << dis2lb << " " << dis2ub << std::endl;
      Vec const cen(is->global2local(cen0));
      if( cen.x() < 0 )  {TR<<"cen x < 0"<<std::endl;continue;}
      //Real const cenl2 = cen.x()*cen.x()+cen.y()*cen.y()+cen.z()*cen.z();
      //Real cendot = cen.x()*cen.x() / cenl2;
      //if( cendotlb > cendot || cendot > cendotub ) continue; // asn specific check

      Vec const axs(is->M.transposed()*axs0);
      Vec const cbcg(Vec( 0.4397778*sqrt(cen.y()*cen.y()+cen.z()*cen.z()),cen.y(),cen.z()).normalized());
      //      TR << cbcg << std::endl;
      Real cendotcg = cbcg.dot(axs); //0.4397778==tan(ang(ca,cb,cg))
      if( cendotcg < dotth )  {TR<<"cendot fail"<<std::endl;continue;}
      // compute chi1
      Real const ch1tmp = acos(cen.y()/sqrt(cen.y()*cen.y()+cen.z()*cen.z()));
      Real const ch1 = (cen.z()>0)?ch1tmp:pi_2-ch1tmp;
      //      Real const ch2 =
      core::conformation::ResidueOP rsd = core::conformation::ResidueFactory::create_residue(frs_->name_map("ASN"),pose_.residue(*ip),pose_.conformation());
      rsd->set_chi(1,numeric::conversions::degrees(ch1));
      if(cen0.distance_squared( rsd->xyz(6)+2.5*(rsd->xyz(6)-rsd->xyz(5))) > disth*disth )  {TR<<"real dis fail"<<std::endl;continue;}
      hits.push_back(rsd);


    }
  }
};


core::pack::rotamer_set::RotamerSetOP get_rotset(Pose & pose, Size icys) {
  core::pack::rotamer_set::RotamerSetOP rotset;
  core::scoring::ScoreFunction dummy_sfxn;
  dummy_sfxn( pose );
  core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
  dummy_task->initialize_from_command_line();
  dummy_task->nonconst_residue_task( icys ).and_extrachi_cutoff(10);
  dummy_task->nonconst_residue_task( icys ).restrict_to_repacking();
  dummy_task->nonconst_residue_task( icys ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
  dummy_task->nonconst_residue_task( icys ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
  core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
  core::pack::rotamer_set::RotamerSetFactory rsf;
  rotset = rsf.create_rotamer_set( pose.residue( icys ) );
  rotset->set_resid( icys );
  rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
  return rotset;
}


void run(std::string fname) {
  using namespace std;
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using namespace core;
  using namespace chemical;
  using namespace id;
  using namespace pose;
  using namespace scoring;

  Pose gln;
  make_pose_from_sequence(gln,"Q","fa_standard",false);
  if(gln.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(gln,1);
  if(gln.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(gln,1);

  // setup stuff
  ResidueTypeSetCAP frs=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  ResidueTypeSetCAP crs=ChemicalManager::get_instance()->residue_type_set( "centroid" );
  ScoreFunctionOP sfstd=getScoreFunction();
  ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");

  Pose pose;
  pose_from_pdb(pose,option[in::file::s]()[1]);
  pose = alapose(pose);
  RotSel & krs(*(new   KinRotSel(pose)));
  RotSel & brs(*(new BruteRotSel(pose)));
  TR << "start" << std::endl;

		//  for(Size igln=1;igln<=pose.n_residue();++igln) {
  for(Size igln=6;igln<=6;++igln) {
    TR << igln << std::endl;
    core::conformation::Residue rtmp(pose.residue(igln));
    pose.replace_residue(igln,gln.residue(1),true);
    core::pack::rotamer_set::RotamerSetOP rots( get_rotset(pose,igln) );
    //for(Size irot = 1; irot <= rots->num_rotamers(); ++irot) {
    for(Size irot = 10; irot <= 10; ++irot) {
      //TR << igln << " " << irot << std::endl;
      pose.set_chi(1,igln,rots->rotamer(irot)->chi(1));
      pose.set_chi(2,igln,rots->rotamer(irot)->chi(2));
      pose.set_chi(3,igln,rots->rotamer(irot)->chi(3));
      //      if(igln==13) pose.dump_pdb("test.pdb");
      Vec cg = pose.residue(igln).xyz(6);
      Vec cd = pose.residue(igln).xyz(7);
      Vec axs = (cg-cd).normalized();
      Vec cen = cd - 1.7*axs;
      vector1<core::conformation::ResidueOP> hits;

      krs.place_n( cen, axs, Vec(0,1,0), hits, 0.7, 0.9 );
      for(Size ih=1;ih<=hits.size();++ih) {
        core::conformation::Residue rtmp2(pose.residue(hits[ih]->seqpos()));
        pose.replace_residue(hits[ih]->seqpos(),*hits[ih],true);
        pose.dump_pdb("KIN_"+lzs(igln,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(hits[ih]->seqpos(),3)+".pdb");
        pose.replace_residue(hits[ih]->seqpos(),rtmp2,false);
      }

      brs.place_n( cen, axs, Vec(0,1,0), hits, 0.7, 0.9 );
      for(Size ih=1;ih<=hits.size();++ih) {
        core::conformation::Residue rtmp2(pose.residue(hits[ih]->seqpos()));
        pose.replace_residue(hits[ih]->seqpos(),*hits[ih],true);
        pose.dump_pdb("BRT_"+lzs(igln,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(hits[ih]->seqpos(),3)+".pdb");
        pose.replace_residue(hits[ih]->seqpos(),rtmp2,false);
      }

    }
    pose.replace_residue(igln,rtmp,false);
  }

  TR << "stop" << std::endl;



  // utility::io::ozstream out("test.pdb");
  // Pose res;
  // make_pose_from_sequence(res,"N","fa_standard",false);
  // core::pose::remove_lower_terminus_type_from_pose_residue(res,1);
  // core::pose::remove_upper_terminus_type_from_pose_residue(res,1);
  // Stub s(res.xyz(AtomID(5,1)),res.xyz(AtomID(2,1)),res.xyz(AtomID(1,1)));
  // xform_pose_rev(res,s);
  // for(Real i1=0; i1<360; i1+= 10.0) {
  //   res.set_chi(1,1,i1);
  //   for(Real i2=0; i2<360; i2+= 10.0) {
  //     res.set_chi(2,1,i2);
  //     Vec fgcen = 2.5*res.xyz(AtomID(6,1)); // cen mid-HB, tilted 4.11416 dg. towards OC
  //     res.dump_pdb(out);
  //     out<<"HETATM"<<I(5,9997)<<' '<<"XXXX"<<' '<< "XXX"<<' '<<"A"<<I(4,995)<<"    "<<F(8,3, fgcen.x())<<F(8,3, fgcen.y())<<F(8,3, fgcen.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
  //     out.close();      utility_exit_with_message("narsiotnra");
  //   }
  // }
  // out.close();
}


int main (int argc, char *argv[]) {

	try {



  devel::init(argc,argv);

  using basic::options::option;
  using namespace basic::options::OptionKeys;

  run(option[in::file::s]()[1]);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



