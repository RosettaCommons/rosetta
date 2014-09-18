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
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
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


// AUTO-REMOVED #include <time.h>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <apps/pilot/will/will_util.ihh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

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
using numeric::constants::d::pi;
using numeric::conversions::degrees;
using utility::io::ozstream;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

typedef utility::vector1<core::Real> Reals;
typedef utility::vector1<core::Size> Sizes;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef utility::vector1<Vec> Vecs;

static thread_local basic::Tracer TR( "IKFGDB" );

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


void dumpcgo(Vec v, string l) {
  std::cerr << "cmd.load_cgo(Vec( " << v.x() << "," << v.y() << "," << v.z() << ").cgo(), '"+l+"')"  << std::endl;
}

void xform_rsd_gl2(Stub const & s, Residue & rsd) {
  for(Size i = 1; i <= rsd.natoms(); ++i) rsd.set_xyz(i, s.global2local(rsd.xyz(i)) );
}

enum KRSQueryType {
  CEN,
  CEN_ANG,
  CEN_AXS,
  CEN_AXS_ORI
};

struct KRSQuery {
  KRSQueryType type;
  Vec  cen;
  Vec  axs;
  Vec  ori;
  Real disth;
  Real angth;
  Real clash;
  KRSQuery(KRSQueryType  typ                                                                            ) : type(typ),cen(0,0,0),axs(1,0,0),ori(0,1,0),disth(0.5),angth(0.0175),clash( 2.8) {}
  KRSQuery(KRSQueryType  typ, Vec  c, Vec  a, Vec  o, Real  dt=0.5, Real  at=0.175, Real  clsh = 2.8*2.8) : type(typ),cen(  c  ),axs(  a  ),ori(  o  ),disth( dt),angth(  at  ),clash(clsh) {}
  KRSQuery(KRSQueryType  typ, Vec  c, Vec  a,         Real  dt=0.5, Real  at=0.175, Real  clsh = 2.8*2.8) : type(typ),cen(  c  ),axs(  a  ),ori(0,0,0),disth( dt),angth(  at  ),clash(clsh) {}
};

class FunGroupTK : public utility::pointer::ReferenceCount {
protected:
  Pose const pose_;
  vector1<Size> const pos_;
  ImplicitFastClashCheck const ifc_;
  std::map<std::string,vector1<Stub> > stb_;
		std::map<std::string,vector1<ResidueOP> > rsd_;
  core::chemical::ResidueTypeSetCAP frs_;
public:
  FunGroupTK( Pose p_in, vector1<Size> & pos )
    : pose_(alapose(p_in)), pos_(allifnone(pos,pose_.n_residue())), ifc_(pose_,2.2)
  {
    frs_ = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
    vector1<string> res_types;
    res_types.push_back("ASN");
    res_types.push_back("CYS");
    res_types.push_back("HIS");
    res_types.push_back("GLN");
				for(vector1<string>::const_iterator it = res_types.begin(); it != res_types.end(); ++it) {
						rsd_[*it].resize(pose_.n_residue());
						for(vector1<Size>::const_iterator i=pos_.begin(); i!=pos_.end(); ++i) {
        ResidueOP rsd = core::conformation::ResidueFactory::create_residue(frs_->name_map(*it),pose_.residue(*i),pose_.conformation());
        stb_[*it].push_back(Stub(rsd->xyz(5),rsd->xyz(2),rsd->xyz(1)));
								rsd_[*it][*i] = rsd;
      }
    }
  }
  ImplicitFastClashCheck const & ifc() const { return ifc_; }
  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const = 0;
  virtual void place_c(KRSQuery const & q, Residue const & qrsd,vector1<core::conformation::ResidueOP> & hits ) const = 0;
};

class BruteFunGroupTK : public FunGroupTK {
public:
  BruteFunGroupTK( Pose p_in, vector1<Size> pos=vector1<Size>() ) : FunGroupTK(p_in,pos) {}
  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const {
  }
  virtual void place_c(KRSQuery const & q, Residue const & qrsd, vector1<core::conformation::ResidueOP> & hits ) const {
    vector1<Size>::const_iterator ip = pos_.begin();
    vector1<Stub> const & stb(stb_.find("CYS")->second);
    for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
      Real cbd2 = q.cen.distance_squared(is->v);
      if(*ip==qrsd.seqpos()) continue;
      if( cbd2 > 14.0 ) continue;
      ResidueOP bestrsd;
      Real bestsc = 9e9;
      ResidueOP rsd = rsd_.find("CYS")->second[*ip];
      rsd->set_xyz(11, rsd->xyz("SG") + 2.1*(rsd->xyz("HG")-rsd->xyz("SG")).normalized() );
      for(Real ch1 = 0; ch1 <= 360; ch1+=1.0) {
        rsd->set_chi(1,ch1);
        Vec sg = rsd->xyz("SG");
        if( sg.distance_squared(q.cen) > 8.0 ) continue;
        for(Real ch2 = 0; ch2 <= 360; ch2+=1.0) {
          rsd->set_chi(2,ch2);
          Vec hg = rsd->xyz("HG");
          Vec cen = sg+2.1*((hg-sg).normalized());
          Real const dsq = cen.distance_squared(q.cen);
          if( dsq > q.disth*q.disth) continue;
          Vec axs = ((sg-q.cen).normalized());
          Real const da = numeric::angle_radians( axs, Vec(0,0,0), q.axs);

          // Vec qc = is->global2local(q.cen);
          // qc = qc * 2.911171 / qc.length();
          // ozstream out("test_"+str(qc.x())+"_"+lzs(qrsd.seqpos(),3)+"_"+lzs(*ip,3)+"_"+lzs((Size)ch1,3)+"_"+lzs((Size)ch2,3)+".pdb");
          // Size ano=0;
          // core::io::pdb::dump_pdb_residue(*rsd,ano,out);
          // out.close();

          if( fabs(da-q.ori.x()) > q.angth ) continue;
          bool clash = false;
          for(Size ia = 6; ia <= rsd->nheavyatoms(); ++ia) {
            if(!ifc_.clash_check(rsd->xyz(ia),*ip) ) { clash = true; break; }
            for(Size ja = 1; ja <= qrsd.nheavyatoms(); ++ja) if( qrsd.xyz(ja).distance_squared(rsd->xyz(ia)) < q.clash ) { clash = true; break; }
            if(clash) break;
          }
          if(clash) continue;
          if( da < bestsc ) {
            //TR << da << " " << q.ori.x() << " " << q.angth << std::endl;
            bestrsd = rsd->clone();
            bestsc = da;
          }
        }
      }
      if(bestsc < 9e8) {
        hits.push_back(bestrsd);
      }
    }
  }
};

class KinFunGroupTK : public FunGroupTK {
public:
  KinFunGroupTK( Pose p_in, vector1<Size> pos=vector1<Size>() ) : FunGroupTK(p_in,pos) {}

  virtual void place_n(Vec const & cen0, Vec const & axs0, Vec const & ori0, vector1<core::conformation::ResidueOP> & hits, Real const disth=0.5, Real const dotth=0.9848078) const {
  }
  virtual void place_c(KRSQuery const & q, Residue const & qrsd, vector1<core::conformation::ResidueOP> & hits_out ) const {
#define CYS_CB_HG_DIS      2.911171 // 3.388379 // 2.771698
#define CYS_SG_CB_H        0.7494478 // height of SG "above" CB
#define CYS_HG_SG_PRJLEN   2.088496  //1.818653 // sin(CB-SG-HG)*len(HG-SG)
#define CYS_1oSIN_CB_SG_HG 1.09875   //
#define CYS_HG_SG_X_DROP   0.0998    // h drop due measured, 2.1*tan((90-84.011803)/180*pi)*tan((90-65.522644)/180*pi)
#define CYS_CB_SG_BLEN     1.808803
#define CYS_CB_SG_PERP     1.646237
    if( fabs(q.axs.length()-1.0) > 0.0000001 ) utility_exit_with_message("place_c query axs not kormalized()!");
    Real const r3o2 = sqrt(3.0)/2.0;
    Real const dis2ub = sqr(        CYS_CB_HG_DIS+q.disth );
    Real const dis2lb = sqr(max(0.0,CYS_CB_HG_DIS-q.disth));
    vector1<Size>::const_iterator ip = pos_.begin();
    vector1<Stub> const & stb(stb_.find("CYS")->second);
    vector1<ResidueOP> const & rsdlst(rsd_.find("CYS")->second);
    for(vector1<Stub>::const_iterator is = stb.begin(); is!=stb.end(); ++is,++ip) {
      if(*ip==qrsd.seqpos()) continue;
						core::conformation::ResidueOP rtmp = rsdlst[*ip];
      Real const cbd2 = q.cen.distance_squared(is->v);
      if( dis2ub < cbd2 || cbd2 < dis2lb ) continue;
      Vec  const qcen0(is->global2local(q.cen));
      Real const pdis = sqrt(sqr(q.disth)-sqr(sqrt(cbd2)-CYS_CB_HG_DIS)) * qcen0.length() / CYS_CB_HG_DIS;
						Size const NPOS = (pdis > 0.2) ? 7 : 1;
      Vec  const Y = Vec(0,1,0).cross(qcen0).normalized();
      Vec  const Z =          Y.cross(qcen0).normalized();
						vector1<core::conformation::ResidueOP> hits;
      for(int flp = -1; flp <= 1; flp+=2) {
        Real best_angerr = 9e9, best_ch1=0.0, best_ch2=0.0, mn_ang=9e9, mx_ang=-9e9;
        for(Size ipos = 0; ipos < NPOS; ++ipos) {
          Real y=0,z=0;
          if     (ipos==1) { y =  1.0; z =   0.0; }
          else if(ipos==2) { y =  0.5; z =  r3o2; }
          else if(ipos==3) { y = -0.5; z =  r3o2; }
          else if(ipos==4) { y = -1.0; z =   0.0; }
          else if(ipos==5) { y = -0.5; z = -r3o2; }
          else if(ipos==6) { y =  0.5; z = -r3o2; }
          Vec  const qcen = qcen0 + y*Y*pdis + z*Z*pdis;
          // calc chi2
          Real const lqcen = qcen.length();
          Real const qcenxsphere = qcen.x() * CYS_CB_HG_DIS / lqcen;
          Real const h = ( qcenxsphere - CYS_SG_CB_H) * CYS_1oSIN_CB_SG_HG - CYS_HG_SG_X_DROP;
          if( CYS_HG_SG_PRJLEN < fabs(h) ) continue;
          Real const ch2_0 = acos( h/CYS_HG_SG_PRJLEN ); // chi2 = pi +- ch2_0
          //calc chi1
          Real const lqcenpx = sqrt(qcen.y()*qcen.y()+qcen.z()*qcen.z());
          Real const ch1_0 = (qcen.y()>0)? asin(qcen.z()/lqcenpx) : pi-asin(qcen.z()/lqcenpx);
          Real const hgz = sin(pi-ch2_0)*CYS_HG_SG_PRJLEN;
          Real const hgdzy = lqcenpx * CYS_CB_HG_DIS / lqcen;
          Real const ch1_ofst = asin( hgz/hgdzy );
          Real const ch1 = ch1_0 + flp*ch1_ofst;
          Real const ch2 = pi    + flp*ch2_0   ;
										// check angle
          Vec  const sg  = is->local2global(Vec(CYS_SG_CB_H, cos(ch1)*CYS_CB_SG_PERP, sin(ch1)*CYS_CB_SG_PERP));
										if(!ifc_.clash_check(sg,*ip) ) continue;

										Real const bang = acos( (sg-q.cen).normalized().dot(q.axs) );
        // // make rsd
        // core::conformation::ResidueOP rtmp = core::conformation::ResidueFactory::create_residue(frs_->name_map("CYS"),pose_.residue(*ip),pose_.conformation());
        // rtmp->set_xyz(11, rtmp->xyz(6) + 2.1*(rtmp->xyz(11)-rtmp->xyz(6)).normalized() );
        // rtmp->set_chi(1,degrees(ch1));
        // rtmp->set_chi(2,degrees(ch2));
								// //        ResidueOP hrsd = qrsd.clone();
								// //        xform_rsd_gl2(*is,*rtmp);
								// //        xform_rsd_gl2(*is,*hrsd);
        // //hrsd->set_xyz("HE2", hrsd->xyz("HE2") * CYS_CB_HG_DIS / hrsd->xyz("HE2").length() );
        // Size tmp = 1;
        // ozstream out1("cys_"+lzs(*ip,3)+"_"+str(ipos)+"_"+str(flp+1)+".pdb");
        // core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
        // //core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
        // out1.close();
								// //								utility_exit_with_message(str(degrees(bang)));
										mn_ang = min(bang,mn_ang);
										mx_ang = max(bang,mx_ang);
          Real const angerr = fabs( bang - q.ori.x() );
          if( angerr < best_angerr ) {
            best_angerr = angerr;
            best_ch1 = ch1;
            best_ch2 = ch2;
          }
        }
								if( q.ori.x() < mn_ang-q.angth || mx_ang+q.angth < q.ori.x() ) continue;

        // make rsd
        if(!rtmp) {
										rtmp = core::conformation::ResidueFactory::create_residue(frs_->name_map("CYS"),pose_.residue(*ip),pose_.conformation());
										rtmp->set_xyz(11, rtmp->xyz(6) + 2.1*(rtmp->xyz(11)-rtmp->xyz(6)).normalized() );
								}
        rtmp->set_chi(1,degrees(best_ch1));
        rtmp->set_chi(2,degrees(best_ch2));

								// if( hits.size() == 1 ) {
								// 		Real tmp = fabs( acos((hits[1]->xyz("SG")-q.cen).normalized().dot(q.axs)) - q.ori.x() );
								// 		if( tmp < best_angerr ) hits[1] = rtmp;
								// } else {
								hits.push_back(rtmp->clone());
										//								}
								//if(hits.size()) TR << *ip << " " << hits.size() << std::endl;

        // ResidueOP hrsd = qrsd.clone();
        // xform_rsd_gl2(*is,*rtmp);
        // xform_rsd_gl2(*is,*hrsd);
        // hrsd->set_xyz("HE2", hrsd->xyz("HE2") * CYS_CB_HG_DIS / hrsd->xyz("HE2").length() );
        // Size tmp = 1;
        // ozstream out1("cys"+str(ipos)+"_"+str(flp+1)+".pdb");
        // core::io::pdb::dump_pdb_residue(*rtmp,tmp,out1);
        // core::io::pdb::dump_pdb_residue(*hrsd,tmp,out1);
        // out1.close();

      }
						hits_out.insert(hits_out.end(),hits.begin(),hits.end());
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

  Pose gln,his;
  make_pose_from_sequence(gln,"Q","fa_standard",false);
  make_pose_from_sequence(his,"H","fa_standard",false);
  if(gln.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(gln,1);
  if(gln.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(gln,1);
  if(his.residue(1).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(his,1);
  if(his.residue(1).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(his,1);

  // setup stuff
  ResidueTypeSetCAP frs=ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
  ResidueTypeSetCAP crs=ChemicalManager::get_instance()->residue_type_set( "centroid" );
  ScoreFunctionOP sfstd=get_score_function();
  ScoreFunctionOP sfcen=ScoreFunctionFactory::create_score_function("score3");

  Pose pose;
  pose_from_pdb(pose,option[in::file::s]()[1]);
  pose = alapose(pose);
  FunGroupTK & krs(*(new   KinFunGroupTK(pose)));
  FunGroupTK & brs(*(new BruteFunGroupTK(pose)));
  TR << "start" << std::endl;

		vector1<core::pack::rotamer_set::RotamerSetOP> rots(pose.n_residue());
		for(Size i=1; i<=pose.n_residue(); ++i) {
				core::conformation::Residue rtmp(pose.residue(i));
				pose.replace_residue(i,his.residue(1),true);
				rots[i] = get_rotset(pose,i);
				pose.replace_residue(i,rtmp,false);
		}

  time_t ktime=0,btime=0;
  Size nkh=0,nbh=0;

  for(int i=0;i<1;i++) {
    nkh=0;
    nbh=0;
    // his / cys
    for(Size ir=1;ir<=pose.n_residue();++ir) {
      core::conformation::Residue rtmp(pose.residue(ir));
      pose.replace_residue(ir,his.residue(1),true);
      for(Size irot = 1; irot <= rots[ir]->num_rotamers(); ++irot) {
        pose.set_chi(1,ir,rots[ir]->rotamer(irot)->chi(1));
        pose.set_chi(2,ir,rots[ir]->rotamer(irot)->chi(2));
        bool clash = false;
        for(Size ia = 6; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
          if( !krs.ifc().clash_check(pose.residue(ir).xyz(ia),ir) ) { clash=true; break; }
        }
        if(clash) continue;
        pose.set_dof(core::id::DOF_ID(AtomID(pose.residue(ir).atom_index("HE2"),ir),D ),2.1);

        Vec ne = pose.residue(ir).xyz("NE2");
        Vec he = pose.residue(ir).xyz("HE2");
        Vec axs = (ne-he).normalized();
        Vec cen = ne + 2.1*(-axs);
        Vec ori = Vec(1.911136,0,0);
        KRSQuery q(CEN_ANG,cen,axs,ori,0.5,0.175);
        vector1<core::conformation::ResidueOP> bhits,khits;

        time_t tmp = clock();
        for(int i=0;i<1;++i) krs.place_c( q, pose.residue(ir), khits );
        ktime += clock()-tmp;
        tmp = clock();
        brs.place_c( q, pose.residue(ir), bhits );
        btime += clock()-tmp;
        nkh += khits.size();
        nbh += bhits.size();

        if( khits.size() != bhits.size() ) {
          for(Size ih=1;ih<=khits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(khits[ih]->seqpos()));
            pose.replace_residue(khits[ih]->seqpos(),*khits[ih],true);
            pose.dump_pdb("KIN_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(khits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(khits[ih]->seqpos(),rtmp2,false);
          }
          for(Size ih=1;ih<=bhits.size();++ih) {
            core::conformation::Residue rtmp2(pose.residue(bhits[ih]->seqpos()));
            pose.replace_residue(bhits[ih]->seqpos(),*bhits[ih],false);
            pose.dump_pdb("BRT_"+lzs(ir,3)+"_"+lzs(irot,3)+"_"+lzs(ih,2)+"_"+lzs(bhits[ih]->seqpos(),3)+".pdb");
            pose.replace_residue(bhits[ih]->seqpos(),rtmp2,false);
          }
          /*utility_exit_with_message*/TR << ("kin vs. brute mismatch!!! rsd: "+str(ir)+" rot: "+str(irot)+" brute: "+str(bhits.size())+" kin: "+str(khits.size())) << std::endl;
        }

      }
      pose.replace_residue(ir,rtmp,false);
    }
  }
  TR << "BRT hits: " << nbh << " BRT time: " << btime << std::endl;
  TR << "KIN hits: " << nkh << " KIN time: " << ktime << std::endl;
  TR << "K/B ratio " << (Real)btime / (Real)ktime << std::endl;

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




