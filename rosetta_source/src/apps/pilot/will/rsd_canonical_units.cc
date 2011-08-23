#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/gpu.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/orbitals/OrbitalType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/RT.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include "will_util.hh"
#include "mynamespaces.hh"

static basic::Tracer TR("rsd_canonical_units");

int main(int argc, char *argv[]) {
  using namespace core::chemical;

  core::init(argc,argv);
  cout<<std::setprecision(8)<<std::fixed;

  std::map<float,vector1<string> > crdgrp;

  for(int aa = 1; aa <= num_canonical_aas; ++aa) {
    //if( aa != aa_asp ) continue;
    Pose p;
    std::string seq; seq += oneletter_code_from_aa((AA)aa);
    make_pose_from_sequence(p,seq,"fa_standard",false);
    if(p.residue(1).nheavyatoms() < 6) continue;
    remove_termini(p);
    ResidueType const & t(p.residue(1).type());

    vector1<AtomIndices> const & catom = t.chi_atoms();
    for(Size i = 1; i <= t.nchi(); ++i) p.set_chi(i,1,0.0);

    Size p1=1,p2=2,p3=5;
    core::kinematics::Stub stb = core::kinematics::Stub(safe_xyz(p,p3,1),safe_xyz(p,p2,1),safe_xyz(p,p1,1));
    xform_pose_rev(p,stb);
    p.dump_pdb(seq+"0.pdb");

    //cout<<"// coord info "<<seq<<endl;
    std::ostringstream sout;
    for(Size lcc=1; lcc <= t.nchi(); ++lcc) {
      AtomIndices const & sa = catom[lcc];
      stb = core::kinematics::Stub(safe_xyz(p,sa[3],1),safe_xyz(p,sa[2],1),safe_xyz(p,sa[1],1));
      xform_pose_rev(p,stb);
      p.dump_pdb(seq+str(lcc)+".pdb");

      Size hvycnt = 0;
      for(Size ai = 1; ai <= t.nheavyatoms(); ++ai) {
        if(lcc!=t.last_controlling_chi(ai)) continue;

        // heavy atom
        {
          Vec X = safe_xyz(p,ai,1);
          std::ostringstream otmp;
          otmp<<"#define "
              <<LJ(15,seq+""+str(t.last_controlling_chi(ai))+"_"+strip(t.atom_name(ai)))
              <<" vec("<<F(11,8,X.x())<<","<<F(11,8,X.y())<<","<<F(11,8,X.z())<<")"<<std::endl;
          cout<<otmp.str();
          crdgrp[X.x()+X.y()+X.z()].push_back( otmp.str() );
          hvycnt++;
        }

        // polar hydrogens
        for(Size hi = t.attached_H_begin(ai); hi <= t.attached_H_end(ai); ++hi) {
          if(25==t.atom_type_index(hi)) continue;
          {
            Vec XH = safe_xyz(p,hi,1);
            std::ostringstream otmp;
            otmp<<"#define "
                <<LJ(15,seq+""+str(t.last_controlling_chi(ai))+"_"+strip(t.atom_name(ai))+"_"+strip(t.atom_name(hi)))
                <<" vec("<<F(11,8,XH.x())<<","<<F(11,8,XH.y())<<","<<F(11,8,XH.z())<<") " <<std::endl;
            cout<<otmp.str();
            crdgrp[XH.x()+XH.y()+XH.z()].push_back(otmp.str());
          }
        }

        // orbitals
        vector1<Size> orb = t.bonded_orbitals(ai);
        if(t.atom_type(ai).is_acceptor()){
          for(Size oi = 1; oi <= orb.size(); ++oi) {
            using namespace orbitals;
            orbital_type_enum ot = t.orbital_type(orb[oi]).orbital_enum();
            if(ot==N_p_sp2||ot==O_p_sp2||ot==O_p_sp3|ot==S_p_sp3) {
              Vec XO = orb_xyz(p,orb[oi],1);
              std::ostringstream otmp;
              otmp<<"#define "
                  <<LJ(15,seq+""+str(t.last_controlling_chi(ai))+"_"+strip(t.atom_name(ai))+"_orb"+t.orbital_name(orb[oi]))
                  <<" vec("<<F(11,8,XO.x())<<","<<F(11,8,XO.y())<<","<<F(11,8,XO.z())<<")"<<std::endl;
              cout<<otmp.str();
              crdgrp[XO.x()+XO.y()+XO.z()].push_back(otmp.str());
              //cout<<"  orbl "<<t.orbital_name(orb[oi])<<" "<<t.orbital_type(orb[oi]).orbital_enum()<<std::endl;
            }
          }
        }
      }
      if(1==lcc || 0==hvycnt) continue;
      core::kinematics::Stub toprev(safe_xyz(p,p3,1),safe_xyz(p,p2,1),safe_xyz(p,p1,1));
      p1 = sa[1];  p2 = sa[2];  p3 = sa[3];
      Vec ct = toprev.v;
      Mat cr = toprev.M;
      Real cw;
      Vec ca = rotation_axis(cr,cw);
      if( ca.normalized().dot(Vec(0,0,-1).normalized()) < 0.99999 ) {
        cout<<ca.x()<<" "<<ca.y()<<" "<<ca.z()<<" "<<cw<<" <- SHOULD BE  0, 0,-1!"<<std::endl;
        utility_exit_with_message("chi "+str(lcc)+" to "+str(lcc-1)+" axis not 0,0,-1!!!");
      }
      if( fabs(ct.y()) > 0.00001 || fabs(ct.z()) > 0.00001 || ct.x() > -1.0 ) {
        cout<<ct.x()<<" "<<ct.y()<<" "<<ct.z()<<"<- SHOULD BE -X, 0, 0!"<<std::endl;
        utility_exit_with_message("chi "+str(lcc)+" to "+str(lcc-1)+" trans  not -R,0,0!!!");
      }
      sout<<"#define "+seq+""+str(lcc)+""+str(lcc-1)+"_rot001_COS "<<F(12,9,cos(-cw))+"f"<<std::endl;
      sout<<"#define "+seq+""+str(lcc)+""+str(lcc-1)+"_rot001_SIN "<<F(12,9,sin(-cw))+"f"<<std::endl;
      sout<<"#define "+seq+""+str(lcc)+""+str(lcc-1)+"_mov100     "<<F(12,9, ct.x() )+"f"<<std::endl;
    }
    cout << sout.str();

  }

  cout << endl << "COORD MULTIPLES" << endl;
  for(std::map<float,vector1<string> >::iterator mi = crdgrp.begin(); mi != crdgrp.end(); ++mi) {
    if(mi->second.size()==1) continue;
    for(Size ci = 1; ci <= mi->second.size(); ++ci) {
      cout << mi->second[ci];
    }
    cout << endl;
  }


}
