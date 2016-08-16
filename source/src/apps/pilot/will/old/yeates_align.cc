// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#include <core/id/AtomID_Map.hh>
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <fstream>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

#include <sstream>
#include <utility/io/ozstream.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


static THREAD_LOCAL basic::Tracer TR( "yeates_align" );


using core::pose::Pose;
using std::string;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using numeric::xyzVector;
using numeric::min;
using namespace ObjexxFCL::format;

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;

int
main (int argc, char *argv[]){

	try {

  using namespace core;
  using basic::options::option;
  using namespace basic::options::OptionKeys;

  devel::init( argc, argv );
  pose::Pose ccpp,ccap;
  import_pose::pose_from_file(ccpp,"input/ccpp.pdb"); remove_termini(ccpp, core::import_pose::PDB_file);
  import_pose::pose_from_file(ccap,"input/ccap.pdb"); remove_termini(ccap, core::import_pose::PDB_file);
  ccpp.set_xyz(AtomID(ccpp.residue(1).atom_index("H"),1),Vec(0,0,0));
  ccpp.set_xyz(AtomID(ccpp.residue(2).atom_index("H"),2),Vec(0,0,1));
  ccap.set_xyz(AtomID(ccap.residue(1).atom_index("H"),1),Vec(0,0,0));
  ccap.set_xyz(AtomID(ccap.residue(2).atom_index("H"),2),Vec(0,0,1));

  for(Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile) {
    string fname = option[in::file::s]()[ifile];
    //TR << fname << std::endl;
    pose::Pose pose;
    import_pose::pose_from_file(pose,fname); remove_termini(pose, core::import_pose::PDB_file);

    core::scoring::dssp::Dssp dssp(pose);
    Size nstart=0,nstop=0,cstart=0,cstop=0;
    for(Size i = 1; i <= pose.n_residue(); ++i) {
      if( !nstart && dssp.get_dssp_secstruct(i) == 'H' ) {
        nstart = i;
      } else if(  nstart && dssp.get_dssp_secstruct(i) != 'H' ) {
        nstop = i;
        break;
      }
    }
    for(Size i = pose.n_residue(); i >= 1; --i) {
      if( !cstart && dssp.get_dssp_secstruct(i) == 'H' ) {
        cstart = i;
      } else if(  cstart && dssp.get_dssp_secstruct(i) != 'H' ) {
        cstop = i;
        break;
      }
    }
    //TR << fname << " " << nstart << " " << nstop << " " << pose.n_residue()-cstart << " " << pose.n_residue()-cstop << std::endl;

    Vec trax = Vec(0,0,1);

    for(Size ap = 0; ap <= 1; ap++) {
      for(Size nc = 0; nc <= 1; nc++) {

        Pose h( ap?ccap:ccpp );
        core::id::AtomID_Map<AtomID> amap;
        core::pose::initialize_atomid_map(amap,h);
        // get sup atom map
        if( nc && nstart < 10 && nstop-nstart > 6 ) {
          Size naln = min(nstop-nstart+1,(Size)9);
          for(Size i = 1; i <= naln; ++i) {
            Size pr = nstart+i-1;
            Size hr = h.n_residue()-naln+i;
            for(Size j = 1; j <= min((Size)5,min(pose.residue(pr).nheavyatoms(),h.residue(hr).nheavyatoms())); ++j) amap[ AtomID(j,hr) ] = AtomID(j,pr);
          }
        } else if( !nc && cstart < 10 && cstart-cstop > 6 ) {
          Size naln = min(cstart-cstop+1,(Size)9);
          for(Size i = 1; i <= naln; ++i) {
            Size pr = cstop+i-1;
            Size hr = i;
            for(Size j = 1; j <= min((Size)5,min(pose.residue(pr).nheavyatoms(),h.residue(hr).nheavyatoms())); ++j) amap[ AtomID(j,hr) ] = AtomID(j,pr);
          }
        } else {
          continue;
        }

        core::scoring::superimpose_pose(h,pose,amap);
        Vec ccax =  h.residue(2).xyz("H") - h.residue(1).xyz("H");
        Vec c = Vec(0,0,0)-h.residue(1).xyz("H");
        Vec x = trax.cross(ccax);
        Real d = fabs(c.dot(x)) / x.length();
        if( d > 8.0 ) { /*TR << "d fail " << d << std::endl;*/ continue; }
        Real ang = numeric::angle_degrees(ccax,Vec(0,0,0),Vec(0,0,1));
        if( fabs(ang-54.7356563997) > 10.0 && fabs((180.0-ang)-54.7356563997) > 15.0 ) { /*TR << "ang fail " << ang << std::endl*/; continue; }

        Pose t2(pose); rot_pose(t2,Vec(0,0,1),120.0);
        Pose t3(pose); rot_pose(t3,Vec(0,0,1),240.0);
        Pose d2(pose); rot_pose(d2,h.residue(2).xyz("H"),180.0);

        Vec c3f1 = Vec(0,0,0);
        Vec c2f1 = h.residue(1).xyz("H");
        Vec isct = - trax * ((c.cross(ccax).dot(trax.cross(ccax)))/((trax.cross(ccax).length()*trax.cross(ccax).length())));

        if( isct.distance(c2f1) < 25.0 || isct.distance(c3f1) < 25.0 || isct.distance(c2f1)/isct.distance(c3f1) > 2.0 || isct.distance(c2f1) / isct.distance(c3f1) < 0.5 ) {
          TR << "fail: too close to c3/c2: " << isct.distance(c2f1) << " " << isct.distance(c3f1) <<  std::endl;
          continue;
        }


        std::string ofn = utility::file_basename(fname)+(ap?"_AP_":"_PR_")+(nc?"N_":"C_");
        TR << "HIT! " << ofn << " " << isct << std::endl;

        pose.dump_pdb("test.pdb");
        h.dump_pdb("test2.pdb");
        utility::io::ozstream out("test3.pdb");
        Vec viz(0,0,0);
        viz = viz; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        viz = isct; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        viz = c2f1; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        viz = h.residue(2).xyz("H"); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        out.close();
        utility_exit_with_message("test.pdb");


        Vec sax1 = Vec( 0.816496579408716,0, 0.57735027133783);
        Vec sax2 = Vec(-0.816496579408716,0,-0.57735027133783);
        Real da = numeric::dihedral(c2f1,Vec(0,0,0),Vec(0,0,1),(ang>90)?sax1:sax2);
        Mat Rsymm = numeric::rotation_matrix(Vec(0,0,1),-da);

        Pose symm = pose;
        trans_pose(symm,-isct);
        rot_pose  (symm,Rsymm);
        trans_pose(h,-isct);
        rot_pose  (h,Rsymm);
        core::pose::symmetry::make_symmetric_pose(symm);
        core::pose::symmetry::make_symmetric_pose(h);

        h   .dump_pdb(ofn+"_0.pdb");
        symm.dump_pdb(ofn+"_1.pdb");

        //utility_exit_with_message("dsgf");
      }

    }
  }

  return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
