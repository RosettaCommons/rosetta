// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#include <core/id/AtomID_Map.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <fstream>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

// AUTO-REMOVED #include <protocols/moves/Mover.hh>
#include <sstream>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>
#include <utility/vector1.hh>

static basic::Tracer TR("symdock_hybrid");


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

bool strip_termini(core::pose::Pose & pose) {
  core::scoring::dssp::Dssp dssp(pose);
  dssp.insert_ss_into_pose(pose);
  int ROOT;
  for(ROOT=1; ROOT <= (int)pose.n_residue(); ++ROOT) if(pose.secstruct(ROOT)!='L') break;
  if(ROOT>=(int)pose.n_residue() || pose.n_residue() < 20) return false;
  core::kinematics::FoldTree ft;
  ObjexxFCL::FArray2D<int> jp(2,1); jp(1,1) = ROOT; jp(2,1) = ROOT+1;
  ObjexxFCL::FArray1D<int> cp(1,ROOT);
  ft.tree_from_jumps_and_cuts( (int)pose.n_residue(), (int)1, jp, cp, ROOT );
  pose.fold_tree(ft);
  while(pose.secstruct(       1        )=='L') pose.delete_polymer_residue(       1        );
  while(pose.secstruct(pose.n_residue())=='L') pose.delete_polymer_residue(pose.n_residue());
  for(Size i = 1; i <= pose.n_residue(); ++i) {
    if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
    if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
  }
  return true;
}

int
main (int argc, char *argv[]){

	try {

  using namespace core;
  using basic::options::option;
  using namespace basic::options::OptionKeys;

  devel::init( argc, argv );
  pose::Pose cc3;
  import_pose::pose_from_pdb(cc3,"input/cc3.pdb"); strip_termini(cc3);
  cc3.set_xyz(AtomID(cc3.residue(1).atom_index("H"),1),Vec(0,0,0));
  cc3.set_xyz(AtomID(cc3.residue(2).atom_index("H"),2),Vec(0,0,1));

  for(Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile) {
    string fname = option[in::file::s]()[ifile];
    TR << fname << std::endl;
    pose::Pose pose;
    import_pose::pose_from_pdb(pose,fname);
    if(!strip_termini(pose)) continue;
    if(pose.n_residue() < 40) continue;

    Vec trax = Vec(0,0,1);

    for(Size nc = 0; nc <= 1; nc++) {
      for(Size of = 0; of <= 7; of++) {

        Pose h(cc3);
        core::id::AtomID_Map<AtomID> amap;
        core::pose::initialize_atomid_map(amap,h,core::id::BOGUS_ATOM_ID);
        // get sup atom map
        Size h1, h2, p1/*, p2*/;  // p2 is unused ~Labonte
        if(nc) {
          p1=1+of; /*p2=7+of;*/ h1=   h.n_residue()-8; h2=h.n_residue();
        } else {
          h1=1; h2=7; p1=pose.n_residue()-8-of; /*p2=pose.n_residue()-of;*/
        }
        for(Size i = h1; i <= h2; ++i) {
          Size pr = p1+(i-h1);
          //for(Size j = 1; j <= min((Size)5,min(pose.residue(pr).nheavyatoms(),h.residue(i).nheavyatoms())); ++j) amap[ AtomID(j,i) ] = AtomID(j,pr);
          amap[ AtomID(1,i) ] = AtomID(1,pr);
          amap[ AtomID(2,i) ] = AtomID(2,pr);
          amap[ AtomID(3,i) ] = AtomID(3,pr);
        }
        core::scoring::superimpose_pose(h,pose,amap);
        Real rms=0;
        for(Size i = h1; i <= h2; ++i) {
          Size pr = p1+(i-h1);
          rms += pose.residue(pr).xyz(1).distance_squared(h.residue(i).xyz(1));
          rms += pose.residue(pr).xyz(2).distance_squared(h.residue(i).xyz(2));
          rms += pose.residue(pr).xyz(3).distance_squared(h.residue(i).xyz(3));
        }
        rms = sqrt( rms/27.0 );
        //TR << fname << " " << rms << std::endl;
        if(rms > 1.0) continue;

        //pose.dump_pdb("pose.pdb");
        //h.dump_pdb("h.pdb");
        //utility_exit_with_message("klfg");

        Vec ccax =  h.residue(2).xyz("H") - h.residue(1).xyz("H");
        Vec c = Vec(0,0,0)-h.residue(1).xyz("H");
        Vec x = trax.cross(ccax);
        Real d = fabs(c.dot(x)) / x.length();
        if( d > 5.0 ) { /*TR << "d fail " << d << std::endl;*/ continue; }
        Real ang = numeric::angle_degrees(ccax,Vec(0,0,0),Vec(0,0,1));
        if( fabs(ang-54.7356563997) > 7.0 && fabs((180.0-ang)-54.7356563997) > 7.0 ) { /*TR << "ang fail " << ang << std::endl*/; continue; }

        Vec c3f1 = Vec(0,0,0);
        Vec c2f1 = h.residue(1).xyz("H");
        Vec isct = - trax * ((c.cross(ccax).dot(trax.cross(ccax)))/((trax.cross(ccax).length()*trax.cross(ccax).length())));

        if( isct.distance(c2f1) < 25.0 || isct.distance(c3f1) < 25.0 || isct.distance(c2f1)/isct.distance(c3f1) > 2.0 || isct.distance(c2f1) / isct.distance(c3f1) < 0.5 ) {
          //TR << "fail on "+fname+": too close to c3/c2: " << isct.distance(c2f1) << " " << isct.distance(c3f1) <<  std::endl;
          continue;
        }

        std::string ofn = utility::file_basename(fname)+"cc3_"+(nc?"N_":"C_")+ObjexxFCL::string_of(of)+"_";
        TR << "HIT! " << ofn << " " << isct << std::endl;

        // pose.dump_pdb("test.pdb");
        // h.dump_pdb("test2.pdb");
        // Utility::io::ozstream out("test3.pdb");
        // Vec viz(0,0,0);
        // viz = viz; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        // viz = isct; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        // viz = c2f1; out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        // viz = h.residue(2).xyz("H"); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        // out.close();
        // utility_exit_with_message("test.pdb");

				{
					Vec  ax = ((c3f1-isct).normalized() + (c2f1-isct).normalized()).normalized();
					Mat  ra = numeric::rotation_matrix_degrees(ax,180.0);
					h.dump_pdb("testh.pdb");
					pose.dump_pdb("test.pdb");
					rot_pose(h,ra,isct);
					rot_pose(pose,ra,isct);
					h.dump_pdb("2testh.pdb");
					pose.dump_pdb("2test.pdb");
				}

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

        h   .dump_pdb(ofn+"_coil.pdb");
        symm.dump_pdb(ofn+"_main.pdb");

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
