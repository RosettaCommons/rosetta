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


#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
//#include <core/conformation/symmetry/SymDof.hh>
//#include <core/conformation/symmetry/SymmetricConformation.hh>
//#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/docking/util.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
//#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
//#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <protocols/protein_interface_design/movers/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>
#include <core/pack/task/ResfileReader.hh>
//Auto Headers

#include <apps/pilot/will/will_util.ihh>


static THREAD_LOCAL basic::Tracer TR( "mkresfile_i213" );

using core::Size;
using core::Real;
using core::pose::Pose;
using std::string;
using utility::vector1;
using std::endl;
typedef vector1<Size> Sizes;


vector1<Size> get_des_pos(core::pose::Pose & pose_for_design) {
  // Find out which positions are near the inter-subunit interfaces
  // These will be further screened below, then passed to design() and minimize()
  Real const contact_dist = 10.0;
  Real const contact_dist_sq = contact_dist * contact_dist;

  core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function();

  ////////////////////////////////////
  vector1<bool> indy_resis(pose_for_design.n_residue(),false);
  for(Size i = 1; i <= pose_for_design.n_residue()/2; ++i) indy_resis[i] = true;
  vector1<bool> subunit_index(pose_for_design.n_residue());
  for(Size i = 1; i <= pose_for_design.n_residue(); ++i) subunit_index[i] = (6*(i-1))/pose_for_design.n_residue();
  vector1<Size> intra_subs; intra_subs.push_back(1); intra_subs.push_back(2); intra_subs.push_back(3);


  Sizes interface_pos;
  for (Size ir=1; ir<= pose_for_design.n_residue()/2; ir++) {
    std::string atom_i = "";
    if (pose_for_design.residue(ir).name3() == "GLY") {
      atom_i = "CA";
    } else {
      atom_i = "CB";
    }
    for (Size jr=1+pose_for_design.n_residue()/2; jr<= pose_for_design.n_residue(); jr++) {
      std::string atom_j = "";
      if (pose_for_design.residue(jr).name3() == "GLY") {
        atom_j = "CA";
      } else {
        atom_j = "CB";
      }
      if (pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
        interface_pos.push_back(ir);
        break;
      }
    }
  }

  return interface_pos;

}

Size argmax(vector1<Size> x) {
  Size amx=1; for(Size i = 2; i <= x.size(); ++i) if(x[i] > x[amx]) amx = i;
  return amx;
}


void pymol(Vec v, Vec c) {
  static int N = 0;
  std::cout << "veccgofrompoint(Vec("<<10*v.x()<<","<<10*v.y()<<","<<10*v.z()<<"), Vec("<<c.x()<<","<<c.y()<<","<<c.z()<<"),'"<< "cgo"+str(N) <<"' )" << std::endl;
  N++;
}

int
main (int argc, char *argv[])
{

	try {

  devel::init(argc,argv);
  using namespace basic::options;

  for(Size ifile = 1; ifile <= option[OptionKeys::in::file::s]().size(); ++ifile) {
    string fname = option[OptionKeys::in::file::s]()[ifile];
    Pose p;
    core::import_pose::pose_from_file(p,fname, core::import_pose::PDB_file);
    Size nsub = p.n_residue()/6;
    Size nsym = 3;

    utility::io::ozstream out("mutalyze_"+fname+".resfile");
    out << "AUTO\nNATRO\n\nstart\n";
    vector1<Size> resis = get_des_pos(p);
    vector1<Size> subs(6,0);
    for(vector1<Size>::iterator i = resis.begin(); i != resis.end(); ++i) {
      ++subs[ (*i-1)/nsub+1 ];
      out << ObjexxFCL::format::I(4,(*i-1)%nsub+1) << " A PIKAA " << p.residue(*i).name1() << std::endl;
    }
    out.close();

    Size psub = argmax(subs);

    Mat R2 = rotation_matrix_degrees(Vec(0,0,1),120.0);
    Mat R3 = rotation_matrix_degrees(Vec(0,0,1),240.0);

    Vec c0,c1,c2,c3;
    Vec a0,a1,a2,a3;
    Vec c00,c10,c20,c30;
    Vec c01,c11,c21,c31;
    Vec c02,c12,c22,c32;
    Vec b00,b10,b20,b30;
    Vec b01,b11,b21,b31;
    Vec b02,b12,b22,b32;

    c0 = com(p,0*nsub*nsym+1,1*nsub*nsym);
    c1 = com(p,1*nsub*nsym+1,2*nsub*nsym);
    c2 = R2 * c1;
    c3 = R3 * c1;

    a0 = symaxis(p,nsub,3,0);
    a1 = symaxis(p,nsub,3,1);
    a2 = R2*a1;
    a3 = R3*a1;

    c00 = com(p,0*nsub*nsym+1,0*nsub*nsym+nsub);
    c01 = rotation_matrix_degrees(a0,120.0)*(c00-c0)+c0;
    c02 = rotation_matrix_degrees(a0,240.0)*(c00-c0)+c0;
    c10 = com(p,1*nsub*nsym+1,1*nsub*nsym+nsub);
    c11 = rotation_matrix_degrees(a1,120.0)*(c10-c1)+c1;
    c12 = rotation_matrix_degrees(a1,240.0)*(c10-c1)+c1;
    c20 = R2*c10;
    c21 = R2*c11;
    c22 = R2*c12;
    c30 = R3*c10;
    c31 = R3*c11;
    c32 = R3*c12;


    b00 = projperp(a0,c00-c0).normalized();
    b10 = projperp(a1,c10-c1).normalized();
    b20 = projperp(a2,c20-c2).normalized();
    b30 = projperp(a3,c30-c3).normalized();

    b01 = projperp(a0,c01-c0).normalized();
    b11 = projperp(a1,c11-c1).normalized();
    b21 = projperp(a2,c21-c2).normalized();
    b31 = projperp(a3,c31-c3).normalized();

    b02 = projperp(a0,c02-c0).normalized();
    b12 = projperp(a1,c12-c1).normalized();
    b22 = projperp(a2,c22-c2).normalized();
    b32 = projperp(a3,c32-c3).normalized();

    // pymol(a0,c0);
    // pymol(b00,c0);
    // pymol(b10,c1);
    // pymol(b20,c2);
    // pymol(b30,c3);

    // pymol(a1,c1);
    // pymol(b01,c0);
    // pymol(b11,c1);
    // pymol(b21,c2);
    // pymol(b31,c3);

    // pymol(a2,c2);
    // pymol(b02,c0);
    // pymol(b12,c1);
    // pymol(b22,c2);
    // pymol(b32,c3);

    utility::io::ozstream sout("mutalyze_"+fname+".sym");
    //std::ostream & sout(std::cout);
    sout << std::fixed << std::setprecision(6);
    sout << "symmetry_name i213" << endl;
    sout << "subunits 12" << endl;
    sout << "number_of_interfaces 11" << endl;
    sout << "E = 12*S00+6*(S00:S01)+6*(S00:S02)+6*(S00:S10)+6*(S00:S11)+6*(S00:S12)+6*(S00:S20)+6*(S00:S21)+6*(S00:S22)+6*(S00:S30)+6*(S00:S31)+6*(S00:S32) " << endl;
    sout << "anchor_residue 1" << endl;
    sout << "virtual_coordinates_start" << endl;
    sout << "xyz C00 " <<a0.x()<<","<<a0.y()<<","<<a0.z() << " "<<b00.x()<<","<<b00.y()<<","<<b00.z() << " "<<c0.x()<<","<<c0.y()<<","<<c0.z() << endl;
    sout << "#" << endl;
    sout << "xyz P00 " <<a0.x()<<","<<a0.y()<<","<<a0.z() << " "<<b00.x()<<","<<b00.y()<<","<<b00.z() << " "<<c0.x()<<","<<c0.y()<<","<<c0.z() << endl;
    sout << "#" << endl;
    sout << "xyz S00 " <<a0.x()<<","<<a0.y()<<","<<a0.z() << " "<<b00.x()<<","<<b00.y()<<","<<b00.z() << " "<<c0.x()<<","<<c0.y()<<","<<c0.z() << endl;
    sout << "xyz S01 " <<a0.x()<<","<<a0.y()<<","<<a0.z() << " "<<b01.x()<<","<<b01.y()<<","<<b01.z() << " "<<c0.x()<<","<<c0.y()<<","<<c0.z() << endl;
    sout << "xyz S02 " <<a0.x()<<","<<a0.y()<<","<<a0.z() << " "<<b02.x()<<","<<b02.y()<<","<<b02.z() << " "<<c0.x()<<","<<c0.y()<<","<<c0.z() << endl;
    sout << "#" << endl;
    sout << "xyz S10 " <<a1.x()<<","<<a1.y()<<","<<a1.z() << " "<<b10.x()<<","<<b10.y()<<","<<b10.z() << " "<<c1.x()<<","<<c1.y()<<","<<c1.z() << endl;
    sout << "xyz S11 " <<a1.x()<<","<<a1.y()<<","<<a1.z() << " "<<b11.x()<<","<<b11.y()<<","<<b11.z() << " "<<c1.x()<<","<<c1.y()<<","<<c1.z() << endl;
    sout << "xyz S12 " <<a1.x()<<","<<a1.y()<<","<<a1.z() << " "<<b12.x()<<","<<b12.y()<<","<<b12.z() << " "<<c1.x()<<","<<c1.y()<<","<<c1.z() << endl;
    sout << "#" << endl;
    sout << "xyz S20 " <<a2.x()<<","<<a2.y()<<","<<a2.z() << " "<<b20.x()<<","<<b20.y()<<","<<b20.z() << " "<<c2.x()<<","<<c2.y()<<","<<c2.z() << endl;
    sout << "xyz S21 " <<a2.x()<<","<<a2.y()<<","<<a2.z() << " "<<b21.x()<<","<<b21.y()<<","<<b21.z() << " "<<c2.x()<<","<<c2.y()<<","<<c2.z() << endl;
    sout << "xyz S22 " <<a2.x()<<","<<a2.y()<<","<<a2.z() << " "<<b22.x()<<","<<b22.y()<<","<<b22.z() << " "<<c2.x()<<","<<c2.y()<<","<<c2.z() << endl;
    sout << "#" << endl;
    sout << "xyz S30 " <<a3.x()<<","<<a3.y()<<","<<a3.z() << " "<<b30.x()<<","<<b30.y()<<","<<b30.z() << " "<<c3.x()<<","<<c3.y()<<","<<c3.z() << endl;
    sout << "xyz S31 " <<a3.x()<<","<<a3.y()<<","<<a3.z() << " "<<b31.x()<<","<<b31.y()<<","<<b31.z() << " "<<c3.x()<<","<<c3.y()<<","<<c3.z() << endl;
    sout << "xyz S32 " <<a3.x()<<","<<a3.y()<<","<<a3.z() << " "<<b32.x()<<","<<b32.y()<<","<<b32.z() << " "<<c3.x()<<","<<c3.y()<<","<<c3.z() << endl;
    sout << "#" << endl;
    sout << "virtual_coordinates_stop" << endl;
    sout << "connect_virtual JP   C00  P00" << endl;
    sout << "connect_virtual J0   P00  S00" << endl;
    sout << "connect_virtual J1   P00  S01" << endl;
    sout << "connect_virtual J2   P00  S02" << endl;
    sout << "connect_virtual J3   C00  S10" << endl;
    sout << "connect_virtual J4   C00  S11" << endl;
    sout << "connect_virtual J5   C00  S12" << endl;
    sout << "connect_virtual J6   C00  S20" << endl;
    sout << "connect_virtual J7   C00  S21" << endl;
    sout << "connect_virtual J8   C00  S22" << endl;
    sout << "connect_virtual J9   C00  S30" << endl;
    sout << "connect_virtual J10  C00  S31" << endl;
    sout << "connect_virtual J11  C00  S32" << endl;
    sout << "connect_virtual JS0  S00 SUBUNIT" << endl;
    sout << "connect_virtual JS1  S01 SUBUNIT" << endl;
    sout << "connect_virtual JS2  S02 SUBUNIT" << endl;
    sout << "connect_virtual JS3  S10 SUBUNIT" << endl;
    sout << "connect_virtual JS4  S11 SUBUNIT" << endl;
    sout << "connect_virtual JS5  S12 SUBUNIT" << endl;
    sout << "connect_virtual JS6  S20 SUBUNIT" << endl;
    sout << "connect_virtual JS7  S21 SUBUNIT" << endl;
    sout << "connect_virtual JS8  S22 SUBUNIT" << endl;
    sout << "connect_virtual JS9  S30 SUBUNIT" << endl;
    sout << "connect_virtual JS10 S31 SUBUNIT" << endl;
    sout << "connect_virtual JS11 S32 SUBUNIT" << endl;
    sout << "set_dof JP x angle_x" << endl;
    sout << "set_jump_group JGS JS0 JS1 JS2 JS3 JS4 JS5 JS6 JS7 JS8 JS9 JS10 JS11" << endl;
    sout.close();

    Pose tmp;
    core::kinematics::FoldTree f(nsub);
    vector1<Size> pos(nsub);
    for(Size i = 1; i <= nsub; ++i) pos[i] = i;
    core::pose::create_subpose(p,pos,f,tmp);
    tmp.dump_pdb("mutalyze_"+fname+"_sub1.pdb");

  }


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

