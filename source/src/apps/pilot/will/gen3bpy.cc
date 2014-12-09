// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/pack/optimizeH.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SymMinMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>



using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static thread_local basic::Tracer TR( "gen3bpy" );
static core::io::silent::SilentFileData sfd;


inline Real sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
  if( sqdist > 9.*9. ) {
    return 0.0;
  } else if( sqdist < 6.*6. ) {
    return 1.0;
  } else {
    Real dist = sqrt( sqdist );
    return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
  }
}


std::pair<Size,Size> makesplitwork_3bpy(Size total) {
  using namespace basic::options::OptionKeys;
  Size size1 = 1;
  Size size2 = total;
  if( option[willmatch::splitwork].user() ) {
    Size part   = option[willmatch::splitwork]()[1];
    Size nparts = option[willmatch::splitwork]()[2];
    size1 = (part-1)*(Size)(std::ceil(((Real)total)/(Real)nparts))+1;
    size2 = (part  )*(Size)(std::ceil(((Real)total)/(Real)nparts));
    if( option[in::file::s]().size() == 1 ) {
      Real frac1 = ((Real)part-1)/(Real)nparts;
      Real frac2 = ((Real)part  )/(Real)nparts;
      frac1 = 1.0 - sqrt(1.0-frac1);
      frac2 = 1.0 - sqrt(1.0-frac2);
      size1 = (Size)(std::ceil(frac1*(Real)total))+1;
      size2 = (Size)(std::ceil(frac2*(Real)total));
    }
  }
  TR << "SIZE " << size1 << " " << size2 << std::endl;
  return std::pair<Size,Size>(size1,size2);
}



void run() {

  using namespace basic::options::OptionKeys;

  bool USEPOI(false);
  Vec POI(0,0,0);
  if( option[willmatch::poi].user() ) {
    POI.x() = option[willmatch::poi]()[1];
    POI.y() = option[willmatch::poi]()[2];
    POI.z() = option[willmatch::poi]()[3];
    USEPOI = true;
  }

  core::chemical::ResidueTypeSetCAP  rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  Pose bpy1,bpy2;
  core::import_pose::pose_from_pdb(bpy1,*rs,"input/bpy_aa_align_trimer_in.pdb");
  core::import_pose::pose_from_pdb(bpy2,*rs,"input/bpy_aa_align_trimer_in2.pdb");
  for(int ibpy = 1; ibpy <= 2; ++ibpy) {
    Pose & bpy( ibpy==1 ? bpy1 : bpy2 );
    core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
    core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);
    bpy.append_residue_by_jump( *core::conformation::ResidueFactory::create_residue(rs->name_map("VRT")),1,"NN1","ORIG");
    core::kinematics::FoldTree ft = bpy.fold_tree();
    ft.reorder(2);
    bpy.fold_tree(ft);
  }
  // ScoreFunctionOP sf = core::scoring::get_score_function();
  ScoreFunctionOP sf = new core::scoring::symmetry::SymmetricScoreFunction(core::scoring::get_score_function());

  for(Size ifile = 1; ifile <= option[in::file::s]().size(); ++ifile) {
    string infile = utility::file_basename(option[in::file::s]()[ifile]);
    Pose init;
    core::import_pose::pose_from_pdb(init,*rs,option[in::file::s]()[ifile]);
    Size nres = init.n_residue();

    vector1<int> fres = option[willmatch::forbid_residues].user() ? option[willmatch::forbid_residues] : vector1<int>();
    vector1<numeric::xyzTriple<Real> > chis;
    for(Size i = 1; i <= init.n_residue(); ++i) {
			if(init.residue(i).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(init,i);
			if(init.residue(i).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(init,i);
      if( std::find(fres.begin(),fres.end(),i) != fres.end() ) continue;
      for(Real chi1 = -180.0; chi1 < 180.0; chi1+=option[willmatch::chi1_increment]()) {
        for(Real chi2 = -180.0; chi2 < 180.0; chi2+=option[willmatch::chi2_increment]()) {
          chis.push_back(numeric::xyzTriple<Real>(i,chi1,chi2));
        }
      }
    }
    std::pair<Size,Size> split = makesplitwork_3bpy(chis.size());

    for(int ibpy = 1; ibpy <= 2; ++ibpy) {
      Pose & bpy( ibpy==1 ? bpy1 : bpy2 );

      Pose pose = init;
      vector1<Stub> bpystub(pose.n_residue());
      for(Size ir = 1; ir <= init.n_residue(); ++ir) {
				core::conformation::ResidueOP rsd = pose.residue(ir).clone();
        pose.replace_residue(ir,bpy.residue(1),true);
        bpystub[ir] = Stub(pose.xyz(AtomID(2,ir)),pose.xyz(AtomID(1,ir)),pose.xyz(AtomID(3,ir)));
				rsd->chain(1);
        pose.replace_residue(ir,*rsd,false);
      }
			//			pose.dump_pdb("test.pdb");

      Mat const R1 = rotation_matrix_degrees(Vec(0,0,1), 120.0);
      Mat const R2 = rotation_matrix_degrees(Vec(0,0,1),-120.0);

      Size Nbpy = bpy.residue(1).nheavyatoms();
      ImplicitFastClashCheck ifc(pose,option[willmatch::clash_dis]());

      for(Size ichi = max((Size)1,split.first); ichi <= min(chis.size(),split.second); ichi++) {
        Size irsd = (Size)chis[ichi].x();
        if( ichi % (chis.size()/10) == 0 ) TR << "progress: " << irsd << std::endl;
        if( pose.residue(irsd).name3()=="GLY" || pose.residue(irsd).name3()=="PRO" || pose.residue(irsd).name3()=="DPR" ) continue;
        Real chi1 = chis[ichi].y();
        Real chi2 = chis[ichi].z();
        bpy.set_phi  (1,pose.phi  (irsd));
        bpy.set_psi  (1,pose.psi  (irsd));
        bpy.set_omega(1,pose.omega(irsd));
        bpy.set_chi(1,1,chi1);
        bpy.set_chi(2,1,chi2);

        Stub stmp( bpy.xyz(AtomID(2,1)), bpy.xyz(AtomID(1,1)), bpy.xyz(AtomID(3,1)));
        Stub & sp( bpystub[irsd] );
        Mat const spMt = sp.M.transposed();
        stmp.v = stmp.v - stmp.M*spMt*sp.v;
        stmp.M = stmp.M * spMt;
        Stub const & sb(stmp);

        // bpy.dump_pdb("bpy.pdb");
        // xform_pose(pose,sb);
        // pose.dump_pdb("pose.pdb");
        // TR << "cmd.load_cgo( Vec(" << sb.local2global(POI) << ").cgo(),'arst')" << std::endl;
        // utility_exit_with_message("aiorstn");

        if( USEPOI && sb.local2global(POI).length() > option[willmatch::poidis]() ) continue;

        for(Size ia = 6; ia <= Nbpy; ++ia) if(!ifc.clash_check(sb.global2local(    bpy.xyz(AtomID(ia,1))),irsd) ) goto cont1;
        for(Size ia = 1; ia <= Nbpy; ++ia) if(!ifc.clash_check(sb.global2local( R1*bpy.xyz(AtomID(ia,1)))     ) ) goto cont1;
        for(Size ia = 1; ia <= Nbpy; ++ia) if(!ifc.clash_check(sb.global2local( R2*bpy.xyz(AtomID(ia,1)))     ) ) goto cont1;
        goto done1; cont1: continue; done1:

        for(Size ir = 1; ir <= nres; ir++) {
          Size N = (pose.residue(ir).name3()=="GLY") ? 4 : 5;
          for(Size ia = 1; ia <= N; ++ia) {
            if(!ifc.clash_check( sb.global2local( R1*sb.local2global( pose.xyz(AtomID(ia,ir))) )) ) goto cont2; // inefficient!!
            if(!ifc.clash_check( sb.global2local( R2*sb.local2global( pose.xyz(AtomID(ia,ir))) )) ) goto cont2;
          }
        }
        goto done2; cont2: continue; done2:

        Real NB = 0;
        for(Size ir = 1; ir <= nres; ir++) {
          if(!pose.residue(ir).is_protein()) continue;
          Vec const X1 = pose.residue(ir).xyz(5);
          Vec const X1A = sb.global2local( R1*sb.local2global( X1 ) );
          Vec const X1B = sb.global2local( R2*sb.local2global( X1 ) );
          for(Size jr = 1; jr <= nres; jr++) {
            if(!pose.residue(jr).is_protein()) continue;
            Vec const X2 = pose.residue(jr).xyz(5);
            NB += sigmoidish_neighbor(X1A.distance_squared(X2));
            NB += sigmoidish_neighbor(X1B.distance_squared(X2));
          }
        }
        if(NB < option[willmatch::interface_size]()) continue;

        string fname = infile+"_"+lead_zero_string_of(irsd,3)+"_bpy"+lead_zero_string_of(ibpy,1)+"_"+lead_zero_string_of((Size)(chi1+180.0),3)+"_"+lead_zero_string_of((Size)(chi2+180.0),3)+".pdb";
        Pose tmp(pose);
        xform_pose(tmp,sb);
        tmp.replace_residue(irsd,bpy.residue(1),false);
        core::pose::symmetry::make_symmetric_pose( tmp );
				ozstream out(fname);
        tmp.dump_pdb(out);
				Vec viz;
				viz =    sb.local2global(POI); out<<"HETATM"<<I(5,9999)<<' '<<"POI "<<' '<<"POI"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				viz = R1*sb.local2global(POI); out<<"HETATM"<<I(5,9999)<<' '<<"POI "<<' '<<"POI"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				viz = R2*sb.local2global(POI); out<<"HETATM"<<I(5,9999)<<' '<<"POI "<<' '<<"POI"<<' '<<"Z"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				out.close();
        //utility_exit_with_message("clash");

      }
    }
  }
}

int main (int argc, char *argv[]) {

	try {

  // Vec p(-6.456746, 5.922204, -0.982538);
  // Vec d(0.393718,  0.677101,  0.621707);
  // Vec v(0.000000,  0.000000,  0.000000);
  // Vec a(0.998233, -0.003844, -0.059301);
  // Real t = 45;
  // vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
  // for(Size i = 1; i <= X.size(); ++i) {
  //  TR << "X " << i << " " << X[i] << std::endl;
  // }
  // utility_exit_with_message("debug line_cone_intersection");

  devel::init(argc,argv);
  run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}




//
//








