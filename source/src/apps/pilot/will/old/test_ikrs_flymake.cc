// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/willmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include <apps/pilot/will/will_util.ihh>

using core::Real;
using core::Size;
using core::pose::Pose;
using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using std::string;
using utility::vector1;
using numeric::min;
using core::import_pose::pose_from_file;
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

static THREAD_LOCAL basic::Tracer TR( "test_ikrs" );


vector1<Reals> vecs2vv(Vecs const & v) {
  vector1<Reals> vv;
  for(Vecs::const_iterator i = v.begin(); i != v.end(); ++i) {
    Reals r(3);
    r[1] = i->x();
    r[2] = i->y();
    r[3] = i->z();
    vv.push_back(r);
  }
  return vv;
}
vector1<Vec> vv2vecs(vector1<Reals> const & vv) {
  vector1<Vec> r;
  for(vector1<Reals>::const_iterator i = vv.begin(); i != vv.end(); ++i) {
    r.push_back(Vec( (*i)[1], (*i)[2], (*i)[3] ));
  }
  return r;
}

struct Hit : utility::pointer::ReferenceCount {
  Hit(core::conformation::Residue const & r1, core::conformation::Residue const & r2) : rsd1(r1),rsd2(r2) {}
  core::conformation::Residue const rsd1,rsd2;
  Vec cen,axs,ori;
};
typedef utility::pointer::owning_ptr<Hit> HitOP;

void ik_arg_glu_frnt(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits) {
  using namespace basic::options::OptionKeys;
  using namespace core::id;

  utility::vector1<Size> pivots (3), order (3);
  utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
  int nsol=0;
  // for eliminating identical solutions
  utility::vector1<utility::vector1<Real> > Q0 (3);
  utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

  // N CA C N CA CB CG CD NE
  Vecs atoms(18);
  atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
  atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
  atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
  atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
  atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
  atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
  atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
  atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
  atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
  atoms[10] = pose.residue(rsd2  ).xyz( "CG" );
  atoms[11] = pose.residue(rsd2  ).xyz( "CB" );
  atoms[12] = pose.residue(rsd2  ).xyz( "CA" );
  atoms[13] = pose.residue(rsd2  ).xyz( "N"  );
  atoms[14] = pose.residue(rsd2-1).xyz( "C"  );
  atoms[15] = pose.residue(rsd2-1).xyz( "CA" );
  atoms[16] = pose.residue(rsd2-1).xyz( "N"  );
  atoms[17] = pose.residue(rsd2-2).xyz( "C"  );
  atoms[18] = pose.residue(rsd2-2).xyz( "CA" );

  order [1]=1; order [2]=2; order [3]=3;
  pivots[1]=5, pivots[2]=8, pivots[3]=11;

  using namespace numeric::kinematic_closure;
  chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
  core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

  db_len[ 9] = 6.9;
  db_ang[ 9] = numeric::angle_degrees(pose.residue(rsd1).xyz("CD"),pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"));
  db_ang[10] = numeric::angle_degrees(pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"),pose.residue(rsd2).xyz("CD"));
  Real phitgt = dt_ang[4];
  Size count = 0;
  for(Size idh = 3; idh <= 360; idh += 10) {
    dt_ang[9] = (Real)idh;
    for(Size ichi2 = 3; ichi2 <= 360; ichi2 += 10) {
      dt_ang[6] = (Real)ichi2;
      bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
      if(nsol==0) break;
      for(int isol = 1; isol <= nsol; isol++) {
        Real phidiff = phitgt-t_ang[isol][4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
        if( fabs(phidiff) > 10.0 ) continue;

        utility::vector1<utility::vector1<core::Real> > vv_atm_out;
        numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
        Vecs apos = vv2vecs(vv_atm_out);

        bool clash = false;
        for( Size i = 1; i <= atoms.size(); ++i ) {
          for( Size j = i+3; j <= atoms.size(); ++j ) {
            if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
          }
          if(clash) break;
        }
        if(clash) break;

        pose.set_chi(1,rsd1,t_ang[isol][ 5]);
        pose.set_chi(2,rsd1,t_ang[isol][ 6]);
        pose.set_chi(3,rsd1,t_ang[isol][ 7]);
        pose.set_chi(4,rsd1,t_ang[isol][ 8]);
        pose.set_chi(3,rsd2,t_ang[isol][ 9]); // this is tricky...
        pose.set_chi(2,rsd2,t_ang[isol][10]);
        pose.set_chi(1,rsd2,t_ang[isol][11]);
        Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
        Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
        if( dun1 > 16.0 ) continue;
        if( dun2 > 12.0 ) continue;

        Vec CA = pose.xyz(AtomID(2,rsd1));
        Mat M = rotation_matrix_degrees( CA-pose.xyz(AtomID(1,rsd1)) , -phidiff );
        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M * (pose.xyz(AtomID(i,rsd1))-CA) + CA );

        for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
        for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

        count++;
        if(clash/*||count!=20*/){if(clash)count--;for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M.transposed()*(pose.xyz(AtomID(i,rsd1))-CA)+CA);continue;}

        HitOP hitop = new Hit( pose.residue(rsd1), pose.residue(rsd2));
        Hit & hit(*hitop);
        hit.cen = hit.rsd1.xyz("CZ");
        hit.axs = ((hit.rsd1.xyz("NH2")+hit.rsd1.xyz("NE"))/2.0 - hit.cen).normalized();
        hit.ori = hit.axs.cross(hit.rsd1.xyz("NE")-hit.cen).normalized();
        hits.push_back(hitop);

        // {
        //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_frnt_glu_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+"_res.pdb");
        //   Size ano = 0;
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd2),ano,out);
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd1),ano,out);
        //   out.close();
        // }
        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M.transposed() * (pose.xyz(AtomID(i,rsd1))-CA) + CA );
      }
    }
    if(nsol==0) break;
  }
}

//*************************************************************************************************************************************************************
//*************************************************************************************************************************************************************
//*************************************************************************************************************************************************************
void ik_arg_glu_side(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits) {
  using namespace basic::options::OptionKeys;
  using namespace core::id;

  utility::vector1<Size> pivots (3), order (3);
  utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
  int nsol=0;
  // for eliminating identical solutions
  utility::vector1<utility::vector1<Real> > Q0 (3);
  utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

  // N CA C N CA CB CG CD NE
  Vecs atoms(18);
  atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
  atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
  atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
  atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
  atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
  atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
  atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
  atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
  atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
  atoms[10] = pose.residue(rsd1  ).xyz( "CZ" );

  atoms[11] = pose.residue(rsd2  ).xyz( "CG" );
  atoms[12] = pose.residue(rsd2  ).xyz( "CB" );
  atoms[13] = pose.residue(rsd2  ).xyz( "CA" );
  atoms[14] = pose.residue(rsd2  ).xyz( "N"  );
  atoms[15] = pose.residue(rsd2-1).xyz( "C"  );
  atoms[16] = pose.residue(rsd2-1).xyz( "CA" );
  atoms[17] = pose.residue(rsd2-1).xyz( "N"  );
  atoms[18] = pose.residue(rsd2-2).xyz( "C"  );

  order [1]=1; order [2]=2; order [3]=3;
  pivots[1]=5, pivots[2]=8, pivots[3]=11;

  using namespace numeric::kinematic_closure;
  chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
  core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

  dt_ang[ 9] = 180.0;
  db_len[10] = 5.6;
  db_ang[10] = 180.0 - numeric::angle_degrees(pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"),pose.residue(rsd1).xyz("NH1"));
  db_ang[11] =         numeric::angle_degrees(pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"),pose.residue(rsd2).xyz("CD"));
  Real phitgt = dt_ang[4];
  Size count = 0;
  for(Size idh = 3; idh <= 360; idh += 10) {
    dt_ang[12] = (Real)idh;
    for(Size ichi2 = 3; ichi2 <= 360; ichi2 += 10) {
      dt_ang[6] = (Real)ichi2;
      bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
      if(nsol==0) break;
      for(int isol = 1; isol <= nsol; isol++) {
        Real phidiff = phitgt-t_ang[isol][4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
        if( fabs(phidiff) > 10.0 ) continue;

        utility::vector1<utility::vector1<core::Real> > vv_atm_out;
        numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
        Vecs apos = vv2vecs(vv_atm_out);

        bool clash = false;
        for( Size i = 1; i <= atoms.size(); ++i ) {
          for( Size j = i+3; j <= atoms.size(); ++j ) {
            if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
          }
          if(clash) break;
        }
        if(clash) break;

        pose.set_chi(1,rsd1,t_ang[isol][ 5]);
        pose.set_chi(2,rsd1,t_ang[isol][ 6]);
        pose.set_chi(3,rsd1,t_ang[isol][ 7]);
        pose.set_chi(4,rsd1,t_ang[isol][ 8]);
        pose.set_chi(3,rsd2,t_ang[isol][10]); // this is tricky...
        pose.set_chi(2,rsd2,t_ang[isol][11]);
        // don't set chi1 rsd2 -- should stay as orig because we've explicitly moved dt_ang[12]
        Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
        Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
        if( dun1 > 16.0 ) continue;
        if( dun2 > 12.0 ) continue;

        Vec CA = pose.xyz(AtomID(2,rsd1));
        Mat M = rotation_matrix_degrees( CA-pose.xyz(AtomID(1,rsd1)) , -phidiff );
        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M * (pose.xyz(AtomID(i,rsd1))-CA) + CA );

        for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
        for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

        count++;
        if(clash){//||count!=10){
          if(clash)count--;
          for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M.transposed()*(pose.xyz(AtomID(i,rsd1))-CA)+CA);
          continue;
        }

        HitOP hitop = new Hit( pose.residue(rsd1), pose.residue(rsd2));
        Hit & hit(*hitop);
        hit.cen = hit.rsd1.xyz("CZ");
        hit.axs = ((hit.rsd1.xyz("NH2")+hit.rsd1.xyz("NH1"))/2.0 - hit.cen).normalized();
        hit.ori = hit.axs.cross(hit.rsd1.xyz("NH1")-hit.cen).normalized();
        hits.push_back(hitop);

        // {
        //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_side_glu_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+"_res.pdb");
        //   Size ano = 0;
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd2),ano,out);
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd1),ano,out);
        //   out.close();
        // }

        // {
        //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_side_glu_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+".pdb");
        //   out<<"ATOM  "<<I(5, 1)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 1].x())<<F(8,3,apos[ 1].y())<<F(8,3,apos[ 1].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 2)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 2].x())<<F(8,3,apos[ 2].y())<<F(8,3,apos[ 2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 3)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 3].x())<<F(8,3,apos[ 3].y())<<F(8,3,apos[ 3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out<<"ATOM  "<<I(5, 4)<<' '<<" N  "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 4].x())<<F(8,3,apos[ 4].y())<<F(8,3,apos[ 4].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 5)<<' '<<" CA "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 5].x())<<F(8,3,apos[ 5].y())<<F(8,3,apos[ 5].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 6)<<' '<<" CB "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 6].x())<<F(8,3,apos[ 6].y())<<F(8,3,apos[ 6].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 7)<<' '<<" CG "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 7].x())<<F(8,3,apos[ 7].y())<<F(8,3,apos[ 7].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 8)<<' '<<" CD "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 8].x())<<F(8,3,apos[ 8].y())<<F(8,3,apos[ 8].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 9)<<' '<<" NE "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 9].x())<<F(8,3,apos[ 9].y())<<F(8,3,apos[ 9].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,10)<<' '<<" CZ "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[10].x())<<F(8,3,apos[10].y())<<F(8,3,apos[10]. z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out<<"ATOM  "<<I(5,17)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[17].x())<<F(8,3,apos[17].y())<<F(8,3,apos[17].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,16)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[16].x())<<F(8,3,apos[16].y())<<F(8,3,apos[16].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,15)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[15].x())<<F(8,3,apos[15].y())<<F(8,3,apos[15].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out<<"ATOM  "<<I(5,14)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[14].x())<<F(8,3,apos[14].y())<<F(8,3,apos[14].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,13)<<' '<<" CA "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[13].x())<<F(8,3,apos[13].y())<<F(8,3,apos[13].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,12)<<' '<<" CB "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[12].x())<<F(8,3,apos[12].y())<<F(8,3,apos[12].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,11)<<' '<<" CG "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[11].x())<<F(8,3,apos[11].y())<<F(8,3,apos[11].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out.close();
        // }
        //        utility_exit_with_message("test");

        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M.transposed() * (pose.xyz(AtomID(i,rsd1))-CA) + CA );
      }
    }
    if(nsol==0) break;
  }
}


//**********************************************************************************************************************************************************8
//**********************************************************************************************************************************************************8
//**********************************************************************************************************************************************************8
void ik_arg_asp_frnt(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits) {
  using namespace basic::options::OptionKeys;
  using namespace core::id;

  utility::vector1<Size> pivots (3), order (3);
  utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
  int nsol=0;
  // for eliminating identical solutions
  utility::vector1<utility::vector1<Real> > Q0 (3);
  utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

  // N CA C N CA CB CG CD NE
  Vecs atoms(15);
  atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
  atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
  atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
  atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
  atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
  atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
  atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
  atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
  atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
  atoms[10] = pose.residue(rsd2  ).xyz( "CB" );
  atoms[11] = pose.residue(rsd2  ).xyz( "CA" );
  atoms[12] = pose.residue(rsd2  ).xyz( "N"  );
  atoms[13] = pose.residue(rsd2-1).xyz( "C"  );
  atoms[14] = pose.residue(rsd2-1).xyz( "CA" );
  atoms[15] = pose.residue(rsd2-1).xyz( "N"  );

  order [1]=1; order [2]=2; order [3]=3;
  pivots[1]=5, pivots[2]=8, pivots[3]=11;

  using namespace numeric::kinematic_closure;
  chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
  core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

  db_len[ 9] = 6.9;
  db_ang[ 9] = numeric::angle_degrees(pose.residue(rsd1).xyz("CD"),pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"));
  db_ang[10] = numeric::angle_degrees(pose.residue(rsd2).xyz("CA"),pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"));

  Real phitgt = dt_ang[ 4];
  Real psitgt = dt_ang[11];
  Size count = 0;
  for(Size idh = 3; idh <= 360; idh += 10) {
    dt_ang[ 9] = (Real)idh;
    for(Size ichi2 = 3; ichi2 <= 360; ichi2 += 10) {
      dt_ang[6] = (Real)ichi2;
      bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
      if(nsol==0) break;
      for(int isol = 1; isol <= nsol; isol++) {
        Real phidiff = phitgt-t_ang[isol][ 4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
        Real ph2diff = psitgt-t_ang[isol][11]; while(ph2diff < -180.0) ph2diff+=360.0; while(ph2diff > 180.0) ph2diff-=360.0;
        if( fabs(phidiff) > 13.0 ) continue;
        if( fabs(ph2diff) > 13.0 ) continue;

        utility::vector1<utility::vector1<core::Real> > vv_atm_out;
        numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
        Vecs apos = vv2vecs(vv_atm_out);

        bool clash = false;
        for( Size i = 1; i <= atoms.size(); ++i ) {
          for( Size j = i+3; j <= atoms.size(); ++j ) {
            if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
          }
          if(clash) break;
        }
        if(clash) break;

        //for(Size i = 1; i <= atoms.size(); ++i) TR << "ICHI2 " << ichi2 << " SOL " << isol << " CHI " << i << " " << t_ang[isol][i] << std::endl;
        pose.set_chi(1,rsd1,t_ang[isol][ 5]);
        pose.set_chi(2,rsd1,t_ang[isol][ 6]);
        pose.set_chi(3,rsd1,t_ang[isol][ 7]);
        pose.set_chi(4,rsd1,t_ang[isol][ 8]);
        pose.set_chi(2,rsd2,t_ang[isol][ 9]); // this is tricky...
        pose.set_chi(1,rsd2,t_ang[isol][10]);
        Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
        Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
        if( dun1 > 16.0 ) continue;
        if( dun2 > 12.0 ) continue;
        Vec CA1 = pose.xyz(AtomID(2,rsd1));
				Mat M1 = rotation_matrix_degrees( CA1-pose.xyz(AtomID(1,rsd1)) , -phidiff );
        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1 * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
        Vec CA2 = pose.xyz(AtomID(2,rsd2));
        Mat M2 = rotation_matrix_degrees( CA2-pose.xyz(AtomID(1,rsd2)) , -ph2diff );
        for(Size i = 5; i <= pose.residue_type(rsd2).natoms(); ++i) pose.set_xyz( AtomID(i,rsd2), M2 * (pose.xyz(AtomID(i,rsd2))-CA2) + CA2 );

        for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
        for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

        count++;
        if(clash){//||count!=1){
          if(clash) count--;
          for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M1.transposed()*(pose.xyz(AtomID(i,rsd1))-CA1)+CA1);
          for(Size i=5; i<=pose.residue_type(rsd2).natoms();++i) pose.set_xyz(AtomID(i,rsd2),M2.transposed()*(pose.xyz(AtomID(i,rsd2))-CA2)+CA2);
          continue;
        }

        HitOP hitop = new Hit( pose.residue(rsd1), pose.residue(rsd2));
        Hit & hit(*hitop);
        hit.cen = hit.rsd1.xyz("CZ");
        hit.axs = ((hit.rsd1.xyz("NH2")+hit.rsd1.xyz("NE"))/2.0 - hit.cen).normalized();
        hit.ori = hit.axs.cross(hit.rsd1.xyz("NE")-hit.cen).normalized();
        hits.push_back(hitop);

        // {
        //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_frnt_asp_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+"_res.pdb");
        //   Size ano = 0;
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd2),ano,out);
        //   core::io::pdb::dump_pdb_residue(pose.residue(rsd1),ano,out);
        //   out.close();
        // }
        // {
        //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_frnt_asp_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+".pdb");
        //   out<<"ATOM  "<<I(5, 1)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 1].x())<<F(8,3,apos[ 1].y())<<F(8,3,apos[ 1].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 2)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 2].x())<<F(8,3,apos[ 2].y())<<F(8,3,apos[ 2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 3)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 3].x())<<F(8,3,apos[ 3].y())<<F(8,3,apos[ 3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out<<"ATOM  "<<I(5, 4)<<' '<<" N  "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 4].x())<<F(8,3,apos[ 4].y())<<F(8,3,apos[ 4].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 5)<<' '<<" CA "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 5].x())<<F(8,3,apos[ 5].y())<<F(8,3,apos[ 5].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 6)<<' '<<" CB "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 6].x())<<F(8,3,apos[ 6].y())<<F(8,3,apos[ 6].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 7)<<' '<<" CG "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 7].x())<<F(8,3,apos[ 7].y())<<F(8,3,apos[ 7].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 8)<<' '<<" CD "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 8].x())<<F(8,3,apos[ 8].y())<<F(8,3,apos[ 8].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5, 9)<<' '<<" NE "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 9].x())<<F(8,3,apos[ 9].y())<<F(8,3,apos[ 9].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out<<"ATOM  "<<I(5,15)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[15].x())<<F(8,3,apos[15].y())<<F(8,3,apos[15].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,14)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[14].x())<<F(8,3,apos[14].y())<<F(8,3,apos[14].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,13)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[13].x())<<F(8,3,apos[13].y())<<F(8,3,apos[13].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,12)<<' '<<" N  "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[12].x())<<F(8,3,apos[12].y())<<F(8,3,apos[12].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,11)<<' '<<" CA "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[11].x())<<F(8,3,apos[11].y())<<F(8,3,apos[11].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        //   out<<"ATOM  "<<I(5,10)<<' '<<" CB "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[10].x())<<F(8,3,apos[10].y())<<F(8,3,apos[10].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

        //   out.close();
        // }
        // utility_exit_with_message("test");
        for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1.transposed() * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
        for(Size i = 5; i <= pose.residue_type(rsd2).natoms(); ++i) pose.set_xyz( AtomID(i,rsd2), M2.transposed() * (pose.xyz(AtomID(i,rsd2))-CA2) + CA2 );
      }
    }
    if(nsol==0) break;
  }
}

//**********************************************************************************************************************************************************8
//**********************************************************************************************************************************************************8
//**********************************************************************************************************************************************************8
void ik_arg_asp_side(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits) {
  using namespace basic::options::OptionKeys;
  using namespace core::id;

  utility::vector1<Size> pivots (3), order (3);
  utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
  int nsol=0;
  utility::vector1<utility::vector1<Real> > Q0 (3);
  utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

  Vecs atoms(15);
  atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
  atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
  atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
  atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
  atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
  atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
  atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
  atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
  atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
  atoms[10] = pose.residue(rsd1  ).xyz( "CZ" );
  atoms[11] = pose.residue(rsd2  ).xyz( "CB" );
  atoms[12] = pose.residue(rsd2  ).xyz( "CA" );
  atoms[13] = pose.residue(rsd2  ).xyz( "N"  );
  atoms[14] = pose.residue(rsd2-1).xyz( "C"  );
  atoms[15] = pose.residue(rsd2-1).xyz( "CA" );

  order [1]=1; order [2]=2; order [3]=3;
  pivots[1]=5, pivots[2]=8, pivots[3]=11;
  numeric::kinematic_closure::chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);
  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
  core::pack::rotamers::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
  core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

  dt_ang[ 9] = 180.0;
  db_len[10] = 5.6;
  db_ang[10] = 180.0 - numeric::angle_degrees(pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"),pose.residue(rsd1).xyz("NH1"));
  db_ang[11] =         numeric::angle_degrees(pose.residue(rsd2).xyz("CA"),pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"));

  Real phitgt = dt_ang[ 4];
  Size count = 0;
  for(Size ichi2 = 3; ichi2 <= 360; ichi2 += 10) {
    dt_ang[6] = (Real)ichi2;
    numeric::kinematic_closure::bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
    if(nsol==0) break;
    for(int isol = 1; isol <= nsol; isol++) {
      Real phidiff = phitgt-t_ang[isol][ 4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
      if( fabs(phidiff) > 13.0 ) continue;

      utility::vector1<utility::vector1<core::Real> > vv_atm_out;
      numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
      Vecs apos = vv2vecs(vv_atm_out);

      bool clash = false;
      for( Size i = 1; i <= atoms.size(); ++i ) {
        for( Size j = i+3; j <= atoms.size(); ++j ) {
          if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
        }
        if(clash) break;
      }
      if(clash) break;

      pose.set_chi(1,rsd1,t_ang[isol][ 5]);
      pose.set_chi(2,rsd1,t_ang[isol][ 6]);
      pose.set_chi(3,rsd1,t_ang[isol][ 7]);
      pose.set_chi(4,rsd1,t_ang[isol][ 8]);
      pose.set_chi(2,rsd2,t_ang[isol][10]); // this is tricky...
      pose.set_chi(1,rsd2,t_ang[isol][11]);
      Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
      Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
      if( dun1 > 16.0 ) continue;
      if( dun2 > 12.0 ) continue;
      Vec CA1 = pose.xyz(AtomID(2,rsd1));
      Mat M1 = rotation_matrix_degrees( CA1-pose.xyz(AtomID(1,rsd1)) , -phidiff );
      for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1 * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );

      for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
      for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

      count++;
      if(clash){//||count!=1){
        if(clash) count--;
        for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M1.transposed()*(pose.xyz(AtomID(i,rsd1))-CA1)+CA1);
        continue;
      }

      HitOP hitop = new Hit( pose.residue(rsd1), pose.residue(rsd2));
      Hit & hit(*hitop);
      hit.cen = hit.rsd1.xyz("CZ");
      hit.axs = ((hit.rsd1.xyz("NH2")+hit.rsd1.xyz("NH1"))/2.0 - hit.cen).normalized();
      hit.ori = hit.axs.cross(hit.rsd1.xyz("NH1")-hit.cen).normalized();
      hits.push_back(hitop);

      // {
      //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_side_asp_"+lzs(ichi2,3)+"_"+lzs(isol,2)+"_res.pdb");
      //   Size ano = 0;
      //   core::io::pdb::dump_pdb_residue(pose.residue(rsd2),ano,out);
      //   core::io::pdb::dump_pdb_residue(pose.residue(rsd1),ano,out);
      //   out.close();
      // }
      // {
      //   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_side_asp_"+lzs(ichi2,3)+"_"+lzs(idh,3)+"_"+lzs(isol,2)+".pdb");
      //   out<<"ATOM  "<<I(5, 1)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 1].x())<<F(8,3,apos[ 1].y())<<F(8,3,apos[ 1].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 2)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 2].x())<<F(8,3,apos[ 2].y())<<F(8,3,apos[ 2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 3)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 3].x())<<F(8,3,apos[ 3].y())<<F(8,3,apos[ 3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

      //   out<<"ATOM  "<<I(5, 4)<<' '<<" N  "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 4].x())<<F(8,3,apos[ 4].y())<<F(8,3,apos[ 4].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 5)<<' '<<" CA "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 5].x())<<F(8,3,apos[ 5].y())<<F(8,3,apos[ 5].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 6)<<' '<<" CB "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 6].x())<<F(8,3,apos[ 6].y())<<F(8,3,apos[ 6].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 7)<<' '<<" CG "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 7].x())<<F(8,3,apos[ 7].y())<<F(8,3,apos[ 7].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 8)<<' '<<" CD "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 8].x())<<F(8,3,apos[ 8].y())<<F(8,3,apos[ 8].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5, 9)<<' '<<" NE "<<' '<<"ARG"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 9].x())<<F(8,3,apos[ 9].y())<<F(8,3,apos[ 9].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

      //   out<<"ATOM  "<<I(5,15)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[15].x())<<F(8,3,apos[15].y())<<F(8,3,apos[15].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5,14)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[14].x())<<F(8,3,apos[14].y())<<F(8,3,apos[14].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5,13)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,5)<<"    "<<F(8,3,apos[13].x())<<F(8,3,apos[13].y())<<F(8,3,apos[13].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5,12)<<' '<<" N  "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[12].x())<<F(8,3,apos[12].y())<<F(8,3,apos[12].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5,11)<<' '<<" CA "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[11].x())<<F(8,3,apos[11].y())<<F(8,3,apos[11].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      //   out<<"ATOM  "<<I(5,10)<<' '<<" CB "<<' '<<"GLU"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,apos[10].x())<<F(8,3,apos[10].y())<<F(8,3,apos[10].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

      //   out.close();
      // }
      //utility_exit_with_message("test");
      for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1.transposed() * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
    }
    if(nsol==0) break;
  }
}

void repack(Pose & arg) {
  ScoreFunctionOP sf = core::scoring::get_score_function();
  core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(arg);
  task->restrict_to_repacking();
  protocols::simple_moves::PackRotamersMover repack( sf, task );
  repack.apply(arg);
}

int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);


  Pose pose,arg,asp,glu,lys;
  core::import_pose::pose_from_file(pose,"input/2vdf_nohet_1.pdb", core::import_pose::PDB_file);
  for(Size i = 1; i <= pose.size(); ++i) {
    if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
    if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
  }

  ImplicitFastClashCheck ifc(pose,2.5);

  vector1<Real> sasa; { core::id::AtomID_Map<Real> atom_sasa; core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 5.0, false ); }

  core::chemical::ResidueTypeSetCAP frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
  core::pose::make_pose_from_sequence(arg,"R",*frs,false); repack(arg);
  core::pose::make_pose_from_sequence(asp,"D",*frs,false); repack(asp);
  core::pose::make_pose_from_sequence(glu,"E",*frs,false); repack(glu);
  core::pose::make_pose_from_sequence(lys,"K",*frs,false); repack(lys);

  vector1<HitOP> hits;
  for(Size ir = 3; ir <= pose.size()-2; ++ir) {
    if(sasa[ir] > 0) continue;
    pose.replace_residue(ir,arg.residue(1),true);
    TR << ir << " " << hits.size() << std::endl;
    for(Size jr = 3; jr <= pose.size()-2; ++jr) {
      if(ir==jr) continue;
      if(sasa[jr] > 0) continue;
      pose.replace_residue(jr,asp.residue(1),true);
			ik_arg_asp_frnt(pose,ir,jr,ifc,hits);
			ik_arg_asp_side(pose,ir,jr,ifc,hits);
      pose.replace_residue(jr,glu.residue(1),true);
			ik_arg_glu_frnt(pose,ir,jr,ifc,hits);
			ik_arg_glu_side(pose,ir,jr,ifc,hits);
	// if(hits.size()){
	// 	Hit & h(*hits[hits.size()]);
  //   ozstream out("test.pdb");
  //   Size ano = 0;
  //   core::io::pdb::dump_pdb_residue(h.rsd2,ano,out);
  //   core::io::pdb::dump_pdb_residue(h.rsd1,ano,out);
	// 	{Vec V=h.cen      ;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
	// 	{Vec V=h.axs+h.cen;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
	// 	{Vec V=h.ori+h.cen;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
  //   out.close();
	// 	utility_exit_with_message("asl;dfj");
  // }
    }
  }

  TR << "TOTAL: " << hits.size() << std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
