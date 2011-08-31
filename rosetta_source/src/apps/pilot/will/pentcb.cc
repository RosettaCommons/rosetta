// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
//#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/symmetry/SymMinMover.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
#include "apps/pilot/will/will_util.hh"
#include "mynamespaces.hh"

#include <sys/stat.h>


#define CONTACT_D2 20.25
#define CONTACT_TH 12
#define NSS 15872

typedef numeric::xyzVector<Real> Vecf;
typedef numeric::xyzMatrix<Real> Matf;



using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;
using core::pose::Pose;
using core::conformation::ResidueOP;

static basic::Tracer TR("pentcb");
static core::io::silent::SilentFileData sfd;



double
sicfast(
vector1<Vecf> & pa,
  vector1<Vecf> & pb,
  vector1<Vecf> & cba,
  vector1<Vecf> & cbb,
//Vec ori,
  int & cbcount,
  bool debug = false
  ) {
double BIN = 1.8;

// get points, rotated ro ori is 0,0,1, might already be done
// Mat rot = Mat::identity();
// if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
// else if( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
// if( rot != Mat::identity() ) {
//   for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
//   for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
// }

// get bounds for plane hashes
double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
}
for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
 ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
 }
 xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
 ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);


 int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
 int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

 // TR << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
 // TR << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

 // insert points into hashes
 int const xsize = xub-xlb+1;
 int const ysize = yub-ylb+1;
 ObjexxFCL::FArray2D<Vecf> ha(xsize,ysize,Vecf(0,0,-9e9)),hb(xsize,ysize,Vecf(0,0,9e9));
 for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
   // int const ix = min(xsize,max(1,(int)ceil(ia->x()/BIN)-xlb));
   // int const iy = min(ysize,max(1,(int)ceil(ia->y()/BIN)-ylb));
   int const ix = (int)ceil(ia->x()/BIN)-xlb;
   int const iy = (int)ceil(ia->y()/BIN)-ylb;
   if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
   if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
 }
 for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
   // int const ix = min(xsize,max(1,(int)ceil(ib->x()/BIN)-xlb));
   // int const iy = min(ysize,max(1,(int)ceil(ib->y()/BIN)-ylb));
   int const ix = (int)ceil(ib->x()/BIN)-xlb;
   int const iy = (int)ceil(ib->y()/BIN)-ylb;
   if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
   if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
 }

 // check hashes for min dis
 int imna=0,jmna=0,imnb=0,jmnb=0;
 double mindis = 9e9;
 for(int i = 1; i <= xsize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
   for(int j = 1; j <= ysize; ++j) {
     for(int k = -2; k <= 2; ++k) {
       if(i+k < 1 || i+k > xsize) continue;
       for(int l = -2; l <= 2; ++l) {
         if(j+l < 1 || j+l > ysize) continue;
         double const xa = ha(i  ,j  ).x();
         double const ya = ha(i  ,j  ).y();
         double const xb = hb(i+k,j+l).x();
         double const yb = hb(i+k,j+l).y();
         double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);

         if( d2 < BIN*BIN*4.0 ) {
           double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(BIN*BIN*4.0-d2);
           if( dz < mindis ) {
             mindis = dz;
             imna = i;
             jmna = j;
             imnb = i+k;
             jmnb = j+l;
           }
         }
       }
     }
   }
 }

 // {
 //  utility::io::ozstream out("cba.pdb");
 //  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
 //    Vec viz = (*ia) + (mindis*ori);
 //    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //  }
 //  out.close();
 // }
 // {
 //  utility::io::ozstream out("cbb.pdb");
 //  for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
 //    Vec viz = (*ib);
 //    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //  }
 //  out.close();
 // }

 cbcount = 0;
 // utility::io::ozstream out("cb8.pdb");
 // TR << "CB0 " << cbcount << std::endl;
 for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
   for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
     if( ib->distance_squared( (*ia) + (mindis*Vecf(0,0,1)) ) < CONTACT_D2 ) {
       cbcount++;
       // Vec viz = (*ia) + (mindis*ori);
       // out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
       // viz = *ib;
       // out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
     }
   }
 }
 // out.close();
 // TR << "CB1 " << cbcount << std::endl;

 // // rotate points back -- needed iff pa/pb come by reference
 //rot = rot.transposed();
 // if( rot != Mat::identity() ) {
 //  for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
 //  for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
 // }

 // uncomment this to get hashes in local space
 // rot = Mat::identity();
 // ori = Vec(0,0,1);

 // if(debug){
 //   {
 //     utility::io::ozstream out("hasha.pdb");
 //     for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
 //       for(int j = 2; j <= ysize-1; ++j) {
 //         Vec viz = rot*ha(i,j) + mindis*ori;
 //         if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
 //         out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //       }
 //     }
 //     Vec viz = rot*ha(imna,jmna) + mindis*ori;
 //     out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"B"<<I(4,100+jmna)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //     out.close();
 //   }
 //   {
 //     utility::io::ozstream out("hashb.pdb");
 //     for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
 //       for(int j = 2; j <= ysize-1; ++j) {
 //         Vec viz = rot*hb(i,j);
 //         if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
 //         out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"C"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //       }
 //     }
 //     Vec viz = rot*hb(imnb,jmnb);
 //     out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"C"<<I(4,100+jmnb)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
 //     out.close();
 //   }
 // }

 return mindis;
}

// double sicfast(
//                Pose const & a,
//                Pose const & b,
//                Vec ori_in,
//                int & cbcount
//                ) {
//   // get points, rotated ro ori is 0,0,1
//   vector1<Vec> pa,pb;
//   vector1<Vec> cba,cbb;
//   Vec ori = ori_in.normalized();
//   Mat rot = Mat::identity();
//   if     ( ori.dot(Vec(0,0,1)) < -0.999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
//   else if( ori.dot(Vec(0,0,1)) <  0.999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
//   for(int i = 1; i <= (int)a.n_residue(); ++i) {
//     if(a.residue(i).name3()=="GLY") {
//       cba.push_back(rot*Vec(a.residue(i).xyz(2)));
//       for(int j = 1; j <= 4; ++j) pa.push_back(rot*Vec(a.residue(i).xyz(j)));
//     } else {
//       cba.push_back(rot*Vec(a.residue(i).xyz(5)));
//       for(int j = 1; j <= 5; ++j) pa.push_back(rot*Vec(a.residue(i).xyz(j)));
//     }
//   }
//   for(int i = 1; i <= (int)b.n_residue(); ++i) {
//     if(b.residue(i).name3()=="GLY") {
//       cbb.push_back(rot*Vec(b.residue(i).xyz(2)));
//       for(int j = 1; j <= 4; ++j) pb.push_back(rot*Vec(b.residue(i).xyz(j)));
//     } else {
//       cbb.push_back(rot*Vec(b.residue(i).xyz(5)));
//       for(int j = 1; j <= 5; ++j) pb.push_back(rot*Vec(b.residue(i).xyz(j)));
//     }
//   }
//   return sicfast( pa, pb, cba, cbb, Vec(0,0,1), cbcount );
// }

struct Hit {
  int iss,irt,cbc;
	core::kinematics::Stub s1,s2;
  Hit(int is, int ir, int cb) : iss(is),irt(ir),cbc(cb) {}
};
bool cmp(Hit i,Hit j) { return i.cbc > j.cbc; }


void dock(Pose const init, std::string const & fn, vector1<xyzVector<double> > const & ssamp) {
  vector1<double> asamp; asamp.reserve(115);
  // asamp.push_back(180.00000);asamp.push_back(179.00354);asamp.push_back(178.00700);asamp.push_back(177.01032);asamp.push_back(176.01340);asamp.push_back(175.01619);
  // asamp.push_back(174.01859);asamp.push_back(173.02055);asamp.push_back(172.02197);asamp.push_back(171.02279);asamp.push_back(170.02292);asamp.push_back(169.02228);
  // asamp.push_back(168.02081);asamp.push_back(167.01842);asamp.push_back(166.01503);asamp.push_back(165.01057);asamp.push_back(164.00494);asamp.push_back(162.99807);
  // asamp.push_back(161.98987);asamp.push_back(160.98027);asamp.push_back(159.96918);asamp.push_back(158.95651);asamp.push_back(157.94217);asamp.push_back(156.92608);
  // asamp.push_back(155.90815);asamp.push_back(154.88828);asamp.push_back(153.86639);asamp.push_back(152.84238);asamp.push_back(151.81616);asamp.push_back(150.78762);
  // asamp.push_back(149.75667);asamp.push_back(148.72321);asamp.push_back(147.68713);asamp.push_back(146.64833);asamp.push_back(145.60670);asamp.push_back(144.56214);
  // asamp.push_back(143.51452);asamp.push_back(142.46373);asamp.push_back(141.40967);asamp.push_back(140.35219);asamp.push_back(139.29119);asamp.push_back(138.22653);
  // asamp.push_back(137.15808);asamp.push_back(136.08570);asamp.push_back(135.00927);asamp.push_back(133.92863);asamp.push_back(132.84364);asamp.push_back(131.75415);
  // asamp.push_back(130.66000);asamp.push_back(129.56103);asamp.push_back(128.45708);asamp.push_back(127.34796);asamp.push_back(126.23351);asamp.push_back(125.11353);
  // asamp.push_back(123.98784);asamp.push_back(122.85624);asamp.push_back(121.71852);asamp.push_back(120.57447);asamp.push_back(119.42386);asamp.push_back(118.26646);
  // asamp.push_back(117.10204);asamp.push_back(115.93033);asamp.push_back(114.75107);asamp.push_back(113.56400);asamp.push_back(112.36882);asamp.push_back(111.16522);
  // asamp.push_back(109.95291);asamp.push_back(108.73153);asamp.push_back(107.50075);asamp.push_back(106.26020);asamp.push_back(105.00950);asamp.push_back(103.74823);
  // asamp.push_back(102.47598);asamp.push_back(101.19227);asamp.push_back(99.89665) ;asamp.push_back(98.58859) ;asamp.push_back(97.26755) ;asamp.push_back(95.93297 );
  // asamp.push_back(94.58422) ;asamp.push_back(93.22066) ;asamp.push_back(91.84158) ;asamp.push_back(90.44624) ;asamp.push_back(89.03383) ;asamp.push_back(87.60349 );
  // asamp.push_back(86.15429) ;asamp.push_back(84.68521) ;asamp.push_back(83.19517) ;asamp.push_back(81.68297) ;asamp.push_back(80.14733) ;asamp.push_back(78.58682 );
  // asamp.push_back(76.99990) ;asamp.push_back(75.38486) ;asamp.push_back(73.73980) ;asamp.push_back(72.06262) ;asamp.push_back(70.35100) ;asamp.push_back(68.60231 );
  // asamp.push_back(66.81358) ;asamp.push_back(64.98148) ;asamp.push_back(63.10218) ;asamp.push_back(61.17128) ;asamp.push_back(59.18369) ;asamp.push_back(57.13344 );
  // asamp.push_back(55.01348) ;asamp.push_back(52.81535) ;asamp.push_back(50.52882) ;asamp.push_back(48.14123) ;asamp.push_back(45.63666) ;asamp.push_back(42.99462 );
  // asamp.push_back(40.18793) ;asamp.push_back(37.17927) ;asamp.push_back(33.91485) ;asamp.push_back(30.31207) ;asamp.push_back(26.23179) ;asamp.push_back(21.40253 );
  // asamp.push_back(15.12286) ;asamp.push_back(0.00000 ) ;
	for(Real i = 0; i < 180; ++i) asamp.push_back(i);

  vector1<Vecf> bb0tmp,cb0tmp;
  for(int ir = 1; ir <= init.n_residue(); ++ir) {
    for(int ia = 1; ia <= ((init.residue(ir).has("CB"))?5:4); ++ia) {
      bb0tmp.push_back(init.xyz(AtomID(ia,ir)));
    }
    if(init.secstruct(ir)=='H') {
      if(init.residue(ir).has("CB")) cb0tmp.push_back(init.xyz(AtomID(5,ir)));
      else                           cb0tmp.push_back(init.xyz(AtomID(4,ir)));
    }
  }
  vector1<Vecf> const bb0(bb0tmp);
  vector1<Vecf> const cb0(cb0tmp);

  Matf const R2 = rotation_matrix_degrees(Vec(1,0,0),180.0);
  Matf const R3 = rotation_matrix_degrees(Vec(1,0,0),120.0);
  Matf const R4 = rotation_matrix_degrees(Vec(1,0,0), 90.0);
  Matf const R5 = rotation_matrix_degrees(Vec(1,0,0), 72.0);
  Matf const R6 = rotation_matrix_degrees(Vec(1,0,0), 60.0);

  vector1<Hit> hits;

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
  for(int iss = 1; iss <= ssamp.size(); ++iss) {
    if(iss%1000==0) TR << iss << " of " << NSS << std::endl;
    Vec axs = ssamp[iss];
    for(int irt = 1; irt <= asamp.size(); ++irt) {
      //if(axs.x() < 0.0) continue;
	//Real const u = uniform();
      Real const rot = asamp[irt];
      Matf const R = rotation_matrix_degrees(axs,rot);

      vector1<Vecf> bb1 = bb0;
      vector1<Vecf> cb1 = cb0;
      for(vector1<Vecf>::iterator i = bb1.begin(); i != bb1.end(); ++i) *i = R*(*i);
      for(vector1<Vecf>::iterator i = cb1.begin(); i != cb1.end(); ++i) *i = R*(*i);
      vector1<Vecf> bb2 = bb1;
      vector1<Vecf> cb2 = cb1;
      for(vector1<Vecf>::iterator i = bb2.begin(); i != bb2.end(); ++i) *i = R5*(*i);
      for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = R5*(*i);

      int cbc;
      Real t = sicfast(bb1,bb2,cb1,cb2,cbc);
      if(cbc >= CONTACT_TH) {
				Hit h(iss,irt,cbc);
				h.s1.from_four_points(bb1[1],bb1[1],bb1[2],bb1[3]);
				h.s2.from_four_points(bb2[1],bb2[1],bb2[2],bb2[3]);
				h.s1.v += t*Vec(0,0,1);
#ifdef USE_OPENMP
#pragma omp critical
#endif
        hits.push_back(h);
      }
    }
    //break;
  }

  std::sort(hits.begin(),hits.end(),cmp);

	if(hits.size()==0) return;
  TR << hits[1].cbc << std::endl;
  Pose p(init),q(init);
	core::kinematics::Stub s(init.xyz(AtomID(1,1)),init.xyz(AtomID(2,1)),init.xyz(AtomID(3,1)));
	xform_pose_rev(p,s); xform_pose(p,hits[1].s1);
	xform_pose_rev(q,s); xform_pose(q,hits[1].s2);
  p.dump_pdb(utility::file_basename(fn)+"_1_A.pdb");
  q.dump_pdb(utility::file_basename(fn)+"_1_B.pdb");
}


int main(int argc, char *argv[]) {
  core::init(argc,argv);
  using namespace basic::options::OptionKeys;

  vector1<xyzVector<double> > ssamp(NSS);
  {
    izstream is;
    basic::database::open(is,"geometry/sphere_"+str(NSS)+".dat");
    for(int i = 1; i <= NSS; ++i) {
      double x,y,z;
      is >> x >> y >> z;
      ssamp[i] = xyzVector<double>(x,y,z);
    }
    is.close();
  }
  // cout <<"checking geom" << endl;
  // double mn = 9e9;
  // int mi,mj;
  // for(int i = 1; i <= ssamp.size(); ++i) {
  //  for(int j = i+1; j <= ssamp.size(); ++j) {
  //    double d = ssamp[i].distance_squared(ssamp[j]);
  //    if( mn > d ) {
  //      mn = d;
  //      mi = i;
  //      mj = j;
  //    }
  //  }
  // }
  // cout << degrees(asin(sqrt(mn))) << " " << mi << " " << mj << endl;
  // utility_exit_with_message("sphere spacing: "+str(sqrt(mn)));

  for(Size ifn = 1; ifn <= option[in::file::s]().size(); ++ifn) {
    string fn = option[in::file::s]()[ifn];
    Pose pnat;
    TR << "searching " << fn << std::endl;
    core::import_pose::pose_from_pdb(pnat,fn);
    trans_pose(pnat,-com(pnat,1,pnat.n_residue()));
    core::scoring::dssp::Dssp dssp(pnat);
    dssp.insert_ss_into_pose(pnat);
    //if( pnat.n_residue() > 150 ) continue;
    Size cyscnt=0, nhelix=0;
    for(Size ir = 2; ir <= pnat.n_residue()-1; ++ir) {
			if(pnat.secstruct(ir) == 'H') nhelix++;
      //if(!pnat.residue(ir).is_protein()) goto cont1;
      if(pnat.residue(ir).is_lower_terminus()) remove_lower_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).is_upper_terminus()) remove_upper_terminus_type_from_pose_residue(pnat,ir);//goto cont1;
      if(pnat.residue(ir).name3()=="CYS") { if(++cyscnt > 3) goto cont1; }
    } goto done1; cont1: TR << "skipping " << fn << std::endl; continue; done1:
		if( nhelix < 20 ) continue;
    Pose pala(pnat);
    dock(pala,fn,ssamp);
  }
}

