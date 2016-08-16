// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilat/will/symdock_enum.cc
/// @brief docks trimer center and pentamer center together

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/will_util.ihh>


static THREAD_LOCAL basic::Tracer TR( "isct_test" );

using core::Size;
using core::Real;
using core::pose::Pose;
using utility::vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using numeric::min;
using numeric::max;
using numeric::random::uniform;
using numeric::random::gaussian;
using core::id::AtomID;

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

int cbcount_vec(vector1<Vec> & cba, vector1<Vec> & cbb) {
  int cbcount = 0;
  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia)
    for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib)
      if( ib->distance_squared(*ia) < 100.0 ) cbcount++;
  return cbcount;
}

void set_cb_pairs(vector1<Vec> & cba, vector1<Vec> & cbb) {
  vector1<Vec> a,b;
  for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
    for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
      if( ib->distance_squared(*ia) < 100.0 ) {
        a.push_back(*ia);
        b.push_back(*ib);
      }
    }
  }
  cba = a;
  cbb = b;
}

Vec randvec() {
  Vec v(uniform(),uniform(),uniform());
  while( v.length() > 1.0 ) v = Vec(uniform(),uniform(),uniform());
  return v.normalized();
}

int pose_cbcount(Pose const & a, Pose const & b) {
  int count = 0;
  for(Size i = 1; i <= a.n_residue(); ++i) {
    for(Size j = 1; j <= b.n_residue(); ++j) {
      if(a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0) {
        count++;
      }
    }
  }
  return count;
}

#define BIN 2.0

struct IsctFast {
  IsctFast(core::pose::Pose const & pose, Real clash_dis){

		for(Size i = 1; i <= pose.n_residue(); ++i) {
			Size const natom = (pose.residue(i).name3()=="GLY") ? 4 : 5;
			for(Size j = 1; j <= natom; ++j) pnt.push_back(Vec(pose.residue(i).xyz(j)));
		}
		// get bounds for plane hashes
		Real xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9,zmx=-9e9,zmn=9e9;
		for(vector1<Vec>::const_iterator ib = pnt.begin(); ib != pnt.end(); ++ib) {
			xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
			ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
			zmx = max(zmx,ib->z()); zmn = min(zmn,ib->z());
		}

		xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,pppp
		ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms
		zlb = (int)floor(zmn/BIN)-2; int zub = (int)ceil(zmx/BIN)+2; // and one extra for outside atoms

		xsz = xub-xlb+1;
		ysz = yub-ylb+1;
		zsz = zub-zlb+1;
		gxl.dimension((Size)ysz,(Size)zsz); qxl.dimension((Size)ysz,(Size)zsz);
		gxu.dimension((Size)ysz,(Size)zsz);	qxu.dimension((Size)ysz,(Size)zsz);
		gyl.dimension((Size)zsz,(Size)xsz);	qyl.dimension((Size)zsz,(Size)xsz);
		gyu.dimension((Size)zsz,(Size)xsz);	qyu.dimension((Size)zsz,(Size)xsz);
		gzl.dimension((Size)xsz,(Size)ysz);	qzl.dimension((Size)xsz,(Size)ysz);
		gzu.dimension((Size)xsz,(Size)ysz);	qzu.dimension((Size)xsz,(Size)ysz);
		for(vector1<Vec>::const_iterator ia = pnt.begin(); ia != pnt.end(); ++ia) {
			int const ix = (int)ceil(ia->x()/BIN)-(int)xlb;
			int const iy = (int)ceil(ia->y()/BIN)-(int)ylb;
			int const iz = (int)ceil(ia->z()/BIN)-(int)zlb;
			if( 1 < iy && iy < ysz && 1 < iz && iz < zsz ) {
				if( gxl(iy,iz).x() > ia->x() ) gxl(iy,iz) = *ia;
				if( gxu(iy,iz).x() < ia->x() ) gxu(iy,iz) = *ia;
			}
			if( 1 < iz && iz < zsz && 1 < ix && ix < xsz ) {
				if( gyl(iz,ix).y() > ia->y() ) gyl(iz,ix) = *ia;
				if( gyu(iz,ix).y() < ia->y() ) gyu(iz,ix) = *ia;
			}
			if( 1 < ix && ix < xsz && 1 < iy && iy < ysz ) {
				if( gzl(ix,iy).z() > ia->z() ) gzl(ix,iy) = *ia;
				if( gzu(ix,iy).z() < ia->z() ) gzu(ix,iy) = *ia;
			}
		}

  }

	bool test(core::kinematics::Stub const & s) {
		for(vector1<Vec>::const_iterator ia = pnt.begin(); ia != pnt.end(); ++ia) {
			Vec const V = s.local2global(*ia);
			int const ix = (int)ceil(V.x()/BIN)-(int)xlb;
			int const iy = (int)ceil(V.y()/BIN)-(int)ylb;
			int const iz = (int)ceil(V.z()/BIN)-(int)zlb;
			if( 1 < iy && iy < ysz && 1 < iz && iz < zsz ) {
				if( gxl(iy,iz).x() > V.x() ) gxl(iy,iz) = V;
				if( gxu(iy,iz).x() < V.x() ) gxu(iy,iz) = V;
			}
			if( 1 < iz && iz < zsz && 1 < ix && ix < xsz ) {
				if( gyl(iz,ix).y() > V.y() ) gyl(iz,ix) = V;
				if( gyu(iz,ix).y() < V.y() ) gyu(iz,ix) = V;
			}
			if( 1 < ix && ix < xsz && 1 < iy && iy < ysz ) {
				if( gzl(ix,iy).z() > V.z() ) gzl(ix,iy) = V;
				if( gzu(ix,iy).z() < V.z() ) gzu(ix,iy) = V;
			}
		}
	}

  vector1<Vec> pnt;
  ObjexxFCL::FArray2D<Vec> gxl,gxu,gyl,gyu,gzl,gzu;
	mutable ObjexxFCL::FArray2D<Vec> qxl,qxu,qyl,qyu,qzl,qzu;
  Real xlb,ylb,zlb,xsz,ysz,zsz;
};

double
sicfast(
        vector1<Vec>  pa,
        vector1<Vec>  pb,
        Vec ori
        ) {

  // get points, rotated ro ori is 0,0,1, might already be done
  Mat rot = Mat::identity();
  if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  else if( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  if( rot != Mat::identity() ) {
    for(vector1<Vec>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
    for(vector1<Vec>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
  }

  // get bounds for plane hashes
  double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
    ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
    ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
  }
  xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
  ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);


  int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,pppp
  int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

  // TR << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
  // TR << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

  // insert points into hashes
  int const xsize = xub-xlb+1;
  int const ysize = yub-ylb+1;
  ObjexxFCL::FArray2D<Vec> ha(xsize,ysize,Vec(0,0,-9e9)),hb(xsize,ysize,Vec(0,0,9e9));
  for(vector1<Vec>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    // int const ix = min(xsize,max(1,(int)ceil(ia->x()/BIN)-xlb));
    // int const iy = min(ysize,max(1,(int)ceil(ia->y()/BIN)-ylb));
    int const ix = (int)ceil(ia->x()/BIN)-xlb;
    int const iy = (int)ceil(ia->y()/BIN)-ylb;
    if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
    if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
  }
  for(vector1<Vec>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
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

          if( d2 < 16.0 ) {
            double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(16.0-d2);
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

  return mindis;
}


bool grd_clash(Pose const & pose, core::kinematics::Stub const & stub, IsctFast & ifst) {
	ifst.test( stub );
}

bool ifc_clash(Pose const & pose, core::kinematics::Stub const & stub, protocols::scoring::ImplicitFastClashCheck const & ifc) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    for(Size ia = 1; ia <= 5; ++ia) {
      Vec A = stub.local2global(pose.xyz(AtomID(ia,ir)));
      if(!ifc.clash_check( A ) ) return true;
    }
  }
  return false;
}

bool brt_clash(Pose const & pose, core::kinematics::Stub const & stub, Real clash_dis2 = 16.0) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    Vec const CA1 = stub.local2global(pose.xyz(AtomID(2,ir)));
    for(Size jr = 1; jr <= pose.n_residue(); ++jr) {
      Vec const CA2 = pose.xyz(AtomID(2,jr));
      if(CA1.distance_squared(CA2) > 81.0) continue;
      for(Size ia = 1; ia <= 5; ++ia) {
        Vec A = stub.local2global(pose.xyz(AtomID(ia,ir)));
        for(Size ja = 1; ja <= 5; ++ja) {
          Vec B = pose.xyz(AtomID(ja,jr));
          if(A.distance_squared(B) < clash_dis2) {
            //            TR << ir << " " << ia << " " << jr << " " << ja << std::endl;
            return true;
          }
        }
      }
    }
  }
  return false;
}


int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);
  using namespace basic::options::OptionKeys;
  std::string fn = basic::options::option[in::file::s]()[1];
  Pose pose;
  core::import_pose::pose_from_file(pose,fn, core::import_pose::PDB_file);
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    core::pose::replace_pose_residue_copying_existing_coordinates(pose,ir,pose.residue(ir).residue_type_set().name_map("ALA"));
  }


  protocols::scoring::ImplicitFastClashCheck ifc(pose,4.0);
	IsctFast ifast(pose,4.0);
  core::kinematics::Stub stub;
  Real cnt = 0.0;
  Real N = 10000;
  time_t bt=0,it=0,gt=0;
  for(Size iter = 1; iter <= N; ++iter) {

    if(iter%1000==0) TR << iter << std::endl;

    stub.M = rotation_matrix(randvec(),uniform()*2.0*numeric::constants::d::pi);
    stub.v = randvec()*gaussian()*40.0;

    time_t tmp = clock();
    bool c1 = brt_clash(pose,stub,4.0*4.0);
    time_t tmp2 = clock();
    bt += tmp2-tmp;
    bool c2 = ifc_clash(pose,stub,ifc);
		time_t tmp3 = clock();
    it += tmp3-tmp2;
		bool c3 = grd_clash(pose,stub,ifast);
		gt += clock()-tmp3;

    if(c1) cnt += 1.0;

    if( c1 != c2 ) {
      TR << "brt " << c1 << " ifc " << c2 << std::endl;
      pose.dump_pdb("test1.pdb");
      xform_pose(pose,stub);
      pose.dump_pdb("test2.pdb");
      utility_exit_with_message("ARST");
    }

  }
  TR << cnt / N << " " << Real(bt) / Real(it) << " " << Real(it) / Real(gt) << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
