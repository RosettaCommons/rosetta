// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/symdock_enum.cc
/// @brief docks trimer center and pentamer center together

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

static THREAD_LOCAL basic::Tracer TR( "symdock_enum" );

using core::Size;
using core::Real;
using core::pose::Pose;
using utility::vector1;
using namespace numeric;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

typedef xyzVector<core::Real> Vec;
typedef xyzMatrix<core::Real> Mat;
typedef xyzVector<double> Vecf;
typedef xyzMatrix<double> Matf;

void trans_pose( Pose & pose, Vec const & trans ) {
	for(core::Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + trans );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot ) {
	for(core::Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot, Vec const & cen ) {
	trans_pose(pose,-cen);
	rot_pose(pose,rot);
	trans_pose(pose,cen);
}

void rot_pose( Pose & pose, Vec const & axis, core::Real const & ang ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

void rot_pose( Pose & pose, Vec const & axis, core::Real const & ang, Vec const & cen ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
}


void alignaxis(Pose & pose, Vec newaxis, Vec oldaxis, Vec cen = Vec(0,0,0) ) {
	newaxis.normalize();
	oldaxis.normalize();
	Vec axis = newaxis.cross(oldaxis).normalized();
	core::Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
	rot_pose(pose,axis,ang,cen);
}

int cbcount_vec(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
	int cbcount = 0;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia)
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib)
			if( ib->distance_squared(*ia) < 100.0 ) cbcount++;
	return cbcount;
}

void set_cb_pairs(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
	vector1<Vecf> a,b;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			if( ib->distance_squared(*ia) < 100.0 ) {
				a.push_back(*ia);
				b.push_back(*ib);
			}
		}
	}
	cba = a;
	cbb = b;
}

int pose_cbcount(Pose const & a, Pose const & b) {
	int count = 0;
	for(core::Size i = 1; i <= a.n_residue(); ++i) {
		for(core::Size j = 1; j <= b.n_residue(); ++j) {
			if(a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0) {
				count++;
			}
		}
	}
	return count;
}

core::Real
angle_degrees(Vec u, Vec v, Vec w) {
	return numeric::conversions::degrees(numeric::angle_radians(u,v,w));
}

double
sicfast(
	vector1<Vecf>  pa,
	vector1<Vecf>  pb,
	vector1<Vecf> & cba,
	vector1<Vecf> & cbb,
	Vecf ori,
	int & cbcount,
	bool debug = false
) {
	double BIN = 2.0;

	// get points, rotated ro ori is 0,0,1, might already be done
	Matf rot = Matf::identity();
	if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	if( rot != Matf::identity() ) {
		for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
		for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	}

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


	int xlb = floor(xmn/BIN)-2; int xub = ceil(xmx/BIN)+2; // one extra on each side for correctness,
	int ylb = floor(ymn/BIN)-2; int yub = ceil(ymx/BIN)+2; // and one extra for outside atoms

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

	// {
	// 	utility::io::ozstream out("cba.pdb");
	// 	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
	// 		Vec viz = (*ia) + (mindis*ori);
	// 		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// 	}
	// 	out.close();
	// }
	// {
	// 	utility::io::ozstream out("cbb.pdb");
	// 	for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
	// 		Vec viz = (*ib);
	// 		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	// 	}
	// 	out.close();
	// }

	cbcount = 0;
	// utility::io::ozstream out("cb8.pdb");
	// TR << "CB0 " << cbcount << std::endl;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			if( ib->distance_squared( (*ia) + (mindis*ori) ) < 100.0 ) {
				cbcount++;
				// Vec viz = (*ia) + (mindis*ori);
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// viz = *ib;
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
	}
	// out.close();
	// TR << "CB1 " << cbcount << std::endl;

	// // rotate points back -- needed iff pa/pb come by reference
	rot = rot.transposed();
	// if( rot != Matf::identity() ) {
	// 	for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
	// 	for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	// }

	// uncomment this to get hashes in local space
	// rot = Matf::identity();
	// ori = Vec(0,0,1);

	if(debug){
	{
		utility::io::ozstream out("hasha.pdb");
		for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
			for(int j = 2; j <= ysize-1; ++j) {
				Vecf viz = rot*ha(i,j) + mindis*ori;
				if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
				out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"B"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
		Vecf viz = rot*ha(imna,jmna) + mindis*ori;
	    out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"B"<<I(4,100+jmna)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		out.close();
	}
	{
		utility::io::ozstream out("hashb.pdb");
		for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
			for(int j = 2; j <= ysize-1; ++j) {
				Vecf viz = rot*hb(i,j);
				if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
				out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"C"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
		Vecf viz = rot*hb(imnb,jmnb);
		out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"C"<<I(4,100+jmnb)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
		out.close();
	}
	}

	return mindis;
}

double sicfast(
	Pose const & a,
	Pose const & b,
	Vecf ori_in,
	int & cbcount
) {
	// get points, rotated ro ori is 0,0,1
	vector1<Vecf> pa,pb;
	vector1<Vecf> cba,cbb;
	Vecf ori = ori_in.normalized();
	Matf rot = Matf::identity();
	if     ( ori.dot(Vec(0,0,1)) < -0.999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if( ori.dot(Vec(0,0,1)) <  0.999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	for(int i = 1; i <= (int)a.n_residue(); ++i) {
		cba.push_back(rot*Vecf(a.residue(i).xyz(2)));
		int const natom = (a.residue(i).name3()=="GLY") ? 4 : 5;
		for(int j = 1; j <= natom; ++j) pa.push_back(rot*Vecf(a.residue(i).xyz(j)));
	}
	for(int i = 1; i <= (int)b.n_residue(); ++i) {
		cbb.push_back(rot*Vecf(b.residue(i).xyz(2)));
		int const natom = (b.residue(i).name3()=="GLY") ? 4 : 5;
		for(int j = 1; j <= natom; ++j) pb.push_back(rot*Vecf(b.residue(i).xyz(j)));
	}
	return sicfast( pa, pb, cba, cbb, Vec(0,0,1), cbcount );
}

// find 2 HIS that can chelate a tetrahedral metal
void run(  ) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	using numeric::conversions::radians;

	core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);


	Pose t_in,p_in;
	core::import_pose::pose_from_file(t_in,*crs,option[in::file::s]()[1], core::import_pose::PDB_file);
	core::import_pose::pose_from_file(p_in,*crs,option[in::file::s]()[2], core::import_pose::PDB_file);

	// set up geometry
	Vecf taxs = Vec( 0.000000, 0.000000,1.000000).normalized();
	Vecf tax2 = Vec(-0.333333,-0.577350,0.745356).normalized(); // 33.4458470159, 10.42594
	Vecf paxs = Vec(-0.607226, 0.000000,0.794529).normalized();
	Vecf pax2 = Vec(-0.491123,-0.850651,0.187593).normalized(); // 63.4311873349, 5.706642
	core::Real alpha = angle_degrees(taxs,Vec(0,0,0),paxs);
	//core::Real alpha = numeric::conversions::degrees(angle(taxs,Vec(0,0,0),paxs));
	rot_pose(p_in,Vec(0,1,0),-alpha,Vec(0,0,0));

	// Vecf tcom(0,0,0),pcom(0,0,0);
	// for(core::Size i = 1; i <= t_in.n_residue()/3; ++i) for(core::Size j = 1; j <= 5; ++j) tcom += t_in.residue(i).xyz(j);
	// tcom /= double(5*t_in.n_residue()/3);
	// for(core::Size i = 1; i <= p_in.n_residue()/5; ++i) for(core::Size j = 1; j <= 5; ++j) pcom += p_in.residue(i).xyz(j);
	// pcom /= double(5*t_in.n_residue()/5);
	// rot_pose(t_in,taxs,dihedral_degrees(paxs,Vec(0,0,0),taxs,tcom));
	// rot_pose(p_in,paxs,dihedral_degrees(taxs,Vec(0,0,0),paxs,pcom));
	// t_in.dump_pdb("t0.pdb");
	// p_in.dump_pdb("p0.pdb");
	// utility_exit_with_message("debug");

	Pose const tinit(t_in);
	Pose const pinit(p_in);


	int ANGLE_INCR = 6;

	// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
	TR << "precomputing pent-pent and tri-tri interactions" << std::endl;
	vector1<double>            pntmnpos(72/ANGLE_INCR,0),   trimnpos(120/ANGLE_INCR,0);
	vector1<double>            pntmnneg(72/ANGLE_INCR,0),   trimnneg(120/ANGLE_INCR,0);
	ObjexxFCL::FArray2D<int>   pntcbpos(72/ANGLE_INCR,97,0),tricbpos(120/ANGLE_INCR,145,0);
	ObjexxFCL::FArray2D<int>   pntcbneg(72/ANGLE_INCR,97,0),tricbneg(120/ANGLE_INCR,145,0);
	{
		Pose p = pinit;
		for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
			rot_pose(p,paxs,ANGLE_INCR);
			vector1<Vecf> ppnt,cbp; // precompute these
			for(int i = 1; i <= (int)p.n_residue(); ++i) {
				cbp.push_back(Vecf(p.residue(i).xyz(2)));
				int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
				for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
			}
			vector1<Vecf> ppn2(ppnt),cb2(cbp);
			Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
			for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
			int cbcount = 0;
			double const d = sicfast(ppnt,ppn2,cbp,cb2,(pax2-paxs).normalized(),cbcount,false);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pax2-paxs).normalized();
			pntmnpos[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
			pntcbpos(ipnt/ANGLE_INCR+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			set_cb_pairs(cbp,cb2);
			for(int i = 2; i <= 97; ++i) {
				for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1*paxs;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*pax2;
				int cbc = 0; for(core::Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
				pntcbpos(ipnt/ANGLE_INCR+1,i) = cbc;
				if(cbc==0) break;
			}
		}
		for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
			// Pose p = pinit;
			rot_pose(p,paxs,ANGLE_INCR);
			vector1<Vecf> ppnt,cbp; // precompute these
			for(int i = 1; i <= (int)p.n_residue(); ++i) {
				cbp.push_back(Vecf(p.residue(i).xyz(2)));
				int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
				for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
			}
			vector1<Vecf> ppn2(ppnt),cb2(cbp);
			Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
			for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
			int cbcount = 0;
			double const d = sicfast(ppnt,ppn2,cbp,cb2,(paxs-pax2).normalized(),cbcount,false);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(paxs-pax2).normalized();
			pntmnneg[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
			pntcbneg(ipnt/ANGLE_INCR+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			set_cb_pairs(cbp,cb2);
			for(int i = 2; i <= 97; ++i) {
				for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1*paxs;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*pax2;
				int cbc = 0; for(core::Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
				pntcbneg(ipnt/ANGLE_INCR+1,i) = cbc;
				if(cbc==0) break;
			}
		}


		for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
			Pose t = tinit;
			rot_pose(t,taxs,(core::Real)itri);
			vector1<Vecf> ptri,cbt; // precompute these
			for(int i = 1; i <= (int)t.n_residue(); ++i) {
				cbt.push_back(Vecf(t.residue(i).xyz(2)));
				int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
				for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
			}
			vector1<Vecf> ptr2(ptri),cb2(cbt);
			Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
			for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
			int cbcount = 0;
			double const d = sicfast(ptri,ptr2,cbt,cb2,(tax2-taxs).normalized(),cbcount,false);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(tax2-taxs).normalized();
			trimnpos[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
			tricbpos(itri/ANGLE_INCR+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			set_cb_pairs(cbt,cb2);
			for(int i = 2; i <= 145; ++i) {
				for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1*taxs;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*tax2;
				int cbc = 0; for(core::Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
				tricbpos(itri/ANGLE_INCR+1,i) = cbc;
				if(cbc==0) break;
			}
			// if(itri != 60) continue;
			// Pose p1(t),p2(t);
			// rot_pose(p2,r);
			// trans_pose(p1,taxs*trimnpos[itri/ANGLE_INCR+1]);
			// trans_pose(p2,tax2*trimnpos[itri/ANGLE_INCR+1]);
			// p1.dump_pdb("tritri1.pdb");
			// p2.dump_pdb("tritri2.pdb");
			// TR << "PNT NEG D " << -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
			// for(int i = 1; i <= 145; ++i) {
			// 	TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbpos(itri/ANGLE_INCR+1,i) << std::endl;
			// 	trans_pose(p1,taxs*0.1);
			// 	trans_pose(p2,tax2*0.1);
			// }
			// TR << "D " << trimnpos[itri/ANGLE_INCR+1] << std::endl;
			// // utility_exit_with_message("test");
		}
		for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
			Pose t = tinit;
			rot_pose(t,taxs,(core::Real)itri);
			vector1<Vecf> ptri,cbt; // precompute these
			for(int i = 1; i <= (int)t.n_residue(); ++i) {
				cbt.push_back(Vecf(t.residue(i).xyz(2)));
				int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
				for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
			}
			vector1<Vecf> ptr2(ptri),cb2(cbt);
			Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
			for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
			int cbcount = 0;
			double const d = sicfast(ptri,ptr2,cbt,cb2,(taxs-tax2).normalized(),cbcount,false);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(taxs-tax2).normalized();
			trimnneg[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
			tricbneg(itri/ANGLE_INCR+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			set_cb_pairs(cbt,cb2);
			for(int i = 2; i <= 145; ++i) {
				for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1*taxs;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*tax2;
				int cbc = 0; for(core::Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
				tricbneg(itri/ANGLE_INCR+1,i) = cbc;
				if(cbc==0) break;
			}
			// if(itri < 70) continue;
			// Pose p1(t),p2(t);
			// rot_pose(p2,r);
			// trans_pose(p1,taxs*trimnneg[itri/ANGLE_INCR+1]);
			// trans_pose(p2,tax2*trimnneg[itri/ANGLE_INCR+1]);
			// p1.dump_pdb("test1.pdb");
			// p2.dump_pdb("test2.pdb");
			// TR << "PNT NEG D " << d << " " << sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
			// for(int i = 1; i <= 145; ++i) {
			// 	TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbneg(itri/ANGLE_INCR+1,i) << std::endl;
			// 	trans_pose(p1,-taxs*0.1);
			// 	trans_pose(p2,-tax2*0.1);
			// }
			// TR << "D " << trimnneg[itri/ANGLE_INCR+1] << std::endl;
			// utility_exit_with_message("test");
		}
	}

	// utility_exit_with_message("testing");
	TR << "done precomputing pentamer / trimer nb counts" << std::endl;

	ObjexxFCL::FArray3D<int>   cbcount((core::Size)floor(72.0/ANGLE_INCR),(core::Size)floor(120.0/ANGLE_INCR),(core::Size)floor(360.0/ANGLE_INCR),0);
	ObjexxFCL::FArray3D<float> surfdis((core::Size)floor(72.0/ANGLE_INCR),(core::Size)floor(120.0/ANGLE_INCR),(core::Size)floor(360.0/ANGLE_INCR),0.0);
	{
		double mxd = 0;
		int cbmax = 0, mxiori = 0, mxipnt = 0, mxitri = 0;
		for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
			TR << "pent rot: " << ipnt << " cbmax: " << cbmax << std::endl;
			Pose p = pinit;
			rot_pose(p,paxs,(core::Real)ipnt);
			vector1<Vecf> pb,cbb; // precompute these
			for(int i = 1; i <= (int)p.n_residue(); ++i) {
				cbb.push_back(Vecf(p.residue(i).xyz(2)));
				int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
				for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
			}
			for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
				Pose t = tinit;
				rot_pose(t,taxs,(core::Real)itri);
				vector1<Vecf> pa,cba; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cba.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
				}
				int iori = -1, ori_stage = 1;
				bool newstage = true;
				// for(iori = 0; iori < 360; iori+=ANGLE_INCR)
				while(ori_stage < 5) {
					if(newstage) {
						if( ori_stage == 1 || ori_stage == 2 ) iori = ( 90.0+double(ANGLE_INCR)/2.0+angle_degrees(taxs,Vecf(0,0,0),paxs));
						if( ori_stage == 3 || ori_stage == 4 ) iori = (270.0+double(ANGLE_INCR)/2.0+angle_degrees(taxs,Vecf(0,0,0),paxs));
						iori = (iori / ANGLE_INCR) * ANGLE_INCR; // round to closest multiple of angle incr
						if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
						newstage = false;
					} else {
						if( ori_stage == 1 || ori_stage == 3 ) iori += ANGLE_INCR;
						if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
					}
					// TR << "IORI " << iori << std::endl;
					Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(core::Real)iori) * Vec(0,0,1)).normalized();
					int tmpcbc;
					double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
					double theta = iori;
					double gamma = theta-alpha;
					double x = d * sin(numeric::conversions::radians(gamma));
					double y = d * cos(numeric::conversions::radians(gamma));
					double w = x / sin(numeric::conversions::radians(alpha));
					double z = x / tan(numeric::conversions::radians(alpha));
					double dpnt = y+z;
					double dtri = w;
					double pntmn,trimn;
					if( w > 0 ) {
						pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
						trimn = trimnpos[itri/ANGLE_INCR+1];
						int dp = (dpnt-pntmn)*10+1;
						int dt = (dtri-trimn)*10+1;
						if( dp < 1 )  { ori_stage++; newstage=true; continue; };
						if( dt < 1 )  { ori_stage++; newstage=true; continue; };
						// if(ipnt==18 && itri==72 && iori==276)
						// 	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
						// 	  << "    " << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
						if( dp <= 97  ) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
						if( dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
						// TR << "CHK " << dpnt << " " << pntmn << "    " << dtri << " " << trimn << std::endl;
					} else {
						pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
						trimn = trimnneg[itri/ANGLE_INCR+1];
						int dp = (-dpnt+pntmn)*10+1;
						int dt = (-dtri+trimn)*10+1;
						if( dp < 1 )  { ori_stage++; newstage=true; continue; };
						if( dt < 1 )  { ori_stage++; newstage=true; continue; };
						// if(ipnt==18 && itri==72 && iori==276)
						// 	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
						// 	  << "    " << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
						if( dp <= 97  ) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
						if( dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
					}

					surfdis(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = d;
					cbcount(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = tmpcbc;
					// d = sicfast(t,p,sicaxis,cbcount);
					// TR << "trial " << ipnt << " " << itri << " " << iori << " " << d << " " << cbcount << std::endl;
					if(tmpcbc > cbmax) {
						cbmax = tmpcbc;
						mxiori = iori;
						mxipnt = ipnt;
						mxitri = itri;
						mxd = d;
					}
				}
			}
		}
		TR << "MAX " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << " " << mxd << std::endl;
		// TR << "MAX " << mxipnt << " " << mxitri << " " << mxiori << " " << cbcount(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << " " << surfdis(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << std::endl;
	}

	ObjexxFCL::FArray3D<int> cbcount_local((core::Size)floor(72.0/ANGLE_INCR),(core::Size)floor(120.0/ANGLE_INCR),(core::Size)floor(360.0/ANGLE_INCR),0);
	int delta = ceil(5.0/(core::Real)ANGLE_INCR);
	TR << "scanning results for local maxima +- " << delta*ANGLE_INCR << "°" << std::endl;
	for(int ipnt = 1; ipnt <=  72/ANGLE_INCR; ++ipnt) {
		for(int itri = 1; itri <= 120/ANGLE_INCR; ++itri) {
			for(int iori = 1; iori <= 360/ANGLE_INCR; ++iori) {
				int lmaxcb = 0;
				for(int dipnt = -delta; dipnt <= delta; ++dipnt) {
					for(int ditri = -delta; ditri <= delta; ++ditri) {
						for(int diori = -delta; diori <= delta; ++diori) {
							int i = (ipnt+dipnt-1) % ( 72/ANGLE_INCR) + 1;
							int j = (itri+ditri-1) % (120/ANGLE_INCR) + 1;
							int k = (iori+diori-1) % (360/ANGLE_INCR) + 1;
							if( cbcount(i,j,k) > lmaxcb ) lmaxcb = cbcount(i,j,k);
						}
					}
				}
				if( cbcount(ipnt,itri,iori) >= lmaxcb ) cbcount_local(ipnt,itri,iori) = cbcount(ipnt,itri,iori);
			}
		}
	}

	TR << "cbcount_local size " << cbcount_local.size() << std::endl;
	vector1<int> cbtmp;
	for(core::Size i = 0; i < cbcount_local.size(); ++i) {
		cbtmp.push_back(cbcount_local[i]);
	}
	std::sort(cbtmp.begin(),cbtmp.end());
	int top10 = cbtmp[cbtmp.size()-9];
	TR << "outputting top10 with cb count >= " << top10 << std::endl;

	for(int ipnt = 1; ipnt <=  72/ANGLE_INCR; ++ipnt) {
		for(int itri = 1; itri <= 120/ANGLE_INCR; ++itri) {
			for(int iori = 1; iori <= 360/ANGLE_INCR; ++iori) {
				if(cbcount_local(ipnt,itri,iori) < top10) continue;

				std::string fname;
				fname +=       utility::file_basename(option[in::file::s]()[1]);
				fname += "_" + utility::file_basename(option[in::file::s]()[2]);
				fname += "_" + ObjexxFCL::lead_zero_string_of(cbcount_local(ipnt,itri,iori),4);
				fname += "_" + ObjexxFCL::lead_zero_string_of((ipnt-1)*ANGLE_INCR,3);
				fname += "_" + ObjexxFCL::lead_zero_string_of((itri-1)*ANGLE_INCR,3);
				fname += "_" + ObjexxFCL::lead_zero_string_of((iori-1)*ANGLE_INCR,3);
				fname += ".pdb.gz";

				TR << "dumping file " << fname << std::endl;

				Pose p = pinit;
				Pose t = tinit;
				rot_pose(p,paxs,(core::Real)((ipnt-1)*ANGLE_INCR));
				rot_pose(t,taxs,(core::Real)((itri-1)*ANGLE_INCR));

				Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(core::Real)((iori-1)*ANGLE_INCR)) * Vec(0,0,1)).normalized();

				core::Real theta = ((iori-1)*ANGLE_INCR);
				core::Real gamma = theta-alpha;
				core::Real x = surfdis(ipnt,itri,iori) * sin(numeric::conversions::radians(gamma));
				core::Real y = surfdis(ipnt,itri,iori) * cos(numeric::conversions::radians(gamma));
				core::Real w = x / sin(numeric::conversions::radians(alpha));
				core::Real z = x / tan(numeric::conversions::radians(alpha));
				trans_pose(p,(y+z)*paxs);
				trans_pose(t,  w  *taxs);

				Pose symm;
				symm.append_residue_by_jump(t.residue(1),1);
				for(core::Size i = 2; i <= t.n_residue()/3; ++i) {
					if(symm.residue(i-1).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(symm,i-1);
					if(   t.residue(i  ).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(t,i);
					symm.append_residue_by_bond(t.residue(i));
				}
				symm.append_residue_by_jump(p.residue(1),1);
				for(core::Size i = 2; i <= p.n_residue()/5; ++i) symm.append_residue_by_bond(p.residue(i));

				core::pose::symmetry::make_symmetric_pose(symm);
				core::io::pdb::old_dump_pdb(symm,fname);

			}
		}
	}

}


int main (int argc, char *argv[]) {
	try{
	devel::init(argc,argv);
	run();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


