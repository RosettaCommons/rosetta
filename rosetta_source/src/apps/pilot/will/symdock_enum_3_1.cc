#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <utility/vector1.hh>
static basic::Tracer TR("symdock_enum");

// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif



OPT_1GRP_KEY( Real , tcdock, clash_dis		)
OPT_1GRP_KEY( Real , tcdock, contact_dis	)
OPT_1GRP_KEY( Real , tcdock, intra	)
OPT_1GRP_KEY( Integer , tcdock, topx	)
OPT_1GRP_KEY( Boolean , tcdock, reverse	)
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( tcdock::clash_dis	 ,"max acceptable clash dis"	,	3.5 );
	NEW_OPT( tcdock::contact_dis ,"max acceptable contact dis", 10.0 );
	NEW_OPT( tcdock::intra       ,"include intra 3-3 5-5 contacts", 1.0 );
	NEW_OPT( tcdock::topx        ,"output top X hits", 10 );
	NEW_OPT( tcdock::reverse     ,"rev.", false );	
}


using core::Size;
using core::pose::Pose;
using utility::vector1;
using ObjexxFCL::fmt::I;
using ObjexxFCL::fmt::F;
using numeric::min;
using numeric::max;

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzVector<double> Vecf;
typedef numeric::xyzMatrix<double> Matf;

template<typename T> inline T sqr(T x) { return x*x; }

inline double sigmoid( double const & sqdist, double const & start, double const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		double dist = sqrt( sqdist );
		return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

void trans_pose( Pose & pose, Vecf const & trans ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot, Vecf const & cen ) {
	trans_pose(pose,-cen);
	rot_pose(pose,rot);
	trans_pose(pose,cen);
}

void rot_pose( Pose & pose, Vecf const & axis, double const & ang ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

void rot_pose( Pose & pose, Vecf const & axis, double const & ang, Vecf const & cen ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
}


void alignaxis(Pose & pose, Vecf newaxis, Vecf oldaxis, Vecf cen = Vecf(0,0,0) ) {
	newaxis.normalize();
	oldaxis.normalize();
	Vecf axis = newaxis.cross(oldaxis).normalized();
	double ang = (double)-acos(numeric::max((double)-1.0,numeric::min((double)1.0,newaxis.dot(oldaxis))))*(double)180.0/numeric::constants::f::pi;
	rot_pose(pose,axis,ang,cen);
}

int cbcount_vec(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
	int cbcount = 0;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia)
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib)
			if( ib->distance_squared(*ia) < 100.0 ) cbcount++;
	return cbcount;
}

void prune_cb_pairs_dis10(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
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
	for(Size i = 1; i <= a.n_residue(); ++i) {
		for(Size j = 1; j <= b.n_residue(); ++j) {
			if(a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0) {
				count++;
			}
		}
	}
	return count;
}



double
sicsafe(
				vector1<Vecf> const & pa,
				vector1<Vecf> const & pb,
				vector1<Vecf> const & cba,
				vector1<Vecf> const & cbb,
				Vecf ori,
				double & cbcount,
				bool debug = false
){
	double const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
	double const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
	double const CONTACT_D2 = sqr(CONTACT_D);
	double const CLASH_D2	 = sqr(CLASH_D);
	ori.normalize();
	Vecf axis = Vecf(0,0,1).cross(ori);
	Matf R = rotation_matrix( axis, (double)-acos(numeric::max((double)-1.0,numeric::min((double)1.0,Vecf(0,0,1).dot(ori)))) );
	double mindis = 9e9;
	int mni,mnj;
	// for(vector1<Vecf>::const_iterator i = pa.begin(); i != pa.end(); ++i) {
	// for(vector1<Vecf>::const_iterator j = pb.begin(); j != pb.end(); ++j) {
	for(int ii = 1; ii <= pa.size(); ++ii) {
		Vecf const i = R * pa[ii];
		for(int jj = 1; jj <= pb.size(); ++jj) {		
			Vecf const j = R * pb[jj];
			double const dxy2 = (i.x()-j.x())*(i.x()-j.x()) + (i.y()-j.y())*(i.y()-j.y());
			if( dxy2 >= CLASH_D2 ) continue;
			double const dz = j.z() - i.z() - sqrt(CLASH_D2-dxy2);
			if( dz < mindis) {
				mindis = dz;
				mni = ii;
				mnj = jj;
			}
		}		
	}
	
	cbcount = 0.0;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		Vecf const & va = *ia;
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			Vecf const & vb = *ib;
			double d2 = vb.distance_squared( (va + mindis*ori ) );
			if( d2 < CONTACT_D2 ) {
				cbcount += sigmoid(d2, CLASH_D, CONTACT_D );
			}
		}
	}
	return mindis;
	
}

double
sicfast(
				vector1<Vecf>	pa,
				vector1<Vecf>	pb,
				vector1<Vecf> const & cba,
				vector1<Vecf> const & cbb,
				Vecf ori,
				double & cbcount,
				bool debug = false
){
	double const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
	double const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
	double const CONTACT_D2 = sqr(CONTACT_D);
	double const CLASH_D2	 = sqr(CLASH_D);	
	double const BIN = CLASH_D / 2.0;

	// get points, rotated ro ori is 0,0,1, might already be done
	Matf rot = Matf::identity();
	if		 ( ori.dot(Vecf(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vecf(1,0,0).cross(ori), (double)-acos(Vecf(0,0,1).dot(ori)) );
	else if( ori.dot(Vecf(0,0,1)) <	 0.99999 ) rot = rotation_matrix( Vecf(0,0,1).cross(ori), (double)-acos(Vecf(0,0,1).dot(ori)) );
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


	int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
	int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

	// std::cerr << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
	// std::cerr << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

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
					double const xa = ha(i	,j	).x();
					double const ya = ha(i	,j	).y();
					double const xb = hb(i+k,j+l).x();
					double const yb = hb(i+k,j+l).y();
					double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);

					if( d2 < CLASH_D2 ) {
						double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(CLASH_D2-d2);
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
	//	utility::io::ozstream out("cba.pdb");
	//	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
	//		Vecf viz = (*ia) + (mindis*ori);
	//		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"A"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//	}
	//	out.close();
	// }
	// {
	//	utility::io::ozstream out("cbb.pdb");
	//	for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
	//		Vecf viz = (*ib);
	//		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//	}
	//	out.close();
	// }

	cbcount = 0;
	// utility::io::ozstream out("cb8.pdb");
	// std::cerr << "CB0 " << cbcount << std::endl;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			double d2 = ib->distance_squared( (*ia) + (mindis*ori) );
			if( d2 < CONTACT_D2 ) {
				cbcount += sigmoid(d2, CLASH_D, CONTACT_D );
				// Vecf viz = (*ia) + (mindis*ori);
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"A"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// viz = *ib;
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"B"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
	}
	// out.close();
	// std::cerr << "CB1 " << cbcount << std::endl;

	// // rotate points back -- needed iff pa/pb come by reference
	rot = rot.transposed();
	// if( rot != Matf::identity() ) {
	//	for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
	//	for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	// }

	// uncomment this to get hashes in local space
	// rot = Matf::identity();
	// ori = Vecf(0,0,1);

	if(debug){
		{
			utility::io::ozstream out("hasha.pdb");
			for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for(int j = 2; j <= ysize-1; ++j) {
					Vecf viz = rot*ha(i,j) + mindis*ori;
					if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100+j)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vecf viz = rot*ha(imna,jmna) + mindis*ori;
			out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"B"<<I(4,100+jmna)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
		{
			utility::io::ozstream out("hashb.pdb");
			for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for(int j = 2; j <= ysize-1; ++j) {
					Vecf viz = rot*hb(i,j);
					if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"C"<<I(4,100+j)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vecf viz = rot*hb(imnb,jmnb);
			out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"C"<<I(4,100+jmnb)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
	}

	return mindis;
}

double sicfast(
							 Pose const & a,
							 Pose const & b,
							 Vecf ori_in,
							 double & cbcount
){
	// get points, rotated ro ori is 0,0,1
	vector1<Vecf> pa,pb;
	vector1<Vecf> cba,cbb;
	Vecf ori = ori_in.normalized();
	Matf rot = Matf::identity();
	if		 ( ori.dot(Vecf(0,0,1)) < -0.999 ) rot = rotation_matrix( Vecf(1,0,0).cross(ori), (double)-acos(Vecf(0,0,1).dot(ori)) );
	else if( ori.dot(Vecf(0,0,1)) <	0.999 ) rot = rotation_matrix( Vecf(0,0,1).cross(ori), (double)-acos(Vecf(0,0,1).dot(ori)) );
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
	return sicfast( pa, pb, cba, cbb, Vecf(0,0,1), cbcount );
}

// void tri1_to_tri2_pose(core::pose::Pose & pose, Vecf triaxs_, Vecf triax2_) {
// 	Matf r = rotation_matrix_degrees( triaxs_.cross(triax2_), angle_degrees(triaxs_,Vecf(0,0,0),triax2_) );
// 	r = r * rotation_matrix_degrees(triaxs_,60.0);
// 	rot_pose(pose,r);
// }
void tri1_to_tri2(vector1<Vecf> & pts, Vecf triaxs_, Vecf triax2_) {
	Matf r = rotation_matrix_degrees( triaxs_.cross(triax2_), angle_degrees(triaxs_,Vecf(0,0,0),triax2_) );
	r = r * rotation_matrix_degrees(triaxs_,(double)60.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}
// void pnt1_to_pnt2_pose(core::pose::Pose & pose, Vecf pntaxs_, Vecf pntax2_) {
// 	Matf r = rotation_matrix_degrees( pntaxs_.cross(pntax2_), angle_degrees(pntaxs_,Vecf(0,0,0),pntax2_) );
// 	r = r * rotation_matrix_degrees(pntaxs_,36.0);
// 	rot_pose(pose,r);
// }
void pnt1_to_pnt2(vector1<Vecf> & pts, Vecf pntaxs_, Vecf pntax2_) {
	Matf r = rotation_matrix_degrees( pntaxs_.cross(pntax2_), angle_degrees(pntaxs_,Vecf(0,0,0),pntax2_) );
	r = r * rotation_matrix_degrees(pntaxs_,(double)36.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}





struct TCDock {

	vector1<double>						  pntmnpos3_,trimnpos3_;
	vector1<double>						  pntmnneg3_,trimnneg3_;
	ObjexxFCL::FArray2D<double>	pntcbpos3_,tricbpos3_;
	ObjexxFCL::FArray2D<double>	pntcbneg3_,tricbneg3_;
	vector1<double>						  pntmnpos1_,trimnpos1_;
	vector1<double>						  pntmnneg1_,trimnneg1_;
	vector1<double>						  pntdspos1_,tridspos1_;
	vector1<double>						  pntdsneg1_,tridsneg1_;
	ObjexxFCL::FArray2D<double>	pntcbpos1_,tricbpos1_;
	ObjexxFCL::FArray2D<double>	pntcbneg1_,tricbneg1_;

	ObjexxFCL::FArray3D<double> gdis;
	ObjexxFCL::FArray3D<double> gcbc;	

	Vecf triaxs_;
	Vecf triax2_;
	Vecf pntaxs_;
	Vecf pntax2_;
	double alpha_;

	core::pose::Pose tri_in_,pnt_in_;
	vector1<Vecf> tri_pts_,tri_cbs_,pnt_pts_,pnt_cbs_;
	std::string trifile,pntfile;
	
	TCDock( Size itrifile, Size ipntfile ) : 
		pntmnpos3_(24,0.0)   ,trimnpos3_(40,0.0),
		pntmnneg3_(24,0.0)   ,trimnneg3_(40,0.0),
		pntcbpos3_(24,97,0.0),tricbpos3_(40,145,0.0),
		pntcbneg3_(24,97,0.0),tricbneg3_(40,145,0.0),
		pntmnpos1_(72,0.0)	 ,trimnpos1_(120,0.0),
		pntmnneg1_(72,0.0)	 ,trimnneg1_(120,0.0),
		pntdspos1_(72,0.0)	 ,tridspos1_(120,0.0),
		pntdsneg1_(72,0.0)	 ,tridsneg1_(120,0.0),
		pntcbpos1_(72,97,0.0),tricbpos1_(120,145,0.0), // axes moving out
		pntcbneg1_(72,97,0.0),tricbneg1_(120,145,0.0),
		gdis(72,120,360,-1.0),
		gcbc(72,120,360,-1.0)
	{
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);
		core::import_pose::pose_from_pdb(tri_in_,*crs,option[in::file::s]()[itrifile]);
		core::import_pose::pose_from_pdb(pnt_in_,*crs,option[in::file::s]()[ipntfile]);
		trifile = utility::file_basename(option[in::file::s]()[itrifile]);
		pntfile = utility::file_basename(option[in::file::s]()[ipntfile]);
		make_trimer(tri_in_);
		make_pentamer(pnt_in_);		
		// set up geometry
		triaxs_ = Vecf( 0.000000, 0.000000,1.000000).normalized();
		triax2_ = Vecf(-0.333333,-0.577350,0.745356).normalized(); // 33.4458470159, 10.42594
		pntaxs_ = Vecf(-0.607226, 0.000000,0.794529).normalized();
		pntax2_ = Vecf(-0.491123,-0.850651,0.187593).normalized(); // 63.4311873349, 5.706642
		alpha_ = angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_);
		rot_pose(pnt_in_,Vecf(0,1,0),-alpha_,Vecf(0,0,0));
		if(option[tcdock::reverse]()) rot_pose(tri_in_,Vecf(0,1,0),180.0);
		
		// init points
		// for(Size i = 1; i <= tri_in_.n_residue(); ++i) {
		// 	if(tri_in_.residue(i).has("CB")) tri_cbs_.push_back(Vecf(tri_in_.residue(i).xyz("CB")));
		// 	int const natom = (tri_in_.residue(i).name3()=="GLY") ? 4 : 5;
		// 	for(int j = 1; j <= natom; ++j) tri_pts_.push_back(Vecf(tri_in_.residue(i).xyz(j)));
		// }
		// for(Size i = 1; i <= pnt_in_.n_residue(); ++i) {
		// 	if(pnt_in_.residue(i).has("CB")) pnt_cbs_.push_back(Vecf(pnt_in_.residue(i).xyz("CB")));
		// 	int const natom = (pnt_in_.residue(i).name3()=="GLY") ? 4 : 5;
		// 	for(int j = 1; j <= natom; ++j) pnt_pts_.push_back(Vecf(pnt_in_.residue(i).xyz(j)));
		// }		
		for(Size i = 1; i <= tri_in_.n_residue(); ++i) {
			tri_cbs_.push_back(Vecf(tri_in_.residue(i).xyz(2)));
			int const natom = (tri_in_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) tri_pts_.push_back(Vecf(tri_in_.residue(i).xyz(j)));
		}
		for(Size i = 1; i <= pnt_in_.n_residue(); ++i) {
			pnt_cbs_.push_back(Vecf(pnt_in_.residue(i).xyz(2)));
			int const natom = (pnt_in_.residue(i).name3()=="GLY") ? 4 : 5;
			for(int j = 1; j <= natom; ++j) pnt_pts_.push_back(Vecf(pnt_in_.residue(i).xyz(j)));
		}		
		
	}

	void get_tri(double atri, vector1<Vecf> & pts, vector1<Vecf> & cbs ) {
		Matf R = rotation_matrix_degrees( triaxs_, atri );
		pts = tri_pts_;
		cbs = tri_cbs_;		
		for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) (*i) = R*(*i);
		for(vector1<Vecf>::iterator i = cbs.begin(); i != cbs.end(); ++i) (*i) = R*(*i);
	}
	void get_pnt(double apnt, vector1<Vecf> & pts, vector1<Vecf> & cbs ) {
		Matf R = rotation_matrix_degrees( pntaxs_, apnt );
		pts = pnt_pts_;
		cbs = pnt_cbs_;		
		for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) (*i) = R*(*i);
		for(vector1<Vecf>::iterator i = cbs.begin(); i != cbs.end(); ++i) (*i) = R*(*i);
	}
	
	

	double
	calc_iface(
					vector1<Vecf>	pb,
					vector1<Vecf>	pa,
					vector1<Vecf> const & cbb,
					vector1<Vecf> const & cba,
					Vecf sicaxis,
					int & tmpcbc,
					double & dout				
	){
		return true;
	}

	void make_dimer(core::pose::Pose & pose) {
		core::pose::Pose t2(pose);
		rot_pose(t2,Vecf(0,0,1),180.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
	}
	void make_trimer(core::pose::Pose & pose) {
		core::pose::Pose t2(pose),t3(pose);
		rot_pose(t2,Vecf(0,0,1),120.0);
		rot_pose(t3,Vecf(0,0,1),240.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
	}
	void make_pentamer(core::pose::Pose & pose) {
		core::pose::Pose t2(pose),t3(pose),t4(pose),t5(pose);
		rot_pose(t2,Vecf(0,0,1), 72.0);
		rot_pose(t3,Vecf(0,0,1),144.0);
		rot_pose(t4,Vecf(0,0,1),216.0);
		rot_pose(t5,Vecf(0,0,1),288.0);
		for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
		for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
		for(Size i = 1; i <= t4.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t4.residue(i),1); else pose.append_residue_by_bond(t4.residue(i));
		for(Size i = 1; i <= t5.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t5.residue(i),1); else pose.append_residue_by_bond(t5.residue(i));
	}
		
	void dump_pdb(double apnt, double atri, double aori, std::string fname, bool sym=true) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		vector1<Vec> pb,pa,cbb,cba;
		get_tri(atri,pa,cba);
		get_pnt(apnt,pb,cbb);

		Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)aori) * Vecf(0,0,1)).normalized();
		double tmpcbc;
		double d = sicsafe(pb,pa,cbb,cba,sicaxis,tmpcbc);
		double theta = aori;
		double gamma = theta-alpha_;
		double x = d * sin(numeric::conversions::radians(gamma));
		double y = d * cos(numeric::conversions::radians(gamma));
		double w = x / sin(numeric::conversions::radians(alpha_));
		double z = x / tan(numeric::conversions::radians(alpha_));
		double dpnt = y+z;
		double dtri = w;
		double pntmn,trimn;

		Pose t,p,symm;
		#ifdef USE_OPENMP
		#pragma omp critical
		#endif
		{
			t = tri_in_;
			p = pnt_in_;
		}
		rot_pose  (p,pntaxs_,apnt);
		rot_pose  (t,triaxs_,atri);
		trans_pose(p,dpnt*pntaxs_);
		trans_pose(t,dtri*triaxs_);
		{
			symm.append_residue_by_jump(t.residue(1),1);
			for(Size i = 2; i <= t.n_residue()/3; ++i) {
				if(symm.residue(i-1).is_terminus()) symm.append_residue_by_jump(t.residue(i),1);
				else                                symm.append_residue_by_bond(t.residue(i));
			}
			symm.append_residue_by_jump(p.residue(1),1);
			for(Size i = 2; i <= p.n_residue()/5; ++i) {
				if(symm.residue(symm.n_residue()).is_terminus()) symm.append_residue_by_jump(p.residue(i),1);
				else                                             symm.append_residue_by_bond(p.residue(i));
			}
			// core::util::switch_to_residue_type_set(symm,"fa_standard");
			// for(Size i = 1; i <= symm.n_residue(); ++i) core::pose::replace_pose_residue_copying_existing_coordinates(symm,i,rs->name_map("ALA"));
			//std::cerr << "making symm" << std::endl;					
		}
		if(sym)	core::pose::symmetry::make_symmetric_pose(symm);				
		core::io::pdb::dump_pdb(symm,option[out::file::o]()+"/"+fname);
		
	}

	double calc_cbc(int ipnt, int itri, int iori, int ANGLE_INCR, double & dpnt, double & dtri, bool safe=true) {
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		vector1<Vecf> pb,cbb; // precompute these
		vector1<Vecf> pa,cba; // precompute these
		get_tri(itri,pa,cba);
		get_pnt(ipnt,pb,cbb);		
		Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
		double tmpcbc,d;
		if(safe) d = sicsafe(pb,pa,cbb,cba,sicaxis,tmpcbc);
		else     d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);

		double theta = iori;
		double gamma = theta-alpha_;
		double x = d * sin(numeric::conversions::radians(gamma));
		double y = d * cos(numeric::conversions::radians(gamma));
		double w = x / sin(numeric::conversions::radians(alpha_));
		double z = x / tan(numeric::conversions::radians(alpha_));
		dpnt = y+z;
		dtri = w;
		double pntmn,trimn;
		if( w > 0 ) {
			pntmn = pntmnpos1_[ipnt/ANGLE_INCR+1];
			trimn = trimnpos1_[itri/ANGLE_INCR+1];
			if( dpnt < pntmn ) utility_exit_with_message("BAD CONFIG");
			if( dtri < trimn ) utility_exit_with_message("BAD CONFIG");
			int dp = (int)(dpnt-pntmn)*10+1;
			int dt = (int)(dtri-trimn)*10+1;
			// if(ipnt==18 && itri==72 && iori==276)
			//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos1_(ipnt/ANGLE_INCR+1,dp)
			//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos1_(itri/ANGLE_INCR+1,dt) << std::endl;
			if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbpos1_(ipnt/ANGLE_INCR+1,dp);
			if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos1_(itri/ANGLE_INCR+1,dt);
			// std::cerr << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
		} else {
			pntmn = pntmnneg1_[ipnt/ANGLE_INCR+1];
			trimn = trimnneg1_[itri/ANGLE_INCR+1];
			if( dpnt > pntmn ) utility_exit_with_message("BAD CONFIG");
			if( dtri > trimn ) utility_exit_with_message("BAD CONFIG");
			int dp = (int)(-dpnt+pntmn)*10+1;
			int dt = (int)(-dtri+trimn)*10+1;
			// if(ipnt==18 && itri==72 && iori==276)
			//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg1_(ipnt/ANGLE_INCR+1,dp)
			//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg1_(itri/ANGLE_INCR+1,dt) << std::endl;
			if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbneg1_(ipnt/ANGLE_INCR+1,dp);
			if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg1_(itri/ANGLE_INCR+1,dt);
		}
		return tmpcbc;
	}

	void run() {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using numeric::conversions::radians;

		double const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
		double const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
		double const CONTACT_D2 = sqr(CONTACT_D);
		double const CLASH_D2	 = sqr(CLASH_D);

		Pose const triinit(tri_in_);
		Pose const pntinit(pnt_in_);

		int ANGLE_INCR = 3;
		if(ANGLE_INCR != 3) utility_exit_with_message("first ANGLE_INCR must be 3!!!");
		{
			// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
			std::cerr << "precomputing pent-pent and tri-tri interactions every 3°" << std::endl;
			{
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
					vector1<Vecf> ppnt,cbp; // precompute these
					get_pnt(ipnt,ppnt,cbp);
					vector1<Vecf> ppn2(ppnt),cb2(cbp);
					pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
					pnt1_to_pnt2( cb2,pntaxs_,pntax2_);				
					double cbcount = 0.0;
					double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntax2_-pntaxs_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntax2_-pntaxs_).normalized();
					pntmnpos3_[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
					pntcbpos3_(ipnt/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbp,cb2);
					for(int i = 2; i <= 97; ++i) {
						for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1f*pntaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*pntax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						pntcbpos3_(ipnt/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
					vector1<Vecf> ppnt,cbp; // precompute these
					get_pnt(ipnt,ppnt,cbp);
					vector1<Vecf> ppn2(ppnt),cb2(cbp);
					pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
					pnt1_to_pnt2( cb2,pntaxs_,pntax2_);
					double cbcount = 0;
					double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntaxs_-pntax2_).normalized(),cbcount,false);				
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntaxs_-pntax2_).normalized();
					pntmnneg3_[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
					pntcbneg3_(ipnt/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbp,cb2);
					for(int i = 2; i <= 97; ++i) {
						for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1f*pntaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*pntax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						pntcbneg3_(ipnt/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}								
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
					vector1<Vecf> ptri,cbt; // precompute these
					get_tri(itri,ptri,cbt);					
					vector1<Vecf> ptr2(ptri),cb2(cbt);
					tri1_to_tri2(ptr2,triaxs_,triax2_);
					tri1_to_tri2( cb2,triaxs_,triax2_);				
					double cbcount = 0;
					double const d = sicsafe(ptri,ptr2,cbt,cb2,(triax2_-triaxs_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triax2_-triaxs_).normalized();
					trimnpos3_[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
					tricbpos3_(itri/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbt,cb2);
					for(int i = 2; i <= 145; ++i) {
						for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1f*triaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*triax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						tricbpos3_(itri/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
					vector1<Vecf> ptri,cbt; // precompute these
					get_tri(itri,ptri,cbt);
					vector1<Vecf> ptr2(ptri),cb2(cbt);
					tri1_to_tri2(ptr2,triaxs_,triax2_);
					tri1_to_tri2( cb2,triaxs_,triax2_);				
					double cbcount = 0;
					double const d = sicsafe(ptri,ptr2,cbt,cb2,(triaxs_-triax2_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triaxs_-triax2_).normalized();
					trimnneg3_[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
					tricbneg3_(itri/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbt,cb2);
					for(int i = 2; i <= 145; ++i) {
						for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1f*triaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*triax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						tricbneg3_(itri/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
			}
		
		
			std::cerr << "main loop 1 over ipnt, itri, iori every 3 degrees" << std::endl;
			// ObjexxFCL::FArray3D<int>	 cbcount((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
			// ObjexxFCL::FArray3D<double> surfdis3((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0.0);
			{
				double mxd = 0;
				double cbmax = 0; int mxiori = 0, mxipnt = 0, mxitri = 0;
				for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
					if(ipnt%9==0 && ipnt!=0) std::cerr << "	 loop1 " << trifile << " " << (100*ipnt)/72 << "\% done, cbmax: " << cbmax << std::endl;
					vector1<Vecf> pb,cbb; // precompute these
					get_pnt(ipnt,pb,cbb);
					#ifdef USE_OPENMP
					#pragma omp parallel for schedule(dynamic,1)
					#endif
					for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
						vector1<Vecf> pa,cba; // precompute these
						get_tri(itri,pa,cba);

						int iori = -1, ori_stage = 1;
						bool newstage = true;
						// for(iori = 0; iori < 360; iori+=ANGLE_INCR)
						while(ori_stage < 5) {
							if(newstage) {
								if( ori_stage == 1 || ori_stage == 2 ) iori = (int)( 90.0+double(ANGLE_INCR)/2.0+angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_));
								if( ori_stage == 3 || ori_stage == 4 ) iori = (int)(270.0+double(ANGLE_INCR)/2.0+angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_));
								iori = (iori / ANGLE_INCR) * ANGLE_INCR; // round to closest multiple of angle incr
								if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
								newstage = false;
							} else {
								if( ori_stage == 1 || ori_stage == 3 ) iori += ANGLE_INCR;
								if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
							}
							// std::cerr << "IORI " << iori << std::endl;
							Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
							double tmpcbc;
							double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);

							double theta = iori;
							double gamma = theta-alpha_;
							double x = d * sin(numeric::conversions::radians(gamma));
							double y = d * cos(numeric::conversions::radians(gamma));
							double w = x / sin(numeric::conversions::radians(alpha_));
							double z = x / tan(numeric::conversions::radians(alpha_));
							double dpnt = y+z;
							double dtri = w;
							double pntmn,trimn;
							if( w > 0 ) {
								pntmn = pntmnpos3_[ipnt/ANGLE_INCR+1];
								trimn = trimnpos3_[itri/ANGLE_INCR+1];
								if( dtri < trimn ) { ori_stage++; newstage=true; continue; };
								if( dpnt < pntmn ) { ori_stage++; newstage=true; continue; };
								int dp = (int)(dpnt-pntmn)*10+1;
								int dt = (int)(dtri-trimn)*10+1;
								// if(ipnt==18 && itri==72 && iori==276)
								//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos3_(ipnt/ANGLE_INCR+1,dp)
								//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos3_(itri/ANGLE_INCR+1,dt) << std::endl;
								if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbpos3_(ipnt/ANGLE_INCR+1,dp);
								if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos3_(itri/ANGLE_INCR+1,dt);
								// std::cerr << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
							} else {
								pntmn = pntmnneg3_[ipnt/ANGLE_INCR+1];
								trimn = trimnneg3_[itri/ANGLE_INCR+1];
								if( dtri > trimn ) { ori_stage++; newstage=true; continue; };
								if( dpnt > pntmn ) { ori_stage++; newstage=true; continue; };
								int dp = (int)(-dpnt+pntmn)*10+1;
								int dt = (int)(-dtri+trimn)*10+1;
								// if(ipnt==18 && itri==72 && iori==276)
								//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg3_(ipnt/ANGLE_INCR+1,dp)
								//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg3_(itri/ANGLE_INCR+1,dt) << std::endl;
								if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbneg3_(ipnt/ANGLE_INCR+1,dp);
								if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg3_(itri/ANGLE_INCR+1,dt);
							}

							// surfdis3(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = d;
							// cbcount(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = tmpcbc;
							gdis(ipnt+1,itri+1,iori+1) = d;
							gcbc(ipnt+1,itri+1,iori+1) = tmpcbc;
							// d = sicsafe(t,p,sicaxis,cbcount);
							//std::cerr << "trial " << ipnt << " " << itri << " " << iori << " " << d << " " << tmpcbc << std::endl;
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							if(tmpcbc > cbmax) {
								cbmax = tmpcbc;
								mxiori = iori;
								mxipnt = ipnt;
								mxitri = itri;
								mxd = d;
							}
						}
						// utility_exit_with_message("testing");
					}
				}
				if(cbmax<0.00001) utility_exit_with_message("tri or pnt too large, no contacts!");
				std::cerr << "MAX3 " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << " " << mxd << std::endl;
				// std::cerr << "MAX " << mxipnt << " " << mxitri << " " << mxiori << " " << cbcount(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << " " << surfdis3(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << std::endl;
			}

		}

		//std::cerr << "main loop 2 over top 1000 3deg hits at 1deg" << std::endl;
		vector1<int> hitpnt,hittri,hitori;
		vector1<double> hitcbc;
		ANGLE_INCR = 1;
		if(ANGLE_INCR != 1) utility_exit_with_message("second ANGLE_INCR must be 1!!!");
		{
			// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
			std::cerr << "precomputing pent-pent and tri-tri interactions every 1°" << std::endl;
			{
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
					vector1<Vecf> ppnt,cbp; // precompute these
					get_pnt(ipnt,ppnt,cbp);
					vector1<Vecf> ppn2(ppnt),cb2(cbp);
					pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
					pnt1_to_pnt2( cb2,pntaxs_,pntax2_);				
					double cbcount = 0;
					double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntax2_-pntaxs_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntax2_-pntaxs_).normalized();
					pntmnpos1_[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
					pntdspos1_[ipnt/ANGLE_INCR+1] = -d;
					pntcbpos1_(ipnt/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbp,cb2);
					for(int i = 2; i <= 97; ++i) {
						for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1f*pntaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*pntax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						pntcbpos1_(ipnt/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
					vector1<Vecf> ppnt,cbp; // precompute these
					get_pnt(ipnt,ppnt,cbp);					
					vector1<Vecf> ppn2(ppnt),cb2(cbp);
					pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
					pnt1_to_pnt2( cb2,pntaxs_,pntax2_);				
					double cbcount = 0;
					double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntaxs_-pntax2_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntaxs_-pntax2_).normalized();
					pntmnneg1_[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
					pntdsneg1_[ipnt/ANGLE_INCR+1] = d;
					pntcbneg1_(ipnt/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbp,cb2);
					for(int i = 2; i <= 97; ++i) {
						for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1f*pntaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*pntax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						pntcbneg1_(ipnt/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif			
				for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
					vector1<Vecf> ptri,cbt; // precompute these
					get_tri(itri,ptri,cbt);
					vector1<Vecf> ptr2(ptri),cb2(cbt);
					tri1_to_tri2(ptr2,triaxs_,triax2_);
					tri1_to_tri2( cb2,triaxs_,triax2_);				
					double cbcount = 0;
					double const d = sicsafe(ptri,ptr2,cbt,cb2,(triax2_-triaxs_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triax2_-triaxs_).normalized();
					trimnpos1_[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
					tridspos1_[itri/ANGLE_INCR+1] = -d;
					tricbpos1_(itri/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbt,cb2);
					for(int i = 2; i <= 145; ++i) {
						for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1f*triaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*triax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						tricbpos1_(itri/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif			
				for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
					vector1<Vecf> ptri,cbt; // precompute these
					get_tri(itri,ptri,cbt);
					vector1<Vecf> ptr2(ptri),cb2(cbt);
					tri1_to_tri2(ptr2,triaxs_,triax2_);
					tri1_to_tri2( cb2,triaxs_,triax2_);				
					double cbcount = 0;
					double const d = sicsafe(ptri,ptr2,cbt,cb2,(triaxs_-triax2_).normalized(),cbcount,false);
					if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
					for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triaxs_-triax2_).normalized();
					trimnneg1_[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
					tridsneg1_[itri/ANGLE_INCR+1] = d;
					tricbneg1_(itri/ANGLE_INCR+1,1) = cbcount;
					// std::cerr << "compute CB" << std::endl;
					prune_cb_pairs_dis10(cbt,cb2);
					for(int i = 2; i <= 145; ++i) {
						for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1f*triaxs_;
						for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*triax2_;
						double cbc = 0.0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100.0) cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
						tricbneg1_(itri/ANGLE_INCR+1,i) = cbc;
						if(cbc==0) break;
					}
				}
				// {
				// 	using utility::file_basename;
				// 	std::cerr << "dumping 1D stats: " << option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_POS_1D.dat" << std::endl;
				// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_POS_1D.dat");
				// 		for(int i = 1; i <=	72; ++i) out << i << " " << pntdspos1_[i] << " " << pntcbpos1_(i,1) << std::endl;
				// 		out.close(); }
				// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_NEG_1D.dat");
				// 		for(int i = 1; i <=	72; ++i) out << i << " " << pntdsneg1_[i] << " " << pntcbneg1_(i,1) << std::endl;
				// 		out.close(); }
				// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[itrifile])+"_POS_1D.dat");
				// 		for(int i = 1; i <= 120; ++i) out << i << " " << tridspos1_[i] << " " << tricbpos1_(i,1) << std::endl;
				// 		out.close(); }
				// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[itrifile])+"_NEG_1D.dat");
				// 		for(int i = 1; i <= 120; ++i) out << i << " " << tridsneg1_[i] << " " << tricbneg1_(i,1) << std::endl;
				// 		out.close(); }
				// }
			
			}

			double top1000_3 = 0;
			{
				vector1<double> cbtmp;
				for(Size i = 0; i < gcbc.size(); ++i) if(gcbc[i] > 0) cbtmp.push_back(gcbc[i]);
				std::sort(cbtmp.begin(),cbtmp.end());
				top1000_3 = cbtmp[max(1,(int)cbtmp.size()-999)];
				std::cerr << "scanning top1000 with cbcount3 >= " << top1000_3 << std::endl;
			}
			// assuming ANGLE_INCR = 1!!!!!!!!!!!
			utility::vector1<vector1<int> > pntlmx,trilmx,orilmx;
			for(int ipnt = 0; ipnt < 72; ipnt+=3) {
				for(int itri = 0; itri < 120; itri+=3) {
					for(int iori = 0; iori < 360; iori+=3) {
						if( gcbc(ipnt+1,itri+1,iori+1) >= top1000_3) {
							vector1<int> pnt,tri,ori;
							for(int i = -1; i <= 1; ++i) pnt.push_back( (ipnt+i+ 72)% 72 );
							for(int j = -1; j <= 1; ++j) tri.push_back( (itri+j+120)%120 );
							for(int k = -1; k <= 1; ++k) ori.push_back( (iori+k+360)%360 );
							pntlmx.push_back(pnt);
							trilmx.push_back(tri);
							orilmx.push_back(ori);
						}
					}
				}
			}

			{ // loop 1deg lmx
				double max1 = 0;
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(Size ilmx = 1; ilmx <= pntlmx.size(); ++ilmx) {
					if( (ilmx-1)%100==0 && ilmx!=1) std::cerr << "	loop2 " << trifile << " " << (double(ilmx-1)/10) << "\% done, cbmax: " << max1 << std::endl;
					// std::cerr << "checking around local max # " << ilmx << std::endl;
					double cbmax = 0; int mxiori = 0, mxipnt = 0, mxitri = 0;
					for(vector1<int>::const_iterator pipnt = pntlmx[ilmx].begin(); pipnt != pntlmx[ilmx].end(); ++pipnt) {
						int ipnt = *pipnt;
						vector1<Vecf> pb,cbb; // precompute these
						get_pnt(ipnt,pb,cbb);
						for(vector1<int>::const_iterator pitri = trilmx[ilmx].begin(); pitri != trilmx[ilmx].end(); ++pitri) {
							int itri = *pitri;
							vector1<Vecf> pa,cba; // precompute these
							get_tri(itri,pa,cba);
							for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
								int iori = *piori;								
								if( gdis(ipnt+1,itri+1,iori+1) < -0.5 ) {								
									Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
									double tmpcbc;
									double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
									double theta = iori;
									double gamma = theta-alpha_;
									double x = d * sin(numeric::conversions::radians(gamma));
									double y = d * cos(numeric::conversions::radians(gamma));
									double w = x / sin(numeric::conversions::radians(alpha_));
									double z = x / tan(numeric::conversions::radians(alpha_));
									double dpnt = y+z;
									double dtri = w;
									double pntmn,trimn;
									if( w > 0 ) {
										pntmn = pntmnpos1_[ipnt/ANGLE_INCR+1];
										trimn = trimnpos1_[itri/ANGLE_INCR+1];
										if( dpnt < pntmn ) continue;
										if( dtri < trimn ) continue;								
										int dp = (int)(dpnt-pntmn)*10+1;
										int dt = (int)(dtri-trimn)*10+1;
										// if(ipnt==18 && itri==72 && iori==276)
										//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos1_(ipnt/ANGLE_INCR+1,dp)
										//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos1_(itri/ANGLE_INCR+1,dt) << std::endl;
										if( dp <= 97	) tmpcbc += option[tcdock::intra]()*pntcbpos1_(ipnt/ANGLE_INCR+1,dp);
										if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos1_(itri/ANGLE_INCR+1,dt);
										// std::cerr << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
									} else {
										pntmn = pntmnneg1_[ipnt/ANGLE_INCR+1];
										trimn = trimnneg1_[itri/ANGLE_INCR+1];
										if( dpnt > pntmn ) continue;
										if( dtri > trimn ) continue;
										int dp = (int)(-dpnt+pntmn)*10+1;
										int dt = (int)(-dtri+trimn)*10+1;
										// if(ipnt==18 && itri==72 && iori==276)
										//	std::cerr << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg1_(ipnt/ANGLE_INCR+1,dp)
										//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg1_(itri/ANGLE_INCR+1,dt) << std::endl;
										if( dp <= 97	) tmpcbc += option[tcdock::intra]()*pntcbneg1_(ipnt/ANGLE_INCR+1,dp);
										if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg1_(itri/ANGLE_INCR+1,dt);
									}
									gdis(ipnt+1,itri+1,iori+1) = d;
									gcbc(ipnt+1,itri+1,iori+1) = tmpcbc;
								}

								#ifdef USE_OPENMP
								#pragma omp critical
								#endif
								if(gcbc(ipnt+1,itri+1,iori+1) > cbmax) {
									cbmax = gcbc(ipnt+1,itri+1,iori+1);
									mxiori = iori;
									mxipnt = ipnt;
									mxitri = itri;
								}
							}

						}
					}
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					{
						hitpnt.push_back(mxipnt);
						hittri.push_back(mxitri);
						hitori.push_back(mxiori);
						hitcbc.push_back(cbmax);
						if(cbmax > max1) max1 = cbmax;
					}
					// std::cerr << "HIT " << ilmx << " " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << std::endl;
				}
			}

			double topX,max1;
			{
				vector1<double> hittmp = hitcbc;
				std::sort(hittmp.begin(),hittmp.end());
				max1	= hittmp[hittmp.size()  ];
				topX = hittmp[hittmp.size()-(option[tcdock::topx]()-1)];
				std::cerr << "topX " << topX << std::endl;
			}
			for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
				if(hitcbc[ihit]==max1)	std::cerr << "MAX1 " << hitpnt[ihit] << " " << hittri[ihit] << " " << hitori[ihit] << " " << hitcbc[ihit] << std::endl;
			}

			core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	// #ifdef USE_OPENMP
	// #pragma omp parallel for schedule(dynamic,1)
	// #endif
			vector1<std::pair<double,Size> > hits;
			for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
				if(hitcbc[ihit] >= topX ) {
					hits.push_back(std::pair<double,Size>(hitcbc[ihit],ihit));
				}
			}
			std::sort(hits.begin(),hits.end());
			for(Size ihiti = 1; ihiti <= hits.size(); ++ihiti) {
				Size ihit = hits[ihiti].second;
				using namespace ObjexxFCL::fmt;
				double dpnt,dtri,apnt=hitpnt[ihit],atri=hittri[ihit];
				double cbcsafe = calc_cbc(hitpnt[ihit],hittri[ihit],hitori[ihit],1,dpnt,dtri,true);
				if( cbcsafe < topX ) continue;
				std::cerr << "HIT " << F(7,3,hitcbc[ihit]) << " " << F(7,3,cbcsafe) << " " 
					<< utility::file_basename(pntfile) << " " << pnt_in_.n_residue() << " " << apnt << " " << dpnt << " "
				  << utility::file_basename(trifile) << " " << tri_in_.n_residue() << " " << atri << " " << dtri << " "
					<< std::endl;
				



				// std::cerr << "dumping grid 1 " << hitcbc[ihit] << std::endl;
				// vector1<numeric::xyzVector<int> > tocheck;
				// {
				// 	int ipnt = hitpnt[ihit];
				// 	int itri = hittri[ihit];
				// 	int iori = hitori[ihit];
				// 	for(int i = -8; i <= 8; ++i)
				// 	for(int j = -8; j <= 8; ++j)
				// 	for(int k = -2; k <= 2; ++k) {
				// 		tocheck.push_back( numeric::xyzVector<int>((ipnt+i+72)%72,(itri+j+120)%120,(iori+k+360)%360) );
				// 	}
				// }
				// ObjexxFCL::FArray3D<double> g11d(17,17,5,0.0),g11c(17,17,5,0.0);
				// // #ifdef USE_OPENMP
				// // #pragma omp parallel for schedule(dynamic,1)
				// // #endif					
				// for(int ic = 0; ic < tocheck.size(); ++ic)	{
				// 	int ipnt2 = tocheck[ic+1].x();
				// 	int itri2 = tocheck[ic+1].y();
				// 	int iori2 = tocheck[ic+1].z();
				// 	//std::cerr << "dumping grid " << ipnt2 << " " << itri2 << " " << iori2 << " " << ic << " " << g11d.size() << std::endl;
				// 	Pose p,t;
				// 	#ifdef USE_OPENMP
				// 	#pragma omp critical
				// 	#endif
				// 	{
				// 		p = pntinit;
				// 		t = triinit;
				// 	}
				// 	rot_pose(p,pntaxs_,(double)ipnt2);
				// 	rot_pose(t,triaxs_,(double)itri2);
				// 
				// 	vector1<Vecf> pb,cbb; // precompute these
				// 	for(int i = 1; i <= (int)p.n_residue(); ++i) {
				// 		cbb.push_back(Vecf(p.residue(i).xyz(2)));
				// 		int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
				// 		for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
				// 	}
				// 	vector1<Vecf> pa,cba; // precompute these
				// 	for(int i = 1; i <= (int)t.n_residue(); ++i) {
				// 		cba.push_back(Vecf(t.residue(i).xyz(2)));
				// 		int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
				// 		for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
				// 	}
				// 	Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori2) * Vecf(0,0,1)).normalized();
				// 	double tmpcbc;
				// 	double d = sicsafe(pb,pa,cbb,cba,sicaxis,tmpcbc);
				// 	double theta = iori2;
				// 	double gamma = theta-alpha_;
				// 	double x = d * sin(numeric::conversions::radians(gamma));
				// 	double y = d * cos(numeric::conversions::radians(gamma));
				// 	double w = x / sin(numeric::conversions::radians(alpha_));
				// 	double z = x / tan(numeric::conversions::radians(alpha_));
				// 	double dpnt = y+z;
				// 	double dtri = w;
				// 	double pntmn,trimn;
				// 	int dp,dt;
				// 	if( w > 0 ) {
				// 		pntmn = pntmnpos1_[ipnt2+1];
				// 		trimn = trimnpos1_[itri2+1];
				// 		dp = (int)((dpnt-pntmn)*10.0)+1;
				// 		dt = (int)((dtri-trimn)*10.0)+1;
				// 		if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbpos1_(ipnt2+1,dp);
				// 		if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos1_(itri2+1,dt);
				// 		if( dpnt < pntmn || dtri < trimn ) {
				// 			g11d[ic] = 1.0/0.0; g11c[ic] = 1.0/0.0;
				// 		} else {
				// 			g11d[ic] = d; g11c[ic] = tmpcbc;
				// 		}
				// 	} else {
				// 		pntmn = pntmnneg1_[ipnt2+1];
				// 		trimn = trimnneg1_[itri2+1];
				// 		dp = (int)((-dpnt+pntmn)*10.0)+1;
				// 		dt = (int)((-dtri+trimn)*10.0)+1;
				// 		if( dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbneg1_(ipnt2+1,dp);
				// 		if( dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg1_(itri2+1,dt);
				// 		if( dpnt > pntmn || dtri > trimn ) {
				// 			g11d[ic] = 1.0/0.0; g11c[ic] = 1.0/0.0;
				// 		} else {
				// 			g11d[ic] = d; g11c[ic] = tmpcbc;
				// 		}
				// 	}
				// 	if(ipnt2==hitpnt[ihit] && itri2==hittri[ihit]) {
				// 		if(iori2-hitori[ihit] == 2 || iori2-hitori[ihit] == -2) dump_pdb(ipnt2,itri2,iori2,"test_ori"+ObjexxFCL::string_of(iori2)+".pdb.gz",false);
				// 		if( iori2-hitori[ihit] == 0 ) dump_pdb(ipnt2,itri2,iori2,"test_ori"+ObjexxFCL::string_of(iori2)+".pdb.gz",true);								
				// 	}
				// 	// trans_pose(p,dpnt*pntaxs_);
				// 	// trans_pose(t,dtri*triaxs_);
				// 	// std::cerr << ipnt2 << " " << itri2 << " " << iori2 << " " << tmpcbc << " " << d << " " << w << " " << dp << " " << pntcbneg1_(ipnt2+1,dp) << " " << dt << " " << tricbneg1_(itri2+1,dt) << std::endl;
				// }
				// 
				// for(int i = 0; i < tocheck.size(); ++i) {
				// 	std::cout << g11d[i] << " " << g11c[i] << std::endl;
				// }
				// utility_exit_with_message("oarftn");
									
				// int ipnt = hitpnt[ihit];
				// int itri = hittri[ihit];
				// int iori = hitori[ihit];
				// double sizefac = tri_in_.n_residue() + pnt_in_.n_residue();
				// std::string fname;
				// fname	= ObjexxFCL::lead_zero_string_of(Size(hitcbc[ihit]/sizefac*1000000.0),7);
				// fname += "_" + utility::file_basename(pntfile);
				// fname += "_" + utility::file_basename(trifile);
				// fname += "_" + ObjexxFCL::lead_zero_string_of(hitcbc[ihit],4);
				// fname += "_" + ObjexxFCL::lead_zero_string_of(ipnt,3);
				// fname += "_" + ObjexxFCL::lead_zero_string_of(itri,3);
				// fname += "_" + ObjexxFCL::lead_zero_string_of(iori,3);
				// fname += ".pdb.gz";
				// 
				// std::cerr << "HIT " << " p" << ipnt << " t" << itri << " o" << iori << " c" << hitcbc[ihit] << " c" << hitcbc[ihit]/sizefac << std::endl;
				// dump_pdb(ipnt,itri,iori,fname);

			}
		}
	}

};

int main (int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	for(Size i = 2; i <= basic::options::option[in::file::s]().size(); ++i) {
		TCDock tcd(i,1);
		tcd.run();
	}
	std::cerr << "DONE testing refactoring 1-main" << std::endl;
}




//
//








