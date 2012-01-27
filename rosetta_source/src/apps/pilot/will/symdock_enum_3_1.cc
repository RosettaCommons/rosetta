#include <basic/options/keys/in.OptionKeys.gen.hh> 
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <devel/init.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
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
static basic::Tracer TR("symdock_enum");
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
		NEW_OPT( tcdock::contact_dis ,"max acceptable contact dis", 12 );
		NEW_OPT( tcdock::intra       ,"include intra 3-3 5-5 contacts", 1.0 );
		NEW_OPT( tcdock::topx        ,"output top X hits", 10 );
		NEW_OPT( tcdock::reverse     ,"rev.", false );	
	}
template<typename T> inline T sqr(T x) { return x*x; }
inline double sigmoid( double const & sqdist, double const & start, double const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		double dist = sqrt( sqdist );
		return (stop-dist)/(stop-start);
		//return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
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
struct Vecf2 {
	Vecf a,b;
	Vecf2() {}
	Vecf2(Vecf _a, Vecf _b) : a(_a),b(_b) {}
};
double sicfast( vector1<Vecf> pa, vector1<Vecf> pb, vector1<Vecf> const & cba, vector1<Vecf> const & cbb, Vecf ori, double & cbcount){ 
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

	// insert points into hashes
	int const xsize = xub-xlb+1;
	int const ysize = yub-ylb+1;
	ObjexxFCL::FArray2D<Vecf2> ha(xsize,ysize,Vecf2(Vecf(0,0,-9e9),Vecf(0,0,-9e9))),hb(xsize,ysize,Vecf2(Vecf(0,0,9e9),Vecf(0,0,9e9)));
	for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
		int const ix = (int)ceil(ia->x()/BIN)-xlb;
		int const iy = (int)ceil(ia->y()/BIN)-ylb;
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( ha(ix,iy).a.z() < ia->z() ) {
			ha(ix,iy).b = ha(ix,iy).a;
			ha(ix,iy).a = *ia;
		} else
		if( ha(ix,iy).b.z() < ia->z() ) {
			ha(ix,iy).b = *ia;			
		}
	}
	for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
		int const ix = (int)ceil(ib->x()/BIN)-xlb;
		int const iy = (int)ceil(ib->y()/BIN)-ylb;
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( hb(ix,iy).a.z() > ib->z() ) {
			hb(ix,iy).b = hb(ix,iy).a;
			hb(ix,iy).a = *ib;			
		} else 
		if( hb(ix,iy).b.z() > ib->z() ) {
			hb(ix,iy).b = *ib;
		}
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
					{
						double const xa = ha(i	,j	).a.x();
						double const ya = ha(i	,j	).a.y();
						double const xb = hb(i+k,j+l).a.x();
						double const yb = hb(i+k,j+l).a.y();
						double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
						if( d2 < CLASH_D2 ) {
							double dz = hb(i+k,j+l).a.z() - ha(i,j).a.z() - sqrt(CLASH_D2-d2);
							if( dz < mindis ) {
								mindis = dz;
								imna = i;
								jmna = j;
								imnb = i+k;
								jmnb = j+l;
							}
						}
					}
					{
						double const xa = ha(i	,j	).a.x();
						double const ya = ha(i	,j	).a.y();
						double const xb = hb(i+k,j+l).b.x();
						double const yb = hb(i+k,j+l).b.y();
						double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
						if( d2 < CLASH_D2 ) {
							double dz = hb(i+k,j+l).b.z() - ha(i,j).a.z() - sqrt(CLASH_D2-d2);
							if( dz < mindis ) {
								mindis = dz;
								imna = i;
								jmna = j;
								imnb = i+k;
								jmnb = j+l;
							}
						}
					}
					{
						double const xa = ha(i	,j	).b.x();
						double const ya = ha(i	,j	).b.y();
						double const xb = hb(i+k,j+l).a.x();
						double const yb = hb(i+k,j+l).a.y();
						double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
						if( d2 < CLASH_D2 ) {
							double dz = hb(i+k,j+l).a.z() - ha(i,j).b.z() - sqrt(CLASH_D2-d2);
							if( dz < mindis ) {
								mindis = dz;
								imna = i;
								jmna = j;
								imnb = i+k;
								jmnb = j+l;
							}
						}
					}
					{
						double const xa = ha(i	,j	).b.x();
						double const ya = ha(i	,j	).b.y();
						double const xb = hb(i+k,j+l).b.x();
						double const yb = hb(i+k,j+l).b.y();
						double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);
						if( d2 < CLASH_D2 ) {
							double dz = hb(i+k,j+l).b.z() - ha(i,j).b.z() - sqrt(CLASH_D2-d2);
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
	}


	cbcount = 0;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			double d2 = ib->distance_squared( (*ia) + (mindis*ori) );
			if( d2 < CONTACT_D2 ) {
				cbcount += sigmoid(d2, CLASH_D, CONTACT_D );
			}
		}
	}
	return mindis;
}
void tri1_to_tri2(vector1<Vecf> & pts, Vecf triaxs_, Vecf triax2_) { 
	Matf r = rotation_matrix_degrees( triaxs_.cross(triax2_), angle_degrees(triaxs_,Vecf(0,0,0),triax2_) );
	r = r * rotation_matrix_degrees(triaxs_,(double)60.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}
void pnt1_to_pnt2(vector1<Vecf> & pts, Vecf pntaxs_, Vecf pntax2_) { 
	Matf r = rotation_matrix_degrees( pntaxs_.cross(pntax2_), angle_degrees(pntaxs_,Vecf(0,0,0),pntax2_) );
	r = r * rotation_matrix_degrees(pntaxs_,(double)36.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}
struct LMAX {
	double score,radius;
	int ipnt,itri,iori;
	LMAX() :	score(0),radius(0),ipnt(0),itri(0),iori(0) {}
	LMAX(double _score, double _radius, int _ipnt, int _itri, int _iori) :
		score(_score),radius(_radius),ipnt(_ipnt),itri(_itri),iori(_iori) {}
};
int compareLMAX(const LMAX a,const LMAX b) {
	return a.score > b.score;
}
struct TCDock {
	vector1<double> pntmnpos_,trimnpos_,pntmnneg_,trimnneg_,pntdspos_,tridspos_,pntdsneg_,tridsneg_;
	ObjexxFCL::FArray2D<double> pntcbpos_,tricbpos_,pntcbneg_,tricbneg_;
	ObjexxFCL::FArray3D<double> gradii,gscore;	
	Vecf triaxs_,triax2_,pntaxs_,pntax2_;
	double alpha_,sin_alpha_,tan_alpha_;
	core::pose::Pose tri_in_,pnt_in_;
	vector1<Vecf> tri_pts_,tri_cbs_,pnt_pts_,pnt_cbs_;
	std::string trifile_,pntfile_;
	TCDock( Size itrifile, Size ipntfile ) : 
		pntmnpos_(72,0.0),trimnpos_(120,0.0),pntmnneg_(72,0.0),trimnneg_(120,0.0),
		pntdspos_(72,0.0),tridspos_(120,0.0),pntdsneg_(72,0.0),tridsneg_(120,0.0),
		pntcbpos_(72,97,0.0),tricbpos_(120,145,0.0),pntcbneg_(72,97,0.0),tricbneg_(120,145,0.0),
		gradii(72,120,360,-9e9),gscore(72,120,360,-9e9)
	{
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);
		core::import_pose::pose_from_pdb(tri_in_,*crs,option[in::file::s]()[itrifile]);
		core::import_pose::pose_from_pdb(pnt_in_,*crs,option[in::file::s]()[ipntfile]);
		trifile_ = utility::file_basename(option[in::file::s]()[itrifile]);
		pntfile_ = utility::file_basename(option[in::file::s]()[ipntfile]);
		make_trimer(tri_in_);
		make_pentamer(pnt_in_);		
		// set up geometry
		triaxs_ = Vecf( 0.000000, 0.000000,1.000000).normalized();
		triax2_ = Vecf(-0.333333,-0.577350,0.745356).normalized(); // 33.4458470159, 10.42594
		pntaxs_ = Vecf(-0.607226, 0.000000,0.794529).normalized();
		pntax2_ = Vecf(-0.491123,-0.850651,0.187593).normalized(); // 63.4311873349, 5.706642
		alpha_ = angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_);
		sin_alpha_ = sin(numeric::conversions::radians(alpha_));
		tan_alpha_ = tan(numeric::conversions::radians(alpha_));
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
	void precompute_intra() {
		double const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
		double const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
		double const CONTACT_D2 = sqr(CONTACT_D);
		double const CLASH_D2	 = sqr(CLASH_D);
		// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
		std::cerr << "precomputing pent-pent and tri-tri interactions every 1°" << std::endl;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(int ipnt = 0; ipnt < 72; ipnt+=1) {
			vector1<Vecf> ppnt,cbp; // precompute these
			get_pnt(ipnt,ppnt,cbp);
			vector1<Vecf> ppn2(ppnt),cb2(cbp);
			pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
			pnt1_to_pnt2( cb2,pntaxs_,pntax2_);				
			double cbcount = 0;
			double const d = sicfast(ppnt,ppn2,cbp,cb2,(pntax2_-pntaxs_).normalized(),cbcount);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntax2_-pntaxs_).normalized();
			pntmnpos_[ipnt+1] = -d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
			pntdspos_[ipnt+1] = -d;
			pntcbpos_(ipnt+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			prune_cb_pairs_dis10(cbp,cb2);
			for(int i = 2; i <= 97; ++i) {
				for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1f*pntaxs_;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*pntax2_;
				double cbc = 0.0;
				for(Size j = 1; j <= cbp.size(); ++j) 
					if(cbp[j].distance_squared(cb2[j]) < 100.0) 
						cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
				pntcbpos_(ipnt+1,i) = cbc;
				if(cbc==0) break;
			}
		}
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(int ipnt = 0; ipnt < 72; ipnt+=1) {
			vector1<Vecf> ppnt,cbp; // precompute these
			get_pnt(ipnt,ppnt,cbp);					
			vector1<Vecf> ppn2(ppnt),cb2(cbp);
			pnt1_to_pnt2(ppn2,pntaxs_,pntax2_);
			pnt1_to_pnt2( cb2,pntaxs_,pntax2_);				
			double cbcount = 0;
			double const d = sicfast(ppnt,ppn2,cbp,cb2,(pntaxs_-pntax2_).normalized(),cbcount);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntaxs_-pntax2_).normalized();
			pntmnneg_[ipnt+1] = d/2.0/sin( angle_radians(pntax2_,Vecf(0,0,0),pntaxs_)/2.0 );
			pntdsneg_[ipnt+1] = d;
			pntcbneg_(ipnt+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			prune_cb_pairs_dis10(cbp,cb2);
			for(int i = 2; i <= 97; ++i) {
				for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1f*pntaxs_;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*pntax2_;
				double cbc = 0.0; 
				for(Size j = 1; j <= cbp.size(); ++j) 
					if(cbp[j].distance_squared(cb2[j]) < 100.0) 
						cbc += sigmoid(cbp[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
				pntcbneg_(ipnt+1,i) = cbc;
				if(cbc==0) break;
			}
		}
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif			
		for(int itri = 0; itri < 120; itri+=1) {
			vector1<Vecf> ptri,cbt; // precompute these
			get_tri(itri,ptri,cbt);
			vector1<Vecf> ptr2(ptri),cb2(cbt);
			tri1_to_tri2(ptr2,triaxs_,triax2_);
			tri1_to_tri2( cb2,triaxs_,triax2_);				
			double cbcount = 0;
			double const d = sicfast(ptri,ptr2,cbt,cb2,(triax2_-triaxs_).normalized(),cbcount);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triax2_-triaxs_).normalized();
			trimnpos_[itri+1] = -d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
			tridspos_[itri+1] = -d;
			tricbpos_(itri+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			prune_cb_pairs_dis10(cbt,cb2);
			for(int i = 2; i <= 145; ++i) {
				for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1f*triaxs_;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1f*triax2_;
				double cbc = 0.0; 
				for(Size j = 1; j <= cbt.size(); ++j) 
					if(cbt[j].distance_squared(cb2[j]) < 100.0) 
						cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
				tricbpos_(itri+1,i) = cbc;
				if(cbc==0) break;
			}
		}
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif			
		for(int itri = 0; itri < 120; itri+=1) {
			vector1<Vecf> ptri,cbt; // precompute these
			get_tri(itri,ptri,cbt);
			vector1<Vecf> ptr2(ptri),cb2(cbt);
			tri1_to_tri2(ptr2,triaxs_,triax2_);
			tri1_to_tri2( cb2,triaxs_,triax2_);				
			double cbcount = 0;
			double const d = sicfast(ptri,ptr2,cbt,cb2,(triaxs_-triax2_).normalized(),cbcount);
			if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
			for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triaxs_-triax2_).normalized();
			trimnneg_[itri+1] = d/2.0/sin( angle_radians(triax2_,Vecf(0,0,0),triaxs_)/2.0 );
			tridsneg_[itri+1] = d;
			tricbneg_(itri+1,1) = cbcount;
			// std::cerr << "compute CB" << std::endl;
			prune_cb_pairs_dis10(cbt,cb2);
			for(int i = 2; i <= 145; ++i) {
				for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1f*triaxs_;
				for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1f*triax2_;
				double cbc = 0.0; 
				for(Size j = 1; j <= cbt.size(); ++j) 
					if(cbt[j].distance_squared(cb2[j]) < 100.0) 
						cbc += sigmoid(cbt[j].distance_squared(cb2[j]), CLASH_D, CONTACT_D );
				tricbneg_(itri+1,i) = cbc;
				if(cbc==0) break;
			}
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
		double const d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
		if(d > 0) utility_exit_with_message("ZERO!!");
		double const theta=aori, gamma=numeric::conversions::radians(theta-alpha_);
		double const sin_gamma=sin(gamma), cos_gamma=cos(gamma), x=d*sin_gamma, y=d*cos_gamma, w=x/sin_alpha_, z=x/tan_alpha_;
		double const dpnt = y+z, dtri = w;
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
		}
		if(sym) core::pose::symmetry::make_symmetric_pose(symm);				
		core::io::pdb::dump_pdb(symm,option[out::file::o]()+"/"+fname);
		
	}
	void dump_onecomp() {
		using	namespace	basic::options;
		using	namespace	basic::options::OptionKeys;		
		using utility::file_basename;
		std::cerr << "dumping 1D stats: " << option[out::file::o]()+"/"+pntfile_+"_POS_1D.dat" << std::endl;
		{ utility::io::ozstream out(option[out::file::o]()+"/"+pntfile_+"_POS_1D.dat");
			for(int i = 1; i <=	72; ++i) out << i << " " << pntdspos_[i] << " " << pntcbpos_(i,1) << std::endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+pntfile_+"_NEG_1D.dat");
			for(int i = 1; i <=	72; ++i) out << i << " " << pntdsneg_[i] << " " << pntcbneg_(i,1) << std::endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+trifile_+"_POS_1D.dat");
			for(int i = 1; i <= 120; ++i) out << i << " " << tridspos_[i] << " " << tricbpos_(i,1) << std::endl;
			out.close(); }
		{ utility::io::ozstream out(option[out::file::o]()+"/"+trifile_+"_NEG_1D.dat");
			for(int i = 1; i <= 120; ++i) out << i << " " << tridsneg_[i] << " " << tricbneg_(i,1) << std::endl;
			out.close(); }
	}
	double calc_cbc(int ipnt, int itri, int iori, double & dpnt, double & dtri, double & icbc, double & pntcbc, double & tricbc, bool cache=false) {
		if(!cache || gscore(ipnt+1,itri+1,iori+1) == -9e9 ) {
			using basic::options::option;
			using namespace basic::options::OptionKeys;
			icbc=0; pntcbc=0; tricbc=0;
			vector1<Vecf> pb,cbb, pa,cba; get_tri(itri,pa,cba); get_pnt(ipnt,pb,cbb);		
			Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
			double const d = sicfast(pb,pa,cbb,cba,sicaxis,icbc);
			if(d > 0) utility_exit_with_message("ZERO!!");
			double const theta=iori, gamma=numeric::conversions::radians(theta-alpha_);
			double const sin_gamma=sin(gamma), cos_gamma=cos(gamma), x=d*sin_gamma, y=d*cos_gamma, w=x/sin_alpha_, z=x/tan_alpha_;
			dpnt = y+z;  dtri = w;
			if( w > 0 ) {
				double const pntmn = pntmnpos_[ipnt+1];
				double const trimn = trimnpos_[itri+1];
				if( dtri < trimn || dpnt < pntmn ) { 
					double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
					double const dmntri = trimn*sin_alpha_/sin_gamma;
					gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
					gscore(ipnt+1,itri+1,iori+1) = 0;
					return 0;
				}
				int dp = (int)(dpnt-pntmn)*10+1;
				int dt = (int)(dtri-trimn)*10+1;
				if( 0 < dp && dp <=  97 ) tricbc = pntcbpos_(ipnt+1,dp);
				if( 0 < dt && dt <= 145 ) pntcbc = tricbpos_(itri+1,dt);
			} else {
				double const pntmn = pntmnneg_[ipnt+1];
				double const trimn = trimnneg_[itri+1];
				if( dtri > trimn || dpnt > pntmn ) { 
					double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
					double const dmntri = trimn*sin_alpha_/sin_gamma;
					gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
					gscore(ipnt+1,itri+1,iori+1) = 0;
					return 0;
				}
				int dp = (int)(-dpnt+pntmn)*10+1;
				int dt = (int)(-dtri+trimn)*10+1;
				if( 0 < dp && dp <=  97 ) tricbc = pntcbneg_(ipnt+1,dp);
				if( 0 < dt && dt <= 145 ) pntcbc = tricbneg_(itri+1,dt);
			}
			// #ifdef USE_OPENMP
			// #pragma omp critical
			// #endif
			// {
				gscore(ipnt+1,itri+1,iori+1) = icbc+option[tcdock::intra]()*(pntcbc+tricbc);
				gradii(ipnt+1,itri+1,iori+1) = d;
			// }
		}
		return gscore(ipnt+1,itri+1,iori+1);
	}
	void run() {
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		using namespace core::id;
		using numeric::conversions::radians;
		double const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
		double const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
		double const CONTACT_D2 = sqr(CONTACT_D);
		double const CLASH_D2	= sqr(CLASH_D);
		Pose const triinit(tri_in_);
		Pose const pntinit(pnt_in_);
		precompute_intra();
		std::cerr << "main loop 1 over ipnt, itri, iori every 3 degrees" << std::endl; { // 3deg loop 
			double mxd = 0;
			double cbmax = 0; int mxiori = 0, mxipnt = 0, mxitri = 0;
			for(int ipnt = 0; ipnt < 72; ipnt+=3) {
				if(ipnt%9==0 && ipnt!=0) std::cerr << "	 loop1 " << trifile_ << " " << (100*ipnt)/72 << "\% done, cbmax: " << cbmax << std::endl;
				vector1<Vecf> pb,cbb; // precompute these
				get_pnt(ipnt,pb,cbb);
				#ifdef USE_OPENMP
				#pragma omp parallel for schedule(dynamic,1)
				#endif
				for(int itri = 0; itri < 120; itri+=3) {
					vector1<Vecf> pa,cba; // precompute these
					get_tri(itri,pa,cba);
					int iori = -1, ori_stage = 1;
					bool newstage = true;
					while(ori_stage < 5) {
						if(newstage) {
							if( ori_stage == 1 || ori_stage == 2 ) iori = (int)( 90.0+double(3)/2.0+angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_));
							if( ori_stage == 3 || ori_stage == 4 ) iori = (int)(270.0+double(3)/2.0+angle_degrees(triaxs_,Vecf(0,0,0),pntaxs_));
							iori = (iori / 3) * 3; // round to closest multiple of angle incr
							if( ori_stage == 2 || ori_stage == 4 ) iori -= 3;
							newstage = false;
						} else {
							if( ori_stage == 1 || ori_stage == 3 ) iori += 3;
							if( ori_stage == 2 || ori_stage == 4 ) iori -= 3;
						}
						Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
						double tmpcbc;
						double const d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
						if(d > 0) utility_exit_with_message("ZERO!!");
						double const theta=iori, gamma=numeric::conversions::radians(theta-alpha_);
						double const sin_gamma=sin(gamma), cos_gamma=cos(gamma), x=d*sin_gamma, y=d*cos_gamma, w=x/sin_alpha_, z=x/tan_alpha_;
						double const dpnt = y+z, dtri = w;
						if( w > 0 ) {							
							double const pntmn = pntmnpos_[ipnt+1], trimn = trimnpos_[itri+1];
							if( dtri < trimn || dpnt < pntmn ) { 
								double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
								double const dmntri = trimn*sin_alpha_/sin_gamma;
								gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
								gscore(ipnt+1,itri+1,iori+1) = 0;
								ori_stage++;
								newstage=true;
								continue;
							}
							int dp = (int)(dpnt-pntmn)*10+1, dt = (int)(dtri-trimn)*10+1;
							if( 0 < dp && dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbpos_(ipnt+1,dp);
							if( 0 < dt && dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos_(itri+1,dt);
						} else {
							double const pntmn = pntmnneg_[ipnt+1], trimn = trimnneg_[itri+1];
							if( dtri > trimn || dpnt > pntmn ) {
								double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
								double const dmntri = trimn*sin_alpha_/sin_gamma;
								gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
								gscore(ipnt+1,itri+1,iori+1) = 0;
								ori_stage++;
								newstage=true;
								continue;
							}
							int dp = (int)(-dpnt+pntmn)*10+1, dt = (int)(-dtri+trimn)*10+1;
							if( 0 < dp && dp <=  97 ) tmpcbc += option[tcdock::intra]()*pntcbneg_(ipnt+1,dp);
							if( 0 < dt && dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg_(itri+1,dt);
						}
						gradii(ipnt+1,itri+1,iori+1) = d;
						gscore(ipnt+1,itri+1,iori+1) = tmpcbc;
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
				}
			}
			if(cbmax<0.00001) utility_exit_with_message("tri or pnt too large, no contacts!");
			std::cerr << "MAX3 " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << " " << mxd << std::endl;
		}
		double topX=0,max1=0;
		utility::vector1<vector1<int> > pntlmx,trilmx,orilmx; { // set up work for main loop 2 
			double top1000_3 = 0;
			vector1<double> cbtmp;
			for(Size i = 0; i < gscore.size(); ++i) if(gscore[i] > 0) cbtmp.push_back(gscore[i]);
			std::sort(cbtmp.begin(),cbtmp.end());
			top1000_3 = cbtmp[max(1,(int)cbtmp.size()-999)];
			std::cerr << "scanning top1000 with cbcount3 >= " << top1000_3 << std::endl;
			for(int ipnt = 0; ipnt < 72; ipnt+=3) {
				for(int itri = 0; itri < 120; itri+=3) {
					for(int iori = 0; iori < 360; iori+=3) {
						if( gscore(ipnt+1,itri+1,iori+1) >= top1000_3) {
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
		}
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(Size ilmx = 1; ilmx <= pntlmx.size(); ++ilmx)  {       //  MAIN LOOP 2                                      
			if( (ilmx-1)%100==0 && ilmx!=1) std::cerr << "	loop2 " << trifile_ 
			         << " " << (double(ilmx-1)/10) << "\% done, cbmax: " << max1 << std::endl;
			double cbmax = 0; int mxiori = 0, mxipnt = 0, mxitri = 0;
			for(vector1<int>::const_iterator pipnt = pntlmx[ilmx].begin(); pipnt != pntlmx[ilmx].end(); ++pipnt) {
				int ipnt = *pipnt; vector1<Vecf> pb,cbb; get_pnt(ipnt,pb,cbb);
				for(vector1<int>::const_iterator pitri = trilmx[ilmx].begin(); pitri != trilmx[ilmx].end(); ++pitri) {
					int itri = *pitri; vector1<Vecf> pa,cba; get_tri(itri,pa,cba);
					for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
						int iori = *piori;								
						// if( gradii(ipnt+1,itri+1,iori+1) == -9e9 ) {
							Vecf sicaxis = (rotation_matrix_degrees(-Vecf(0,1,0),(double)iori) * Vecf(0,0,1)).normalized();
							double tmpcbc;
							double const d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
							if(d > 0) utility_exit_with_message("ZERO!!");
							double const theta=iori, gamma=numeric::conversions::radians(theta-alpha_);
							double const sin_gamma=sin(gamma), cos_gamma=cos(gamma), x=d*sin_gamma, y=d*cos_gamma, w=x/sin_alpha_, z=x/tan_alpha_;
							double const dpnt = y+z, dtri = w;
							if( w > 0 ) {
								double const pntmn = pntmnpos_[ipnt+1];
								double const trimn = trimnpos_[itri+1];
								if( dpnt < pntmn || dtri < trimn ) {
									double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
									double const dmntri = trimn*sin_alpha_/sin_gamma;
									gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
									gscore(ipnt+1,itri+1,iori+1) = 0;
									continue;								
								}
								int dp = (int)(dpnt-pntmn)*10+1;
								int dt = (int)(dtri-trimn)*10+1;
								if( 0 < dp && dp <= 97  ) tmpcbc += option[tcdock::intra]()*pntcbpos_(ipnt+1,dp);
								if( 0 < dt && dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbpos_(itri+1,dt);
							} else {
								double const pntmn = pntmnneg_[ipnt+1];
								double const trimn = trimnneg_[itri+1];
								if( dpnt > pntmn || dtri > trimn ){
									double const dmnpnt = pntmn/(cos_gamma+sin_gamma/tan_alpha_);
									double const dmntri = trimn*sin_alpha_/sin_gamma;
									gradii(ipnt+1,itri+1,iori+1) = min(dmnpnt,dmntri);
									gscore(ipnt+1,itri+1,iori+1) = 0;
									continue;
								}
								int dp = (int)(-dpnt+pntmn)*10+1;
								int dt = (int)(-dtri+trimn)*10+1;
								if( 0 < dp && dp <= 97  ) tmpcbc += option[tcdock::intra]()*pntcbneg_(ipnt+1,dp);
								if( 0 < dt && dt <= 145 ) tmpcbc += option[tcdock::intra]()*tricbneg_(itri+1,dt);
							}
							gradii(ipnt+1,itri+1,iori+1) = d;
							gscore(ipnt+1,itri+1,iori+1) = tmpcbc;
						//}

						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						if(gscore(ipnt+1,itri+1,iori+1) > cbmax) {
							cbmax = gscore(ipnt+1,itri+1,iori+1);
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
			if(cbmax > max1) max1 = cbmax;
		}
				vector1<LMAX> local_maxima; {                        // get local radial disp maxima (minima, really)     
			double mn=9e9,mx=-9e9;
			LMAX highscore; highscore.score = -9e9;
			// vector1<int> ehist(26,0);
			for(int i = 1; i <= gradii.size1(); ++i){
				for(int j = 1; j <= gradii.size2(); ++j){
					for(int k = 1; k <= gradii.size3(); ++k){
						double const val = gradii(i,j,k);
						if( val < -9e8 ) continue;
						double nbmax = -9e9;
						int nedge = 0;
						for(int di = -1; di <= 1; ++di){
							for(int dj = -1; dj <= 1; ++dj){
								for(int dk = -1; dk <= 1; ++dk){
									if(di==0 && dj==0 && dk==0) continue;
									int i2 = (i+di+gradii.size1()-1)%gradii.size1()+1;
									int j2 = (j+dj+gradii.size2()-1)%gradii.size2()+1;
									int k2 = (k+dk+gradii.size3()-1)%gradii.size3()+1;
									double const nbval = gradii(i2,j2,k2);
									nbmax = max(nbmax,nbval);
								}
							}
						}
						if( nbmax != -9e9 && val >= nbmax ) {
							local_maxima.push_back( LMAX(gscore(i,j,k),gradii(i,j,k),i-1,j-1,k-1) );
							mx = max(val,mx);
							mn = min(val,mn);
						}
					}
				}
			}
			std::sort(local_maxima.begin(),local_maxima.end(),compareLMAX);
			std::cerr << "N local max (tri-pent disp.) " << local_maxima.size() << std::endl;					
		}
		for(Size ilm = 1; ilm <= min(local_maxima.size(),(Size)option[tcdock::topx]()); ++ilm) { // dump top hit info 
			LMAX const & h(local_maxima[ilm]);
			double dpnt,dtri,icbc,tricbc,pntcbc;
			ObjexxFCL::FArray3D<double> grid(35,35,35,0.0);
			#ifdef USE_OPENMP
			#pragma omp parallel for schedule(dynamic,1)
			#endif			
			for(int di = -17; di <= 17; ++di) {
				for(int dj = -17; dj <= 17; ++dj) {
					for(int dk = -17; dk <= 17; ++dk) {
						int i = (h.ipnt+di+gscore.size1())%gscore.size1();
						int j = (h.itri+dj+gscore.size2())%gscore.size2();
						int k = (h.iori+dk+gscore.size3())%gscore.size3();
						calc_cbc(i,j,k,dpnt,dtri,icbc,pntcbc,tricbc,true);
						grid(di+18,dj+18,dk+18) = gradii(i+1,j+1,k+1);
					}
				}
			}
			std::ofstream o(("out/grid_"+ObjexxFCL::string_of(ilm)+".dat").c_str());
			for(int i = 1; i <= grid.size1(); ++i) {
				for(int j = 1; j <= grid.size2(); ++j) {
					for(int k = 1; k <= grid.size3(); ++k) {
						o << grid(i,j,k) << std::endl;
					}
				}
			}
			o.close();
			double cbc = calc_cbc(h.ipnt,h.itri,h.iori,dpnt,dtri,icbc,pntcbc,tricbc,false);
			std::cerr << "HIT " 
			          << F(7,3,h.score) << " "
			          << F(7,3,cbc) << " " 
			          << F(7,3,icbc) << " " 
			          << F(7,3,pntcbc) << " " 
			          << F(7,3,tricbc) << " " 
			          << pntfile_ << " " 
			          << I(4,pnt_in_.n_residue()) << " " 
			          << I(2,h.ipnt) << " " 
			          << F(7,3,dpnt) << " "
			          << trifile_ << " " 
			          << I(4,tri_in_.n_residue()) << " " 
			          << I(3,h.itri) << " " 
			          << F(7,3,dtri) << " "
			          << I(3,h.iori) << " " 
			          << std::endl;
			dump_pdb(h.ipnt,h.itri,h.iori,"test_"+ObjexxFCL::string_of(cbc)+".pdb");
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
	std::cerr << "DONE testing grid" << std::endl;
}
