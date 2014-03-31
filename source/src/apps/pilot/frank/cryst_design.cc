/// @file
/// @brief


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>
#include <queue>


using namespace basic;
using namespace core;
using namespace core::pose;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

static basic::Tracer TR("cryst.design");

OPT_1GRP_KEY(String, crystdock, spacegroup)
OPT_1GRP_KEY(Real, crystdock, A)
OPT_1GRP_KEY(Real, crystdock, B)
OPT_1GRP_KEY(Real, crystdock, C)
OPT_1GRP_KEY(Real, crystdock, alpha)
OPT_1GRP_KEY(Real, crystdock, beta)
OPT_1GRP_KEY(Real, crystdock, gamma)
OPT_1GRP_KEY(Real, crystdock, maxclash)
OPT_1GRP_KEY(Real, crystdock, mininterface)
OPT_1GRP_KEY(Real, crystdock, trans_step)
OPT_1GRP_KEY(Real, crystdock, rot_step)
OPT_1GRP_KEY(Integer, crystdock, nmodels)
OPT_1GRP_KEY(Integer, crystdock, rotnum)
OPT_1GRP_KEY(Integer, crystdock, symnum)
OPT_1GRP_KEY(Boolean, crystdock, debug)
OPT_1GRP_KEY(Boolean, crystdock, debug_exact)

////////////////////////////////////////////////
// helper functions
inline int pos_mod(int x,int y) {
	int r=x%y; if (r<0) r+=y;
	return r;
}
inline Real pos_mod(Real x,Real y) {
	Real r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y; if (r<-y/2) r+=y;if (r>=y/2) r-=y;
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if (r<-0.5*y) r+=y;if (r>=0.5*y) r-=y;
	return r;
}

inline double sign(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

/// backward interpolation (unrotated->rotated)
///    used to generate grids
template <class S>
core::Real interp_linear(
		ObjexxFCL::FArray3D< S > const & data ,
		numeric::xyzVector< core::Real > const & idxX ) {

	numeric::xyzVector<int> pt000, pt111;
	numeric::xyzVector< core::Real > fpart,neg_fpart;
	numeric::xyzVector<int> srcgrid(data.u1(),data.u2(),data.u3());

	// 0->1-based indexing
	pt000[0] = (int)(std::floor(idxX[0]+1e-6))+1;
	pt000[1] = (int)(std::floor(idxX[1]+1e-6))+1;
	pt000[2] = (int)(std::floor(idxX[2]+1e-6))+1;

	// bound check
	if (pt000[0] < -srcgrid[0]/2+1 || pt000[0] >= srcgrid[0]/2+1) return 0.0;
	if (pt000[1] < -srcgrid[1]/2+1 || pt000[1] >= srcgrid[1]/2+1) return 0.0;
	if (pt000[2] < -srcgrid[2]/2+1 || pt000[2] >= srcgrid[2]/2+1) return 0.0;
	if (pt000[0] <= 0 ) pt000[0] += srcgrid[0];
	if (pt000[1] <= 0 ) pt000[1] += srcgrid[1];
	if (pt000[2] <= 0 ) pt000[2] += srcgrid[2];

	// interpolation coeffs
	fpart[0] = idxX[0]-std::floor(idxX[0]); neg_fpart[0] = 1-fpart[0];
	fpart[1] = idxX[1]-std::floor(idxX[1]); neg_fpart[1] = 1-fpart[1];
	fpart[2] = idxX[2]-std::floor(idxX[2]); neg_fpart[2] = 1-fpart[2];

	S retval = (S)0.0;

	pt111[0] = (pt000[0]+1); if (pt111[0]>srcgrid[0]) pt111[0]=1;
	pt111[1] = (pt000[1]+1); if (pt111[1]>srcgrid[1]) pt111[1]=1;
	pt111[2] = (pt000[2]+1); if (pt111[2]>srcgrid[2]) pt111[2]=1;

	retval += neg_fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt000[0],pt000[1],pt000[2]);
	retval += neg_fpart[0]*neg_fpart[1]*    fpart[2] * data(pt000[0],pt000[1],pt111[2]);
	retval += neg_fpart[0]*    fpart[1]*neg_fpart[2] * data(pt000[0],pt111[1],pt000[2]);
	retval += neg_fpart[0]*    fpart[1]*    fpart[2] * data(pt000[0],pt111[1],pt111[2]);
	retval +=     fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt111[0],pt000[1],pt000[2]);
	retval +=     fpart[0]*neg_fpart[1]*    fpart[2] * data(pt111[0],pt000[1],pt111[2]);
	retval +=     fpart[0]*    fpart[1]*neg_fpart[2] * data(pt111[0],pt111[1],pt000[2]);
	retval +=     fpart[0]*    fpart[1]*    fpart[2] * data(pt111[0],pt111[1],pt111[2]);

	return (core::Real)retval;
}


// rotation stuff
void
euler2rot( core::Real a, core::Real b, core::Real g, numeric::xyzMatrix<Real> &R) {
	R.xx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + cos(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*sin(DEG2RAD*g) + sin(DEG2RAD*a)*cos(DEG2RAD*g);
	R.xz() =  sin(DEG2RAD*b)*sin(DEG2RAD*g);
	R.yx() = -sin(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - cos(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yy() =  cos(DEG2RAD*a)*cos(DEG2RAD*b)*cos(DEG2RAD*g) - sin(DEG2RAD*a)*sin(DEG2RAD*g);
	R.yz() =  sin(DEG2RAD*b)*cos(DEG2RAD*g);
	R.zx() =  sin(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zy() = -cos(DEG2RAD*a)*sin(DEG2RAD*b);
	R.zz() =  cos(DEG2RAD*b);
}

// angle of rotation from rotation matrix
// an abomination of 0- and 1- indexing
core::Real
R2ang(numeric::xyzMatrix<Real> R) {
	Real q0 = ( R(1,1) + R(2,2) + R(3,3) + 1.0) / 4.0;
	Real q1 = ( R(1,1) - R(2,2) - R(3,3) + 1.0) / 4.0;
	Real q2 = (-R(1,1) + R(2,2) - R(3,3) + 1.0) / 4.0;
	Real q3 = (-R(1,1) - R(2,2) + R(3,3) + 1.0) / 4.0;

	if (q0 < 0.0) q0 = 0.0f; else q0 = sqrt(q0);
	if (q1 < 0.0) q1 = 0.0f; else q1 = sqrt(q1);
	if (q2 < 0.0) q2 = 0.0f; else q2 = sqrt(q2);
	if (q3 < 0.0) q3 = 0.0f; else q3 = sqrt(q3);

	if(q0 >= q1 && q0 >= q2 && q0 >= q3) {
		q1 *= sign(R(3,2) - R(2,3));
		q2 *= sign(R(1,3) - R(3,1));
		q3 *= sign(R(2,1) - R(1,2));
	} else if(q1 >= q0 && q1 >= q2 && q1 >= q3) {
		q0 *= sign(R(3,2) - R(2,3));
		q2 *= sign(R(2,1) + R(1,2));
		q3 *= sign(R(1,3) + R(3,1));
	} else if(q2 >= q0 && q2 >= q1 && q2 >= q3) {
		q0 *= sign(R(1,3) - R(3,1));
		q1 *= sign(R(2,1) + R(1,2));
		q3 *= sign(R(3,2) + R(2,3));
	} else if(q3 >= q0 && q3 >= q1 && q3 >= q2) {
		q0 *= sign(R(2,1) - R(1,2));
		q1 *= sign(R(3,1) + R(1,3));
		q2 *= sign(R(3,2) + R(2,3));
	}
	Real r = norm4(q0, q1, q2, q3);
	q0 /= r; q1 /= r; q2 /= r; q3 /= r;
	Real angle = 2*RAD2DEG*std::fabs(acos(q0));

	return (angle);
}


////////////////////////////////////////////////
// Uniformly sample rotation space
//   * use idea from A. Yershova, et al, International Journal of Robotics Research, 2009.
//       to get SO3 samples from S1+S2 uniform sampling
//   * use icosahedral embedding to generate approximately uniform coverage of S2
struct Quat {
	Real x_,y_,z_,w_;
	Quat() {
		x_=0; y_=0; z_=0; w_=1;
	}
	Quat( Real x_in,Real y_in,Real z_in,Real w_in) {
		x_=x_in; y_=y_in; z_=z_in; w_=w_in;
	}

	// construct from a rotation matrix
	// should be stable
	Quat( numeric::xyzMatrix<Real> const& R) {
		Real S;
		if (R.xx() + R.yy() + R.zz() > 0) {
			S = 0.5 / sqrt(R.xx() + R.yy() + R.zz());
			w_ = 0.25 / S;
			x_ = ( R.zy() - R.yz() ) * S;
			y_ = ( R.xz() - R.zx() ) * S;
			z_ = ( R.yx() - R.xy() ) * S;
		} else if ( R.xx() > R.yy() && R.xx() > R.zz() )  {
			S  = sqrt( 1.0 + R.xx() - R.yy() - R.zz() ) * 2;
			x_ = 0.25 * S;
			y_ = (R.yx() + R.xy() ) / S;
			z_ = (R.zx() + R.xz() ) / S;
			w_ = (R.zy() - R.yz() ) / S;
		} else if ( R.yy() > R.zz() ) {
			S  = sqrt( 1.0 - R.xx() + R.yy() - R.zz() ) * 2;
			x_ = (R.yx() + R.xy() ) / S;
			y_ = 0.25 * S;
			z_ = (R.zy() + R.yz() ) / S;
			w_ = (R.xz() - R.zx() ) / S;
		} else {
			S  = sqrt( 1.0 - R.xx() - R.yy() + R.zz() ) * 2;
			x_ = (R.xz() + R.zx() ) / S;
			y_ = (R.zy() + R.yz() ) / S;
			z_ = 0.25 * S;
			w_ = (R.yx() - R.xy() ) / S;
		}
	}

	numeric::xyzMatrix<Real> asR() const {
		if (1-w_*w_ < 1e-6)
			return numeric::xyzMatrix<core::Real>::rows(1,0,0, 0,1,0, 0,0,1);

		Real xx = x_*x_, xy = x_*y_, xz = x_*z_, xw = x_*w_;
		Real yy = y_*y_, yz = y_*z_, yw = y_*w_;
		Real zz = z_*z_, zw = z_*w_;
		//Real ww = w_*w_;

		return numeric::xyzMatrix<core::Real>::rows(
			1 - 2 * ( yy+zz ) ,     2 * ( xy-zw ) ,     2 * ( xz+yw ) ,
			    2 * ( xy+zw ) , 1 - 2 * ( xx+zz ) ,     2 * ( yz-xw ) ,
			    2 * ( xz-yw ) ,     2 * ( yz+xw ) , 1 - 2 * ( xx+yy ) );
	}
};


class UniformRotationSampler {
private:
	utility::vector1< Quat > rotlist_;
	Real theta_;

public:
	Size
	nrots() const { return rotlist_.size(); }

	void get(Size ii, numeric::xyzMatrix<Real> &R) const {
		if (ii==0) {
			R = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
		} else {
			R = rotlist_[ii].asR();
		}
	}

	void
	generateIcosahedralSamples(utility::vector1<numeric::xyzVector<Real> > &ico, Size nsub) {
		using numeric::constants::d::pi;
		Real theta = 26.56505117707799 * pi / 180.0;
		Real stheta = sin(theta);
		Real ctheta = cos(theta);
		ico.clear();

		ico.push_back(numeric::xyzVector<Real>(0.0,0.0,-1.0));
		Real phi = pi / 5.0;
		for (int i=0; i<5; ++i) {
		  ico.push_back(numeric::xyzVector<Real>( ctheta * cos(phi), ctheta * sin(phi), -stheta ));
		  phi += 2.0 * pi / 5.0;
		}
		phi = 0.0;
		for (int i=0; i<5; ++i) {
		  ico.push_back(numeric::xyzVector<Real>( ctheta * cos(phi), ctheta * sin(phi), stheta ));
		  phi += 2.0 * pi / 5.0;
		}
		ico.push_back(numeric::xyzVector<Real>(0.0,0.0,1.0));
		Size TRIS[][3] = {
			{0,2,1}, {0,3,2}, {0,4,3}, {0,5,4}, {0,1,5}, {1,2,7}, {2,3,8}, {3,4,9}, {4,5,10}, {5,1,6},
			{1,7,6}, {2,8,7}, {3,9,8}, {4,10,9}, {5,6,10}, {6,7,11}, {7,8,11}, {8,9,11}, {9,10,11}, {10,6,11}
		};
		Size EDGES[][2] =  {
			{0,1}, {0,2}, {0,3}, {0,4}, {0,5}, {1,2}, {1,5}, {1,6}, {1,7}, {2,3}, {2,7}, {2,8}, {3,4}, {3,8}, {3,9},
			{4,5}, {4,9}, {4,10}, {5,6}, {5,10}, {6,7}, {6,10}, {6,11}, {7,8}, {7,11}, {8,9}, {8,11}, {9,10}, {9,11}, {10,11}
		};

		// subdivide edges
		for (int i=0; i<30; ++i) {
			numeric::xyzVector<Real> a = ico[EDGES[i][0]+1];
			numeric::xyzVector<Real> b = ico[EDGES[i][1]+1];
			numeric::xyzVector<Real> ab = b-a;
			for (int j=1; j<(int)nsub; ++j) {
				numeric::xyzVector<Real> new_j;
				new_j = a+(((Real)j)/((Real)nsub))*ab;
				new_j = new_j/new_j.length();
		  		ico.push_back(new_j);
			}
		}

		// subdivide faces
		for (int i=0; i<20; ++i) {
			numeric::xyzVector<Real> a = ico[TRIS[i][0]+1];
			numeric::xyzVector<Real> b = ico[TRIS[i][1]+1];
			numeric::xyzVector<Real> c = ico[TRIS[i][2]+1];
			numeric::xyzVector<Real> ab = b-a;
			numeric::xyzVector<Real> ac = c-a;
			for (int j=1; j<(int)nsub-1; ++j)
			for (int k=j; k<(int)nsub-1; ++k) {
				numeric::xyzVector<Real> new_j;
				new_j = a + ((Real)j)/((Real)nsub)*ac + ((Real)(nsub-1-k))/((Real)nsub)*ab;
				new_j = new_j/new_j.length();
		  		ico.push_back(new_j);
			}
		}
	}

	UniformRotationSampler(Real theta) {
		theta_ = theta;
		if (theta<=1e-6) {
			rotlist_.resize(1);
			return;
		}

		using numeric::constants::d::pi;
		int S1subs = (int)std::floor(360.0 / theta + 0.5);
		int S2subs = (int)std::floor(60.0 / theta + 0.5);

		utility::vector1<numeric::xyzVector<Real> > ico;
		generateIcosahedralSamples(ico, S2subs);

		int nsamples = S1subs*ico.size();
		rotlist_.resize( nsamples );
		int counter=1;

		// precompute cos/sin (psi/2)for S1
		utility::vector1<core::Real> cos_psi_over_2(S1subs), sin_psi_over_2(S1subs);
		for (int j=0; j<S1subs; ++j) {
			Real psi_over_2 = j*pi/((Real)S1subs);
			cos_psi_over_2[j+1] = cos(psi_over_2);
			sin_psi_over_2[j+1] = sin(psi_over_2);
		}
		for (int i=1; i<=(int)ico.size(); ++i) {
			// convert to spherical ASSUMES NORMALIZED
			Real theta = pi-acos( ico[i][2] );
			Real phi = atan2( ico[i][1], ico[i][0] ) + pi;
			for (int j=0; j<S1subs; ++j) {
				Real psi = 2.0*j*pi/((Real)S1subs);
				Real w = sin(theta/2)*sin(phi+psi/2);
				Real x = cos(theta/2)*cos(psi/2);
				Real y = cos(theta/2)*sin(psi/2);
				Real z = sin(theta/2)*cos(phi+psi/2);
				rotlist_[counter++] = Quat(x,y,z,w);
			}
		}
	}

	void
	remove_redundant(
			utility::vector1<core::kinematics::RT> const &symmops,
			numeric::xyzMatrix<Real> const &i2c, numeric::xyzMatrix<Real> const &c2i ) {

		utility::vector1<bool> tokeep(rotlist_.size(), true);
		utility::vector1<numeric::xyzMatrix<Real> > Rs(rotlist_.size());
		utility::vector1<numeric::xyzMatrix<Real> > Rinvs(rotlist_.size());

		for (int i=1; i<=(int)rotlist_.size(); ++i) {
			Rs[i] = rotlist_[i].asR();
			Rinvs[i] = numeric::inverse( Rs[i] );
		}


		for (int s=2; s<=(int)symmops.size(); ++s) {
			numeric::xyzMatrix<Real> S = symmops[s].get_rotation();
			if (S.xy()==0 && S.xz()==0 && S.yz()==0 && S.yx()==0 && S.zx()==0 && S.zy()==0 && S.xx()==1 && S.yy()==1 && S.zz()==1)
				continue; // identity

			// cartesian-space S
			numeric::xyzMatrix<Real> Scart = i2c*S*c2i;

			for (int i=1; i<=(int)rotlist_.size(); ++i) {
				if (!tokeep[i]) { continue; }
				numeric::xyzMatrix<Real> SR = Scart*Rs[i];
				for (int j=i+1; j<=(int)rotlist_.size(); ++j) {
					if (!tokeep[j]) { continue; }

					Real ang_ij = R2ang( SR*Rinvs[j] );
					if (ang_ij<0.65*theta_) {
						tokeep[j]=false;
					}
				}
			}
		}

		utility::vector1< Quat > rotlist_old = rotlist_;
		rotlist_.clear();
		for (int i=1; i<=(int)rotlist_old.size(); ++i) {
			if (tokeep[i]) rotlist_.push_back(rotlist_old[i]);
		}

		TR << "Trimmed rot list has " << rotlist_.size() << " rotations (was " << rotlist_old.size() << ")\n";
	}
};


///////////////////////////////////////////////////////////////////////////////
// fpd this class stores information on the top-N interfaces
//    underlying heap makes this reasonably fast
//    N is 5 by default (covers up to 3 unique interfaces)
struct SingleInterface {
	SingleInterface( Size index_in, Real cb_overlap_in ) {
		index = index_in;
		cb_overlap = cb_overlap_in;
	}
    Size index;
	Real cb_overlap;
};

class SingleInterfaceComparitor {
public:
	bool operator()(SingleInterface& t1, SingleInterface& t2) {
		return (t1.cb_overlap > t2.cb_overlap);
	}
};

class InterfaceInfo {
private:
	core::Size N_;
    std::priority_queue<SingleInterface, std::vector<SingleInterface>, SingleInterfaceComparitor> queue_;

public:
	InterfaceInfo() { N_ = 5; }

	void add_interface( SingleInterface interface ) {
		if (queue_.size() <N_) {
			queue_.push( interface );
		} else if ( interface.cb_overlap > queue_.top().cb_overlap ) {
			queue_.pop();
			queue_.push( interface );
		}
	}

	SingleInterface pop() {
		SingleInterface retval = queue_.top();
		queue_.pop();
		return retval;
	}

	Size size() { return queue_.size(); }
};

///////////////////////////////////////////////////////////////////////////////

// fpd this class is used for storing hits over all rotations
//    same top-N priority_queue but much larger store
struct InterfaceHit {
	InterfaceHit( Real score_in, Real x_in, Real y_in, Real z_in, Size rot_index_in ) {
		score=score_in;
		x=x_in; y=y_in; z=z_in;
		rot_index=rot_index_in;
	}
	Real score;
	Real x,y,z;
	Size rot_index;
};

class InterfaceHitComparitor {
public:
	bool operator()(InterfaceHit& t1, InterfaceHit& t2) {
		return (t1.score > t2.score);
	}
};

class InterfaceHitDatabase {
private:
	core::Size N_;
    std::priority_queue<InterfaceHit, std::vector<InterfaceHit>, InterfaceHitComparitor> queue_;

public:
	InterfaceHitDatabase() { N_ = option[crystdock::nmodels]; }

	void add_interface( InterfaceHit interface ) {
		if (queue_.size() <N_) {
			queue_.push( interface );
		} else if ( interface.score > queue_.top().score ) {
			queue_.pop();
			queue_.push( interface );
		}
	}

	InterfaceHit pop() {
		InterfaceHit retval = queue_.top();
		queue_.pop();
		return retval;
	}

	InterfaceHit top() {
		InterfaceHit retval = queue_.top();
		return retval;
	}

	Size size() { return queue_.size(); }
};


///////////////////////////////////////////////////////////////////////////////

enum SpacegroupSetting {
	TRICLINIC,
	MONOCLINIC,
	ORTHORHOMBIC,
	TETRAGONAL,
	HEXAGONAL,  // includes trigonal
	CUBIC,
	UNDEFINED
};

struct CheshireCell {
	CheshireCell() {};
	CheshireCell( numeric::xyzVector<Real> low_in, numeric::xyzVector<Real> high_in ) {
		low=low_in;
		high=high_in;
	}
	numeric::xyzVector<Real> low, high;
};

class Spacegroup {
private:
	std::string name_;
	SpacegroupSetting setting_;

	// cryst stuff
	numeric::xyzMatrix<Real> f2c_, c2f_;
	Real a_, b_, c_, alpha_, beta_, gamma_, V_;
public:
	Spacegroup() {
		a_ = b_ = c_ = alpha_ = beta_ = gamma_ = V_ = 0;
		setting_ = UNDEFINED;
	}

	Spacegroup(std::string name_in) {
		init( name_in );
	}

	void init( std::string name_in ) {
 		name_ = name_in;

		// setting
		if ( name_ == "P1" )
			setting_ = TRICLINIC;
		else if ( (name_ == "P2") ||  (name_ == "P21") ||  (name_ == "C2") )
			setting_ = MONOCLINIC;
		else if ( (name_ == "P23") ||  (name_ == "F23") ||  (name_ == "I23") ||
			      (name_ == "P213") ||  (name_ == "I213") ||  (name_ == "P432") ||
			      (name_ == "P4232") ||  (name_ == "F432") ||  (name_ == "F4132") ||
			      (name_ == "I432") ||  (name_ == "P4332") ||  (name_ == "P4132") ||
			      (name_ == "I4132") )
			setting_ = CUBIC;
		else if ( (name_ == "P222") ||  (name_ == "P2221") ||  (name_ == "P21212") ||
			      (name_ == "P212121") ||  (name_ == "C2221") ||  (name_ == "C222") ||
			      (name_ == "F222") ||  (name_ == "I222") ||  (name_ == "I212121") )
			setting_ = ORTHORHOMBIC;
		else if ( (name_ == "P4") ||  (name_ == "P41") ||  (name_ == "P42") ||
			      (name_ == "P43") ||  (name_ == "I4") ||  (name_ == "I41") ||
			      (name_ == "P422") ||  (name_ == "P4212") ||  (name_ == "P4122") ||
			      (name_ == "P41212") ||  (name_ == "P4222") ||  (name_ == "P42212") ||
			      (name_ == "P4322") ||  (name_ == "P43212") ||  (name_ == "I422") ||
			      (name_ == "I4122") )
			setting_ = TETRAGONAL;
		else if ( (name_ == "P3") ||  (name_ == "P31") ||  (name_ == "P32") ||
			      (name_ == "H3") ||  (name_ == "P312") ||  (name_ == "P321") ||
			      (name_ == "P3112") ||  (name_ == "P3121") ||  (name_ == "P3212") ||
			      (name_ == "P3221") ||  (name_ == "H32") ||  (name_ == "P6") ||
			      (name_ == "P61") ||  (name_ == "P65") ||  (name_ == "P62") ||
			      (name_ == "P64") ||  (name_ == "P63") ||  (name_ == "P622") ||
			      (name_ == "P6122") ||  (name_ == "P6522") ||  (name_ == "P6222") ||
			      (name_ == "P6422") ||  (name_ == "P6322") )
			setting_ = HEXAGONAL;
		else
			utility_exit_with_message("Unknown spacegroup!");
	}

	numeric::xyzMatrix<Real> const &f2c() const { return f2c_; };
	numeric::xyzMatrix<Real> const &c2f() const { return c2f_; };
	Real A() const { return a_; }
	Real B() const { return b_; }
	Real C() const { return c_; }
	Real alpha() const { return alpha_; }
	Real beta() const { return beta_; }
	Real gamma() const { return gamma_; }
	Real volume() const { return V_; }

	// grid spacing must be a multiple of this number
	Size minmult() const {
		if (setting_ == TRICLINIC)   return 2;
		if (setting_ == MONOCLINIC)   return 4;
		if (setting_ == CUBIC)        return 4;
		if (setting_ == ORTHORHOMBIC) return 4;
		if (setting_ == TETRAGONAL)   return 8;
		if (setting_ == HEXAGONAL)    return 6;
		return 2;  // should not reach here ever
	}

	// sets and validates input parameters
	void set_parameters(Real a_in, Real b_in, Real c_in, Real alpha_in, Real beta_in, Real gamma_in) {
		if (setting_ == TRICLINIC) {
			a_=a_in; b_=b_in; c_=c_in; alpha_=alpha_in; beta_=beta_in; gamma_=gamma_in;
		}
		else if (setting_ == MONOCLINIC) {
			a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=beta_in; gamma_=90.0;
		}
		else if (setting_ == CUBIC) {
			a_=a_in; b_=a_in; c_=a_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
		}
		else if (setting_ == ORTHORHOMBIC) {
			a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
		}
		else if (setting_ == TETRAGONAL) {
			a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
		}
		else if (setting_ == HEXAGONAL) {
			a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=120.0;
		}

		// transformation matrices
		core::Real ca = cos(DEG2RAD*alpha_), cb = cos(DEG2RAD*beta_), cg = cos(DEG2RAD*gamma_);
		core::Real sb = sin(DEG2RAD*beta_), sg = sin(DEG2RAD*gamma_);
		core::Real f33 = (cb*cg - ca)/(sb*sg);
		f2c_ = numeric::xyzMatrix<core::Real>::rows(
			a_  , b_ * cg , c_ * cb,
			0.0 , b_ * sg , c_ * (ca - cb*cg) / sg,
			0.0 , 0.0     , c_ * sb * sqrt(1.0 - (f33*f33))
		);
		c2f_ = numeric::inverse(f2c_);
		V_ = a_*b_*c_* sqrt(1-(ca*ca)-(cb*cb)-(cg*cg)+2*ca*cb*cg);

		// report
		if (a_!=a_in || b_!=b_in || c_!=c_in || alpha_!=alpha_in || beta_!=beta_in || gamma_!=gamma_in) {
			// might want to die on this instead
			TR << "Overriding input crystal parameters with [ "
			          << a_ << "," << b_ << "," << c_ << " , " << alpha_ << ","  << beta_ << ","  << gamma_ << " ]" << std::endl;
		}
	}

	// get symmops
	void get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const;

	// get symmops
	std::string pdbname() const;

	// TO DO get this working for all settings
	Real get_interface_score(InterfaceInfo iinfo) {
		if (iinfo.size() < 3) return 0;

		iinfo.pop();
		iinfo.pop();
		Real r3 = iinfo.pop().cb_overlap;

		return 2*r3;
	}
};


///////////////////////////////////////////////////////////////////////////////

class CrystDock : public protocols::moves::Mover {
private:
	core::Real ca_clashdist_, cb_clashdist_, n_clashdist_, c_clashdist_, o_clashdist_;
	core::Real interfacedist_;
	bool ss_only_;

	numeric::xyzVector<int> grid_, oversamplegrid_;

	numeric::xyzMatrix<Real> i2c_, c2i_;

	Spacegroup sg_;

public:
	CrystDock() {
		// to do: make these options
		n_clashdist_  = 1.75;   // ros VDW=1.75
		ca_clashdist_ = 2.00;   // ros VDW=2.0
		c_clashdist_  = 2.00;   // ros VDW=2.0
		o_clashdist_  = 1.55;   // ros VDW=1.55
		cb_clashdist_ = 1.70;   // ros VDW=2.0    .. make this a bit smaller to encourage better packing
		interfacedist_ = 3.0;
		ss_only_ = true;
	}

	virtual std::string get_name() const { return "CrystDock"; }

	// write density grids in MRC --  debugging only for now
	void
	writeMRC(FArray3D<Real> density, std::string mapfilename, bool is_oversampled=false) {
		const int CCP4HDSIZE = 1024;  // size of CCP4/MRC header
		std::fstream outx( (mapfilename).c_str() , std::ios::binary | std::ios::out );

		float buff_f, buff_vf[3];
		int buff_i, buff_vi[3], symBytes = 0;

		if (!outx ) {
			TR << "Error opening " << mapfilename << " for writing." << std::endl;
			return;
		}

		// extent
		buff_vi[0] = density.u1(); buff_vi[1] = density.u2(); buff_vi[2] = density.u3();
		outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

		// mode
		buff_i = 2;
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

		// origin
		int ori_int[3] = {0,0,0};
		outx.write(reinterpret_cast <char*>(ori_int), sizeof(int)*3);

		// grid
		int grid[3] = {density.u1(),density.u2(),density.u3()};
		outx.write(reinterpret_cast <char*>(&grid[0]), sizeof(int)*3);

		// cell params
		Real Aeff=sg_.A(), Beff=sg_.B(), Ceff=sg_.C();
		if (is_oversampled) {
			Aeff *= ((Real)oversamplegrid_[0]) / ((Real)grid_[0]);
			Beff *= ((Real)oversamplegrid_[1]) / ((Real)grid_[1]);
			Ceff *= ((Real)oversamplegrid_[2]) / ((Real)grid_[2]);
		}

		float cellDimensions[3] = {Aeff,Beff,Ceff};
		float cellAngles[3] = {sg_.alpha(),sg_.beta(),sg_.gamma()};
		outx.write(reinterpret_cast <char*>(&cellDimensions), sizeof(float)*3);
		outx.write(reinterpret_cast <char*>(&cellAngles), sizeof(float)*3);

		// crs2xyz
		buff_vi[0] = 1; buff_vi[1] = 2; buff_vi[2] = 3;
		outx.write(reinterpret_cast <char*>(buff_vi), sizeof(int)*3);

		// min, max, mean dens
		buff_vf[0] = -100.0; buff_vf[1] = 100.0; buff_vf[2] = 0.0;
		outx.write(reinterpret_cast <char*>(buff_vf), sizeof(float)*3);

		// 4 bytes junk
		buff_i = 0;
		outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));

		// symmops
		outx.write(reinterpret_cast <char*>(&symBytes), sizeof(int));

		// 104 bytes junk
		buff_i = 0;
		for (int i=0; i<25; ++i) {
			outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
		}

		// alt origin (MRC)
		float ori_float[3]={0,0,0};
		outx.write(reinterpret_cast <char*>(ori_float), sizeof(float)*3);

		// Write "MAP" at byte 208, indicating a CCP4 file.
		char buff_s[80]; strcpy(buff_s, "MAP DD");
		outx.write(reinterpret_cast <char*>(buff_s), 8);

		// fill remainder of head with junk
		int nJunkWords = (CCP4HDSIZE - 216) /4;
		buff_i = 0;
		for (int i=0; i<nJunkWords; ++i) {
			outx.write(reinterpret_cast <char*>(&buff_i), sizeof(int));
		}

		// data
		int coord[3];
		for (coord[2] = 1; coord[2] <= density.u3(); coord[2]++) {
			for (coord[1] = 1; coord[1] <= density.u2(); coord[1]++) {
				for (coord[0] = 1; coord[0] <= density.u1(); coord[0]++) {
					buff_f = (float) density(coord[0],coord[1],coord[2]);
					outx.write(reinterpret_cast <char*>(&buff_f), sizeof(float));
				}
			}
		}
	}

	// build occupancy and interaction masks from pose
	// since we will be sampling rotations from this map, we sample over a larger volume than the unit cell
	void setup_maps( Pose & pose, FArray3D<Real> &rho_ca, FArray3D<Real> &rho_cb, Real trans_step) {
		Real ATOM_MASK_PADDING = 2.0;
		Real UNIT_CELL_PADDING = 4.0;  // IN GRID POINTS!

		Size minmult = sg_.minmult();

		// find true grid
		grid_[0] = minmult*(Size)std::floor(sg_.A()/(minmult*trans_step) + 0.5);
		grid_[1] = minmult*(Size)std::floor(sg_.B()/(minmult*trans_step) + 0.5);
		grid_[2] = minmult*(Size)std::floor(sg_.C()/(minmult*trans_step) + 0.5);

		numeric::xyzMatrix<Real> i2f, f2i;
		f2i = numeric::xyzMatrix<Real>::rows( (Real)grid_[0],0,0,  0,(Real)grid_[1],0,  0,0,(Real)grid_[2] );
		i2f = numeric::xyzMatrix<Real>::rows( 1.0/((Real)grid_[0]),0,0,  0,1.0/((Real)grid_[1]),0,  0,0,1.0/((Real)grid_[2]) );

		i2c_ = sg_.f2c()*i2f;
		c2i_ = f2i*sg_.c2f();

		// find oversampled grid
		oversamplegrid_ = numeric::xyzVector<int>(0,0,0);
		for (int i=1 ; i<=(int)pose.total_residue(); ++i) {
			if (!pose.residue(i).is_protein()) continue;
			for (int j=1; j<=4; ++j) {
				numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
				numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
				for (int k=0; k<3; ++k) oversamplegrid_[k] = std::max(oversamplegrid_[k], 2*(int)std::floor( (atm_idx[k]+UNIT_CELL_PADDING) )+1 );
			}
		}
		TR << "Base grid = [" <<  oversamplegrid_[0] << " , " <<  oversamplegrid_[1]<< " , " << oversamplegrid_[2] << "]" << std::endl;

		rho_ca.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_ca = 1;
		rho_cb.dimension( oversamplegrid_[0], oversamplegrid_[1], oversamplegrid_[2] ); rho_cb = 1;

		// loop over bb heavyatoms
		for (int i=1 ; i<=(int)pose.total_residue(); ++i) {
			if (!pose.residue(i).is_protein()) continue;

			for (int j=1; j<=4; ++j) {
				core::Real clashdist=0.0;
				if (j==1) clashdist = n_clashdist_;
				if (j==2) clashdist = ca_clashdist_;
				if (j==3) clashdist = c_clashdist_;
				if (j==4) clashdist = o_clashdist_;

				numeric::xyzVector< core::Real> xyz_j = pose.residue(i).atom(j).xyz();
				numeric::xyzVector< core::Real> atm_idx = c2i_*xyz_j;
				numeric::xyzVector< core::Real> atm_j, del_ij;

				for (int z=1; z<=oversamplegrid_[2]; ++z) {
					atm_j[2] = z-1;
					del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
					del_ij[0] = del_ij[1] = 0.0;
					if ((i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING)) continue;
					for (int y=1; y<=oversamplegrid_[1]; ++y) {
						atm_j[1] = y-1;
						del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
						del_ij[0] = 0.0;
						if ((i2c_*del_ij).length_squared() > (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING)) continue;
						for (int x=1; x<=oversamplegrid_[0]; ++x) {
							atm_j[0] = x-1;
							del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
							numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
							core::Real d2 = (cart_del_ij).length_squared();
							if (d2 <= (clashdist+ATOM_MASK_PADDING)*(clashdist+ATOM_MASK_PADDING)) {
								core::Real doff = sqrt(d2) - clashdist;
								core::Real sig = 1 / ( 1 + exp ( -6*doff ) );   // '6' is sigmoid dropoff ... (reso dependent?)
								rho_ca(x,y,z) *= sig;
							}
						}
					}
				}
			}
		}

		// loop over CBs
		for (int i=1 ; i<=(int)pose.total_residue(); ++i) {
			if (pose.residue(i).aa() == core::chemical::aa_gly) continue;
			numeric::xyzVector< core::Real> CB = pose.residue(i).atom(5).xyz();
			numeric::xyzVector< core::Real> atm_idx = c2i_*CB;
			numeric::xyzVector< core::Real> atm_j, del_ij;

			for (int z=1; z<=oversamplegrid_[2]; ++z) {
				atm_j[2] = z-1;
				del_ij[2] = min_mod(atm_idx[2] - atm_j[2], (Real)oversamplegrid_[2]);
				del_ij[0] = del_ij[1] = 0.0;
				if ((i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING)) continue;
				for (int y=1; y<=oversamplegrid_[1]; ++y) {
					atm_j[1] = y-1;
					del_ij[1] = min_mod(atm_idx[1] - atm_j[1], (Real)oversamplegrid_[1]);
					del_ij[0] = 0.0;
					if ((i2c_*del_ij).length_squared() > (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING)) continue;
					for (int x=1; x<=oversamplegrid_[0]; ++x) {
						atm_j[0] = x-1;
						del_ij[0] = min_mod(atm_idx[0] - atm_j[0], (Real)oversamplegrid_[0]);
						numeric::xyzVector< core::Real > cart_del_ij = i2c_*del_ij;
						core::Real d2 = (cart_del_ij).length_squared();
						if (d2 <= (interfacedist_+ATOM_MASK_PADDING)*(interfacedist_+ATOM_MASK_PADDING)) {
							core::Real doff = sqrt(d2) - cb_clashdist_;
							core::Real sig = 1 / ( 1 + exp ( -6*doff ) );   // '6' gives sigmoid dropoff
							rho_ca(x,y,z) *= sig;

							doff = sqrt(d2) - interfacedist_;
							sig = 1 / ( 1 + exp ( -6*doff ) );   // '6' gives sigmoid dropoff
							rho_cb(x,y,z) *= sig;
						}
					}
				}
			}
		}

		// factor CA mask out of CB mask; invert both
		for (int i=0 ; i<oversamplegrid_[0]*oversamplegrid_[1]*oversamplegrid_[2]; ++i) {
			rho_cb[i] = (1-rho_cb[i])*rho_ca[i];
			rho_ca[i] = (1-rho_ca[i]);
		}
	}

	// trilinear interpolation of map subject to real-space rotation
	void resample_maps(
			FArray3D<Real> const &rho_ca, FArray3D<Real> const &rho_cb,
			numeric::xyzMatrix<Real> R,
			FArray3D<Real> &r_rho_ca, FArray3D<Real> &r_rho_cb ) {
		r_rho_ca.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_ca=0;
		r_rho_cb.dimension( grid_[0], grid_[1], grid_[2] ); r_rho_cb=0;

		numeric::xyzMatrix<Real> Ri = numeric::inverse(R);

		// smarter logic here could possibly save some time
		for (int a=-1; a<=1; ++a)
		for (int b=-1; b<=1; ++b)
		for (int c=-1; c<=1; ++c) {
			for (int z=1; z<=grid_[2]; ++z)
			for (int y=1; y<=grid_[1]; ++y)
			for (int x=1; x<=grid_[0]; ++x) {
				int cx=x-1,cy=y-1,cz=z-1;
				if (cx>grid_[0]/2) cx -= grid_[0];
				if (cy>grid_[1]/2) cy -= grid_[1];
				if (cz>grid_[2]/2) cz -= grid_[2];

				cx += a*grid_[0];
				cy += b*grid_[1];
				cz += c*grid_[2];

				numeric::xyzMatrix<Real> Rgridspace = c2i_*Ri*i2c_;
				numeric::xyzVector<Real> rx = Rgridspace*numeric::xyzVector<Real>(cx,cy,cz);

				Real rho_ca_rx = interp_linear( rho_ca, rx );
				Real rho_cb_rx = interp_linear( rho_cb, rx );
				r_rho_ca(x,y,z) += rho_ca_rx;
				r_rho_cb(x,y,z) += rho_cb_rx;
			}

		}
	}

	void get_self_interactions(
			FArray3D<Real> const &r_rho_ca, FArray3D<Real> const &r_rho_cb,
			Real &ca_overlap, Real &cb_overlap ) {
		int nx=grid_[0], ny=grid_[1], nz=grid_[2];
		ca_overlap = cb_overlap = 0;

		// to compute these quantities exactly, we need the product of overlaps; we have the sum
		// we therefor assume:
		//   1 ) values <1 all come from 1 component (underestimation)
		//   2 ) values >1 are equally split between the two (overestimation)
		// add a buffer to handle interpolation effects
		// hopefully the over- and under- estimate cancel out
		for (int z=1; z<=nz; ++z)
		for (int y=1; y<=ny; ++y)
		for (int x=1; x<=nx; ++x) {
			if (r_rho_ca(x,y,z) > 1.000) {
				//if (r_rho_ca(x,y,z) > 1.0) TR << x<<","<<y<<","<<z<<": "<<r_rho_ca(x,y,z)<<std::endl;
				ca_overlap += (r_rho_ca(x,y,z)/2.0)*(r_rho_ca(x,y,z)/2.0);
			}
			if (r_rho_cb(x,y,z) > 1.000)
				cb_overlap += (r_rho_cb(x,y,z)/2.0)*(r_rho_cb(x,y,z)/2.0);
		}

		if (option[ crystdock::debug ]) {
			TR << "ca_overlap = " << ca_overlap << std::endl;
			TR << "cb_overlap = " << cb_overlap << std::endl;
		}
	}

	void center_pose_at_origin( Pose & pose ) {
		numeric::xyzVector<Real> com(0,0,0);
		int count=0;
		for (int i=1; i<=(int)pose.total_residue(); ++i) {
			if (!pose.residue(i).is_protein()) continue;
			com += pose.residue(i).atom(2).xyz();
			count++;
		}
		com /= count;

		for (int i=1; i<=(int)pose.total_residue(); ++i) {
			for (int j=1; j<=(int)pose.residue(i).natoms(); ++j) {
				pose.set_xyz(id::AtomID(j,i), pose.residue(i).atom(j).xyz()-com );
			}
		}
	}

	void add_crystinfo_to_pose( Pose & pose ) {
		CrystInfo ci;

		ci.A(sg_.A()); ci.B(sg_.B()); ci.C(sg_.C()); ci.alpha(sg_.alpha()); ci.beta(sg_.beta()); ci.gamma(sg_.gamma());
		ci.spacegroup(sg_.pdbname());

		pose.pdb_info()->set_crystinfo(ci);
	}

	void dump_transformed_pdb( Pose pose, InterfaceHit ih, UniformRotationSampler const &urs, Real score, std::string outname ) {
		numeric::xyzMatrix<Real> R;
		urs.get( ih.rot_index, R );
		numeric::xyzVector<Real> T (ih.x, ih.y, ih.z);

		// add score to header
		core::pose::RemarkInfo remark;
		std::ostringstream oss;
		oss << "  score = " << score << std::endl;
		remark.num = 1;	remark.value = oss.str();
		pose.pdb_info()->remarks().push_back( remark );

		pose.apply_transform_Rx_plus_v( R,T );
		pose.dump_pdb( outname );
	}

	// nearest-neighbor interpolation subject to grid-space transform
	void transform_map(
			FArray3D<Real> const &rho,
			numeric::xyzMatrix<Real> S, numeric::xyzVector<Real> T,
			FArray3D<Real> &Srho) {
		Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
		for (int z=1; z<=rho.u3(); ++z)
		for (int y=1; y<=rho.u2(); ++y)
		for (int x=1; x<=rho.u1(); ++x) {
			int cx=x-1,cy=y-1,cz=z-1;
			int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) + (T[0]*rho.u1()) + 0.5 ) , rho.u1());
			int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) + (T[1]*rho.u2()) + 0.5 ) , rho.u2());
			int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) + (T[2]*rho.u3()) + 0.5 ) , rho.u3());
			Real rho_sx = rho(rx,ry,rz);
			Srho(x,y,z) = rho_sx;
		}
	}

	// same as previous function, but no T and applies an offset of 1 (necessary for symm xform)
	void transform_map_offset1(
			FArray3D<Real> const &rho,
			numeric::xyzMatrix<Real> S,
			FArray3D<Real> &Srho) {
		Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
		for (int z=1; z<=rho.u3(); ++z)
		for (int y=1; y<=rho.u2(); ++y)
		for (int x=1; x<=rho.u1(); ++x) {
			int cx=x,cy=y,cz=z;
			int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) ) , rho.u1());
			int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) ) , rho.u2());
			int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) ) , rho.u3());
			Real rho_sx = rho(rx,ry,rz);
			Srho(x,y,z) = rho_sx;
		}
	}


	// same as previous function, applies an offset of 0
	void transform_map_offset0(
			FArray3D<Real> const &rho,
			numeric::xyzMatrix<Real> S,
			FArray3D<Real> &Srho) {
		Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
		for (int z=1; z<=rho.u3(); ++z)
		for (int y=1; y<=rho.u2(); ++y)
		for (int x=1; x<=rho.u1(); ++x) {
			int cx=x-1,cy=y-1,cz=z-1;
			int rx = 1 + pos_mod( (int)std::floor( (S.xx()*cx) + (S.xy()*cy) + (S.xz()*cz) ) , rho.u1());
			int ry = 1 + pos_mod( (int)std::floor( (S.yx()*cx) + (S.yy()*cy) + (S.yz()*cz) ) , rho.u2());
			int rz = 1 + pos_mod( (int)std::floor( (S.zx()*cx) + (S.zy()*cy) + (S.zz()*cz) ) , rho.u3());
			Real rho_sx = rho(rx,ry,rz);
			Srho(x,y,z) = rho_sx;
		}
	}

	// same as previous function, applies an offset of 0
	void transform_map_offset0_subsample(
			FArray3D<Real> const &rho,
			numeric::xyzMatrix<Real> S,
			FArray3D<Real> &Srho) {
		Srho.dimension( rho.u1(), rho.u2(), rho.u3() );
		for (int z=1; z<=rho.u3(); ++z)
		for (int y=1; y<=rho.u2(); ++y)
		for (int x=1; x<=rho.u1(); ++x) {
			int cx=x-1,cy=y-1,cz=z-1;
			numeric::xyzVector<Real> rx = S*numeric::xyzVector<Real>(cx,cy,cz);
			rx[0] = min_mod( rx[0], (Real)grid_[0]);
			rx[1] = min_mod( rx[1], (Real)grid_[1]);
			rx[2] = min_mod( rx[2], (Real)grid_[2]);
			Real rho_sx = interp_linear( rho, rx );
			Srho(x,y,z) = rho_sx;
		}
	}


	// wrapper does a series of 2d ffts
	//    this could be an option to the low-level fft code to possibly save time
	void
	fft2dslice( FArray3D<Real> const &rho, FArray3D< std::complex<Real> > &Frho, int axisSlice ) {
		FArray2D<Real> rhoSlice;
		FArray2D< std::complex<Real> > FrhoSlice;
		int xi = rho.u1(), yi = rho.u2(), zi = rho.u3();
		Frho.dimension(xi,yi,zi);

		if (axisSlice == 1) {
			rhoSlice.dimension(yi,zi);
			for (int ii=1; ii<=xi; ii++) {
				for (int jj=1; jj<=yi; jj++)
				for (int kk=1; kk<=zi; kk++) {
					rhoSlice(jj,kk) = rho(ii,jj,kk);
				}
				numeric::fourier::fft2(rhoSlice, FrhoSlice);
				for (int jj=1; jj<=yi; jj++)
				for (int kk=1; kk<=zi; kk++) {
					Frho(ii,jj,kk) = FrhoSlice(jj,kk);
				}
			}
		} else if (axisSlice == 2) {
			rhoSlice.dimension(xi,zi);
			for (int jj=1; jj<=yi; jj++) {
				for (int ii=1; ii<=xi; ii++)
				for (int kk=1; kk<=zi; kk++) {
					rhoSlice(ii,kk) = rho(ii,jj,kk);
				}
				numeric::fourier::fft2(rhoSlice, FrhoSlice);
				for (int ii=1; ii<=xi; ii++)
				for (int kk=1; kk<=zi; kk++) {
					Frho(ii,jj,kk) = FrhoSlice(ii,kk);
				}
			}
		} else if (axisSlice == 3) {
			rhoSlice.dimension(xi,yi);
			for (int kk=1; kk<=zi; kk++) {
				for (int ii=1; ii<=xi; ii++)
				for (int jj=1; jj<=yi; jj++) {
					rhoSlice(ii,jj) = rho(ii,jj,kk);
				}
				numeric::fourier::fft2(rhoSlice, FrhoSlice);
				for (int ii=1; ii<=xi; ii++)
				for (int jj=1; jj<=yi; jj++) {
					Frho(ii,jj,kk) = FrhoSlice(ii,jj);
				}
			}
		} else {
			utility_exit_with_message( "ERROR! Bad axis specified!");
		}
	}

	// wrapper does a series of 2d iffts
	void
	ifft2dslice( FArray3D< std::complex<Real> > const &Frho, FArray3D<Real> &rho, int axisSlice ) {
		FArray2D<Real> rhoSlice;
		FArray2D< std::complex<Real> > FrhoSlice;
		int xi = Frho.u1(), yi = Frho.u2(), zi = Frho.u3();
		rho.dimension(xi,yi,zi);
		rho = 0;

		if (axisSlice == 1) {
			FrhoSlice.dimension(yi,zi);
			for (int ii=1; ii<=xi; ii++) {
				for (int jj=1; jj<=yi; jj++)
				for (int kk=1; kk<=zi; kk++) {
					FrhoSlice(jj,kk) = Frho(ii,jj,kk);
				}
				numeric::fourier::ifft2(FrhoSlice, rhoSlice);
				for (int jj=1; jj<=yi; jj++)
				for (int kk=1; kk<=zi; kk++) {
					rho(ii,jj,kk) = rhoSlice(jj,kk);
				}
			}
		} else if (axisSlice == 2) {
			FrhoSlice.dimension(xi,zi);
			for (int jj=1; jj<=yi; jj++) {
				for (int ii=1; ii<=xi; ii++)
				for (int kk=1; kk<=zi; kk++) {
					FrhoSlice(ii,kk) = Frho(ii,jj,kk);
				}
				numeric::fourier::ifft2(FrhoSlice, rhoSlice);
				for (int ii=1; ii<=xi; ii++)
				for (int kk=1; kk<=zi; kk++) {
					rho(ii,jj,kk) = rhoSlice(ii,kk);
				}
			}
		} else if (axisSlice == 3) {
			FrhoSlice.dimension(xi,yi);
			for (int kk=1; kk<=zi; kk++) {
				for (int ii=1; ii<=xi; ii++)
				for (int jj=1; jj<=yi; jj++) {
					FrhoSlice(ii,jj) = Frho(ii,jj,kk);
				}
				numeric::fourier::ifft2(FrhoSlice, rhoSlice);
				for (int ii=1; ii<=xi; ii++)
				for (int jj=1; jj<=yi; jj++) {
					rho(ii,jj,kk) = rhoSlice(ii,jj);
				}
			}
		} else {
			utility_exit_with_message( "ERROR! Bad axis specified!");
		}
	}

	void
	project_along_axis( FArray3D<Real> &rho, numeric::xyzVector<Real> axis ) {
		FArray3D<Real> rho_input = rho;

		Real scalefact = std::max(std::fabs(axis[0]), std::max(std::fabs(axis[1]),std::fabs(axis[2])));
		axis /= scalefact;

		int ngrid = (std::fabs(axis[0]>0.999))? grid_[0] : ((std::fabs(axis[1])>0.999)? grid_[1] : grid_[2]);

		// for pp = 1:D-1; mfft = mfft+circshift( mfft_orig, [pp*slide(1,1),pp*sl ide(2,2),pp*slide(3,3)]);
 		for (int P=1; P<ngrid; ++P) {  // one less than full cycle
 			for (int z=1; z<=grid_[2]; ++z)
 			for (int y=1; y<=grid_[1]; ++y)
 			for (int x=1; x<=grid_[0]; ++x) {
 				int xt = pos_mod((int)floor( x+P*axis[0]-0.5 ), grid_[0]) + 1;
 				int yt = pos_mod((int)floor( y+P*axis[1]-0.5 ), grid_[1]) + 1;
 				int zt = pos_mod((int)floor( z+P*axis[2]-0.5 ), grid_[2]) + 1;
 				rho(x,y,z) += rho_input(xt,yt,zt);
 			}
 		}
	}

	// the main convolution operation used for both CA and CB maps
	void
	do_convolution( FArray3D<Real> const &rho, FArray3D<Real> const &Srho, numeric::xyzMatrix<Real> S, FArray3D<Real> &conv_out) {
		FArray3D<Real> working_s, working_trans, working_strans;
		FArray3D< std::complex<Real> > Fworking_trans, Fworking_strans;
		Size Npoints = rho.u1()*rho.u2()*rho.u3();

		// find out what plane we are slicing along
		Quat Sinv_Q(S);

		// many different possibilities explain the convoluted logic below
		//    - this may return (1,0,0), (sqrt(1/2),sqrt(1/2),0) or (1/2,1/2,1/2) for X, XY, and XYZ symm axes
		//    - also in some spacegps( e.g. P3112[5] ) it may be (1,1/4,0) -- since it's not a rotation but a skew
		//    - in others ( e.g. P6122[8] ) it may also be (1,1/4,0) but we need to pre-skew
		//    - in these cases, we want to process as 'x' rotation with the skew handled in the projection
		int slice_x = (Sinv_Q.x_>=0.5) ? 1 : ((Sinv_Q.x_<=-0.5) ? -1:0);
		int slice_y = (Sinv_Q.y_>=0.5) ? 1 : ((Sinv_Q.y_<=-0.5) ? -1:0);
		int slice_z = (Sinv_Q.z_>=0.5) ? 1 : ((Sinv_Q.z_<=-0.5) ? -1:0);

		numeric::xyzVector<Real> rotaxis(
			std::sqrt(std::fabs(Sinv_Q.x_)),
			std::sqrt(std::fabs(Sinv_Q.y_)),
			std::sqrt(std::fabs(Sinv_Q.z_)));
		if (Sinv_Q.x_<0) rotaxis[0] *= -1;
		if (Sinv_Q.y_<0) rotaxis[1] *= -1;
		if (Sinv_Q.z_<0) rotaxis[2] *= -1;

		// do we need to do fft convolution at all?
		if (slice_x==0 && slice_y==0 && slice_z==0) {

			// translation only, just take dot product
			if (option[ crystdock::debug ]|| option[ crystdock::debug_exact ]) {
				TR << "no slice"  << std::endl;
				TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
			}
			core::Real dotProd=0;
			for (int i=0; i<(int)Npoints; ++i) dotProd += rho[i]*Srho[i];
			conv_out.dimension( rho.u1(),rho.u2(),rho.u3() );
			conv_out = dotProd;
		} else {

			// do the fft convolution

			// Rinv maps xyz->pij (our reindexed coordinates with p is perpendicular to the symmaxis)
			// R maps pij->xyz
			numeric::xyzMatrix<Real> Rinv;
			int slice_axis = 0;
			bool R_is_identity=false;
			if (slice_x != 0) {
				slice_axis = 1;
				Rinv.row_x( Vector( slice_x, slice_y, slice_z ) );
				R_is_identity = (slice_y==0 && slice_z==0);
			} else {
				Rinv.row_x( Vector( 1, 0, 0 ) );
			}

			if (slice_x == 0 && slice_y != 0) {
				slice_axis = 2;
				Rinv.row_y( Vector( slice_x, slice_y, slice_z ) );
				R_is_identity = (slice_z==0);
			} else {
				Rinv.row_y( Vector( 0, 1, 0 ) );
			}
			if (slice_x == 0 && slice_y == 0 && slice_z != 0) {
				slice_axis = 3;
				Rinv.row_z( Vector( slice_x, slice_y, slice_z ) );
				R_is_identity = true;
			} else {
				Rinv.row_z( Vector( 0, 0, 1 ) );
			}

			// oddball case (all rotations of it too)
			if (S == numeric::xyzMatrix<Real>::rows(1,-1,0,   0,-1,0,   0,0,-1) ) {
				Rinv = numeric::xyzMatrix<Real>::rows( 2,-1,0,  1,0,0,  0,0,1 );
				R_is_identity = false;
				rotaxis = Vector(1,0,0);
			}
			if (S == numeric::xyzMatrix<Real>::rows(-1,0,0, -1,1,0, 0,0,-1) ) {
				Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,2,0, 0,0,1 );
				R_is_identity = false;
				rotaxis = Vector(0,1,0);
			}
			if (S == numeric::xyzMatrix<Real>::rows(-1,0,0, 0,-1,0, 0,-1,1) ) {
				Rinv = numeric::xyzMatrix<Real>::rows(1,0,0, 0,0,1, 0,-1,2);
				R_is_identity = false;
				rotaxis = Vector(0,0,1);
			}
			if (S == numeric::xyzMatrix<Real>::rows(1,0,-1, 0,-1,0, 0,0,-1) ) {
				Rinv = numeric::xyzMatrix<Real>::rows( 2,0,-1, 0,1,0, 1,0,0 );
				R_is_identity = false;
				rotaxis = Vector(1,0,0);
			}
			if (S == numeric::xyzMatrix<Real>::rows(0,0,-1, -1,0,0, 0,1,-1) ) {
				Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, 1,0,0, 0,2,-1 );
				R_is_identity = false;
				rotaxis = Vector(0,0,1);
			}
			if (S == numeric::xyzMatrix<Real>::rows(0,-1,0, -1,0,1, -1,0,0) ) {
				Rinv = numeric::xyzMatrix<Real>::rows( 0,1,0, -1,0,2, 0,0,1 );
				R_is_identity = false;
				rotaxis = Vector(0,1,0);
			}

			numeric::xyzMatrix<Real> R = numeric::inverse(Rinv);

			// Q is our symmetric transformation in our new coordinates
			numeric::xyzMatrix<Real> Q = Rinv*S*R;

			// Qstar is the FFT transformation
			//   Q-I on diagonal elements not on the slice axis
			//   partially apply the elements in the column of the slice axis (only applied to symm-sampled map)
			//        to compensate for point group sliding as we move up the symm axis in resampled groups
			numeric::xyzMatrix<Real> QstarF, QstarG;
			if (slice_axis == 1) {
				QstarF = numeric::xyzMatrix<Real>::rows( 1,0,0,   Q.yx(),Q.yy()-1,Q.yz(),   Q.zx(),Q.zy(),Q.zz()-1 );
				QstarG = numeric::xyzMatrix<Real>::rows( 1,0,0,        0,Q.yy()-1,Q.yz(),        0,Q.zy(),Q.zz()-1 );
			} else if (slice_axis == 2) {
				QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,  Q.xy(),Q.xz(), 0,1,0, Q.zx(),  Q.zy(),Q.zz()-1 );
				QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,       0,Q.xz(), 0,1,0, Q.zx(),       0,Q.zz()-1 );
			} else if (slice_axis == 3) {
				QstarF = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),  Q.xz(), Q.yx(),Q.yy()-1,  Q.yz(), 0,0,1 );
				QstarG = numeric::xyzMatrix<Real>::rows( Q.xx()-1,Q.xy(),       0, Q.yx(),Q.yy()-1,       0, 0,0,1 );
			}

			if (option[ crystdock::debug ]|| option[ crystdock::debug_exact ]) {
				TR << "  slice on axis :    " << slice_axis << std::endl;
				TR << "slice along axis: " << slice_x << " " << slice_y << " " << slice_z << std::endl;
				TR << "                : " << Sinv_Q.x_ << " " << Sinv_Q.y_ << " " << Sinv_Q.z_ << std::endl;
				TR << "                : " << rotaxis[0] << " " << rotaxis[1] << " " << rotaxis[2] << std::endl;
				TR << "         S = [" << S.xx() << "," <<  S.xy() << "," <<  S.xz() << ";  "
				                       << S.yx() << "," <<  S.yy() << "," <<  S.yz() << ";  "
				                       << S.zx() << "," <<  S.zy() << "," <<  S.zz() << "]" << std::endl;
				if (!R_is_identity) {
					TR << "      Rinv = [" << Rinv.xx() << "," << Rinv.xy() << "," << Rinv.xz() << ";  "
										   << Rinv.yx() << "," << Rinv.yy() << "," << Rinv.yz() << ";  "
										   << Rinv.zx() << "," << Rinv.zy() << "," << Rinv.zz() << "]" << std::endl;
					TR << "         R = [" << R.xx() << "," <<  R.xy() << "," <<  R.xz() << ";  "
										   << R.yx() << "," <<  R.yy() << "," <<  R.yz() << ";  "
										   << R.zx() << "," <<  R.zy() << "," <<  R.zz() << "]" << std::endl;
				}
				TR << "         Q = [" << Q.xx() << "," << Q.xy() << "," << Q.xz() << ";  "
									   << Q.yx() << "," << Q.yy() << "," << Q.yz() << ";  "
									   << Q.zx() << "," << Q.zy() << "," << Q.zz() << "]" << std::endl;
				TR << "        Qf = [" << QstarF.xx() << "," <<  QstarF.xy() << "," <<  QstarF.xz() << ";  "
									   << QstarF.yx() << "," <<  QstarF.yy() << "," <<  QstarF.yz() << ";  "
									   << QstarF.zx() << "," <<  QstarF.zy() << "," <<  QstarF.zz() << "]" << std::endl;
				TR << "        Qg = [" << QstarG.xx() << "," <<  QstarG.xy() << "," <<  QstarG.xz() << ";  "
									   << QstarG.yx() << "," <<  QstarG.yy() << "," <<  QstarG.yz() << ";  "
									   << QstarG.zx() << "," <<  QstarG.zy() << "," <<  QstarG.zz() << "]" << std::endl;
			}

			// transform
			transform_map_offset0( rho, R, working_trans);
			transform_map_offset0( Srho, R, working_strans);

			// transform in each plane
			// fft in each plane
			FArray3D<Real> rworking_trans, rworking_strans;
			transform_map_offset0( working_trans, QstarF, rworking_trans);
			transform_map_offset0( working_strans, QstarG, rworking_strans);

			fft2dslice( rworking_trans, Fworking_trans, slice_axis );
			fft2dslice( rworking_strans, Fworking_strans, slice_axis );
			for (int i=0; i<(int)Npoints; ++i) Fworking_trans[i] *= std::conj(Fworking_strans[i]);


			if (R_is_identity) {
				ifft2dslice(Fworking_trans, conv_out, slice_axis);
			} else {
				ifft2dslice(Fworking_trans, working_trans, slice_axis);
				transform_map_offset0( working_trans, Rinv, conv_out);
			}
			project_along_axis( conv_out, rotaxis );

			//project_along_axis( conv_out, numeric::xyzVector<int>( slice_x, slice_y, slice_z ) );
		}
	}

	void apply( Pose & pose) {
		// set up crystal info
		sg_.init( option[ crystdock::spacegroup ] );
		sg_.set_parameters(
			option[ crystdock::A ],option[ crystdock::B ],option[ crystdock::C ],
			option[ crystdock::alpha ],option[ crystdock::beta ],option[ crystdock::gamma ] );

		Real rotstep = option[ crystdock::rot_step ];
		Real maxclash = option[ crystdock::maxclash ];

		// center pose at origin
		center_pose_at_origin( pose );
		add_crystinfo_to_pose( pose );

		// get SS
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );

		// lookup symmops
		utility::vector1<core::kinematics::RT> rts;
		CheshireCell cc;
		sg_.get_symmops( rts, cc );

		if (option[ crystdock::debug ] || option[ crystdock::debug_exact ] )
			TR << "search [ " << cc.low[0] <<","<< cc.low[1] <<","<< cc.low[2] << " : "
			          << cc.high[0] <<","<< cc.high[1] <<","<< cc.high[2] << "]" << std::endl;

		// compute rho_ca, rho_cb
		FArray3D<Real> rho_ca, rho_cb, r_rho_ca, r_rho_cb;

		setup_maps( pose, rho_ca, rho_cb, option[ crystdock::trans_step ]);

		Size Npoints = grid_[0]*grid_[1]*grid_[2];
		Real voxel_volume = sg_.volume() / Npoints;

		numeric::xyzVector<Size> ccIndexLow(
			(Size)std::floor(cc.low[0]*grid_[0]+1.5),
			(Size)std::floor(cc.low[1]*grid_[1]+1.5),
			(Size)std::floor(cc.low[2]*grid_[2]+1.5));
		numeric::xyzVector<Size> ccIndexHigh(
			(Size)std::floor(cc.high[0]*grid_[0]+0.5),
			(Size)std::floor(cc.high[1]*grid_[1]+0.5),
			(Size)std::floor(cc.high[2]*grid_[2]+0.5));

		for (int i=0; i<3; i++)
			if (ccIndexHigh[i]<ccIndexLow[i])
				ccIndexHigh[i]=ccIndexLow[i];

		if (option[ crystdock::debug ]|| option[ crystdock::debug_exact ]) {
			writeMRC( rho_ca, "ca_mask.mrc", true );
			writeMRC( rho_cb, "cb_mask.mrc", true );
		}

		// space for intermediate results
		FArray3D<Real> working_s, conv_out;

		numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
		numeric::xyzMatrix<Real> r_local;

		// the collection of hits
		InterfaceHitDatabase IDB;
		core::Size nnonclashing = 0;

		// foreach rotation
		UniformRotationSampler urs( rotstep );
		urs.remove_redundant( rts, i2c_, c2i_ );

		Size rot_lb = 1, rot_ub = urs.nrots();
		if ( option[ crystdock::rotnum ].user() ) {
			rot_lb = rot_ub = option[ crystdock::rotnum ]();
		}
		Size sym_lb = 2, sym_ub = rts.size();
		if ( option[ crystdock::symnum ].user() ) {
			sym_lb = sym_ub = option[ crystdock::symnum ]();
		}

		for (int ctr=(int)rot_lb; ctr<=(int)rot_ub ; ++ctr) {
			TR << "Rotation " << ctr << " of " << urs.nrots() << std::endl;
			urs.get(ctr, r_local);

			resample_maps( rho_ca, rho_cb, r_local, r_rho_ca, r_rho_cb );

			if (option[ crystdock::debug ]|| option[ crystdock::debug_exact ] ) {
				std::ostringstream oss1; oss1 << "rot"<<ctr<<".mrc";
				writeMRC( r_rho_ca, oss1.str() );
				std::ostringstream oss2; oss2 << "rot"<<ctr<<".pdb";
				dump_transformed_pdb( pose, InterfaceHit( 0,0,0,0, ctr ), urs, 0.0, oss2.str() );
			}

			// store info
			FArray3D<Real> collision_map, ex_collision_map;
			collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
			collision_map=0;

			if (option[ crystdock::debug_exact ] ) {
				ex_collision_map.dimension( grid_[0] , grid_[1] , grid_[2] );
				ex_collision_map=0;
			}

			FArray3D<InterfaceInfo> interface_map;
			interface_map.dimension( grid_[0] , grid_[1] , grid_[2] );

			// get "P1" interactions
			Real self_ca, self_cb;
			get_self_interactions( r_rho_ca, r_rho_cb , self_ca, self_cb );

			// for exact debugging always run
			if (!option[ crystdock::debug_exact ] ) {
				if (self_ca >= maxclash) {
					TR << "   self clashing!" << std::endl;
					continue;   // no need to continue
				}
			}

			// foreach symmop
			for (int s=sym_lb; s<=(int)sym_ub; ++s) {   // s==1 is always identity
				numeric::xyzMatrix<Real> s_i = rts[s].get_rotation();
				numeric::xyzMatrix<Real> s_inv = numeric::inverse(s_i);
				numeric::xyzVector<Real> t_inv = s_inv*(-rts[s].get_translation());

				//transform_map( r_rho_ca, rts[s].get_rotation(), rts[s].get_translation(), working_s);
				transform_map( r_rho_ca, s_inv, t_inv, working_s);

				// all the magic is in here
				do_convolution( r_rho_ca, working_s, s_i, conv_out);

				for (int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z)
				for (int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y)
				for (int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x) {
					 collision_map(x,y,z) += conv_out(x,y,z);
				}

				// debug: exact CA fft
				if (option[ crystdock::debug_exact ] ) {
					FArray3D<Real> shift_rho_ca, s_shift_rho_ca;
					shift_rho_ca.dimension( grid_[0] , grid_[1] , grid_[2] );
					for (int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz)
					for (int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy)
					for (int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx) {
						for (int z=1; z<=grid_[2]; ++z)
						for (int y=1; y<=grid_[1]; ++y)
						for (int x=1; x<=grid_[0]; ++x) {
							int cx = pos_mod( x-sx-1 , grid_[0] ) + 1;
							int cy = pos_mod( y-sy-1 , grid_[1] ) + 1;
							int cz = pos_mod( z-sz-1 , grid_[2] ) + 1;
							shift_rho_ca(x,y,z) = r_rho_ca(cx,cy,cz);
						}
						transform_map( shift_rho_ca, rts[s].get_rotation(), rts[s].get_translation(), s_shift_rho_ca);
						for (int i=0; i<(int)Npoints; ++i) {
							ex_collision_map(sx+1,sy+1,sz+1) += shift_rho_ca[i] * s_shift_rho_ca[i];
						}
					}
				}


				// do CB fft
				transform_map( r_rho_cb, s_inv, t_inv, working_s);
				do_convolution( r_rho_cb, working_s, s_i, conv_out);
				for (int i=0; i<(int)Npoints; ++i) interface_map[i].add_interface( SingleInterface(s, conv_out[i]) );
			}

			if (option[ crystdock::debug ] || option[ crystdock::debug_exact ]) {
				std::ostringstream oss; oss << "collisionmap_"<<ctr<<".mrc";
				FArray3D<Real> collision_map_dump = collision_map;
				for (int i=0; i<(int)Npoints; ++i) collision_map_dump[i] = /*maxclash-*/collision_map[i];
				writeMRC( collision_map_dump, oss.str() );

				if( option[ crystdock::debug_exact ] ) {
					std::ostringstream oss2; oss2 << "ex_collisionmap_"<<ctr<<".mrc";
					for (int i=0; i<(int)Npoints; ++i) collision_map_dump[i] = /*maxclash-*/ex_collision_map[i];
					writeMRC( collision_map_dump, oss2.str() );

					// get correl
					Real x2=0,y2=0,xy=0,x=0,y=0, Ncc=0;
					for (int sz=(int)ccIndexLow[2]-1; sz<=(int)ccIndexHigh[2]-1; ++sz)
					for (int sy=(int)ccIndexLow[1]-1; sy<=(int)ccIndexHigh[1]-1; ++sy)
					for (int sx=(int)ccIndexLow[0]-1; sx<=(int)ccIndexHigh[0]-1; ++sx) {
						x2 += collision_map(sx+1,sy+1,sz+1) * collision_map(sx+1,sy+1,sz+1);
						y2 += ex_collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
						xy += collision_map(sx+1,sy+1,sz+1) * ex_collision_map(sx+1,sy+1,sz+1);
						x  += collision_map(sx+1,sy+1,sz+1);
						y  += ex_collision_map(sx+1,sy+1,sz+1);
						Ncc++;
					}
					Real correl = 1, slope = 0;
					Real sx = (Ncc*x2-x*x);
					Real sy = (Ncc*y2-y*y);
					if (sx*sy > 0) {
						correl = (Ncc*xy - x*y) / std::sqrt( (sx) * (sy) );
						slope = correl * sx/sy;
					}
					TR << "correl = " << correl << "   scale = " << slope << std::endl;

					collision_map = ex_collision_map;
				}
			}

			// now add nonclashing interfaces to the DB
			for (int z=(int)ccIndexLow[2]; z<=(int)ccIndexHigh[2]; ++z)
			for (int y=(int)ccIndexLow[1]; y<=(int)ccIndexHigh[1]; ++y)
			for (int x=(int)ccIndexLow[0]; x<=(int)ccIndexHigh[0]; ++x) {
				if (collision_map(x,y,z) < maxclash) {
					nnonclashing++;
					core::Real score_xyz = voxel_volume * sg_.get_interface_score( interface_map(x,y,z) );
					numeric::xyzVector<Real> xyz((Real)x-1,(Real)y-1,(Real)z-1);
					xyz = i2c_*xyz;
					IDB.add_interface( InterfaceHit( score_xyz, xyz[0],xyz[1],xyz[2], ctr ) );
				}
			}
			if (IDB.size()>0)
				TR << IDB.size() << " of " << nnonclashing << " nonclashing configurations; min_score = " << IDB.top().score << std::endl;
			else
				TR << IDB.size() << " configurations" << std::endl;
		}

		// DONE!  dump hits to stdout
		int nhits = IDB.size();
		for (int i=1; i<=nhits; ++i) {
			InterfaceHit ih = IDB.pop();
			TR << i << ": " << ih.score << " " << ih.x << " "  << ih.y << " "  << ih.z << " " << ih.rot_index  << " " << std::endl;

			// Treat tags as file names so that we put the number before the extension.
			std::string base_name = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
			utility::vector1< std::string > temp_out_names= utility::split( base_name );
			utility::file::FileName out_name = utility::file::combine_names( temp_out_names );
			base_name = out_name.base();
			std::string outname = base_name+"_"+option[ out::suffix ]()+right_string_of( i, 8, '0' )+".pdb";
			dump_transformed_pdb( pose, ih, urs, ih.score, outname );
		}
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( new CrystDock() );

	// force some options
	option[ out::nooutput ].value(true);
	option[ in::preserve_crystinfo ].value(true);

	// main loop
	protocols::jd2::JobDistributor::get_instance()->go( seq );

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
try {
	NEW_OPT(crystdock::spacegroup, "spacegroup (no spaces)", "P23");
	NEW_OPT(crystdock::A, "unit cell A", 50);
	NEW_OPT(crystdock::B, "unit cell B", 50);
	NEW_OPT(crystdock::C, "unit cell C", 50);
	NEW_OPT(crystdock::alpha, "unit cell alpha", 90);
	NEW_OPT(crystdock::beta, "unit cell beta", 90);
	NEW_OPT(crystdock::gamma, "unit cell gamma", 90);
	NEW_OPT(crystdock::maxclash, "max allowed clashscore", 5);
	NEW_OPT(crystdock::trans_step, "translational stepsize (A)", 1);
	NEW_OPT(crystdock::rot_step, "rotational stepsize (degrees) ([debug] 0 searches input rotation only)", 10);
	NEW_OPT(crystdock::nmodels, "number of models to output", 1000);
	NEW_OPT(crystdock::rotnum, "[debug] only run a single rotation", 0);
	NEW_OPT(crystdock::symnum, "[debug] only run a single symmop", 0);
	NEW_OPT(crystdock::debug, "[debug] dump intermediate info", false);
	NEW_OPT(crystdock::debug_exact, "[debug] debug mode with exact (non-FFT) calculations (slow!)", false);

	devel::init( argc, argv );
	protocols::viewer::viewer_main( my_main );

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}
	return 0;
}


///
///
///
void Spacegroup::get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const {
	if ( name_ == "P1" ) {
		rt_out.resize(1);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(0 ,0 ,0 ) );
	}
	if ( name_ == "P121" || name_ == "P2") {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>(0,0,0),numeric::xyzVector<Real>(0.5,0,0.5) );
	}
	else if ( name_ == "P1211" || name_ == "P21" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0.5,0) );
		//rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0,0) );
		//rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(0.5 ,0 ,0.5 ) );
	}
	else if ( name_ == "C121" || name_ == "C2" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(0.5 ,0 ,0.5 ) );
	}
	else if ( name_ == "P4" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "P41" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "P42" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "P43" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "I4" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "I41" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
 		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0 ) );
	}
	else if ( name_ == "P422" ) {
 		rt_out.resize(8);
 		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
 		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P4212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P4122" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P41212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P4222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P42212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P4322" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.75) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P43212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "I422" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=8; ++ii) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "I4122" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<Real>(0.5,0,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=8; ++ii) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,0.5 ,0.5 ) );
	}
	else if ( name_ == "P23" ) {
		rt_out.resize(12);
  		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
  		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>(0,0,0),numeric::xyzVector<Real>(1,1,1) );
	}
	else if ( name_ == "F23" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[24+ii] = rt_out[ii];
			rt_out[36+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(0.5,0.5,0.5 ) );
	}
	else if ( name_ == "I23" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "P213" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "I213" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		for (int ii=1; ii<=12; ++ii) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "P432" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "P4232" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "F432" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(0.5,0.5,0.5 ) );
	}
	else if ( name_ == "F4132" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		// cenops
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(0.5,0.5,0.5 ) );
	}
	else if ( name_ == "I432" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "P4332" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "P4132" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.75) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if ( name_ == "I4132" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<Real>(0.5,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<Real>(0.75,0.25,0.75) );
		// cenops
		for (int ii=1; ii<=24; ++ii) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0,0,0),numeric::xyzVector<Real>(1,1,1 ) );
	}
	else if (name_ == "P3") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	}
	else if (name_ == "P31") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	}
	else if (name_ == "P32") {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	}
	else if (name_ == "H3") {
		rt_out.resize(9);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=3; ++ii) {
			rt_out[3+ii] = rt_out[ii];
			rt_out[3+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	}
	else if (name_ == "P312") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	}
	else if (name_ == "P321") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P3112") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	}
	else if (name_ == "P3121") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P3212") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	}
	else if (name_ == "P3221") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "H32") {
		rt_out.resize(18);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=6; ++ii) {
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P61") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P65") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P62") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P64") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P63") {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0 ) );
	}
	else if (name_ == "P622") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6122") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6522") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.166666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.833333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6222") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6422") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P6322") {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>( 0, 0, 0),numeric::xyzVector<Real>(1 ,1 ,0.5 ) );
	}
	else if (name_ == "P222") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "P2221") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "P21212") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "P212121") {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "C2221") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "C222") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "F222") {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[8+ii] = rt_out[ii];
			rt_out[12+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0,0.5,0.5) );
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0,0.5) );
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "I222") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
	else if (name_ == "I212121") {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<Real>(0,0.5,0.5) );
		// cenops
		for (int ii=1; ii<=4; ++ii) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<Real>(0, 0, 0),numeric::xyzVector<Real>(0.5 ,0.5 ,0.5) );
	}
}

// name output in pdbheader
std::string Spacegroup::pdbname() const {
	if ( name_ == "P1" ) return "P 1";
	if ( name_ == "P121" || name_ == "P2") return "P 1 2 1";
	if ( name_ == "P1211" || name_ == "P21" ) return "P 1 21 1";
	if ( name_ == "C121" || name_ == "C2" ) return "C 1 2 1";
	if ( name_ == "P4" ) return "P 4";
	if ( name_ == "P41" ) return "P 41";
	if ( name_ == "P42" ) return "P 42";
	if ( name_ == "P43" ) return "P 43";
	if ( name_ == "I4" ) return "I 4";
	if ( name_ == "I41" )return "I 41";
	if ( name_ == "P422" ) return "P 4 2 2";
	if ( name_ == "P4212" ) return "P 4 21 2";
	if ( name_ == "P4122" ) return "P 41 2 2";
	if ( name_ == "P41212" ) return "P 41 21 2";
	if ( name_ == "P4222" ) return "P 42 2 2";
	if ( name_ == "P42212" ) return "P 42 21 2";
	if ( name_ == "P4322" ) return "P 43 2 2";
	if ( name_ == "P43212" ) return "P 43 21 2";
	if ( name_ == "I422" ) return "I 4 2 2";
	if ( name_ == "I4122" ) return "I 41 2 2";
	if ( name_ == "P23" ) return "P 2 3";
	if ( name_ == "F23" ) return "F 2 3";
	if ( name_ == "I23" ) return "I 2 3";
	if ( name_ == "P213" ) return "P 21 3";
	if ( name_ == "I213" ) return "I 21 3";
	if ( name_ == "P432" ) return "P 4 3 2";
	if ( name_ == "P4232" ) return "P 42 3 2";
	if ( name_ == "F432" ) return "F 4 3 2";
	if ( name_ == "F4132" ) return "F 41 3 2";
	if ( name_ == "I432" ) return "I 4 3 2";
	if ( name_ == "P4332" ) return "P 43 3 2";
	if ( name_ == "P4132" ) return "P 41 3 2";
	if ( name_ == "I4132" ) return "I 41 3 2";
	if ( name_ == "P3") return "P 3";
	if ( name_ == "P31") return "P 31";
	if ( name_ == "P32") return "P 32";
	if ( name_ == "H3") return "H 3";
	if ( name_ == "P312") return "P 3 1 2";
	if ( name_ == "P321") return "P 3 2 1";
	if ( name_ == "P3112") return "P 31 1 2";
	if ( name_ == "P3121") return "P 31 2 1";
	if ( name_ == "P3212") return "P 3 21 2";
	if ( name_ == "P3221") return "P 3 2 21";
	if ( name_ == "H32") return "H 3 2";
	if ( name_ == "P6") return "P 6";
	if ( name_ == "P61") return "P 61";
	if ( name_ == "P65") return "P 65";
	if ( name_ == "P62") return "P 62";
	if ( name_ == "P64") return "P 64";
	if ( name_ == "P63") return "P 63";
	if ( name_ == "P622") return "P 6 2 2";
	if ( name_ == "P6122") return "P 61 2 2";
	if ( name_ == "P6522") return "P 65 2 2";
	if ( name_ == "P6222") return "P 62 2 2";
	if ( name_ == "P6422") return "P 64 2 2";
	if ( name_ == "P6322") return "P 63 2 2";
	if ( name_ == "P222") return "P 2 2 2";
	if ( name_ == "P2221") return "P 2 2 21";
	if ( name_ == "P21212") return "P 21 21 2";
	if ( name_ == "P212121") return "P 21 21 21";
	if ( name_ == "C2221") return "C 2 2 21";
	if ( name_ == "C222") return "C 2 2 2";
	if ( name_ == "F222") return "F 2 2 2";
	if ( name_ == "I222") return "I 2 2 2";
	if ( name_ == "I212121") return "I 21 21 21";
	return "?";
}
