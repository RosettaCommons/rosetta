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
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/cryst/spacegroup.hh>
#include <protocols/cryst/refinable_lattice.hh>

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
#include <cmath>

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

static basic::Tracer TR("cryst.gen");

////////////////////////////////////////////////
// helper functions
inline int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
inline Real pos_mod(Real x,Real y) {
	Real r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y;
	if ( r<-y/2 ) r+=y;
	if ( r>=y/2 ) r-=y;
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y);
	if ( r<-0.5*y ) r+=y;
	if ( r>=0.5*y ) r-=y;
	return r;
}

inline double sign(double x) {return (x >= 0.0) ? 1.0 : -1.0;}
inline double norm4(double a, double b, double c, double d) {return sqrt(a * a + b * b + c * c + d * d);}

inline numeric::xyzMatrix<Real>
axis_angle( numeric::xyzVector<Real> u, Real theta) {
	u.normalize();
	Real ct = cos(theta);
	Real st = sin(theta);

	return numeric::xyzMatrix<Real>::rows(
		ct+u[0]*u[0]*(1-ct)     , u[0]*u[1]*(1-ct)-u[2]*st, u[0]*u[2]*(1-ct)+u[1]*st ,
		u[0]*u[1]*(1-ct)+u[2]*st, ct+u[1]*u[1]*(1-ct)     , u[1]*u[2]*(1-ct)-u[0]*st ,
		u[0]*u[2]*(1-ct)-u[1]*st, u[1]*u[2]*(1-ct)+u[0]*st, ct+u[2]*u[2]*(1-ct)
	);
}

inline bool transforms_equiv_mod1(
	numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
	numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() ) +
		std::fabs( min_mod( T1[0] - T2[0] , 1.0)) +
		std::fabs( min_mod( T1[1] - T2[1] , 1.0)) +
		std::fabs( min_mod( T1[2] - T2[2] , 1.0));
	return (err <= 1e-6);
}

inline bool transforms_equiv(
	numeric::xyzMatrix<Real> const &S1, numeric::xyzVector<Real> const &T1,
	numeric::xyzMatrix<Real> const &S2, numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() ) +
		std::fabs( ( T1[0] - T2[0] )) +
		std::fabs( ( T1[1] - T2[1] ))  +
		std::fabs( ( T1[2] - T2[2] ));
	return (err <= 1e-6);
}

inline bool transforms_equiv(
	numeric::xyzMatrix<Real> const &S1,
	numeric::xyzMatrix<Real> const &S2
) {
	Real err =
		std::fabs( S1.xx() - S2.xx() ) +  std::fabs( S1.xy() - S2.xy() ) +  std::fabs( S1.xz() - S2.xz() ) +
		std::fabs( S1.yx() - S2.yx() ) +  std::fabs( S1.yy() - S2.yy() ) +  std::fabs( S1.yz() - S2.yz() ) +
		std::fabs( S1.zx() - S2.zx() ) +  std::fabs( S1.zy() - S2.zy() ) +  std::fabs( S1.zz() - S2.zz() );
	return (err <= 1e-6);
}

inline bool transforms_equiv(
	numeric::xyzVector<Real> const &T1,
	numeric::xyzVector<Real> const &T2
) {
	Real err =
		std::fabs( ( T1[0] - T2[0] )) +
		std::fabs( ( T1[1] - T2[1] )) +
		std::fabs( ( T1[2] - T2[2] ));
	return (err <= 1e-6);
}

void
disp_symmop( core::kinematics::RT const &rt ) {
	numeric::xyzMatrix<Real> const &S=rt.get_rotation();
	numeric::xyzVector<Real> const &T=rt.get_translation();

	std::cout << "S=["<<S.xx()<<","<<S.xy()<<","<<S.xz()<<";"
		<<S.yx()<<","<<S.yy()<<","<<S.yz()<<";"
		<<S.zx()<<","<<S.zy()<<","<<S.zz()<<"]  T=["
		<<T[0]<<","<<T[1]<<","<<T[2]<<"]" << std::endl;
}


// linalg utility
numeric::xyzVector<Real>
get_rotation_axis( numeric::xyzMatrix<Real> const& R) {
	numeric::xyzVector<Real> axis;

	Real S;
	if ( R.xx() + R.yy() + R.zz() > 0 ) {
		S = sqrt(1 + R.xx() + R.yy() + R.zz()) * 2;
		axis[0] = ( R.zy() - R.yz() ) * S;
		axis[1] = ( R.xz() - R.zx() ) * S;
		axis[2] = ( R.yx() - R.xy() ) * S;
	} else if ( R.xx() > R.yy() && R.xx() > R.zz() )  {
		S  = sqrt( 1.0 + R.xx() - R.yy() - R.zz() ) * 2;
		axis[0] = 0.25 * S;
		axis[1] = (R.yx() + R.xy() ) / S;
		axis[2] = (R.zx() + R.xz() ) / S;
	} else if ( R.yy() > R.zz() ) {
		S  = sqrt( 1.0 - R.xx() + R.yy() - R.zz() ) * 2;
		axis[0] = (R.yx() + R.xy() ) / S;
		axis[1] = 0.25 * S;
		axis[2] = (R.zy() + R.yz() ) / S;
	} else {
		S  = sqrt( 1.0 - R.xx() - R.yy() + R.zz() ) * 2;
		axis[0] = (R.xz() + R.zx() ) / S;
		axis[1] = (R.zy() + R.yz() ) / S;
		axis[2] = 0.25 * S;
	}
	axis.normalize();
	return axis;
}

// linalg utility
numeric::xyzVector<Real>
get_reflection_axis( numeric::xyzMatrix<Real> const& R) {
	numeric::xyzVector<Real> axis;

	numeric::xyzVector< core::Real > eigval;
	numeric::xyzMatrix< core::Real > eigvec;
	eigval = numeric::eigenvector_jacobi( R, 1e-6, eigvec );

	// special case: Ci symmetry
	if ( std::fabs( eigval[0] + 1 ) < 1e-6 && std::fabs( eigval[1] + 1 ) < 1e-6 && std::fabs( eigval[2] + 1 ) < 1e-6 ) {
		// inversion about origin, we'll set axis to (0,0,0)
		return numeric::xyzVector<Real>(0,0,0);
	}

	// otherwise, mirror plane ...  find the eigvec corresponding to eigval == -1
	if ( std::fabs( eigval[0] + 1 ) < 1e-6 ) {
		axis = eigvec.col(1); // 1 indexed (?)
	} else if ( std::fabs( eigval[1] + 1 ) < 1e-6 ) {
		axis = eigvec.col(2);
	} else if ( std::fabs( eigval[2] + 1 ) < 1e-6 ) {
		axis = eigvec.col(3);
	} else {
		std::cout << "eigvals: " << eigval[0] << "," << eigval[1] << "," << eigval[2] << std::endl;
		utility_exit_with_message("error in get_reflection_axis");
	}
	axis.normalize();
	return axis;
}

class pointGroupHit {
public:
	friend void
	get_angle_and_offset(
		pointGroupHit hit_i, pointGroupHit hit_j,
		core::Real &angle, core::Real &offset, core::Real &shift
	);

public:
	pointGroupHit() {}

	pointGroupHit(
		std::string name_in,
		numeric::xyzVector<Real> origin_in,
		numeric::xyzVector<Real> axis_in
	) {
		name_=name_in;
		origin_=origin_in;
		axes.push_back(axis_in);
		resolve_name();
	}

	pointGroupHit(
		std::string name_in,
		numeric::xyzVector<Real> origin_in,
		utility::vector1< numeric::xyzVector<Real> > axes_in
	) {
		name_=name_in;
		origin_=origin_in;
		axes = axes_in;
		resolve_name();
	}

	pointGroupHit(
		pointGroupHit pg1,
		pointGroupHit pg2,
		utility::vector1<core::kinematics::RT> const &rts
	);

	void
	construct_from_basis(
		pointGroupHit hit_i,
		pointGroupHit hit_j,
		utility::vector1<core::kinematics::RT> const &rts,
		core::Size expected_count
	);


	void
	add_symmop( int symm, numeric::xyzVector<int> ls) {
		symmops.push_back(symm);
		lattice_shifts.push_back(ls);
	}

	void append_symmops( utility::vector1<core::kinematics::RT> &rts_symm, utility::vector1<core::kinematics::RT> const &rts ) {
		for ( core::Size i=1; i<=symmops.size(); ++i ) {
			rts_symm.push_back( core::kinematics::RT(
				rts[symmops[i]].get_rotation(),
				rts[symmops[i]].get_translation() +
				numeric::xyzVector<Real> ( lattice_shifts[i][0],  lattice_shifts[i][1],  lattice_shifts[i][2] ) )
			);
		}
	}

	bool
	expand(utility::vector1<core::kinematics::RT> const &rts) {
		bool done=false;
		int rounds=0;
		int nsymmops0=symmops.size();
		while ( !done && ++rounds<16 ) {
			done=true;
			int nsymmops=symmops.size();
			for ( int i=1; i<=nsymmops; ++i ) {
				for ( int j=1; j<=nsymmops0; ++j ) {
					numeric::xyzMatrix<Real> const &Si=rts[symmops[i]].get_rotation();
					numeric::xyzVector<Real> Ti = rts[symmops[i]].get_translation()
						+ numeric::xyzVector<Real>(lattice_shifts[i][0],
						lattice_shifts[i][1],
						lattice_shifts[i][2]) ;
					numeric::xyzMatrix<Real> const &Sj=rts[symmops[j]].get_rotation();
					numeric::xyzVector<Real> Tj = rts[symmops[j]].get_translation()
						+ numeric::xyzVector<Real>(lattice_shifts[j][0],
						lattice_shifts[j][1],
						lattice_shifts[j][2]) ;

					numeric::xyzMatrix<Real> Sij = Sj*Si;
					numeric::xyzVector<Real> Tij = Sj*Ti + Tj;

					// check if unique
					bool unique = true;
					for ( int k=1; k<=(int)symmops.size() && unique; ++k ) {
						if ( transforms_equiv(Sij, Tij, rts[symmops[k]].get_rotation(),
								rts[symmops[k]].get_translation()
								+ numeric::xyzVector<Real>(lattice_shifts[k][0],
								lattice_shifts[k][1],
								lattice_shifts[k][2]) ) ) {
							unique = false;
						}
					}

					// look up spacegroup
					if ( unique ) {
						bool found = false;
						for ( int k=1; k<=(int)rts.size() && !found; ++k ) {
							if ( transforms_equiv_mod1(Sij, Tij, rts[k].get_rotation(), rts[k].get_translation()) ) {
								symmops.push_back(k);
								numeric::xyzVector<int> ls(
									(int)std::floor( Tij[0]-rts[k].get_translation()[0]+0.5 ),
									(int)std::floor( Tij[1]-rts[k].get_translation()[1]+0.5 ),
									(int)std::floor( Tij[2]-rts[k].get_translation()[2]+0.5 )
								);
								lattice_shifts.push_back(ls);
								found = true;
								done = false;
							}
						}

						if ( !found ) {
							std::cerr << "In " << name_ << std::endl;
							utility_exit_with_message("error in expand (symmop not found)");
						}
					}
				}
			}
		}

		return done;
	}

	void
	show(std::ostream & oss) const {
		oss << "  SYMMOP " << name_ << ": nsymm=" << symmops.size();
		oss << " origin=[" << origin_[0] << "," << origin_[1] << "," << origin_[2] << "] ";
		for ( int i=1; i<=std::min((int)axes.size(),2); ++i ) {
			oss << "axis_" << i << "=[" << axes[i][0] << "," << axes[i][1] << "," << axes[i][2] << "] ";
		}
		oss << "symmops=" << symmops[1] ;
		for ( int i=2; i<=(int)symmops.size(); ++i ) {
			oss << "," << symmops[i];
		}
		oss << std::endl;
	}

	void
	show() const {
		show(std::cout);
	}

	pointGroupHit
	transform( utility::vector1<core::kinematics::RT> const &rts, core::kinematics::RT const &rt ) {
		std::string name_new = name_;
		numeric::xyzVector<Real> origin_new = rt.get_rotation()*origin() + rt.get_translation();
		utility::vector1< numeric::xyzVector<Real> > axes_new = axes;

		for ( int i=1; i<=(int)axes.size(); ++i ) {
			axes_new[i] = rt.get_rotation()*axes[i];
		}

		pointGroupHit tX(name_new, origin_new, axes_new);

		for ( int i=1; i<=(int)symmops.size(); ++i ) {
			numeric::xyzMatrix<Real> const &Si=
				rt.get_rotation()*rts[symmops[i]].get_rotation();
			numeric::xyzVector<Real> Ti =
				rt.get_rotation() * rts[symmops[i]].get_translation()
				+ rt.get_translation()
				+ numeric::xyzVector<Real>(lattice_shifts[i][0],
				lattice_shifts[i][1],
				lattice_shifts[i][2]) ;

			bool found=false;
			for ( int k=1; k<=(int)rts.size() && !found; ++k ) {
				if ( transforms_equiv_mod1(Si, Ti, rts[k].get_rotation(), rts[k].get_translation()) ) {
					tX.symmops.push_back(k);
					numeric::xyzVector<int> ls(
						(int)std::floor( Ti[0]-rts[k].get_translation()[0]+0.5 ),
						(int)std::floor( Ti[1]-rts[k].get_translation()[1]+0.5 ),
						(int)std::floor( Ti[2]-rts[k].get_translation()[2]+0.5 )
					);
					tX.lattice_shifts.push_back(ls);
					found=true;
				}
			}
		}
		runtime_assert( tX.symmops.size() == symmops.size() );
		return tX;
	}

	int intersects( pointGroupHit const & pg ) const {
		int retval=0;
		int nsymmops0=symmops.size();
		int nsymmops1=pg.symmops.size();
		for ( int i=1; i<=nsymmops0; ++i ) {
			for ( int j=1; j<=nsymmops1; ++j ) {
				if ( symmops[i] == pg.symmops[j] && lattice_shifts[i] == pg.lattice_shifts[j] ) {
					retval++;
				}
			}
		}
		return retval;
	}

	bool is_subset_of( pointGroupHit const & pg ) const {
		int nsymmops0=symmops.size();
		int nsymmops1=pg.symmops.size();
		for ( int i=1; i<=nsymmops0; ++i ) {
			bool contained_in = false;
			for ( int j=1; j<=nsymmops1 && !contained_in; ++j ) {
				contained_in = (symmops[i] == pg.symmops[j] && lattice_shifts[i] == pg.lattice_shifts[j]);
			}
			if ( !contained_in ) return false;
		}
		return true;
	}

	bool equals( pointGroupHit const & pg ) const {
		return (is_subset_of(pg) && pg.is_subset_of(*this));
	}

	bool is_subset_of_ignore_shift( pointGroupHit const & pg ) const {
		int nsymmops0=symmops.size();
		int nsymmops1=pg.symmops.size();
		for ( int i=1; i<=nsymmops0; ++i ) {
			bool contained_in = false;
			for ( int j=1; j<=nsymmops1 && !contained_in; ++j ) {
				contained_in = (symmops[i] == pg.symmops[j]);
			}
			if ( !contained_in ) return false;
		}
		return true;
	}

	bool equals_ignore_shift( pointGroupHit const & pg ) const {
		return (this->is_subset_of_ignore_shift(pg) && pg.is_subset_of_ignore_shift(*this));
	}

	int order() const {
		return (symmops.size());
	}

	std::string name() const { return name_; }

	Real
	score() {
		Real score=0;
		score += std::fabs(origin_[0]) + std::fabs(origin_[1]) + std::fabs(origin_[2]);
		return score;
	}

	std::string
	as_string( ) {
		std::ostringstream oss;
		oss << axes[1][0] <<","<< axes[1][1] <<","<< axes[1][2] << " ";
		if ( axes.size() == 2 ) {
			oss << axes[2][0] <<","<< axes[2][1] <<","<< axes[2][2] << " ";
		} else {
			oss << "- ";
		}
		oss << origin()[0] <<","<< origin()[1] <<","<< origin()[2] <<" ";
		return (oss.str());
	}

	bool isC() const { return isC_; }
	bool isD() const { return isD_; }
	bool isT() const { return isT_; }
	bool isO() const { return isO_; }
	bool isCs() const { return isCs_; }
	bool isCi() const { return isCi_; }
	bool isS() const { return isS_; }
	bool isCh() const { return isCh_; }
	bool isCv() const { return isCv_; }
	bool isDd() const { return isDd_; }
	bool isDh() const { return isDh_; }

	// getters
	numeric::xyzVector<Real> & origin() { return origin_; }
	numeric::xyzVector<Real> const & origin() const { return origin_; }

private:
	std::string name_;

	void resolve_name() {
		isC_=isD_=isT_=isO_=isCs_=isCi_=isS_=isCh_=isCv_=isDd_=isDh_=false;
		if ( name_ == "C1" || name_ == "C2" || name_ == "C3" || name_ == "C4" || name_ == "C6" ) {
			isC_ = true;
		} else if ( name_ == "D2" || name_ == "D3" || name_ == "D4" || name_ == "D6" ) {
			isD_ = true;
		} else if ( name_ == "T" ) {
			isT_ = true;
		} else if ( name_ == "O" ) {
			isO_ = true;
		} else if ( name_ == "Cs" ) {
			isCs_ = true; // Cs == C1v == C1h
		} else if ( name_ == "Ci" ) {
			isCi_ = true; // Ci == S2
		} else if ( name_ == "S2" || name_ == "S4" || name_ == "S6" || name_ == "S8" || name_ == "S12" ) {
			isS_ = true;
		} else if ( name_ == "C2h" || name_ == "C3h" || name_ == "C4h" || name_ == "C6h" ) {
			isCh_ = true;
		} else if ( name_ == "C2v" || name_ == "C3v" || name_ == "C4v" || name_ == "C6v" ) {
			isCv_ = true;
		} else if ( name_ == "D2d" || name_ == "D3d" || name_ == "D4d" || name_ == "D6d" ) {
			isDd_ = true;
		} else if ( name_ == "D2h" || name_ == "D3h" || name_ == "D4h" || name_ == "D6h" ) {
			isDh_ = true;
		}
	}

	bool isC_, isD_, isT_, isO_; // chiral
	bool isCs_, isCi_; // simple achiral
	bool isS_, isCh_, isCv_, isDd_, isDh_; // more complicated achiral

	utility::vector1<int> symmops;
	utility::vector1< numeric::xyzVector<int> > lattice_shifts;
	utility::vector1< numeric::xyzVector<Real> > axes;
	numeric::xyzVector<Real> origin_;
};


struct latticeHit {
	latticeHit( pointGroupHit hit1_in, pointGroupHit hit2_in, Real angle_in, Real offset_in, Real shift_in ) {
		hit1=hit1_in; hit2=hit2_in; angle=angle_in; offset=offset_in; shift=shift_in;
	}

	pointGroupHit hit1,hit2;
	Real angle;
	Real offset;
	Real shift;
};

void
get_angle_and_offset(
	pointGroupHit hit_i, pointGroupHit hit_j,
	core::Real &angle, core::Real &offset, core::Real &shift
) {
	numeric::xyzVector<core::Real> Iaxis=hit_i.axes[1], Jaxis=hit_j.axes[1];
	angle = RAD2DEG * acos( std::max( std::min( Iaxis.dot( Jaxis ), 1.0 ), -1.0) );
	angle = std::min( 180-angle, angle );
	offset=0;
	shift=0;

	if ( angle <= 1e-6 ) {
		offset = ((hit_j.origin()-hit_i.origin()) - (hit_j.origin()-hit_i.origin()).dot( Iaxis )* Iaxis).length();
		shift = (hit_j.origin()-hit_i.origin()).dot( Iaxis );
	} else {
		offset = std::fabs (
			(hit_j.origin()-hit_i.origin()).dot( Iaxis.cross( Jaxis ) ) / Iaxis.cross( Jaxis ).length() );
	}
}


void
pointGroupHit::construct_from_basis(
	pointGroupHit hit_i,
	pointGroupHit hit_j,
	utility::vector1<core::kinematics::RT> const &rts,
	core::Size expected_count
) {
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	numeric::xyzVector<core::Real> zerovec(0,0,0);

	axes.push_back(hit_j.axes[1]);
	axes.push_back(hit_i.axes[1]);

	for ( int i=1; i<=hit_i.order()+hit_j.order(); ++i ) {
		pointGroupHit const &h = (i<=hit_i.order() ? hit_i : hit_j);
		int idx = (i<=hit_i.order() ? i : i-hit_i.order());

		numeric::xyzMatrix<Real> Si = rts[h.symmops[idx]].get_rotation();
		numeric::xyzVector<Real> Ti = rts[h.symmops[idx]].get_translation()
			+ numeric::xyzVector<Real>(h.lattice_shifts[idx][0], h.lattice_shifts[idx][1], h.lattice_shifts[idx][2]) ;

		if ( !transforms_equiv( Si, Ti, identity, zerovec ) ) {
			symmops.push_back(h.symmops[idx]);
			lattice_shifts.push_back(h.lattice_shifts[idx]);
		}
	}
	// now add in identity
	symmops.push_back( 1 ); // op 1 always identity
	lattice_shifts.push_back( numeric::xyzVector<int>(0,0,0) );

	expand( rts );
	if ( symmops.size() != expected_count ) {
		utility_exit_with_message("error in expand (wrong # of symmops)");
	}
}


// expand pairs of point-group hits
pointGroupHit::pointGroupHit(
	pointGroupHit pg1,
	pointGroupHit pg2,
	utility::vector1<core::kinematics::RT> const &rts
) {
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	name_="";

	// possible expansions:
	//   a) C2+Cn -> Dn
	//   b) Cs+Cn -> Cnv or Cnh
	//   c) Cs+Sn -> Dnd
	//   d) Cs+Cnv -> Dnh (?)
	//   e) C2+C3 -> T or O
	//   *** f) S4+C3 -> Td (?)   TO DO
	//   *** g) S6+Ci -> Th or Oh (?)    TO DO

	// cases a & e
	if ( pg1.name() == "C2" || pg2.name() == "C2" ) {
		pointGroupHit const &hit_i = (pg1.name() == "C2") ? pg1 : pg2;
		pointGroupHit const &hit_j = (pg1.name() == "C2") ? pg2 : pg1;

		if ( hit_j.isC() ) {
			numeric::xyzMatrix<Real> Si = rts[hit_i.symmops[1]].get_rotation();
			numeric::xyzVector<Real> Ti = rts[hit_i.symmops[1]].get_translation()
				+ numeric::xyzVector<Real>(hit_i.lattice_shifts[1][0], hit_i.lattice_shifts[1][1], hit_i.lattice_shifts[1][2]) ;
			numeric::xyzMatrix<Real> Sj = rts[hit_j.symmops[1]].get_rotation();
			numeric::xyzVector<Real> Tj = rts[hit_j.symmops[1]].get_translation()
				+ numeric::xyzVector<Real>(hit_j.lattice_shifts[1][0], hit_j.lattice_shifts[1][1], hit_j.lattice_shifts[1][2]) ;

			// *** Dn SYMMETRY ***
			Real dotprod = hit_i.axes[1].dot(hit_j.axes[1] );
			if ( std::fabs(dotprod) <=1e-4 ) {
				numeric::xyzVector<Real> del = Si*Tj+Ti-(Sj*Ti+Tj);
				if ( transforms_equiv(del,numeric::xyzVector<core::Real>(0,0,0) ) ) {
					name_ = "D"+hit_j.name_.substr(1,1);
					origin_ = (hit_j.origin() + Si*hit_j.origin() + Ti)/2.0;
					construct_from_basis( hit_i,hit_j, rts, 2*hit_j.order() );
				}
			}

			// *** T/O SYMMETRY ***
			// if (hit_j.name() == "C3") {
			//  Real b = hit_i.axes[1].dot(hit_j.axes[1] );
			//
			//  bool is_potential_T = (std::abs( b-0.57735026919) <= 1e-4);
			//  bool is_potential_O = (std::abs( b-0.81649658093) <= 1e-4);
			//
			//  if (is_potential_T || is_potential_O) {
			//   Real angle = RAD2DEG*acos( b );
			//
			//   numeric::xyzVector<Real> w0=hit_i.origin()-hit_j.origin();
			//   Real d=dot(hit_j.axes[1],w0);
			//   Real e=dot(hit_i.axes[1],w0);
			//   numeric::xyzVector<Real> x = (hit_i.origin()-hit_j.origin())+((b*e-d)*hit_i.axes[1]-(e-b*d)*hit_j.axes[1])/(1-b*b);
			//   Real offset = x.length();
			//   Real sc=(b*e-d)/(1-b*b);
			//
			//   bool intersects = ( offset<=1e-4 );  // do symmaxes intersect?
			//
			//   numeric::xyzVector<Real> next3fold = axis_angle(hit_i.axes[1], 180*DEG2RAD )*hit_j.axes[1];
			//   if ( intersects && is_potential_T ) {
			//    name_ = "T";
			//    origin_ = hit_i.origin()+hit_i.axes[1]*sc;
			//
			//    // rotate by 120 around axis 1
			//    numeric::xyzMatrix<Real> axisgen = axis_angle( hit_j.axes[1], 120*DEG2RAD );
			//    axes.push_back(hit_j.axes[1]);
			//    axes.push_back(next3fold);
			//    axes.push_back(axisgen*hit_j.axes[1]);
			//    axes.push_back(axisgen*axisgen*next3fold);
			//    construct_from_basis( hit_i,hit_j, rts, 12 );
			//   } else if ( intersects && is_potential_O ) {
			//    name_ = "O"; //targetct=24;
			//    origin_ = hit_i.origin()+hit_i.axes[1]*sc;
			//
			//    // rot 90 about axis 1
			//    numeric::xyzMatrix<Real> axisgen = axis_angle( hit_j.axes[1], 120*DEG2RAD );
			//    axes.push_back(hit_j.axes[1]);
			//    axes.push_back(next3fold);
			//    axes.push_back(-hit_j.axes[1]);
			//    axes.push_back(-next3fold);
			//    axes.push_back(axisgen*hit_j.axes[1]);
			//    axes.push_back(-(axisgen*next3fold));
			//    axes.push_back(axisgen*axisgen*hit_j.axes[1]);
			//    axes.push_back(-(axisgen*axisgen*next3fold));
			//    construct_from_basis( hit_i,hit_j, rts, 24 );
			//
			//   }
			//  }
			// }
		}
	}

	// cases b & c & e
	if ( pg1.name() == "Cs" || pg2.name() == "Cs" ) {
		pointGroupHit const &hit_i = (pg1.name() == "Cs") ? pg1 : pg2;
		pointGroupHit const &hit_j = (pg1.name() == "Cs") ? pg2 : pg1;

		numeric::xyzMatrix<Real> Si = rts[hit_i.symmops[1]].get_rotation();
		numeric::xyzVector<Real> Ti = rts[hit_i.symmops[1]].get_translation()
			+ numeric::xyzVector<Real>(hit_i.lattice_shifts[1][0], hit_i.lattice_shifts[1][1], hit_i.lattice_shifts[1][2]) ;
		numeric::xyzMatrix<Real> Sj = rts[hit_j.symmops[1]].get_rotation();
		numeric::xyzVector<Real> Tj = rts[hit_j.symmops[1]].get_translation()
			+ numeric::xyzVector<Real>(hit_j.lattice_shifts[1][0], hit_j.lattice_shifts[1][1], hit_j.lattice_shifts[1][2]) ;

		if ( hit_j.isC() ) {
			// *** Cnv/Cnh SYMMETRY ***
			Real dotprod = hit_i.axes[1].dot( hit_j.axes[1] );
			if ( std::fabs(dotprod-1) <=1e-4 ) {
				// mirror plane normal parallel to symm axis ==> Cnh
				// no translation check needed
				name_ = "C"+hit_j.name_.substr(1,1)+"h";
				// origin is intersection of Cn axis and mirror plane
				core::Real d = (hit_i.origin()-hit_j.origin()).dot( hit_i.axes[1] );
				origin_ = d*hit_j.axes[1] + hit_j.origin();
				construct_from_basis( hit_i,hit_j, rts, 2*hit_j.order() );
			} else if ( std::fabs(dotprod) <=1e-4 ) {
				// mirror plane normal perpendicular to symm axis ==> Cnv
				// ensure mirror plane intersects symm axis
				numeric::xyzVector<Real> Tij = hit_i.origin()-hit_j.origin();
				numeric::xyzVector<Real> del = Si*Tj+Ti-(Sj*Ti+Tj);

				if ( std::fabs( Tij.dot(hit_i.axes[1]) ) < 1e-6 && transforms_equiv(del,numeric::xyzVector<core::Real>(0,0,0)) ) {
					name_ = "C"+hit_j.name_.substr(1,1)+"v";
					// origin is origin of the Cn
					origin_ = (hit_j.origin() + Si*hit_j.origin() + Ti)/2.0;
					construct_from_basis( hit_i,hit_j, rts, 2*hit_j.order() );
				}
			}
		} else if ( hit_j.isS() ) {
			// *** Dnd SYMMETRY ***
			// check mirror plane normal perpendicular to symm axis
			Real dotprod = hit_i.axes[1].dot( hit_j.axes[1] );
			if ( std::fabs(dotprod) <= 1e-4 ) {
				// mirror plane normal perpendicular to symm axis ==> Dnd
				// ensure mirror plane intersects symm axis
				numeric::xyzVector<Real> Tij = hit_i.origin()-hit_j.origin();
				if ( std::fabs( Tij.dot(hit_i.axes[1]) ) < 1e-6 ) {
					name_ = "D"+hit_j.name_.substr(1,1)+"d";
					// origin is origin of the S
					origin_ = (hit_j.origin() + Si*hit_j.origin() + Ti)/2.0;
					construct_from_basis( hit_i,hit_j, rts, 2*hit_j.order() );
				}
			}
		} else if ( hit_j.isCv() ) {
			// *** Dnh SYMMETRY ***
			// a bit hairy since this is composed from 3 primative ops so there are several
			//     different ways we can generate this...
			// we'll just try this one for now (Cs+Cnv)
			//
			// hit_j, the Cnv symm group: axis[1] is Cn, axis2 is Cs
			Real dotprod = hit_i.axes[1].dot( hit_j.axes[1] );
			if ( std::fabs(dotprod-1) <=1e-4 ) {
				// mirror plane normal is parallel to the Cn axis
				name_ = "D"+hit_j.name_.substr(1,1)+"h";
				origin_ = (hit_j.origin() + Si*hit_j.origin() + Ti)/2.0;
				construct_from_basis( hit_i,hit_j, rts, 2*hit_j.order() );
			}
		}
	}

	// cases f & g ::: TO DO


	resolve_name(); // final call
}


void
dedup( utility::vector1< pointGroupHit > &hits ) {
	utility::vector1< pointGroupHit > hits_dedup;

	// remove duplicates
	int nhits = (int)hits.size();
	for ( int i=1; i<=nhits; ++i ) {
		bool uniq = true;
		int j;
		for ( j=1; j<=(int)hits_dedup.size() && uniq; ++j ) {
			if ( (hits[i].origin() - hits_dedup[j].origin()).length() < 1e-4 && hits[i].equals_ignore_shift( hits_dedup[j] ) ) {
				uniq=false;
			}

		}
		if ( uniq ) {
			hits_dedup.push_back(hits[i]);
		} else {
			if ( hits_dedup[j-1].score() > hits[i].score() ) {
				hits_dedup[j-1] = hits[i];
			}
		}
	}

	hits = hits_dedup;
	hits_dedup.clear();

	// [ii] subsets
	nhits = (int)hits.size();
	for ( int i=1; i<=nhits; ++i ) {
		bool uniq = true;
		for ( int j=1; j<=nhits && uniq; ++j ) {
			if ( i==j ) continue;

			// don't generate embeddings of same point symmetry
			if ( hits[i].isC() && hits[j].isC() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;
			if ( hits[i].isCh() && hits[j].isCh() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;
			if ( hits[i].isCv() && hits[j].isCv() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;
			if ( hits[i].isD() && hits[j].isD() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;
			if ( hits[i].isDd() && hits[j].isDd() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;
			if ( hits[i].isDh() && hits[j].isDh() && hits[i].order() < hits[j].order() && hits[i].is_subset_of( hits[j] ) ) uniq=false;

			// don't generate T embedded in O
			if ( hits[i].isT() && hits[j].isO() && hits[i].is_subset_of( hits[j] ) ) uniq=false;

			// we could further reduce this, but we might not want to
			// e.g. could kill C embedded in D
		}
		if ( uniq ) hits_dedup.push_back(hits[i]);
	}
	hits = hits_dedup;
}

void
remove_drift( utility::vector1< pointGroupHit > &hits ) {
	utility::vector1< pointGroupHit > hits_dedup;
	int nhits = (int)hits.size();
	for ( int i=1; i<=nhits; ++i ) {
		if (
				hits[i].origin()[0] <= 1 || hits[i].origin()[0] > -1 ||
				hits[i].origin()[1] <= 1 || hits[i].origin()[1] > -1 ||
				hits[i].origin()[2] <= 1 || hits[i].origin()[2] > -1 ) {
			hits_dedup.push_back (hits[i]);
		}
	}
	hits = hits_dedup;
}

// get a point symmetries from a single symmop
pointGroupHit
get_primary_point_group(
	numeric::xyzMatrix<Real> const & Si,
	numeric::xyzVector<Real> const & Ti,
	numeric::xyzMatrix<core::Real> const &skewM
) {
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);

	// 1: are we a mirror symm?
	bool is_inversion = (Si.det()<0);

	// keeping track of variables:
	//   (Ti,Si) = base transform
	//   SS,TT = running totals of Si, Ti
	//   Tsum = running total of n*Ti+(n-1)*SiTi+(n-2)*Si^2Ti+...
	numeric::xyzMatrix<Real> ST=Si;
	numeric::xyzVector<Real> TT=Ti, Tsum=Ti;

	if ( is_inversion ) {
		ST = Si*Si;
		TT = Ti + Si*Ti;

		if ( transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0)) ) {
			numeric::xyzVector<Real> invaxis = get_reflection_axis( skewM*Si*numeric::inverse(skewM) );
			if ( invaxis[0] == 0.0 && invaxis[1] == 0.0 && invaxis[2] == 0.0 ) {
				return (pointGroupHit("Ci", Ti/2.0, invaxis ));
			} else {
				return (pointGroupHit("Cs", Ti/2.0, invaxis ));
			}
		}
		if ( transforms_equiv(ST,identity) ) return (pointGroupHit());  // shift only

		Tsum += TT;
	}

	// get rotation axis
	numeric::xyzVector<Real> axis = get_rotation_axis( skewM*ST*numeric::inverse(skewM) );

	// C2/3/4/6 or S4/6/8/12
	ST = Si*ST;
	TT = Ti + Si*Ti;
	Tsum += TT;
	if ( is_inversion ) {
		ST = Si*ST;
		TT = Ti + Si*Ti;
		Tsum += TT;
	}

	if ( transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0)) ) {
		if ( is_inversion ) {
			return (pointGroupHit("S4", Tsum/4.0, axis ));
		} else {
			return (pointGroupHit("C2", Tsum/2.0, axis ));
		}
	}

	ST = Si*ST;
	TT = Ti + Si*Ti;
	Tsum += TT;
	if ( is_inversion ) {
		ST = Si*ST;
		TT = Ti + Si*Ti;
		Tsum += TT;
	}
	if ( transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0)) ) {
		if ( is_inversion ) {
			return (pointGroupHit("S6", Tsum/6.0, axis ));
		} else {
			return (pointGroupHit("C3", Tsum/3.0, axis ));
		}
	}

	ST = Si*ST;
	TT = Ti + Si*Ti;
	Tsum += TT;
	if ( transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0)) ) {
		return (pointGroupHit("C4", Tsum/4.0, axis ));
	}

	ST = Si*ST;
	TT = Ti + Si*Ti;
	Tsum += TT;
	// dont care about 5-folds

	ST = Si*ST;
	TT = Ti + Si*Ti;
	Tsum += TT;
	if ( transforms_equiv(ST,TT,identity,numeric::xyzVector<core::Real>(0,0,0)) ) {
		return (pointGroupHit("C6", Tsum/6.0, axis ));
	}
	return (pointGroupHit());
}


// get all point symmetries from symmops
utility::vector1< pointGroupHit >
get_point_groups(utility::vector1<core::kinematics::RT> const &rts, numeric::xyzMatrix<core::Real> const &skewM ) {
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
	utility::vector1< pointGroupHit > hits;
	int nxforms = (int)rts.size();

	// [1] point symm around (0,0,0)
	int MAXS = 1;
	for ( int x=-MAXS; x<=MAXS; ++x ) {
		for ( int y=-MAXS; y<=MAXS; ++y ) {
			for ( int z=-MAXS; z<=MAXS; ++z ) {
				for ( int i=2; i<=nxforms; ++i ) {
					numeric::xyzMatrix<Real> const &Si=rts[i].get_rotation();
					numeric::xyzVector<Real> lattice_shift = numeric::xyzVector<Real>(x,y,z);
					numeric::xyzVector<Real> Ti = rts[i].get_translation() + lattice_shift;

					if ( Ti[0] > 1 || Ti[0] < -1 ) continue;
					if ( Ti[1] > 1 || Ti[1] < -1 ) continue;
					if ( Ti[2] > 1 || Ti[2] < -1 ) continue;

					if ( transforms_equiv(Si,identity) ) continue;

					pointGroupHit pg = get_primary_point_group(Si, Ti, skewM);
					if ( pg.name() != "" ) {
						pg.add_symmop( i,lattice_shift );
						pg.expand( rts );
						hits.push_back( pg);
					}
				}
			}
		}
	}

	// remove drift
	remove_drift(hits);
	dedup(hits);

	// detect higher order
	int nhits = (int)hits.size();
	for ( int i=1; i<=nhits; ++i ) {
		for ( int j=i+1; j<=nhits; ++j ) {
			pointGroupHit const &hit_i = hits[i];
			pointGroupHit const &hit_j = hits[j];

			if ( hit_i.intersects(hit_j) > 1 ) continue;

			pointGroupHit hit_ij(hit_i, hit_j, rts);
			if ( hit_ij.name() != "" ) {
				hits.push_back( hit_ij );
			}
		}
	}

	// need to do this again... (3rd order interactions)
	//   however, we only need to consider Cs+Cnv
	nhits = (int)hits.size();
	for ( int i=1; i<=nhits; ++i ) {
		for ( int j=i+1; j<=nhits; ++j ) {
			pointGroupHit const &hit_i = hits[i];
			pointGroupHit const &hit_j = hits[j];

			if ( hit_i.intersects(hit_j) > 1 ) continue;

			if ( (hit_i.isCs() && hit_j.isCv()) || (hit_j.isCs() && hit_i.isCv()) ) {
				pointGroupHit hit_ij(hit_i, hit_j, rts);
				if ( hit_ij.name() != "" ) {
					hits.push_back( hit_ij );
				}
			}
		}
	}


	remove_drift(hits);
	dedup(hits);

	return hits;
}


bool
check_if_forms_lattice(
	utility::vector1<core::kinematics::RT> const &rts,
	utility::vector1<core::kinematics::RT> const &rts_all,
	bool verbose=false
) {
	numeric::xyzMatrix<Real> identity = numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);

	utility::vector1<core::kinematics::RT> rts_exp = rts;
	int nxformsOrig = (int)rts_exp.size();

	if ( verbose ) {
		for ( int i=1; i<=(int)rts_exp.size(); ++i ) {
			disp_symmop(rts_exp[i]);
		}
	}

	bool done = false;
	int rounds=0;
	int startXform = 1;
	while ( !done && ++rounds<16 ) {
		int nxforms = (int)rts_exp.size();
		done = true;
		for ( int i=startXform; i<=nxforms; ++i ) {
			for ( int j=1; j<=nxformsOrig; ++j ) {
				numeric::xyzMatrix<Real> const &Si=rts_exp[i].get_rotation();
				numeric::xyzVector<Real> const &Ti=rts_exp[i].get_translation();
				numeric::xyzMatrix<Real> const &Sj=rts_exp[j].get_rotation();
				numeric::xyzVector<Real> const &Tj=rts_exp[j].get_translation();

				numeric::xyzMatrix<Real> Sij = Sj*Si;
				numeric::xyzVector<Real> Tij = Tj+Sj*Ti;

				//if ( std::fabs(Tij[0])>2 || std::fabs(Tij[1])>2 || std::fabs(Tij[2])>2 ) continue;
				if ( Tij[0]<-1 || Tij[1]<-1 || Tij[2]<-1 ) continue;
				if ( Tij[0]>2 || Tij[1]>2 || Tij[2]>2 ) continue;

				// check if unique
				bool unique = true;
				for ( int k=1; k<=(int)rts_exp.size() && unique; ++k ) {
					if ( transforms_equiv( Sij, Tij, rts_exp[k].get_rotation(), rts_exp[k].get_translation() ) ) {
						unique = false;
					}
				}

				if ( unique ) {
					done = false;
					rts_exp.push_back( core::kinematics::RT( Sij, Tij ) );
				}
			}
		}
		startXform = nxforms+1; // we can start where we left off

		// check this at each round to see if we can exit early. Check if:
		//    1 - all generated transformations are in the lattice group
		//    2 - all lattice transforms are generated
		bool connected = (rts_exp.size() > rts_all.size());
		if ( connected ) {
			for ( int i=1; i<=(int)rts_all.size() && connected; ++i ) {
				bool contains_j = false;
				for ( int j=1; j<=(int)rts_exp.size(); ++j ) {
					bool i_equals_j = transforms_equiv_mod1(
						rts_exp[j].get_rotation(), rts_exp[j].get_translation(),
						rts_all[i].get_rotation(), rts_all[i].get_translation());
					contains_j |= i_equals_j;
				}
				connected &= contains_j;
			}
			if ( connected ) {
				bool has_001=false, has_010=false, has_100=false;
				for ( int i=1; i<=(int)rts_exp.size(); ++i ) {
					numeric::xyzMatrix<Real> const &Si=rts_exp[i].get_rotation();
					numeric::xyzVector<Real> const &Ti=rts_exp[i].get_translation();
					if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,0,1)) ) has_001 = true;
					if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(0,1,0)) ) has_010 = true;
					if ( transforms_equiv(Si, Ti, identity, numeric::xyzVector<Real>(1,0,0)) ) has_100 = true;
				}
				connected = (has_001 && has_010 && has_100);
			}
		}

		if ( connected ) return true;
	}

	return false;
}

//void get_symmops(std::string name_, utility::vector1<core::kinematics::RT> &rt_out, utility::vector1<numeric::xyzVector<core::Real> > &cenop);
void
set_all_spacegroups( utility::vector1<std::string> &spacegroups ) {
	spacegroups = utility::vector1<std::string>({
		"P1", "P-1", "P121", "P1211", "C121", "P1m1", "P1c1", "C1m1", "C1c1", "P12/m1", "P121/m1", "C12/m1", "P12/c1", "P121/c1", "C12/c1",
		"P222", "P2221", "P21212", "P212121", "C2221", "C222", "F222", "I222", "I212121", "Pmm2", "Pmc21", "Pcc2", "Pma2", "Pca21", "Pnc2",
		"Pmn21", "Pba2", "Pna21", "Pnn2", "Cmm2", "Cmc21", "Ccc2", "Amm2", "Abm2", "Ama2", "Aba2", "Fmm2", "Fdd2", "Imm2", "Iba2",
		"Ima2", "Pmmm", "Pnnn:2", "Pccm", "Pban:2", "Pmma", "Pnna", "Pmna", "Pcca", "Pbam", "Pccn", "Pbcm", "Pnnm", "Pmmn:2", "Pbcn",
		"Pbca", "Pnma", "Cmcm", "Cmca", "Cmmm", "Cccm", "Cmma", "Ccca:2", "Fmmm", "Fddd:2", "Immm", "Ibam", "Ibca", "Imma",
		"P4", "P41", "P42", "P43", "I4", "I41", "P-4", "I-4", "P4/m", "P42/m", "P4/n:2", "P42/n:2", "I4/m", "I41/a:2", "P422", "P4212",
		"P4122", "P41212", "P4222", "P42212", "P4322", "P43212", "I422", "I4122", "P4mm", "P4bm", "P42cm", "P42nm", "P4cc", "P4nc", "P42mc",
		"P42bc", "I4mm", "I4cm", "I41md", "I41cd", "P-42m", "P-42c", "P-421m", "P-421c", "P-4m2", "P-4c2", "P-4b2", "P-4n2", "I-4m2", "I-4c2",
		"I-42m", "I-42d", "P4/mmm", "P4/mcc", "P4/nbm:2", "P4/nnc:2", "P4/mbm", "P4/mnc", "P4/nmm:2", "P4/ncc:2", "P42/mmc", "P42/mcm",
		"P42/nbc:2", "P42/nnm:2", "P42/mbc", "P42/mnm", "P42/nmc:2", "P42/ncm:2", "I4/mmm", "I4/mcm", "I41/amd:2", "I41/acd:2",
		"P3", "P31", "P32", "R3:H", "P-3", "R-3:H", "P312", "P321", "P3112", "P3121", "P3212", "P3221", "R32:H", "P3m1", "P31m", "P3c1",
		"P31c", "R3m:H", "R3c:H", "P-31m", "P-31c", "P-3m1", "P-3c1", "R-3m:H", "R-3c:H", "P6", "P61", "P65", "P62", "P64", "P63", "P-6",
		"P6/m", "P63/m", "P622", "P6122", "P6522", "P6222", "P6422", "P6322", "P6mm", "P6cc", "P63cm", "P63mc", "P-6m2", "P-6c2", "P-62m",
		"P-62c", "P6/mmm", "P6/mcc", "P63/mcm", "P63/mmc", "P23", "F23", "I23", "P213", "I213", "Pm-3", "Pn-3:2", "Fm-3", "Fd-3:2", "Im-3", "Pa-3",
		"Ia-3", "P432", "P4232", "F432", "F4132", "I432", "P4332", "P4132", "I4132", "P-43m", "F-43m", "I-43m", "P-43n", "F-43c", "I-43d",
		"Pm-3m", "Pn-3n:2", "Pm-3n", "Pn-3m:2", "Fm-3m", "Fm-3c", "Fd-3m:2", "Fd-3c:2", "Im-3m", "Ia-3d"
		});
}

int
main( int argc, char * argv [] ) {
	try {
		// all spacegroups
		utility::vector1<std::string> spacegroups;

		std::string mode;
		if ( argc<2 ) {
			mode = "pt";
		} else {
			mode = argv[1];
		}

		if ( mode!="pt" && mode!="gen1" && mode!="gen2" ) {
			std::cerr << "USAGE: " << argv[0] << " (pt|gen1|gen2) [spacegroup]*" << std::endl;
			exit(0);
		}

		bool user_spacegroups=false;
		if ( argc>2 ) {
			for ( int i=2; i<argc; ++i ) {
				char chk=argv[i][0];
				if ( chk=='P' || chk == 'F' || chk=='I' || chk=='H' || chk=='C' ) user_spacegroups = true;
			}
		}

		if ( !user_spacegroups ) {
			set_all_spacegroups( spacegroups );
		} else {
			for ( int i=2; i<argc; ++i ) {
				char chk=argv[i][0];
				if ( chk=='P' || chk == 'F' || chk=='I' || chk=='H' || chk=='C' ) {
					spacegroups.push_back(std::string(argv[i]));
				}
			}
		}

		for ( int sx=1; sx<=(int) spacegroups.size(); ++sx ) {
			std::string name=spacegroups[sx];

			numeric::xyzMatrix<core::Real> skewM=numeric::xyzMatrix<Real>::rows(  1,0,0,   0,1,0,   0,0,1);
			if (
					name == "P3" || name == "P31" || name == "P32" || name == "R3:H" || name == "P-3" || name == "R-3:H"
					|| name == "P312" || name == "P321" || name == "P3112" || name == "P3121" || name == "P3212" || name == "P3221"
					|| name == "R32:H" || name == "P3m1" || name == "P31m" || name == "P3c1" || name == "P31c" || name == "R3m:H"
					|| name == "R3c:H" || name == "P-31m" || name == "P-31c" || name == "P-3m1" || name == "P-3c1"
					|| name == "R-3m:H" || name == "R-3c:H" || name == "P6" || name == "P61" || name == "P65" || name == "P62"
					|| name == "P64" || name == "P63" || name == "P-6" || name == "P6/m" || name == "P63/m" || name == "P622"
					|| name == "P6122" || name == "P6522" || name == "P6222" || name == "P6422" || name == "P6322"
					|| name == "P6mm" || name == "P6cc" || name == "P63cm" || name == "P63mc" || name == "P-6m2"
					|| name == "P-6c2" || name == "P-62m" || name == "P-62c" || name == "P6/mmm" || name == "P6/mcc"
					|| name == "P63/mcm" || name == "P63/mmc"
					) {
				skewM = numeric::xyzMatrix<core::Real>::rows( 1,-0.5,0, 0,sqrt(3)/2,0, 0,0,1);
			}

			protocols::cryst::Spacegroup sg;
			sg.set_spacegroup( spacegroups[sx] );
			utility::vector1<core::kinematics::RT> rts = sg.symmops();
			utility::vector1< pointGroupHit > pgs = get_point_groups( rts, skewM );

			if ( mode == "pt" ) {
				for ( int i=1; i<=(int)pgs.size(); ++i ) {
					pointGroupHit hit_i = pgs[i];
					std::cout << spacegroups[sx] << " ";
					hit_i.show(std::cout);
				}
			} else {
				utility::vector1< latticeHit > finalHits;
				std::cout << "Found " << pgs.size() << " point groups" << std::endl;
				if ( mode == "gen2" ) {
					// sample all pairs of point groups
					for ( int i=1; i<=(int)pgs.size(); ++i ) {
						for ( int j=i+1; j<=(int)pgs.size(); ++j ) {
							pointGroupHit hit_i = pgs[i];
							pointGroupHit hit_j = pgs[j];

							// get RTs
							utility::vector1<core::kinematics::RT> rts_symm;
							rts_symm.reserve( hit_i.order() + hit_j.order() );
							hit_i.append_symmops( rts_symm, rts );
							hit_j.append_symmops( rts_symm, rts );

							// expand
							bool is_lattice = check_if_forms_lattice( rts_symm, rts );   // check if we are fully connected

							if ( is_lattice ) {
								core::Real angle,offset,shift;
								get_angle_and_offset( hit_i, hit_j, angle,offset,shift );

								bool dup = false;
								for ( int q = 1; q<= (int)finalHits.size() && !dup; ++q ) {
									bool name_match = (
										(finalHits[q].hit1.name() == hit_i.name() && finalHits[q].hit2.name() == hit_j.name()) ||
										(finalHits[q].hit2.name() == hit_i.name() && finalHits[q].hit1.name() == hit_j.name())
									);
									bool angle_match = std::fabs( angle - finalHits[q].angle ) < 1e-6;
									if ( name_match && angle_match ) {
										if ( (offset==0 && finalHits[q].offset == 0) || (offset!=0 && finalHits[q].offset != 0) ) {
											dup = true;
											if ( offset<finalHits[q].offset ) {
												finalHits[q] = latticeHit( hit_i, hit_j, angle, offset, shift );
											}
										}
									}
								}
								if ( !dup ) {
									finalHits.push_back( latticeHit( hit_i, hit_j, angle, offset, shift ));
								}
							}
						}
					}
				}

				for ( int i=1; i<=(int)finalHits.size(); ++i ) {
					latticeHit lh = finalHits[i];

					// one line
					std::cout << name << " " << lh.hit1.name() << " " << lh.hit2.name() << " " << lh.angle << " " << lh.offset << " ";
					if ( lh.angle == 0 ) {
						std::cout << lh.shift << " ";
					} else {
						std::cout << "0 ";
					}
					std::cout << lh.hit1.as_string() << " ";
					std::cout << lh.hit2.as_string() << std::endl;
				}
			}
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;
}

