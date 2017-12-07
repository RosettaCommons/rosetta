/// @file
/// @brief

#ifndef INCLUDED_protocols_cryst_spacegroup_hh
#define INCLUDED_protocols_cryst_spacegroup_hh

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/RT.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

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

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

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

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>


#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

static basic::Tracer TSG( "spacegroup" );

core::Real absfpart( core::Real x) {
	return ( std::fabs( x-std::floor( x+0.5) ) );
}

core::Size denom( core::Real x ) {
	if ( absfpart(x) <= 1e-4 ) { return 1; }
	for ( core::Size i=2; i<=6; ++i ) {
		if ( absfpart(i*x) <= 1e-4 ) { return i;}
	}

	utility_exit_with_message( "error in denom()");
	return 0;
}

enum SpacegroupSetting {
	TRICLINIC,
	MONOCLINIC,
	ORTHORHOMBIC,
	TETRAGONAL,
	HEXAGONAL,  // includes trigonal
	CUBIC,
	NOSETTING
};

struct CheshireCell {
	CheshireCell() : low(0.0,0.0,0.0),high(0.0,0.0,0.0) {};
	CheshireCell( numeric::xyzVector<core::Real> low_in, numeric::xyzVector<core::Real> high_in ) {
		low=low_in;
		high=high_in;
	}
	numeric::xyzVector<core::Real> low, high;
};

class Spacegroup {
private:
	std::string name_;
	SpacegroupSetting setting_;

	// cryst stuff
	numeric::xyzMatrix<core::Real> f2c_, c2f_;
	core::Real a_, b_, c_, alpha_, beta_, gamma_, V_;
	core::Size ncopies_;

	utility::vector1<core::kinematics::RT> symmops_;
	CheshireCell cc_;

public:
	Spacegroup();
	Spacegroup(std::string name_in);

	void
	set_spacegroup( std::string name_in);

	numeric::xyzMatrix<core::Real> const &f2c() const { return f2c_; };
	numeric::xyzMatrix<core::Real> const &c2f() const { return c2f_; };
	core::Real A() const { return a_; }
	core::Real B() const { return b_; }
	core::Real C() const { return c_; }
	core::Real alpha() const { return alpha_; }
	core::Real beta() const { return beta_; }
	core::Real gamma() const { return gamma_; }
	core::Real volume() const { return V_; }

	SpacegroupSetting setting() const { return setting_; }
	utility::vector1<core::kinematics::RT> const &symmops() const { return symmops_; }
	core::kinematics::RT symmop(core::Size i) const { return symmops_[i]; }
	core::Size nsymmops() const { return symmops_.size(); }
	CheshireCell cheshire_cell() const { return cc_; }

	// grid spacing must be a multiple of this number
	core::Size minmult() const {
		if ( setting_ == TRICLINIC )   return 2;
		if ( setting_ == MONOCLINIC )   return 4;
		if ( setting_ == CUBIC )        return 4;
		if ( setting_ == ORTHORHOMBIC ) return 4;
		if ( setting_ == TETRAGONAL )   return 8;
		if ( setting_ == HEXAGONAL )    return 6;
		return 0;
	}

	// sets AND VALIDATES input parameters
	void set_parameters(core::Real a_in, core::Real b_in, core::Real c_in, core::Real alpha_in, core::Real beta_in, core::Real gamma_in);

	std::string
	get_moveable_dofs() const {
		utility::vector1<core::kinematics::RT> rt;
		CheshireCell  cc;
		get_symmops (rt,cc);
		std::string retval;
		if ( cc.high[0]>cc.low[0] ) retval += "x ";
		if ( cc.high[1]>cc.low[1] ) retval += "y ";
		if ( cc.high[2]>cc.low[2] ) retval += "z ";
		return retval;
	}

	numeric::xyzVector<core::Size>
	get_nsubdivisions() {
		numeric::xyzVector<core::Size> retval(1,1,1);
		for ( int i=1; i<=(int)nsymmops(); ++i ) {
			numeric::xyzVector<core::Real> const &T = symmops_[i].get_translation();
			retval[0] = std::max( retval[0], denom( T[0] ) );
			retval[1] = std::max( retval[1], denom( T[1] ) );
			retval[2] = std::max( retval[2], denom( T[2] ) );
		}
		return retval;
	}

	numeric::xyzVector<core::Size>
	get_trans_dofs() {
		if ( setting_ == TRICLINIC ) {
			return numeric::xyzVector<core::Size>(1,2,3);
		} else if ( setting_ == MONOCLINIC ) {
			return numeric::xyzVector<core::Size>(1,2,3);
		} else if ( setting_ == CUBIC ) {
			return numeric::xyzVector<core::Size>(1,1,1);
		} else if ( setting_ == ORTHORHOMBIC ) {
			return numeric::xyzVector<core::Size>(1,2,3);
		} else if ( setting_ == TETRAGONAL ) {
			return numeric::xyzVector<core::Size>(1,1,3);
		} else if ( setting_ == HEXAGONAL ) {
			return numeric::xyzVector<core::Size>(1,1,3);
		}

		return numeric::xyzVector<core::Size>(1,2,3);  // no warnings
	}

	core::Size
	get_nrot_dofs() {
		if ( setting_ == TRICLINIC ) {
			return 3;
		} else if ( setting_ == MONOCLINIC ) {
			return 1;
		}
		return 0;
	}

	// get 'CRYST1' name of spacegroup
	std::string pdbname() const;

private:
	// get symmops
	void get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const;
};

///
Spacegroup::Spacegroup() {
	name_="";
	setting_=NOSETTING;
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;
}

Spacegroup::Spacegroup(std::string name_in) {
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;
	set_spacegroup( name_in );
}

void
Spacegroup::set_spacegroup( std::string name_in) {
	// reset lattice
	a_=b_=c_=0;
	alpha_=beta_=gamma_=90.0;
	V_=0;
	ncopies_=0;

	name_ = name_in;
	name_.erase( std::remove( name_.begin(), name_.end(), ' ' ), name_.end() );

	// setting + validation
	if ( name_ == "P1" ) {
		setting_ = TRICLINIC;
	} else if ( (name_ == "P2") ||  (name_ == "P21") ||  (name_ == "C2") || (name_ == "P121") ||  (name_ == "P1211") ||  (name_ == "C121") ) {
		setting_ = MONOCLINIC;
	} else if ( (name_ == "P23") ||  (name_ == "F23") ||  (name_ == "I23") ||
			(name_ == "P213") ||  (name_ == "I213") ||  (name_ == "P432") ||
			(name_ == "P4232") ||  (name_ == "F432") ||  (name_ == "F4132") ||
			(name_ == "I432") ||  (name_ == "P4332") ||  (name_ == "P4132") ||
			(name_ == "I4132") ) {
		setting_ = CUBIC;
	} else if ( (name_ == "P222") ||  (name_ == "P2221") ||  (name_ == "P21212") ||
			(name_ == "P212121") ||  (name_ == "C2221") ||  (name_ == "C222") ||
			(name_ == "F222") ||  (name_ == "I222") ||  (name_ == "I212121") ) {
		setting_ = ORTHORHOMBIC;
	} else if ( (name_ == "P4") ||  (name_ == "P41") ||  (name_ == "P42") ||
			(name_ == "P43") ||  (name_ == "I4") ||  (name_ == "I41") ||
			(name_ == "P422") ||  (name_ == "P4212") ||  (name_ == "P4122") ||
			(name_ == "P41212") ||  (name_ == "P4222") ||  (name_ == "P42212") ||
			(name_ == "P4322") ||  (name_ == "P43212") ||  (name_ == "I422") ||
			(name_ == "I4122") ) {
		setting_ = TETRAGONAL;
	} else if ( (name_ == "P3") ||  (name_ == "P31") ||  (name_ == "P32") ||
			(name_ == "H3") ||  (name_ == "P312") ||  (name_ == "P321") ||
			(name_ == "P3112") ||  (name_ == "P3121") ||  (name_ == "P3212") ||
			(name_ == "P3221") ||  (name_ == "H32") ||  (name_ == "P6") ||
			(name_ == "P61") ||  (name_ == "P65") ||  (name_ == "P62") ||
			(name_ == "P64") ||  (name_ == "P63") ||  (name_ == "P622") ||
			(name_ == "P6122") ||  (name_ == "P6522") ||  (name_ == "P6222") ||
			(name_ == "P6422") ||  (name_ == "P6322") ) {
		setting_ = HEXAGONAL;
	} else {
		utility_exit_with_message("Unknown spacegroup! "+name_);
	}

	// lookup
	symmops_.clear();
	get_symmops(symmops_, cc_);
}


// sets AND VALIDATES input parameters
void Spacegroup::set_parameters(core::Real a_in, core::Real b_in, core::Real c_in, core::Real alpha_in, core::Real beta_in, core::Real gamma_in) {
	if ( setting_ == TRICLINIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=alpha_in; beta_=beta_in; gamma_=gamma_in;
	} else if ( setting_ == MONOCLINIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=beta_in; gamma_=90.0;
	} else if ( setting_ == CUBIC ) {
		a_=a_in; b_=a_in; c_=a_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == ORTHORHOMBIC ) {
		a_=a_in; b_=b_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == TETRAGONAL ) {
		a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=90.0;
	} else if ( setting_ == HEXAGONAL ) {
		a_=a_in; b_=a_in; c_=c_in; alpha_=90.0; beta_=90.0; gamma_=120.0;
	}

	// transformation matrices
	core::Real ca = cos(DEG2RAD*alpha_), cb = cos(DEG2RAD*beta_), cg = cos(DEG2RAD*gamma_);
	core::Real /*sa = sin(DEG2RAD*alpha_),*/ sb = sin(DEG2RAD*beta_), sg = sin(DEG2RAD*gamma_);
	f2c_ = numeric::xyzMatrix<core::Real>::rows(
		a_  , b_ * cg , c_ * cb,
		0.0 , b_ * sg , c_ * (ca - cb*cg) / sg,
		0.0 , 0.0     , c_ * sb * sqrt(1.0 - numeric::square((cb*cg - ca)/(sb*sg)))
	);
	c2f_ = numeric::inverse(f2c_);
	V_ = a_*b_*c_* sqrt(1-numeric::square(ca)-numeric::square(cb)-numeric::square(cg)+2*ca*cb*cg);

	// report
	if ( a_!=a_in || b_!=b_in || c_!=c_in || alpha_!=alpha_in || beta_!=beta_in || gamma_!=gamma_in ) {
		TSG << "Overriding input crystal parameters with [ "
			<< a_ << "," << b_ << "," << c_ << " , " << alpha_ << ","  << beta_ << ","  << gamma_ << " ]" << std::endl;
	}
}

///
void Spacegroup::get_symmops(utility::vector1<core::kinematics::RT> &rt_out, CheshireCell &cc) const {
	if ( name_ == "P1" ) {
		rt_out.resize(1);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0 ,0 ,0 ) );
	}
	if ( name_ == "P121" || name_ == "P2" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0,0,0),numeric::xyzVector<core::Real>(0.5,0,0.5) );
	} else if ( name_ == "P1211" || name_ == "P21" ) {
		rt_out.resize(2);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0 ,0.5 ) );
	} else if ( name_ == "C121" || name_ == "C2" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0 ,0.5 ) );
	} else if ( name_ == "P4" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "P41" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "P42" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "P43" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "I4" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "I41" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0 ) );
	} else if ( name_ == "P422" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P4212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P4122" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P41212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P4222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P42212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P4322" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.75) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P43212" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "I422" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=8; ++ii ) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "I4122" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,-1,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   -1,0,0,   0,0,1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0.5,0.25) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,1,0,   1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   -1,0,0,   0,1,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0.5,0,0.75) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   0,-1,0,   -1,0,0,   0,0,-1 )  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=8; ++ii ) {
			rt_out[8+ii] = rt_out[ii];
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,0.5 ,0.5 ) );
	} else if ( name_ == "P23" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0,0,0),numeric::xyzVector<core::Real>(1,1,1) );
	} else if ( name_ == "F23" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[24+ii] = rt_out[ii];
			rt_out[36+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[36+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(0.5,0.5,0.5 ) );
	} else if ( name_ == "I23" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P213" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "I213" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		for ( int ii=1; ii<=12; ++ii ) {
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P432" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P4232" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "F432" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(0.5,0.5,0.5 ) );
	} else if ( name_ == "F4132" ) {
		rt_out.resize(96);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		// cenops
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[48+ii] = rt_out[ii];
			rt_out[72+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[48+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[72+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(0.5,0.5,0.5 ) );
	} else if ( name_ == "I432" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P4332" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P4132" ) {
		rt_out.resize(24);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.75) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "I4132" ) {
		rt_out.resize(48);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,1,   0,1,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.25) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   -1,0,0,   0,1,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,-1,   0,1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.75) );
		rt_out[13] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   -1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[14] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,0,1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.25,0.75,0.75) );
		rt_out[15] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   1,0,0,   0,-1,0)  , numeric::xyzVector<core::Real>(0.5,0,0) );
		rt_out[16] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,0,-1,   0,-1,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[17] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,1,   1,0,0)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[18] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   0,0,-1,   -1,0,0)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[19] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.25) );
		rt_out[20] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,1,   -1,0,0)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		rt_out[21] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,-1,0,   -1,0,0)  , numeric::xyzVector<core::Real>(0.25,0.25,0.25) );
		rt_out[22] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   0,0,-1,   1,0,0)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[23] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,1,   0,-1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.75,0.25) );
		rt_out[24] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,0,-1,   0,1,0,   1,0,0)  , numeric::xyzVector<core::Real>(0.75,0.25,0.75) );
		// cenops
		for ( int ii=1; ii<=24; ++ii ) {
			rt_out[24+ii] = rt_out[ii];
			rt_out[24+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0,0,0),numeric::xyzVector<core::Real>(1,1,1 ) );
	} else if ( name_ == "P3" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	} else if ( name_ == "P31" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	} else if ( name_ == "P32" ) {
		rt_out.resize(3);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	} else if ( name_ == "H3" ) {
		rt_out.resize(9);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=3; ++ii ) {
			rt_out[3+ii] = rt_out[ii];
			rt_out[3+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0 ) );
	} else if ( name_ == "P312" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	} else if ( name_ == "P321" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P3112" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	} else if ( name_ == "P3121" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P3212" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(2.0/3.0 ,2.0/3.0 ,0.5 ) );
	} else if ( name_ == "P3221" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "H32" ) {
		rt_out.resize(18);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=6; ++ii ) {
			rt_out[6+ii] = rt_out[ii];
			rt_out[6+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.666666666666667,0.333333333333333,0.333333333333333) );
			rt_out[12+ii] = rt_out[ii];
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.333333333333333,0.666666666666667,0.666666666666667) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P61" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P65" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P62" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P64" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P63" ) {
		rt_out.resize(6);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0 ) );
	} else if ( name_ == "P622" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6122" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6522" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.166666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.833333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6222" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6422" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.333333333333333) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.666666666666667) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P6322" ) {
		rt_out.resize(12);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   0,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   1,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[5] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   -1,0,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[6] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   -1,1,0,   0,0,1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[7] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,-1,0,   -1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[8] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,-1,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[9] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  1,0,0,   1,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[10] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  0,1,0,   1,0,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[11] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,1,0,   0,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[12] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   -1,1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>( 0, 0, 0),numeric::xyzVector<core::Real>(1 ,1 ,0.5 ) );
	} else if ( name_ == "P222" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P2221" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P21212" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "P212121" ) {
		rt_out.resize(4);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0.5,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0.5,0.5,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "C2221" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "C222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "F222" ) {
		rt_out.resize(16);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[8+ii] = rt_out[ii];
			rt_out[12+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0,0.5,0.5) );
			rt_out[8+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0,0.5) );
			rt_out[12+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "I222" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	} else if ( name_ == "I212121" ) {
		rt_out.resize(8);
		rt_out[1] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0, 1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0,0) );
		rt_out[2] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0,-1,0,   0,0, 1)  , numeric::xyzVector<core::Real>(0,0.5,0) );
		rt_out[3] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(   1,0,0,   0,-1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0,0.5) );
		rt_out[4] = core::kinematics::RT( numeric::xyzMatrix<core::Real>::rows(  -1,0,0,   0, 1,0,   0,0,-1)  , numeric::xyzVector<core::Real>(0,0.5,0.5) );
		// cenops
		for ( int ii=1; ii<=4; ++ii ) {
			rt_out[4+ii] = rt_out[ii];
			rt_out[4+ii].set_translation( rt_out[ii].get_translation() + numeric::xyzVector<core::Real>(0.5,0.5,0.5) );
		}
		cc = CheshireCell( numeric::xyzVector<core::Real>(0, 0, 0),numeric::xyzVector<core::Real>(0.5 ,0.5 ,0.5) );
	}

	// cenops can make translation outside of [0,1], correct this
	for ( int ii=1; ii<=(int)rt_out.size(); ++ii ) {
		core::Vector T_i = rt_out[ii].get_translation();
		rt_out[ii].set_translation( core::Vector(std::fmod(T_i[0],1), std::fmod(T_i[1],1), std::fmod(T_i[2],1) ) );
	}
}

// name output in pdbheader
std::string Spacegroup::pdbname() const {
	if ( name_ == "P1" ) return "P 1";
	if ( name_ == "P121" || name_ == "P2" ) return "P 1 2 1";
	if ( name_ == "P1211" || name_ == "P21" ) return "P 1 21 1";
	if ( name_ == "C121" || name_ == "C2" ) return "C 1 2 1";
	if ( name_ == "P4" ) return "P 4";
	if ( name_ == "P41" ) return "P 41";
	if ( name_ == "P42" ) return "P 42";
	if ( name_ == "P43" ) return "P 43";
	if ( name_ == "I4" ) return "I 4";
	if ( name_ == "I41" ) return "I 41";
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
	if ( name_ == "P3" ) return "P 3";
	if ( name_ == "P31" ) return "P 31";
	if ( name_ == "P32" ) return "P 32";
	if ( name_ == "H3" ) return "H 3";
	if ( name_ == "P312" ) return "P 3 1 2";
	if ( name_ == "P321" ) return "P 3 2 1";
	if ( name_ == "P3112" ) return "P 31 1 2";
	if ( name_ == "P3121" ) return "P 31 2 1";
	if ( name_ == "P3212" ) return "P 32 1 2";
	if ( name_ == "P3221" ) return "P 32 2 1";
	if ( name_ == "H32" ) return "H 3 2";
	if ( name_ == "P6" ) return "P 6";
	if ( name_ == "P61" ) return "P 61";
	if ( name_ == "P65" ) return "P 65";
	if ( name_ == "P62" ) return "P 62";
	if ( name_ == "P64" ) return "P 64";
	if ( name_ == "P63" ) return "P 63";
	if ( name_ == "P622" ) return "P 6 2 2";
	if ( name_ == "P6122" ) return "P 61 2 2";
	if ( name_ == "P6522" ) return "P 65 2 2";
	if ( name_ == "P6222" ) return "P 62 2 2";
	if ( name_ == "P6422" ) return "P 64 2 2";
	if ( name_ == "P6322" ) return "P 63 2 2";
	if ( name_ == "P222" ) return "P 2 2 2";
	if ( name_ == "P2221" ) return "P 2 2 21";
	if ( name_ == "P21212" ) return "P 21 21 2";
	if ( name_ == "P212121" ) return "P 21 21 21";
	if ( name_ == "C2221" ) return "C 2 2 21";
	if ( name_ == "C222" ) return "C 2 2 2";
	if ( name_ == "F222" ) return "F 2 2 2";
	if ( name_ == "I222" ) return "I 2 2 2";
	if ( name_ == "I212121" ) return "I 21 21 21";
	return "X";
}

#endif
