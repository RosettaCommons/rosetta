/// @file
/// @brief

#ifndef INCLUDED_protocols_cryst_spacegroup_hh
#define INCLUDED_protocols_cryst_spacegroup_hh

#include <protocols/cryst/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
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



namespace protocols {
namespace cryst {

core::Real absfpart( core::Real x);

core::Size denom( core::Real x );

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

	std::string
	name() const { return name_; }

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


}
}

#endif
