// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <basic/Tracer.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/abinitio/FragmentMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>


using numeric::conversions::radians;

static basic::Tracer TR("dubois");

using core::Size;
using core::Real;
using core::id::AtomID;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;
using numeric::random::uniform;


Vec helix_axis(core::pose::Pose const & pose) {
	Vec axis(0,0,0);
	for(Size i = 2; i <= pose.n_residue()-4; ++i) {
		axis += ( pose.residue(i+4).xyz(1) - pose.residue(i).xyz(3) );
	}
	axis.normalize();
	return axis;
}

inline Vec center_of_mass( core::pose::Pose const & pose, Size nres = 0 ) {
	if( 0 == nres ) nres = pose.n_residue();
	Vec com(0,0,0);
	Size count = 0;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			com += pose.xyz(core::id::AtomID(ia,ir));
			count++;
		}
	}
	return com/(Real)count;
}

inline void trans_pose( core::pose::Pose & pose, Vec const & trans ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + trans );
		}
	}
}

inline void rot_pose( core::pose::Pose & pose, Mat const & rot ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}

inline void rot_pose( core::pose::Pose & pose, Mat const & rot, Vec const & cen ) {
	trans_pose(pose,-cen);
	rot_pose(pose,rot);
	trans_pose(pose,cen);
}

inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

inline void rot_pose( core::pose::Pose & pose, Vec const & axis, Real const & ang, Vec const & cen ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
}

core::pose::Pose make_helix(std::string seq) {
	using namespace core::chemical;
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose,"input/start_HC.pdb");
	core::pose::remove_lower_terminus_type_from_pose_residue(pose,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(pose,1);
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );
	Size nres = seq.size();
	for( Size i = 2; i <= nres; ++i ){
		std::string name3 = name_from_aa(aa_from_oneletter_code(seq[i-1]));
		pose.append_residue_by_bond(*core::conformation::ResidueFactory::create_residue(residue_set->name_map(name3)),true);
	}
	core::pose::add_upper_terminus_type_to_pose_residue(pose,pose.n_residue());
	for( Size i = 1; i <= seq.size(); i++ ) {
		pose.set_phi  (i,-60.16731);
		pose.set_psi  (i,-45.19451);
		pose.set_omega(i,180.00000);
	}

	trans_pose(pose,-center_of_mass(pose));
	Vec axis = helix_axis(pose);
	Vec z = Vec(0,0,1);
	Vec rot_axis = z.cross(axis);
	Real rot_ang = acos(z.dot(axis));
	Mat rot = rotation_matrix_radians( rot_axis, -rot_ang );
	std::cerr << helix_axis(pose) << std::endl;
	rot_pose(pose,rot);
	std::cerr << helix_axis(pose) << std::endl;
	return pose;
}

core::conformation::symmetry::SymmData
make_symm_data(
	core::pose::Pose const & pose,
	core::Real rot,
	core::Real trans,
	core::Size n
) {
	using namespace ObjexxFCL;
	Size anchor = pose.n_residue() / 2;
	std::string s = "";
	s += "symmetry_name c12345\nsubunits "+string_of(n)+"\nnumber_of_interfaces "+string_of(n-1)+"\n";
	s += "E = 1.0*VRT0001";
	for(Size i = 2; i<=n; i++) s+= " + 1*(VRT1:VRT000"+string_of(i)+")";
	s += "\nanchor_residue " + string_of(anchor) + "\n";
	s += "virtual_transforms_start consecutive\nstart -1,0,0 0,1,0 0,0,0\n";
	s += "rot Rz_angle " + string_of(rot) + "\n";
	s += "trans 0,0," + string_of(trans);
	s += "\nvirtual_transforms_stop\n";
	for(Size i = 2; i<=n; i++) s+= "connect_virtual J"+string_of(i)+" VRT000"+string_of(i-1)+" VRT000"+string_of(i)+"\n";
	s += "set_dof BASEJUMP x angle_x angle_y angle_z\n";
	s += "set_dof J2 z angle_z\n";
	// TR << "================= symm dat ==================" << std::endl;
	// TR << s << std::endl;
	// TR << "================= symm dat ==================" << std::endl;
	std::istringstream iss(s);
	core::conformation::symmetry::SymmData symdat( pose.n_residue(), pose.num_jump() );
	symdat.read_symmetry_data_from_stream(iss);
	return symdat;
}

struct CCParam {
	CCParam() {
		randomize();
	}
	void randomize() {
		//nres  = 24;
		//nsub  = 6;
		//rot   = -179.8721277314733;
		//trans = 9.975812165562381;
		//x     = 4.441085794869124;
		//rh    = 269.1484789308233;
		//rhx   = -29.57581636464426;
		//rhy   = -6.035161566019054;
		using namespace numeric::random;
		nres = 15.0 + 15.0*(uniform());
		rot = 120.0 + 30.0*gaussian();
		x = 6.0 + 2.0*gaussian();
		rh = 360.0*uniform();
		rhx = 7.0*gaussian()+20.0*32.0/(Real)nres;
		if(uniform()<0.5) rot *= -1.0;
		if(uniform()<0.5) rhx *= -1.0;
		rhy =  6.0*gaussian();
		trans = 9.0 + 3*uniform();
		nsub = 9;
	}
	Size nres,nsub;
	Real rot, trans, x, rh, rhx, rhy;
	void show() {
		TR << "nres "  << nres << " ";
		TR << "nsub "	 << nsub << " ";
		TR << "rot "   << rot << " ";
		TR << "trans " << trans << " ";
		TR << "x "     << x << " ";
		TR << "rh "    << rh << " ";
		TR << "rhx "   << rhx << " ";
		TR << "rhy "   << rhy;
	}
	std::string str() {
		return string_of(nres)+" "+string_of(nsub)+" "+string_of(rot)+" "+string_of(trans)+" "+string_of(x)+" "
		         +string_of(rh)+" "+string_of(rhx)+" "+string_of(rhy);
	}
	utility::vector1<Size> cys_pos;
};

core::pose::Pose make_coiled_coil(CCParam & p) {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace pose;
	using namespace conformation::symmetry;

	std::string seq = "";
	for(Size i = 1; i <= p.nres; ++i) seq += "A";
	for(Size i = 1; i <= p.cys_pos.size(); ++i) seq[p.cys_pos[i]-1] = 'C';
	Pose helix = make_helix(seq);
	rot_pose(  helix,Vec(0,0,1), p.rh );
	rot_pose(  helix,Vec(1,0,0), p.rhx);
	rot_pose(  helix,Vec(0,1,0), p.rhy);
	trans_pose(helix,Vec(p.x,0,0));

	Real mnz=99999,mxz=-99999;
	for(Size ir = 1; ir <= helix.n_residue(); ++ir) {
		for(Size ia = 1; ia <= helix.residue_type(ir).natoms(); ++ia) {
			Real z = helix.xyz(core::id::AtomID(ia,ir)).z();
			if(z<mnz) mnz = z;
			if(z>mxz) mxz = z;
		}
	}
	p.nsub = std::ceil(( mxz-mnz+10 ) / p.trans) + 1;

	Pose cc = helix;
	SymmData symmdata( make_symm_data(cc,p.rot,p.trans,p.nsub) );
	core::pose::symmetry::make_symmetric_pose( cc, symmdata );

	return cc;
}



/////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace core;
	using namespace chemical;
	using namespace conformation;
	using namespace pose;
	using namespace protocols;
	using namespace moves;
	using namespace ObjexxFCL::format;
	using numeric::random::uniform;

	devel::init(argc,argv);

	CCParam cp;
	Pose pose = make_coiled_coil(cp);
	pose.dump_pdb("test.pdb");


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
