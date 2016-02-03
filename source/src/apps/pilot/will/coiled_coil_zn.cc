// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilat/will/coiled_coil.cc
/// @brief samples coiled coils

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/fragment/FragmentIO.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <devel/init.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
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
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <protocols/toolbox/SwitchResidueTypeSet.hh>


static THREAD_LOCAL basic::Tracer TR( "coiled_coil" );

using core::Size;
using core::Real;
typedef numeric::xyzVector<Real> Vec;
typedef utility::vector1<Vec>    Vecs;
typedef numeric::xyzMatrix<Real> Mat;
using protocols::moves::MoverOP;
using core::scoring::ScoreFunctionOP;

Vec helix_axis(core::pose::Pose const & pose) {
	Vec axis(0,0,0);
	for(Size i = 1; i <= pose.n_residue()-4; ++i) {
		axis += ( pose.residue_type(i+4).xyz(1) - pose.residue(i).xyz(3) );
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
	core::pose::Pose pose;
	core::pose::make_pose_from_sequence(pose,seq,core::chemical::CENTROID,true);

	utility::vector1<Real> tmpphi(98,0.0),tmppsi(98,0.0),tmpomg(98,0.0);
	tmpphi[1 ] = -69.3868; tmppsi[1 ] = -37.2318; tmpomg[1 ] = 171.228;
	tmpphi[2 ] = -74.9636; tmppsi[2 ] = -49.2507; tmpomg[2 ] = 174.801;
	tmpphi[3 ] = -50.566 ; tmppsi[3 ] = -31.1118; tmpomg[3 ] = 179.936;
	tmpphi[4 ] = -71.1152; tmppsi[4 ] = -77.8957; tmpomg[4 ] = -177.588;
	tmpphi[5 ] = -39.9142; tmppsi[5 ] = -27.5922; tmpomg[5 ] = 178.027;
	tmpphi[6 ] = -63.6357; tmppsi[6 ] = -68.5906; tmpomg[6 ] = -173.661;
	tmpphi[7 ] = -58.1097; tmppsi[7 ] = -20.1639; tmpomg[7 ] = 166.724;
	tmpphi[8 ] = -66.8294; tmppsi[8 ] = -35.5721; tmpomg[8 ] = 175.247;
	tmpphi[9 ] = -66.9206; tmppsi[9 ] = -48.3089; tmpomg[9 ] = 179.32;
	tmpphi[10] = -55.5917; tmppsi[10] = -60.8831; tmpomg[10] = 178.377;
	tmpphi[11] = -49.067 ; tmppsi[11] = -38.5048; tmpomg[11] = -178.697;
	tmpphi[12] = -54.9044; tmppsi[12] = -65.6644; tmpomg[12] = 179.703;
	tmpphi[13] = -50.5591; tmppsi[13] = -37.7453; tmpomg[13] = -166.71;
	tmpphi[14] = -75.4586; tmppsi[14] = -31.3726; tmpomg[14] = 167.743;
	tmpphi[15] = -67.4081; tmppsi[15] = -39.1921; tmpomg[15] = 172.28;
	tmpphi[16] = -58.2457; tmppsi[16] = -47.3504; tmpomg[16] = 176.534;
	tmpphi[17] = -58.0694; tmppsi[17] = -47.2142; tmpomg[17] = -178.599;
	tmpphi[18] = -57.3874; tmppsi[18] = -42.3596; tmpomg[18] = 178.427;
	tmpphi[19] = -57.2778; tmppsi[19] = -37.1843; tmpomg[19] = -177.884;
	tmpphi[20] = -71.8025; tmppsi[20] = -42.2116; tmpomg[20] = -179.249;
	tmpphi[21] = -62.6323; tmppsi[21] = -42.4217; tmpomg[21] = 169.996;
	tmpphi[22] = -57.6275; tmppsi[22] = -37.9848; tmpomg[22] = -176.473;
	tmpphi[23] = -69.3893; tmppsi[23] = -45.9868; tmpomg[23] = 176.756;
	tmpphi[24] = -64.1536; tmppsi[24] = -43.659 ; tmpomg[24] = -178.45;
	tmpphi[25] = -58.063 ; tmppsi[25] = -49.091 ; tmpomg[25] = 178.833;
	tmpphi[26] = -54.1827; tmppsi[26] = -45.4343; tmpomg[26] = 179.171;
	tmpphi[27] = -61.254 ; tmppsi[27] = -39.9318; tmpomg[27] = 178.987;
	tmpphi[28] = -63.9572; tmppsi[28] = -48.9621; tmpomg[28] = 179.806;
	tmpphi[29] = -59.0385; tmppsi[29] = -40.852 ; tmpomg[29] = 176.567;
	tmpphi[30] = -65.5365; tmppsi[30] = -46.0753; tmpomg[30] = 171.496;
	tmpphi[31] = -51.2758; tmppsi[31] = -55.4363; tmpomg[31] = -179.938;
	tmpphi[32] = -49.7449; tmppsi[32] = -50.1227; tmpomg[32] = -178.027;
	tmpphi[33] = -54.6838; tmppsi[33] = -51.7226; tmpomg[33] = 177.936;
	tmpphi[34] = -54.0262; tmppsi[34] = -58.297 ; tmpomg[34] = -175.144;
	tmpphi[35] = -55.5292; tmppsi[35] = -37.5977; tmpomg[35] = -177.94;
	tmpphi[36] = -70.9469; tmppsi[36] = -38.6609; tmpomg[36] = 171.637;
	tmpphi[37] = -63.1542; tmppsi[37] = -48.1853; tmpomg[37] = 177.163;
	tmpphi[38] = -57.3797; tmppsi[38] = -45.2332; tmpomg[38] = 179.523;
	tmpphi[39] = -63.1354; tmppsi[39] = -36.7263; tmpomg[39] = 175.482;
	tmpphi[40] = -61.251 ; tmppsi[40] = -47.5537; tmpomg[40] = 175.511;
	tmpphi[41] = -55.7829; tmppsi[41] = -53.0446; tmpomg[41] = -175.078;
	tmpphi[42] = -58.7367; tmppsi[42] = -36.4854; tmpomg[42] = -179.347;
	tmpphi[43] = -74.9232; tmppsi[43] = -28.9601; tmpomg[43] = 171.245;
	tmpphi[44] = -67.8677; tmppsi[44] = -49.0195; tmpomg[44] = 173.924;
	tmpphi[45] = -61.6278; tmppsi[45] = -38.1132; tmpomg[45] = 176.816;
	tmpphi[46] = -62.8141; tmppsi[46] = -50.6933; tmpomg[46] = 174.209;
	tmpphi[47] = -49.6604; tmppsi[47] = -46.9487; tmpomg[47] = -175.874;
	tmpphi[48] = -63.1361; tmppsi[48] = -50.1631; tmpomg[48] = -175.802;
	tmpphi[49] = -58.9914; tmppsi[49] = -38.8711; tmpomg[49] = 179.107;
	tmpphi[50] = -61.1823; tmppsi[50] = -46.9917; tmpomg[50] = 176.944;
	tmpphi[51] = -55.8222; tmppsi[51] = -52.0404; tmpomg[51] = -179.515;
	tmpphi[52] = -58.1787; tmppsi[52] = -42.9166; tmpomg[52] = -179.104;
	tmpphi[53] = -63.5307; tmppsi[53] = -33.8841; tmpomg[53] = 177.535;
	tmpphi[54] = -67.4866; tmppsi[54] = -50.5217; tmpomg[54] = 177.468;
	tmpphi[55] = -54.0615; tmppsi[55] = -63.3318; tmpomg[55] = -179.227;
	tmpphi[56] = -49.6363; tmppsi[56] = -44.352 ; tmpomg[56] = 179.875;
	tmpphi[57] = -58.303;  tmppsi[57] = -43.9049; tmpomg[57] = 168.371;
	tmpphi[58] = -53.8322; tmppsi[58] = -61.1025; tmpomg[58] = -179.348;
	tmpphi[59] = -47.9876; tmppsi[59] = -48.9279; tmpomg[59] = -177.232;
	tmpphi[60] = -56.8017; tmppsi[60] = -42.6614; tmpomg[60] = -179.848;
	tmpphi[61] = -64.1835; tmppsi[61] = -52.2563; tmpomg[61] = 179.246;
	tmpphi[62] = -55.37;   tmppsi[62] = -46.1485; tmpomg[62] = -179.665;
	tmpphi[63] = -63.891;  tmppsi[63] = -31.0802; tmpomg[63] = 178.843;
	tmpphi[64] = -71.39;   tmppsi[64] = -43.0953; tmpomg[64] = 174.516;
	tmpphi[65] = -55.3363; tmppsi[65] = -54.6028; tmpomg[65] = 179.936;
	tmpphi[66] = -53.6256; tmppsi[66] = -45.5313; tmpomg[66] = -175.104;
	tmpphi[67] = -63.5452; tmppsi[67] = -41.9229; tmpomg[67] = 178.466;
	tmpphi[68] = -66.5179; tmppsi[68] = -37.4687; tmpomg[68] = 175.584;
	tmpphi[69] = -60.6738; tmppsi[69] = -44.5559; tmpomg[69] = -179.889;
	tmpphi[70] = -65.6707; tmppsi[70] = -40.7223; tmpomg[70] = 177.846;
	tmpphi[71] = -58.901;  tmppsi[71] = -39.7334; tmpomg[71] = 174.322;
	tmpphi[72] = -67.4215; tmppsi[72] = -49.0965; tmpomg[72] = -179.667;
	tmpphi[73] = -57.7398; tmppsi[73] = -46.0866; tmpomg[73] = -178.983;
	tmpphi[74] = -60.9306; tmppsi[74] = -47.2598; tmpomg[74] = 178.843;
	tmpphi[75] = -55.8155; tmppsi[75] = -47.8979; tmpomg[75] = -179.639;
	tmpphi[76] = -59.9858; tmppsi[76] = -45.5085; tmpomg[76] = -178.519;
	tmpphi[77] = -64.2232; tmppsi[77] = -44.4745; tmpomg[77] = 179.473;
	tmpphi[78] = -61.4602; tmppsi[78] = -37.5473; tmpomg[78] = 178.294;
	tmpphi[79] = -64.8219; tmppsi[79] = -47.4759; tmpomg[79] = 177.284;
	tmpphi[80] = -63.4321; tmppsi[80] = -38.5366; tmpomg[80] = 176.909;
	tmpphi[81] = -59.3389; tmppsi[81] = -46.5287; tmpomg[81] = 179.374;
	tmpphi[82] = -62.3524; tmppsi[82] = -44.0625; tmpomg[82] = 178.957;
	tmpphi[83] = -63.7316; tmppsi[83] = -33.8172; tmpomg[83] = 178.284;
	tmpphi[84] = -71.527;  tmppsi[84] = -35.4306; tmpomg[84] = 177.693;
	tmpphi[85] = -63.6935; tmppsi[85] = -41.2806; tmpomg[85] = 174.469;
	tmpphi[86] = -64.5189; tmppsi[86] = -51.5024; tmpomg[86] = 176.156;
	tmpphi[87] = -51.4985; tmppsi[87] = -48.6718; tmpomg[87] = 176.922;
	tmpphi[88] = -52.4326; tmppsi[88] = -63.4492; tmpomg[88] = 178.796;
	tmpphi[89] = -42.2522; tmppsi[89] = -56.5929; tmpomg[89] = 176.315;
	tmpphi[90] = -46.0111; tmppsi[90] = -42.5002; tmpomg[90] = -171.783;
	tmpphi[91] = -71.787;  tmppsi[91] = -46.7328; tmpomg[91] = -178.402;
	tmpphi[92] = -69.1741; tmppsi[92] = -19.614 ; tmpomg[92] = 174.636;
	tmpphi[93] = -82.3169; tmppsi[93] = -36.7285; tmpomg[93] = 173.756;
	tmpphi[94] = -65.6821; tmppsi[94] = -31.1765; tmpomg[94] = 168.989;
	tmpphi[95] = -65.3523; tmppsi[95] = -58.1832; tmpomg[95] = 176.046;
	tmpphi[96] = -41.9;    tmppsi[96] = -69.0843; tmpomg[96] = -177.093;
	tmpphi[97] = -48.6064; tmppsi[97] = -50.5125; tmpomg[97] = -179.32;
	tmpphi[98] = -44.7473; tmppsi[98] = -53.924 ; tmpomg[98] = 170.504;

	Size pos = (Size)floor(numeric::random::uniform() * (tmpphi.size()-pose.n_residue()+0.99999999));
	for ( core::Size i = 1; i <= pose.n_residue(); i++ ) {
		pose.set_phi  ( i,tmpphi[i+pos] );
		pose.set_psi  ( i,tmppsi[i+pos] );
		pose.set_omega( i,tmpomg[i+pos] );
	}
	trans_pose(pose,-center_of_mass(pose));
	Vec axis = helix_axis(pose);
	Vec z = Vec(0,0,1);
	Vec rot_axis = z.cross(axis);
	Real rot_ang = acos(z.dot(axis));
	Mat rot = rotation_matrix_radians( rot_axis, -rot_ang );
	rot_pose(pose,rot);
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
	s += "E = 1.0*VRT1";
	for(Size i = 2; i<=n; i++) s+= " + 1*(VRT1:VRT"+string_of(i)+")";
	s += "\nanchor_residue " + string_of(anchor) + "\n";
	s += "virtual_transforms_start consecutive\nstart -1,0,0 0,1,0 0,0,0\n";
	s += "rot Rz_angle " + string_of(rot) + "\n";
	s += "trans 0,0," + string_of(trans);
	s += "\nvirtual_transforms_stop\n";
	for(Size i = 2; i<=n; i++) s+= "connect_virtual J"+string_of(i)+" VRT"+string_of(i-1)+" VRT"+string_of(i)+"\n";
	s += "set_dof BASEJUMP x angle_x angle_y angle_z\n";
	s += "set_dof J2 z angle_z\n";
	// TR << "================= symm dat ==================" << std::endl;
	// TR << s << std::endl;
	// TR << "================= symm dat ==================" << std::endl;
	std::istringstream iss(s);
	core::conformation::symmetry::SymmData symdat( pose.n_residue(), pose.num_jump() );
//	symdat.read_symmetry_data_from_stream(iss);
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
		nres = 15.0 + 35.0*(uniform());
		rot = 120.0 + 30.0*gaussian();
		x = 6.0 + 2.0*gaussian();
		rh = 360.0*uniform();
		rhx = 7.0*gaussian()+20.0*32.0/(Real)nres;
		if(uniform()<0.5) rot *= -1.0;
		if(uniform()<0.5) rhx *= -1.0;
		rhy =  6.0*gaussian();
		trans = 7.0 + 4*uniform();
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
	utility::vector1<Size> vcb_pos;
};

core::pose::Pose make_coiled_coil(CCParam & p) {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace pose;
	using namespace conformation::symmetry;

	std::string seq = "";
	for(Size i = 1; i <= p.nres; ++i) seq += "A";
	for(Size i = 1; i <= p.vcb_pos.size(); ++i) seq[p.vcb_pos[i]-1] = 'C';
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


core::kinematics::MoveMapOP make_move_map( core::pose::Pose & pose ) {
	using namespace core;
	using namespace core::conformation::symmetry;

	kinematics::MoveMapOP movemap = new kinematics::MoveMap();
	movemap->set_bb( true );
	movemap->set_chi( true );
	movemap->set_jump( false );

	assert( is_symmetric( pose ) );
	SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
	SymmetryInfo const & symm_info( symm_conf.Symmetry_Info() );
	std::map< Size, SymDof > dofs ( symm_info.get_dofs() );
	std::map< Size, SymDof >::iterator it;
	std::map< Size, SymDof >::iterator it_begin = dofs.begin();
	std::map< Size, SymDof >::iterator it_end   = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first );
		SymDof dof( (*it).second );
		if(dof.allow_dof(X_DOF      )) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,1)),true);
		if(dof.allow_dof(Y_DOF      )) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,2)),true);
		if(dof.allow_dof(Z_DOF      )) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,3)),true);
		if(dof.allow_dof(X_ANGLE_DOF)) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,4)),true);
		if(dof.allow_dof(Y_ANGLE_DOF)) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,5)),true);
		if(dof.allow_dof(Z_ANGLE_DOF)) movemap->set(pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,6)),true);
	}

	// if jump is between independent, non-virtual residues, move it
	for(Size i = 1; i <= pose.num_jump(); ++i) {
   	Size res1 = pose.fold_tree().upstream_jump_residue( i );
   	Size res2 = pose.fold_tree().downstream_jump_residue( i );
		if ( symm_info.fa_is_independent(res1) && !symm_info.is_virtual(res1) &&
		     symm_info.fa_is_independent(res2) && !symm_info.is_virtual(res2)  ) {
			movemap->set_jump(i,true);
		}
	}

	return movemap;
}

Real rg2d(core::pose::Pose & pose) {
	Real rg = 0;
	Size count = 0;
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= numeric::min(pose.residue_type(ir).natoms(),(Size)5); ++ia) {
			Vec v = pose.xyz(core::id::AtomID(ia,ir));
			rg += sqrt(sqrt(sqrt(v.x()*v.x()+v.y()*v.y())));
			count++;
		}
	}
	rg /= count;
	return rg*rg*rg*rg;
}

void addcc(core::pose::Pose & pose, core::id::AtomID aid, core::id::AtomID anchor, core::Real mult = 1.0 ) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::CoordinateConstraint(
		aid, anchor, pose.xyz(aid), new core::scoring::constraints::HarmonicFunc(0,mult) );
	pose.add_constraint(cc);
}

void add_apc(core::pose::Pose & pose, core::id::AtomID aid1, core::id::AtomID aid2, core::Real mean, core::Real sd) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AtomPairConstraint(
		aid1, aid2, new core::scoring::constraints::HarmonicFunc(mean,sd) );
	pose.add_constraint(cc);
}

void
add_agc(
	core::pose::Pose & pose,
	core::id::AtomID aid1,
	core::id::AtomID aid2,
	core::id::AtomID aid3,
	core::Real mean,
	core::Real sd
) {
	core::scoring::constraints::ConstraintOP cc = new core::scoring::constraints::AngleConstraint(
		aid1, aid2, aid3, new core::scoring::constraints::HarmonicFunc(mean,sd) );
	pose.add_constraint(cc);
}

void
add_dhc(
	core::pose::Pose & pose,
	core::id::AtomID aid1,
	core::id::AtomID aid2,
	core::id::AtomID aid3,
	core::id::AtomID aid4,
	core::Real mean,
	core::Real sd
) {
	using namespace core::scoring::constraints;
	core::scoring::constraints::ConstraintOP dh1 = new DihedralConstraint(
		aid1, aid2, aid3, aid4, new CircularHarmonicFunc(mean,sd) );
	core::scoring::constraints::ConstraintOP dh2 = new DihedralConstraint(
		aid1, aid2, aid3, aid4, new CircularHarmonicFunc(-mean,sd) );
	ConstraintCOPs csts;
	csts.push_back(dh1);
	csts.push_back(dh2);
	pose.add_constraint( new AmbiguousConstraint(csts) );
}


void repack(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	task->restrict_to_repacking();
	task->or_include_current(true);
	// for(Size i = 1; i <= task->total_residue(); ++i) {
	// 	// if(cc.residue(i).name3().substr(0,2)=="CY") {
	// 		// task->nonconst_residue_task(i).or_ex1(true);
	// 		// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
	// 	}
	// }
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_target(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	utility::vector1<Size> aalist = basic::options::option[basic::options::OptionKeys::in::target_residues]();
	for(Size i = 1; i <= aalist.size(); ++i) {
		aas[aalist[i]] = true;
	}
	utility::vector1< bool > cys(20,false);
	cys[core::chemical::aa_cys] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_all(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	task->or_include_current(true);
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		}
	}
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_FILV(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	aas[core::chemical::aa_phe] = true;
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_val] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_AFILV(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	aas[core::chemical::aa_ala] = true;
	aas[core::chemical::aa_phe] = true;
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_val] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_AL(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	aas[core::chemical::aa_ala] = true;
	aas[core::chemical::aa_leu] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_FILVEK(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	aas[core::chemical::aa_phe] = true;
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_val] = true;
	aas[core::chemical::aa_lys] = true;
	aas[core::chemical::aa_glu] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}

void design_AFILVEK(core::pose::Pose & cc, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(cc);
	utility::vector1< bool > aas(20,false);
	aas[core::chemical::aa_ala] = true;
	aas[core::chemical::aa_phe] = true;
	aas[core::chemical::aa_ile] = true;
	aas[core::chemical::aa_leu] = true;
	aas[core::chemical::aa_val] = true;
	aas[core::chemical::aa_lys] = true;
	aas[core::chemical::aa_glu] = true;
	for(Size i = 1; i <= task->total_residue(); ++i) {
		if(cc.residue(i).name3().substr(0,2)=="CY") {
			task->nonconst_residue_task(i).restrict_to_repacking();
			// task->nonconst_residue_task(i).or_ex1(true);
			// task->nonconst_residue_task(i).or_ex1_sample_level(EX_FOUR_HALF_STEP_STDDEVS);
		} else {
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
		}
	}
	task->or_include_current(true);
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(cc);
}


void minimize(core::pose::Pose & cc, ScoreFunctionOP sf) {
	core::kinematics::MoveMapOP movemap = make_move_map(cc);
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-2, true );
	m.apply(cc);
}

Size natoms(core::pose::Pose & cc, Size nres) {
	Size count = 0;
	for(Size ir = 1; ir <= nres; ++ir) count += cc.residue_type(ir).nheavyatoms();
	return count;
}

Real zcyl_score(core::pose::Pose & cc, Size nres, Real trans) {
	Real score = 0;
	for(Size ir = 1; ir <= nres; ++ir) {
		for(Size ia = 1; ia <= cc.residue_type(ir).nheavyatoms(); ++ia) {
			Vec v = cc.xyz(core::id::AtomID(ia,ir));
			v.z(0);
			score += exp(-v.length()*v.length()/16.0);
		}
	}
	return score / trans;
}

void add_lig_cst( core::pose::Pose & pose, Size r1, Size r2, Size r3, Size r4, Size l1, Real bondwt = 1.0 ) {
	using core::id::AtomID;
	using numeric::conversions::radians;
	// // sf4
	// core::Real cb_sg      = 1.804057, cb_sg_sd      = 0.036217;
	// core::Real cb_sg_m    = 108.2215, cb_sg_m       = 9.0377 ;
	// core::Real ca_cb_sg   = 114.1297, ca_cb_sg      = 3.232815;
	// core::Real ca_cb_sg_m = 86.13981, ca_cb_sg_m_sd = 26.35887;
	// ZN
	Real    cb_sg   =           1.8113 ,    cb_sg_sd   =          0.03128  * 1.0;
	Real    cb_sg_m = radians(106.7646),    cb_sg_m_sd = radians(13.43894) * 1.0;
	Real ca_cb_sg   = radians(114.3113), ca_cb_sg_sd   = radians( 3.26554) * 1.0;
	// Real ca_cb_sg_m = radians( 95.6460), ca_cb_sg_m_sd = radians(25.57046) * 1.0;
	Real    cb_cg   =           1.49743,    cb_cg_sd   =       0.006614549 * 5.0;
	Real    cb_cg_n = radians(113.5352),    cb_cg_n_sd = radians(0.968983) * 3.0;
	Real ca_cb_cg   = radians(123.0981), ca_cb_cg_sd   = radians(0.892374) * 3.0;
	add_apc(pose,               AtomID(5,r1), AtomID( 6,l1),                   cb_cg  ,    cb_cg_sd / bondwt );
	add_apc(pose,               AtomID(5,r2), AtomID(11,l1),                   cb_cg  ,    cb_cg_sd / bondwt );
	add_apc(pose,               AtomID(5,r3), AtomID( 2,l1),                   cb_sg  ,    cb_sg_sd / bondwt );
	add_apc(pose,               AtomID(5,r4), AtomID( 3,l1),                   cb_sg  ,    cb_sg_sd / bondwt );
	add_agc(pose,               AtomID(5,r1), AtomID( 6,l1), AtomID( 7,l1),    cb_cg_n,    cb_cg_n_sd        );
	add_agc(pose,               AtomID(5,r2), AtomID(11,l1), AtomID(12,l1),    cb_cg_n,    cb_cg_n_sd        );
	add_agc(pose,               AtomID(5,r1), AtomID( 6,l1), AtomID( 5,l1),    cb_cg_n,    cb_cg_n_sd        );
	add_agc(pose,               AtomID(5,r2), AtomID(11,l1), AtomID(10,l1),    cb_cg_n,    cb_cg_n_sd        );
	add_agc(pose,               AtomID(5,r3), AtomID( 2,l1), AtomID( 1,l1),    cb_sg_m,    cb_sg_m_sd        );
	add_agc(pose,               AtomID(5,r4), AtomID( 3,l1), AtomID( 1,l1),    cb_sg_m,    cb_sg_m_sd        );
	add_agc(pose, AtomID(2,r1), AtomID(5,r1), AtomID( 6,l1),                ca_cb_cg  , ca_cb_cg_sd          );
	add_agc(pose, AtomID(2,r2), AtomID(5,r2), AtomID(11,l1),                ca_cb_cg  , ca_cb_cg_sd          );
	add_agc(pose, AtomID(2,r3), AtomID(5,r3), AtomID( 2,l1),                ca_cb_sg  , ca_cb_sg_sd          );
	add_agc(pose, AtomID(2,r4), AtomID(5,r4), AtomID( 3,l1),                ca_cb_sg  , ca_cb_sg_sd          );
	// // add_dhc(pose, AtomID(2,r1), AtomID(5,r1), AtomID(1,l1), AtomID(5,l1), ca_cb_sg_m, ca_cb_sg_m_sd        );
	// // add_dhc(pose, AtomID(2,r2), AtomID(5,r2), AtomID(2,l1), AtomID(5,l1), ca_cb_sg_m, ca_cb_sg_m_sd        );
	// add_dhc(pose, AtomID(2,r3), AtomID(5,r3), AtomID( 2,l1), AtomID( 1,l1), ca_cb_sg_m, ca_cb_sg_m_sd        );
	// add_dhc(pose, AtomID(2,r4), AtomID(5,r4), AtomID( 3,l1), AtomID( 1,l1), ca_cb_sg_m, ca_cb_sg_m_sd        );
}

void add_fa_cst(core::pose::Pose & cc, CCParam & p) {
	using core::id::AtomID;
	cc.remove_constraints();
	for(Size j = 1; j <= p.nres; ++j) addcc(cc,AtomID(2,j),AtomID(1,cc.n_residue()+1-p.nsub),2.0);

	//std::cerr << p.vcb_pos.size() << std::endl;
	Size c1 = p.vcb_pos[1]; assert("CY"==cc.residue(c1).name3().substr(0,2));
	Size c2 = p.vcb_pos[2]; assert("CY"==cc.residue(c2).name3().substr(0,2));
	Size c3 = p.vcb_pos[3]; assert("CY"==cc.residue(c3).name3().substr(0,2));
	Size c4 = p.vcb_pos[4]; assert("CY"==cc.residue(c4).name3().substr(0,2));

	add_lig_cst(cc,c1,c2,c3,c4,p.nres+1,0.3);
}


typedef ScoreFunctionOP SFOP;
void make_sf(SFOP& sfc, SFOP& sfd, SFOP& sf1, SFOP& sf2, SFOP& sf3, SFOP& sf4, SFOP& sf5 ) {
	using namespace core::scoring;

	sfc = ScoreFunctionFactory::create_score_function( CENTROID_WTS ) ;
	sfd = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS ) ;
	sf1 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf2 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf3 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf4 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );
	sf5 = get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS );

	// sfc->set_weight(metal_placement,1.0);

	Real gbl_cst_wt = 0.6;
	sfc->set_weight(coordinate_constraint,gbl_cst_wt*1.0);
	sfd->set_weight(coordinate_constraint,gbl_cst_wt*0.1);
	sf1->set_weight(coordinate_constraint,gbl_cst_wt*0.2);
	sf2->set_weight(coordinate_constraint,gbl_cst_wt*0.3);
	sf3->set_weight(coordinate_constraint,gbl_cst_wt*0.5);
	sf4->set_weight(coordinate_constraint,gbl_cst_wt*0.7);
	sf5->set_weight(coordinate_constraint,gbl_cst_wt*1.0);

	sfc->set_weight(atom_pair_constraint,gbl_cst_wt*1.0);
	sfd->set_weight(atom_pair_constraint,gbl_cst_wt*0.8);
	sf1->set_weight(atom_pair_constraint,gbl_cst_wt*0.9);
	sf2->set_weight(atom_pair_constraint,gbl_cst_wt*1.1);
	sf3->set_weight(atom_pair_constraint,gbl_cst_wt*1.3);
	sf4->set_weight(atom_pair_constraint,gbl_cst_wt*1.6);
	sf5->set_weight(atom_pair_constraint,gbl_cst_wt*2.0);

	sfc->set_weight(angle_constraint,gbl_cst_wt*1.0);
	sfd->set_weight(angle_constraint,gbl_cst_wt*0.1);
	sf1->set_weight(angle_constraint,gbl_cst_wt*0.2);
	sf2->set_weight(angle_constraint,gbl_cst_wt*0.3);
	sf3->set_weight(angle_constraint,gbl_cst_wt*0.4);
	sf4->set_weight(angle_constraint,gbl_cst_wt*0.7);
	sf5->set_weight(angle_constraint,gbl_cst_wt*1.0);

	sfc->set_weight(dihedral_constraint,gbl_cst_wt*1.0);
	sfd->set_weight(dihedral_constraint,gbl_cst_wt*0.1);
	sf1->set_weight(dihedral_constraint,gbl_cst_wt*0.1);
	sf2->set_weight(dihedral_constraint,gbl_cst_wt*0.2);
	sf3->set_weight(dihedral_constraint,gbl_cst_wt*0.3);
	sf4->set_weight(dihedral_constraint,gbl_cst_wt*0.4);
	sf5->set_weight(dihedral_constraint,gbl_cst_wt*1.0);


	// Real orig_srbb = sf5->get_weight(hbond_sr_bb);
	// sfd->set_weight(hbond_sr_bb,2.0*orig_srbb);
	// sf1->set_weight(hbond_sr_bb,2.0*orig_srbb);
	// sf2->set_weight(hbond_sr_bb,2.0*orig_srbb);
	// sf3->set_weight(hbond_sr_bb,2.0*orig_srbb);
	// sf4->set_weight(hbond_sr_bb,2.0*orig_srbb);
	// sf5->set_weight(hbond_sr_bb,2.0*orig_srbb);


//  	Real dun_orig = sf5->get_weight(fa_dun);
//  	sfd->set_weight(fa_dun,0.1*dun_orig);
//  	sf1->set_weight(fa_dun,0.2*dun_orig);
//  	sf2->set_weight(fa_dun,0.3*dun_orig);
//  	sf3->set_weight(fa_dun,0.4*dun_orig);
//  	sf4->set_weight(fa_dun,0.5*dun_orig);
//  	sf5->set_weight(fa_dun,0.6*dun_orig);

//  	Real atr_orig = sf5->get_weight(fa_atr);
//  	sfd->set_weight(fa_atr,1.2*atr_orig);
//  	sf1->set_weight(fa_atr,1.2*atr_orig);
//  	sf2->set_weight(fa_atr,1.2*atr_orig);
//  	sf3->set_weight(fa_atr,1.2*atr_orig);
//  	sf4->set_weight(fa_atr,1.2*atr_orig);
//  	sf5->set_weight(fa_atr,1.2*atr_orig);

//  	Real sol_orig = sf5->get_weight(fa_sol);
//  	sfd->set_weight(fa_sol,0.7*sol_orig);
//  	sf1->set_weight(fa_sol,0.7*sol_orig);
//  	sf2->set_weight(fa_sol,0.7*sol_orig);
//  	sf3->set_weight(fa_sol,0.7*sol_orig);
//  	sf4->set_weight(fa_sol,0.7*sol_orig);
//  	sf5->set_weight(fa_sol,0.7*sol_orig);

	sfd->set_weight(fa_rep,1.00);
	sf1->set_weight(fa_rep,0.01);
	sf2->set_weight(fa_rep,0.025);
	sf3->set_weight(fa_rep,0.050);
	sf4->set_weight(fa_rep,0.100);
	sf4->set_weight(fa_rep,0.200);

	sfc = new symmetry::SymmetricScoreFunction(*sfc);
	sfd = new symmetry::SymmetricScoreFunction(*sfd);
	sf1 = new symmetry::SymmetricScoreFunction(*sf1);
	sf2 = new symmetry::SymmetricScoreFunction(*sf2);
	sf3 = new symmetry::SymmetricScoreFunction(*sf3);
	sf4 = new symmetry::SymmetricScoreFunction(*sf4);
	sf5 = new symmetry::SymmetricScoreFunction(*sf5);

}

inline Real sq(Real x) { return x*x; }

core::Real
align_zns(
	core::pose::Pose const & cc,
	CCParam & p,
	core::pose::Pose & tp,
	Size c1, Size c2, Size c3, Size c4
) {
	using core::id::AtomID;
	using namespace core::conformation;
	using namespace core::chemical;
	// core::pose::Pose cc = inpose;
	core::scoring::ScoreFunctionOP sf = new core::scoring::ScoreFunction;
	sf->set_weight(core::scoring::atom_pair_constraint,1.0);
	sf->set_weight(core::scoring::angle_constraint    ,1.0);
	sf->set_weight(core::scoring::dihedral_constraint ,1.0);
	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_chi(false); mm->set_bb(false);	mm->set_jump(false);
	mm->set_jump(4,true);
	// mm->set(core::id::RB1,false);   mm->set(core::id::RB2,false);   mm->set(core::id::RB3,false);
	protocols::simple_moves::MinMover mnm( mm, sf, "dfpmin_armijo_nonmonotone", 1e-4, true );
	ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
	ResidueOP cys = ResidueFactory::create_residue( residue_set->name_map("CYV") );
	ResidueOP zns = ResidueFactory::create_residue( residue_set->name_map("ZHC") );
	// cc.replace_residue(c1,*cys,true);
	// cc.replace_residue(c2,*cys,true);
	// cc.replace_residue(c3,*cys,true);
	// cc.replace_residue(c4,*cys,true);
	tp.append_residue_by_bond(cc.residue(c1)  );
	tp.append_residue_by_jump(cc.residue(c2),1);
	tp.append_residue_by_jump(cc.residue(c3),1);
	tp.append_residue_by_jump(cc.residue(c4),1);
	tp.append_residue_by_jump(*zns,1);

	// for(Size i=1; i<= tp.fold_tree().num_jump(); ++i)
		// std::cerr << tp.fold_tree().jump_edge(i) << std::endl;

	add_lig_cst(tp,1,2,3,4,5,0.3);
	mnm.apply(tp);

	core::pose::Pose tmp = tp;
	tp.remove_constraints();
	add_lig_cst(tp,1,2,4,3,5,0.3);
	mnm.apply(tp);

	// std::cerr << "DBG " << (*sf)(tp) << " " << (*sf)(tmp) << std::endl;
	bool minus = true;
	if( (*sf)(tp) > (*sf)(tmp) ) {
		minus = false;
		tp = tmp;
	}
	Vec com(0,0,0);
	for(Size i = 1; i <= tp.residue_type(5).nheavyatoms(); ++i) com += tp.xyz(AtomID(i,5));
	com.z(0);
	using numeric::min;
	using numeric::max;
	// !!!!!!!!!!! assumes no ZHC in yet!!!
	Size nt =          min((c1-1)%p.nres,min((c2-1)%p.nres,min((c3-1)%p.nres,(c4-1)%p.nres)))+1;
	Size ct = p.nres - max((c1-1)%p.nres,max((c2-1)%p.nres,max((c3-1)%p.nres,(c4-1)%p.nres)));
	Real score = (*sf)(tp) + sq(com.length());
	Real tailbonus = + 50*( max(0.0,4.0-nt) + max(0.0,4.0-ct) );
	// std::cerr << "tailbonus " << tailbonus << std::endl;
	// if( nt <=4 ||  ct <= 4 )
	// 	std::cerr << "ntct " << nt << " " << ct << " " << score << " " << tailbonus << " " << std::endl;
	if(minus) return -max(1.0,score+tailbonus);
	return            max(1.0,score+tailbonus);
}


utility::vector1< utility::vector1<Size> >
place_zns(
	core::pose::Pose & cc,
	CCParam & p,
	core::Real cstthresh
) {
	utility::vector1< utility::vector1<Size> > idxs;
	Real bestscore=9e9;
	using core::id::AtomID;
	for(Size i =   (p.nres*p.nsub/2)+1; i <= (p.nres*(p.nsub/2+1))-4     ; ++i ) {
	Size j = i+4;/*for(Size j = i+1; j <= p.nres*p.nsub; ++j )*/ {
		// if( i%p.nres == j%p.nres ) continue;
		// if( THRESH < cc.xyz(AtomID(5,i)).distance_squared(cc.xyz(AtomID(5,j))) ) continue;
		// if( NEAR   > cc.xyz(AtomID(5,i)).distance_squared(cc.xyz(AtomID(5,j))) ) continue;
	for(Size k = 1; k <= p.nres*p.nsub; ++k ) {
		if( i%p.nres == k%p.nres ) continue; if( j%p.nres == k%p.nres ) continue;
		Real const dik = cc.xyz(AtomID(5,i)).distance_squared(cc.xyz(AtomID(5,k)));
		Real const djk = cc.xyz(AtomID(5,j)).distance_squared(cc.xyz(AtomID(5,k)));
		if( 81.0 < dik ) continue;	if( 81.0 < djk ) continue;
	for(Size l = k+1; l <= p.nres*p.nsub; ++l ) {
		if( i%p.nres == l%p.nres ) continue; if( j%p.nres == l%p.nres ) continue; if( k%p.nres == l%p.nres ) continue;
		Real const dil = cc.xyz(AtomID(5,i)).distance_squared(cc.xyz(AtomID(5,l)));
		Real const djl = cc.xyz(AtomID(5,j)).distance_squared(cc.xyz(AtomID(5,l)));
		if( 81.0 < dil ) continue;	if( 81.0 < djl ) continue;
		Real const dimn = numeric::min(dik,dil), dimx = numeric::max(dik,dil);
		Real const djmn = numeric::min(djk,djl), djmx = numeric::max(djk,djl);
		if( 49.0 < dimn )	continue; if( 81.0 < dimx ) continue;
		if( 56.3 < djmn ) continue; if( 72.3 < djmx ) continue;

		// check for space for ZHC
		Vec const a = cc.xyz(AtomID(5,i));
		Vec const b = cc.xyz(AtomID(5,j));
		Vec const c = cc.xyz(AtomID(5,k));
		Vec const d = cc.xyz(AtomID(5,l));
		Vec const aa = cc.xyz(AtomID(2,i));
		Vec const ba = cc.xyz(AtomID(2,j));
		Vec const ca = cc.xyz(AtomID(2,k));
		Vec const da = cc.xyz(AtomID(2,l));
		Vec com = (a+b+c+d)/4.0;

		// check CA-CB pointing roughly right way
		if( (com-a).dot(a-aa) < 0.0 ) continue;
		if( (com-b).dot(b-ba) < 0.0 ) continue;
		if( (com-c).dot(c-ca) < 0.0 ) continue;
		if( (com-d).dot(d-da) < 0.0 ) continue;

		// std::cout << "checking zns spot " << i << " " << j << " " << k << " " << l;

		bool collision = false;
		for(Size ir = 1; ir <= p.nres*p.nsub; ++ir) {
			for(Size ia = 1; ia <= cc.residue_type(ir).nheavyatoms(); ++ia) {
				if( com.distance_squared(cc.xyz(AtomID(ia,ir))) < 4.0 ) {
					collision = true;
				}
			}
		}
		if( collision ) {
			// std::cerr << " collision" << std::endl;
			continue;
		}

		core::pose::Pose tmp;
		Real znsscore = fabs(align_zns(cc,p,tmp,i,j,k,l));
		// std::cout << "align_zns score " << znsscore << std::endl;
		// std::cerr << "znsscore " << znsscore << " "
		//           << i+floor((i-1)/(p.nres)) << " "
		//           << j+floor((j-1)/(p.nres)) << " "
		//           << k+floor((k-1)/(p.nres)) << " "
		//           << l+floor((l-1)/(p.nres)) << std::endl;
		if(znsscore < bestscore) bestscore = znsscore;
		if(znsscore < cstthresh) {
			utility::vector1<Size> idx;
			idx.push_back(i);
			idx.push_back(j);
			idx.push_back(k);
			idx.push_back(l);
			idxs.push_back(idx);
		}
		// std::cerr << " " << bestm << " " << bestatpos << std::endl;
	} // l
	} // k
	} // j
	} // i
	/*if(idxs.size()==0) */TR << "place ZHC score " << bestscore << std::endl;
	return idxs;
}


bool
add_symm_zns(
	core::pose::Pose & cc,
	CCParam & p,
	SFOP sf,
	core::Real cstthresh
) {
	using core::id::AtomID;
	using namespace core::conformation;
	using namespace core::chemical;
	core::Real bestscore = 9e9;
	core::pose::Pose bestpose;
	Size bestii = 0;
	utility::vector1<utility::vector1<Size> > idxs( place_zns(cc,p,cstthresh) );
	utility::vector1<Size> cpos(4);
	if( 0 == idxs.size() ) return false;
	for( Size ii = 1; ii <= idxs.size(); ++ii ) {
		core::pose::Pose tmppose = cc;
		utility::vector1<Size> t = idxs[ii];
		Size c1=t[1],c2=t[2],c3=t[3],c4=t[4];

		core::pose::Pose tp;
		// std::cerr << "place " << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
		Real minus = align_zns(tmppose,p,tp,c1,c2,c3,c4);

		// make symm pose w/ZHC
		ResidueTypeSetCAP residue_set( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );
		ResidueOP zns = ResidueFactory::create_residue( residue_set->name_map("ZHC") );
		core::pose::replace_pose_residue_copying_existing_coordinates(tmppose,(c1-1)%p.nres+1,residue_set->name_map("CYV"));
		core::pose::replace_pose_residue_copying_existing_coordinates(tmppose,(c2-1)%p.nres+1,residue_set->name_map("CYV"));
		core::pose::replace_pose_residue_copying_existing_coordinates(tmppose,(c3-1)%p.nres+1,residue_set->name_map("CYV"));
		core::pose::replace_pose_residue_copying_existing_coordinates(tmppose,(c4-1)%p.nres+1,residue_set->name_map("CYV"));
		std::string seq = tmppose.sequence().substr(0,p.nres);

		core::pose::Pose pose;
		core::pose::make_pose_from_sequence(pose,seq,core::chemical::CENTROID,true);
		pose.copy_segment(p.nres,tmppose,1,1);
		for(Size i = 1; i <= zns->nheavyatoms(); ++i ) {
			zns->set_xyz(i,tp.xyz(AtomID(i,5)));
		}
		using numeric::min;
		pose.append_residue_by_jump( *zns, min(min(min(p.nres,c1),c3),c4) );
		core::conformation::symmetry::SymmData sd = make_symm_data(pose,p.rot,p.trans,p.nsub);
		core::pose::symmetry::make_symmetric_pose( pose, sd );

		// for(Size i = 1; i <= pose.num_jump(); ++i){
		// 	TR << "fold tree jump " << pose.fold_tree().jump_edge(i) << std::endl;
		// }

		// if jump is between independent, non-virtual residues, move it
		using namespace core::conformation::symmetry;
		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
		SymmetryInfo & symm_info( symm_conf.Symmetry_Info() );
		std::map< Size, SymDof > dofs ( symm_info.get_dofs() );
		for(Size i = 1; i <= pose.num_jump(); ++i) {
	   	Size res1 = pose.fold_tree().upstream_jump_residue( i );
	   	Size res2 = pose.fold_tree().downstream_jump_residue( i );
			// TR << "checking jump " << i << " " << res1 << " " << res2 << std::endl;
			if ( symm_info.fa_is_independent(res1) && !symm_info.is_virtual(res1) &&
			     symm_info.fa_is_independent(res2) && !symm_info.is_virtual(res2)  ) {
				SymDof d;
				// TR << "adding SymDOF for jump " << i << std::endl;
				d.read("x y z angle_x angle_y angle_z");
				dofs.insert(std::make_pair( i, d ));
			}
		}
		symm_info.set_dofs( dofs );


		tmppose = pose;
		tmppose.remove_constraints();
		c1 = c1+floor((c1-1)/(p.nres)); assert("CY"==tmppose.residue(c1).name3().substr(0,2));
		c2 = c2+floor((c2-1)/(p.nres)); assert("CY"==tmppose.residue(c2).name3().substr(0,2));
		c3 = c3+floor((c3-1)/(p.nres)); assert("CY"==tmppose.residue(c3).name3().substr(0,2));
		c4 = c4+floor((c4-1)/(p.nres)); assert("CY"==tmppose.residue(c4).name3().substr(0,2));
		if( minus < 0 ) { Size tmp = c3; c3 = c4; c4 = tmp; }
		cpos[1] = c1; cpos[2] = c2; cpos[3] = c3; cpos[4] = c4;

		add_lig_cst(tmppose,c1,c2,c3,c4,p.nres+1,0.1);

		core::Real score = (*sf)(tmppose);
		// std::cout << ii << std::endl;
		// sf->show(tmppose);
		// tmppose.dump_pdb("tmppose.pdb");
		// std::cerr << "zns censcore " << score/p.nres << " " << c1 << " " << c2 << " " << c3 << " " << c4 << std::endl;
		if(score < bestscore) {
			bestscore = score;
			bestpose  = tmppose;
			bestii    = ii;
			p.vcb_pos = cpos;
		}
	}
	cc = bestpose;

	// std::cerr << "bestii " << bestii << std::endl;
	return true;
}


void read_frags() {
	utility::io::izstream in( basic::database::full_name("sampling/ss_fragfiles/HHH.fragfile") );
	Size n;
	utility::vector1< utility::pointer::owning_ptr<core::fragment::BBTorsionSRFD> > fds;
	while( in >> n ) {
	 	std::string pdb;

		std::cerr << n << std::endl;
		for( Size i = 1; i <= n; ++i ) {
			fds.push_back( new core::fragment::BBTorsionSRFD );
			in >> pdb >> *fds.back();
			std::cerr << pdb << " " << *fds.back() << std::endl;
		}
	}
	in.close();
}

void do_centroid_stuff(core::pose::Pose & cenpose, ScoreFunctionOP sfc ) {
	using namespace protocols;
	using namespace abinitio;
	using namespace moves;
	using namespace core;
	using namespace fragment;
	core::kinematics::MoveMapOP movemap = make_move_map(cenpose);
	MoverOP cendock = new protocols::symmetric_docking::SymDockingLowRes(sfc);
	// FragmentIO io;
	// FragSetOP hhh = io.read( basic::database::full_name("sampling/ss_fragfiles/HHH.fragfile") );
	// std::cerr << "frags " << hhh << std::endl;
	// MoverOP fragins( new protocols::abinitio::ClassicFragmentMover(hhh,movemap) );
	// MonteCarloOP mc( new MonteCarlo( cenpose, *sfc, 2000000.0 ) );
	// TrialMoverOP trials( new TrialMover( fragins, mc ) );
	//
	// for(Size i = 1; i <= 5; ++i) {
	// 	trials->keep_stats_type( all_stats );
	// 	for ( Size ii = 1; ii <= 100; ++ii ) {
	// 		trials->apply(cenpose);
	// 		mc->show_state();
	// 	}
	// 	mc->recover_low(cenpose);
	// 	std::cout << "frag accept rate " << trials->acceptance_rate() << std::endl;
	// }
	cendock->apply(cenpose);

}

/////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace core;
	using namespace scoring;
	using namespace protocols;
	using namespace moves;
	using namespace simple_moves::symmetry;
	using namespace core::pack::task;
	using namespace ObjexxFCL::format;
	using namespace id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::ostringstream sslog;

	devel::init(argc,argv);

	// read_frags(); std::exit(0);

	ScoreFunctionOP sfc,sfd,sf1,sf2,sf3,sf4,sf5;
	make_sf(sfc,sfd,sf1,sf2,sf3,sf4,sf5);

	Size NPACK=3,iter = 0;
	Real zcyl,censcore,censcore2,cstscore,znsscore;
	while(true) {
		iter++;
		//for( Size iter = 1; iter <= NSTRUCT; iter++ ) {
		core::pose::Pose cc,init,cenpose;
		std::string tag = string_of(numeric::random::uniform());
		CCParam p;
		cenpose = make_coiled_coil(p) ;

		zcyl = zcyl_score(cenpose,p.nres,p.trans); if( zcyl     <  1.7 ) { continue; }
		censcore = (*sfc)(cenpose)/p.nres;         if( censcore >= 1.0 ) {
			// TR << "censcore 1 fail " << censcore << std::endl;
			continue;
		}
		if( !add_symm_zns(cenpose,p,sfc,250.0) ) { continue; }

				// cenpose.dump_pdb("cen0.pdb");
		do_centroid_stuff(cenpose,sfc);
				// cenpose.dump_pdb("cen1.pdb");
		// std::exit(0);

		censcore2 = (*sfc)(cenpose);
		cstscore = sfc->get_weight( atom_pair_constraint) * cenpose.energies().total_energies()[ atom_pair_constraint]+
                 sfc->get_weight(     angle_constraint) * cenpose.energies().total_energies()[     angle_constraint]+
                 sfc->get_weight(coordinate_constraint) * cenpose.energies().total_energies()[coordinate_constraint]+
		           sfc->get_weight(  dihedral_constraint) * cenpose.energies().total_energies()[  dihedral_constraint];
		censcore2 = (censcore2-cstscore)/(p.nres+1);
		znsscore = cenpose.energies().residue_total_energy(p.nres+1);
		// TR << "cstscore " << cstscore << std::endl;
		//std::cerr << "CENSCORE AFTER MIN ZHC " << censcore << " " << censcore2 << " " << cstscore << std::endl;
		if( censcore2 >= 3.0 ) {
			TR << "ZHC censcore fail " << censcore2 << std::endl;
			continue;
		}

		TR << LJ(35,"cenp_cen_"+tag+".pdb.gz") << " " << F(10,5,censcore) << " " << F(10,5,censcore2) << " " << F(10,5,cstscore) << " " << F(10,5,znsscore) << " " << p.nres << " " << p.trans << " " << p.rot << std::endl;
		TR << p.str() << std::endl;

		cc = cenpose;
		cenpose = cc;

		cc.remove_constraints();
		protocols::toolbox::switch_to_residue_type_set( cc, core::chemical::FA_STANDARD );
		// cc.dump_pdb("fa0.pdb");
		// std::cerr << "fa0 " << (*sf2)(cc) << std::endl;
		add_fa_cst(cc,p);
		repack(cc,sf2);
		// std::cerr << "rpk " << (*sf2)(cc) << std::endl;
		// cc.dump_pdb("fa_rpk.pdb");
		minimize(cc,sf2);
		//		cc.dump_pdb("fa_min.pdb");
		// sf2->show(cc);
		// if( (*sf2)(cc) > 100 ) continue;
		init = cc;
		// continue;
		// std::exit(-1);

		core::pose::Pose best = cc;
		for(Size k = 1; k <= NPACK; ++k) {
			// design_AL(cc,sfd); minimize(cc,sfd); //sfd->show(cc);
			// design_AL(cc,sf1); minimize(cc,sf1); //sf1->show(cc);
			// design_AL(cc,sf2); minimize(cc,sf2); //sf2->show(cc);
			design_AL(cc,sf3); minimize(cc,sf3); //sf3->show(cc);
			design_AL(cc,sf4); minimize(cc,sf4); //sf4->show(cc);
			design_AL(cc,sf5); minimize(cc,sf5); //sf5->show(cc);
			TR << "design_AL/min12345    " << I(2,k) << " " << F(10,3,(*sf5)(cc)/p.nres) << " " << cc.sequence().substr(0,p.nres) << std::endl;
			if( (*sf5)(best) >= (*sf5)(cc) ) best = cc;
		}
		cc = best;
		for(Size k = 1; k <= NPACK; ++k) {
			// design_FILV(cc,sfd); minimize(cc,sfd);
			// design_FILV(cc,sf1); minimize(cc,sf1);
			// design_FILV(cc,sf2); minimize(cc,sf2);
			design_FILV(cc,sf3); minimize(cc,sf3);
			design_FILV(cc,sf4); minimize(cc,sf4);
			design_FILV(cc,sf5); minimize(cc,sf5);
			TR << "deisgnFILV/min12345   " << I(2,k) << " " << F(10,3,(*sf5)(cc)/p.nres) << " " << cc.sequence().substr(0,p.nres) << std::endl;
			if( (*sf5)(best) >= (*sf5)(cc) ) best = cc;
		}
		cc = best;
		for(Size k = 1; k <= NPACK; ++k) {
			// design_AFILV(cc,sfd); minimize(cc,sfd);
			// design_AFILV(cc,sf1); minimize(cc,sf1);
			// design_AFILV(cc,sf2); minimize(cc,sf2);
			design_AFILV(cc,sf3); minimize(cc,sf3);
			design_AFILV(cc,sf4); minimize(cc,sf4);
			design_AFILV(cc,sf5); minimize(cc,sf5);
			TR << "deisgnAFILV/min12345   " << I(2,k) << " " << F(10,3,(*sf5)(cc)/p.nres) << " " << cc.sequence().substr(0,p.nres) << std::endl;
			if( (*sf5)(best) >= (*sf5)(cc) ) best = cc;
		}
		cc = best;
		for(Size k = 1; k <= NPACK; ++k) {
			// design_AFILVEK(cc,sfd); minimize(cc,sfd);
			// design_AFILVEK(cc,sf1); minimize(cc,sf1);
			// design_AFILVEK(cc,sf2); minimize(cc,sf2);
			design_AFILVEK(cc,sf3); minimize(cc,sf3);
			design_AFILVEK(cc,sf4); minimize(cc,sf4);
			design_AFILVEK(cc,sf5); minimize(cc,sf5);
			TR << "deisgnAFILVEK/min12345   " << I(2,k) << " " << F(10,3,(*sf5)(cc)/p.nres) << " " << cc.sequence().substr(0,p.nres) << std::endl;
			if( (*sf5)(best) >= (*sf5)(cc) ) best = cc;
		}
		cc = best;

		Real score = (*sf5)(cc);
		Real ap  = sf5->get_weight(atom_pair_constraint)  * cc.energies().total_energies()[atom_pair_constraint ];
		Real ang = sf5->get_weight(angle_constraint)      * cc.energies().total_energies()[angle_constraint];
		Real dhc = sf5->get_weight(dihedral_constraint)   * cc.energies().total_energies()[dihedral_constraint];
		Real coc = sf5->get_weight(coordinate_constraint) * cc.energies().total_energies()[coordinate_constraint];
		score = score - ap - ang - coc - dhc;
		//		sf5->show(cc);

		sslog << LJ(25,("cc"+tag+".pdb").c_str()) << " fa: " << F(10,5,score) << " " << F(10,5,ap) << " " << F(10,5,ang) << " " << F(10,5,score/p.nres) << " " << F(10,5,score/natoms(cc,p.nres))
		          << " zcyl " << F(8,5,zcyl_score(cc,p.nres,p.trans))
		          << " srms " << F(8,5,CA_rmsd(cc,cenpose)) << " " << I(7,iter) << " "
					 << cc.sequence().substr(0,p.nres) << " 		";
		sslog << p.str() << std::endl;

		if( basic::options::option[in::file::silent_energy_cut]() > score/p.nres ) {
			cenpose.dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/cc"+tag+"_cen.pdb.gz");
			cc     .dump_pdb(std::string(option[OptionKeys::out::file::o]())+"/cc"+tag+"_fa.pdb.gz");
			std::cout << sslog.str();
			std::cout.flush();
			sslog.clear();
			sslog.str("");
		}

	}

	std::cout << sslog.str();
	std::cout.flush();
	sslog.clear();
	sslog.str("");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
