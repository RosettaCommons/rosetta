#include <core/scoring/motif/motif_hash_stuff.hh>

	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <basic/Tracer.hh>
	#include <basic/database/open.hh>
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/mh.OptionKeys.gen.hh>
	#include <basic/options/option_macros.hh>
	#include <basic/pymol_chains.hh>
	#include <core/chemical/AtomType.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/conformation/symmetry/util.hh>
	#include <core/conformation/Residue.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/motif/reference_frames.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/util.hh>
	#include <core/pose/symmetry/util.hh>
	#include <core/conformation/symmetry/SymmetryInfo.hh>
	#include <core/pose/xyzStripeHashPose.hh>
	#include <core/io/pdb/pose_io.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/EnergyGraph.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoreTypeManager.hh>
	#include <core/scoring/dssp/Dssp.hh>
	#include <core/scoring/dssp/StrandPairing.hh>
	#include <core/scoring/hbonds/HBondOptions.hh>
	#include <core/scoring/methods/EnergyMethodOptions.hh>
	#include <core/scoring/packing/compute_holes_score.hh>
	#include <core/scoring/rms_util.hh>
	#include <core/scoring/sasa.hh>
	#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/file/file_sys_util.hh>
	#include <utility/fixedsizearray1.hh>
	#include <numeric/xyzTransform.hh>

	#include <numeric/geometry/hashing/SixDHasher.hh>

	#include <boost/unordered_set.hpp>

	#include <boost/foreach.hpp>
	#include <bitset>
    #ifndef _WIN32  
	#include <pthread.h>
    #endif
   


	// #include <core/pack/task/PackerTask.hh>
	// #include <core/pack/packer_neighbors.hh>
	// #include <core/pack/rotamer_set/RotamerSet.hh>
	// #include <core/pack/rotamer_set/RotamerSetFactory.hh>
	// #include <core/pack/task/TaskFactory.hh>

	#define MAX_UINT16 65535
	#define MAX_UINT8    255
	#define XFORM_SCORE_FILE_VERSION 1

namespace core {
namespace scoring {
namespace motif {

	using numeric::Xforms;


/************************************************* types ************************************************/
	static basic::Tracer TR("core.scoring.motif.util");

	using core::pose::PoseCoordPickMode_BB;
	using core::pose::PoseCoordPickMode_N_C_O;

	typedef utility::fixedsizearray1<float,20> float20;
	typedef utility::fixedsizearray1<float,9> float9;
	using std::make_pair;
	using core::chemical::AA;
	using core::id::AtomID;
	using basic::options::option;
	namespace mh = basic::options::OptionKeys::mh;
	using core::pose::Pose;
	using core::Real;
	using core::scoring::ScoreFunctionOP;
	using core::Size;
	using numeric::max;
	using numeric::min;
	using numeric::random::gaussian;
	using numeric::random::uniform;
	using numeric::rotation_matrix_degrees;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;
	using namespace ObjexxFCL::format;
	using ObjexxFCL::string_of;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::ostream;
	using std::string;
	using utility::io::izstream;
	using utility::io::ozstream;
	using utility::file_basename;
	using utility::vector1;
	using std::endl;
	using core::import_pose::pose_from_pdb;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	using core::pose::xyzStripeHashPoseCOP;
	using core::pose::initialize_atomid_map;
	using core::pose::motif::get_nterminal_peptide_bond_reference_frame;
	using core::pose::motif::get_cterminal_peptide_bond_reference_frame;
	using core::pose::motif::get_backbone_reference_frame;
	using core::pose::motif::get_sidechain_reference_frame;
	using core::pose::motif::get_nterminal_peptide_bond_reference_frame_atomids;
	using core::pose::motif::get_cterminal_peptide_bond_reference_frame_atomids;
	using core::pose::motif::get_backbone_reference_frame_atomids;
	using core::pose::motif::get_sidechain_reference_frame_atomids;
	using core::pose::motif::get_sidechain_reference_frame_atomids_with_downstream;
	using core::pose::motif::get_backbone_reference_frame_atomids_with_downstream;

	using core::chemical::aa_gly;

void xform_pose( core::pose::Pose & pose, Xform const & s, Size sres=1, Size eres=0 ) {
  if(eres==0) eres = pose.n_residue();
  for(Size ir = sres; ir <= eres; ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
 }
Mat random_rotation(){ // from quaternion
	Real a = uniform() * 2.0 - 1.0; Real b = uniform() * 2.0 - 1.0;
	Real c = uniform() * 2.0 - 1.0; Real d = uniform() * 2.0 - 1.0;
	while( a*a + b*b + c*c + d*d > 1.0 ){
		a=uniform()*2.0-1.0; b=uniform()*2.0-1.0;
		c=uniform()*2.0-1.0; d=uniform()*2.0-1.0;
	 }
	Real len = sqrt(a*a+b*b+c*c+d*d);
	a /= len; b /= len; c /= len; d /= len;
	return Mat::cols(
		a*a + b*b - c*c - d*d,        2.0*(b*c - a*d),         2.0*(b*d + a*c),
				2.0*(b*c + a*d),      a*a - b*b + c*c - d*d,       2.0*(c*d - a*b),
				2.0*(b*d - a*c),          2.0*(c*d + a*b),     a*a - b*b - c*c + d*d
	);
 }

std::string tag_from_pdb_fname(string const & fname0){
	string fname = utility::file_basename(fname0);
	if     (fname.substr(fname.size()-3,3)==".gz" ) fname = fname.substr(0,fname.size()-3);
	if     (fname.substr(fname.size()-4,4)==".pdb") fname = fname.substr(0,fname.size()-4);
	else if(fname.substr(fname.size()-5,4)==".pdb") fname = fname.substr(0,fname.size()-5) + fname[fname.size()-1];
	else utility_exit_with_message("file is not .pdb or .pdb#");
	if(fname.size()==4) fname = fname+"_";
	runtime_assert(fname.size()==5);
	return fname;
 }

ostream & operator<<(ostream & out, Real6 const & r6){
	out << r6[1] <<' ' << r6[2] <<' ' << r6[3] <<' ' << r6[4] <<' ' << r6[5] <<' ' << r6[6];
	return out;
 }

Real rt6_rt6_bb_dis2_explicit_stupid(
	Real6 const & x1,
	Real6 const & x2
){
	utility_exit_with_message("being refactored");
	Xform ht1,ht2;
	ht1.from_euler_angles_deg(numeric::xyzVector<Real>(x1[4],x1[5],x1[6]));
	ht2.from_euler_angles_deg(numeric::xyzVector<Real>(x2[4],x2[5],x2[6]));
	ht1.t = (Vec(x1[1],x1[2],x1[3]));
	ht2.t = (Vec(x2[1],x2[2],x2[3]));
	Vec const   N(      1.454298456301180 ,     0.000000000000000 ,     0.000000000000000 );
	Vec const  CA(     -0.000000000000000 ,     0.000000000000000 ,     0.000000000000000 );
	Vec const   C(    -0.5232708573008347 ,     1.456155420928565 ,     0.000000000000000 );
	Vec const   O(   -0.04476385140748412 ,     2.289714637034945 ,    -0.762673638298255 );
	Vec const  CB(    -0.5411598950614691 ,    -0.849773440692597 ,     1.16326783995368  );
	Vec const CG1(    -0.6475754656270934 ,    -2.320037628202743 ,     0.732861119204208 );
	Vec const CG2(     -1.882609438342826 ,    -0.327829466486401 ,     1.63164013911332  );
	Vec const CD1(     -1.738720129313219 ,    -3.190855526746606 ,     1.428714918852484 );
	Real d2 = 0.0;
	d2 += (ht1*  N).distance_squared(ht2*  N);
	d2 += (ht1* CA).distance_squared(ht2* CA);
	d2 += (ht1*  C).distance_squared(ht2*  C);
	d2 += (ht1*  O).distance_squared(ht2*  O);
	d2 += (ht1* CB).distance_squared(ht2* CB);
	d2 += (ht1*CG1).distance_squared(ht2*CG1);
	d2 += (ht1*CG2).distance_squared(ht2*CG2);
	d2 += (ht1*CD1).distance_squared(ht2*CD1);
	d2 /= 8.0;
	return d2;
 }

Real rt6_rt6_dis2(Real6 const & x1, Real6 const & x2, Real const & lever){
	utility_exit_with_message("being refactored");
	Xform ht1,ht2;
	ht1.from_euler_angles_deg(numeric::xyzVector<Real>(x1[4],x1[5],x1[6]));
	ht2.from_euler_angles_deg(numeric::xyzVector<Real>(x2[4],x2[5],x2[6]));
	Xform r = ht1 * ht2.inverse();
	Real cos_theta = (r.xx()+r.yy()+r.zz()-1.0)/2.0;
	Real sin2_theta = 1.0-cos_theta*cos_theta;
	Real dis2 = Vec(x1[1],x1[2],x1[3]).distance_squared(Vec(x2[1],x2[2],x2[3])) + lever*lever*sin2_theta;
	return dis2;
 }
Real rt6_rt6_bb_dis2(Real6 const & x1, Real6 const & x2){
	utility_exit_with_message("being refactored");
	return rt6_rt6_dis2(x1,x2,std::sqrt(6.0));
 }

Real6 inverse_rt6(Real6 const & rt){
	return Xform(rt).inverse().rt6();
 }

Real6
rt_to_real6(core::kinematics::RT const & rt){
	utility_exit_with_message("being refactored");
	Xform ht( rt.get_rotation() , rt.get_translation() );
	numeric::xyzVector < Real > euler_angles =  ht.euler_angles_deg();
	Real6 rt6;
	rt6[1] = rt.get_translation().x();
	rt6[2] = rt.get_translation().y();
	rt6[3] = rt.get_translation().z();
	rt6[4] = fmod(euler_angles.x(),360.0);
	rt6[5] = fmod(euler_angles.y(),360.0);
	rt6[6] = fmod(euler_angles.z(),360.0);
	rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
	rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
	rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
	return rt6;
 }
core::kinematics::RT
real6_to_rt(Real6 const & rt6){
	utility_exit_with_message("being refactored");
	Xform x;
	x.from_euler_angles_deg(numeric::xyzVector<Real>(rt6[4],rt6[5],rt6[6]));
	numeric::xyzVector<Real> t(rt6[1],rt6[2],rt6[3]);
	return core::kinematics::RT(x.R,t);
 }

	RM_Type rpm_type1(RPM_Type const & type) {
		switch(type){
			case SC_SC: return RM_SC;
			case SC_BB: return RM_SC;
			case SC_PH: return RM_SC;
			case SC_PO: return RM_SC;
			case BB_BB: return RM_BB;
			case BB_PH: return RM_BB;
			case BB_PO: return RM_BB;
			case PH_PO: return RM_PH;
			case RPM_Type_NONE: utility_exit_with_message("RPM_Type_NONE not valid!");
			default: utility_exit_with_message("arst");
		}
	}
	RM_Type rpm_type2(RPM_Type const & type) {
		switch(type){
			case SC_SC: return RM_SC;
			case SC_BB: return RM_BB;
			case SC_PH: return RM_PH;
			case SC_PO: return RM_PO;
			case BB_BB: return RM_BB;
			case BB_PH: return RM_PH;
			case BB_PO: return RM_PO;
			case PH_PO: return RM_PO;
			case RPM_Type_NONE: utility_exit_with_message("RPM_Type_NONE not valid!");
			default: utility_exit_with_message("arst");
		}
	}

Xform get_residue_pair_xform(Pose const & pose1, Size ir, Pose const & pose2, Size jr, RPM_Type const & type){
	// Xform const bbr1 = core::pose::motif::get_backbone_reference_frame(pose1,ir);
	// Xform const bbr2 = core::pose::motif::get_backbone_reference_frame(pose2,jr);
	// return (~bbr1*bbr2).rt6();
	// //!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Xform frame1,frame2;
	switch(rpm_type1(type)){
		case RM_SC: frame1 =              get_sidechain_reference_frame(pose1,ir); break;
		case RM_BB: frame1 =               get_backbone_reference_frame(pose1,ir); break;
		case RM_PH: frame1 = get_nterminal_peptide_bond_reference_frame(pose1,ir); break;
		case RM_PO: frame1 = get_cterminal_peptide_bond_reference_frame(pose1,ir); break;
		default: utility_exit_with_message("unknown motif type");
	}
	switch(rpm_type2(type)){
		case RM_SC: frame2 =              get_sidechain_reference_frame(pose2,jr); break;
		case RM_BB: frame2 =               get_backbone_reference_frame(pose2,jr); break;
		case RM_PH: frame2 = get_nterminal_peptide_bond_reference_frame(pose2,jr); break;
		case RM_PO: frame2 = get_cterminal_peptide_bond_reference_frame(pose2,jr); break;
		default: utility_exit_with_message("unknown motif type");
	}
	if( frame1.bad() || frame2.bad() ) return Xform::BAD_XFORM();
	return ~frame1*frame2;
 }
Xform get_residue_pair_xform(Pose const & pose, Size ir,  Size jr, RPM_Type const & type){
	return get_residue_pair_xform(pose,ir,pose,jr,type);
 }
Real6 get_residue_pair_rt6(Pose const & pose1, Size ir, Pose const & pose2, Size jr, RPM_Type const & type){
	return get_residue_pair_xform(pose1,ir,pose2,jr,type).rt6();
}
Real6 get_residue_pair_rt6(Pose const & pose, Size ir, Size jr, RPM_Type const & type){
	return get_residue_pair_rt6(pose,ir,pose,jr,type);
 }

void set_residue_pair_xform(Xform const & x, Pose & pose, Size ir, Size jr, RPM_Type const & type){
	Xform frame1,frame2;
	switch(rpm_type1(type)){
		case RM_SC: frame1 =              get_sidechain_reference_frame(pose,ir); break;
		case RM_BB: frame1 =               get_backbone_reference_frame(pose,ir); break;
		case RM_PH: frame1 = get_nterminal_peptide_bond_reference_frame(pose,ir); break;
		case RM_PO: frame1 = get_cterminal_peptide_bond_reference_frame(pose,ir); break;
		default: utility_exit_with_message("unknown motif type");
	}
	switch(rpm_type2(type)){
		case RM_SC: frame2 =              get_sidechain_reference_frame(pose,jr); break;
		case RM_BB: frame2 =               get_backbone_reference_frame(pose,jr); break;
		case RM_PH: frame2 = get_nterminal_peptide_bond_reference_frame(pose,jr); break;
		case RM_PO: frame2 = get_cterminal_peptide_bond_reference_frame(pose,jr); break;
		default: utility_exit_with_message("unknown motif type");
	}
	xform_pose(pose,  ~frame1,1,1);
	xform_pose(pose,x*~frame2,2,2);

	// core::id::StubID id1( AtomID(1,ir), AtomID(2,ir), AtomID(3,ir) );
	// core::id::StubID id2( AtomID(1,jr), AtomID(2,jr), AtomID(3,jr) );
	// Xform const xref = get_residue_pair_xform(pose,ir,jr,type);
	// Xform const delta = x * ~xref;
	// core::kinematics::RT rt0 = pose.conformation().get_stub_transform(id1,id2);
	// Xform const x0(rt0.get_rotation(),rt0.get_translation());
	// Xform xmov = delta*x0;
	// core::kinematics::RT rt(xmov.R,xmov.t);
	// pose.conformation().set_stub_transform(id1,id2,rt);
 }

Reals get_sasa(Pose const & pose, Real const & probesize){
	core::id::AtomID_Map<Real> atom_sasa;
	Reals rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probesize, false );
	if(rsd_sasa.size() != pose.n_residue()) utility_exit_with_message("bad sasa!");
	return rsd_sasa;
 }
Reals get_nbrs(Pose const & pose){
	Reals nbrs(pose.n_residue(),0);
	for(Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(!pose.residue(ir).is_protein()) continue;
		Vec CB1 = pose.residue(ir).xyz(2);
		for(Size jr = 1; jr <= pose.n_residue(); ++jr){
			if(!pose.residue(jr).is_protein()) continue;
			if(CB1.distance_squared(pose.residue(jr).xyz(2)) < 100.0){
				nbrs[ir] += 1.0;
			}
		}
	 }
	return nbrs;
 }
// Real get_sc_sasa(Pose const & pose, Size const & ir){
	// 	core::id::AtomID_Map<Real> atom_sasa;
	// 	Reals rsd_sasa;
	// 	core::id::AtomID_Map<bool> atom_subset;
	// 	initialize_atomid_map( atom_subset, pose, false );
	// 	for(Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia) atom_subset[AtomID(ia,ir)] = true;
	// 	return core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 1.4, false, atom_subset );
	// }
	// Real get_sc_sasa(Pose const & pose, Size const & ir, Size const & jr){
	// 	core::id::AtomID_Map<Real> atom_sasa;
	// 	Reals rsd_sasa;
	// 	core::id::AtomID_Map<bool> atom_subset;
	// 	initialize_atomid_map( atom_subset, pose, false );
	// 	for(Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia) atom_subset[AtomID(ia,ir)] = true;
	// 	for(Size ia = 5; ia <= pose.residue(jr).nheavyatoms(); ++ia) atom_subset[AtomID(ia,jr)] = true;
	// 	return core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 1.4, false, atom_subset );
	// }
	// Real get_sc_sasa_2rsd(Pose const & pose, Size ir, Size jr){
	// 	Pose tmp;
	// 	tmp.append_residue_by_jump(pose.residue(ir),1);
	// 	tmp.append_residue_by_jump(pose.residue(jr),1);
	// 	Real sasa0 = get_sc_sasa(tmp,1,2);
	// 	Real sasa1 = get_sc_sasa(tmp,1);
	// 	Real sasa2 = get_sc_sasa(tmp,2);
	// 	Real sciface = sasa1+sasa2-sasa0;
	// 	// if( pose.residue(ir).aa() == core::chemical::aa_ser ){
	// 	//  cerr << "get_sc_sasa_2rsd " << ir << " " << jr << " " << sasa0 <<' '<< sasa1 <<' '<< sasa2 <<' '<< sasa1+sasa2-sasa0 << endl;
	// 	// }
	// 	// tmp.dump_pdb("sasa_test.pdb");
	// 	// utility_exit_with_message("sasa test");
	// 	return sciface;
	// }

void HACK_dump_helix(Pose const & pose, string fname, int beg, int end){
	cout << "dump helix to " << fname << endl;
	utility::io::ozstream out("helices/"+fname);
	Size ano = 1;
	for(int i = beg; i <= end; ++i){
		core::io::pdb::dump_pdb_residue(pose.residue(i),ano,out);
	 }
	out.close();

 }
int HACK_dump_helices(Pose const & pose, string tag, int nres, int minlen=10){
	int hbeg = 'H'==pose.secstruct(1);
	int hend = 'H'==pose.secstruct(1);
	int nhelix = 0;
	for(int ir = 2; ir <= nres; ++ir){
		cout << "HACK_dump_helix " << ir << " " << pose.secstruct(ir) << " " << hbeg << " " << hend << " " << tag << endl;
		bool newc =	pose.chain(ir-1)!=pose.chain(ir);
		if(newc){
			if(hbeg&&hend && hend-hbeg >= minlen-1) ++nhelix;
			if(hbeg&&hend && hend-hbeg >= minlen-1) HACK_dump_helix(pose,tag+((nhelix>8)?"_":"__")+string_of(nhelix)+".pdb",hbeg,hend);
			hbeg = hend = 'H'==pose.secstruct(ir) ? ir : 0;
		} else {
			bool newh = 'H'!=pose.secstruct(ir-1) && 'H'==pose.secstruct(ir);
			bool endh = 'H'==pose.secstruct(ir-1) && 'H'!=pose.secstruct(ir);
			if('H'==pose.secstruct(ir) && hbeg) hend = ir;
			if(newh) hbeg = hend = ir;
			if(endh && hend-hbeg >= minlen-1) ++nhelix;
			if(endh && hend-hbeg >= minlen-1) HACK_dump_helix(pose,tag+(nhelix>8?"_":"__")+string_of(nhelix)+".pdb",hbeg,hend);
			if(endh) hbeg = hend = 0;
		}
	 }
	return nhelix;
 }

Real6 get_bins(Real c, Real a){
	utility::fixedsizearray1<Real,6> bins(c);
	bins[4] = a;
	bins[5] = a;
	bins[6] = a;
	return bins;
 }



core::id::AtomID_Mask get_motif_atom_mask( Pose const & motif_pose, RPM_Type const & type, bool with_Hpol ){
	core::id::AtomID_Mask mask; {
		core::pose::initialize_atomid_map(mask,motif_pose,false);
		switch(rpm_type1(type)){
			case RM_PH:
				BOOST_FOREACH(AtomID const & aid, get_nterminal_peptide_bond_reference_frame_atomids(motif_pose,1) )  mask[aid] = true;
				break;
			case RM_PO:
				BOOST_FOREACH(AtomID const & aid, get_cterminal_peptide_bond_reference_frame_atomids(motif_pose,1) )  mask[aid] = true;
				break;
			case RM_BB:
				BOOST_FOREACH(AtomID const & aid, get_backbone_reference_frame_atomids_with_downstream(motif_pose,1) )  mask[aid] = true;
				break;
			case RM_SC:
				if(with_Hpol) BOOST_FOREACH(AtomID const & aid, get_sidechain_reference_frame_atomids_with_downstream(motif_pose,1) )  mask[aid] = true;
				else          BOOST_FOREACH(AtomID const & aid, get_sidechain_reference_frame_atomids                (motif_pose,1) )  mask[aid] = true;
				break;
			default: utility_exit_with_message("unknown motif type");
		}
		switch(rpm_type2(type)){
			case RM_PH:
				BOOST_FOREACH(AtomID const & aid, get_nterminal_peptide_bond_reference_frame_atomids(motif_pose,2) )  mask[aid] = true;
				break;
			case RM_PO:
				BOOST_FOREACH(AtomID const & aid, get_cterminal_peptide_bond_reference_frame_atomids(motif_pose,2) )  mask[aid] = true;
				break;
			case RM_BB:
				BOOST_FOREACH(AtomID const & aid, get_backbone_reference_frame_atomids_with_downstream(motif_pose,2) )  mask[aid] = true;
				break;
			case RM_SC:
				if(with_Hpol) BOOST_FOREACH(AtomID const & aid, get_sidechain_reference_frame_atomids_with_downstream(motif_pose,2) )  mask[aid] = true;
				else          BOOST_FOREACH(AtomID const & aid, get_sidechain_reference_frame_atomids                (motif_pose,2) )  mask[aid] = true;
				break;
			default: utility_exit_with_message("unknown motif type");
		}
	}
	return mask;
 }

Real align_motif_pose_NCAC_super( Pose & motif_pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type ){
	using namespace core::id;
	if(type!=BB_BB) utility_exit_with_message("not impl, NCAC_super and not BB_BB");
	// using rms super
	Pose ref;
	ref.append_residue_by_jump(paln1.residue(ir),1);
	ref.append_residue_by_jump(paln2.residue(jr),1);
	AtomID_Map<AtomID> atommap;
	initialize_atomid_map(atommap,motif_pose,BOGUS_ATOM_ID);
	atommap[AtomID(1,1)] = AtomID(1,1);
	atommap[AtomID(2,1)] = AtomID(2,1);
	atommap[AtomID(3,1)] = AtomID(3,1);
	atommap[AtomID(1,2)] = AtomID(1,2);
	atommap[AtomID(2,2)] = AtomID(2,2);
	atommap[AtomID(3,2)] = AtomID(3,2);
	core::scoring::superimpose_pose(motif_pose,ref,atommap);
	Real myrms = sqrt( ref.xyz(AtomID(1,1)).distance_squared(motif_pose.xyz(AtomID(1,1))) +
	                   ref.xyz(AtomID(2,1)).distance_squared(motif_pose.xyz(AtomID(2,1))) +
	                   ref.xyz(AtomID(3,1)).distance_squared(motif_pose.xyz(AtomID(3,1))) +
	                   ref.xyz(AtomID(1,2)).distance_squared(motif_pose.xyz(AtomID(1,2))) +
	                   ref.xyz(AtomID(2,2)).distance_squared(motif_pose.xyz(AtomID(2,2))) +
	                   ref.xyz(AtomID(3,2)).distance_squared(motif_pose.xyz(AtomID(3,2))) / 6.0 );
	return myrms;
 }

Real align_motif_pose_break( Pose & /*motif_pose*/, Pose const & /*paln1*/, Size const & /*ir*/, Pose const & /*paln2*/, Size const & /*jr*/, RPM_Type const & /*type*/ ){
	utility_exit_with_message("NOT IMPL");
    using namespace core::id;
    return(1);//python build was complaining.
	// if(type!=BB_BB) utility_exit_with_message("not impl");
	// Pose start(motif_pose);
	// //align separately to bb's
	// Xform const stub1  = get(paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(1,ir)),paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(3,ir)));
	// Xform const stub2  = get(paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(1,jr)),paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(3,jr)));
	// Xform const mstub1 = Xform().from_four_points(motif_pose.xyz(AtomID(2,1)),motif_pose.xyz(AtomID(1,1)),motif_pose.xyz(AtomID(2,1)),motif_pose.xyz(AtomID(3,1)));
	// Xform const mstub2 = Xform().from_four_points(motif_pose.xyz(AtomID(2,2)),motif_pose.xyz(AtomID(1,2)),motif_pose.xyz(AtomID(2,2)),motif_pose.xyz(AtomID(3,2)));
	// xform_pose(motif_pose,stub1*~mstub1,1,1);
	// xform_pose(motif_pose,stub2*~mstub2,2,2);

	// return core::scoring::bb_rmsd(motif_pose,start);

 }

Real align_motif_pose_by_one_frame( Pose & /*motif_pose*/, Pose const & /*paln1*/, Size const & /*ir*/, Pose const & /*paln2*/, Size const & /*jr*/, RPM_Type const & /*type*/ ){
	utility_exit_with_message("align_motif_pose_by_one_frame not implemented");
	using namespace core::id;
 	// Size mres1=0,mres2=0;
	// get_motif_resnums(motif_pose,mres1,mres2,type);

	// Xform frame1, frame2, mframe1, mframe2;
	// switch(type){
	// 	case RM_PO: mframe1 = get_cterminal_peptide_bond_reference_frame(motif_pose,mres1); break;
	// 	case RM_BB: mframe1 =               get_backbone_reference_frame(motif_pose,mres1); break;
	// 	case RM_SC: mframe1 =              get_sidechain_reference_frame(motif_pose,mres1); break;
	// 	default: utility_exit_with_message("unknown motif type");
	// }
	// switch(type){
	// 	case RM_PO: mframe2 = get_cterminal_peptide_bond_reference_frame(motif_pose,mres2); break;
	// 	case RM_BB: mframe2 =               get_backbone_reference_frame(motif_pose,mres2); break;
	// 	case RM_SC: mframe2 =              get_sidechain_reference_frame(motif_pose,mres2); break;
	// 	default: utility_exit_with_message("unknown motif type");
	// }
	// switch(type){
	// 	case RM_PO: frame1 = get_cterminal_peptide_bond_reference_frame(paln1,ir); break;
	// 	case RM_BB: frame1 =               get_backbone_reference_frame(paln1,ir); break;
	// 	case RM_SC: frame1 =              get_sidechain_reference_frame(paln1,ir); break;
	// 	default: utility_exit_with_message("unknown motif type");
	// }
	// switch(type){
	// 	case RM_PO: frame2 = get_cterminal_peptide_bond_reference_frame(paln2,jr); break;
	// 	case RM_BB: frame2 =               get_backbone_reference_frame(paln2,jr); break;
	// 	case RM_SC: frame2 =              get_sidechain_reference_frame(paln2,jr); break;
	// 	default: utility_exit_with_message("unknown motif type");
	// }

	// // motif_pose.dump_pdb("test0.pdb");

	// switch(type){
	// 	case RM_BB: case PB_SC: xform_pose(motif_pose,frame2*~mframe2); break;
	// 	default:                                        xform_pose(motif_pose,frame1*~mframe1); break;
	// }

	// // motif_pose.dump_pdb("test1.pdb");
	// // utility_exit_with_message("aoirts");

	// //align separately to bb's
	// // xform_pose(motif_pose,frame1*~mframe1);
	// // xform_pose(motif_pose,frame1*~mframe1,1,1);
	// // xform_pose(motif_pose,frame2*~mframe2,2,2);

	return -1;

 }

Real align_motif_pose_super( Pose & motif_pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type ){
	using namespace core::id;
	vector1<AtomID> modids1,refids1,modids2,refids2;
	Pose refer_pose;
	AtomID_Map<AtomID> atommap;
	initialize_atomid_map(atommap,motif_pose,BOGUS_ATOM_ID);
	refer_pose.append_residue_by_jump(paln1.residue(ir),1);
	switch(rpm_type1(type)){
		case RM_PH:
			modids1 = get_nterminal_peptide_bond_reference_frame_atomids(motif_pose,1);
			refids1 = get_nterminal_peptide_bond_reference_frame_atomids(refer_pose,1);
			runtime_assert(modids1.size()==refids1.size());
			for(Size i=1; i <= modids1.size(); ++i)
				if(modids1[i]!=BOGUS_ATOM_ID)
					atommap[modids1[i]] = refids1[i];
			break;
		case RM_PO:
			modids1 = get_cterminal_peptide_bond_reference_frame_atomids(motif_pose,1);
			refids1 = get_cterminal_peptide_bond_reference_frame_atomids(refer_pose,1);
			runtime_assert(modids1.size()==refids1.size());
			for(Size i=1; i <= modids1.size(); ++i)
				if(modids1[i]!=BOGUS_ATOM_ID)
					atommap[modids1[i]] = refids1[i];
			break;
		case RM_BB:
			// cout << "case1 RM_BB" << endl;
			modids1 = get_backbone_reference_frame_atomids(motif_pose,1);
			refids1 = get_backbone_reference_frame_atomids(refer_pose,1);
			runtime_assert(modids1.size()==refids1.size());
			for(Size i=1; i <= modids1.size(); ++i)
				if(modids1[i]!=BOGUS_ATOM_ID && refids1[i]!=BOGUS_ATOM_ID)
					atommap[modids1[i]] = refids1[i];
			break;
		case RM_SC:
			// cout << "case1 RM_SC" << endl;
			modids1 = get_sidechain_reference_frame_atomids(motif_pose,1);
			refids1 = get_sidechain_reference_frame_atomids(refer_pose,1);
			runtime_assert(modids1.size()==refids1.size());
			for(Size i=1; i <= modids1.size(); ++i)
				if(modids1[i]!=BOGUS_ATOM_ID)
					atommap[modids1[i]] = refids1[i];
			break;
		default: utility_exit_with_message("unknown motif type");
	}
	refer_pose.append_residue_by_jump(paln2.residue(jr),1);
	switch(rpm_type2(type)){
		case RM_PH:
			modids2 = get_nterminal_peptide_bond_reference_frame_atomids(motif_pose,2);
			refids2 = get_nterminal_peptide_bond_reference_frame_atomids(refer_pose,2);
			runtime_assert(modids2.size()==refids2.size());
			for(Size i=1; i <= modids2.size(); ++i)
				if(modids2[i]!=BOGUS_ATOM_ID)
					atommap[modids2[i]] = refids2[i];
			break;
		case RM_PO:
			modids2 = get_cterminal_peptide_bond_reference_frame_atomids(motif_pose,2);
			refids2 = get_cterminal_peptide_bond_reference_frame_atomids(refer_pose,2);
			runtime_assert(modids2.size()==refids2.size());
			for(Size i=1; i <= modids2.size(); ++i)
				if(modids2[i]!=BOGUS_ATOM_ID)
					atommap[modids2[i]] = refids2[i];
			break;
		case RM_BB:
			// cout << "case2 RM_BB" << endl;
			modids2 = get_backbone_reference_frame_atomids(motif_pose,2);
			refids2 = get_backbone_reference_frame_atomids(refer_pose,2);
			runtime_assert(modids2.size()==refids2.size());
			for(Size i=1; i <= modids2.size(); ++i)
				if(modids2[i]!=BOGUS_ATOM_ID && refids2[i]!=BOGUS_ATOM_ID)
					atommap[modids2[i]] = refids2[i];
			break;
		case RM_SC:
			// cout << "case2 RM_SC" << endl;
			modids2 = get_sidechain_reference_frame_atomids(motif_pose,2);
			refids2 = get_sidechain_reference_frame_atomids(refer_pose,2);
			runtime_assert(modids2.size()==refids2.size());
			for(Size i=1; i <= modids2.size(); ++i)
				if(modids2[i]!=BOGUS_ATOM_ID)
					atommap[modids2[i]] = refids2[i];
			break;
		default: utility_exit_with_message("unknown motif type");
	}
	core::scoring::superimpose_pose(motif_pose,refer_pose,atommap);
	Real rms = 0.0;
	for(Size i=1; i <= modids1.size(); ++i) if(modids1[i]!=BOGUS_ATOM_ID&&refids1[i]!=BOGUS_ATOM_ID) rms += motif_pose.xyz(modids1[i]).distance_squared(refer_pose.xyz(refids1[i]));
	for(Size i=1; i <= modids2.size(); ++i) if(modids2[i]!=BOGUS_ATOM_ID&&refids2[i]!=BOGUS_ATOM_ID) rms += motif_pose.xyz(modids2[i]).distance_squared(refer_pose.xyz(refids2[i]));
	if(rpm_type1(type)==RM_SC) rms *= 0.8;
	if(rpm_type2(type)==RM_SC) rms *= 0.8;
	// if(rms > 0.1){
	// 	for(Size i=1; i <= modids1.size(); ++i)
	// 		if(modids1[i]!=BOGUS_ATOM_ID)
	// 			cout << motif_pose.residue(modids1[i].rsd()).name() << "/" << motif_pose.residue(modids1[i].rsd()).atom_name(modids1[i].atomno()) << " -- "
	// 			     << refer_pose.residue(refids1[i].rsd()).name() << "/" << refer_pose.residue(refids1[i].rsd()).atom_name(refids1[i].atomno()) << endl;
	// 	for(Size i=1; i <= modids2.size(); ++i)
	// 		if(modids2[i]!=BOGUS_ATOM_ID)
	// 			cout << motif_pose.residue(modids2[i].rsd()).name() << "/" << motif_pose.residue(modids2[i].rsd()).atom_name(modids2[i].atomno()) << " -- "
	// 			     << refer_pose.residue(refids2[i].rsd()).name() << "/" << refer_pose.residue(refids2[i].rsd()).atom_name(refids2[i].atomno()) << endl;
	// 	cout << endl;
	// }
	return rms;
 }

Real align_motif_pose( Pose & motif_pose, Pose const & paln1, Size const & ir, Pose const & paln2, Size const & jr, RPM_Type const & type ){
	// return align_motif_pose_by_one_frame(motif_pose,paln1,ir,paln2,jr,type);
	return align_motif_pose_super(motif_pose,paln1,ir,paln2,jr,type);
 }



core::Real aa_trustworthiness(char aa){
	switch(aa){
		case 'A': return 1.0;
		case 'C': return 0.8;
		case 'D': return 0.8;
		case 'E': return 0.6;
		case 'F': return 0.7;
		case 'G': return 1.0;
		case 'H': return 0.7;
		case 'I': return 0.9;
		case 'K': return 0.5;
		case 'L': return 0.9;
		case 'M': return 0.7;
		case 'N': return 0.8;
		case 'P': return 1.0;
		case 'Q': return 0.6;
		case 'R': return 0.5;
		case 'S': return 0.8;
		case 'T': return 0.8;
		case 'V': return 0.9;
		case 'W': return 0.6;
		case 'Y': return 0.7;
		default: return 1.0;
	}
	return 1.0;
}


}}}

