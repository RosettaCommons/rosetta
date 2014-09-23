// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:

// headers
	#define MAX_UINT16 65535
	#define MAX_UINT8    255

	#include <protocols/motif_hash/motif_hash_stuff.hh>

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
	#include <core/conformation/ResidueFactory.hh>
	#include <core/import_pose/import_pose.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pack/rotamer_set/RotamerSet.hh>
	#include <core/pack/task/TaskFactory.hh>
	#include <core/pose/Pose.fwd.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/util.hh>
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
	#include <devel/init.hh>
	#include <numeric/conversions.hh>
	#include <numeric/model_quality/rms.hh>
	#include <numeric/random/random.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <numeric/xyzVector.hh>
	#include <protocols/sic_dock/RigidScore.hh>
	#include <protocols/sic_dock/SICFast.hh>
	#include <protocols/sic_dock/xyzStripeHashPose.hh>
	#include <protocols/sic_dock/util.hh>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <utility/file/file_sys_util.hh>
	#include <utility/fixedsizearray1.hh>

	#include <numeric/geometry/hashing/SixDHasher.hh>
	#include <numeric/HomogeneousTransform.hh>

	#include <protocols/sic_dock/read_biounit.hh>

	#include <boost/unordered_set.hpp>
	#include <boost/cstdint.hpp>
	#include <boost/foreach.hpp>


namespace protocols {
namespace motif_hash {

// static data
static thread_local basic::Tracer TR( "protocols.motif_hash" );

// types

	typedef numeric::xyzVector<core::Real> Vec;
	typedef numeric::xyzMatrix<core::Real> Mat;
	typedef numeric::xyzTransform<core::Real> Xform;
	typedef utility::fixedsizearray1<float,20> float20;
	typedef utility::fixedsizearray1<float,9> float9;
	using std::make_pair;
	using core::chemical::AA;
	using numeric::HomogeneousTransform;
	using core::id::AtomID;
	using basic::options::option;
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
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::I;
	using ObjexxFCL::string_of;
	using std::cerr;
	using std::cout;
	using std::endl;
	using std::string;
	using utility::io::izstream;
	using utility::io::ozstream;
	using utility::file_basename;
	using utility::vector1;
	using std::endl;
	using core::import_pose::pose_from_pdb;
	using numeric::geometry::hashing::Real3;
	using numeric::geometry::hashing::Real6;
	using protocols::sic_dock::sqr;

// helpers
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

	std::string tag_from_pdb_fname(std::string const & fname0){
		string fname = utility::file_basename(fname0);
		if(fname.substr(fname.size()-3,3)==".gz" ) fname = fname.substr(0,fname.size()-3);
		if     (fname.substr(fname.size()-4,4)==".pdb") fname = fname.substr(0,fname.size()-4);
		else if(fname.substr(fname.size()-5,4)==".pdb") fname = fname.substr(0,fname.size()-5) + fname[fname.size()-1];
		if(fname.size()!=5){
			fname = fname.substr(0,4)+"S";
		}
		return fname;
	}

	std::ostream & operator<<(std::ostream & out, Real6 const & r6){
		out << r6[1] <<' ' << r6[2] <<' ' << r6[3] <<' ' << r6[4] <<' ' << r6[5] <<' ' << r6[6];
		return out;
	}

	Real
	rt6_rt6_bb_dis2_explicit_stupid(
		Real6 const & x1,
		Real6 const & x2
	){
		numeric::HomogeneousTransform<float> ht1,ht2;
		ht1.from_euler_angles_deg(numeric::xyzVector<Real>(x1[4],x1[5],x1[6]));
		ht2.from_euler_angles_deg(numeric::xyzVector<Real>(x2[4],x2[5],x2[6]));
		ht1.set_point(Vec(x1[1],x1[2],x1[3]));
		ht2.set_point(Vec(x2[1],x2[2],x2[3]));
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

	Real rt6_rt6_dis2(Real6 const & x1, Real6 const & x2, Real const lever){
		numeric::HomogeneousTransform<Real> ht1,ht2;
		ht1.from_euler_angles_deg(numeric::xyzVector<Real>(x1[4],x1[5],x1[6]));
		ht2.from_euler_angles_deg(numeric::xyzVector<Real>(x2[4],x2[5],x2[6]));
		HomogeneousTransform<Real> r = ht1 * ht2.inverse();
		Real cos_theta = (r.xx()+r.yy()+r.zz()-1.0)/2.0;
		Real sin2_theta = 1.0-cos_theta*cos_theta;
		Real dis2 = Vec(x1[1],x1[2],x1[3]).distance_squared(Vec(x2[1],x2[2],x2[3])) + lever*lever*sin2_theta;
		return dis2;
	}
	Real rt6_rt6_bb_dis2(Real6 const & x1, Real6 const & x2){
		return rt6_rt6_dis2(x1,x2,std::sqrt(6.0));
	}

	Real6 inverse_rt6(Real6 const & rt){
		numeric::HomogeneousTransform<Real> x;
		x.from_euler_angles_deg(numeric::xyzVector<Real>(rt[4],rt[5],rt[6]));
		x.set_point(Vec(rt[1],rt[2],rt[3]));
		x = x.inverse();
		numeric::xyzVector<Real> euler_angles =  x.euler_angles_deg();
		Real6 rt6;
		rt6[1] = x.px();
		rt6[2] = x.py();
		rt6[3] = x.pz();
		rt6[4] = fmod(euler_angles.x(),360.0);
		rt6[5] = fmod(euler_angles.y(),360.0);
		rt6[6] = fmod(euler_angles.z(),360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		return rt6;
	}

	Real6
	rt_to_real6(core::kinematics::RT const & rt){
		numeric::HomogeneousTransform<Real> ht( rt.get_rotation() , rt.get_translation() );
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
		numeric::HomogeneousTransform<Real> x;
		x.from_euler_angles_deg(numeric::xyzVector<Real>(rt6[4],rt6[5],rt6[6]));
		numeric::xyzMatrix<Real> R( numeric::xyzMatrix<Real>::cols_constructor(x.xaxis(),x.yaxis(),x.zaxis()) );
		numeric::xyzVector<Real> t(rt6[1],rt6[2],rt6[3]);
		return core::kinematics::RT(R,t);
	}

	Real6
	get_residue_pair_rt6(Pose const & pose1, Size ir, Pose const & pose2, Size jr){
		Xform x1(pose1.xyz(AtomID(2,ir)),pose1.xyz(AtomID(1,ir)),pose1.xyz(AtomID(2,ir)),pose1.xyz(AtomID(3,ir)));
		Xform x2(pose2.xyz(AtomID(2,jr)),pose2.xyz(AtomID(1,jr)),pose2.xyz(AtomID(2,jr)),pose2.xyz(AtomID(3,jr)));
		Xform x ( ~x1 * x2 );
		Vec e = x.euler_angles_deg();
		if(e.x()<0||e.x()>360.0||e.y()<0||e.y()>360.0){
			cout << e << endl;
			cout << ir << endl;
			cout << jr << endl;
		}
		Real6 rt;
		rt[1] = x.t.x();
		rt[2] = x.t.y();
		rt[3] = x.t.z();
		rt[4] = e.x();
		rt[5] = e.y();
		rt[6] = e.z();
		return rt;
	}


	Real6
	get_residue_pair_rt6(Pose const & pose, Size ir, Size jr){
		///////////// TODO ///////////
		// does it matter much how this is defined? best would be(?):
		//      cen = CAs midpoint
		//        X = CA1 -> CA2
		//        Y = ?
		core::id::StubID id1( AtomID(2,ir), AtomID(1,ir), AtomID(2,ir), AtomID(3,ir) );
		core::id::StubID id2( AtomID(2,jr), AtomID(1,jr), AtomID(2,jr), AtomID(3,jr) );
		core::kinematics::RT rt = pose.conformation().get_stub_transform(id1,id2);
		return rt_to_real6(rt);
	}

	void
	set_residue_pair_rt(Real6 const & rt6, Pose & pose, Size ir, Size jr){
		core::id::StubID id1( AtomID(2,ir), AtomID(1,ir), AtomID(2,ir), AtomID(3,ir) );
		core::id::StubID id2( AtomID(2,jr), AtomID(1,jr), AtomID(2,jr), AtomID(3,jr) );
		pose.conformation().set_stub_transform(id1,id2,real6_to_rt(rt6));
	}


	Real uint16_to_real(uint16_t const & val, Real const & lb, Real const & ub){
		return min(ub,max(lb,((Real)val)/((Real)MAX_UINT16) * (ub-lb) + lb));
	}
	uint16_t real_to_uint16(Real const & val, Real const & lb, Real const & ub){
		return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT16 );
	}

	Real uint8_to_real(uint8_t const & val, Real const & lb, Real const & ub){
		return min(ub,max(lb,((Real)val)/((Real)MAX_UINT8) * (ub-lb) + lb));
	}
	uint8_t real_to_uint8(Real const & val, Real const & lb, Real const & ub){
		return (uint16_t)( min(1.0,max(0.0,(val-lb)/(ub-lb))) * (Real)MAX_UINT8 );
	}

	vector1<Real> get_sasa(Pose const & pose, Real const & probesize){
		core::id::AtomID_Map<Real> atom_sasa;
		vector1<Real> rsd_sasa;
		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probesize, false );
		if(rsd_sasa.size() != pose.n_residue()) utility_exit_with_message("bad sasa!");
		return rsd_sasa;
	}
	vector1<Real> get_nbrs(Pose const & pose){
		vector1<Real> nbrs(pose.n_residue(),0);
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
	// 	vector1<Real> rsd_sasa;
	// 	core::id::AtomID_Map<bool> atom_subset;
	// 	core::pose::initialize_atomid_map( atom_subset, pose, false );
	// 	for(Size ia = 5; ia <= pose.residue(ir).nheavyatoms(); ++ia) atom_subset[AtomID(ia,ir)] = true;
	// 	return core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, 1.4, false, atom_subset );
	// }
	// Real get_sc_sasa(Pose const & pose, Size const & ir, Size const & jr){
	// 	core::id::AtomID_Map<Real> atom_sasa;
	// 	vector1<Real> rsd_sasa;
	// 	core::id::AtomID_Map<bool> atom_subset;
	// 	core::pose::initialize_atomid_map( atom_subset, pose, false );
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

	void HACK_dump_helix(core::pose::Pose const & pose, string fname, int beg, int end){
		cout << "dump helix to " << fname << endl;
		utility::io::ozstream out("helices/"+fname);
		Size ano = 1;
		for(int i = beg; i <= end; ++i){
			core::io::pdb::dump_pdb_residue(pose.residue(i),ano,out);
		}
		out.close();

	}
	int HACK_dump_helices(core::pose::Pose const & pose, string tag, int nres, int minlen=10){
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

	utility::fixedsizearray1<Real,6> get_bins(Real c, Real a){
		utility::fixedsizearray1<Real,6> bins(c);
		bins[4] = a;
		bins[5] = a;
		bins[6] = a;
		return bins;
	}

	Real
	align_motif_pose_super(
		Pose & pose,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr
	){
		using namespace core::id;

		// using rms super
		Pose ref;
		ref.append_residue_by_jump(paln1.residue(ir),1);
		ref.append_residue_by_jump(paln2.residue(jr),1);
		AtomID_Map<AtomID> atommap;
		core::pose::initialize_atomid_map(atommap,pose,BOGUS_ATOM_ID);
		atommap[AtomID(1,1)] = AtomID(1,1);
		atommap[AtomID(2,1)] = AtomID(2,1);
		atommap[AtomID(3,1)] = AtomID(3,1);
		atommap[AtomID(1,2)] = AtomID(1,2);
		atommap[AtomID(2,2)] = AtomID(2,2);
		atommap[AtomID(3,2)] = AtomID(3,2);
		return core::scoring::superimpose_pose(pose,ref,atommap);

		// //align separately to bb's
		// Xform stub1(paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(1,ir)),paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(3,ir)));
		// Xform stub2(paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(1,jr)),paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(3,jr)));
		// Xform mstub1(pose.xyz(AtomID(2,1)),pose.xyz(AtomID(1,1)),pose.xyz(AtomID(2,1)),pose.xyz(AtomID(3,1)));
		// Xform mstub2(pose.xyz(AtomID(2,2)),pose.xyz(AtomID(1,2)),pose.xyz(AtomID(2,2)),pose.xyz(AtomID(3,2)));
		// protocols::sic_dock::xform_pose(pose,stub1*~mstub1,1,1);
		// protocols::sic_dock::xform_pose(pose,stub2*~mstub2,2,2);

	}
	Real
	align_motif_pose_break(
		Pose & pose,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr
	){
		using namespace core::id;
		Pose start(pose);
		//align separately to bb's
		Xform stub1(paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(1,ir)),paln1.xyz(AtomID(2,ir)),paln1.xyz(AtomID(3,ir)));
		Xform stub2(paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(1,jr)),paln2.xyz(AtomID(2,jr)),paln2.xyz(AtomID(3,jr)));
		Xform mstub1(pose.xyz(AtomID(2,1)),pose.xyz(AtomID(1,1)),pose.xyz(AtomID(2,1)),pose.xyz(AtomID(3,1)));
		Xform mstub2(pose.xyz(AtomID(2,2)),pose.xyz(AtomID(1,2)),pose.xyz(AtomID(2,2)),pose.xyz(AtomID(3,2)));
		protocols::sic_dock::xform_pose(pose,stub1*~mstub1,1,1);
		protocols::sic_dock::xform_pose(pose,stub2*~mstub2,2,2);

		return core::scoring::bb_rmsd(pose,start);

	}


// ResPairMotif
	void ResPairMotif::reset(){
		pdb_=0;
		rt6_=0;
		resi1_=0;
		resi2_=0;
		chain1_=0;
		chain2_=0;
		ss1_=0;
		ss2_=0;
		aa1_=0;
		aa2_=0;
		chi1_[1]=0;
		chi1_[2]=0;
		chi1_[3]=0;
		chi1_[4]=0;
		chi2_[1]=0;
		chi2_[2]=0;
		chi2_[3]=0;
		chi2_[4]=0;
		fa_atr_=0;
		hb_sc_=0;
		hb_bb_sc_=0;
		hb_bb_=0;
		sasa_=0;
		nbrs_=0;
		count_=1;
		misc1_=0;
		misc2_=0;
		misc3_=0;
		misc4_=0;
		misc5_=0;
		bfac1_=0;
		bfac2_=0;
		lg1x_=0;
		lg1y_=0;
		lg1z_=0;
		lg2x_=0;
		lg2y_=0;
		lg2z_=0;
	}
	ResPairMotif::ResPairMotif(){
		reset();
	}
	ResPairMotif::ResPairMotif(
		std::string const & tag,
		core::pose::Pose const & _pose,
		Real6 const & _xform,
		Size const & _resi1,   Size const & _resi2,
		Real const & _sasa ,   Real const & _nbrs, Real const & _fa_atr,
		Real const & _hb_sc,   Real const & _hb_bb_sc, Real const & _hb_bb,
		Real const & _bfac1,   Real const & _bfac2, bool is_sspaired
	){
		reset();
		if( tag.size()!=5 && tag.size()!=8 ) utility_exit_with_message("bad tag for ResPairMotif not nchar 5 or 8");
		pdb_[1]=tag[0]; pdb_[2]=tag[1]; pdb_[3]=tag[2]; pdb_[4]=tag[3]; pdb_[5]=tag[4];
		if(tag.size()==8) {
			pdb_[6]=tag[5]; pdb_[7]=tag[6]; pdb_[8]=tag[7];
		} else {
			pdb_[6]='_'; pdb_[7]='_'; pdb_[8]='_';
		}

		this->rt(_xform);

		resi1_  = _resi1;
		resi2_  = _resi2;
		if( is_sspaired ){
			ss1_ = 'P';
			ss2_ = 'P';
		} else {
			ss1_ = _pose.secstruct(_resi1);
			ss2_ = _pose.secstruct(_resi2);
		}
		aa1_ = _pose.residue(_resi1).name1();
		aa2_ = _pose.residue(_resi2).name1();
		chain1_ = basic::get_pymol_chain(_pose.chain(_resi1));
		chain2_ = basic::get_pymol_chain(_pose.chain(_resi2));

		this->bfac1   (_bfac1);
		this->bfac2   (_bfac2);
		this->sasa    (_sasa);
		this->nbrs    (_nbrs);
		this->fa_atr  (_fa_atr);
		this->hb_sc   (_hb_sc);
		this->hb_bb_sc(_hb_bb_sc);
		this->hb_bb   (_hb_bb);

		for(Size i = 1; i <= _pose.residue(_resi1).nchi(); ++i) chi1_[i] = real_to_uint8(_pose.residue(_resi1).chi(i),-180.0,180.0);
		for(Size i = 1; i <= _pose.residue(_resi2).nchi(); ++i) chi2_[i] = real_to_uint8(_pose.residue(_resi2).chi(i),-180.0,180.0);
	}
	char ResPairMotif::aa1() const { return aa1_; }
	char ResPairMotif::aa2() const { return aa2_; }
	char ResPairMotif::ss1() const { if(ss1_=='H'||ss1_=='I'||ss1_=='G') return 'H'; else if(ss1_=='P'||ss1_=='E'||ss1_=='U'||ss1_=='B') return 'E'; else if(ss1_=='S'||ss1_=='T'||ss1_==' '||ss1_=='L') return 'L'; else utility_exit_with_message("unknown SS "+string_of(ss1_)); }
	char ResPairMotif::ss2() const { if(ss2_=='H'||ss2_=='I'||ss2_=='G') return 'H'; else if(ss2_=='P'||ss2_=='E'||ss2_=='U'||ss2_=='B') return 'E'; else if(ss2_=='S'||ss2_=='T'||ss2_==' '||ss2_=='L') return 'L'; else utility_exit_with_message("unknown SS "+string_of(ss2_)); }
	char ResPairMotif::dssp1() const { return ss1_; }
	char ResPairMotif::dssp2() const { return ss2_; }
	void ResPairMotif::reverse_in_place(){
		ht( ht().inverse() );
		std::swap(resi1_  ,resi2_);
		std::swap(chain1_ ,chain2_);
		std::swap(ss1_    ,ss2_);
		std::swap(aa1_    ,aa2_);
		std::swap(chi1_   ,chi2_);
		std::swap(bfac1_  ,bfac2_);
		// lg1x_,lg1y_,lg1z_; // dummy for lig / metal / water ???????????
		// lg2x_,lg2y_,lg2z_; // euler ang for lig? ???????????
	}

	Real6 ResPairMotif::rt() const {
		Real6 rt6;
		for(Size i = 1; i <= 3; ++i) rt6[i] = uint16_to_real(rt6_[i],-16.0,  16.0 );
		for(Size i = 4; i <= 5; ++i) rt6[i] = uint16_to_real(rt6_[i],  0.0, 360.0 );
		for(Size i = 6; i <= 6; ++i) rt6[i] = uint16_to_real(rt6_[i],  0.0, 180.0 );
		return rt6;
	}
	void  ResPairMotif::rt(Real6 const & rt_in) {
		Real6 rt6 = rt_in;
		rt6[4] = fmod(rt_in[4],360.0);
		rt6[5] = fmod(rt_in[5],360.0);
		rt6[6] = fmod(rt_in[6],360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		for(Size i = 1; i <= 3; ++i) rt6_[i] = real_to_uint16(rt_in[i],-16.0,  16.0 );
		for(Size i = 4; i <= 5; ++i) rt6_[i] = real_to_uint16(rt_in[i],  0.0, 360.0 );
		for(Size i = 6; i <= 6; ++i) rt6_[i] = real_to_uint16(rt_in[i],  0.0, 180.0 );
	}
	HomogeneousTransform<Real> ResPairMotif::ht() const{
		Real6 rt6 = rt();
		numeric::HomogeneousTransform<Real> x;
		x.from_euler_angles_deg(numeric::xyzVector<Real>(rt6[4],rt6[5],rt6[6]));
		x.set_point(Vec(rt6[1],rt6[2],rt6[3]));
		return x;
	}
	void ResPairMotif::ht(HomogeneousTransform<Real> const & x){
		numeric::xyzVector < Real > euler_angles =  x.euler_angles_deg();
		Real6 rt6;
		rt6[1] = x.px();
		rt6[2] = x.py();
		rt6[3] = x.pz();
		rt6[4] = fmod(euler_angles.x(),360.0);
		rt6[5] = fmod(euler_angles.y(),360.0);
		rt6[6] = fmod(euler_angles.z(),360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		if(rt6[6] > 180.0) utility_exit_with_message("bad rt6[6]");
		rt(rt6);
	}

	Real ResPairMotif::dist2() const {
		Real sqdis = 0.0;
		for(Size i = 1; i <= 3; ++i){
			Real d = uint16_to_real(rt6_[i],
			-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](),
			 basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]());
			sqdis += d*d;
		}
		return sqdis;
	}

	// char ss_pair_to_code(char const & ss1, char const & ss2){
	// 	char code(0);
	// 	code = ('H'==ss1 && 'H'==ss2 ) ? (char)1 : code;
	// 	code = ('H'==ss1 && 'E'==ss2 ) ? (char)2 : code;
	// 	code = ('H'==ss1 && 'L'==ss2 ) ? (char)3 : code;
	// 	code = ('E'==ss1 && 'H'==ss2 ) ? (char)4 : code;
	// 	code = ('E'==ss1 && 'E'==ss2 ) ? (char)5 : code;
	// 	code = ('E'==ss1 && 'L'==ss2 ) ? (char)6 : code;
	// 	code = ('L'==ss1 && 'H'==ss2 ) ? (char)7 : code;
	// 	code = ('L'==ss1 && 'E'==ss2 ) ? (char)8 : code;
	// 	code = ('L'==ss1 && 'L'==ss2 ) ? (char)9 : code;
	// 	return code;
	// }

	Real ResPairMotif::bfac1   (                      ) const { return      uint8_to_real(  bfac1_    ,  0.0,  2.55  ); }
	void ResPairMotif::bfac1   (Real const & _bfac1   )       { bfac1_    = real_to_uint8( _bfac1     ,  0.0,  2.55  ); }
	Real ResPairMotif::bfac2   (                      ) const { return      uint8_to_real(  bfac2_    ,  0.0,  2.55  ); }
	void ResPairMotif::bfac2   (Real const & _bfac2   )       { bfac2_    = real_to_uint8( _bfac2     ,  0.0,  2.55  ); }

	Real ResPairMotif::sasa    (                      ) const { return      uint8_to_real(  sasa_     ,  0.0,  1.0   ); }
	void ResPairMotif::sasa    (Real const & _sasa    )       { sasa_     = real_to_uint8( _sasa      ,  0.0,  1.0   ); }
	Real ResPairMotif::nbrs    (                      ) const { return      uint8_to_real(  nbrs_     ,  0.0, 63.75  ); }
	void ResPairMotif::nbrs    (Real const & _nbrs    )       { nbrs_     = real_to_uint8( _nbrs      ,  0.0, 63.75  ); }

	Real ResPairMotif::fa_atr  (                      ) const { return      uint8_to_real(  fa_atr_   , -7.96875,   0.0  ); }
	void ResPairMotif::fa_atr  (Real const & _fa_atr  )       { fa_atr_   = real_to_uint8( _fa_atr    , -7.96875,   0.0  ); }
	Real ResPairMotif::hb_sc   (                      ) const { return      uint8_to_real(  hb_sc_    , -7.96875,   0.0  ); }
	void ResPairMotif::hb_sc   (Real const & _hb_sc   )       { hb_sc_    = real_to_uint8( _hb_sc     , -7.96875,   0.0  ); }
	Real ResPairMotif::hb_bb_sc(                      ) const { return      uint8_to_real(  hb_bb_sc_ , -7.96875,   0.0  ); }
	void ResPairMotif::hb_bb_sc(Real const & _hb_bb_sc)       { hb_bb_sc_ = real_to_uint8( _hb_bb_sc  , -7.96875,   0.0  ); }
	Real ResPairMotif::hb_bb   (                      ) const { return      uint8_to_real(  hb_bb_    , -7.96875,   0.0  ); }
	void ResPairMotif::hb_bb   (Real const & _hb_bb   )       { hb_bb_    = real_to_uint8( _hb_bb     , -7.96875,   0.0  ); }

	Real ResPairMotif::chi11() const { return uint8_to_real(chi1_[1],-180.0,180.0); }
	Real ResPairMotif::chi12() const { return uint8_to_real(chi1_[2],-180.0,180.0); }
	Real ResPairMotif::chi13() const { return uint8_to_real(chi1_[3],-180.0,180.0); }
	Real ResPairMotif::chi14() const { return uint8_to_real(chi1_[4],-180.0,180.0); }
	Real ResPairMotif::chi21() const { return uint8_to_real(chi2_[1],-180.0,180.0); }
	Real ResPairMotif::chi22() const { return uint8_to_real(chi2_[2],-180.0,180.0); }
	Real ResPairMotif::chi23() const { return uint8_to_real(chi2_[3],-180.0,180.0); }
	Real ResPairMotif::chi24() const { return uint8_to_real(chi2_[4],-180.0,180.0); }

	bool ResPairMotif::filter() const {
		using std::string;
		namespace f = basic::options::OptionKeys::mh::filter;
		if( option[f::seqsep   ].user() && option[f::seqsep]() >= abs(resi1_-resi2_) && chain1_==chain2_ ) return false;
		if( option[f::no_hb_bb ].user() && option[f::no_hb_bb ]() && hb_bb() != 0.0)  return false;
		if( option[f::no_hb_bb ].user() && option[f::no_hb_bb ]() && (ss1_=='P'||ss2_=='P') ) return false;
		if( option[f::ss1      ].user() && option[f::ss1      ]()[0] != ss1()    ) return false;
		if( option[f::ss2      ].user() && option[f::ss2      ]()[0] != ss2()    ) return false;
		if( option[f::dssp1    ].user() && option[f::dssp1    ]()[0] != dssp1()  ) return false;
		if( option[f::dssp2    ].user() && option[f::dssp2    ]()[0] != dssp2()  ) return false;
		if( option[f::aa1      ].user() && option[f::aa1      ]()[0] != aa1_     ) return false;
		if( option[f::aa2      ].user() && option[f::aa2      ]()[0] != aa2_     ) return false;
		if( option[f::sasa     ].user() && option[f::sasa     ]() >  sasa()      ) return false;
		if( option[f::faatr    ].user() && option[f::faatr    ]() <  fa_atr()    ) return false;
		if( option[f::hb_sc    ].user() && option[f::hb_sc    ]() <  hb_sc()     ) return false;
		if( option[f::hb_bb_sc ].user() && option[f::hb_bb_sc ]() <  hb_bb_sc()  ) return false;
		if( option[f::hb_bb    ].user() && option[f::hb_bb    ]() <  hb_bb()     ) return false;
		if( option[f::coorderr ].user() && option[f::coorderr ]() <  bfac1()     ) return false;
		if( option[f::coorderr ].user() && option[f::coorderr ]() <  bfac2()     ) return false;
		Real x = uint16_to_real(rt6_[1],-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](), basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]() );
		Real y = uint16_to_real(rt6_[2],-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](), basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]() );
		Real z = uint16_to_real(rt6_[3],-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](), basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]() );
		Real dis2 = x*x+y*y+z*z;
		if( option[f::mindist2 ].user() && option[f::mindist2 ]() > dis2 ) return false;
		if( option[f::maxdist2 ].user() && option[f::maxdist2 ]() < dis2 ) return false;

		if( option[f::faatr_or_hbbb].user()       &&
			option[f::faatr_or_hbbb] < fa_atr()   &&
			option[f::faatr_or_hbbb] < hb_sc()    &&
			option[f::faatr_or_hbbb] < hb_bb_sc() &&
			option[f::faatr_or_hbbb] < hb_bb()	  ) return false;

		if( option[f::faatr_or_hb].user()       &&
			option[f::faatr_or_hb] < fa_atr()   &&
			option[f::faatr_or_hb] < hb_sc()	  ) return false;
		// if( option[mh::lg1   ].user() && option[mh::lg1   ]() != lg1x_     ) return false;
		// if( option[mh::lg2   ].user() && option[mh::lg2   ]() != lg2x_     ) return false;

		if( option[f::noloops ]() && (ss1()=='L' || ss2() =='L')  ) return false;
		if( option[f::oneloop ]() && (ss1()=='L' && ss2() =='L')  ) return false;
		if( option[f::nodisulf]() && (aa1_=='C' && aa2_=='C')  ) return false;

		if( option[f::    restype    ].user() && (option[f::    restype    ]().find(aa1_)==string::npos || option[f::    restype    ]().find(aa2_)==string::npos) ) return false;
		if( option[f::    restype_one].user() && (option[f::    restype_one]().find(aa1_)==string::npos && option[f::    restype_one]().find(aa2_)==string::npos) ) return false;
		if( option[f::not_restype    ].user() && (option[f::not_restype    ]().find(aa1_)!=string::npos && option[f::not_restype    ]().find(aa2_)!=string::npos) ) return false;
		if( option[f::not_restype_one].user() && (option[f::not_restype_one]().find(aa1_)!=string::npos || option[f::not_restype_one]().find(aa2_)!=string::npos) ) return false;
		return true;
	}
	void
	ResPairMotif::dump_pdb(
		std::string const & fname,
		std::string const & tag
	) const {
		utility::io::ozstream out(fname);
		dump_pdb(out,Xform(),tag);
		out.close();
	}
	Real
	ResPairMotif::dump_aligned_motif_super(
		std::ostream & out,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr,
		int const & tag
	) const {
		using namespace core::id;
		Pose pose;
		fill_pose_with_motif(pose,ir,paln1.n_residue()+jr);
		Real rms = align_motif_pose_super(pose,paln1,ir,paln2,jr);
		out << "MODEL " << tag << endl;
		core::io::pdb::dump_pdb(pose,out);
		out << "ENDMDL" << endl;
		return rms;
	}
	Real
	ResPairMotif::dump_aligned_motif_break(
		std::ostream & out,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr,
		int const & tag
	) const {
		using namespace core::id;
		Pose pose;
		fill_pose_with_motif(pose,ir,paln1.n_residue()+jr);
		Real rms = align_motif_pose_break(pose,paln1,ir,paln2,jr);
		out << "MODEL " << tag << endl;
		core::io::pdb::dump_pdb(pose,out);
		out << "ENDMDL" << endl;
		return rms;
	}
	Real
	ResPairMotif::dump_aligned_motif(
		std::ostream & out,
		Pose const & paln1, Size const & ir,
		Pose const & paln2, Size const & jr,
		int const & tag
	) const {
		return dump_aligned_motif_super(out,paln1,ir,paln2,jr,tag);
	}
	Real
	ResPairMotif::dump_aligned_motif(
		std::ostream & out,
		Pose const & paln,
		Size const & ir,
		Size const & jr,
		int const & tag
	) const {
		return dump_aligned_motif(out,paln,ir,paln,jr,tag);
	}

	void
	ResPairMotif::fill_pose_with_motif(
		Pose & pose,
		int ir, int jr
	) const {
		string seq = std::string("")+aa1_+aa2_;
		// cerr << "make pose " << seq << endl;
		core::pose::make_pose_from_sequence(pose,seq,"fa_standard");
		core::pose::remove_lower_terminus_type_from_pose_residue(pose,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(pose,2);
		core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( pose, true ) );
		pdb_info->number(1,ir);
		pdb_info->number(2,jr);
		pose.pdb_info( pdb_info );
		// pose.dump_pdb("testempty.pdb");
		for(Size i=1; i<=pose.residue(1).nchi(); ++i){ pose.set_chi(i,1,uint8_to_real(chi1_[i],-180.0,180.0)); }
		for(Size i=1; i<=pose.residue(2).nchi(); ++i){ pose.set_chi(i,2,uint8_to_real(chi2_[i],-180.0,180.0)); }
		// cerr << "make fold tree" << endl;
		core::kinematics::FoldTree ft(pose.fold_tree());
		ObjexxFCL::FArray2D<int> jumps(2,1); jumps(1,1)=1;jumps(2,1)=2;
		ObjexxFCL::FArray1D<int> cuts(1,1); cuts(1)=1;
		ft.tree_from_jumps_and_cuts(2,1,jumps,cuts);
		pose.fold_tree(ft);
		// cerr << "set rt" << endl;
		set_residue_pair_rt(this->rt(),pose,1,2);
	}

	void
	ResPairMotif::dump_pdb(
		std::ostream & out,
		Xform const & x = Xform(),
		std::string const & tag
	) const {
		Pose pose;
		fill_pose_with_motif(pose);
		Xform xl;
		xl.from_four_points( pose.residue(1).xyz(2),pose.residue(1).xyz(1),pose.residue(1).xyz(2),pose.residue(1).xyz(3) );
		protocols::sic_dock::xform_pose(pose,x * ~xl );
		// cerr << "dump_pdb" << endl;
		out << "MODEL " << tag << endl;
		core::io::pdb::dump_pdb(pose,out);
		out << "ENDMDL" << endl;
	}
	void ResPairMotif::print_header(std::ostream & out){
	//	out << "ResPairMotif pdbid lig   -0.132   -2.504   -3.917   -3.738 -171.551  107.124   36 E E A 0.360 -175.76 -180.00    0.71 -180.00  157 F E A 0.300  -71.29   76.94 -180.00 -180.00  99.205  -2.896   0.000   0.000   0.000" << endl;
		out << "ResPairMotif pdbid lig       dx       dy       dz       ex       ey       ez res1 A S C crder    chi1    chi2    chi3    chi4 res2 A S C crder    chi1    chi2    chi3    chi4  sasa   faatr   hb_sc  hbscbb   hb_bb" << endl;
	}
		 // 8  utility::fixedsizearray1<unsigned char,8> pdb_; // simplify this // last3 are lig resn
		// 12  utility::fixedsizearray1<uint16_t,6> rt6_;
		 // 4  uint16_t resi1_,resi2_;
		 // 2  char chain1_,chain2_;
		 // 2  char ss1_,ss2_;
		 // 2  uint8_t aa1_,aa2_;
		 // 2  uint8_t  bfac1_,bfac2_;
		 // 8  utility::fixedsizearray1<uint8_t,4> chi1_,chi2_;
		// 12  uint16_t sasa_,fa_atr_,hb_sc_,hb_bb_sc_,hb_bb_,MISC;
		 // 6  uint16_t lg1x_,lg1y_,lg1z_; // dummy for lig / metal / water
		 // 6  uint16_t lg2x_,lg2y_,lg2z_; // euler ang for lig?
	std::ostream & operator<<(std::ostream & out, ResPairMotif const & x){
		Real6 const rt = x.rt();
		out << x.pdb_[1]<<x.pdb_[2]<<x.pdb_[3]<<x.pdb_[4]<<x.pdb_[5]<<' '<<x.pdb_[6]<<x.pdb_[7]<<x.pdb_[8];
		out << I(4,x.resi1_) <<" ";
		out << I(4,x.resi2_) <<" ";
		out << x.aa1_ <<"";
		out << x.aa2_ <<" ";
		out << x.ss1_ <<"";
		out << x.ss2_ <<" ";
		out << x.chain1_ <<"";
		out << x.chain2_ <<" ";
		out << F(5,3,x.bfac1())    <<" ";
		out << F(5,3,x.bfac2())    <<" ";
		out << F(5,3,x.sasa())  <<" ";
		out << F(4,1,x.nbrs())  <<" ";
		out << "e ";
		out << I(3,x.count())   <<" ";
		out << F(7,3,x.fa_atr())   <<" ";
		out << F(7,3,x.hb_sc())    <<" ";
		out << F(7,3,x.hb_bb_sc()) <<" ";
		out << F(7,3,x.hb_bb())    <<" ";
		out << "rt " << F(8,3,rt[1])<<' '<<F(8,3,rt[2])<<' '<<F(8,3,rt[3])<<' '<<F(8,3,rt[4])<<' '<<F(8,3,rt[5])<<' '<<F(8,3,rt[6]) <<" ";
		out << "chi " << I(4,(int)uint8_to_real(x.chi1_[1],-180.0,180.0))<<" "<<I(4,(int)uint8_to_real(x.chi1_[2],-180.0,180.0))<<' '<<I(4,(int)uint8_to_real(x.chi1_[3],-180.0,180.0))<<' '<<I(4,(int)uint8_to_real(x.chi1_[4],-180.0,180.0))<<' ';
		out << I(4,(int)uint8_to_real(x.chi2_[1],-180.0,180.0))<<' '<<I(4,(int)uint8_to_real(x.chi2_[2],-180.0,180.0))<<' '<<I(4,(int)uint8_to_real(x.chi2_[3],-180.0,180.0))<<' '<<I(4,(int)uint8_to_real(x.chi2_[4],-180.0,180.0))     ;
		return out;
	}
	bool write_motifs_binary(std::ostream & out, ResPairMotifs const & motifs){
		boost::uint64_t n = motifs.size();
		out.write((char*)&n,sizeof(boost::uint64_t));
		for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
			ResPairMotif const & sm( *i );
			out.write((char*)&sm,sizeof(ResPairMotif));
			if(!out.good()) return false;
		}
		return true;
	}
	bool write_motifs_binary(std::string const & fname, ResPairMotifs const & motifs){
		utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
		if(!out.good()){
			TR.Error << "write_motifs_binary(fname): problem opening motif output file " << fname << endl;
			return false;
		}
		if(!write_motifs_binary(out,motifs)){
			TR.Error << "write_motifs_binary(fname): problem while writing motif file " << fname << endl;
			return false;
		}
		if(!out.good()){
			TR.Error << "write_motifs_binary(fname): problem after writing motif file " << fname << endl;
			out.close();
			return false;
		}
		out.close();
		return true;
	}
	bool read_motifs_binary(std::istream & in, ResPairMotifs & motifs){
		while(!in.eof()){
			boost::uint64_t n0=motifs.size(),n=0;
			in.read((char*)&n,sizeof(boost::uint64_t));
			if(option[basic::options::OptionKeys::mh::generate_reverse_motifs]()){
				motifs.resize(motifs.size()+2*n);
				for(Size i = 1; i <= n; ++i){
					in.read((char*)&motifs[n0+2*i-1],sizeof(ResPairMotif));
					motifs[n0+2*i] = motifs[n0+2*i-1];
					motifs[n0+2*i].reverse_in_place();
				}
			} else {
				motifs.resize(motifs.size()+n);
				for(Size i = 1; i <= n; ++i){
					in.read((char*)&motifs[n0+i],sizeof(ResPairMotif));
				}
			}
		}
		return true;
	}
	bool read_motifs_binary(std::string const & fname, ResPairMotifs & motifs){
		utility::io::izstream in(fname,std::ios::in|std::ios::binary);
		if(!in.good()){
			TR.Error << "read_motifs_binary(fname): problem opening motif input file " << fname << endl;
			return false;
		}
		if(!read_motifs_binary(in,motifs)){
			TR.Error << "read_motifs_binary(fname): problem while reading motif file " << fname << endl;
			return false;
		}
		// if(!in.good()){
		// 	TR.Error << "problem after reading motif file " << fname << endl;
		// 	in.close();
		// 	return false;
		// }
		// in.close();
		return true;
	}
	bool read_motifs_binary(vector1<string> const & fnames, ResPairMotifs & motifs){
		TR << "Nfiles: " << fnames.size() << endl;
		for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
			// if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
			TR << "reading binary ResPairMotifs " << *i << endl;
			if(!read_motifs_binary(*i,motifs)){
				TR.Error << "read_motifs_binary(fnames): error reading file "+*i << endl;
				if(!option[basic::options::OptionKeys::mh::ignore_io_errors]()) utility_exit_with_message("");
			}
		}
		TR << endl;
		return true;
	}

	void load_motifs(vector1<string> const & fnames, ResPairMotifs & motifs){
		ResPairMotifs raw_motifs;
		if(!read_motifs_binary(fnames,raw_motifs)) utility_exit_with_message("error reading motif files");
		TR << "read " << ((Real)raw_motifs.size())/1000000.0 << " million binary motifs from " << fnames.size() << " files" << endl;
		motifs.clear();
		if(option[basic::options::OptionKeys::mh::filter::filter_io]()){
			filter_motifs(raw_motifs,motifs);
			TR << "filter: keep " << ((Real)motifs.size())/1000000.0 << " million binary motifs from " << fnames.size() << " files" << endl;
		}
		if(option[basic::options::OptionKeys::mh::merge_similar_motifs].user()){
			Real cartsize = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[1];
			Real cartresl = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[2];
			Real anglresl = option[basic::options::OptionKeys::mh::merge_similar_motifs]()[3];
			motifs.filter_structurally_similar_motifs(cartsize,cartresl,anglresl);
			TR << "merge:  keep " << ((Real)motifs.size())/1000000.0 << " million binary motifs from " << fnames.size() << " files" << endl;
		}
	}

	void load_motifs(string const & fname, ResPairMotifs & motifs){
		load_motifs(utility::vector1<string>(1,fname),motifs);
	}



	struct MyHash {
		boost::uint64_t operator()(ResPairMotif const & rpm) const {
			boost::uint64_t const *x = (boost::uint64_t const *)(&rpm);
			return x[0] ^ x[1] ^ x[2] ^ x[3] ^ x[4] ^ x[5];
		}
	};
	struct MyPred {
		bool operator()(ResPairMotif const & a, ResPairMotif const & b) const {
			return 0 == memcmp(&a,&b,48);
		}
	};
	typedef boost::unordered_set<ResPairMotif,MyHash,MyPred> MotifSet;
	void ResPairMotifs::filter_structurally_identical_motifs(){
		MotifSet motifset;
		BOOST_FOREACH(ResPairMotif const & rpm, *this) motifset.insert(rpm);
		clear();
		reserve(motifset.size());
		BOOST_FOREACH(ResPairMotif const & rpm, motifset) push_back(rpm);
	}
	void ResPairMotifs::filter_structurally_similar_motifs(Real cart_size, Real cart_resl, Real ang_resl){
		for(int iofst = 0; iofst < 30; ++iofst){
			Real6 shift(0);
			if(iofst>0){ shift[1]=uniform(); shift[2]=uniform(); shift[3]=uniform(); shift[4]=uniform(); shift[5]=uniform(); shift[6]=0; }
			MotifHash mh(cart_size + (iofst? cart_resl/2.0 : 0.0)  ,cart_resl,ang_resl);
			ResPairMotifs all_keepers;
			BOOST_FOREACH(ResPairMotif const & rpm, *this) mh.add_motif_shift(rpm,shift);
			for(MotifHash::KeySet::const_iterator i = mh.begin(); i != mh.end(); ++i){
				ResPairMotifs this_bin,keepers;
				mh.find_motifs(*i,this_bin);
				BOOST_FOREACH(ResPairMotif const & m, this_bin){
					bool seenit = false;
					BOOST_FOREACH(ResPairMotif & n, keepers){
						bool seenthis = ( m.aa1()==n.aa1() && m.aa2()==n.aa2() && m.ss1()==n.ss1() && m.ss2()==n.ss2() );
						if( min(fabs(m.chi11()-n.chi11()),360.0-fabs(m.chi11()-n.chi11())) > ang_resl*1.5 ) seenthis = false;
						if( min(fabs(m.chi21()-n.chi21()),360.0-fabs(m.chi21()-n.chi21())) > ang_resl*1.5 ) seenthis = false;
						if( min(fabs(m.chi12()-n.chi12()),360.0-fabs(m.chi12()-n.chi12())) > ang_resl*2.0 ) seenthis = false;
						if( min(fabs(m.chi22()-n.chi22()),360.0-fabs(m.chi22()-n.chi22())) > ang_resl*2.0 ) seenthis = false;
						if( min(fabs(m.chi13()-n.chi13()),360.0-fabs(m.chi13()-n.chi13())) > ang_resl*2.5 ) seenthis = false;
						if( min(fabs(m.chi23()-n.chi23()),360.0-fabs(m.chi23()-n.chi23())) > ang_resl*2.5 ) seenthis = false;
						if( min(fabs(m.chi14()-n.chi14()),360.0-fabs(m.chi14()-n.chi14())) > ang_resl*3.0 ) seenthis = false;
						if( min(fabs(m.chi24()-n.chi24()),360.0-fabs(m.chi24()-n.chi24())) > ang_resl*3.0 ) seenthis = false;
						if(seenthis){
							seenit = true;
							n.addcount();
							break;
						}
					}
					if(!seenit) keepers.push_back(m);
				}
				all_keepers.insert(all_keepers.end(),keepers.begin(),keepers.end());
			}
			*this = all_keepers;
			cout << iofst << " " << size() << endl;
		}
	}

// Xfrags
	Xfres::Xfres(char _aa, char _ss, Real _phi, Real _psi, Real _omg, Real _chi1, Real _chi2, Real _chi3, Real _chi4){
		aa  (_aa  );
		ss  (_ss  );
		phi (_phi );
		psi (_psi );
		omg (_omg );
		chi1(_chi1);
		chi2(_chi2);
		chi3(_chi3);
		chi4(_chi4);
	}
	Real Xfres::phi () const { return uint16_to_real(phi_  ,-180.0,180.0); }
	Real Xfres::psi () const { return uint16_to_real(psi_  ,-180.0,180.0); }
	Real Xfres::omg () const { return uint16_to_real(omg_  ,-180.0,180.0); }
	Real Xfres::chi1() const { return uint8_to_real(chi_[1],-180.0,180.0); }
	Real Xfres::chi2() const { return uint8_to_real(chi_[2],-180.0,180.0); }
	Real Xfres::chi3() const { return uint8_to_real(chi_[3],-180.0,180.0); }
	Real Xfres::chi4() const { return uint8_to_real(chi_[4],-180.0,180.0); }
	void Xfres::phi (Real _phi) { phi_    = real_to_uint16(_phi,-180.0,180.0); }
	void Xfres::psi (Real _psi) { psi_    = real_to_uint16(_psi,-180.0,180.0); }
	void Xfres::omg (Real _omg) { omg_    = real_to_uint16(_omg,-180.0,180.0); }
	void Xfres::chi1(Real _chi1){ chi_[1] = real_to_uint8(_chi1,-180.0,180.0); }
	void Xfres::chi2(Real _chi2){ chi_[2] = real_to_uint8(_chi2,-180.0,180.0); }
	void Xfres::chi3(Real _chi3){ chi_[3] = real_to_uint8(_chi3,-180.0,180.0); }
	void Xfres::chi4(Real _chi4){ chi_[4] = real_to_uint8(_chi4,-180.0,180.0); }
	void Xfres::place_sidechain_in_pose(Pose & pose, int ir) const{
		string const & name( core::chemical::name_from_aa(core::chemical::aa_from_oneletter_code(aa())) );
		core::conformation::ResidueOP newres = core::conformation::ResidueFactory::create_residue(pose.residue(ir).residue_type_set().name_map(name));
		pose.replace_residue(ir,*newres,true);
		if(pose.residue(ir).nchi() >= 1) pose.set_chi(1,ir,chi1());
		if(pose.residue(ir).nchi() >= 2) pose.set_chi(2,ir,chi2());
		if(pose.residue(ir).nchi() >= 3) pose.set_chi(3,ir,chi3());
		if(pose.residue(ir).nchi() >= 4) pose.set_chi(4,ir,chi4());
	}
	std::ostream & operator << (std::ostream & out, Xfres const & x){
		out<<"Xfres "<<x.aa()<<" "<<x.ss()<<" "<<F(8,3,x.phi())<<" "<<F(8,3,x.psi())<<" "<<F(8,3,x.omg())<<" "<<F(8,3,x.chi1())<<" "<<F(8,3,x.chi2())<<" "<<F(8,3,x.chi3())<<" "<<F(8,3,x.chi4());
		return out;
	}
	std::istream & operator >> (std::istream & in , Xfres & x){
		string tag; in >> tag; if("Xfres"!=tag) utility_exit_with_message("bad tag for Xfres");
		char _aa,_ss; Real _phi,_psi,_omg,_chi1,_chi2,_chi3,_chi4;
		in>>_aa>>_ss>>_phi>>_psi>>_omg>>_chi1>>_chi2>>_chi3>>_chi4;
		x.aa  (_aa  );
		x.ss  (_ss  );
		x.phi (_phi );
		x.psi (_psi );
		x.omg (_omg );
		x.chi1(_chi1);
		x.chi2(_chi2);
		x.chi3(_chi3);
		x.chi4(_chi4);
		return in;
	}
	bool read_xfres_binary( string  const & fname , vector1<Xfres> & xfres){
		utility::io::izstream in(fname,std::ios::in|std::ios::binary);
		if(!in.good()){
			TR.Error << "read_xfres_binary(fname): problem opening xfres input file " << fname << endl;
			return false;
		}
		if(!read_xfres_binary(in,xfres)){
			TR.Error << "read_xfres_binary(fname): problem while reading xfres file " << fname << endl;
			return false;
		}
		return true;
	}
	bool read_xfres_binary( std::istream & in, vector1<Xfres> & xfres){
		boost::uint64_t n0=xfres.size(),n=0;
		in.read((char*)&n,sizeof(boost::uint64_t));
		xfres.resize(xfres.size()+n);
		for(Size i = 1; i <= n; ++i){
			in.read((char*)&xfres[n0+i],sizeof(Xfres));
		}
		return true;
	}
	bool write_xfres_binary( std::ostream & out , vector1<Xfres> const & xfres){
		boost::uint64_t n = xfres.size();
		out.write((char*)&n,sizeof(boost::uint64_t));
		cout << "write_xfres_binary " << n << endl;
		for(vector1<Xfres>::const_iterator i = xfres.begin(); i != xfres.end(); ++i){
			Xfres const & sm( *i );
			out.write((char*)&sm,sizeof(Xfres));
			if(!out.good()) return false;
		}
		return true;
	}
	bool write_xfres_binary( string  const & fname , vector1<Xfres> const & xfres){
		utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
		if(!out.good()){
			TR.Error << "write_xfres_binary(fname): problem opening xfres output file " << fname << endl;
			return false;
		}
		if(!write_xfres_binary(out,xfres)){
			TR.Error << "write_xfres_binary(fname): problem while writing xfres file " << fname << endl;
			return false;
		}
		if(!out.good()){
			TR.Error << "write_xfres_binary(fname): problem after writing xfres file " << fname << endl;
			out.close();
			return false;
		}
		out.close();
		return true;
	}
	bool read_xfres_binary( vector1<string> const & fnames, vector1<Xfres> & xfres){
		Size count = 0;
		TR << "Nfiles: " << fnames.size() << endl;
		for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
			if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
			if(!read_xfres_binary(*i,xfres)){
				TR.Error << "read_xfres_binary(fnames): error reading file "+*i << endl;
				return false;
			}
		}
		TR << endl;
		return true;
	}
	Xfrag::Xfrag(Real6 _rt, Size _position, uint8_t _size, char _sscomp, Real _bfac, Real _burial){
		rt6(_rt);
		position(_position);
		sscomp(_sscomp);
		size(_size);
		bfac(_bfac);
		burial(_burial);
	}
	Real6 Xfrag::rt6() const {
		Real6 rt6;
		for(Size i = 1; i <= 3; ++i) rt6[i] = uint16_to_real(rt6_[i],
			-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](),
			 basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]());
		for(Size i = 4; i <= 5; ++i) rt6[i] = uint16_to_real(rt6_[i], 0.0, 360.0 );
		for(Size i = 6; i <= 6; ++i) rt6[i] = uint16_to_real(rt6_[i], 0.0, 180.0 );
		return rt6;
	}
	void  Xfrag::rt6(Real6 const & rt_in) {
		Real6 rt6 = rt_in;
		rt6[4] = fmod(rt_in[4],360.0);
		rt6[5] = fmod(rt_in[5],360.0);
		rt6[6] = fmod(rt_in[6],360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		for(Size i = 1; i <= 3; ++i) rt6_[i] = real_to_uint16(rt_in[i],
			-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size](),
			 basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]());
		for(Size i = 4; i <= 5; ++i) rt6_[i] = real_to_uint16(rt_in[i], 0.0, 360.0 );
		for(Size i = 6; i <= 6; ++i) rt6_[i] = real_to_uint16(rt_in[i], 0.0, 180.0 );
	}
	Size    Xfrag::position() const { return position_; }
	Size    Xfrag::size    () const { return size_; }
	char    Xfrag::sscomp  () const { return sscomp_; }
	Real    Xfrag::bfac    () const { return uint8_to_real(bfac_,  0.0, 2.55); }
	Real    Xfrag::burial  () const { return uint8_to_real(burial_,0.0, 1.0 ); }
	Real    Xfrag::ex1     () const { return uint8_to_real(ex1_,   0.0, 1.0 ); }
	Real    Xfrag::ex2     () const { return uint8_to_real(ex2_,   0.0, 1.0 ); }
	Real    Xfrag::ex3     () const { return uint8_to_real(ex3_,   0.0, 1.0 ); }
	Real    Xfrag::ex4     () const { return uint8_to_real(ex4_,   0.0, 1.0 ); }
	void    Xfrag::position(Size    _position) { position_ = _position; }
	void    Xfrag::size    (Size    _size    ) { size_     = _size; }
	void    Xfrag::sscomp  (char    _sscomp  ) { sscomp_   = _sscomp; }
	void    Xfrag::bfac    (Real    _bfac    ) { bfac_     = real_to_uint8(_bfac,  0.0, 2.55); }
	void    Xfrag::burial  (Real    _burial  ) { burial_   = real_to_uint8(_burial,0.0, 1.0 ); }
	void    Xfrag::ex1     (Real    _ex1     ) { ex1_      = real_to_uint8(_ex1,   0.0, 1.0 ); }
	void    Xfrag::ex2     (Real    _ex2     ) { ex2_      = real_to_uint8(_ex2,   0.0, 1.0 ); }
	void    Xfrag::ex3     (Real    _ex3     ) { ex3_      = real_to_uint8(_ex3,   0.0, 1.0 ); }
	void    Xfrag::ex4     (Real    _ex4     ) { ex4_      = real_to_uint8(_ex4,   0.0, 1.0 ); }
	void    Xfrag::insert (Pose & pose, vector1<Xfres> const & xfres, int lowres, int highres) const {
		int cutres=0;
		for(int ic=1; ic <= pose.fold_tree().num_cutpoint(); ++ic){
			int c = pose.fold_tree().cutpoint(ic);
			if(c < lowres || c >= highres) continue;
			if(cutres) utility_exit_with_message("more than one cutpoint in insertion");
			cutres = c;
		}
		if(!cutres) utility_exit_with_message("must be a cutpoint in inertion");
		core::kinematics::FoldTree ft( pose.fold_tree() );
		for(int ij=1; ij <= (int)pose.num_jump(); ++ij){
			int ur = pose.fold_tree().upstream_jump_residue(ij);
			int dr = pose.fold_tree().downstream_jump_residue(ij);
			if( lowres < ur && ur <= cutres  ) { ft.slide_jump(ij, lowres, dr);      cout << "ft.slide_jump("<<ij<<", "<<lowres   <<", "<<dr<<");" << endl; }
			if( lowres < dr && dr <= cutres  ) { ft.slide_jump(ij, ur, lowres);      cout << "ft.slide_jump("<<ij<<", "<<ur       <<", "<<lowres<<");" << endl; }
			if( cutres < ur && ur <= highres ) { ft.slide_jump(ij, highres+1, dr);   cout << "ft.slide_jump("<<ij<<", "<<highres+1<<", "<<dr<<");" << endl; }
			if( cutres < dr && dr <= highres ) { ft.slide_jump(ij, ur, highres+1);   cout << "ft.slide_jump("<<ij<<", "<<ur       <<", "<<highres+1<<");" << endl; }
			cout << ft;
			for(int i=1;i<=(int)ft.nres();++i) cout << i << " " << ft.is_jump_point(i) << endl;
		}
		ft.slide_cutpoint(cutres, highres);
		ft.reorder(1);
		cout << ft;
		for(int i=1;i<=(int)ft.nres();++i) cout << i << " " << ft.is_jump_point(i) << endl;
		cout << "OLD " << pose.fold_tree();
		pose.fold_tree(ft);
		cout << "NEW " << pose.fold_tree();

		for(int ir=lowres+1; ir <= highres; ++ir){
			cout << pose.fold_tree();
			cout << "delete " << lowres+1 << endl;
			pose.delete_polymer_residue(lowres+1);
		}
		core::conformation::ResidueOP dummyres = core::conformation::ResidueFactory::create_residue(pose.residue(1).residue_type_set().name_map("GLY"));
		for(int ir=2; ir <= (int)size(); ++ir) pose.append_polymer_residue_after_seqpos(*dummyres,lowres-2+ir,true);
		vector1<Xfres>::const_iterator ifrag = xfres.begin() + position();
		for(int ir=1; ir <= (int)size(); ++ir,++ifrag){
			int resno = lowres+ir-1;
			pose.set_phi  (resno,ifrag->phi());
			pose.set_psi  (resno,ifrag->psi());
			pose.set_omega(resno,ifrag->omg());
			ifrag->place_sidechain_in_pose(pose,resno);
		}
	}

	std::ostream & operator << (std::ostream & out, Xfrag const & x){
		out<<"Xfrag "<<I(9,x.position())<<" "<<I(2,x.size())<<" "<<x.sscomp()<<" "<<F(7,3,x.bfac())<<" "<<F(7,3,x.burial())<<" "<<F(7,3,x.ex1())<<" "<<F(7,3,x.ex2())<<" "<<F(7,3,x.ex3())<<" "<<F(7,3,x.ex4());
		out<<F(7,3,x.rt6()[1])<<" "<<F(7,3,x.rt6()[2])<<" "<<F(7,3,x.rt6()[3])<<" "<<F(7,3,x.rt6()[4])<<" "<<F(7,3,x.rt6()[5])<<" "<<F(7,3,x.rt6()[6]);
		return out;
	}
	std::istream & operator >> (std::istream & in , Xfrag & x){
		string tag; in >> tag; if("Xfrag"!=tag) utility_exit_with_message("bad tag for Xfrag");
		char _position,_size,_sscomp; Real _bfac,_burial,_ex1,_ex2,_ex3,_ex4;
		in>>_position>>_size>>_sscomp>>_bfac>>_burial>>_ex1>>_ex2>>_ex3>>_ex4;
		x.position  (_position  );
		x.size  (_size  );
		x.sscomp (_sscomp );
		x.bfac (_bfac );
		x.burial (_burial );
		x.ex1(_ex1);
		x.ex2(_ex2);
		x.ex3(_ex3);
		x.ex4(_ex4);
		return in;
	}
	bool read_xfrag_binary( string  const & fname , vector1<Xfrag> & xfrag){
		utility::io::izstream in(fname,std::ios::in|std::ios::binary);
		if(!in.good()){
			TR.Error << "read_xfrag_binary(fname): problem opening xfrag input file " << fname << endl;
			return false;
		}
		if(!read_xfrag_binary(in,xfrag)){
			TR.Error << "read_xfrag_binary(fname): problem while reading xfrag file " << fname << endl;
			return false;
		}
		return true;
	}
	bool read_xfrag_binary( std::istream & in, vector1<Xfrag> & xfrag){
		boost::uint64_t n0=xfrag.size(),n=0;
		in.read((char*)&n,sizeof(boost::uint64_t));
		xfrag.resize(xfrag.size()+n);
		for(Size i = 1; i <= n; ++i){
			in.read((char*)&xfrag[n0+i],sizeof(Xfrag));
		}
		return true;
	}
	bool write_xfrag_binary( std::ostream & out , vector1<Xfrag> const & xfrag){
		boost::uint64_t n = xfrag.size();
		out.write((char*)&n,sizeof(boost::uint64_t));
		cout << "write_xfrag_binary " << n << endl;
		for(vector1<Xfrag>::const_iterator i = xfrag.begin(); i != xfrag.end(); ++i){
			Xfrag const & sm( *i );
			out.write((char*)&sm,sizeof(Xfrag));
			if(!out.good()) return false;
		}
		return true;
	}
	bool write_xfrag_binary( string  const & fname , vector1<Xfrag> const & xfrag){
		utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
		if(!out.good()){
			TR.Error << "write_xfrag_binary(fname): problem opening xfrag output file " << fname << endl;
			return false;
		}
		if(!write_xfrag_binary(out,xfrag)){
			TR.Error << "write_xfrag_binary(fname): problem while writing xfrag file " << fname << endl;
			return false;
		}
		if(!out.good()){
			TR.Error << "write_xfrag_binary(fname): problem after writing xfrag file " << fname << endl;
			out.close();
			return false;
		}
		out.close();
		return true;
	}
	bool read_xfrag_binary( vector1<string> const & fnames, vector1<Xfrag> & xfrag){
		Size count = 0;
		TR << "Nfiles: " << fnames.size() << endl;
		for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
			if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
			if(!read_xfrag_binary(*i,xfrag)){
				TR.Error << "read_xfrag_binary(fnames): error reading file "+*i << endl;
				return false;
			}
		}
		TR << endl;
		return true;
	}




	bool read_xfrag_binary (        string  const & fname , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
		utility::io::izstream in(fname,std::ios::in|std::ios::binary);
		if(!in.good()){
			TR.Error << "read_xfrag_binary(fname): problem opening xfrag input file " << fname << endl;
			return false;
		}
		if(!read_xfrag_binary(in,xfrag,xfres)){
			TR.Error << "read_xfrag_binary(fname): problem while reading xfrag file " << fname << endl;
			return false;
		}
		return true;

	}
	bool read_xfrag_binary (  std::istream        & in    , vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
		if(!read_xfres_binary(in,xfres)) return false;
		if(!read_xfrag_binary(in,xfrag)) return false;
		return true;
	}
	bool write_xfrag_binary(  std::ostream        & out   , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres){
		if(!write_xfres_binary(out,xfres)) return false;
		if(!write_xfrag_binary(out,xfrag)) return false;
		return true;
	}
	bool write_xfrag_binary(        string  const & fname , vector1<Xfrag> const & xfrag, vector1<Xfres> const & xfres){
		utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
		if(!out.good()){
			TR.Error << "write_xfrag_binary(fname): problem opening xfrag output file " << fname << endl;
			return false;
		}
		if(!write_xfrag_binary(out,xfrag,xfres)){
			TR.Error << "write_xfrag_binary(fname): problem while writing xfrag file " << fname << endl;
			return false;
		}
		if(!out.good()){
			TR.Error << "write_xfrag_binary(fname): problem after writing xfrag file " << fname << endl;
			out.close();
			return false;
		}
		out.close();
		return true;
	}
	bool read_xfrag_binary (vector1<string> const & fnames, vector1<Xfrag>       & xfrag, vector1<Xfres>       & xfres){
		Size count = 0;
		TR << "Nfiles: " << fnames.size() << endl;
		for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
			if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
			if(!read_xfrag_binary(*i,xfrag,xfres)){
				TR.Error << "read_xfrag_binary(fnames): error reading file "+*i << endl;
				return false;
			}
		}
		TR << endl;
		return true;

	}


	XfragSet::XfragSet(core::Real /*cartsize*/, core::Real cartresl, core::Real angleresl)
	 : cart_resl_(cartresl), angle_resl_(angleresl),
			hasher_(
			numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(
				Vec(-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]()),
				Vec( basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]())
			),
			utility::fixedsizearray1<Size,3>(0),
			get_bins(cart_resl_,angle_resl_)
		)
 {}

// MotifHash
	MotifHash::MotifHash():
		cart_size_ (basic::options::option[basic::options::OptionKeys::mh::hash_cart_size ]()),
		cart_resl_ (basic::options::option[basic::options::OptionKeys::mh::hash_cart_resl ]()),
		angle_resl_(basic::options::option[basic::options::OptionKeys::mh::hash_angle_resl]()),
		hasher_(
			numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec(cart_size_)),
			utility::fixedsizearray1<Size,3>(0),
			get_bins(cart_resl_,angle_resl_)
		)
	{}
	MotifHash::MotifHash(Real cart_size, Real cart_resl, Real ang_resl):
		cart_size_ (cart_size ),
		cart_resl_ (cart_resl ),
		angle_resl_( ang_resl),
		hasher_(
			numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec( cart_size_)),
			utility::fixedsizearray1<Size,3>(0),
			get_bins(cart_resl_,angle_resl_)
		)
	{}
	MotifHash::MotifHash( ResPairMotifs const & motifs ):
		cart_size_ (basic::options::option[basic::options::OptionKeys::mh::hash_cart_size ]()),
		cart_resl_ (basic::options::option[basic::options::OptionKeys::mh::hash_cart_resl ]()),
		angle_resl_(basic::options::option[basic::options::OptionKeys::mh::hash_angle_resl]()),
		hasher_(
			numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(Vec(-cart_size_),Vec( cart_size_)),
			utility::fixedsizearray1<Size,3>(0),
			get_bins(cart_resl_,angle_resl_)
		)
	{
		for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
			add_motif(*i);
		}
	}

	MotifHash::Key MotifHash::bin_index(Real6 const & rt) const {
		Real6 rt6 = rt;
		rt6[4] = fmod(rt6[4],360.0);
		rt6[5] = fmod(rt6[5],360.0);
		rt6[6] = fmod(rt6[6],360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		return hasher_.bin_index(rt6);
	}
	MotifHash::Key MotifHash::bin_index(Motif const & m) const {
		return bin_index(m.rt());
	}
	void MotifHash::add_motif(Motif const & d, Key const & key){
		motif_hash_.insert( std::make_pair( key , d ) );
		key_set_.insert(key);
		// Motif tmp(d);         // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should I do this???
		// tmp.reverse_in_place();
		// motif_hash_.insert( std::make_pair( bin_index(tmp.rt()) , tmp ) );
	}
	void MotifHash::add_motif(Motif const & d){
		Real6 const rt( d.rt() );
		if( fabs(rt[1]) > cart_size_ ) return;
		if( fabs(rt[2]) > cart_size_ ) return;
		if( fabs(rt[3]) > cart_size_ ) return;
		Key key = bin_index(rt);
		add_motif(d,key);
	}
	void MotifHash::add_motif_shift(Motif const & d, Real6 const & shift){
		Real6 rt( d.rt() );
		if( fabs(rt[1]) > cart_size_ ) return;
		if( fabs(rt[2]) > cart_size_ ) return;
		if( fabs(rt[3]) > cart_size_ ) return;
		rt[1] = rt[1] + cart_resl_ * shift[1];
		rt[2] = rt[2] + cart_resl_ * shift[2];
		rt[3] = rt[3] + cart_resl_ * shift[3];
		rt[4] = fmod(rt[4] + angle_resl_*shift[4],360.0);
		rt[5] = fmod(rt[5] + angle_resl_*shift[5],360.0);
		rt[6] = fmod(rt[6] + angle_resl_*shift[6],180.0);
		Key key = bin_index(rt);
		add_motif(d,key);
	}
	void MotifHash::set_score(Key const & k, float const & score){
		score_hash_[k] = score;
	}
	void MotifHash::find_motifs(Key const & k, ResPairMotifs & results ) const {
		std::pair<MotifMap::const_iterator,MotifMap::const_iterator> range = motif_hash_.equal_range(k);
		for( MotifMap::const_iterator it = range.first; it != range.second; ++it) {
			results.push_back( it->second );
		}
	}
	void MotifHash::find_motifs(Real6 const & rt6, ResPairMotifs & results ) const {
		find_motifs(bin_index(rt6),results);
	}
	int MotifHash::count_motifs(Key const & k ) const {
		std::pair<MotifMap::const_iterator,MotifMap::const_iterator> range = motif_hash_.equal_range(k);
		int c = 0;
		for( MotifMap::const_iterator it = range.first; it != range.second; ++it) ++c;
		return c;
	}
	int MotifHash::count_motifs(Real6 const & rt6) const {
		return count_motifs(bin_index(rt6));
	}


	void MotifHash::find_motifs_with_radius(Real6 const & rt, Real radius, vector1<Motif> & results) const {
		Real cart_frac = radius / cart_resl_;
		Real ang2dis = cart_resl_/angle_resl_;
		for(int lookuprad = 0; lookuprad <= (int)(cart_frac+1.9999); ++lookuprad){
			if(lookuprad>5) continue;
			Xform x1,x2;
			x1.from_euler_angles_deg(Vec(rt[4],rt[5],rt[6]));
			std::vector<MotifHash::Key> surrounding_cells = hasher_.radial_bin_index(lookuprad,rt);
			// cout << surrounding_cells.size() << endl;
			// utility_exit_with_message("foo");
			BOOST_FOREACH(MotifHash::Key key,surrounding_cells){
				MotifHash::MotifMap::const_iterator i = motif_hash_.equal_range(key).first;
				MotifHash::MotifMap::const_iterator e = motif_hash_.equal_range(key).second;
				for(; i != e; ++i){
					ResPairMotif const & sm(i->second);
					Real6 ot = sm.rt();
					Real dis2 = sqr(rt[1]-ot[1])+sqr(rt[2]-ot[2])+sqr(rt[3]-ot[3]);
					x2.from_euler_angles_deg(Vec(ot[4],ot[5],ot[6]));
					Mat r = x1.R * x2.R.transposed();
					Real cos_theta = (r.xx()+r.yy()+r.zz()-1.0)/2.0;
					// Real sin2_theta = 1.0-cos_theta*cos_theta;
					Real theta = degrees(acos(cos_theta));
					dis2 += ang2dis*ang2dis*theta*theta;
					if( dis2 <= radius*radius ) results.push_back(sm);
				}
			}
		}
	}

	float MotifHash::motif_score(Real6 const & rt) const {
		return motif_score( bin_index(rt) );
	}
	void MotifHash::generate_scoring_hash(){
		utility_exit_with_message("generate_scoring_hash");
		// score all cells
	}
	// template<class Archive>
	// void
	// MotifHash::serialize(Archive & ar, const unsigned int /*version*/){
	// 	// ar & hasher_;
	// 	// ar & motif_hash_;
	// 	ar & score_hash_;
	// }
	std::ostream & operator << (std::ostream & out, MotifHash const & x){
		numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> > const & bb(x.hasher_.bounding_box() );
		Real3 const & eo( x.hasher_.euler_offsets() );
		Real6 const & bw( x.hasher_.bin_widths() );
		out << bb.lower() << ' ' << bb.upper() << endl;
		out << eo[1] << ' ' << eo[2] << ' ' << eo[3] << endl;
		out << bw[1] << ' ' << bw[2] << ' ' << bw[3] << ' ' << bw[4] << ' ' << bw[5] << ' ' << bw[6] << endl;
		out << x.motif_hash_.size() << endl;
		out << x.motif_hash_.begin()->second << endl;
		// for(vector1<MotifHash::Motif>::const_iterator i = x.motifs_.begin(); i != x.motifs_.end(); ++i){
		// 	out << *i << endl;
		// }
		// MotifMap motif_hash_;
		// ScoreMap score_hash_;
		return out;
	}
	// std::istream & operator >> (std::istream & in , MotifHash & x){
	// 	if(x.motifs_.size()) utility_exit_with_message("try to input to non-empty motif hash");
	// 	numeric::xyzVector<numeric::Real> lb,ub;
	// 	Real eo1,eo2,eo3,bw1,bw2,bw3,bw4,bw5,bw6;
	// 	in >> lb >> ub;
	// 	in >> eo1 >> eo2 >> eo3 >> bw1 >> bw2 >> bw3 >> bw4 >> bw5 >> bw6;

	// 	numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> > const & bb(x.hasher_.bounding_box() );
	// 	Real3 const & eo( x.hasher_.euler_offsets() );
	// 	Real6 const & bw( x.hasher_.bin_widths() );
	// 	if( bb.lower() != lb || bb.upper() != ub ||
	// 	    eo[1]!=eo1 || eo[2]!=eo2 || eo[3]!=eo3 ||
	// 	    bw[1]!=bw1 || bw[2]!=bw2 || bw[3]!=bw3 || bw[4]!=bw4 || bw[5]!=bw5 || bw[6]!=bw6 ) {
	// 		utility_exit_with_message("MotifHash params differ, old version of stored?");
	// 	}


	// 	Size nmotif; in >> nmotif;
	// 	x.motifs_.resize(nmotif);
	// 	for(Size i = 1; i <= nmotif; ++i){
	// 		in >> x.motifs_[i];
	// 	}
	// 	return in;
	// }

	////////////////////////////// motif /////////////////////////////////////////////////

// protocols

	void print_motifs(std::ostream & out){
		using namespace basic::options::OptionKeys;
		ResPairMotifs motifs;
		load_motifs( option[mh::print_motifs](), motifs );
		for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
			out << *i << endl;
		}
	}

	void merge_motifs(){
		using namespace basic::options::OptionKeys;
		if( ! basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user() ){
			utility_exit_with_message("must sepcify name for merged file -motif_out_file");
		}
		vector1<string> const & fnames ( option[mh::merge_motifs]() );
		string          const & outfile( option[mh::motif_out_file ]()+".smtf.bin.gz" );
		ResPairMotifs motifs; load_motifs( fnames, motifs );
		TR << motifs.size() << endl;
		motifs.filter_structurally_identical_motifs();
		TR << motifs.size() << endl;

		// !!!!!!!!!!!!!


		TR << "writing them to binary file " << outfile << endl;
		write_motifs_binary(outfile,motifs);
	}

	void dump_motif_pdbs(){
		using namespace basic::options::OptionKeys;

		TR << "dump_motif_pdbs: read motifs" << endl;
		vector1<string> const & fnames( option[mh::dump_motif_pdbs]() );
		ResPairMotifs motifs; load_motifs( fnames, motifs );

		TR << "dump_motif_pdbs: hash motifs" << endl;
		MotifHash mh;
		for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
			// TR << "ADD_KEY " << mh.bin_index(*i) << " " <<  i->rt() << endl;
			mh.add_motif(*i);
		}

		TR << "dump_motif_pdbs: find largest bucket / check counts" << endl;
		boost::uint64_t maxnum=0;
		vector1<boost::uint64_t> keys;
		if( !option[mh::dump_motif_pdbs_min_counts].user() ) keys.resize(1);
		for(MotifHash::MotifMap::const_iterator i = mh.motif_hash_.cbegin(); i != mh.motif_hash_.cend(); ++i){
			boost::uint64_t n = mh.motif_hash_.count(i->first);
			if( n > maxnum ){
				maxnum = n;
				if( option[mh::dump_motif_pdbs_min_counts].user() ){
					if((boost::uint64_t)option[mh::dump_motif_pdbs_min_counts]() <= n ){
						keys.push_back(i->first);
					}
				} else {
					keys[1] = i->first;
				}
			}
			// ResPairMotif const & sm(i->second);
			// string fn = string(option[mh::motif_out_file]())+
			// 	"_motif_"+
			// 	string_of(sm.aa1())    +
			// 	string_of(sm.aa2())    + "_"+
			// 	string_of(sm.pdb_[1]) + "_"+
			// 	string_of(sm.resi1_ ) + "_"+
			// 	string_of(sm.resi2_ ) + ".pdb";
			// TR << "DUMPING: " << fn << " " << sm << endl;
			// i->second.dump_pdb(fn);
		}

		TR << "dump_motif_pdbs: max motifs in bucket, dumping that one bucket: " << maxnum << " " << (Real)maxnum/(Real)motifs.size() << endl;
		string tag = string(option[mh::motif_out_file]());
		if(tag=="") tag = "NOTAG";
		ResPairMotif::print_header(TR);
		for(Size k = 1; k <= keys.size(); ++k){
			TR << "dump motif bin " << k << " " << endl;
			boost::uint64_t key = keys[k];
			MotifHash::MotifMap::const_iterator l = mh.motif_hash_.equal_range(key).first;
			MotifHash::MotifMap::const_iterator e = mh.motif_hash_.equal_range(key).second;
			for(; l != e; ++l){
				ResPairMotif const & sm(l->second);
				string fn = tag+
					+"_"+string_of(key)+
					"_largest_bucket_"    +
					string_of(sm.ss1())    +
					string_of(sm.ss2())    + "_"+
					string_of(sm.aa1())    +
					string_of(sm.aa2())    + "_"+
					string_of(sm.pdb_[1]) +string_of(sm.pdb_[2]) +string_of(sm.pdb_[3]) +string_of(sm.pdb_[4]) +string_of(sm.pdb_[5]) + "_"+
					string_of(sm.resi1_ ) + "_"+
					string_of(sm.resi2_ ) + ".pdb";
				TR << sm << " " << fn << endl;
				l->second.dump_pdb(fn);
			}
		}
	}

	void harvest_scores(){
		using namespace basic::options::OptionKeys;

		TR << "harvest_scores: read motifs" << endl;
		vector1<string> const & fnames( option[mh::harvest_scores]() );
		ResPairMotifs motifs; load_motifs( fnames, motifs );


		TR << "computing counts with requested bin sizes" << endl;
		XformScore xh;
		for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
			xh.add_xform_count(i->rt());
		}

		TR << xh << endl;

		xh.write_binary(option[mh::motif_out_file]()+".xh.bin.gz");

		// xh.print_scores(cout);
	}

	void print_scores(){
		using namespace basic::options::OptionKeys;
		XformScore xh;
		xh.read_binary(option[mh::print_scores]());
		TR << xh << endl;
		xh.print_scores(cout);
	}

	void dump_matching_motifs(){
		using namespace basic::options::OptionKeys;

		vector1<string> const & fnames ( option[mh::input_motifs]() );
		ResPairMotifs motifs; load_motifs( fnames, motifs );

		TR << "create motif hash" << endl;
		MotifHash mh(motifs);
		mh.hasher().tree_init(5);

		vector1<string> const & pdbfiles ( option[mh::dump_matching_motifs]() );

		for(int ifile = 1; ifile <= (int)pdbfiles.size(); ++ifile ){
			// cout << "READ BEGIN" << endl;
			string fname = pdbfiles[ifile];
			Pose pose;
			core::import_pose::centroid_pose_from_pdb(pose,fname);
			core::scoring::dssp::Dssp dssp(pose);
			dssp.insert_edge_ss_into_pose(pose);

			// cout << "SCORE BEGIN" << endl;
			for(Size ir = 1; ir <= pose.total_residue(); ++ir){
				for(Size jr = ir+1; jr <= pose.total_residue(); ++jr){
					// cout << "check res pair " << ir<<' '<<jr << endl;
					std::ostringstream oss;
					int nhits = 0;
					if(option[mh::score_across_chains_only]() && pose.chain(ir)==pose.chain(jr)) continue;
					float dsq = pose.xyz(AtomID(2,ir)).distance_squared(pose.xyz(AtomID(2,jr)));
					if(dsq > mh.cart_size()*mh.cart_size()) continue;
					Real6 rt = get_residue_pair_rt6(pose,ir,jr);
					std::vector<MotifHash::Key> surrounding_cells = mh.hasher().radial_bin_index(2,rt);
					BOOST_FOREACH(MotifHash::Key key,surrounding_cells){
						MotifHash::MotifMap::const_iterator i = mh.motif_hash_.equal_range(key).first;
						MotifHash::MotifMap::const_iterator e = mh.motif_hash_.equal_range(key).second;
						for(; i != e; ++i){
							ResPairMotif const & sm(i->second);
							float err = rt6_rt6_bb_dis2(rt,sm.rt());
							if( err <= option[mh::dump_matching_motifs_cutoff]()*option[mh::dump_matching_motifs_cutoff]() ){
								sm.dump_aligned_motif(oss,pose,ir,jr,ir+jr);
								++nhits;
							}
						}
					}
					if(nhits){
						cout << "hits " << ir << " " << jr << " " << nhits << endl;
						utility::io::ozstream o(utility::file_basename(fname)+"_motifs_"+string_of(ir)+"_"+string_of(jr)+".pdb");
						o << oss.str();
						o.close();
					}
				}
			}

		}
	}
	void filter_motifs(
		ResPairMotifs const & motifs_in,
		ResPairMotifs       & motifs_out
	){
		for(ResPairMotifs::const_iterator i = motifs_in.begin(); i != motifs_in.end(); ++i){
			if(i->filter()) motifs_out.push_back(*i);
		}
	}

	Real6 & reframe_rt_centered(Real6 & rt){
		return rt;
	}


// XformScore
	XformScore::XformScore():
		hasher_(
			numeric::geometry::BoundingBox<numeric::xyzVector<numeric::Real> >(
				Vec(-basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]()),
				Vec( basic::options::option[basic::options::OptionKeys::mh::hash_cart_size]())
			),
			utility::fixedsizearray1<Size,3>(0),
			get_bins(
				basic::options::option[basic::options::OptionKeys::mh::hash_cart_resl ](),
				basic::options::option[basic::options::OptionKeys::mh::hash_angle_resl]()
			)
		)
	{
		using namespace basic::options::OptionKeys;
		cart_size_  = option[mh::hash_cart_size  ]();
		cart_resl_  = option[mh::hash_cart_resl  ]();
		angle_resl_ = option[mh::hash_angle_resl ]();
		num_samples_ = 0;
	}
	void XformScore::clear() {
		total_count_.clear();
		num_samples_ = 0;
	}
	Size XformScore::num_bins() const { return total_count_.size(); }
	Size XformScore::num_samples() const { return num_samples_; }
	void XformScore::print_quantiles(std::ostream & out, int num) const {
		utility::vector1<Count> vals;
		for(CountMap::const_iterator i=total_count_.begin(); i != total_count_.end(); ++i){
			vals.push_back(i->second);
		}
		std::sort(vals.begin(),vals.end());
		int N = vals.size()/num;
		for(int i = 1; i <= (int)vals.size(); i += N){
			out <<" " <<vals[i];

		}
		out << " " << vals.back();
	}
	XformScore::Key XformScore::bin_index(Real6 const & rt) const {
		Real6 rt6 = rt;
		rt6[4] = fmod(rt6[4],360.0);
		rt6[5] = fmod(rt6[5],360.0);
		rt6[6] = fmod(rt6[6],360.0);
		rt6[4] = rt6[4]<0.0 ? rt6[4]+360.0 : rt6[4];
		rt6[5] = rt6[5]<0.0 ? rt6[5]+360.0 : rt6[5];
		rt6[6] = rt6[6]<0.0 ? rt6[6]+360.0 : rt6[6];
		return hasher_.bin_index(rt6);
	}
	void XformScore::add_xform_count(Real6 const & xform){
		Key const key = bin_index(xform);
		total_count_.emplace( std::make_pair(key,(Count)0) );
		++(total_count_.find(key)->second);
		++num_samples_;
	}
	void XformScore::set_score(Key const & key, Count const & count){ total_count_.insert( std::make_pair( key, count                     ) ); }
	void XformScore::add_score(Key const & key, Count const & count){ total_count_.insert( std::make_pair( key, count + score_of_bin(key) ) ); }

	XformScore::Count XformScore::score_of_bin(Real6 const & xform) const {
		Key key = bin_index(xform);
		CountMap::const_iterator i = total_count_.find(key);
		if( i != total_count_.end() ){
			return i->second;
		} else {
			return 0;
		}
	}
	XformScore::Count XformScore::score_of_bin(Key const & key) const {
		CountMap::const_iterator i = total_count_.find(key);
		if( i != total_count_.end() ){
			return i->second;
		} else {
			return 0;
		}
	}
	Real XformScore::score_of_bin_normalized(Real6 const & xform) const {
		Key key = bin_index(xform);
		Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
		Real count = score_of_bin(key);
		Real freq = count / sin(radians(center[6]));
		return freq;
	}
	bool XformScore::write_binary(std::ostream & out) const {
		numeric::Real  cart_size =  cart_size_;
		numeric::Real  cart_resl =  cart_resl_;
		numeric::Real angle_resl = angle_resl_;
		boost::uint64_t nmotifs = num_samples();
		boost::uint64_t ncounts = num_bins();
		out.write((char*)&   ncounts,sizeof(boost::uint64_t));
		out.write((char*)&   nmotifs,sizeof(boost::uint64_t));
		out.write((char*)& cart_size,sizeof(numeric::Real));
		out.write((char*)& cart_resl,sizeof(numeric::Real));
		out.write((char*)&angle_resl,sizeof(numeric::Real));
		for(CountMap::const_iterator i = total_count_.begin(); i != total_count_.end(); ++i){
			Key key = i->first;
			Count count = i->second;
			if(count < option[basic::options::OptionKeys::mh::harvest_scores_min_count]()) continue;
			out.write((char*)&key  ,sizeof(Key));
			out.write((char*)&count,sizeof(Count));
		}
		return true;
	}
	bool XformScore::write_binary(string const & fname) const {
		utility::io::ozstream out(fname,std::ios::out|std::ios::binary);
		if(!out.good()){
			TR.Error << "XformScore::write_binary(fname): problem opening XformScore output file " << fname << endl;
			return false;
		}
		if(!write_binary(out)){
			TR.Error << "XformScore::write_binary(fname): problem while writing XformScore file " << fname << endl;
			out.close();
			return false;
		}
		out.close();
		return true;
	}
	bool XformScore::read_binary(std::istream & in,bool clearme){
		if(clearme){
			clear();
		}
		numeric::Real  cart_size=0;
		numeric::Real  cart_resl=0;
		numeric::Real angle_resl=0;
		boost::uint64_t ncounts,nmotifs;
		in.read((char*)&   ncounts,sizeof(boost::uint64_t));
		in.read((char*)&   nmotifs,sizeof(boost::uint64_t));
		in.read((char*)& cart_size,sizeof(numeric::Real));
		in.read((char*)& cart_resl,sizeof(numeric::Real));
		in.read((char*)&angle_resl,sizeof(numeric::Real));
		if(cart_size != cart_size_ || cart_resl!= cart_resl_ || angle_resl != angle_resl_){
			cout << "from file " << cart_size  << " " << cart_resl  << " " << angle_resl  << endl;
			cout << "from this " << cart_size_ << " " << cart_resl_ << " " << angle_resl_ << endl;
			utility_exit_with_message("hash param mismatch!");
		}
		for(Size i = 1; i <= ncounts; ++i){
			Key key;
			Count count;
			in.read((char*)&key  ,sizeof(Key));
			in.read((char*)&count,sizeof(Count));
			add_score(key,count);
		}
		return true;
	}
	bool XformScore::read_binary(string const & fname, bool clearme){
		TR.Error << "XformScore::read_binary(fname): " << fname << endl;
		utility::io::izstream in(fname,std::ios::in|std::ios::binary);
		if(!in.good()){
			TR.Error << "XformScore::read_binary(fname): problem opening XformScore input file " << fname << endl;
			return false;
		}
		if(!read_binary(in,clearme)){
			TR.Error << "XformScore::read_binary(fname): problem while writing XformScore file " << fname << endl;
			in.close();
			return false;
		}
		in.close();
		return true;
	}
	bool XformScore::read_binary(vector1<string> const & fnames){
		clear();
		Size count = 0;
		TR << "Nfiles: " << fnames.size() << endl;
		for(vector1<string>::const_iterator i = fnames.begin(); i != fnames.end(); ++i){
			if( ++count %50 == 0 ){ TR << "... read " << count << endl; }
			if(!read_binary(*i,false)){
				TR.Error << "XformScore::read_binary(fnames): error reading file "+*i << endl;
				if(!option[basic::options::OptionKeys::mh::ignore_io_errors]()) utility_exit_with_message("XformScore::read_binary(fnames): see tracer output");
			}
		}
		TR << endl;
		return true;
	}



	void XformScore::print_scores(std::ostream & out, std::string const tag) const {
		for(CountMap::const_iterator i = total_count_.begin(); i != total_count_.end(); ++i){
			Key key = i->first;
			// Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
			Count count = i->second;
			out << tag << key << " " << count << endl;
		}
	}

	void XformScore::print_scores_norm(std::ostream & out, std::string const tag) const {
		for(CountMap::const_iterator i = total_count_.begin(); i != total_count_.end(); ++i){
			Key key = i->first;
			Real6 center = hasher_.bin_center_point(hasher_.bin_from_index(key));
			Real count = i->second;
			Real freq = count / sin(radians(center[6]));
			out << tag << key << " " << freq << endl;
		}
	}

	std::ostream & operator<<(std::ostream & out, XformScore const & xh){
		out << "XformScore grid: "  << xh.cart_size_ << " " << xh.cart_resl_ << " " << xh.angle_resl_;
		out << ", " << xh.num_samples() << " xforms, " << xh.num_bins() << " bins, " << (Real)xh.num_samples()/(Real)xh.num_bins() << "/bin Qtile ";
		xh.print_quantiles(out,100);
		return out;
	}


// gather motifs
	bool gather_simplemotifs(
		ResPairMotifs & motifs,
		string const & tag,
		Pose const & pose,
		core::scoring::dssp::Dssp & dssp,
		Size const & ir,
		core::scoring::EnergyEdge const & edge,
		vector1<Real> const & occupancy,
		vector1<Real> const & bfactors,
		vector1<Real> const & rsd_sasa,
		vector1<Real> const & nbrs,
		int & nbad, int & nbadbfac, int & nbadocc, int & nbadiface, int & nbaddist
	){
		using namespace basic::options::OptionKeys;

		Size const jr( edge.get_second_node_ind() );
		Size const ichain = pose.chain(ir);
		Size const jchain = pose.chain(jr);

		if( abs((int)ir-(int)jr) < option[mh::filter::seqsep]() && ichain==jchain )  { ++nbad; ++nbaddist; return false; } // no same or adjcent res
		if(!pose.residue(jr).is_protein())                                           { ++nbad; ++nbadocc;  return false; }
		if( bfactors [ir] < 0.00 || bfactors [ir] > option[mh::filter::coorderr]() ) { ++nbad; ++nbadbfac; return false; }
		if( bfactors [jr] < 0.00 || bfactors [jr] > option[mh::filter::coorderr]() ) { ++nbad; ++nbadbfac; return false; }
		if( occupancy[ir] < 0.99 || occupancy[ir] >      1.01                      ) { ++nbad; ++nbadocc ; return false; }
		if( occupancy[jr] < 0.99 || occupancy[jr] >      1.01                      ) { ++nbad; ++nbadocc ; return false; }

		Real6 rt = get_residue_pair_rt6(pose,ir,jr);
		if( option[mh::filter::maxdist2]() < rt[1]*rt[1]+rt[2]*rt[2]+rt[3]*rt[3] ) { ++nbad; ++nbaddist; return false; }

		Real const fa_atr   = edge[core::scoring::fa_atr     ];
		Real const sasa0    = rsd_sasa[ir]+rsd_sasa[jr];
		Real const hb_sc    = edge[core::scoring::hbond_sc   ];
		Real const hb_bb_sc = edge[core::scoring::hbond_bb_sc];
		Real hb_bb = edge[core::scoring::hbond_sr_bb]+edge[core::scoring::hbond_lr_bb];
		bool is_sspaired = dssp.in_paired_strands(ir,jr) ;

		if(fa_atr > -0.001 )  { ++nbad; ++nbadiface ; return false; }

		Real sasa; {
			core::id::AtomID_Map<Real> atom_sasa0;
			vector1<Real> rsd_sasa0;
			core::id::AtomID_Map<bool> atom_subset0;
			core::pose::initialize_atomid_map( atom_subset0, pose, false );
			for(Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia) atom_subset0[AtomID(ia,ir)] = true;
			for(Size ia = 1; ia <= pose.residue(jr).natoms(); ++ia) atom_subset0[AtomID(ia,jr)] = true;
			sasa = sasa0 / core::scoring::calc_per_atom_sasa( pose, atom_sasa0, rsd_sasa0, 1.7, false, atom_subset0 );
			// cout << sasa0 << " " << core::scoring::calc_per_atom_sasa( pose, atom_sasa0, rsd_sasa0, 1.7, false, atom_subset0 ) << endl;
		}

		ResPairMotif sm(tag,pose,rt,ir,jr,sasa,(nbrs[ir]+nbrs[jr])/2.0,fa_atr,hb_sc,hb_bb_sc,hb_bb,bfactors[ir],bfactors[jr],is_sspaired);

		if( !option[mh::filter::filter_harvest]() || sm.filter() ) motifs.push_back(sm);
		else { ++nbad; ++nbadiface; return false; }

		return true;
	}

	void harvest_motifs(){
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace ObjexxFCL::format;

		if(64!=sizeof(ResPairMotif)) utility_exit_with_message("ResPairMotif should be 64 bytes in size!!!!!!");
		if(12!=sizeof(Xfres)) utility_exit_with_message("Xfres should be 12 bytes in size!!!!!!");
		if(24!=sizeof(Xfrag)) utility_exit_with_message("Xfrag should be 24 bytes in size!!!!!!");

		bool allonefile = false;
		string motif_out_file = basic::options::option[basic::options::OptionKeys::mh::motif_out_file]()+".smtf.bin.gz";
		utility::io::ozstream allout;
		if(basic::options::option[basic::options::OptionKeys::mh::motif_out_file].user()){
			if(utility::file::file_exists(motif_out_file)) utility_exit_with_message("already done!");
			allonefile = true;
			allout.open(motif_out_file,std::ios::out|std::ios::binary);
			if(!allout.good()){
				allout.close();
				utility_exit_with_message("error opening "+motif_out_file );
			}
		}

		vector1<string> const & fnames( option[mh::harvest_motifs]() );
		for(int ifile = 1; ifile <= (int)fnames.size(); ++ifile ){
			ResPairMotifs motifs;
			string fname = fnames[ifile];
			string const tag = tag_from_pdb_fname(fname);
			if(""==tag){
				TR.Error << "WARNING: skipping non-pdb fname: " << fname << endl;
				continue;
			}

			if(!allonefile && utility::file::file_exists(option[out::file::o]()+"/"+tag.substr(1,2)+"/"+tag+".smtf.bin.gz")) continue;

			Pose pose; vector1<Real> bfactors,occupancy;
			utility::vector1<int> pdbres;
			std::map<int,char> pdbchain;
			int nresmodel1;
			if( !protocols::sic_dock::read_biounit(fname,pose,bfactors,occupancy,pdbres,pdbchain,nresmodel1,99999,false) ){
				TR.Error << "FAIL TO READ " << fname << endl;
				continue;
			}
			ScoreFunctionOP sf = core::scoring::get_score_function();
			core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
			myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
			sf->set_energy_method_options(myopt);
			sf->score(pose);

			core::scoring::dssp::Dssp dssp(pose);
			dssp.insert_edge_ss_into_pose(pose);

			Energies    const & energies     ( pose.energies() );
			EnergyGraph const & energy_graph ( energies.energy_graph() );
			vector1<Real> rsd_sasa = get_sasa(pose,1.7);
			vector1<Real> nbrs = get_nbrs(pose);

			int ssbdy = option[mh::harvest_motifs_min_hh_ends]();
			int min_boundary = ssbdy;

			// scan pairs
			int nbad=0,nbadbfac=0,nbadocc=0,nbadiface=0,nbaddist=0,nbadsselem=0;
			for(int ir = 1+min_boundary; ir <= nresmodel1-min_boundary; ++ir){
				if(!pose.residue(ir).is_protein()) { ++nbad; continue; }
				for ( core::graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(ir)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(ir)->const_upper_edge_list_end();
						iru != irue; ++iru ) {
					EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
					if(!gather_simplemotifs(motifs,tag,pose,dssp,ir,edge,occupancy,bfactors,rsd_sasa,nbrs,nbad,nbadbfac,nbadocc,nbadiface,nbaddist)) continue;
				}
			}

			// do I/O
			TR << fname <<' '<< motifs.size() <<" nbad: "<< I(5,nbad) <<" nbadbfac: "<< I(5,nbadbfac) <<" nbadocc: "<< I(5,nbadocc) <<" nbadiface: "<< I(5,nbadiface) <<" nbaddist: "<< I(5,nbaddist) <<" nbadsselem: "<<I(5,nbadsselem)<< endl;
			if(allonefile){
				TR << "appending " << motifs.size() << " motifs to " << motif_out_file+".smtf.bin.gz" << endl;
				if(!write_motifs_binary(allout,motifs)) utility_exit_with_message("error writing to file "+motif_out_file+".smtf.bin.gz");
			} else {
				string outfile = option[out::file::o]()+"/"+tag.substr(1,2)+"/"+tag+".smtf.bin.gz";
				TR << "dumping " << motifs.size() << " motifs to " << outfile << endl;
				write_motifs_binary(outfile,motifs);
			}

		}
		// DONE: ;
		if(allonefile){
			if(!allout.good()){
				allout.close();
				utility_exit_with_message("error after writing to "+motif_out_file+".smtf.bin.gz" );
			}
			allout.close();
		}
		TR << "DONE HARVESTING" << endl;
	}


// manager
	MotifHashManager * MotifHashManager::instance_( 0 );
	MotifHashManager * MotifHashManager::get_instance() {
		if ( instance_ == 0 ) instance_ = new MotifHashManager();
		return instance_;
	}

	MotifHashCOP MotifHashManager::motif_hash_from_cli(){
		if( cli_motif_hash_ == 0 ){
			ResPairMotifs motifs; load_motifs( option[basic::options::OptionKeys::mh::input_motifs](), motifs );
			TR << "create motif hash" << endl;
			cli_motif_hash_ = MotifHashOP( new MotifHash(motifs) );
			cli_motif_hash_->hasher().tree_init(5);
		}
		return cli_motif_hash_;
	}

	XformScoreOP get_xform_score_from_file(XformScoreOP xs,vector1<string> const & datfiles,XformScoreOP defaultval=NULL){
		if( xs == 0 ){
			if(datfiles.size()){
				xs = XformScoreOP( new XformScore );
				if(!xs->read_binary(datfiles)) utility_exit_with_message("bad file -mh:score_data");
				TR << "created xs from " << datfiles.size() << ":" << endl;
				BOOST_FOREACH(string s,datfiles) TR << s << endl;
				// cout <<"MotifHashManager: " << xs->cart_size_ << " " << xs->cart_resl_ << " " << xs->angle_resl_ << " " << endl;
			 } else {
				xs = defaultval;
			 }
		}
		return xs;
	}

	XformScoreCOP MotifHashManager::xform_score_from_cli()       { return get_xform_score_from_file(cli_xform_score_       ,option[basic::options::OptionKeys::mh::xform_score_data        ]());  }
	XformScoreCOP MotifHashManager::xform_score_ee_from_cli()    { return get_xform_score_from_file(cli_xform_score_ee_    ,option[basic::options::OptionKeys::mh::xform_score_data_ee     ](),cli_xform_score_);  }
	XformScoreCOP MotifHashManager::xform_score_eh_from_cli()    { return get_xform_score_from_file(cli_xform_score_eh_    ,option[basic::options::OptionKeys::mh::xform_score_data_eh     ](),cli_xform_score_);  }
	XformScoreCOP MotifHashManager::xform_score_he_from_cli()    { return get_xform_score_from_file(cli_xform_score_he_    ,option[basic::options::OptionKeys::mh::xform_score_data_he     ](),cli_xform_score_);  }
	XformScoreCOP MotifHashManager::xform_score_hh_from_cli()    { return get_xform_score_from_file(cli_xform_score_hh_    ,option[basic::options::OptionKeys::mh::xform_score_data_hh     ](),cli_xform_score_);  }
	XformScoreCOP MotifHashManager::xform_score_sspair_from_cli(){ return get_xform_score_from_file(cli_xform_score_sspair_,option[basic::options::OptionKeys::mh::xform_score_data_sspair]());  }


	bool MotifHit::operator<(MotifHit const & other) const{
		// !!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!
		return memcmp(this,&other,sizeof(MotifHit)-16) < 0;

	}
	bool MotifHit::operator==(MotifHit const & other) const{
		// !!!!!!!!!!!!!!!!!!! HERE BE DRAGONS !!!!!!!!!!!!!!!!!!!!!!
		return memcmp(this,&other,sizeof(MotifHit)-16) == 0;
	}

	void MotifHits::filter_redundant(){
		std::set<MotifHit> hitset(begin(),end());
		clear();
		std::copy(hitset.begin(),hitset.end(),std::back_insert_iterator<MotifHits>(*this));
	}
	void MotifHits::compute_metrics() const {
		std::set<int> resset;
		std::set<std::pair<int,int> > pairset;
		total_score = 0.0;
		for(MotifHits::const_iterator i = begin(); i != end(); ++i){
			resset.insert(i->residue1);
			resset.insert(i->residue2);
			pairset.insert(std::make_pair(i->residue1,i->residue2));
			total_score += i->score;
		}
		rescover = resset.size();
		paircover = pairset.size();
	}

	std::ostream & operator<<(std::ostream & out, MotifHit  const & h){
		out << I(4,h.residue1) << " " << I(4,h.residue2) << " " << F(7,3,h.score) << " " << h.motif;
		return out;
	}
	std::ostream & operator<<(std::ostream & out, MotifHits const & h){
		h.compute_metrics();
		out << "MotifHits" << " score: " << h.total_score  << ", num: " << h.size() << ", rescover: " << h.rescover << ", paircover: " << h.paircover;
		for(MotifHits::const_iterator i = h.begin(); i != h.end(); ++i){
			out << endl << "    HIT " << *i;
		}
		return out;
	}


// utils for MotifHashRigidScore

	int
	MotifHash::get_matching_motifs(
		Pose const & pose1,
		Pose const & pose2,
		Bools const & useres1,
		Bools const & useres2,
		MotifHits & hits,
		ClashChecker clash_check,
		Real radius
	) const {
		using namespace protocols::sic_dock;

		xyzStripeHashPoseCOP ccheck1bb32,ccheck2bb32,ccheck1nco2,ccheck2nco2;
		if(clash_check){
			ccheck1bb32 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose1,BB ,3.2) );
			ccheck2bb32 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose2,BB ,3.2) );
			ccheck1nco2 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose1,NCO,2.5) );
			ccheck2nco2 = xyzStripeHashPoseCOP( new xyzStripeHashPose(pose2,NCO,2.5) );
		}

		vector1<Xform> stubs1,stubs2;
		for(Size ir = 1; ir <= pose1.total_residue(); ++ir) stubs1.push_back(Xform(pose1.xyz(AtomID(2,ir)),pose1.xyz(AtomID(1,ir)),pose1.xyz(AtomID(2,ir)),pose1.xyz(AtomID(3,ir))));
		for(Size ir = 1; ir <= pose2.total_residue(); ++ir) stubs2.push_back(Xform(pose2.xyz(AtomID(2,ir)),pose2.xyz(AtomID(1,ir)),pose2.xyz(AtomID(2,ir)),pose2.xyz(AtomID(3,ir))));
		Pose motifpose;
		int nbadmotifs = 0;
		for(Size ir = 1; ir <= pose1.total_residue(); ++ir){
			if(!useres1[ir]) continue;
			char /*aa1 = pose1.residue(ir).name1(),*/ ss1=pose1.secstruct(ir);
			for(Size jr = 1; jr <= pose2.total_residue(); ++jr){
				if(!useres2[jr]) continue;
				char /*aa2 = pose2.residue(jr).name1(),*/ ss2=pose2.secstruct(jr);
				float dsq = pose1.xyz(AtomID(2,ir)).distance_squared(pose2.xyz(AtomID(2,jr)));
				if(dsq > cart_size()*cart_size()) continue;
				Real6 rt = get_residue_pair_rt6(pose1,ir,pose2,jr);
				ResPairMotifs motifs;
				if(radius>0.0) find_motifs_with_radius(rt,radius,motifs);
				else           find_motifs(rt,motifs);
				for(ResPairMotifs::const_iterator i = motifs.begin(); i != motifs.end(); ++i){
					if( ss1!='L' && ss2!='L' && ss1!=i->ss1() ) continue;
					if( ss2!='L' && ss1!='L' && ss2!=i->ss2() ) continue;
					i->fill_pose_with_motif(motifpose);
					align_motif_pose_super(motifpose,pose1,ir,pose2,jr);

					// cout << ir << " " << jr << " " << *i << endl;

					if(clash_check){
						bool badmotif = false;
						Size nstart1intra=6,nstart2intra=6,nstart1inter=6,nstart2inter=6;
						if(i->aa1()=='P') nstart1intra = 999;
						if(i->aa2()=='P') nstart2intra = 999;
						for(Size ia = nstart1inter; ia <= motifpose.residue(1).nheavyatoms(); ++ia) if(ccheck2bb32->clash(motifpose.xyz(AtomID(ia,1)))) badmotif = true;
						for(Size ia = nstart1intra; ia <= motifpose.residue(1).nheavyatoms(); ++ia) if(ccheck1nco2->clash(motifpose.xyz(AtomID(ia,1)))) badmotif = true;
						for(Size ia = nstart2inter; ia <= motifpose.residue(2).nheavyatoms(); ++ia) if(ccheck1bb32->clash(motifpose.xyz(AtomID(ia,2)))) badmotif = true;
						for(Size ia = nstart2intra; ia <= motifpose.residue(2).nheavyatoms(); ++ia) if(ccheck2nco2->clash(motifpose.xyz(AtomID(ia,2)))) badmotif = true;
						for(Size ia = nstart1intra; ia <= motifpose.residue(1).nheavyatoms(); ++ia) if(clash_check->clash(motifpose.xyz(AtomID(ia,1)))) badmotif = true;
						for(Size ia = nstart2intra; ia <= motifpose.residue(2).nheavyatoms(); ++ia) if(clash_check->clash(motifpose.xyz(AtomID(ia,2)))) badmotif = true;
						if(badmotif) {
							++nbadmotifs;
							continue;
						}
					}
					// motifpose.dump_pdb("test1.pdb");
					// pose1.dump_pdb("pose1.pdb");
					// pose2.dump_pdb("pose2.pdb");
					// utility_exit_with_message("foo");

					hits.push_back(MotifHit(ir,jr,0,*i));
				}
				// cout << ir <<" "<< jr <<" "<< hits.size() << endl;
			}
		}
		// cout << "NBADMOTIFS: " <<  nbadmotifs << " " << hits.size() << endl;
		return hits.size();
	}

	int
	MotifHash::get_matching_motifs( Pose const & pose1, Pose const & pose2, MotifHits & hits, ClashChecker clash_check, Real radius ) const {
		Bools useres1(pose1.n_residue(),true);
		Bools useres2(pose2.n_residue(),true);
		return get_matching_motifs(pose1,pose2,useres1,useres2,hits,clash_check,radius);
	}

	int
	MotifHash::get_matching_motifs( Pose const & pose, core::pack::task::PackerTaskCOP ptask, MotifHits & hits, ClashChecker clash_check, Real radius ) const {
		Bools useres = ptask->designing_residues();
		return get_matching_motifs(pose,pose,useres,useres,hits,clash_check,radius);
	}


	int
	base_dump_matching_motifs( Pose const & pose1, Pose const & pose2, std::ostream & out, int & count, bool print, MotifHits const & hits ) {
		for(MotifHits::const_iterator i = hits.begin(); i != hits.end(); ++i){
			i->motif.dump_aligned_motif(out,pose1,i->residue1,pose2,i->residue2,++count);
		}
		if(print) cout << hits << endl;
		return hits.size();
	}

	int
	MotifHash::dump_matching_motifs( Pose const & pose1, Pose const & pose2, Real radius, std::ostream & out, int & count, ClashChecker clash_check, bool print ) const {
		MotifHits hits;
		get_matching_motifs(pose1,pose2,hits,clash_check,radius);
		return base_dump_matching_motifs(pose1,pose2,out,count,print,hits);
	}

	int
	MotifHash::dump_matching_motifs( Pose const & pose1, Pose const & pose2, std::ostream & out, int & count, ClashChecker clash_check, bool print ) const {
		MotifHits hits;
		get_matching_motifs(pose1,pose2,hits,clash_check);
		return base_dump_matching_motifs(pose1,pose2,out,count,print,hits);
	}
	int
	MotifHash::print_matching_motifs( Pose const & pose1, Pose const & pose2, std::ostream & out, ClashChecker clash_check ) const {
		MotifHits hits;
		get_matching_motifs(pose1,pose2,hits,clash_check);
		out << hits << endl;
		return hits.size();
	}

	int
	MotifHash::stat_matching_motifs( Pose const & pose1, Pose const & pose2, std::map<std::string,Real> & stats, ClashChecker clash_check, Real radius) const {
		using std::map;
		using std::pair;
		using std::make_pair;
		Real hh=0,he=0,hl=0,ee=0,el=0,ll=0,pp=0,totscore=0,totsqstscore=0,count=0;
		vector1<Real> scores;
		map<pair<int,int>,Real> pairscores;
		map<int,Real> resscores;
		MotifHits hits;
		get_matching_motifs(pose1,pose2,hits,clash_check,radius);
		BOOST_FOREACH(MotifHit & hit,hits){
			Size ir = hit.residue1;
			Size jr = hit.residue2;
			char aa1 = pose1.residue(ir).name1(), ss1=pose1.secstruct(ir);
			char aa2 = pose2.residue(jr).name1(), ss2=pose2.secstruct(jr);
			Real6 rt = get_residue_pair_rt6(pose1,ir,pose2,jr);
			Real score=1.0;
			ResPairMotif const & sm(hit.motif);
			Real err = rt6_rt6_bb_dis2(rt,sm.rt());
			score *= exp( -2.0 * err );
			if( sm.aa1()==aa1 && sm.aa2()==aa2 ) score *= 1.00;
			else                                 score *= 0.50;
			if( sm.ss1()==ss1 && sm.ss2()==ss2 ) score *= 1.00;
			else                                 score *= 0.30;
			totscore += score;
			totsqstscore += sqrt(score);
			scores.push_back(score);
			pairscores[make_pair(ir,jr)] += score;
			resscores [ir        ] += score;
			resscores [jr+1000000] += score;
			count += 1.0;
					 if( (sm.dssp1()=='P' && sm.dssp2()=='P') ) pp += 1.0;
			else if( (sm.ss1()=='H' && sm.ss2()=='H') || (sm.ss1()=='H' && sm.ss2()=='H') ) hh += 1.0;
			else if( (sm.ss1()=='H' && sm.ss2()=='E') || (sm.ss1()=='E' && sm.ss2()=='H') ) he += 1.0;
			else if( (sm.ss1()=='H' && sm.ss2()=='L') || (sm.ss1()=='L' && sm.ss2()=='H') ) hl += 1.0;
			else if( (sm.ss1()=='E' && sm.ss2()=='E') || (sm.ss1()=='E' && sm.ss2()=='E') ) ee += 1.0;
			else if( (sm.ss1()=='E' && sm.ss2()=='L') || (sm.ss1()=='L' && sm.ss2()=='E') ) el += 1.0;
			else if( (sm.ss1()=='L' && sm.ss2()=='L') || (sm.ss1()=='L' && sm.ss2()=='L') ) ll += 1.0;
			else TR.Error << "BAD SS IN MOTIF " << sm.ss1() << " " << sm.ss2() << endl;
		}


		if( stats.find("M_NUM"    )==stats.end()) stats["M_NUM"    ] = 0;
		// if( stats.find("M_xpsc"   )==stats.end()) stats["M_xpsc"   ] = 0;
		// if( stats.find("M_xpscrt" )==stats.end()) stats["M_xpscrt" ] = 0;
		// if( stats.find("M_PP"     )==stats.end()) stats["M_PP"     ] = 0;
		if( stats.find("M_HH"     )==stats.end()) stats["M_HH"     ] = 0;
		if( stats.find("M_HE"     )==stats.end()) stats["M_HE"     ] = 0;
		// if( stats.find("M_HL"     )==stats.end()) stats["M_HL"     ] = 0;
		if( stats.find("M_EE"     )==stats.end()) stats["M_EE"     ] = 0;
		// if( stats.find("M_EL"     )==stats.end()) stats["M_EL"     ] = 0;
		// if( stats.find("M_LL"     )==stats.end()) stats["M_LL"     ] = 0;

		stats["M_NUM"     ] += count;
		// stats["M_xpsc"    ] += totscore;
		// stats["M_xpscrt" ] += totsqstscore;
		// stats["M_PP" ] += pp;
		stats["M_HH"      ] += hh;
		stats["M_HE"      ] += he;
		// stats["M_HL"      ] += hl;
		stats["M_EE"      ] += ee;
		// stats["M_EL"      ] += el;
		// stats["M_LL"      ] += ll;

		// cout <<endl<< "stat_matching_motifs " << count << " " << hits.size() << " " << stats["M_NUM"     ] << endl<<endl;

		return hits.size();
	}



// packing

 MotifRotamerSetOperationCOP
 MotifHash::get_matching_motifs_rotsetop(
	Pose const & pose,
	core::pack::task::PackerTaskCOP ptask,
	Real radius,
	ClashChecker clash_check
 ) const {
	if(!ptask) ptask = core::pack::task::TaskFactory::create_packer_task(pose);
	// cout << *ptask << endl;
	MotifHits hits;
	get_matching_motifs(pose,ptask,hits,clash_check,radius);
	return MotifRotamerSetOperationCOP( new MotifRotamerSetOperation(pose,hits) );
 }

 MotifRotamerSetOperation::MotifRotamerSetOperation(
	Pose const & refpose,
	MotifHits const & motifs
 ):
	motifs_(motifs)
 {
	BOOST_FOREACH(MotifHit const & h, motifs_){
		core::pose::PoseOP pose( new Pose );
		h.motif.fill_pose_with_motif(*pose);
		align_motif_pose_super(*pose,refpose,h.residue1,refpose,h.residue2);
		if(h.residue1>(int)res1_poses_.size()) res1_poses_.resize(h.residue1);
		if(h.residue2>(int)res2_poses_.size()) res2_poses_.resize(h.residue2);
		res1_poses_[h.residue1].push_back(pose);
		res2_poses_[h.residue2].push_back(pose);
	}
 }

 void
 MotifRotamerSetOperation::alter_rotamer_set(
	core::pose::Pose const & /*pose*/,
	core::scoring::ScoreFunction const & /*sfxn*/,
	core::pack::task::PackerTask const & /*ptask*/,
	core::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
 ){
	BOOST_FOREACH(PoseCOP pp, res1_poses_[rotamer_set.resid()]){
		rotamer_set.add_rotamer(pp->residue(1));
	}
	BOOST_FOREACH(PoseCOP pp, res2_poses_[rotamer_set.resid()]){
		rotamer_set.add_rotamer(pp->residue(2));
	}
 }

 void
 MotifRotamerSetOperation::dump_pdb(std::string const & fname) const {
	utility::io::ozstream out(fname);
	Size ano = 0;
	Size nres = res1_poses_.size();
	for(Size ir = 1; ir <= nres; ++ir){
		BOOST_FOREACH(PoseCOP pose, res1_poses_[ir]){
			out << "MODEL" << endl;
			core::io::pdb::dump_pdb_residue(pose->residue(1),ano,out);
			core::io::pdb::dump_pdb_residue(pose->residue(2),ano,out);
			out << "ENDMDL" << endl;
		}
		BOOST_FOREACH(PoseCOP pose, res2_poses_[ir]){
			out << "MODEL" << endl;
			core::io::pdb::dump_pdb_residue(pose->residue(1),ano,out);
			core::io::pdb::dump_pdb_residue(pose->residue(2),ano,out);
			out << "ENDMDL" << endl;
		}
	}
	out.close();
 }

 Real
 MotifRotamerSetOperation::increase_packer_residue_radius(
	core::pose::Pose const & /*pose*/,
	core::pack::task::PackerTaskCOP /*the_task*/,
	core::Size /*residue_in*/
 ) const {
	return 0.0;
 }




}
}
