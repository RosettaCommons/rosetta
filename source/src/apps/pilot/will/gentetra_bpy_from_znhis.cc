// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/genmatch.cc
/// @brief ???

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/smhybrid.OptionKeys.gen.hh>
#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/pack/optimizeH.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <sstream>
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// #include <devel/init.hh>

// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/pose/util.tmpl.hh>
#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>


using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static basic::Tracer TR("gentetra");
static core::io::silent::SilentFileData sfd;


inline Real sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	if( sqdist > 9.*9. ) {
		return 0.0;
	} else if( sqdist < 6.*6. ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0  - sqr( (dist - 6.) / (9. - 6.) ) );
	}
}


Real iface_check_c3(Pose & pose, Size nres, vector1<Size> const & iface_candidates) {
	Real num = 0;
	for(vector1<Size>::const_iterator i=iface_candidates.begin(),ie=iface_candidates.end(); i != ie; ++i) {
		for(vector1<Size>::const_iterator j=iface_candidates.begin(),je=iface_candidates.end(); j != je; ++j) {
			num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+1*nres).xyz(5)));
			num += sigmoidish_neighbor(pose.residue(*i).xyz(5).distance_squared(pose.residue(*j+2*nres).xyz(5)));
		}
	}
	return num;
}

vector1<Size> read_res_list(string fn) {
	vector1<Size> l;
	if(fn=="") return l;
	if(fn=="_") return l;
	if(fn.size()==1 && fn[0]==(char)0) return l;
	izstream in(fn);
	if(!in.good()) {
		utility_exit_with_message("can't open res list file '"+fn+"'");
	}
	Size r;
	while( in >> r ) l.push_back(r);
	return l;
}

void repack(Pose & pose, Size nres, ScoreFunctionOP sf) {
	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	task->initialize_extra_rotamer_flags_from_command_line();
	for(Size i = 1; i <= nres; ++i) {
		if(pose.residue(i).name3()=="BPY") {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
		}
	}
	// TR << *task << std::endl;
	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}

void design(Pose & pose, Size nres, ScoreFunctionOP sf) {
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		atom_map.set(AtomID(2,ir) , true );
		atom_map.set(AtomID(3,ir) , true );
		atom_map.set(AtomID(5,ir) , true );
	}
	core::id::AtomID_Map<Real> atom_sasa; utility::vector1<Real> sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 2.3, false, atom_map );
	for(Size i = 1; i <= sasa.size(); ++i) if( atom_sasa.n_atom(i) > 4 ) sasa[i] = atom_sasa[AtomID(5,i)];

	using namespace core::pack::task;
	PackerTaskOP task = TaskFactory::create_packer_task(pose);
	vector1< bool > aas(20,true);
	aas[core::chemical::aa_cys] = false;
	aas[core::chemical::aa_his] = false;
	// aas[core::chemical::aa_met] = false;
	aas[core::chemical::aa_pro] = false;
	aas[core::chemical::aa_gly] = false;
	aas[core::chemical::aa_gly] = false;
	if(option[basic::options::OptionKeys::willmatch::exclude_ala]()) aas[core::chemical::aa_ala] = false;
	if(option[basic::options::OptionKeys::smhybrid::design_hydrophobic]()) {
		aas[core::chemical::aa_ser] = false;
		aas[core::chemical::aa_thr] = false;
		aas[core::chemical::aa_asp] = false;
		aas[core::chemical::aa_glu] = false;
		aas[core::chemical::aa_lys] = false;
		aas[core::chemical::aa_arg] = false;
		aas[core::chemical::aa_asn] = false;
		aas[core::chemical::aa_gln] = false;
	}

	vector1<Size> fixed;
	if(option[basic::options::OptionKeys::willmatch::fixed_res].user()) {
		utility::io::izstream in(option[basic::options::OptionKeys::willmatch::fixed_res]());
		Size tmp;
		while(in>>tmp) fixed.push_back(tmp);
		in.close();
	}

	vector1<Size> interface;
	if(option[basic::options::OptionKeys::willmatch::design_interface]()) {
		for(Size i = 1; i <= nres; ++i) {
			if( sasa.size() >= i && sasa[i] > 15.0 ) continue;
			AtomID aid(5,i);
			if(pose.residue(i).nheavyatoms() < 5) aid.atomno() = 2;
			for(Size j = nres+1; j <= 3*nres; ++j) {
				AtomID aid2(5,j);
				if(pose.residue(j).nheavyatoms() < 5) aid.atomno() = 2;
				if(pose.xyz(aid).distance_squared(pose.xyz(aid2)) < 49) {
					interface.push_back(i);
				}
			}
		}
		for(Size i = 1; i <= nres; ++i) {
			Vec xyz = pose.xyz(AtomID(5,i));
			xyz.z() = 0;
			if( xyz.length() < 8.0 ) interface.push_back(i);
		}
	}
	for(Size i = 1; i <= nres; ++i) {
		if(pose.residue(i).name3()=="BPY") {
			task->nonconst_residue_task(i).prevent_repacking();
			// std::exit(-1);
		} else if(std::find(fixed.begin(),fixed.end(),i)!=fixed.end()){
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		} else if(std::find(interface.begin(),interface.end(),i)!=interface.end()){
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(aas);
			task->nonconst_residue_task(i).initialize_extra_rotamer_flags_from_command_line();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
			task->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_ONE_STDDEV );
			task->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_ONE_STDDEV );
		}
	}
	for(Size i = nres+1; i <= pose.n_residue(); ++i) {
		task->nonconst_residue_task(i).prevent_repacking();
	}
	TR << *task << std::endl;

	protocols::simple_moves::symmetry::SymPackRotamersMover repack( sf, task );
	repack.apply(pose);
}

vector1<Vec> line_cone_intersection(Vec p, Vec d, Vec v, Vec a, Real t) {
	vector1<Vec> sol;
	t = numeric::conversions::radians(t);
	Mat M = numeric::outer_product(a,a) - cos(t)*cos(t)*Mat::identity();
	Vec D = p-v;
	Real c2 =   d.dot(M*d);
	Real c1 = 2*d.dot(M*D);
	Real c0 =   D.dot(M*D);
	Real disc = c1*c1 - 4*c0*c2;
	if( disc == 0) sol.push_back( p + (-c1)/(2.0*c2)*d );
	else if( disc > 0) {
		disc = sqrt(disc);
		sol.push_back(p+(-c1+disc)/(2.0*c2)*d);
		sol.push_back(p+(-c1-disc)/(2.0*c2)*d);
	}
	return sol;
}

// Vec projperp(Vec u, Vec v) {
// 	return v - projection_matrix(u)*v;
// }

vector1<std::pair<Vec,Vec> > intersecting_bpy_axes(Vec CB, Vec CG, Vec FE, Vec symmaxis, Vec symmcen = Vec(0,0,0)) {
	vector1<std::pair<Vec,Vec> > sol;
	Vec  p = symmcen;
	Vec  d = symmaxis.normalized();
	Vec  a = (CG-CB).normalized();
	Vec  v = CG + projection_matrix(a)*(FE-CG);
	Real t = 35.2643434495;
	Real l = v.distance(FE);
	vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
	for(Size i = 1; i <= X.size(); ++i) {
		Vec x = X[i];
		Vec o = projperp(a,x-v);
		Real L = o.length();
		o = o.normalized() * l;
		Real ang = 90.0 - numeric::conversions::degrees( atan(l/L) );
		Vec o1 = rotation_matrix_degrees(a, ang) * o;
		Vec o2 = rotation_matrix_degrees(a,-ang) * o;
		sol.push_back(std::pair<Vec,Vec>(x,v+o1));
		sol.push_back(std::pair<Vec,Vec>(x,v+o2));
	}
	return sol;
}



void minimize(Pose & pose, Size nres, Size , ScoreFunctionOP sf, int bb=0) {
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	// core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
	movemap->set_chi(true);
	movemap->set_bb(false);
	movemap->set_jump(false);
	if(bb==1) for(Size i = 1; i <= nres; ++i) if(pose.secstruct(i)=='L') movemap->set_bb(i,true);
	if(bb>=2) for(Size i = 1; i <= nres; ++i) movemap->set_bb(i,true);
	// movemap->set_chi(bpyres,false);

	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );

	m.apply(pose);

}

std::pair<Size,Size> makesplitwork_3bpy(Size total) {
	using namespace basic::options::OptionKeys;
	Size size1 = 1;
	Size size2 = total;
	if( option[willmatch::splitwork].user() ) {
		Size part   = option[willmatch::splitwork]()[1];
		Size nparts = option[willmatch::splitwork]()[2];
		size1 = (Size)((part-1)*std::ceil(((Real)total)/(Real)nparts))+1;
		size2 = (Size)((part  )*std::ceil(((Real)total)/(Real)nparts));
		if( option[in::file::s]().size() == 1 ) {
			Real frac1 = ((Real)part-1)/(Real)nparts;
			Real frac2 = ((Real)part  )/(Real)nparts;
			frac1 = 1.0 - sqrt(1.0-frac1);
			frac2 = 1.0 - sqrt(1.0-frac2);
			size1 = (Size)std::ceil(frac1*(Real)total)+1;
			size2 = (Size)std::ceil(frac2*(Real)total);
		}
	}
	TR << "SIZE " << size1 << " " << size2 << std::endl;
	return std::pair<Size,Size>(size1,size2);
}

vector1<Vec> get_zn_axes(core::pose::Pose const & pose) {
	Vec nd1 = pose.residue(1).xyz("ND1");
	Vec ne1 = pose.residue(1).xyz("NE2");
	Vec nd2 = pose.residue(2).xyz("ND1");
	Vec ne2 = pose.residue(2).xyz("NE2");
	Vec zn  = pose.residue(3).xyz("ZN");
	Vec n1 = ( nd1.distance(zn) < ne1.distance(zn) ) ? nd1 : ne1;
	Vec n2 = ( nd2.distance(zn) < ne2.distance(zn) ) ? nd2 : ne2;
	Vec z = (zn - (n1+n2)/2.0).normalized();
	Vec y = projperp(z,n2-n1).normalized();
	Vec x = y.cross(z).normalized();
	vector1<Vec> axes;
	axes.push_back((x+y).normalized()); // 45 degrees off to make tetrahederal site
	axes.push_back((x-y).normalized());
	return axes;
}

void fixH(core::pose::Pose & pose) {
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(!pose.residue(i).has("H")) continue;
		numeric::xyzVector<Real> n  = pose.residue(i).xyz("N");
		numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
		Size in = i-1;
		if(in == 0) in = pose.n_residue();
		numeric::xyzVector<Real> c  = pose.residue(in).xyz("C");
		numeric::xyzVector<Real> h  = n + (n-(ca+c)/2.0).normalized()*1.01;
		pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
	}
}


void run() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ScoreFunctionOP sf = core::scoring::get_score_function();
	ScoreFunctionOP sfrep = new core::scoring::ScoreFunction;
	sfrep->set_weight(core::scoring::fa_rep,1.0);

	core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	Pose bpy;
	core::import_pose::pose_from_pdb(bpy ,*rs,"input/bpy_ideal.pdb");
	core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);

	Pose nat;
	core::import_pose::pose_from_pdb(nat,*rs,option[in::file::s]()[1]);
	std::string infile = utility::file_basename(option[in::file::s]()[1]);
	vector1<Pose> m;
	core::import_pose::pose_from_pdb(m,*rs,option[in::file::s]()[2]);
	Pose base(nat);
	for(Size ilig = 2; ilig <= m.size(); ++ilig) {
		Pose lig(m[ilig]);
		core::pose::remove_lower_terminus_type_from_pose_residue(lig,1);
		core::pose::remove_upper_terminus_type_from_pose_residue(lig,1);
		core::pose::remove_lower_terminus_type_from_pose_residue(lig,2);
		core::pose::remove_upper_terminus_type_from_pose_residue(lig,2);
		vector1<Vec> znaxes = get_zn_axes(lig);
		for(Size iaxs = 1; iaxs <= znaxes.size(); ++iaxs) {
			Vec a2f1 = znaxes[iaxs];
			Vec c2f1 = lig.residue(3).xyz("ZN");
			Mat Rzn = rotation_matrix_degrees(a2f1,180.0);
			Size ihis = lig.pdb_info()->number(1);
			Size jhis = lig.pdb_info()->number(2);
			// clash check
			for(Size ir = 1; ir <= base.n_residue(); ++ir) {
				for(Size ia = 1; ia <= base.residue(ir).nheavyatoms(); ia++) {
					Vec xi = Rzn*(base.xyz(AtomID(ia,ir))-c2f1)+c2f1;
					for(Size ja = 6; ja <= lig.residue(1).nheavyatoms(); ja++) {
						if( xi.distance_squared(lig.xyz(AtomID(ja,1))) < 9.0 ) goto clash1;
						if( xi.distance_squared(lig.xyz(AtomID(ja,2))) < 9.0 ) goto clash1;
					}
					for(Size jr = 1; jr <= base.n_residue(); ++jr) {
						for(Size ja = 1; ja <= base.residue(jr).nheavyatoms(); ja++) {
							if( base.xyz(AtomID(ja,jr)).distance_squared( xi ) < 9.0 ) goto clash1;
						}
					}
				}
			}
			goto noclash1;	clash1: continue; noclash1:
			TR << "============================= found c2f1 dimer ilig " << ilig << " =============================" << std::endl;
			base.replace_residue(ihis,lig.residue(1),true);
			base.replace_residue(jhis,lig.residue(2),true);
			fixH(base);
			// base.dump_pdb("test1.pdb");
			// rot_pose(base,Rzn,c2f1);
			// base.dump_pdb("test2.pdb");

			for(Size ibpy = 1; ibpy <= base.n_residue(); ++ibpy) {
				if(ibpy == ihis || ibpy == jhis ) continue;
				base = nat;
				base.replace_residue(ihis,lig.residue(1),true);
				base.replace_residue(jhis,lig.residue(2),true);
				core::pose::replace_pose_residue_copying_existing_coordinates(base,ibpy,rs->name_map("BPY"));
				Real chi1_incr = option[willmatch::chi1_increment]();
				for(Real bch1 = 0.0; bch1 <= 360; bch1 += chi1_incr) {
					base.set_chi(1,ibpy,bch1);
					for(Size ir = 1; ir <= base.n_residue(); ++ir) {
						Size natom = (ir==ibpy) ? 5 : base.residue(ir).nheavyatoms();
						for(Size ia = 1; ia <= natom; ia++) {
							if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CZ")) < 9.0 ) goto clash2;
							if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CP")) < 9.0 ) goto clash2;
							if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz("CM")) < 9.0 ) goto clash2;
						}
					}
					goto noclash2;	clash2: continue; noclash2:
					Vec CB = base.residue(ibpy).xyz("CB");
					Vec CG = base.residue(ibpy).xyz("CG");
					Vec FE = base.residue(ibpy).xyz("ZN");
					vector1<std::pair<Vec,Vec> > baxes = intersecting_bpy_axes(CB,CG,FE,a2f1,c2f1);
					for(Size jaxs = 1; jaxs <= baxes.size(); ++jaxs) {
						Vec isct = baxes[jaxs].first;
						Vec c3f1 = baxes[jaxs].second;
						if( isct.distance_squared(c3f1) < 25.0 ) continue;
						if( isct.distance_squared(c2f1) < 25.0 ) continue;
						Vec a3f1 = (isct-c3f1).normalized();
						Real orig_ang = angle_degrees( c2f1, isct, c3f1);
						Real ang = (orig_ang > 90.0) ? 180.0-orig_ang : orig_ang;
						if( fabs(ang-54.7356563997) > 1.0 ) continue;
						// aln bpy
						Vec v = CG + projection_matrix(CB-CG)*(FE-CG);
						ang = dihedral_degrees( c3f1,v,CG,FE );
						base.set_chi(2,ibpy, base.chi(2,ibpy) + ang );
						for(Size ir = 1; ir <= base.n_residue(); ++ir) {
							if(ir==ibpy) continue;
							for(Size ia = 1; ia <= base.residue(ir).nheavyatoms(); ia++) {
								for(Size ja = 7; ja <= base.residue(ibpy).nheavyatoms(); ja++) {
									if( base.xyz(AtomID(ia,ir)).distance_squared(base.xyz(AtomID(ja,ibpy))) < 9.0 ) goto clash3;
								}
							}
						}
						goto noclash3;	clash3: continue; noclash3:

						Mat R1 = rotation_matrix_degrees(a3f1,120.0);
						Mat R2 = rotation_matrix_degrees(a3f1,240.0);
						for(Size ir = 1; ir <= base.n_residue(); ++ir) {
							for(Size ia = 1; ia <= base.residue(ir).nheavyatoms(); ia++) {
								if(ir==ihis || ir==jhis) {
									if(base.residue(ir).atom_name(ia)=="NE2" || base.residue(ir).atom_name(ia)=="ND1") continue;
								} else if(ir==ibpy) {
									if(base.residue(ir).atom_name(ia)=="NE1" || base.residue(ir).atom_name(ia)=="NN1") continue;
								}
								Vec x1 = R1*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
								Vec x2 = R2*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
								for(Size jr = 1; jr <= base.n_residue(); ++jr) {
									for(Size ja = 1; ja <= 5; ja++) {
										if( x1.distance_squared(base.xyz(AtomID(ja,jr))) < 9.0 ) goto clash4;
										if( x1.distance_squared( Rzn*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < 9.0 ) goto clash4;
										if( x2.distance_squared( Rzn*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < 9.0 ) goto clash4;

									}
								}
							}
						}
						goto noclash4;	clash4: continue; noclash4:


						Mat                 Rsymm = alignVectorSets( (c3f1-isct).normalized(), (c2f1-isct).normalized(), Vec(0,0,1), Vec( 0.816496579408716,0, 0.57735027133783) );
						if(orig_ang > 90.0) Rsymm = alignVectorSets( (c3f1-isct).normalized(), (c2f1-isct).normalized(), Vec(0,0,1), Vec(-0.816496585484756,0,-0.577350262745012) );

						Real baserep = sfrep->score(base);

						Pose symm = base;
						Pose symm_bare = symm;
						symm.append_residue_by_jump(*core::conformation::ResidueFactory::create_residue(rs->name_map("FE")),ibpy);
						symm.set_xyz(AtomID(1,symm.n_residue()),c3f1);
						symm.append_residue_by_jump(*core::conformation::ResidueFactory::create_residue(rs->name_map("ZN")),ihis);
						symm.set_xyz(AtomID(1,symm.n_residue()),c2f1);

						trans_pose(symm,-isct);
						rot_pose(symm,Rsymm);
						core::pose::symmetry::make_symmetric_pose(symm);
						trans_pose(symm_bare,-isct);
						rot_pose(symm_bare,Rsymm);
						core::pose::symmetry::make_symmetric_pose(symm_bare);

							// utility::io::ozstream out("test.pdb");
							// base.dump_pdb(out);
							// Vec viz(0,0,0);
							// viz = Rsymm*(c2f1        - isct); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"B"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = Rsymm*(c2f1+a2f1*5 - isct); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = Rsymm*(c3f1        - isct); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"C"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = Rsymm*(c3f1+a3f1*5 - isct); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"D"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// viz = Rsymm*(isct        - isct); out<<"HETATM"<<I(5,9999)<<' '<<"ZN  "<<' '<<" ZN"<<' '<<"I"<<I(4,3)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
							// out.close();
							// utility_exit_with_message("debug");

						// symm.conformation().declare_chemical_bond( ihis, "ND1", jhis, "ND1" );
						// symm.conformation().declare_chemical_bond( ihis, "NE2", jhis, "ND1" );
						// symm.conformation().declare_chemical_bond( ihis, "ND1", jhis, "NE2" );
						// symm.conformation().declare_chemical_bond( ihis, "NE2", jhis, "NE2" );


						string fname = infile+"_H"+lead_zero_string_of(ihis,3)+"-"+lead_zero_string_of(jhis,3)+"-"+lead_zero_string_of(ilig,4)+"_B"+lead_zero_string_of(ibpy,3)+"-"+lead_zero_string_of((Size)(bch1),3)+"_"+lead_zero_string_of(iaxs,1)+"_"+lead_zero_string_of(jaxs,1)+".pdb";


						sf->score(symm_bare);
						// TR << "REP " << 2*baserep << " " << symm_bare.energies().total_energies()[core::scoring::fa_rep] << std::endl;
						if( symm_bare.energies().total_energies()[core::scoring::fa_rep] - 2*baserep > 1000.0) continue;

						TR << "HIT " << ihis << "-" << jhis << " " << ibpy << " " << bch1 << " " << iaxs << " " << jaxs << std::endl;

						Real tmp1 = symm.chi(1,ibpy);
						Real tmp2 = symm.chi(2,ibpy);
						core::pose::replace_pose_residue_copying_existing_coordinates(symm,ibpy,rs->name_map("PHE"));
						symm.set_chi(1,ibpy,tmp1);
						symm.set_chi(2,ibpy,tmp2);

						symm.dump_pdb(fname);
						core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
						ss_out->fill_struct(symm_bare,fname);
						sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

						core::pose::replace_pose_residue_copying_existing_coordinates(base,ibpy,rs->name_map("BPY"));
						symm.set_chi(1,ibpy,tmp1);
						symm.set_chi(2,ibpy,tmp2);

					}
				}
			}
		}
	}



/*

							string fname = infile+"_"+lead_zero_string_of(irsd,3)+"_"+lead_zero_string_of((Size)(chi1+180.0),3)+"_"+lead_zero_string_of((Size)(chi2+180.0),3)+"_"+lead_zero_string_of(jrsd,3)+"_"+lead_zero_string_of((Size)(totrot),3)+"_"+lead_zero_string_of(iaxs,3)+".pdb";

							sf->score(symm);
							if(symm.energies().total_energies()[core::scoring::rg] > 14.0) {
								TR << "DANGER DANGER DANGER!!!" << std::endl;
								continue;
							}
							symm.dump_pdb(fname);

							core::pose::replace_pose_residue_copying_existing_coordinates(symm,jrsd,rs->name_map("ALA"));
							sf->score(symm);
							core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
							ss_out->fill_struct(symm,fname);
							sfd.write_silent_struct( *ss_out, option[ basic::options::OptionKeys::out::file::silent ]() );

						}

					}



				}
				// utility_exit_with_message("debug SG");
			} else {
				//matches?;
			}


		}
	}
*/
}

int main (int argc, char *argv[]) {

	try {


	// Vec p(-6.456746, 5.922204, -0.982538);
	// Vec d(0.393718,  0.677101,  0.621707);
	// Vec v(0.000000,  0.000000,  0.000000);
	// Vec a(0.998233, -0.003844, -0.059301);
	// Real t = 45;
	// vector1<Vec> X = line_cone_intersection(p,d,v,a,t);
	// for(Size i = 1; i <= X.size(); ++i) {
	// 	TR << "X " << i << " " << X[i] << std::endl;
	// }
	// utility_exit_with_message("debug line_cone_intersection");

	devel::init(argc,argv);
	run();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}




//
//








