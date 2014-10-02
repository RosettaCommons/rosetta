// includes
	#include <basic/options/keys/in.OptionKeys.gen.hh>
	#include <basic/options/keys/out.OptionKeys.gen.hh>
	#include <basic/options/keys/match.OptionKeys.gen.hh>
	#include <basic/options/keys/willmatch.OptionKeys.gen.hh>
	#include <basic/options/option.hh>
	#include <basic/options/util.hh>
	#include <basic/Tracer.hh>
	#include <core/chemical/ChemicalManager.hh>
	#include <core/chemical/ResidueTypeSet.hh>
	#include <core/chemical/util.hh>
	#include <core/conformation/Residue.hh>
	#include <core/conformation/ResidueFactory.hh>
	#include <core/import_pose/import_pose.hh>
	#include <devel/init.hh>
	#include <core/io/pdb/pose_io.hh>
	#include <core/kinematics/Stub.hh>
	#include <core/kinematics/MoveMap.hh>
	#include <core/graph/Graph.hh>
	#include <core/pack/packer_neighbors.hh>
	#include <core/pack/rotamer_set/RotamerSetFactory.hh>
	#include <core/pack/rotamer_set/RotamerSet.hh>
	#include <core/pack/optimizeH.hh>
	#include <core/pack/dunbrack/RotamerLibrary.hh>
	#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
	#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
	#include <core/pack/task/PackerTask.hh>
	#include <core/pack/task/TaskFactory.hh>
	#include <core/pose/annotated_sequence.hh>
	#include <core/pose/Pose.hh>
	#include <core/pose/PDBInfo.hh>
	#include <core/pose/symmetry/util.hh>
	#include <core/pose/util.hh>
	#include <core/scoring/Energies.hh>
	#include <core/scoring/rms_util.hh>
	#include <core/scoring/sasa.hh>
	#include <core/scoring/ScoreFunction.hh>
	#include <core/scoring/ScoreFunctionFactory.hh>
	#include <core/scoring/ScoringManager.hh>
	#include <core/scoring/constraints/AtomPairConstraint.hh>
	#include <core/scoring/constraints/DihedralConstraint.hh>
	#include <core/scoring/constraints/AngleConstraint.hh>
	#include <core/io/silent/ScoreFileSilentStruct.hh>
	#include <core/io/silent/SilentFileData.hh>
	#include <numeric/xyz.functions.hh>
	#include <numeric/xyz.io.hh>
	#include <ObjexxFCL/FArray2D.hh>
	#include <ObjexxFCL/format.hh>
	#include <ObjexxFCL/string.functions.hh>
	#include <protocols/scoring/ImplicitFastClashCheck.hh>
	#include <protocols/simple_moves/PackRotamersMover.hh>
	#include <protocols/simple_moves/MinMover.hh>
	#include <numeric/kinematic_closure/bridgeObjects.hh>
	#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
	#include <sstream>
	#include <utility/io/izstream.hh>
	#include <utility/io/ozstream.hh>
	#include <numeric/NumericTraits.hh>



	// #include <devel/init.hh>

	// #include <core/scoring/constraints/LocalCoordinateConstraint.hh>
	#include <apps/pilot/will/will_util.ihh>

	using core::Real;
	using core::Size;
	using core::pose::Pose;
	using core::kinematics::Stub;
	using protocols::scoring::ImplicitFastClashCheck;
	using std::string;
	using utility::vector1;
	using numeric::min;
	using core::import_pose::pose_from_pdb;
	using basic::options::option;
	using numeric::min;
	using numeric::max;
	using utility::io::ozstream;
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;

	typedef utility::vector1<core::Real> Reals;
	typedef utility::vector1<core::Size> Sizes;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;

static thread_local basic::Tracer TR( "test_ikrs" );

core::io::silent::SilentFileData sfd;
vector1<Reals> vecs2vv(Vecs const & v) {
	vector1<Reals> vv;
	for(Vecs::const_iterator i = v.begin(); i != v.end(); ++i) {
		Reals r(3);
		r[1] = i->x();
		r[2] = i->y();
		r[3] = i->z();
		vv.push_back(r);
	}
	return vv;
}
vector1<Vec> vv2vecs(vector1<Reals> const & vv) {
	vector1<Vec> r;
	for(vector1<Reals>::const_iterator i = vv.begin(); i != vv.end(); ++i) {
		r.push_back(Vec( (*i)[1], (*i)[2], (*i)[3] ));
	}
	return r;
}
struct Hit : utility::pointer::ReferenceCount {
	Hit(core::conformation::Residue const & r1,
			core::conformation::Residue const & r2,
			core::kinematics::Stub & stb,
			Real pd1 = 0.0,
			Real pd2 = 0.0,
			bool frnt_in=false, bool negori=false) : stub(stb), phid1(pd1), phid2(pd2), frnt(frnt_in), nego(negori) {
		rsd1 = r1.seqpos();
		rsd2 = r2.seqpos();
		k1 = r1.chi(1);
		k2 = r1.chi(2);
		k3 = r1.chi(3);
		k4 = r1.chi(4);
		c1 = r2.chi(1);
		c2 = r2.chi(2);
		if(r2.aa()==core::chemical::aa_glu) c3 = r2.chi(3);
		aa1 = r1.aa();
		aa2 = r2.aa();
	}

	void apply(Pose & wp) const {
		core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd1,wp.residue(rsd1).residue_type_set().name_map("GLY")); // reset all SC
		core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd2,wp.residue(rsd2).residue_type_set().name_map("GLY")); // reset all SC
		/**/ if(aa1==core::chemical::aa_arg) core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd1,wp.residue(rsd1).residue_type_set().name_map("ARG"));
		else if(aa1==core::chemical::aa_lys) core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd1,wp.residue(rsd1).residue_type_set().name_map("LYS"));
		else utility_exit_with_message("unknown aa "+str(aa1)+" for Hit class rsd1");
		/**/ if(aa2==core::chemical::aa_asp) core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd2,wp.residue(rsd2).residue_type_set().name_map("ASP"));
		else if(aa2==core::chemical::aa_glu) core::pose::replace_pose_residue_copying_existing_coordinates(wp,rsd2,wp.residue(rsd2).residue_type_set().name_map("GLU"));
		else utility_exit_with_message("unknown aa "+str(aa2)+" for Hit class rsd2");
		wp.set_chi(   1   ,       rsd1, k1 );
		wp.set_chi(   2   ,       rsd1, k2 );
		wp.set_chi(   3   ,       rsd1, k3 );
		wp.set_chi(   4   ,       rsd1, k4 );
		wp.set_chi(   1   ,       rsd2, c1 );
		wp.set_chi(   2   ,       rsd2, c2 );
		if(wp.residue(rsd2).nchi() > 2) wp.set_chi(3,rsd2, c3 );
		if(0.0 != phid1) {
			Vec C1 = wp.xyz(AtomID(2,rsd1));
			Mat M1 = rotation_matrix_degrees( C1-wp.xyz(AtomID(1,rsd1)) , phid1 );
			for(Size i = 5; i <= wp.residue_type(rsd1).natoms(); ++i) wp.set_xyz( AtomID(i,rsd1), M1 * (wp.xyz(AtomID(i,rsd1))-C1) + C1 );
		}
		if(0.0 != phid2) {
			Vec C2 = wp.xyz(AtomID(2,rsd2));
			Mat M2 = rotation_matrix_degrees( C2-wp.xyz(AtomID(1,rsd2)) , phid2 );
			for(Size i = 5; i <= wp.residue_type(rsd2).natoms(); ++i) wp.set_xyz( AtomID(i,rsd2), M2 * (wp.xyz(AtomID(i,rsd2))-C2) + C2 );
		}
	}

	Size rsd1,rsd2;
	core::kinematics::Stub stub;
	Real phid1,phid2,k1,k2,k3,k4,c1,c2,c3;
	core::chemical::AA aa1,aa2;
	bool frnt,nego;
};
typedef utility::pointer::owning_ptr<Hit> HitOP;

Real ik_arg_glu_frnt(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 144.0 ) return 0;
	Real mxcb = 0.0;
	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//  Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	Size iO2=ctp.residue(1).atom_index("O2"), iO5=ctp.residue(1).atom_index("O5"), iC5=ctp.residue(1).atom_index("C5");
	Size iO1=ctp.residue(1).atom_index("O1"), iC2=ctp.residue(1).atom_index("C2"), iC7=ctp.residue(1).atom_index("C7");
	Size iCZ=pose.residue(rsd1).atom_index("CZ"),iNE=pose.residue(rsd1).atom_index("NE"),iNH2=pose.residue(rsd1).atom_index("NH2");

	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	// N CA C N CA CB CG CD NE
	Vecs atoms(18);
	atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
	atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
	atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
	atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
	atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
	atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
	atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
	atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
	atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
	atoms[10] = pose.residue(rsd2  ).xyz( "CG" );
	atoms[11] = pose.residue(rsd2  ).xyz( "CB" );
	atoms[12] = pose.residue(rsd2  ).xyz( "CA" );
	atoms[13] = pose.residue(rsd2  ).xyz( "N"  );
	atoms[14] = pose.residue(rsd2-1).xyz( "C"  );
	atoms[15] = pose.residue(rsd2-1).xyz( "CA" );
	atoms[16] = pose.residue(rsd2-1).xyz( "N"  );
	atoms[17] = pose.residue(rsd2-2).xyz( "C"  );
	atoms[18] = pose.residue(rsd2-2).xyz( "CA" );

	order [1]=1; order [2]=2; order [3]=3;
	pivots[1]=5, pivots[2]=8, pivots[3]=11;

	using namespace numeric::kinematic_closure;
	chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	db_len[ 9] = 7.1;
	db_ang[ 9] = numeric::angle_degrees(pose.residue(rsd1).xyz("CD"),pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"));
	db_ang[10] = numeric::angle_degrees(pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"),pose.residue(rsd2).xyz("CD"));
	Real phitgt = dt_ang[4];
	Size count = 0;
	for(Size idh = 6; idh <= 360; idh += 12) {
		dt_ang[9] = (Real)idh;
		for(Size ichi2 = 6; ichi2 <= 360; ichi2 += 12) {
			dt_ang[6] = (Real)ichi2;
			bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
			for(int isol = 1; isol <= nsol; isol++) {
				Real phidiff = phitgt-t_ang[isol][4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
				if( fabs(phidiff) > 10.0 ) continue;

				utility::vector1<utility::vector1<core::Real> > vv_atm_out;
				numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
				Vecs apos = vv2vecs(vv_atm_out);

				bool clash = false;
				for( Size i = 1; i <= atoms.size(); ++i ) {
					for( Size j = i+3; j <= atoms.size(); ++j ) {
						if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
					}
					if(clash) break;
				}
				if(clash) continue;

				pose.set_chi(1,rsd1,t_ang[isol][ 5]);
				pose.set_chi(2,rsd1,t_ang[isol][ 6]);
				pose.set_chi(3,rsd1,t_ang[isol][ 7]);
				pose.set_chi(4,rsd1,t_ang[isol][ 8]);
				pose.set_chi(3,rsd2,t_ang[isol][ 9]); // this is tricky...
				pose.set_chi(2,rsd2,t_ang[isol][10]);
				pose.set_chi(1,rsd2,t_ang[isol][11]);
				//Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
				//if( dun1 > arg_dun_th ) continue;
				//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
				//if( dun2 > coo_dun_th ) continue;

				Vec CA = pose.xyz(AtomID(2,rsd1));
				Mat M = rotation_matrix_degrees( CA-pose.xyz(AtomID(1,rsd1)) , -phidiff );
				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M * (pose.xyz(AtomID(i,rsd1))-CA) + CA );

				for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
				for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }
				count++;
				if(clash/*||count!=20*/){if(clash)count--;for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M.transposed()*(pose.xyz(AtomID(i,rsd1))-CA)+CA);continue;}

				Vec const cen = pose.residue(rsd1).xyz(iCZ);
				Vec const axs = ((pose.residue(rsd1).xyz(iNH2)+pose.residue(rsd1).xyz(iNE))/2.0 - cen).normalized();
				Vec const ori = axs.cross(pose.residue(rsd1).xyz(iNE)-cen).normalized();
				Vec const or2 = ori.cross(axs);

				for(Size negori = 0; negori < 2; negori++) {
					Vec ctp_axs = (ctp.residue(1).xyz(iC2)-ctp.residue(1).xyz(iC7)).normalized();
					Vec ctp_ori = ctp_axs.cross( (ctp.residue(1).xyz(iO5) - ctp.residue(1).xyz(iC7)).normalized() );
					if(negori) ctp_ori = -1.0 * ctp_ori;
					Mat R = numeric::alignVectorSets( ctp_axs, ctp_ori, axs, ori );

					rot_pose(ctp,R);
					Vec t = cen + 4.05*axs - ctp.residue(1).xyz(iC7);
					trans_pose(ctp,t);

					// clash check CTP w/ scaffold BB
					//for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) { if( !ifc.clash_check(ctp.residue(1).xyz(i)) ) clash=true; if(clash) break; } if(clash) continue;

					core::kinematics::Stub lcstub(ctp.residue(1).xyz(iC5),ctp.residue(1).xyz(iO1),ctp.residue(1).xyz(iO2));

					Vec ccom(0,0,0); for(Size i=1; i<=ctp.residue(1).nheavyatoms(); ++i) ccom+=ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					hits.push_back( new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub, -phidiff, 0.0, true, negori ) );

					// ozstream out("ikrs_arg_frnt_glu_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(idh,3)+"_"+str(ichi2)+"_"+str(isol)+"_"+str(negori)+"_res.pdb");
					// Size ano = 0;
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
					// core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);
					// //out<<"ATOM  "<<I(5,11)<<' '<<" COM"<<' '<<"COM"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,ccom.x())<<F(8,3,ccom.y())<<F(8,3,ccom.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// //out<<"ATOM  "<<I(5,11)<<' '<<" CEN"<<' '<<"CEN"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,cen.x())<<F(8,3,cen.y())<<F(8,3,cen.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// //out<<"ATOM  "<<I(5,11)<<' '<<" AXS"<<' '<<"AXS"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,(cen+axs).x())<<F(8,3,(cen+axs).y())<<F(8,3,(cen+axs).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out.close();

					// for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M.transposed() * (pose.xyz(AtomID(i,rsd1))-CA) + CA );
					// Hit & rhit( *(hits[hits.size()]) );
					// Pose wp(pose);
					// rhit.apply(wp);
					// //wp.dump_pdb("test2.pdb");
					// Pose tmp;
					// tmp.append_residue_by_jump(wp.residue(rhit.rsd1),1);
					// tmp.append_residue_by_jump(wp.residue(rhit.rsd2),1);
					// tmp.dump_pdb("test.pdb");
					// utility_exit_with_message("test.pdb");

				}
				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M.transposed() * (pose.xyz(AtomID(i,rsd1))-CA) + CA );
			}
		}
	}

	return mxcb;
}
Real ik_arg_glu_side(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 100.0 ) return 0;
	Real mxcb = 0.0;
	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	Size iO2=ctp.residue(1).atom_index("O2"), iO5=ctp.residue(1).atom_index("O5"), iC5=ctp.residue(1).atom_index("C5");
	Size iO1=ctp.residue(1).atom_index("O1"), iC2=ctp.residue(1).atom_index("C2"), iC7=ctp.residue(1).atom_index("C7");
	Size iCZ=pose.residue(rsd1).atom_index("CZ"),iNH1=pose.residue(rsd1).atom_index("NH1"),iNH2=pose.residue(rsd1).atom_index("NH2");

	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	// N CA C N CA CB CG CD NE
	Vecs atoms(18);
	atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
	atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
	atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
	atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
	atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
	atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
	atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
	atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
	atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
	atoms[10] = pose.residue(rsd1  ).xyz( "CZ" );

	atoms[11] = pose.residue(rsd2  ).xyz( "CG" );
	atoms[12] = pose.residue(rsd2  ).xyz( "CB" );
	atoms[13] = pose.residue(rsd2  ).xyz( "CA" );
	atoms[14] = pose.residue(rsd2  ).xyz( "N"  );
	atoms[15] = pose.residue(rsd2-1).xyz( "C"  );
	atoms[16] = pose.residue(rsd2-1).xyz( "CA" );
	atoms[17] = pose.residue(rsd2-1).xyz( "N"  );
	atoms[18] = pose.residue(rsd2-2).xyz( "C"  );

	order [1]=1; order [2]=2; order [3]=3;
	pivots[1]=5, pivots[2]=8, pivots[3]=11;

	using namespace numeric::kinematic_closure;
	chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	dt_ang[ 9] = 180.0;
	db_len[10] = 5.8;
	db_ang[10] = 180.0 - numeric::angle_degrees(pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"),pose.residue(rsd1).xyz("NH1"));
	db_ang[11] =         numeric::angle_degrees(pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"),pose.residue(rsd2).xyz("CD"));
	Real phitgt = dt_ang[4];
	Size count = 0;
	for(Size idh = 6; idh <= 360; idh += 12) {
		dt_ang[12] = (Real)idh;
		for(Size ichi2 = 6; ichi2 <= 360; ichi2 += 12) {
			dt_ang[6] = (Real)ichi2;
			bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
			for(int isol = 1; isol <= nsol; isol++) {
				Real phidiff = phitgt-t_ang[isol][4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
				if( fabs(phidiff) > 10.0 ) continue;

				utility::vector1<utility::vector1<core::Real> > vv_atm_out;
				numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
				Vecs apos = vv2vecs(vv_atm_out);

				bool clash = false;
				for( Size i = 1; i <= atoms.size(); ++i ) {
					for( Size j = i+3; j <= atoms.size(); ++j ) {
						if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
					}
					if(clash) break;
				}
				if(clash) continue;

				pose.set_chi(1,rsd1,t_ang[isol][ 5]);
				pose.set_chi(2,rsd1,t_ang[isol][ 6]);
				pose.set_chi(3,rsd1,t_ang[isol][ 7]);
				pose.set_chi(4,rsd1,t_ang[isol][ 8]);
				pose.set_chi(3,rsd2,t_ang[isol][10]); // this is tricky...
				pose.set_chi(2,rsd2,t_ang[isol][11]);
				// don't set chi1 rsd2 -- should stay as orig because we've explicitly moved dt_ang[12]
				//Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
				//if( dun1 > arg_dun_th ) continue;
				//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
				//if( dun2 > coo_dun_th ) continue;

				Vec CA = pose.xyz(AtomID(2,rsd1));
				Mat M = rotation_matrix_degrees( CA-pose.xyz(AtomID(1,rsd1)) , -phidiff );
				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M * (pose.xyz(AtomID(i,rsd1))-CA) + CA );

				for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
				for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

				count++;
				if(clash){//||count!=10){
					if(clash)count--;
					for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M.transposed()*(pose.xyz(AtomID(i,rsd1))-CA)+CA);
					continue;
				}

				Vec const cen = pose.residue(rsd1).xyz(iCZ);
				Vec const axs = ((pose.residue(rsd1).xyz(iNH2)+pose.residue(rsd1).xyz(iNH1))/2.0 - cen).normalized();
				Vec const ori = axs.cross(pose.residue(rsd1).xyz(iNH1)-cen).normalized();
				Vec const or2 = ori.cross(axs);

				for(Size negori = 0; negori < 2; negori++) {
					Vec ctp_axs = (ctp.residue(1).xyz(iC2)-ctp.residue(1).xyz(iC7)).normalized();
					Vec ctp_ori = ctp_axs.cross( (ctp.residue(1).xyz(iO5) - ctp.residue(1).xyz(iC7)).normalized() );
					if(negori) ctp_ori = -1.0 * ctp_ori;
					Mat R = numeric::alignVectorSets( ctp_axs, ctp_ori, axs, ori );

					rot_pose(ctp,R);
					Vec t = cen + 4.05*axs - ctp.residue(1).xyz(iC7);
					trans_pose(ctp,t);

					// clash check CTP w/ scaffold BB
					//for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) { if( !ifc.clash_check(ctp.residue(1).xyz(i)) ) clash=true; if(clash) break; } if(clash) continue;

					core::kinematics::Stub lcstub(ctp.residue(1).xyz(iC5),ctp.residue(1).xyz(iO1),ctp.residue(1).xyz(iO2));

					Vec ccom(0,0,0); for(Size i=1; i<=ctp.residue(1).nheavyatoms(); ++i) ccom+=ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					hits.push_back( new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub, -phidiff, 0.0, false, negori ) );

					// ozstream out("ikrs_arg_side_glu_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(idh,3)+"_"+str(ichi2)+"_"+str(isol)+"_"+str(negori)+"_res.pdb");
					// Size ano = 0;
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
					// core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);
					// out<<"ATOM  "<<I(5,11)<<' '<<" COM"<<' '<<"COM"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,ccom.x())<<F(8,3,ccom.y())<<F(8,3,ccom.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out<<"ATOM  "<<I(5,11)<<' '<<" CEN"<<' '<<"CEN"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,cen.x())<<F(8,3,cen.y())<<F(8,3,cen.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out<<"ATOM  "<<I(5,11)<<' '<<" AXS"<<' '<<"AXS"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,(cen+axs).x())<<F(8,3,(cen+axs).y())<<F(8,3,(cen+axs).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out.close();
					// utility_exit_with_message("asldfkj");

				}

				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M.transposed() * (pose.xyz(AtomID(i,rsd1))-CA) + CA );
			}
		}
	}
	return mxcb;
}
Real ik_arg_asp_frnt(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 120.0 ) return 0;
	Real mxcb = 0.0;
	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	Size iO2=ctp.residue(1).atom_index("O2"), iO5=ctp.residue(1).atom_index("O5"), iC5=ctp.residue(1).atom_index("C5");
	Size iO1=ctp.residue(1).atom_index("O1"), iC2=ctp.residue(1).atom_index("C2"), iC7=ctp.residue(1).atom_index("C7");
	Size iCZ=pose.residue(rsd1).atom_index("CZ"),iNE=pose.residue(rsd1).atom_index("NE"),iNH2=pose.residue(rsd1).atom_index("NH2");

	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	// N CA C N CA CB CG CD NE
	Vecs atoms(15);
	atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
	atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
	atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
	atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
	atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
	atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
	atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
	atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
	atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
	atoms[10] = pose.residue(rsd2  ).xyz( "CB" );
	atoms[11] = pose.residue(rsd2  ).xyz( "CA" );
	atoms[12] = pose.residue(rsd2  ).xyz( "N"  );
	atoms[13] = pose.residue(rsd2-1).xyz( "C"  );
	atoms[14] = pose.residue(rsd2-1).xyz( "CA" );
	atoms[15] = pose.residue(rsd2-1).xyz( "N"  );

	order [1]=1; order [2]=2; order [3]=3;
	pivots[1]=5, pivots[2]=8, pivots[3]=11;

	using namespace numeric::kinematic_closure;
	chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);

	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	db_len[ 9] = 7.1;
	db_ang[ 9] = numeric::angle_degrees(pose.residue(rsd1).xyz("CD"),pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"));
	db_ang[10] = numeric::angle_degrees(pose.residue(rsd2).xyz("CA"),pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"));

	Real phitgt = dt_ang[ 4];
	Real psitgt = dt_ang[11];
	Size count = 0;
	for(Size idh = 6; idh <= 360; idh += 12) {
		dt_ang[ 9] = (Real)idh;
		for(Size ichi2 = 6; ichi2 <= 360; ichi2 += 12) {
			dt_ang[6] = (Real)ichi2;
			bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
			//      if(nsol==0) break;
			for(int isol = 1; isol <= nsol; isol++) {
				Real phidiff = phitgt-t_ang[isol][ 4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
				Real ph2diff = psitgt-t_ang[isol][11]; while(ph2diff < -180.0) ph2diff+=360.0; while(ph2diff > 180.0) ph2diff-=360.0;
				if( fabs(phidiff) > 13.0 ) continue;
				if( fabs(ph2diff) > 13.0 ) continue;

				utility::vector1<utility::vector1<core::Real> > vv_atm_out;
				numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
				Vecs apos = vv2vecs(vv_atm_out);

				bool clash = false;
				for( Size i = 1; i <= atoms.size(); ++i ) {
					for( Size j = i+3; j <= atoms.size(); ++j ) {
						if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
					}
					if(clash) break;
				}
				if(clash) continue;

				//for(Size i = 1; i <= atoms.size(); ++i) TR << "ICHI2 " << ichi2 << " SOL " << isol << " CHI " << i << " " << t_ang[isol][i] << std::endl;
				pose.set_chi(1,rsd1,t_ang[isol][ 5]);
				pose.set_chi(2,rsd1,t_ang[isol][ 6]);
				pose.set_chi(3,rsd1,t_ang[isol][ 7]);
				pose.set_chi(4,rsd1,t_ang[isol][ 8]);
				pose.set_chi(2,rsd2,t_ang[isol][ 9]); // this is tricky...
				pose.set_chi(1,rsd2,t_ang[isol][10]);
				//Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
				//if( dun1 > arg_dun_th ) continue;
				//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
				//if( dun2 > coo_dun_th ) continue;

				Vec CA1 = pose.xyz(AtomID(2,rsd1));
				Mat M1 = rotation_matrix_degrees( CA1-pose.xyz(AtomID(1,rsd1)) , -phidiff );
				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1 * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
				Vec CA2 = pose.xyz(AtomID(2,rsd2));
				Mat M2 = rotation_matrix_degrees( CA2-pose.xyz(AtomID(1,rsd2)) , -ph2diff );
				for(Size i = 5; i <= pose.residue_type(rsd2).natoms(); ++i) pose.set_xyz( AtomID(i,rsd2), M2 * (pose.xyz(AtomID(i,rsd2))-CA2) + CA2 );

				for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
				for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

				count++;
				if(clash){//||count!=1){
					if(clash) count--;
					for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M1.transposed()*(pose.xyz(AtomID(i,rsd1))-CA1)+CA1);
					for(Size i=5; i<=pose.residue_type(rsd2).natoms();++i) pose.set_xyz(AtomID(i,rsd2),M2.transposed()*(pose.xyz(AtomID(i,rsd2))-CA2)+CA2);
					continue;
				}
				Vec const cen = pose.residue(rsd1).xyz(iCZ);
				Vec const axs = ((pose.residue(rsd1).xyz(iNH2)+pose.residue(rsd1).xyz(iNE))/2.0 - cen).normalized();
				Vec const ori = axs.cross(pose.residue(rsd1).xyz(iNE)-cen).normalized();
				Vec const or2 = ori.cross(axs);

				for(Size negori = 0; negori < 2; negori++) {
					Vec ctp_axs = (ctp.residue(1).xyz(iC2)-ctp.residue(1).xyz(iC7)).normalized();
					Vec ctp_ori = ctp_axs.cross( (ctp.residue(1).xyz(iO5) - ctp.residue(1).xyz(iC7)).normalized() );
					if(negori) ctp_ori = -1.0 * ctp_ori;
					Mat R = numeric::alignVectorSets( ctp_axs, ctp_ori, axs, ori );

					rot_pose(ctp,R);
					Vec t = cen + 4.05*axs - ctp.residue(1).xyz(iC7);
					trans_pose(ctp,t);

					// clash check CTP w/ scaffold BB
					//          for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) { if( !ifc.clash_check(ctp.residue(1).xyz(i)) ) clash=true; if(clash) break; } if(clash) continue;

					Vec ccom(0,0,0); for(Size i=1; i<=ctp.residue(1).nheavyatoms(); ++i) ccom+=ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();
					core::kinematics::Stub lcstub(ctp.residue(1).xyz(iC5),ctp.residue(1).xyz(iO1),ctp.residue(1).xyz(iO2));

					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					hits.push_back( new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub, -phidiff, -ph2diff, true, negori ) );

					// ozstream out("ikrs_arg_frnt_asp_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(idh,3)+"_"+str(ichi2)+"_"+str(isol)+"_"+str(negori)+"_res.pdb");
					// Size ano = 0;
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
					// core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
					// core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);
					// out<<"ATOM  "<<I(5,11)<<' '<<" COM"<<' '<<"COM"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,ccom.x())<<F(8,3,ccom.y())<<F(8,3,ccom.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out<<"ATOM  "<<I(5,11)<<' '<<" CEN"<<' '<<"CEN"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,cen.x())<<F(8,3,cen.y())<<F(8,3,cen.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out<<"ATOM  "<<I(5,11)<<' '<<" AXS"<<' '<<"AXS"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,(cen+axs).x())<<F(8,3,(cen+axs).y())<<F(8,3,(cen+axs).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					// out.close();

				}
				for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1.transposed() * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
				for(Size i = 5; i <= pose.residue_type(rsd2).natoms(); ++i) pose.set_xyz( AtomID(i,rsd2), M2.transposed() * (pose.xyz(AtomID(i,rsd2))-CA2) + CA2 );
			}
		}
		//    if(nsol==0) break;
	}
	return mxcb;
}
Real ik_arg_asp_side(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 100.0 ) return 0;
	Real mxcb = 0.0;
	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//  Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	Size iO2=ctp.residue(1).atom_index("O2"), iO5=ctp.residue(1).atom_index("O5"), iC5=ctp.residue(1).atom_index("C5");
	Size iO1=ctp.residue(1).atom_index("O1"), iC2=ctp.residue(1).atom_index("C2"), iC7=ctp.residue(1).atom_index("C7");
	Size iCZ=pose.residue(rsd1).atom_index("CZ"),iNH1=pose.residue(rsd1).atom_index("NH1"),iNH2=pose.residue(rsd1).atom_index("NH2");

	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	Vecs atoms(15);
	atoms[ 1] = pose.residue(rsd1-1).xyz( "N"  );
	atoms[ 2] = pose.residue(rsd1-1).xyz( "CA" );
	atoms[ 3] = pose.residue(rsd1-1).xyz( "C"  );
	atoms[ 4] = pose.residue(rsd1  ).xyz( "N"  );
	atoms[ 5] = pose.residue(rsd1  ).xyz( "CA" );
	atoms[ 6] = pose.residue(rsd1  ).xyz( "CB" );
	atoms[ 7] = pose.residue(rsd1  ).xyz( "CG" );
	atoms[ 8] = pose.residue(rsd1  ).xyz( "CD" );
	atoms[ 9] = pose.residue(rsd1  ).xyz( "NE" );
	atoms[10] = pose.residue(rsd1  ).xyz( "CZ" );
	atoms[11] = pose.residue(rsd2  ).xyz( "CB" );
	atoms[12] = pose.residue(rsd2  ).xyz( "CA" );
	atoms[13] = pose.residue(rsd2  ).xyz( "N"  );
	atoms[14] = pose.residue(rsd2-1).xyz( "C"  );
	atoms[15] = pose.residue(rsd2-1).xyz( "CA" );

	order [1]=1; order [2]=2; order [3]=3;
	pivots[1]=5, pivots[2]=8, pivots[3]=11;
	numeric::kinematic_closure::chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);
	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd1).type() );
	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	dt_ang[ 9] = 180.0;
	db_len[10] = 5.8;
	db_ang[10] = 180.0 - numeric::angle_degrees(pose.residue(rsd1).xyz("NE"),pose.residue(rsd1).xyz("CZ"),pose.residue(rsd1).xyz("NH1"));
	db_ang[11] =         numeric::angle_degrees(pose.residue(rsd2).xyz("CA"),pose.residue(rsd2).xyz("CB"),pose.residue(rsd2).xyz("CG"));

	Real phitgt = dt_ang[ 4];
	Size count = 0;
	for(Size ichi2 = 6; ichi2 <= 360; ichi2 += 12) {
		dt_ang[6] = (Real)ichi2;
		numeric::kinematic_closure::bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
		//    if(nsol==0) break;
		for(int isol = 1; isol <= nsol; isol++) {
			Real phidiff = phitgt-t_ang[isol][ 4]; while(phidiff < -180.0) phidiff+=360.0; while(phidiff > 180.0) phidiff-=360.0;
			if( fabs(phidiff) > 13.0 ) continue;

			utility::vector1<utility::vector1<core::Real> > vv_atm_out;
			numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
			Vecs apos = vv2vecs(vv_atm_out);

			bool clash = false;
			for( Size i = 1; i <= atoms.size(); ++i ) {
				for( Size j = i+3; j <= atoms.size(); ++j ) {
					if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
				}
				if(clash) break;
			}
			if(clash) continue;

			pose.set_chi(1,rsd1,t_ang[isol][ 5]);
			pose.set_chi(2,rsd1,t_ang[isol][ 6]);
			pose.set_chi(3,rsd1,t_ang[isol][ 7]);
			pose.set_chi(4,rsd1,t_ang[isol][ 8]);
			pose.set_chi(2,rsd2,t_ang[isol][10]); // this is tricky...
			pose.set_chi(1,rsd2,t_ang[isol][11]);
			//Real dun1 = dunlib1->rotamer_energy( pose.residue(rsd1), scratch );
			//if( dun1 > arg_dun_th ) continue;
			//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
			//if( dun2 > coo_dun_th ) continue;

			Vec CA1 = pose.xyz(AtomID(2,rsd1));
			Mat M1 = rotation_matrix_degrees( CA1-pose.xyz(AtomID(1,rsd1)) , -phidiff );
			for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1 * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );

			for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
			for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }

			count++;
			if(clash){//||count!=1){
				if(clash) count--;
				for(Size i=5; i<=pose.residue_type(rsd1).natoms();++i) pose.set_xyz(AtomID(i,rsd1),M1.transposed()*(pose.xyz(AtomID(i,rsd1))-CA1)+CA1);
				continue;
			}

			Vec const cen = pose.residue(rsd1).xyz(iCZ);
			Vec const axs = ((pose.residue(rsd1).xyz(iNH2)+pose.residue(rsd1).xyz(iNH1))/2.0 - cen).normalized();
			Vec const ori = axs.cross(pose.residue(rsd1).xyz(iNH1)-cen).normalized();
			Vec const or2 = ori.cross(axs);

			for(Size negori = 0; negori < 2; negori++) {
				Vec ctp_axs = (ctp.residue(1).xyz(iC2)-ctp.residue(1).xyz(iC7)).normalized();
				Vec ctp_ori = ctp_axs.cross( (ctp.residue(1).xyz(iO5) - ctp.residue(1).xyz(iC7)).normalized() );
				if(negori) ctp_ori = -1.0 * ctp_ori;
				Mat R = numeric::alignVectorSets( ctp_axs, ctp_ori, axs, ori );

				rot_pose(ctp,R);
				Vec t = cen + 4.05*axs - ctp.residue(1).xyz(iC7);
				trans_pose(ctp,t);

				// clash check CTP w/ scaffold BB
				//        for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) { if( !ifc.clash_check(ctp.residue(1).xyz(i)) ) clash=true; if(clash) break; } if(clash) continue;

				Vec ccom(0,0,0); for(Size i=1; i<=ctp.residue(1).nheavyatoms(); ++i) ccom+=ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();
				core::kinematics::Stub lcstub(ctp.residue(1).xyz(iC5),ctp.residue(1).xyz(iO1),ctp.residue(1).xyz(iO2));

				#ifdef USE_OPENMP
				#pragma omp critical
				#endif
				hits.push_back( new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub, -phidiff, 0.0, false, negori ) );

				// ozstream out("ikrs_arg_side_asp_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+str(ichi2)+"_"+str(isol)+"_"+str(negori)+"_res.pdb");
				// Size ano = 0;
				// core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
				// core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
				// core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
				// core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
				// core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);
				// out<<"ATOM  "<<I(5,11)<<' '<<" COM"<<' '<<"COM"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,ccom.x())<<F(8,3,ccom.y())<<F(8,3,ccom.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// out<<"ATOM  "<<I(5,11)<<' '<<" CEN"<<' '<<"CEN"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,cen.x())<<F(8,3,cen.y())<<F(8,3,cen.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// out<<"ATOM  "<<I(5,11)<<' '<<" AXS"<<' '<<"AXS"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,(cen+axs).x())<<F(8,3,(cen+axs).y())<<F(8,3,(cen+axs).z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// out.close();

			}
			for(Size i = 5; i <= pose.residue_type(rsd1).natoms(); ++i) pose.set_xyz( AtomID(i,rsd1), M1.transposed() * (pose.xyz(AtomID(i,rsd1))-CA1) + CA1 );
		}
		//    if(nsol==0) break;
	}
	return mxcb;
}
void ik_lys_ctp_asp(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 225.0 ) return;

	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//  Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	core::kinematics::Stub cstub(ctp.residue(1).xyz("C5"),ctp.residue(1).xyz("O1"),ctp.residue(1).xyz("O2"));
	for(Size i = 1; i <= ctp.residue(1).natoms(); ++i) ctp.set_xyz(AtomID(i,1), cstub.global2local(ctp.xyz(AtomID(i,1))));
	Vec ccom(0,0,0); for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) ccom += ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();

	Size iOD1=pose.residue(rsd2).atom_index("OD1"),iCG=pose.residue(rsd2).atom_index("CG"),iCB=pose.residue(rsd2).atom_index("CB");

	utility::vector1<Size> pivots(3), order(3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0(3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0(3);

	order [1]= 1; order [2]= 2; order [3]= 3;
	pivots[1]= 8, pivots[2]=11, pivots[3]=14;

	//core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	//core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	// get rotamers
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn( pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
	dummy_task->initialize_from_command_line();
	dummy_task->nonconst_residue_task( rsd1 ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( rsd1 ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( rsd1 ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
	core::pack::rotamer_set::RotamerSetFactory rsf;
	core::pack::rotamer_set::RotamerSetOP rotset( rsf.create_rotamer_set( pose.residue( rsd1 ) ) );
	rotset->set_resid( rsd1 );
	rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );

	//TR << "ik_lys_ctp_asp" << rsd1 << " " << rsd2 << " " << rotset->num_rotamers() << std::endl;

	// do IK on each lys rot
	for(Size krot = 1; krot <= rotset->num_rotamers(); ++krot) {
		//TR << "ik_lys_ctp_asp" << rsd1 << " " << rsd2 << " " << krot << std::endl;
		pose.set_chi(1,rsd1,rotset->rotamer(krot)->chi(1));
		pose.set_chi(2,rsd1,rotset->rotamer(krot)->chi(2));
		pose.set_chi(3,rsd1,rotset->rotamer(krot)->chi(3));
		pose.set_chi(4,rsd1,rotset->rotamer(krot)->chi(4));
		Vec const ctp_oc  = pose.residue(rsd1).xyz("NZ" ) + 3.0*(pose.residue(rsd1).xyz("1HZ")-pose.residue(rsd1).xyz("NZ")).normalized();
		Vec const ctp_oh  = pose.residue(rsd2).xyz("OD1") + 3.0*(pose.residue(rsd2).xyz("CG" )-pose.residue(rsd2).xyz("CB")).normalized();
		Vec const ctp_coh = ctp_oh + 1.5*(pose.residue(rsd2).xyz("OD1" )-pose.residue(rsd2).xyz("CG")).normalized();
		Vecs atoms(21);
		atoms[ 1] = pose.residue(rsd1-1).xyz( "C"  );
		atoms[ 2] = pose.residue(rsd1  ).xyz( "N"  );
		atoms[ 3] = pose.residue(rsd1  ).xyz( "CA" );
		atoms[ 4] = pose.residue(rsd1  ).xyz( "CB" );
		atoms[ 5] = pose.residue(rsd1  ).xyz( "CG" );
		atoms[ 6] = pose.residue(rsd1  ).xyz( "CD" );
		atoms[ 7] = pose.residue(rsd1  ).xyz( "CE" );
		atoms[ 8] = pose.residue(rsd1  ).xyz( "NZ" );
		atoms[ 9] = ctp_oc;
		atoms[10] = ctp_coh;
		atoms[11] = ctp_oh;
		atoms[12] = pose.residue(rsd2  ).xyz( "OD1");
		atoms[13] = pose.residue(rsd2  ).xyz( "CG" );
		atoms[14] = pose.residue(rsd2  ).xyz( "CB" );
		atoms[15] = pose.residue(rsd2  ).xyz( "CA" );
		atoms[16] = pose.residue(rsd2  ).xyz( "N"  );
		atoms[17] = pose.residue(rsd2-1).xyz( "C"  );
		atoms[18] = pose.residue(rsd2-1).xyz( "CA" );
		atoms[19] = pose.residue(rsd2-1).xyz( "N"  );
		atoms[20] = pose.residue(rsd2-2).xyz( "C"  );
		atoms[21] = pose.residue(rsd2-2).xyz( "CA" );
		numeric::kinematic_closure::chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);
		db_len[ 9] = 2.7836;
		db_ang[10] = 146.3;
		db_ang[11] = 105.0;
		db_ang[12] = 120.0;

		for(Size ock = 0; ock <= 1; ock++) {
			db_ang[ 9] = ock?98.0: 72.9;
			dt_ang[ 9] = ock?56.0:-62.9;
			for(Size ocd = 0; ocd <= 1; ocd++) {
				dt_ang[12] = ocd?0.0:180.0;

				Size count = 0;
				numeric::kinematic_closure::bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
				for(int isol = 1; isol <= nsol; isol++) {

					utility::vector1<utility::vector1<core::Real> > vv_atm_out;
					numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
					Vecs apos = vv2vecs(vv_atm_out);

					bool clash = false;
					for( Size i = 1; i <= atoms.size(); ++i ) {
						for( Size j = i+3; j <= atoms.size(); ++j ) {
							if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
						}
						if(clash) break;
					}
					if(clash) continue;

					pose.set_chi(2,rsd2,t_ang[isol][13]);
					pose.set_chi(1,rsd2,t_ang[isol][14]);
					//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
					//if( dun2 > coo_dun_th*2.0/3.0 ) continue;

					for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
					for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }
					if(clash) continue;

					core::kinematics::Stub ikstub(apos[12],apos[13],apos[14]);
					core::kinematics::Stub postub(pose.residue(rsd2).xyz(iOD1),pose.residue(rsd2).xyz(iCG),pose.residue(rsd2).xyz(iCB));
					Vec const tmp_ctp_oc  = postub.local2global(ikstub.global2local(apos[ 9]));
					Vec const tmp_ctp_coh = postub.local2global(ikstub.global2local(apos[10]));
					Vec const tmp_ctp_oh  = postub.local2global(ikstub.global2local(apos[11]));

					core::kinematics::Stub lcstub(tmp_ctp_coh,tmp_ctp_oh,tmp_ctp_oc);

					for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) {
						Vec X = lcstub.local2global(ctp.residue(1).xyz(i));
						for(Size j = 1; j <= pose.residue(rsd1).nheavyatoms()-1; ++j) if( X.distance_squared(pose.residue(rsd1).xyz(j)) < 8.0 ) clash=true;
						for(Size j = 1; j <= pose.residue(rsd2).nheavyatoms()-3; ++j) if( X.distance_squared(pose.residue(rsd2).xyz(j)) < 8.0 ) clash=true;
						//if( !ifc.clash_check( X ) ) clash=true;
						if(clash) break;
					}
					if(clash) continue;

					Vec comtmp = lcstub.local2global(ccom);
					HitOP h = new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub );
					#ifdef USE_OPENMP
					#pragma omp critical
					#endif
					hits.push_back(h);
					count++;

					// if( ock!=1 && ocd!=1 ) continue;

					// if(1){
					//   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(krot,3)+"_"+str(ock)+"_"+str(ocd)+"_res.pdb");
					//   Size ano = 0;
					//   core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
					//   core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
					//   core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
					//   core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
					//  for(Size i = 1; i <= ctp.residue(1).natoms(); ++i) ctp.set_xyz(AtomID(i,1), lcstub.local2global(ctp.xyz(AtomID(i,1))));
					//   core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);

					//  out<<"ATOM  "<<I(5, 9)<<' '<<" O2 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,tmp_ctp_oc.x())<<F(8,3,tmp_ctp_oc.y())<<F(8,3,tmp_ctp_oc.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,10)<<' '<<" C5 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,tmp_ctp_coh.x())<<F(8,3,tmp_ctp_coh.y())<<F(8,3,tmp_ctp_coh.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,11)<<' '<<" O1 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,tmp_ctp_oh.x())<<F(8,3,tmp_ctp_oh.y())<<F(8,3,tmp_ctp_oh.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

					//   out.close();
					// }

					// {
					//   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(krot,3)+"_"+str(ock)+"_"+str(ocd)+".pdb");

					//   out<<"ATOM  "<<I(5, 2)<<' '<<" N  "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 2].x())<<F(8,3,apos[ 2].y())<<F(8,3,apos[ 2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 3)<<' '<<" CA "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 3].x())<<F(8,3,apos[ 3].y())<<F(8,3,apos[ 3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 4)<<' '<<" CB "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 4].x())<<F(8,3,apos[ 4].y())<<F(8,3,apos[ 4].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 5)<<' '<<" CG "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 5].x())<<F(8,3,apos[ 5].y())<<F(8,3,apos[ 5].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 6)<<' '<<" CD "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 6].x())<<F(8,3,apos[ 6].y())<<F(8,3,apos[ 6].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 7)<<' '<<" CE "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 7].x())<<F(8,3,apos[ 7].y())<<F(8,3,apos[ 7].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 8)<<' '<<" NZ "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 8].x())<<F(8,3,apos[ 8].y())<<F(8,3,apos[ 8].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5, 9)<<' '<<" O2 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 9].x())<<F(8,3,apos[ 9].y())<<F(8,3,apos[ 9].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,10)<<' '<<" C5 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[10].x())<<F(8,3,apos[10].y())<<F(8,3,apos[10].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,11)<<' '<<" O1 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[11].x())<<F(8,3,apos[11].y())<<F(8,3,apos[11].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

					//   out<<"ATOM  "<<I(5,19)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[19].x())<<F(8,3,apos[19].y())<<F(8,3,apos[19].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,18)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[18].x())<<F(8,3,apos[18].y())<<F(8,3,apos[18].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,17)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[17].x())<<F(8,3,apos[17].y())<<F(8,3,apos[17].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

					//   out<<"ATOM  "<<I(5,16)<<' '<<" N  "<<' '<<"ASP"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[16].x())<<F(8,3,apos[16].y())<<F(8,3,apos[16].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,15)<<' '<<" CA "<<' '<<"ASP"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[15].x())<<F(8,3,apos[15].y())<<F(8,3,apos[15].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,14)<<' '<<" CB "<<' '<<"ASP"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[14].x())<<F(8,3,apos[14].y())<<F(8,3,apos[14].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,13)<<' '<<" CG "<<' '<<"ASP"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[13].x())<<F(8,3,apos[13].y())<<F(8,3,apos[13].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
					//   out<<"ATOM  "<<I(5,12)<<' '<<" OD1"<<' '<<"ASP"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[12].x())<<F(8,3,apos[12].y())<<F(8,3,apos[12].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

					//   out.close();
					// }
				}
			}
		}
	}
}
void ik_lys_ctp_glu(Pose & pose, Size rsd1, Size rsd2, ImplicitFastClashCheck & ifc, vector1<HitOP> & hits, Pose ctp) {
	using namespace basic::options::OptionKeys;
	using namespace core::id;

	if( pose.residue(rsd1).xyz("CB").distance_squared(pose.residue(rsd2).xyz("CB")) > 300.0 ) return;

	//Real arg_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::arg_dun_th]();
	//  Real coo_dun_th = basic::options::option[basic::options::OptionKeys::willmatch::glu_dun_th]();

	core::kinematics::Stub cstub(ctp.residue(1).xyz("C5"),ctp.residue(1).xyz("O1"),ctp.residue(1).xyz("O2"));
	for(Size i = 1; i <= ctp.residue(1).natoms(); ++i) ctp.set_xyz(AtomID(i,1), cstub.global2local(ctp.xyz(AtomID(i,1))));
	Vec ccom(0,0,0); for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) ccom += ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();

	Size iOE1=pose.residue(rsd2).atom_index("OE1"),iCD=pose.residue(rsd2).atom_index("CD"),iCG=pose.residue(rsd2).atom_index("CG");

	utility::vector1<Size> pivots (3), order (3);
	utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
	int nsol=0;
	// for eliminating identical solutions
	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> dt_ang, db_len, db_ang, R0 (3);

	order [1]=1; order [2]= 2; order [3]= 3;
	pivots[1]=8, pivots[2]=11, pivots[3]=14;

	//  core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib2 = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( pose.residue(rsd2).type() );
	///core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	// get rotamers
	core::pack::rotamer_set::RotamerSetOP rotset;
	{
		core::scoring::ScoreFunction dummy_sfxn;
		dummy_sfxn( pose );
		core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( pose );
		dummy_task->initialize_from_command_line();
		dummy_task->nonconst_residue_task( rsd1 ).restrict_to_repacking();
		dummy_task->nonconst_residue_task( rsd1 ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
		dummy_task->nonconst_residue_task( rsd1 ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
		core::graph::GraphOP dummy_png = core::pack::create_packer_graph( pose, dummy_sfxn, dummy_task );
		core::pack::rotamer_set::RotamerSetFactory rsf;
		rotset = rsf.create_rotamer_set( pose.residue( rsd1 ) );
		rotset->set_resid( rsd1 );
		rotset->build_rotamers( pose, dummy_sfxn, *dummy_task, dummy_png );
	}
	vector1<Real> echi1; echi1.push_back(-175); echi1.push_back(-75), echi1.push_back(-65), echi1.push_back(-55), echi1.push_back(175);
	//TR << "ik_lys_ctp_glu" << rsd1y << " " << rsd2 << " " << rotset->num_rotamers() << std::endl;

	// do IK on each lys rot
	for(Size krot = 1; krot <= rotset->num_rotamers(); ++krot) {
		pose.set_chi(1,rsd1,rotset->rotamer(krot)->chi(1));
		pose.set_chi(2,rsd1,rotset->rotamer(krot)->chi(2));
		pose.set_chi(3,rsd1,rotset->rotamer(krot)->chi(3));
		pose.set_chi(4,rsd1,rotset->rotamer(krot)->chi(4));
		Vec const ctp_oc  = pose.residue(rsd1).xyz("NZ" ) + 3.0*(pose.residue(rsd1).xyz("1HZ")-pose.residue(rsd1).xyz("NZ")).normalized();
		for(Size iech = 1; iech <= echi1.size(); ++iech) {
			Real ech1 = echi1[iech];
			pose.set_chi(1,rsd2,ech1);
			Vec const ctp_oh  = pose.residue(rsd2).xyz(iOE1) + 3.0*(pose.residue(rsd2).xyz("CD" )-pose.residue(rsd2).xyz("CG")).normalized();
			Vec const ctp_coh = ctp_oh + 1.5*(pose.residue(rsd2).xyz(iOE1)-pose.residue(rsd2).xyz("CD")).normalized();
			Vecs atoms(21);
			atoms[ 1] = pose.residue(rsd1-1).xyz( "C"  );
			atoms[ 2] = pose.residue(rsd1  ).xyz( "N"  );
			atoms[ 3] = pose.residue(rsd1  ).xyz( "CA" );
			atoms[ 4] = pose.residue(rsd1  ).xyz( "CB" );
			atoms[ 5] = pose.residue(rsd1  ).xyz( "CG" );
			atoms[ 6] = pose.residue(rsd1  ).xyz( "CD" );
			atoms[ 7] = pose.residue(rsd1  ).xyz( "CE" );
			atoms[ 8] = pose.residue(rsd1  ).xyz( "NZ" );
			atoms[ 9] = ctp_oc;
			atoms[10] = ctp_coh;
			atoms[11] = ctp_oh;
			atoms[12] = pose.residue(rsd2  ).xyz( "OE1");
			atoms[13] = pose.residue(rsd2  ).xyz( "CD" );
			atoms[14] = pose.residue(rsd2  ).xyz( "CG" );
			atoms[15] = pose.residue(rsd2  ).xyz( "CB" );
			atoms[16] = pose.residue(rsd2  ).xyz( "CA" );
			atoms[17] = pose.residue(rsd2  ).xyz( "N"  );
			atoms[18] = pose.residue(rsd2-1).xyz( "C"  );
			atoms[19] = pose.residue(rsd2-1).xyz( "CA" );
			atoms[20] = pose.residue(rsd2-1).xyz( "N"  );
			atoms[21] = pose.residue(rsd2-2).xyz( "C"  );
			numeric::kinematic_closure::chainTORS(atoms.size(), vecs2vv(atoms), dt_ang, db_ang, db_len, R0, Q0);
			db_len[ 9] = 2.7836;
			db_ang[10] = 146.3;
			db_ang[11] = 105.0;
			db_ang[12] = 120.0;

			for(Size ock = 0; ock <= 1; ock++) {
				db_ang[ 9] = ock?98.0: 72.9;
				dt_ang[ 9] = ock?56.0:-62.9;
				for(Size ocd = 0; ocd <= 1; ocd++) {
					dt_ang[12] = ocd?0.0:180.0;

					Size count = 0;
					numeric::kinematic_closure::bridgeObjects(vecs2vv(atoms), dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
					for(int isol = 1; isol <= nsol; isol++) {

						utility::vector1<utility::vector1<core::Real> > vv_atm_out;
						numeric::kinematic_closure::chainXYZ(atoms.size(),b_len[isol],b_ang[isol],t_ang[isol],false,R0,Q0,vv_atm_out);
						Vecs apos = vv2vecs(vv_atm_out);

						bool clash = false;
						for( Size i = 1; i <= atoms.size(); ++i ) {
							for( Size j = i+3; j <= atoms.size(); ++j ) {
								if( apos[i].distance_squared(apos[j]) < 8.0 ) { clash=true; break; }
							}
							if(clash) break;
						}
						if(clash) continue;

						pose.set_chi(3,rsd2,t_ang[isol][13]);
						pose.set_chi(2,rsd2,t_ang[isol][14]);
						pose.set_chi(1,rsd2,t_ang[isol][15]);
						//Real dun2 = dunlib2->rotamer_energy( pose.residue(rsd2), scratch );
						//if( dun2 > coo_dun_th ) continue;

						for(Size i = 6; i <= pose.residue(rsd1).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd1)), rsd1 ) ) { clash=true; break; }
						for(Size i = 6; i <= pose.residue(rsd2).nheavyatoms(); ++i) if(! ifc.clash_check( pose.xyz(AtomID(i,rsd2)), rsd2 ) ) { clash=true; break; }
						if(clash) continue;

						core::kinematics::Stub ikstub(apos[12],apos[13],apos[14]);
						core::kinematics::Stub postub(pose.residue(rsd2).xyz(iOE1),pose.residue(rsd2).xyz(iCD),pose.residue(rsd2).xyz(iCG));
						Vec const tmp_ctp_oc  = postub.local2global(ikstub.global2local(apos[ 9]));
						Vec const tmp_ctp_coh = postub.local2global(ikstub.global2local(apos[10]));
						Vec const tmp_ctp_oh  = postub.local2global(ikstub.global2local(apos[11]));

						core::kinematics::Stub lcstub(tmp_ctp_coh,tmp_ctp_oh,tmp_ctp_oc);

						for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) {
							Vec X = lcstub.local2global(ctp.residue(1).xyz(i));
							for(Size j = 1; j <= pose.residue(rsd1).nheavyatoms()-1; ++j) if( X.distance_squared(pose.residue(rsd1).xyz(j)) < 8.0 ) clash=true;
							for(Size j = 1; j <= pose.residue(rsd2).nheavyatoms()-3; ++j) if( X.distance_squared(pose.residue(rsd2).xyz(j)) < 8.0 ) clash=true;
							//if( !ifc.clash_check( X ) ) clash=true;
							if(clash) break;
						}
						if(clash) continue;

						Vec comtmp = lcstub.local2global(ccom);
						HitOP h = new Hit( pose.residue(rsd1), pose.residue(rsd2), lcstub );
						#ifdef USE_OPENMP
						#pragma omp critical
						#endif
						hits.push_back(h);
						count++;

						// if(ock!=1 || ocd!=0) continue;

						// if(1){
						//   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(krot,3)+"_"+str(ock)+"_"+str(ocd)+"_res.pdb");
						//   Size ano = 0;
						//   core::io::pdb::dump_pdb_residue(pose.residue(rsd2-1),ano,out);
						//   core::io::pdb::dump_pdb_residue(pose.residue(rsd2  ),ano,out);
						//   core::io::pdb::dump_pdb_residue(pose.residue(rsd1-1),ano,out);
						//   core::io::pdb::dump_pdb_residue(pose.residue(rsd1  ),ano,out);
						//  for(Size i = 1; i <= ctp.residue(1).natoms(); ++i) ctp.set_xyz(AtomID(i,1), lcstub.local2global(ctp.xyz(AtomID(i,1))));
						//  core::io::pdb::dump_pdb_residue(ctp.residue(1),ano,out);
						//  out<<"ATOM  "<<I(5, 9)<<' '<<" O2 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3, tmp_ctp_oc.x())<<F(8,3, tmp_ctp_oc.y())<<F(8,3, tmp_ctp_oc.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//  out<<"ATOM  "<<I(5,10)<<' '<<" C5 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,tmp_ctp_coh.x())<<F(8,3,tmp_ctp_coh.y())<<F(8,3,tmp_ctp_coh.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//  out<<"ATOM  "<<I(5,11)<<' '<<" O1 "<<' '<<"VIZ"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3, tmp_ctp_oh.x())<<F(8,3, tmp_ctp_oh.y())<<F(8,3, tmp_ctp_oh.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//  out<<"ATOM  "<<I(5,11)<<' '<<" COM"<<' '<<"COM"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,     h->com.x())<<F(8,3,     h->com.y())<<F(8,3,     h->com.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out.close();
						// }

						// {
						//   ozstream out("ikrs_"+lzs(rsd1,3)+"_"+lzs(rsd2,3)+"_"+lzs(krot,3)+"_"+str(ock)+"_"+str(ocd)+".pdb");

						//   out<<"ATOM  "<<I(5, 2)<<' '<<" N  "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 2].x())<<F(8,3,apos[ 2].y())<<F(8,3,apos[ 2].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 3)<<' '<<" CA "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 3].x())<<F(8,3,apos[ 3].y())<<F(8,3,apos[ 3].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 4)<<' '<<" CB "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 4].x())<<F(8,3,apos[ 4].y())<<F(8,3,apos[ 4].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 5)<<' '<<" CG "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 5].x())<<F(8,3,apos[ 5].y())<<F(8,3,apos[ 5].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 6)<<' '<<" CD "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 6].x())<<F(8,3,apos[ 6].y())<<F(8,3,apos[ 6].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 7)<<' '<<" CE "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 7].x())<<F(8,3,apos[ 7].y())<<F(8,3,apos[ 7].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 8)<<' '<<" NZ "<<' '<<"ARG"<<' '<<"A"<<I(4,1)<<"    "<<F(8,3,apos[ 8].x())<<F(8,3,apos[ 8].y())<<F(8,3,apos[ 8].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5, 9)<<' '<<" O2 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[ 9].x())<<F(8,3,apos[ 9].y())<<F(8,3,apos[ 9].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,10)<<' '<<" C5 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[10].x())<<F(8,3,apos[10].y())<<F(8,3,apos[10].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,11)<<' '<<" O1 "<<' '<<"CTP"<<' '<<"A"<<I(4,2)<<"    "<<F(8,3,apos[11].x())<<F(8,3,apos[11].y())<<F(8,3,apos[11].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

						//   out<<"ATOM  "<<I(5,20)<<' '<<" N  "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[20].x())<<F(8,3,apos[20].y())<<F(8,3,apos[20].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,19)<<' '<<" CA "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[19].x())<<F(8,3,apos[19].y())<<F(8,3,apos[19].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,18)<<' '<<" C  "<<' '<<"ALA"<<' '<<"A"<<I(4,3)<<"    "<<F(8,3,apos[18].x())<<F(8,3,apos[18].y())<<F(8,3,apos[18].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

						//   out<<"ATOM  "<<I(5,17)<<' '<<" N  "<<' '<<"GLU"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[17].x())<<F(8,3,apos[17].y())<<F(8,3,apos[17].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,16)<<' '<<" CA "<<' '<<"GLU"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[16].x())<<F(8,3,apos[16].y())<<F(8,3,apos[16].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,15)<<' '<<" CB "<<' '<<"GLU"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[15].x())<<F(8,3,apos[15].y())<<F(8,3,apos[15].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,14)<<' '<<" CG "<<' '<<"GLU"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[14].x())<<F(8,3,apos[14].y())<<F(8,3,apos[14].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,13)<<' '<<" CD "<<' '<<"GLU"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[13].x())<<F(8,3,apos[13].y())<<F(8,3,apos[13].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
						//   out<<"ATOM  "<<I(5,12)<<' '<<" OE1"<<' '<<"Glu"<<' '<<"A"<<I(4,4)<<"    "<<F(8,3,apos[12].x())<<F(8,3,apos[12].y())<<F(8,3,apos[12].z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';

						//   out.close();
						// }
					}
				}
			}
		}
	}
}

void repack(Pose & arg) {
	ScoreFunctionOP sf = core::scoring::get_score_function();
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task(arg);
	task->restrict_to_repacking();
	protocols::simple_moves::PackRotamersMover repack( sf, task );
	repack.apply(arg);
}

int main (int argc, char *argv[]) {

	try {

  devel::init(argc,argv);

	const core::Real PI = numeric::NumericTraits<Real>::pi();
	Pose ctp; pose_from_pdb(ctp,"input/ctp.pdb");
	Pose pose,arg,asp,glu,lys,ala;
	string infile = basic::options::option[basic::options::OptionKeys::in::file::s]()[1];
	core::import_pose::pose_from_pdb(pose,infile);
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
		if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
	}
	Pose pala(pose);
	for(Size ir = 1; ir <= pala.n_residue(); ++ir) {
		pala.replace_residue(ir,ala.residue(1),true);
	}

	ImplicitFastClashCheck ifc (pose,2.0);
	ImplicitFastClashCheck ifc3(pose,2.6);
	vector1<Real> sasa; { core::id::AtomID_Map<Real> atom_sasa; core::scoring::calc_per_atom_sasa( pose, atom_sasa, sasa, 5.0, false ); }
	if( basic::options::option[basic::options::OptionKeys::match::scaffold_active_site_residues].user() ) {
		string posfile = basic::options::option[basic::options::OptionKeys::match::scaffold_active_site_residues]();
		TR << "reading pos file " << posfile << std::endl;
		utility::io::izstream in(posfile);
		vector1<Size> pos;
		Size tmp; while(in >> tmp) pos.push_back(tmp);
		in.close();
		for(Size i = 1; i <= sasa.size(); ++i) {
			if( std::find(pos.begin(),pos.end(),i) == pos.end() ) sasa[i] = 9e9;
			else sasa[i] = 0.0;
		}
	}

	core::chemical::ResidueTypeSetCAP frs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::pose::make_pose_from_sequence(ala,"A",*frs,false); repack(ala);
	core::pose::make_pose_from_sequence(arg,"R",*frs,false); repack(arg);
	core::pose::make_pose_from_sequence(asp,"D",*frs,false); repack(asp);
	core::pose::make_pose_from_sequence(glu,"E",*frs,false); repack(glu);
	core::pose::make_pose_from_sequence(lys,"K",*frs,false); repack(lys);

	TR << "scanning LYS-LG-OOC" << std::endl;
	vector1<HitOP> khits; {
		vector1<vector1<HitOP> > hitsvec(pose.n_residue());
		Size nhit=0;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(int ir = 3; ir <= (int)pose.n_residue()-2; ++ir) {
			if(sasa[ir] > 0.1) continue;
			Pose wp,lg;
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{ wp=pose; lg=ctp; }
			wp.replace_residue(ir,lys.residue(1),true);
			for(Size jr = 3; jr <= pose.n_residue()-2; ++jr) {
				if(ir==jr) continue;
				if(sasa[jr] > 0.1) continue;
				//if( (ir*jr+jr+3*ir+999999999)%10!=0 ) continue;
				wp.replace_residue(jr,asp.residue(1),true);
				ik_lys_ctp_asp(wp,ir,jr,ifc,hitsvec[ir],lg);
				wp.replace_residue(jr,glu.residue(1),true);
				ik_lys_ctp_glu(wp,ir,jr,ifc,hitsvec[ir],lg);
			}
			nhit += hitsvec[ir].size();
			if(ir%3==0) TR << Real(ir-2)/Real(pose.n_residue()-4) * 100.0 << " percent done ik_lys_ctp " << nhit << " hits" << std::endl;
		}
		for(vector1<vector1<HitOP> >::const_iterator i = hitsvec.begin(); i != hitsvec.end(); ++i) khits.insert(khits.end(),i->begin(),i->end());
		hitsvec.clear(); // recover mem
	}

	TR << "scanning ARG-OOC" << std::endl;
	vector1<HitOP> rhits; {
		vector1<vector1<HitOP> > hitsvec(pose.n_residue());
		Size nhit=0;
		#ifdef USE_OPENMP
		#pragma omp parallel for schedule(dynamic,1)
		#endif
		for(int ir = 3; ir <= (int)pose.n_residue()-2; ++ir) {
			if(sasa[ir] > 0.1) continue;
			//TR << ir << std::endl;
			Pose wp,lg;
			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{ wp=pose; lg=ctp; }
			wp.replace_residue(ir,arg.residue(1),true);
			for(Size jr = 3; jr <= pose.n_residue()-2; ++jr) {
				if(ir==jr) continue;
				if(sasa[jr] > 0.1) continue;
				//if( (ir*jr+jr+3*ir+999999999)%5!=0 ) continue;
				wp.replace_residue(jr,asp.residue(1),true);
				ik_arg_asp_frnt(wp,ir,jr,ifc,hitsvec[ir],lg);
				ik_arg_asp_side(wp,ir,jr,ifc,hitsvec[ir],lg);
				wp.replace_residue(jr,glu.residue(1),true);
				ik_arg_glu_frnt(wp,ir,jr,ifc,hitsvec[ir],lg);
				ik_arg_glu_side(wp,ir,jr,ifc,hitsvec[ir],lg);
			}
			nhit += hitsvec[ir].size();
			if(ir%3==0) TR << Real(ir-2)/Real(pose.n_residue()-4) * 100.0 << " percent done ik_arg " << nhit << " hits" << std::endl;
		}
		for(vector1<vector1<HitOP> >::const_iterator i = hitsvec.begin(); i != hitsvec.end(); ++i) rhits.insert(rhits.end(),i->begin(),i->end());
		hitsvec.clear(); // recover mem
	}
	TR << "LYS/OOC HITS: " << khits.size() << std::endl;
	TR << "ARG bkp HITS: " << rhits.size() << std::endl;

	// set up lig position
	core::kinematics::Stub const lgstub(ctp.residue(1).xyz("C5"),ctp.residue(1).xyz("O1"),ctp.residue(1).xyz("O2"));
	for(Size i = 1; i <= ctp.residue(1).natoms(); ++i) ctp.set_xyz(AtomID(i,1), lgstub.global2local(ctp.xyz(AtomID(i,1))));
	Vec ccom(0,0,0); for(Size i = 1; i <= ctp.residue(1).nheavyatoms(); ++i) ccom += ctp.residue(1).xyz(i); ccom /= ctp.residue(1).nheavyatoms();
	Mat LM = numeric::alignVectorSets( ctp.residue(1).xyz("O4")-ctp.residue(1).xyz("C9"), ctp.residue(1).xyz("O3")-ctp.residue(1).xyz("C9"),
																		 ctp.residue(1).xyz("O5")-ctp.residue(1).xyz("C7"), ctp.residue(1).xyz("O6")-ctp.residue(1).xyz("C7") );
	Pose lgB(ctp);
	rot_pose(lgB,LM);
	trans_pose(lgB,ctp.residue(1).xyz("C7")-lgB.residue(1).xyz("C9"));
	Vec ccomB(0,0,0); for(Size i = 1; i <= lgB.residue(1).nheavyatoms(); ++i) ccomB += lgB.residue(1).xyz(i); ccomB /= lgB.residue(1).nheavyatoms();

	// set up lys grid
	Vec lb(9e9,9e9,9e9),ub(-9e9,-9e9,-9e9);
	for(Size i = 1; i <= khits.size(); ++i) {
		Vec com  = khits[i]->stub.local2global(ccom );  lb.min( com -3.0 ); ub.max( com +3.0 );
	}
	if(rhits.size()==0) { lb = 0; ub = 1; }
	ObjexxFCL::FArray3D<vector1<HitOP> > kgrid(std::ceil(ub.x()-lb.x()),std::ceil(ub.y()-lb.y()),std::ceil(ub.z()-lb.z()));
	for(Size i = 1; i <= khits.size(); ++i) {
		Vec com  = khits[i]->stub.local2global(ccom );
		Size ix = std::ceil( com.x() - lb.x() );
		Size iy = std::ceil( com.y() - lb.y() );
		Size iz = std::ceil( com.z() - lb.z() );
		kgrid(ix,iy,iz).push_back(khits[i]);
	}

	Size nhit=0,nrhit=0;

	Size const iO1=ctp.residue(1).atom_index("O1"), iC7=ctp.residue(1).atom_index("C7"), iC9=ctp.residue(1).atom_index("C9");
	Size const iO2=ctp.residue(1).atom_index("O2"), iC5=ctp.residue(1).atom_index("C5");
	Vec  const lgO1 (ctp.residue(1).xyz(iO1)),lgC7 (ctp.residue(1).xyz(iC7)),lgC9 (ctp.residue(1).xyz(iC9));
	Vec  const lgBO1(lgB.residue(1).xyz(iO1)),lgBC7(lgB.residue(1).xyz(iC7)),lgBC9(lgB.residue(1).xyz(iC9));

	#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
	#endif
	for(int ir = 1; ir <= (int)rhits.size(); ++ir) {


		if(ir%100==0) TR << Real(ir)/Real(rhits.size()) * 100.0 << " percent done combine " << nhit << " hits" << std::endl;
		Hit const rhit(*rhits[ir]);
		Vec rcom = rhit.stub.local2global(ccom);
		Pose wp,lg;
		ScoreFunctionOP sf;
		core::kinematics::MoveMapOP movemap;
		#ifdef USE_OPENMP
		#pragma omp critical
		#endif
		{
			sf = new core::scoring::ScoreFunction();
			movemap = new core::kinematics::MoveMap;
			wp=pala; lg=ctp;
		}
		rhit.apply(wp);

		sf->set_weight(core::scoring::fa_atr,1.0);
		sf->set_weight(core::scoring::fa_rep,1.0);
		sf->set_weight(core::scoring::fa_intra_rep,1.0);
		sf->set_weight(core::scoring::hbond_sc,5.0);
		sf->set_weight(core::scoring::atom_pair_constraint,1.1);
		sf->set_weight(core::scoring::    angle_constraint,1.1);
		movemap->set_bb(false);
		movemap->set_chi(true);
		movemap->set_jump(false);
		movemap->set_jump(4,true);
		protocols::simple_moves::MinMover minm( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
		for(Size is = 1; is <= rhits.size(); ++is) {
			Hit const shit(*rhits[is]);
			if(rhit.rsd1==shit.rsd1||rhit.rsd1==shit.rsd2||rhit.rsd2==shit.rsd1||rhit.rsd2==shit.rsd2) continue;
			Vec scom = shit.stub.local2global(ccomB);
			if(rcom.distance_squared(scom) > 9.0) continue;

			if(rhit.stub.local2global(lgO1).distance_squared(shit.stub.local2global(lgBO1)) > 4.0) continue;
			if(rhit.stub.local2global(lgC7).distance_squared(shit.stub.local2global(lgBC7)) > 4.0) continue;
			if(rhit.stub.local2global(lgC9).distance_squared(shit.stub.local2global(lgBC9)) > 4.0) continue;

			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			{ shit.apply(wp); }


			Pose lg1=lg; for(Size i = 1; i <= lg1.residue_type(1).natoms(); ++i) lg1.set_xyz(AtomID(i,1), rhit.stub.local2global(lg1.xyz(AtomID(i,1))));
			Pose tmp;
			tmp.append_residue_by_jump(wp.residue(rhit.rsd1),1);
			tmp.append_residue_by_jump(wp.residue(rhit.rsd2),1);
			tmp.append_residue_by_jump(wp.residue(shit.rsd1),1);
			tmp.append_residue_by_jump(wp.residue(shit.rsd2),1);
			tmp.append_residue_by_jump(lg1.residue(1),1);

			Size n11 = rhit.frnt?8:10, n12 = 11, n1a = rhit.frnt?10:8;
			Size n21 = shit.frnt?8:10, n22 = 11, n2a = shit.frnt?10:8;
			Size c11 = rhit.nego ? ctp.residue(1).atom_index("O6"): ctp.residue(1).atom_index("O5");
			Size c12 = rhit.nego ? ctp.residue(1).atom_index("O5"): ctp.residue(1).atom_index("O6");
			Size c21 = shit.nego ? ctp.residue(1).atom_index("O3"): ctp.residue(1).atom_index("O4");
			Size c22 = shit.nego ? ctp.residue(1).atom_index("O4"): ctp.residue(1).atom_index("O3");
			using namespace core::scoring::constraints;
			tmp.add_constraint( new AtomPairConstraint( AtomID(n11,1), AtomID(c11,5), new HarmonicFunc(2.8,0.2)) );
			tmp.add_constraint( new AtomPairConstraint( AtomID(n12,1), AtomID(c12,5), new HarmonicFunc(2.8,0.2)) );
			tmp.add_constraint( new AtomPairConstraint( AtomID(n21,3), AtomID(c21,5), new HarmonicFunc(2.8,0.2)) );
			tmp.add_constraint( new AtomPairConstraint( AtomID(n22,3), AtomID(c22,5), new HarmonicFunc(2.8,0.2)) );
			// TR << "n11 " << tmp.residue(1).atom_name(n11) << " c11 " << tmp.residue(5).atom_name(c11) << std::endl;
			// TR << "n12 " << tmp.residue(1).atom_name(n12) << " c12 " << tmp.residue(5).atom_name(c12) << std::endl;
			// TR << "n21 " << tmp.residue(3).atom_name(n21) << " c21 " << tmp.residue(5).atom_name(c21) << std::endl;
			// TR << "n22 " << tmp.residue(3).atom_name(n22) << " c22 " << tmp.residue(5).atom_name(c22) << std::endl;
			tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C2"),5),AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),               new HarmonicFunc(PI,0.3)) );
			tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),AtomID(n1a,1), new HarmonicFunc(PI,0.3)) );
			tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C8"),5),AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),               new HarmonicFunc(PI,0.3)) );
			tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),AtomID(n2a,3), new HarmonicFunc(PI,0.3)) );
			// TR << "C2 C7 " << tmp.residue(1).atom_name(9) << " " << tmp.residue(1).atom_name(n1a) << std::endl;
			// TR << "C8 C9 " << tmp.residue(3).atom_name(9) << " " << tmp.residue(3).atom_name(n2a) << std::endl;
			Size bk11 = (tmp.residue(2).aa()==core::chemical::aa_glu)?6:5, bk12 = (tmp.residue(2).aa()==core::chemical::aa_glu)?7:6, bk13 = 9, bk14=n11;
			tmp.add_constraint( new AngleConstraint( AtomID(bk11,2),AtomID(bk12,2),AtomID(bk13,1),                new HarmonicFunc(PI,0.3)) );
			tmp.add_constraint( new AngleConstraint(                AtomID(bk12,2),AtomID(bk13,1),AtomID(bk14,1), new HarmonicFunc(PI,0.3)) );
			Size bk21 = (tmp.residue(4).aa()==core::chemical::aa_glu)?6:5, bk22 = (tmp.residue(4).aa()==core::chemical::aa_glu)?7:6, bk23 = 9, bk24=n21;
			tmp.add_constraint( new AngleConstraint( AtomID(bk21,4),AtomID(bk22,4),AtomID(bk23,3),                new HarmonicFunc(PI,0.3)) );
			tmp.add_constraint( new AngleConstraint(                AtomID(bk22,4),AtomID(bk23,3),AtomID(bk24,3), new HarmonicFunc(PI,0.3)) );

			// #ifdef USE_OPENMP
			// #pragma omp critical
			// #endif
			minm.apply(tmp);

			{
				bool clash = false;
				for(Size i = 1; i <= tmp.residue(5).nheavyatoms(); ++i) {
					if( !ifc.clash_check( tmp.residue(5).xyz(i) ) ) clash = true;
				}
				if(clash) continue;
			}

			#ifdef USE_OPENMP
			#pragma omp critical
			#endif
			nrhit++;

			Vec const lg12O1(tmp.residue(5).xyz(iO1));
			Vec const lg12C7(tmp.residue(5).xyz(iC7));
			Vec const lg12C9(tmp.residue(5).xyz(iC9));
			Vec ccom12(0,0,0); for(Size i=1; i <= tmp.residue(5).nheavyatoms(); ++i) ccom12 += tmp.residue(5).xyz(i); ccom12 /= tmp.residue(5).nheavyatoms();
			Size ixr = std::ceil( ccom12.x() - lb.x() );
			Size iyr = std::ceil( ccom12.y() - lb.y() );
			Size izr = std::ceil( ccom12.z() - lb.z() );

			for(Size ix = max(ixr-1,((Size)1)); ix <= min(ixr+((Size)1),kgrid.size1()); ++ix) {
				for(Size iy = max(iyr-1,((Size)1)); iy <= min(iyr+((Size)1),kgrid.size2()); ++iy) {
					for(Size iz = max(izr-1,((Size)1)); iz <= min(izr+((Size)1),kgrid.size3()); ++iz) {
						vector1<HitOP> const & ktmphit( kgrid(ix,iy,iz) );
						for(Size ikh = 1; ikh <= ktmphit.size(); ++ikh) {
							Hit const khit(*ktmphit[ikh]);
							if(khit.rsd1==rhit.rsd1||khit.rsd1==rhit.rsd2||khit.rsd2==rhit.rsd1||khit.rsd2==rhit.rsd2) continue;
							if(khit.rsd1==shit.rsd1||khit.rsd1==shit.rsd2||khit.rsd2==shit.rsd1||khit.rsd2==shit.rsd2) continue;
							Vec kcom = khit.stub.local2global(ccom);
							//if( kcom.distance_squared(ccom12) > 0.25 ) continue;

							if(lg12O1.distance_squared(khit.stub.local2global(lgO1)) > 4.0) continue;
							if(lg12C7.distance_squared(khit.stub.local2global(lgC7)) > 4.0) continue;
							if(lg12C9.distance_squared(khit.stub.local2global(lgC9)) > 4.0) continue;



							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							khit.apply(wp);

							Pose lgK(lg); for(Size i = 1; i <= lgK.residue_type(1).natoms(); ++i) lgK.set_xyz(AtomID(i,1), khit.stub.local2global(lgK.xyz(AtomID(i,1))));

							Pose tmp2;
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							{ tmp2 = pala; }

							tmp2.append_residue_by_jump(wp.residue(rhit.rsd1),1);
							tmp2.append_residue_by_jump(wp.residue(rhit.rsd2),1);
							tmp2.append_residue_by_jump(wp.residue(shit.rsd1),1);
							tmp2.append_residue_by_jump(wp.residue(shit.rsd2),1);
							tmp2.append_residue_by_jump(lg1.residue(1),1);
							tmp2.append_residue_by_jump(wp.residue(khit.rsd1),1);
							tmp2.append_residue_by_jump(wp.residue(khit.rsd2),1);
							tmp2.add_constraint( new AtomPairConstraint( AtomID(9                              ,6), AtomID(iO2,5), new HarmonicFunc(3.0,0.2)) );
							tmp2.add_constraint( new AtomPairConstraint( AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), new HarmonicFunc(3.0,0.2)) );
							tmp2.add_constraint( new AngleConstraint( AtomID(tmp2.residue(7).nheavyatoms()-2,7),AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), new HarmonicFunc(PI*2.0/3.0,0.2)) );
							tmp2.add_constraint( new AngleConstraint( AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), AtomID(iC5,5), new HarmonicFunc(1.838,0.2)) );
							//tmp2.add_constraint( new DihedralConstraint( AtomID(tmp2.residue(7).nheavyatoms()-2,7),AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), AtomID(iC5,5), new HarmonicFunc(PI,0.1)) );

							//tmp2.append_residue_by_jump(lgK.residue(1),1);
							Pose prev(tmp2);

							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							minm.apply(tmp2);
							{
								bool clash = false;
								for(Size i = 1; i <= tmp.residue(5).nheavyatoms(); ++i) {
									if( !ifc3.clash_check( tmp.residue(5).xyz(i) ) ) clash = true;
								}
								if(clash) continue;
							}
							string fn = lzs(rhit.rsd1,3)+"_"+tmp2.residue(2).name1()+lzs(rhit.rsd2,3)
								+"_"+lzs(shit.rsd1,3)+"_"+tmp2.residue(4).name1()+lzs(shit.rsd2,3)
								+"_"+lzs(khit.rsd1,3)+"_"+tmp2.residue(7).name1()+lzs(khit.rsd2,3);
							//prev.dump_pdb("test_"+fn+".pdb");
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							tmp2.dump_pdb(utility::file_basename(infile)+"_"+fn+"_min.pdb");
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							nhit++;


							Real rd1 = tmp2.xyz(AtomID(n11,1)).distance(tmp2.xyz(AtomID(c11,5)));
							Real rd2 = tmp2.xyz(AtomID(n12,1)).distance(tmp2.xyz(AtomID(c12,5)));
							Real rd3 = tmp2.xyz(AtomID(n21,3)).distance(tmp2.xyz(AtomID(c21,5)));
							Real rd4 = tmp2.xyz(AtomID(n22,3)).distance(tmp2.xyz(AtomID(c22,5)));
							Real ra1 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C2"),5)),tmp2.xyz(AtomID(tmp.residue(5).atom_index("C7"),5)),tmp2.xyz(AtomID(9,1)));
							Real ra2 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C7"),5)),tmp2.xyz(AtomID(9,1)),tmp2.xyz(AtomID(n1a,1)));
							Real ra3 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C8"),5)),tmp2.xyz(AtomID(tmp.residue(5).atom_index("C9"),5)),tmp2.xyz(AtomID(9,3)));
							Real ra4 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C9"),5)),tmp2.xyz(AtomID(9,3)),tmp2.xyz(AtomID(n2a,3)));
							Real ra5 = numeric::angle_degrees(tmp2.xyz(AtomID(bk11,2)),tmp2.xyz(AtomID(bk12,2)),tmp2.xyz(AtomID(bk13,1)));
							Real ra6 = numeric::angle_degrees(tmp2.xyz(AtomID(bk12,2)),tmp2.xyz(AtomID(bk13,1)),tmp2.xyz(AtomID(bk14,1)));
							Real ra7 = numeric::angle_degrees(tmp2.xyz(AtomID(bk21,4)),tmp2.xyz(AtomID(bk22,4)),tmp2.xyz(AtomID(bk23,3)));
							Real ra8 = numeric::angle_degrees(tmp2.xyz(AtomID(bk22,4)),tmp2.xyz(AtomID(bk23,3)),tmp2.xyz(AtomID(bk24,3)));
							Real kd1 = tmp2.xyz(AtomID(9                              ,6)).distance(tmp2.xyz(AtomID(iO2,5)));
							Real kd2 = tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)).distance(tmp2.xyz(AtomID(iO1,5)));
							Real ka1 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-2,7)),tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)), tmp2.xyz(AtomID(iO1,5)));
							Real ka2 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)),tmp2.xyz(AtomID(iO1,5)), tmp2.xyz(AtomID(iC5,5)));





							// minimal pose:
							//               Pose tmp2(tmp);
							//               tmp2.append_residue_by_jump(wp.residue(khit.rsd1),1);
							//               tmp2.append_residue_by_jump(wp.residue(khit.rsd2),1);
							//               tmp2.add_constraint( new AtomPairConstraint( AtomID(9                              ,6), AtomID(iO2,5), new HarmonicFunc(3.0,0.2)) );
							//               tmp2.add_constraint( new AtomPairConstraint( AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), new HarmonicFunc(3.0,0.2)) );
							//               tmp2.add_constraint( new AngleConstraint( AtomID(tmp2.residue(7).nheavyatoms()-2,7),AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), new HarmonicFunc(PI*2.0/3.0,0.2)) );
							//               tmp2.add_constraint( new AngleConstraint( AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), AtomID(iC5,5), new HarmonicFunc(1.838,0.2)) );
							//               //tmp2.add_constraint( new DihedralConstraint( AtomID(tmp2.residue(7).nheavyatoms()-2,7),AtomID(tmp2.residue(7).nheavyatoms()-1,7), AtomID(iO1,5), AtomID(iC5,5), new HarmonicFunc(PI,0.1)) );

							//               //tmp2.append_residue_by_jump(lgK.residue(1),1);
							//               Pose prev(tmp2);

							// #ifdef USE_OPENMP
							// #pragma omp critical
							// #endif
							//               minm.apply(tmp2);
							//               {
							//                 bool clash = false;
							//                 for(Size i = 1; i <= tmp.residue(5).nheavyatoms(); ++i) {
							//                   if( !ifc3.clash_check( tmp.residue(5).xyz(i) ) ) clash = true;
							//                 }
							//                 if(clash) continue;
							//               }
							//               string fn = lzs(rhit.rsd1,3)+"_"+tmp2.residue(2).name1()+lzs(rhit.rsd2,3)
							//                 +"_"+lzs(shit.rsd1,3)+"_"+tmp2.residue(4).name1()+lzs(shit.rsd2,3)
							//                 +"_"+lzs(khit.rsd1,3)+"_"+tmp2.residue(7).name1()+lzs(khit.rsd2,3);
							//               //prev.dump_pdb("test_"+fn+".pdb");
							// #ifdef USE_OPENMP
							// #pragma omp critical
							// #endif
							//               tmp2.dump_pdb(utility::file_basename(infile)+"_"+fn+"_min.pdb");
							// #ifdef USE_OPENMP
							// #pragma omp critical
							// #endif
							//               nhit++;


							//               Real rd1 = tmp2.xyz(AtomID(n11,1)).distance(tmp2.xyz(AtomID(c11,5)));
							//               Real rd2 = tmp2.xyz(AtomID(n12,1)).distance(tmp2.xyz(AtomID(c12,5)));
							//               Real rd3 = tmp2.xyz(AtomID(n21,3)).distance(tmp2.xyz(AtomID(c21,5)));
							//               Real rd4 = tmp2.xyz(AtomID(n22,3)).distance(tmp2.xyz(AtomID(c22,5)));
							//               Real ra1 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C2"),5)),tmp2.xyz(AtomID(tmp.residue(5).atom_index("C7"),5)),tmp2.xyz(AtomID(9,1)));
							//               Real ra2 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C7"),5)),tmp2.xyz(AtomID(9,1)),tmp2.xyz(AtomID(n1a,1)));
							//               Real ra3 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C8"),5)),tmp2.xyz(AtomID(tmp.residue(5).atom_index("C9"),5)),tmp2.xyz(AtomID(9,3)));
							// 							Real ra4 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp.residue(5).atom_index("C9"),5)),tmp2.xyz(AtomID(9,3)),tmp2.xyz(AtomID(n2a,3)));
							//               Real ra5 = numeric::angle_degrees(tmp2.xyz(AtomID(bk11,2)),tmp2.xyz(AtomID(bk12,2)),tmp2.xyz(AtomID(bk13,1)));
							//               Real ra6 = numeric::angle_degrees(tmp2.xyz(AtomID(bk12,2)),tmp2.xyz(AtomID(bk13,1)),tmp2.xyz(AtomID(bk14,1)));
							//               Real ra7 = numeric::angle_degrees(tmp2.xyz(AtomID(bk21,4)),tmp2.xyz(AtomID(bk22,4)),tmp2.xyz(AtomID(bk23,3)));
							//               Real ra8 = numeric::angle_degrees(tmp2.xyz(AtomID(bk22,4)),tmp2.xyz(AtomID(bk23,3)),tmp2.xyz(AtomID(bk24,3)));
							//               Real kd1 = tmp2.xyz(AtomID(9                              ,6)).distance(tmp2.xyz(AtomID(iO2,5)));
							//               Real kd2 = tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)).distance(tmp2.xyz(AtomID(iO1,5)));
							//               Real ka1 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-2,7)),tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)), tmp2.xyz(AtomID(iO1,5)));
							// 							Real ka2 = numeric::angle_degrees(tmp2.xyz(AtomID(tmp2.residue(7).nheavyatoms()-1,7)),tmp2.xyz(AtomID(iO1,5)), tmp2.xyz(AtomID(iC5,5)));

							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							sf->score(tmp2);
							core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
							ss_out_all->fill_struct(tmp2,fn);
							ss_out_all->add_energy( "rd1", rd1 );
							ss_out_all->add_energy( "rd2", rd2 );
							ss_out_all->add_energy( "rd3", rd3 );
							ss_out_all->add_energy( "rd4", rd4 );
							ss_out_all->add_energy( "ra1", ra1 );
							ss_out_all->add_energy( "ra2", ra2 );
							ss_out_all->add_energy( "ra3", ra3 );
							ss_out_all->add_energy( "ra4", ra4 );
							ss_out_all->add_energy( "ra5", ra5 );
							ss_out_all->add_energy( "ra6", ra6 );
							ss_out_all->add_energy( "ra7", ra7 );
							ss_out_all->add_energy( "ra8", ra8 );
							ss_out_all->add_energy( "kd1", kd1 );
							ss_out_all->add_energy( "kd2", kd2 );
							ss_out_all->add_energy( "ka1", ka1 );
							ss_out_all->add_energy( "ka2", ka2 );
							ss_out_all->add_energy( "irep" ,
																			tmp2.energies().residue_total_energies(1)[core::scoring::fa_intra_rep]+
																			tmp2.energies().residue_total_energies(2)[core::scoring::fa_intra_rep]+
																			tmp2.energies().residue_total_energies(3)[core::scoring::fa_intra_rep]+
																			tmp2.energies().residue_total_energies(4)[core::scoring::fa_intra_rep]+
																			tmp2.energies().residue_total_energies(6)[core::scoring::fa_intra_rep]+
																			tmp2.energies().residue_total_energies(7)[core::scoring::fa_intra_rep] );

							ss_out_all->add_energy( "irep1" , tmp2.energies().residue_total_energies(1)[core::scoring::fa_intra_rep] );
							ss_out_all->add_energy( "irep2" , tmp2.energies().residue_total_energies(2)[core::scoring::fa_intra_rep] );
							ss_out_all->add_energy( "irep3" , tmp2.energies().residue_total_energies(3)[core::scoring::fa_intra_rep] );
							ss_out_all->add_energy( "irep4" , tmp2.energies().residue_total_energies(4)[core::scoring::fa_intra_rep] );
							ss_out_all->add_energy( "irep5" , tmp2.energies().residue_total_energies(6)[core::scoring::fa_intra_rep] );
							ss_out_all->add_energy( "irep6" , tmp2.energies().residue_total_energies(7)[core::scoring::fa_intra_rep] );
							#ifdef USE_OPENMP
							#pragma omp critical
							#endif
							sfd.write_silent_struct( *ss_out_all, utility::file_basename(infile)+".sc" );
							//tmp2.dump_pdb(utility::file_basename(infile)+"_"+fn+".pdb");

						}
					}
				}
			}
		}
	}

	TR << "MATCHES: " << nrhit << " " << nhit << std::endl;

	return 0;

	// //#ifdef USE_OPENMP
	// //#pragma omp parallel for schedule(dynamic,1)
	// //#endif
	// for(Size ixr = 1; ixr <= rgrid.size1(); ++ixr) {
	//   ScoreFunctionOP sf = new core::scoring::ScoreFunction();
	//   sf->set_weight(core::scoring::fa_atr,1.0);
	//   sf->set_weight(core::scoring::fa_rep,1.0);
	//   sf->set_weight(core::scoring::fa_dun,1.0);
	//   sf->set_weight(core::scoring::hbond_sc,10.0);
	//   sf->set_weight(core::scoring::atom_pair_constraint,1.0);
	//   sf->set_weight(core::scoring::    angle_constraint,1.0);
	//   core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	//   movemap->set_bb(false);
	//   movemap->set_chi(true);
	//   movemap->set_jump(false);
	//   movemap->set_jump(4,true);
	//   protocols::simple_moves::MinMover minm( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );

	//   TR << ixr << std::endl;
	//   for(Size iyr = 1; iyr <= rgrid.size2(); ++iyr) {
	//     for(Size izr = 1; izr <= rgrid.size3(); ++izr) {
	//       vector1<HitOP> const & rtmphit( rgrid(ixr,iyr,izr) );
	//       for(Size irh = 1; irh <= rtmphit.size(); ++irh) {
	//         Hit const & rhit(*rtmphit[irh]);
	//         Vec rcom = rhit.stub.local2global(ccom);
	//         Pose lg1(lg);
	//         for(Size i = 1; i <= lg1.residue_type(1).natoms(); ++i) lg1.set_xyz(AtomID(i,1), rhit.stub.local2global(lg1.xyz(AtomID(i,1))));
	//         Vec lg1O1(lg1.residue(1).xyz(iO1)),lg1C7(lg1.residue(1).xyz(iC7)),lg1C9(lg1.residue(1).xyz(iC9));

	//         for(Size ixs = max(ixr-1,((Size)1)); ixs <= min(ixr+((Size)1),rgrid.size1()); ++ixs) {
	//           for(Size iys = max(iyr-1,((Size)1)); iys <= min(iyr+((Size)1),rgrid.size2()); ++iys) {
	//             for(Size izs = max(izr-1,((Size)1)); izs <= min(izr+((Size)1),rgrid.size3()); ++izs) {
	//               vector1<HitOP> const & stmphit( sgrid(ixs,iys,izs) );
	//               //TR << ixr << " " << ixs << " " << rtmphit.size() << " " << stmphit.size() << std::endl;
	//               for(Size ish = 1; ish <= stmphit.size(); ++ish) {
	//                 Hit const & shit(*stmphit[ish]);
	//                 if(rhit.rsd1==shit.rsd1||rhit.rsd1==shit.rsd2||rhit.rsd2==shit.rsd1||rhit.rsd2==shit.rsd2) continue;
	//                 Vec scom = shit.stub.local2global(ccomB);

	//                 Pose lg2(lgB);
	//                 for(Size i = 1; i <= lg2.residue_type(1).natoms(); ++i) lg2.set_xyz(AtomID(i,1), shit.stub.local2global(lg2.xyz(AtomID(i,1))));

	//                 if(lg1O1.distance_squared(shit.stub.local2global(lgBO1)) > 4.25) continue;
	//                 if(lg1C7.distance_squared(shit.stub.local2global(lgBC7)) > 4.25) continue;
	//                 nhit++; TR << nhit << std::endl;
	//                 if(lg1C9.distance_squared(shit.stub.local2global(lgBC9)) > 4.25) continue;

	//                 continue;
	//                 rhit.apply(wp);
	//                 shit.apply(wp);
	//                 // check & refine
	//                 Pose tmp;
	//                 tmp.append_residue_by_jump(wp.residue(rhit.rsd1),1);
	//                 tmp.append_residue_by_jump(wp.residue(rhit.rsd2),1);
	//                 tmp.append_residue_by_jump(wp.residue(shit.rsd1),1);
	//                 tmp.append_residue_by_jump(wp.residue(shit.rsd2),1);
	//                 tmp.append_residue_by_jump(lg1.residue(1),1);

	//                 Size n11 = rhit.frnt?8:10, n12 = 11, n1a = rhit.frnt?10:8;
	//                 Size n21 = shit.frnt?8:10, n22 = 11, n2a = shit.frnt?10:8;
	//                 Size c11 = rhit.nego ? ctp.residue(1).atom_index("O6"): ctp.residue(1).atom_index("O5");
	//                 Size c12 = rhit.nego ? ctp.residue(1).atom_index("O5"): ctp.residue(1).atom_index("O6");
	//                 Size c21 = shit.nego ? ctp.residue(1).atom_index("O3"): ctp.residue(1).atom_index("O4");
	//                 Size c22 = shit.nego ? ctp.residue(1).atom_index("O4"): ctp.residue(1).atom_index("O3");
	//                 using namespace core::scoring::constraints;
	//                 tmp.add_constraint( new AtomPairConstraint( AtomID(n11,1), AtomID(c11,5), new HarmonicFunc(2.8,0.1)) );
	//                 tmp.add_constraint( new AtomPairConstraint( AtomID(n12,1), AtomID(c12,5), new HarmonicFunc(2.8,0.1)) );
	//                 tmp.add_constraint( new AtomPairConstraint( AtomID(n21,3), AtomID(c21,5), new HarmonicFunc(2.8,0.1)) );
	//                 tmp.add_constraint( new AtomPairConstraint( AtomID(n22,3), AtomID(c22,5), new HarmonicFunc(2.8,0.1)) );
	//                 // TR << "n11 " << tmp.residue(1).atom_name(n11) << " c11 " << tmp.residue(5).atom_name(c11) << std::endl;
	//                 // TR << "n12 " << tmp.residue(1).atom_name(n12) << " c12 " << tmp.residue(5).atom_name(c12) << std::endl;
	//                 // TR << "n21 " << tmp.residue(3).atom_name(n21) << " c21 " << tmp.residue(5).atom_name(c21) << std::endl;
	//                 // TR << "n22 " << tmp.residue(3).atom_name(n22) << " c22 " << tmp.residue(5).atom_name(c22) << std::endl;
	//                 tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C2"),5),AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),               new HarmonicFunc(PI,0.2)) );
	//                 tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),AtomID(n1a,1), new HarmonicFunc(PI,0.2)) );
	//                 tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C8"),5),AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),               new HarmonicFunc(PI,0.2)) );
	//                 tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),AtomID(n2a,3), new HarmonicFunc(PI,0.2)) );
	//                 // TR << "C2 C7 " << tmp.residue(1).atom_name(9) << " " << tmp.residue(1).atom_name(n1a) << std::endl;
	//                 // TR << "C8 C9 " << tmp.residue(3).atom_name(9) << " " << tmp.residue(3).atom_name(n2a) << std::endl;
	//                 Size bk11 = (tmp.residue(2).aa()==core::chemical::aa_glu)?6:5, bk12 = (tmp.residue(2).aa()==core::chemical::aa_glu)?7:6, bk13 = 9, bk14=n11;
	//                 tmp.add_constraint( new AngleConstraint( AtomID(bk11,2),AtomID(bk12,2),AtomID(bk13,1),                new HarmonicFunc(PI,0.2)) );
	//                 tmp.add_constraint( new AngleConstraint(                AtomID(bk12,2),AtomID(bk13,1),AtomID(bk14,1), new HarmonicFunc(PI,0.2)) );
	//                 Size bk21 = (tmp.residue(4).aa()==core::chemical::aa_glu)?6:5, bk22 = (tmp.residue(4).aa()==core::chemical::aa_glu)?7:6, bk23 = 9, bk24=n21;
	//                 tmp.add_constraint( new AngleConstraint( AtomID(bk21,4),AtomID(bk22,4),AtomID(bk23,3),                new HarmonicFunc(PI,0.2)) );
	//                 tmp.add_constraint( new AngleConstraint(                AtomID(bk22,4),AtomID(bk23,3),AtomID(bk24,3), new HarmonicFunc(PI,0.2)) );

	//                 static vector1<bool> seenit(16,false);
	//                 int type = 1 + 1*(int)rhit.frnt + 2*(int)rhit.nego + 4*(int)shit.frnt + 8*(int)shit.nego;
	//                 if( !seenit[type] ) {
	//                   seenit[type] = true;
	//                   minm.apply(tmp);
	//                   tmp.dump_pdb("test"+lzs(type,2)+".pdb");
	//                   int sum = 0; for(Size i = 1; i <= 16; i++) sum += seenit[i];
	//                   TR << "SUM " << sum << std::endl;
	//                   if(sum==16) utility_exit_with_message("lasdfj");
	//                 }
	//                 continue;
	//                 for(Size ix = max(ixr-1,((Size)1)); ix <= min(ixr+((Size)1),rgrid.size1()); ++ix) {
	//                   for(Size iy = max(iyr-1,((Size)1)); iy <= min(iyr+((Size)1),rgrid.size2()); ++iy) {
	//                     for(Size iz = max(izr-1,((Size)1)); iz <= min(izr+((Size)1),rgrid.size3()); ++iz) {
	//                       vector1<HitOP> const & ktmphit( kgrid(ix,iy,iz) );
	//                       for(Size ikh = 1; ikh <= ktmphit.size(); ++ikh) {
	//                         Hit const & khit(*ktmphit[ikh]);
	//                         if(khit.rsd1==rhit.rsd1||khit.rsd1==rhit.rsd2||khit.rsd2==rhit.rsd1||khit.rsd2==rhit.rsd2) continue;
	//                         if(khit.rsd1==shit.rsd1||khit.rsd1==shit.rsd2||khit.rsd2==shit.rsd1||khit.rsd2==shit.rsd2) continue;
	//                         Vec kcom = shit.stub.local2global(ccom);
	//                         if( rcom.distance_squared(kcom) > 2.25 ) continue;
	//                         if( scom.distance_squared(kcom) > 2.25 ) continue;
	//                         khit.apply(wp);

	//                         nhit++;
	//                         Pose tmp;
	//                         tmp.append_residue_by_jump(wp.residue(rhit.rsd1),1);
	//                         tmp.append_residue_by_jump(wp.residue(rhit.rsd2),1);
	//                         tmp.append_residue_by_jump(wp.residue(shit.rsd1),1);
	//                         tmp.append_residue_by_jump(wp.residue(shit.rsd2),1);
	//                         tmp.append_residue_by_jump(wp.residue(khit.rsd1),1);
	//                         tmp.append_residue_by_jump(wp.residue(khit.rsd2),1);
	//                         tmp.dump_pdb("test.pdb");
	//                         utility_exit_with_message("las;dfj");

	//                       }
	//                     }
	//                   }
	//                 }
	//               }
	//             }
	//           }
	//         }
	//       }
	//     }
	//   }
	// }

	// TR << "matching bkarg pairs & ring COO" << std::endl;
	// Size nhit=0;
	// //#ifdef USE_OPENMP
	// //#pragma omp parallel for schedule(dynamic,1)
	// //#endif
	// for(Size ih = 1; ih <= hits.size(); ih+=1) {
	//   //    ;    findhits_ih(frs, sasa, pose, asp, glu, ifc, hits, ih, dhit);
	//   Real DTH =  1.5;
	//   Real ATH = 25.0;

	//   if(ih%100==0) TR << ih << " " << rhit << " " << dhit << std::endl;
	//   Hit hi(*hits[ih]);

	//   Pose wp(pose);
	//   core::pack::dunbrack::SingleResidueRotamerLibraryCAP dlib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( asp.residue(1).type() );
	//   core::pack::dunbrack::SingleResidueRotamerLibraryCAP elib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( glu.residue(1).type() );
	//   core::pack::dunbrack::RotamerLibraryScratchSpace scratch;
	//   ScoreFunctionOP sf = new core::scoring::ScoreFunction();
	//   sf->set_weight(core::scoring::fa_atr,1.0);
	//   sf->set_weight(core::scoring::fa_rep,1.0);
	//   sf->set_weight(core::scoring::fa_dun,1.0);
	//   sf->set_weight(core::scoring::hbond_sc,10.0);
	//   sf->set_weight(core::scoring::atom_pair_constraint,1.0);
	//   sf->set_weight(core::scoring::    angle_constraint,1.0);
	//   core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	//   movemap->set_bb(false);
	//   movemap->set_chi(true);
	//   movemap->set_jump(false);
	//   movemap->set_jump(4,true);
	//   protocols::simple_moves::MinMover minm( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );

	//   Pose ctp; pose_from_pdb(ctp,"inrut/ctp.pdb");
	//   Size iO2=ctp.residue(1).atom_index("O2"), iC5=ctp.residue(1).atom_index("C5"), iO1=ctp.residue(1).atom_index("O1");

	//   Vec or2 = hi.ori.cross(hi.axs);

	//   for(Size negori = 0; negori < 2; negori++) {
	//     Vec ctp_axs = (ctp.residue(1).xyz("C2")-ctp.residue(1).xyz("C7")).normalized();
	//     Vec ctp_ori = ctp_axs.cross( (ctp.residue(1).xyz("O5") - ctp.residue(1).xyz("C7")).normalized() );
	//     if(negori) ctp_ori = -1.0 * ctp_ori;
	//     Mat R = numeric::alignVectorSets( ctp_axs, ctp_ori, hi.axs, hi.ori );
	//     rot_pose(ctp,R);
	//     trans_pose(ctp, hi.cen + 4.05*hi.axs - ctp.residue(1).xyz("C7") );
	//     // {
	//     //  bool clash = false;
	//     //  for(Size ia = 1; ia <= ctp.residue(1).nheavyatoms(); ia++) {
	//     //    if(!ifc.clash_check( ctp.residue(1).xyz(ia) ) ) clash = true;
	//     //  }
	//     //  if(clash) continue;
	//     // }

	//     for(Size jh = 1; jh <= hits.size(); jh+=1) {
	//       Hit hj(*hits[jh]);

	//       //if( hi.cen.distance_squared(hj.cen) > (13.6)*(13.6) ) continue;
	//       //if( hi.cen.distance_squared(hj.cen) < (10.5)*(10.5) ) continue;
	//       if( fabs( hi.axs.dot(hj.cen-hi.cen) - 10.455 ) > DTH ) continue;
	//       Real z = hi.ori.dot(hj.cen-hi.cen);
	//       Real y = or2   .dot(hj.cen-hi.cen);
	//       bool up = true;
	//       if( fabs(y-3.256) < DTH && fabs(z+4.945) < DTH ) up = false;
	//       if((fabs(y+3.256) > DTH || fabs(z-4.945) > DTH) && up ) continue;

	//       if( ctp.residue(1).xyz("C9").distance_squared(hj.cen+4.05*hj.axs) > DTH*DTH ) continue;

	//       Vec ctp_axs2 = (ctp.residue(1).xyz("C8")-ctp.residue(1).xyz("C9")).normalized();
	//       if( numeric::angle_degrees(ctp_axs2,Vec(0,0,0),hj.axs) > ATH ) continue;
	//       //Vec ctp_ori2 = ctp_axs.cross( (ctp.residue(1).xyz("O3") - ctp.residue(1).xyz("C9")).normalized() );

	//       Vec ctp_ori2 = ctp_axs2.cross( ctp.residue(1).xyz("O3") - ctp.residue(1).xyz("C9") );
	//       Real ooang = numeric::angle_degrees(ctp_ori2,Vec(0,0,0),hj.ori);

	//       if( ooang > ATH && 180.0-ooang > ATH ) continue;
	//       bool negori2 = ooang > ATH;

	//       Pose tmp;
	//       tmp.append_residue_by_jump(hi.rsd1,1);
	//       tmp.append_residue_by_jump(hi.rsd2,1);
	//       tmp.append_residue_by_jump(hj.rsd1,1);
	//       tmp.append_residue_by_jump(hj.rsd2,1);
	//       tmp.append_residue_by_jump(ctp.residue(1),1);
	//       //minimize
	//       Size n11 = hi.frnt?8:10, n12 = 11, n1a = hi.frnt?10:8;
	//       Size n21 = hj.frnt?8:10, n22 = 11, n2a = hj.frnt?10:8;
	//       Size c11 = negori ? ctp.residue(1).atom_index("O6"): ctp.residue(1).atom_index("O5");
	//       Size c12 = negori ? ctp.residue(1).atom_index("O5"): ctp.residue(1).atom_index("O6");
	//       Size c21 = negori2? ctp.residue(1).atom_index("O4"): ctp.residue(1).atom_index("O3");
	//       Size c22 = negori2? ctp.residue(1).atom_index("O3"): ctp.residue(1).atom_index("O4");
	//       using namespace core::scoring::constraints;
	//       tmp.add_constraint( new AtomPairConstraint( AtomID(n11,1), AtomID(c11,5), new HarmonicFunc(2.8,0.1)) );
	//       tmp.add_constraint( new AtomPairConstraint( AtomID(n12,1), AtomID(c12,5), new HarmonicFunc(2.8,0.1)) );
	//       tmp.add_constraint( new AtomPairConstraint( AtomID(n21,3), AtomID(c21,5), new HarmonicFunc(2.8,0.1)) );
	//       tmp.add_constraint( new AtomPairConstraint( AtomID(n22,3), AtomID(c22,5), new HarmonicFunc(2.8,0.1)) );
	//       // TR << "n11 " << tmp.residue(1).atom_name(n11) << " c11 " << tmp.residue(5).atom_name(c11) << std::endl;
	//       // TR << "n12 " << tmp.residue(1).atom_name(n12) << " c12 " << tmp.residue(5).atom_name(c12) << std::endl;
	//       // TR << "n21 " << tmp.residue(3).atom_name(n21) << " c21 " << tmp.residue(5).atom_name(c21) << std::endl;
	//       // TR << "n22 " << tmp.residue(3).atom_name(n22) << " c22 " << tmp.residue(5).atom_name(c22) << std::endl;
	//       tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C2"),5),AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),               new HarmonicFunc(PI,0.2)) );
	//       tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C7"),5),AtomID(9,1),AtomID(n1a,1), new HarmonicFunc(PI,0.2)) );
	//       tmp.add_constraint( new AngleConstraint( AtomID(tmp.residue(5).atom_index("C8"),5),AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),               new HarmonicFunc(PI,0.2)) );
	//       tmp.add_constraint( new AngleConstraint(                                           AtomID(tmp.residue(5).atom_index("C9"),5),AtomID(9,3),AtomID(n2a,3), new HarmonicFunc(PI,0.2)) );
	//       //        TR << "C2 C7 " << tmp.residue(1).atom_name(9) << " " << tmp.residue(1).atom_name(n1a) << std::endl;
	//       //        TR << "C8 C9 " << tmp.residue(3).atom_name(9) << " " << tmp.residue(3).atom_name(n2a) << std::endl;
	//       Size bk11 = (tmp.residue(2).aa()==core::chemical::aa_glu)?6:5, bk12 = (tmp.residue(2).aa()==core::chemical::aa_glu)?7:6, bk13 = 9, bk14=n11;
	//       tmp.add_constraint( new AngleConstraint( AtomID(bk11,2),AtomID(bk12,2),AtomID(bk13,1),                new HarmonicFunc(PI,0.2)) );
	//       tmp.add_constraint( new AngleConstraint(                AtomID(bk12,2),AtomID(bk13,1),AtomID(bk14,1), new HarmonicFunc(PI,0.2)) );
	//       Size bk21 = (tmp.residue(4).aa()==core::chemical::aa_glu)?6:5, bk22 = (tmp.residue(4).aa()==core::chemical::aa_glu)?7:6, bk23 = 9, bk24=n21;
	//       tmp.add_constraint( new AngleConstraint( AtomID(bk21,4),AtomID(bk22,4),AtomID(bk23,3),                new HarmonicFunc(PI,0.2)) );
	//       tmp.add_constraint( new AngleConstraint(                AtomID(bk22,4),AtomID(bk23,3),AtomID(bk24,3), new HarmonicFunc(PI,0.2)) );
	//       minm.apply(tmp);
	//       {
	//         bool clash = false;
	//         for(Size ia = 1; ia <= tmp.residue(5).nheavyatoms(); ia++) {
	//           if(!ifc.clash_check( tmp.residue(5).xyz(ia) ) ) clash = true;
	//         }
	//         if(clash) continue;
	//       }
	//       rhit++;

	//       Vec OCpose = tmp.residue(5).xyz(iO2);
	//       Vec COpose = tmp.residue(5).xyz(iC5);
	//       Vec OHpose = tmp.residue(5).xyz(iO1);
	//       Size ixr = std::ceil( COpose.x() - lb.x() );
	//       Size iyr = std::ceil( COpose.y() - lb.y() );
	//       Size izr = std::ceil( COpose.z() - lb.z() );

	//       for(Size ix = max(ixr-1,((Size)1)); ix <= min(ixr+((Size)1),grid.size1()); ++ix) {
	//         for(Size iy = max(iyr-1,((Size)1)); iy <= min(iyr+((Size)1),grid.size2()); ++iy) {
	//           for(Size iz = max(izr-1,((Size)1)); iz <= min(izr+((Size)1),grid.size3()); ++iz) {
	//             vector1<HitOP> & hitlst(grid(ix,iy,iz));
	//             for(Size ig = 1; ig <= hitlst.size(); ++ig) {
	//               Hit & hk(*hitlst[ig]);
	//               // Vec const & CO(hk.ori); if(COpose.distance_squared(CO) > 1.6) continue;
	//               // Vec const & OC(hk.cen); if(OCpose.distance_squared(OC) > 1.6) continue;
	//               // Vec const & OH(hk.axs); if(OHpose.distance_squared(OH) > 1.6) continue;
	//               // Pose tmp2(tmp);
	//               utility_exit_with_message("lasdfkj");
	//               //                 tmp2.append_residue_by_jump(hk.rsd1,1);
	//               //                 tmp2.append_residue_by_jump(hk.rsd2,1);
	//               //                 string fn = "IKM_"+lzs(hi.rsd1.seqpos(),4)+"_"+lzs(hj.rsd1.seqpos(),4)+"_"+lzs(hk.rsd1.seqpos(),4)+"_";
	//               //                 fn += lzs(hi.rsd2.seqpos(),4)+"_"+lzs(hj.rsd2.seqpos(),4)+"_"+lzs(hk.rsd2.seqpos(),4)+"_";
	//               //                 fn += (hi.frnt?"s":"f");fn += (negori ?"d":"u");fn += "_";
	//               //                 fn += (hj.frnt?"s":"f");fn += (negori2?"d":"u");fn += "_"+lzs(ih,5)+"_"+lzs(jh,5)+"_";
	//               //                 fn += lzs(ix,3)+"-"+lzs(iy,3)+"-"+lzs(iz,3);

	//               //                 dhit++;

	//               //                 {
	//               //                   sf->score(tmp2);
	//               //                   core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
	//               //                   ss_out_all->fill_struct(tmp2,fn);
	//               //                   ss_out_all->add_energy( "frntside1" , hi.frnt );
	//               //                   ss_out_all->add_energy( "frntside2" , hj.frnt );
	//               //                   sfd.write_silent_struct( *ss_out_all, option[ basic::options::OptionKeys::out::file::silent ]() );
	//               //                   tmp2.dump_pdb(fn+".pdb");
	//               //                 }
	//               // tmp.add_constraint( new AtomPairConstraint( AtomID(iO2,5), AtomID(9,6), new HarmonicFunc(3.0,0.1)) );
	//               //                if(pose.residue(7).name3()=="ASP") tmp.add_constraint( new AtomPairConstraint( AtomID(iO1,5), AtomID(7,6), new HarmonicFunc(3.0,0.1)) );
	//               //                if(pose.residue(7).name3()=="GLU") tmp.add_constraint( new AtomPairConstraint( AtomID(iO1,5), AtomID(8,6), new HarmonicFunc(3.0,0.1)) );
	//               //                 minm.apply(tmp2);
	//               //                 {
	//               //                   sf->score(tmp2);
	//               //                   core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
	//               //                   ss_out_all->fill_struct(tmp2,fn+"_min");
	//               //                   ss_out_all->add_energy( "frntside1" , hi.frnt );
	//               //                   ss_out_all->add_energy( "frntside2" , hj.frnt );
	//               //                   sfd.write_silent_struct( *ss_out_all, option[ basic::options::OptionKeys::out::file::silent ]() );
	//               //                   tmp2.dump_pdb(fn+"_min.pdb");
	//               //                }
	//             }
	//           }
	//         }
	//       }
	//     }
	//   }
	//  }
	TR << "MATCHES: " << nhit << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


// old crap
	// if(hits.size()){
	//  Hit & h(*hits[hits.size()]);
	//   ozstream out("test.pdb");
	//   Size ano = 0;
	//   core::io::pdb::dump_pdb_residue(h.rsd2,ano,out);
	//   core::io::pdb::dump_pdb_residue(h.rsd1,ano,out);
	//  {Vec V=h.cen      ;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
	//  {Vec V=h.axs+h.cen;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
	//  {Vec V=h.ori+h.cen;out<<"HETATM"<<I(5,11)<<' '<<" VIZ"<<' '<<"VIZ"<<' '<<"A"<<I(4,6)<<"    "<<F(8,3,V.x())<<F(8,3,V.y())<<F(8,3,V.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';}
	//   out.close();
	//  utility_exit_with_message("asl;dfj");
	// }

	// TR << "precalc glu " << std::endl;
	// vector1<Real> CHI1,CHI2;
	// for(Real x = 4.5; x <=360.0; x += 5.0) CHI1.push_back(x);
	// for(Real x = 4.5; x <=360.0; x += 9.0) CHI2.push_back(x);
	// // precompute acceptable chis
	// ObjexxFCL::FArray1D<Real> eallow0(                        pose.n_residue(),9e9);
	// ObjexxFCL::FArray2D<Real> eallow1(            CHI1.size(),pose.n_residue(),9e9);
	// ObjexxFCL::FArray3D<Real> eallow2(CHI2.size(),CHI1.size(),pose.n_residue(),9e9);
	// {
	//   Pose wp(pose);
	//   core::pack::dunbrack::SingleResidueRotamerLibraryCAP dlib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( asp.residue(1).type() );
	//   core::pack::dunbrack::SingleResidueRotamerLibraryCAP elib = core::pack::dunbrack::RotamerLibrary::get_instance()->get_rsd_library( glu.residue(1).type() );
	//   core::pack::dunbrack::RotamerLibraryScratchSpace scratch;
	//   for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
	//     if(sasa[ir] > 0) continue;
	//     wp.replace_residue(ir,asp.residue(1),true); asp.replace_residue(1,wp.residue(ir),false);
	//     wp.replace_residue(ir,glu.residue(1),true); glu.replace_residue(1,wp.residue(ir),false);
	//     for(Real i1 = 1; i1 <= CHI1.size(); i1++) {
	//       Real kch1 = CHI1[i1];
	//       for(Size i2 = 1; i2 <= CHI2.size(); i2++) {
	//         Real kch2 = CHI2[i2];
	//         glu.set_chi(2,1,i2);
	//         bool clash = false;
	//         for(Size ia = 6; ia <= glu.residue(1).nheavyatoms(); ia++) {
	//           if(!ifc.clash_check(glu.residue(1).xyz(ia),ir)) clash=true;
	//         }
	//         if(clash) {
	//           eallow2(i2,i1,ir) = 9e9;
	//         } else {
	//           Real mndn = 9e9;
	//           for(Size i3 = 1; i3 <= CHI2.size(); i3++) {
	//             Real kch3 = CHI2[i3];
	//             glu.set_chi(3,1,i3);
	//             Real tmp = elib->rotamer_energy(glu.residue(1),scratch);
	//             if( tmp < mndn ) mndn = tmp;
	//           }
	//           eallow2(i2,i1,ir) = mndn;
	//         }
	//         if( eallow1(i1,ir) > eallow2(i2,i1,ir) ) eallow1(i1,ir) = eallow2(i2,i1,ir);
	//       }
	//       if( eallow0(ir) > eallow1(i1,ir) ) eallow0(ir) = eallow1(i1,ir);
	//     }
	//   }
	// }

	// Vec const OH(tmp.residue(5).xyz("O1"));
	//         Vec const C5(tmp.residue(5).xyz("C5"));
	//         for(Size kr = 1; kr <= wp.n_residue(); ++kr) {
	//           if(sasa[kr] > 0) continue;
	//           //if(eallow0(kr) > 10.0) continue;
	//           if( wp.residue(kr).xyz(5).distance_squared(OH) > 36.0 ) continue;
	//           wp.replace_residue(kr,asp.residue(1),true); asp.replace_residue(1,wp.residue(kr),false);
	//           wp.replace_residue(kr,glu.residue(1),true); glu.replace_residue(1,wp.residue(kr),false);
	//           bool asphit=false,gluhit=false;
	//           for(Real ikch1 = 1; ikch1 <= CHI1.size(); ikch1++) {
	//             Real kch1 = CHI1[ikch1];
	//             { // asp check
	//               if(asphit) goto gluchk;
	//               asp.set_chi(1,1,kch1);
	//               if( asp.residue(1).xyz(6).distance_squared(OH) > 21.0 ) goto gluchk;
	//               Vec const CG(asp.residue(1).xyz(6));
	//               Real const gd2 = OH.distance_squared(CG);
	//               if( 11.15 > gd2 || gd2 > 18.9 ) goto gluchk;
	//               asp.set_chi(2,1, numeric::dihedral_degrees(OH,CG,asp.residue(1).xyz(5),asp.residue(1).xyz(2)) );
	//               Vec OD(asp.residue(1).xyz(7));
	//               Real od2 = OH.distance_squared(OD);
	//               if( 2.3*2.3 > od2 || od2 > 3.3*3.3 ) goto gluchk;
	//               Real tmpa1 = numeric::angle_degrees(C5,OH,OD);
	//               if( fabs(tmpa1-105.3) > 12.0 ) goto gluchk;
	//               Real tmpa2 = numeric::angle_degrees(OH,OD,CG);
	//               if( fabs(tmpa2-120.0) > 12.0 ) goto gluchk;
	//               bool clash = false;
	//               for(Size ia = 1; ia <= tmp.residue(5).nheavyatoms(); ia++) {
	//                 if( tmp.residue(5).xyz(ia).distance_squared(asp.residue(1).xyz(5)) < 3.2*3.2 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(asp.residue(1).xyz(6)) < 3.2*3.2 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(asp.residue(1).xyz(7)) < 2.8*2.8 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(asp.residue(1).xyz(8)) < 2.8*2.8 ) clash = true;
	//               }
	//               if(clash) goto gluchk;
	//               if( dlib->rotamer_energy(asp.residue(1),scratch) > 4.0 ) continue;
	//               asphit=true;
	//               dhit++;
	//               Pose tmp2(tmp);
	//               Vec const H5(tmp2.residue(5).xyz("H5"));
	//               Mat Rtmp = numeric::rotation_matrix_degrees( C5-OH, numeric::dihedral_degrees(H5,OH,C5,OD) );
	//               tmp2.set_xyz(AtomID(tmp2.residue(5).atom_index("H5"),5), Rtmp*(H5-OH)+OH );
	//               tmp2.append_residue_by_jump(asp.residue(1),1);
	//               tmp2.dump_pdb(fn+"_d"+lzs(kr,4)+"_"+lzs(Size(kch1),3)+".pdb");
	//               sf->score(tmp2);
	//               core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
	//               ss_out_all->fill_struct(tmp2,fn+"_e"+lzs(kr,4)+"_"+lzs(Size(kch1),3)+".pdb");
	//               ss_out_all->add_energy( "frntside1" , hi.frnt );
	//               ss_out_all->add_energy( "frntside2" , hj.frnt );
	//               sfd.write_silent_struct( *ss_out_all, option[ basic::options::OptionKeys::out::file::silent ]() );
	//             }
	//           gluchk:
	//             //if(eallow1(ikch1,kr) > 6.0) continue;
	//             if(gluhit) continue;
	//             glu.set_chi(1,1,kch1);
	//             for(Size ikch2 = 1; ikch2 <= CHI2.size(); ikch2++) {
	//               //if(eallow2(ikch2,ikch1,kr) > 10.0) continue;
	//               Real kch2 = CHI2[ikch2];
	//               glu.set_chi(2,1,kch2);
	//               // check E
	//               if( glu.residue(1).xyz(7).distance_squared(OH) > 21.0 ) continue;
	//               Vec const CD(glu.residue(1).xyz(7));
	//               Real const dd2 = OH.distance_squared(CD);
	//               if( 11.15 > dd2 || dd2 > 18.9 ) continue;
	//               glu.set_chi(3,1, numeric::dihedral_degrees(OH,CD,glu.residue(1).xyz(6),glu.residue(1).xyz(5)) );
	//               Vec OE(glu.residue(1).xyz(8));
	//               Real od2 = OH.distance_squared(OE);
	//               if( 2.3*2.3 > od2 || od2 > 3.3*3.3 ) continue;
	//               Real tmpa1 = numeric::angle_degrees(C5,OH,OE);
	//               if( fabs(tmpa1-105.3) > 12.0 ) continue;
	//               Real tmpa2 = numeric::angle_degrees(OH,OE,CD);
	//               if( fabs(tmpa2-120.0) > 12.0 ) continue;
	//               if( elib->rotamer_energy(glu.residue(1),scratch) > 6.0 ) continue;
	//               bool clash = false;
	//               for(Size ia = 1; ia <= tmp.residue(5).nheavyatoms(); ia++) {
	//                 if( tmp.residue(5).xyz(ia).distance_squared(glu.residue(1).xyz(5)) < 3.2*3.2 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(glu.residue(1).xyz(6)) < 3.2*3.2 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(glu.residue(1).xyz(7)) < 3.2*3.2 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(glu.residue(1).xyz(8)) < 2.8*2.8 ) clash = true;
	//                 if( tmp.residue(5).xyz(ia).distance_squared(glu.residue(1).xyz(9)) < 2.8*2.8 ) clash = true;
	//               }
	//               if(clash) continue;
	//               gluhit=true;
	//               dhit++;
	//               Pose tmp2(tmp);
	//               Vec const H5(tmp2.residue(5).xyz("H5"));
	//               Mat Rtmp = numeric::rotation_matrix_degrees( C5-OH, numeric::dihedral_degrees(H5,OH,C5,OE) );
	//               tmp2.set_xyz(AtomID(tmp2.residue(5).atom_index("H5"),5), Rtmp*(H5-OH)+OH );
	//               tmp2.append_residue_by_jump(glu.residue(1),1);
	//               tmp2.dump_pdb(fn+"_e"+lzs(kr,4)+"_"+lzs(Size(kch1),3)+"_"+lzs(Size(kch2),3)+".pdb");
	//               sf->score(tmp2);
	//               core::io::silent::SilentStructOP ss_out_all( new core::io::silent::ScoreFileSilentStruct );
	//               ss_out_all->fill_struct(tmp2,fn+"_e"+lzs(kr,4)+"_"+lzs(Size(kch1),3)+"_"+lzs(Size(kch2),3)+".pdb");
	//               ss_out_all->add_energy( "frntside1" , hi.frnt );
	//               ss_out_all->add_energy( "frntside2" , hj.frnt );
	//               sfd.write_silent_struct( *ss_out_all, option[ basic::options::OptionKeys::out::file::silent ]() );
	//             }
	//           }
	//         }
