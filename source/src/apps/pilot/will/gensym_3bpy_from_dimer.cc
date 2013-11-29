// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/graph/Graph.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/optimizeH.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <devel/init.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/scoring/ImplicitFastClashCheck.hh>
#include <protocols/design_opt/GreedyOptMutationMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/toolbox/task_operations/JointSequenceOperation.hh>
#include <sstream>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/mynamespaces.ihh>
#include <apps/pilot/will/will_util.ihh>

#ifdef USE_OPENMP
#include <omp.h>
#endif
int num_threads() {
	#ifdef USE_OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
}
int thread_num() {
	#ifdef USE_OPENMP
		return omp_get_thread_num();
	#else
		return 0;
	#endif
}

// OPT_1GRP_KEY( Integer, bpytoi, min_bpy_cb_mono_nbr_count   );
// OPT_1GRP_KEY( Integer, bpytoi, min_bpy_bb_tri_nbr_count   );
// OPT_1GRP_KEY( Integer, bpytoi, min_bpy_bb_tri_nbr_count   );
OPT_1GRP_KEY( Integer, bpytoi, chi1_increment  )
OPT_1GRP_KEY( Integer, bpytoi, max_nres        )
OPT_1GRP_KEY( Integer, bpytoi, max_cys         )
OPT_1GRP_KEY( Real   , bpytoi, max_bpy_dun     )
OPT_1GRP_KEY( Real   , bpytoi, bpy_clash_dis   )
OPT_1GRP_KEY( Real   , bpytoi, max_sym_error   )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( bpytoi::chi1_increment            ,"", 1   );
	NEW_OPT( bpytoi::max_nres                  ,"", 300 );
	NEW_OPT( bpytoi::max_cys                   ,"", 3   );
	NEW_OPT( bpytoi::max_bpy_dun               ,"", 3.0 );
	NEW_OPT( bpytoi::bpy_clash_dis             ,"", 3.0 );
	NEW_OPT( bpytoi::max_sym_error             ,"", 0.7 );
}

using core::kinematics::Stub;
using protocols::scoring::ImplicitFastClashCheck;

static basic::Tracer TR("gensym_3bpy_from_dimer");
static core::io::silent::SilentFileData sfd;

inline Real const sqr(Real const r) { return r*r; }
inline Real sigmoidish_neighbor( Real const & sqdist ) {
	if( sqdist > 9.*9. ) {
		return 0.0;
	} else if( sqdist < 6.*6. ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0	- sqr( (dist - 6.) / (9. - 6.) ) );
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
vector1<Vec> line_cone_intersection(Vec p, Vec d, Vec v, Vec a, Real t) {
	vector1<Vec> sol;
	t = numeric::conversions::radians(t);
	Mat M = numeric::outer_product(a,a) - cos(t)*cos(t)*Mat::identity();
	Vec D = p-v;
	Real c2 =	 d.dot(M*d);
	Real c1 = 2*d.dot(M*D);
	Real c0 =	 D.dot(M*D);
	Real disc = c1*c1 - 4*c0*c2;
	if( disc == 0) sol.push_back( p + (-c1)/(2.0*c2)*d );
	else if( disc > 0) {
		disc = sqrt(disc);
		sol.push_back(p+(-c1+disc)/(2.0*c2)*d);
		sol.push_back(p+(-c1-disc)/(2.0*c2)*d);
	}
	return sol;
}
vector1<std::pair<Vec,Vec> > intersecting_bpy_axes(Vec CB, Vec CG, Vec FE, Vec symmaxis, Vec symmcen = Vec(0,0,0)) {
	vector1<std::pair<Vec,Vec> > sol;
	Vec	p = symmcen;
	Vec	d = symmaxis.normalized();
	Vec	a = (CG-CB).normalized();
	Vec	v = CG + projection_matrix(a)*(FE-CG);
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
// void fixH(core::pose::Pose & pose) {
// 	for(Size i = 1; i <= pose.n_residue(); ++i) {
// 		if(!pose.residue(i).has("H")) continue;
// 		numeric::xyzVector<Real> n	= pose.residue(i).xyz("N");
// 		numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
// 		Size in = i-1;
// 		if(in == 0) in = pose.n_residue();
// 		numeric::xyzVector<Real> c	= pose.residue(in).xyz("C");
// 		numeric::xyzVector<Real> h	= n + (n-(ca+c)/2.0).normalized()*1.01;
// 		pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
// 	}
// }
bool is_near_C2Z_iface(Pose const & p, Size rsd) {
	Mat Rc2 = rotation_matrix_degrees(Vec(0,0,1),180.0);
	Vec CA = p.xyz(AtomID(2,rsd));
	for (Size ir=1; ir < p.n_residue(); ++ir) {
		if(p.residue(ir).is_protein())
			if( CA.distance_squared( Rc2*p.xyz(AtomID(2,ir))) < 100.0 )
				return true;
	}
	return false;
}

void fixbb_design(Pose & pose, Size ibpy, Size dsub) {
	using namespace core::chemical;
	using namespace core::conformation::symmetry;
	using namespace core::pack::task;
	using namespace core::scoring::constraints;
	ScoreFunctionOP sf = core::scoring::getScoreFunction();
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	vector1<bool> allowed_aas(20,false);
	allowed_aas[aa_ala] = true;
	allowed_aas[aa_phe] = true;
	allowed_aas[aa_ile] = true;
	allowed_aas[aa_leu] = true;
	allowed_aas[aa_met] = true;
	allowed_aas[aa_val] = true;
	allowed_aas[aa_trp] = true;

	vector1<bool> allowed_all(20,true);

	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	std::set<Size> design_pos;
	for( Size i=1; i <= sym_info->num_independent_residues() ; i++) {
		if(i==ibpy) continue;
		if(!sym_info->bb_is_independent(i)) continue;
		if(pose.residue(i).name3()=="GLY"||pose.residue(i).name3()=="PRO") continue;
		for( Size j=sym_info->num_independent_residues()+1; j <= sym_info->num_total_residues_without_pseudo(); j++) {
			if(sym_info->subunit_index(j) == dsub ) continue;
			if(pose.residue(j).name3()=="GLY"||pose.residue(j).name3()=="PRO") continue;
			for(int ia=1; ia <= pose.residue(i).nheavyatoms();++ia){
				for(int ja=1; ja <= pose.residue(j).nheavyatoms();++ja){
					if( pose.residue(i).xyz(ia).distance_squared( pose.residue(j).xyz(ja) ) < 36.0 ) design_pos.insert(i);
				}
			}
		}
	}
	std::set<Size> design_pos_bpy;
	std::set<Size> design_pos_bpy_ex;
	for( Size i=1; i <= sym_info->num_independent_residues() ; i++) {
		if(i==ibpy) continue;
		if(!sym_info->bb_is_independent(i)) continue;
		for( Size j=1; j <= sym_info->num_total_residues_without_pseudo(); j++) {
			if(pose.residue(j).name3()!="BPY") continue;
			for(int ja=6; ja <= pose.residue(j).nheavyatoms();++ja){
				if( pose.residue(i).has("CB") ){
					if(pose.residue(i).xyz("CB").distance_squared( pose.residue(j).xyz(ja) ) < 81.0 ){
						design_pos_bpy.insert(i);
					}
				}
				for(int ia=5; ia <= pose.residue(i).nheavyatoms();++ia){
					if( pose.residue(i).xyz(ia).distance_squared( pose.residue(j).xyz(ja) ) < 36.0 ) design_pos_bpy_ex.insert(i);
				}
			}
		}
	}

	for( Size i=1; i<=pose.n_residue(); i++) {
		if (!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY" || i==ibpy) {
			//TR << "res " << i << " fix" << std::endl;
			// Don't mess with Pros or Glys at the interfaces
			task->nonconst_residue_task(i).prevent_repacking();
		} else if (find(design_pos_bpy_ex.begin(), design_pos_bpy_ex.end(), i) != design_pos_bpy_ex.end()) {
			std::cout << "design_pos_bpy_ex " << i << " " << pose.residue(i).name3() << std::endl;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
		//	EX_ONE_STDDEV,                 //1
		//	EX_ONE_HALF_STEP_STDDEV,       //2
		//	EX_TWO_FULL_STEP_STDDEVS,      //3
		//	EX_TWO_HALF_STEP_STDDEVS,      //4
		//	EX_FOUR_HALF_STEP_STDDEVS,     //5
		//	EX_THREE_THIRD_STEP_STDDEVS,   //6
		//	EX_SIX_QUARTER_STEP_STDDEVS,   //7
			task->nonconst_residue_task(i).or_ex1_sample_level( EX_TWO_FULL_STEP_STDDEVS );
			task->nonconst_residue_task(i).or_ex2_sample_level( EX_TWO_FULL_STEP_STDDEVS );
		} else if (find(design_pos_bpy.begin(), design_pos_bpy.end(), i) != design_pos_bpy.end()) {
			std::cout << "design_pos_bpy " << i << " " << pose.residue(i).name3() << std::endl;
			bool temp = allowed_aas[pose.residue(i).aa()];
			allowed_aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_from_command_line();
			allowed_aas[pose.residue(i).aa()] = temp;
		} else if (find(design_pos.begin(), design_pos.end(), i) != design_pos.end()) {
			std::cout << "design_pos " << i << " " << pose.residue(i).name3() << std::endl;
			bool temp = allowed_aas[pose.residue(i).aa()];
			allowed_aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_all);
			task->nonconst_residue_task(i).or_include_current(true);
			task->nonconst_residue_task(i).initialize_from_command_line();
			allowed_aas[pose.residue(i).aa()] = temp;
		} else {
			//TR << "res " << i << " fix" << std::endl;
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}

	Real worig = sf->get_weight(core::scoring::res_type_constraint);
	if( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,1.0);
	utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);
	pose.add_constraints( res_cst );

	core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

	pose.remove_constraints( res_cst );
	sf->set_weight(core::scoring::res_type_constraint,worig);
}

void refine(Pose & pose, Size ibpy, Size dsub) {
	using namespace core::chemical;
	using namespace core::conformation::symmetry;
	using namespace core::pack::task;
	using namespace core::scoring::constraints;

  const core::Real PI = numeric::NumericTraits<Real>::pi();

	ScoreFunctionOP sf = core::scoring::getScoreFunction();
	sf->set_weight(core::scoring::atom_pair_constraint ,10.0);
	sf->set_weight(core::scoring::angle_constraint		 ,10.0);
	sf->set_weight(core::scoring::coordinate_constraint,10.0);
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	{
		// Set allowed AAs.
		vector1<bool> allowed_aas(20,false);
		allowed_aas[aa_ala] = true;

		PackerTaskOP task( TaskFactory::create_packer_task( pose ));

		vector1<Size> design_pos;
		for( Size i=1; i<=pose.n_residue(); i++) {
			if(i==ibpy) continue;
			if(!sym_info->bb_is_independent(i)) continue;
			for( Size j=1; j<=pose.n_residue(); j++) {
				if(sym_info->bb_is_independent(j)) continue;
				if(sym_info->subunit_index(j) == dsub ) continue;
				if( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) < 100.0 ) {
					//TR << "design_pos " << i << " " << j << std::endl;
					if( find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) design_pos.push_back(i);
				}
			}
		}

		// Set which residues can be designed
		for( Size i=1; i<=pose.n_residue(); i++) {
			if (!sym_info->bb_is_independent(i)) {
				task->nonconst_residue_task(i).prevent_repacking();
			} else if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
				// Don't mess with Pros or Glys at the interfaces
				task->nonconst_residue_task(i).prevent_repacking();
			} else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
				task->nonconst_residue_task(i).prevent_repacking();
			} else {
				bool temp = allowed_aas[pose.residue(i).aa()];
				allowed_aas[pose.residue(i).aa()] = true;
				task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
				task->nonconst_residue_task(i).or_include_current(true);
				task->nonconst_residue_task(i).initialize_from_command_line();
				allowed_aas[pose.residue(i).aa()] = temp;
			}
		}

		Real worig = sf->get_weight(core::scoring::res_type_constraint);
		if( worig == 0.0 ) sf->set_weight(core::scoring::res_type_constraint,10.0);
		utility::vector1< core::scoring::constraints::ConstraintCOP > res_cst = add_favor_native_cst(pose);
		pose.add_constraints( res_cst );

		// Actually perform design.
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
		protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
		packer->apply(pose);

		pose.remove_constraints( res_cst );
		sf->set_weight(core::scoring::res_type_constraint,worig);

	}
	{
		pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+0*sym_info->num_independent_residues()),
																								 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+1*sym_info->num_independent_residues()),
																								 new HarmonicFunc(0,0.02) ) );
		pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+0*sym_info->num_independent_residues()),
																								 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+2*sym_info->num_independent_residues()),
																								 new HarmonicFunc(0,0.02) ) );

		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(PI,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );

		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(PI,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
	}
	AtomID REF(1,sym_info->num_total_residues_without_pseudo()+1);
	for(Size ir = 1; ir <= sym_info->num_independent_residues(); ++ir) {
		if(ir==ibpy) continue;
		for(Size ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
			pose.add_constraint(new CoordinateConstraint( AtomID(ia,ir), REF, pose.xyz(AtomID(ia,ir)), new HarmonicFunc(0.0,1.0) ) );
		}
	}
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(true);
	movemap->set_bb(true);
	movemap->set_chi(true);
	core::pose::symmetry::make_symmetric_movemap(pose,*movemap);
	protocols::simple_moves::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
	pose.remove_constraints();


	fixbb_design(pose,ibpy,dsub);


	{
		pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+0*sym_info->num_independent_residues()),
																								 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+1*sym_info->num_independent_residues()),
																								 new HarmonicFunc(0,0.02) ) );
		pose.add_constraint( new AtomPairConstraint( AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+0*sym_info->num_independent_residues()),
																								 AtomID(pose.residue(ibpy).atom_index("ZN"),ibpy+2*sym_info->num_independent_residues()),
																								 new HarmonicFunc(0,0.02) ) );

		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(PI,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );

		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+1*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NE1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(PI,0.1) ) );
		pose.add_constraint( new AngleConstraint( AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("ZN" ),ibpy+0*sym_info->num_independent_residues()),
																							AtomID(pose.residue(ibpy).atom_index("NN1"),ibpy+2*sym_info->num_independent_residues()),
																							new HarmonicFunc(1.570796,0.1) ) );
	}
	protocols::simple_moves::symmetry::SymMinMover( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false ).apply(pose);
	pose.remove_constraints();

}
void repack(Pose & pose, Size ibpy, Size dsub) {
	using namespace core::chemical;
	using namespace core::conformation::symmetry;
	using namespace core::pack::task;
	using namespace core::scoring::constraints;
	ScoreFunctionOP sf = core::scoring::getScoreFunction();

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	vector1<Size> design_pos;
	for( Size i=1; i<=pose.n_residue(); i++) {
		if(i==ibpy) continue;
		if(!sym_info->bb_is_independent(i)) continue;
		for( Size j=1; j<=pose.n_residue(); j++) {
			if(sym_info->bb_is_independent(j)) continue;
			if(sym_info->subunit_index(j) == dsub ) continue;
			if( pose.residue(i).nbr_atom_xyz().distance_squared( pose.residue(j).nbr_atom_xyz() ) < 100.0 ) {
				//TR << "design_pos " << i << " " << j << std::endl;
				if( find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) design_pos.push_back(i);
			}
		}
	}

	// Set which residues can be designed
	for( Size i=1; i<=pose.n_residue(); i++) {
		if (!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if( pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
			// Don't mess with Pros or Glys at the interfaces
			task->nonconst_residue_task(i).prevent_repacking();
		} else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			task->nonconst_residue_task(i).restrict_to_repacking();
		}
	}

	// Actually perform design.
	core::pack::make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);
}
void new_sc(Pose &pose, utility::vector1<Size> intra_subs, Real& int_area, Real& sc) {

	using namespace core;

	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	core::scoring::sc::ShapeComplementarityCalculator scc;
	scc.Init();

	// Figure out which chains touch chain A, and add the residues from those chains
	// into the sc surface objects
	Size nres_monomer = symm_info->num_independent_residues();
	for (Size i=1; i<=nres_monomer; ++i) {
		scc.AddResidue(0, pose.residue(i));
	}
	for (Size i=1; i<=symm_info->subunits(); ++i) {
		if (std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end()) continue;
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for (Size ir=1; ir<=nres_monomer; ir++) {
			if (pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0) {
				contact = true;
				break;
			}
		}
		if (contact) {
			for (Size ir=1; ir<=nres_monomer; ir++) {
				scc.AddResidue(1, pose.residue(ir+start));
			}
		}
	}
	if (scc.Calc()) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	} else {
		sc = 0.0;
		int_area = 0.0;
	}
}
Size num_trimer_contacts(Pose const & psym, Size nres) {
	Size count = 0;
	for(Size ir = 1; ir <= nres; ++ir) {
		for(Size jr = nres+1; jr <= 3*nres; ++jr) {
			if( psym.residue(ir).nbr_atom_xyz().distance_squared( psym.residue(jr).nbr_atom_xyz() ) < 100 ) count++;
		}
		for(Size jr = 4*nres+1; jr <= 12*nres; ++jr) {
			if( psym.residue(ir).nbr_atom_xyz().distance_squared( psym.residue(jr).nbr_atom_xyz() ) < 100 ) count++;
		}
	}
	return count;
}
// Real get_bare_rep(Pose tmp, ScoreFunctionOP sfrepsym) {
// 	option[basic::options::OptionKeys::symmetry::symmetry_definition]("input/sym/C2.sym");
// 	core::pose::symmetry::make_symmetric_pose(tmp);
// 	sfrepsym->score(tmp);
// 	return tmp.energies().total_energies()[core::scoring::fa_rep];
// }
Real dimer_rms(Pose const & psym, Pose const & pdimer) {
	Pose dimer2;
	for(Size i = 1; i <= pdimer.n_residue()*2; ++i) {
		if(pdimer.n_residue()/2 < i && i <= pdimer.n_residue()*3/2 ) continue;
		if(psym.residue(i).is_lower_terminus()) dimer2.append_residue_by_jump(psym.residue(i),1);
		else                                    dimer2.append_residue_by_bond(psym.residue(i));
	}
	// dimer2.dump_pdb("dimer2.pdb");
	// pdimer.dump_pdb("pdimer.pdb");
	return core::scoring::CA_rmsd(dimer2,pdimer);
}
int neighbor_count(Pose const &pose, int ires, double distance_threshold=10.0) {
	core::conformation::Residue const resi( pose.residue( ires ) );
	if(!resi.has("CB")) return 0;
	Size resi_neighbors( 0 );
	for(Size jres = 1; jres <= pose.n_residue(); ++jres) {
		core::conformation::Residue const resj( pose.residue( jres ) );
		if(!resj.has("CB")) continue;
		double const distance = resi.xyz("CB").distance(resj.xyz("CB"));
		if( distance <= distance_threshold ){
			++resi_neighbors;
		}
	}
	return resi_neighbors;
}



#define ATET 54.735610317245360079 // asin(sr2/sr3)
#define AOCT 35.264389682754668343 // asin(sr1/sr3)
#define AICS 20.89774264557		     // asin(G/2/sr3)

void run() {
	TR << "START RUN!" << std::endl;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ScoreFunctionOP sf = core::scoring::getScoreFunction();
	ScoreFunctionOP sfrep = new core::scoring::ScoreFunction; sfrep->set_weight(core::scoring::fa_rep,1.0);
	ScoreFunctionOP sfdun = new core::scoring::ScoreFunction; sfdun->set_weight(core::scoring::fa_dun,1.0);
	ScoreFunctionOP sfsym = new core::scoring::symmetry::SymmetricScoreFunction(sf);
	sfsym->set_weight(core::scoring::atom_pair_constraint,1.0);
	ScoreFunctionOP sfrepsym = new core::scoring::symmetry::SymmetricScoreFunction(sfrep);


	core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::pack::dunbrack::SingleResidueRotamerLibraryCAP dunlib1 = core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( (rs->name_map("TYR")) );
	core::pack::dunbrack::RotamerLibraryScratchSpace scratch;

	Pose bpy,ala,tyr;
	core::import_pose::pose_from_pdb(bpy ,*rs,"input/bpy_ideal.pdb");
	core::pose::remove_lower_terminus_type_from_pose_residue(bpy,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(bpy,1);
	int const iCZ = bpy.residue(1).atom_index("CZ");
	int const iCP = bpy.residue(1).atom_index("CP");
	int const iCM = bpy.residue(1).atom_index("CM");
	int const iCB = bpy.residue(1).atom_index("CB");
	int const iCG = bpy.residue(1).atom_index("CG");
	int const iZN = bpy.residue(1).atom_index("ZN");
	int const iNN = bpy.residue(1).atom_index("NN1");
	int const iNE = bpy.residue(1).atom_index("NE1");

	make_pose_from_sequence(ala,"A","fa_standard",false);
	make_pose_from_sequence(tyr,"Y","fa_standard",false);
	remove_lower_terminus_type_from_pose_residue(ala,1);
	remove_upper_terminus_type_from_pose_residue(ala,1);

	Vec a2f1 = Vec(0,0,1);
	Vec c2f1 = Vec(0,0,0);
	Mat Rc2 = rotation_matrix_degrees(a2f1,180.0);
	Real bpy_clash_d2 = option[bpytoi::bpy_clash_dis]() * option[bpytoi::bpy_clash_dis]();

	vector1<string> const & fnames(option[OptionKeys::in::file::s]());
	TR << "gensym_3bpy_from_dimer " << fnames.size() << " files" << std::endl;

	for(Size ifile = 1; ifile <= fnames.size(); ++ifile) {
		string fname = fnames[ifile], infile = utility::file_basename(fnames[ifile]);
		cout << "PROGRESS "<< ifile << " " << fname << endl;
		Pose nat; core::import_pose::pose_from_pdb(nat,*rs,fname);
		core::scoring::dssp::Dssp dssp(nat);
		dssp.insert_ss_into_pose(nat);
		if( nat.n_residue() > option[bpytoi::max_nres]() ) { /*SIZE*/ continue; };
		Pose base(nat), pdimer(nat);
		make_dimer(pdimer);
		Size cyscnt = 0; {
			for(Size ir = 1; ir <= base.n_residue(); ++ir) if(base.residue(ir).name3()=="CYS") cyscnt++;
			if(cyscnt > option[bpytoi::max_cys]()) { /*CYS*/ continue; };
		}

		//TR << "gensym_3bpy_from_dimer " << ifile << " " << fnames[ifile] << " " << base.n_residue() << " residues" << " " << cyscnt << std::endl;
		for(Size ibpy = 1; ibpy <= base.n_residue(); ++ibpy) {
			if(!base.residue(ibpy).is_protein()) continue;
			// if(is_near_C2Z_iface(base,ibpy)) continue;
			if(base.residue(ibpy).is_lower_terminus()) continue;
			if(base.residue(ibpy).is_upper_terminus()) continue;
			// if(base.secstruct(ibpy)=='L') continue;
			// if(neighbor_count(base,ibpy) < option[bpytoi::min_bpy_cb_mono_nbr_count]()) continue;

			//cout << fname << " " << ibpy << endl;
			base = nat;
			base.replace_residue(ibpy,bpy.residue(1),true);
			Real chi1_incr = option[bpytoi::chi1_increment]();
			for(Real bch1 = 0.0; bch1 < 360.0; bch1 += chi1_incr) {
				base.set_chi(1,ibpy,bch1);
				// clash check atoms along chi2 vector
				for(Size ir = 1; ir <= base.n_residue(); ++ir) {
					if(ir==ibpy) continue;
					for(Size ia = 1; ia <= (base.residue(ir).name3()=="GLY"?4:5); ia++) {
						if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz(iCZ)) < 6.0 ) goto clash2;
						if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz(iCP)) < 6.0 ) goto clash2;
						if( base.xyz(AtomID(ia,ir)).distance_squared(base.residue(ibpy).xyz(iCM)) < 6.0 ) goto clash2;
					}
				}
				goto noclash2; clash2: { /*CLASH2*/ continue; }; noclash2:
				Vec const CB = base.residue(ibpy).xyz(iCB);
				Vec const CG = base.residue(ibpy).xyz(iCG);
				vector1<std::pair<Vec,Vec> > baxes = intersecting_bpy_axes(CB,CG,base.residue(ibpy).xyz("ZN"),a2f1,c2f1);
				for(Size jaxs = 1; jaxs <= baxes.size(); ++jaxs) {
					Vec const isct = baxes[jaxs].first;
					Vec const c3f1 = baxes[jaxs].second;
					// if( isct.distance_squared(c3f1) < 15.0 ) { /*ISCTDIS*/ continue; };
					// if( isct.distance_squared(c2f1) < 15.0 ) { /*ISCTDIS*/ continue; };
					Vec const a3f1 = (isct-c3f1).normalized();
					Mat const R1 = rotation_matrix_degrees(a3f1,120.0);
					Mat const R2 = rotation_matrix_degrees(a3f1,240.0);
					Real const orig_ang = angle_degrees( c2f1, isct, c3f1);
					Real const ang = (orig_ang > 90.0) ? 180.0-orig_ang : orig_ang;
					if( fabs(ang-ATET)*(isct.distance(c3f1))*0.01745418/3.0 > option[bpytoi::max_sym_error]() && // sin(x) ~ x for small x
					    fabs(ang-AOCT)*(isct.distance(c3f1))*0.01745418/3.0 > option[bpytoi::max_sym_error]() && // sin(x) ~ x for small x
					    fabs(ang-AICS)*(isct.distance(c3f1))*0.01745418/3.0 > option[bpytoi::max_sym_error]() ){ /*SYM_GEOM*/ continue; }; // sin(x) ~ x for small x

					Real const dang = dihedral_degrees( c3f1,CG,CB,base.residue(ibpy).xyz(iZN) );
					base.set_chi(2,ibpy, base.chi(2,ibpy) + dang );

					base.replace_residue(ibpy,tyr.residue(1),true);
					base.set_chi(1,ibpy, bch1 );
					base.set_chi(2,ibpy, base.chi(2,ibpy) + dang );
					Real bpy_dun = dunlib1->rotamer_energy( base.residue(ibpy), scratch );
					base.replace_residue(ibpy,bpy.residue(1),true);
					base.set_chi(1,ibpy, bch1 );
					base.set_chi(2,ibpy, base.chi(2,ibpy) + dang );
					if(bpy_dun > option[bpytoi::max_bpy_dun]()) { /*DUN*/ continue; };

					// clash check BPY vs trimer BB
					int bpy_mono_bb_atom_nbrs = 0;
					int bpy_tri_bb_atom_nbrs = 0;
					for(Size ir = 1; ir <= base.n_residue(); ++ir) {
						if(ir==ibpy) continue;
						for(Size ia = 1; ia <= min(5ul,base.residue(ir).nheavyatoms()); ia++) {
					// foreach bb atoms besides bpy:
							Vec const x0 =     base.xyz(AtomID(ia,ir))           ;
							Vec const x1 = R1*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
							Vec const x2 = R2*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
							for(Size ja = 7; ja <= base.residue(ibpy).nheavyatoms(); ja++) {
								if( x0.distance_squared(base.xyz(AtomID(ja,ibpy))) < bpy_clash_d2 ) goto clash3;
								if( x1.distance_squared(base.xyz(AtomID(ja,ibpy))) < bpy_clash_d2 ) goto clash3;
								if( x2.distance_squared(base.xyz(AtomID(ja,ibpy))) < bpy_clash_d2 ) goto clash3;
								if(ia == 5) {
									if( x0.distance_squared(base.xyz(AtomID(ja,ibpy))) < 64.0 ) bpy_mono_bb_atom_nbrs++;
									if( x0.distance_squared(base.xyz(AtomID(ja,ibpy))) < 64.0 ) bpy_tri_bb_atom_nbrs++;
									if( x1.distance_squared(base.xyz(AtomID(ja,ibpy))) < 64.0 ) bpy_tri_bb_atom_nbrs++;
									if( x2.distance_squared(base.xyz(AtomID(ja,ibpy))) < 64.0 ) bpy_tri_bb_atom_nbrs++;
								}
							}
						}
					}
					goto noclash3;	clash3: { /*CLASH3*/ continue; }; noclash3:
					// if( bpy_mono_bb_atom_nbrs < option[min_bpy_cb_mono_nbr_count]() ) continue;
					// if( bpy_tri_bb_atom_nbrs  < option[min_bpy_cb_tri_nbr_count]()  ) continue;

					// clash check trimer BB
					for(Size ir = 1; ir <= base.n_residue(); ++ir) {
						for(Size ia = 1; ia <= min(5ul,base.residue(ir).nheavyatoms()); ia++) {
							if(ir==ibpy) if(ia==iNE||ia==iNN) continue;
					// foreach BB
							Vec const x1 = R1*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
							Vec const x2 = R2*(base.xyz(AtomID(ia,ir))-c3f1)+c3f1;
							for(Size jr = 1; jr <= base.n_residue(); ++jr) {
								for(Size ja = 1; ja <= min(5ul,base.residue(jr).nheavyatoms()); ja++) {
							// foreach BB
									if( x1.distance_squared(      base.xyz(AtomID(ja,jr)))             < bpy_clash_d2 ) goto clash4;
									if( x1.distance_squared( Rc2*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < bpy_clash_d2 ) goto clash4;
									if( x2.distance_squared( Rc2*(base.xyz(AtomID(ja,jr))-c2f1)+c2f1 ) < bpy_clash_d2 ) goto clash4;
								}
							}
						}
					}
					goto noclash4;	clash4: { /*CLASH4*/ continue; }; noclash4:

					Pose psym = base;
					trans_pose(psym,-isct);
					Pose psym_bare = psym;
					Mat Rsymm;
					string symtag;
					Vec f3=(c3f1-isct).normalized(),f2=(c2f1-isct).normalized(),t1=Vec(0,0,1);
					//Vec newf2,newf3;
					Size dimersub = 0;
					Real angerr = 0.0;
					{
						Real dtet = fabs(ang-ATET);
						Real doct = fabs(ang-AOCT);
						Real dics = fabs(ang-AICS);
						Real dmin = min(dtet,min(doct,dics));
						if( dmin==dtet ) {
							angerr = fabs(ang-ATET);
							symtag = "TET";
							option[OptionKeys::symmetry::symmetry_definition]("input/sym/tetra.sym");
							Rsymm = alignVectorSets(f3,f2,t1, (orig_ang<90.0?1.0:-1.0)*Vec(0.8164965743782284,0.0,0.5773502784520137) );
							//newf2 = (orig_ang<90.0?1.0:-1.0)*Vec(0.8164965743782284,0.0,0.5773502784520137);
							//newf3 = Vec(0,0,1);
							dimersub =	4;
						} else
						if( dmin==doct ) {
							angerr = fabs(ang-AOCT);
							symtag = "OCT";
							option[OptionKeys::symmetry::symmetry_definition]("input/sym/octa.sym");
							Rsymm = alignVectorSets(f3,f2,t1, (orig_ang<90.0?1.0:-1.0)*Vec(0.4082482904638630,0.4082482904638626,0.8164965809277260));
							//newf2 = (orig_ang<90.0?1.0:-1.0)*Vec(0.4082482904638630,0.4082482904638626,0.8164965809277260);
							//newf3 = Vec(0,0,1);
							dimersub =	4;
						} else
						if( dmin==dics ) {
							angerr = fabs(ang-AICS);
							symtag = "ICS";
							option[OptionKeys::symmetry::symmetry_definition]("input/sym/icosa.sym");
							Rsymm = alignVectorSets(f3,f2,t1, (orig_ang<90.0?1.0:-1.0)*rotation_matrix_degrees(Vec(0,0,1),120.0)*Vec(0.35670090519235864157,0.0,0.93421863834701557305));
							//newf2 = (orig_ang<90.0?1.0:-1.0)*rotation_matrix_degrees(Vec(0,0,1),120.0)*Vec(0.35670090519235864157,0.0,0.93421863834701557305);
							//newf3 = Vec(0,0,1);
							dimersub =	4;
						} else {
							utility_exit_with_message("closed symm not yet supported");
						}
					}
					rot_pose(psym,Rsymm);
					//alignaxis(psym,newf2,f2); // f2/3 now invalid!!!!
					//Vec iron = psym.residue(ibpy).xyz("ZN");
					//rot_pose( psym, newf2, -dihedral_degrees(Vec(0,0,1),Vec(0,0,0),newf2,iron) );

					string outfname = symtag+"_"+infile+"_B"+lead_zero_string_of(ibpy,3)+"-"+lead_zero_string_of((Size)(bch1),3)+"_"+lead_zero_string_of(jaxs,1)+".pdb";
					core::pose::symmetry::make_symmetric_pose(psym);

					//psym.dump_pdb("test.pdb");
					//utility_exit_with_message("arisetnario");

					Real bpymdis = psym.residue(ibpy).xyz("ZN").distance(psym.residue(ibpy+base.n_residue()).xyz("ZN"));
					if( bpymdis > option[bpytoi::max_sym_error]() ) { /*BPYGEOM*/ continue; };

					bool clash = false;
					for(Size ir = 4; ir <= base.n_residue()-3; ++ir) {
						Size natom =	(psym.residue(ir).name3()=="BPY") ? 17 : 5;
						if( psym.residue(ir).name3()=="GLY" ) natom = 4;
						if( psym.residue(ir).name3()=="PRO" ) natom = 7;
						for(Size ia = 1; ia <= natom; ++ia) {
							Vec ip = psym.xyz(AtomID(ia,ir));
							for(Size is = 2; is <= 12; ++is) {
								if(is==4) continue;
								for(Size jr = 4 + (is-1)*base.n_residue(); jr <= is*base.n_residue()-3; ++jr) {
									for(Size ja = 1; ja <= 4; ja++) {
										if( ip.distance_squared( psym.xyz(AtomID(ja,jr)) ) < bpy_clash_d2 ) clash = true;
									}
								}
							}
						}
					}
					for(Size ir = 3*base.n_residue()+4; ir <= 4*base.n_residue()-3; ++ir) {
						Size natom =	(psym.residue(ir).name3()=="BPY") ? 17 : 5;
						if( psym.residue(ir).name3()=="GLY" ) natom = 4;
						for(Size ia = 1; ia <= natom; ++ia) {
							Vec ip = psym.xyz(AtomID(ia,ir));
							for(Size is = 2; is <= 12; ++is) {
								if(is==4) continue;
								for(Size jr = 4 + (is-1)*base.n_residue(); jr <= is*base.n_residue()-3; ++jr) {
									for(Size ja = 1; ja <= 4; ja++) {
										if( ip.distance_squared( psym.xyz(AtomID(ja,jr)) ) < bpy_clash_d2 ) clash = true;
									}
								}
							}
						}
					}
					//TR << "!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n!\n" << std::endl;
					if(clash) { /*CLASH5*/ continue; };

					Size ncontact = num_trimer_contacts(psym,base.n_residue());
					Real drms = dimer_rms(psym,pdimer);
					if( drms > option[bpytoi::max_sym_error]() ) { /*DIMER_RMS*/ continue; };

					psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_base.pdb");

					// Real batr = psym.energies().residue_total_energies(ibpy)[core::scoring::fa_atr];
					fixbb_design(psym,ibpy,dimersub);
					int nmut = 0; for(Size i = 1; i <= base.n_residue(); ++i) if(base.residue(i).name3()!=                psym.residue(i).name3()) nmut++;
					int nala = 0; for(Size i = 1; i <= base.n_residue(); ++i) if(base.residue(i).name3()!="ALA" && "ALA"==psym.residue(i).name3()) nala++;
					std::cout << "dump des" << endl;
					psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_des.pdb");

					std::cout << "fixbb done" << endl;

					using namespace core::pack::task;
					using namespace core::pack::task::operation;
					TaskFactoryOP task_factory( new TaskFactory );
					protocols::toolbox::task_operations::JointSequenceOperationOP revert = new protocols::toolbox::task_operations::JointSequenceOperation;
					Pose basesym(base);
					core::pose::symmetry::make_symmetric_pose(basesym);
					revert->add_pose(psym);
					revert->add_pose(basesym);
				 	task_factory->push_back(revert);
					protocols::filters::FilterOP filter = new protocols::simple_filters::ScoreTypeFilter(sfsym,core::scoring::total_score,sfsym->score(psym));
					// filter->apply(psym);
					// double tmp1 = filter->report_sm(psym);
					// double tmp2 =  sfsym->score(psym);
					// std::cout << tmp1 << " " << tmp2 << std::endl;
					// utility_exit_with_message("aoristn");
					using namespace core::scoring::constraints;
					using core::id::AtomID;
					psym.add_constraint(new AtomPairConstraint(AtomID(iZN,ibpy),AtomID(iZN,ibpy+1*base.n_residue()),new HarmonicFunc(0,0.01)));
					psym.add_constraint(new AtomPairConstraint(AtomID(iZN,ibpy),AtomID(iZN,ibpy+2*base.n_residue()),new HarmonicFunc(0,0.01)));
					core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
					movemap->set_jump(true); movemap->set_bb(false); movemap->set_chi(true);
					protocols::moves::MoverOP relax_mover = new protocols::simple_moves::symmetry::SymMinMover( movemap, sfsym, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
					relax_mover->apply(psym);
					// psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_min.pdb");
					sfsym->show(psym);
					protocols::design_opt::GreedyOptMutationMover gomm;
					gomm.task_factory( task_factory );
					gomm.scorefxn( sfsym );
					gomm.filter( filter );
					gomm.relax_mover( relax_mover );
					gomm.apply(psym);
					psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_greedy.pdb");


					// Real rms = core::scoring::CA_rmsd(psym,psym_bare);
					Real sc,int_area;
					vector1<Size> intra_subs; intra_subs.push_back(1); intra_subs.push_back(4);
					sfsym->score(psym);
					new_sc(psym,intra_subs,int_area,sc);
					// if(int_area < 200.0 || sc < 0.4) { /*SC*/ continue; }

					std::cout << "calc buttressing" << endl;

					// BPY buttressing atom count
					int bpy_tri_atom_nbrs  = 0;
					int bpy_mono_atom_nbrs = 0;
					Mat Rz0 = rotation_matrix_degrees(Vec(0,0,1),  0.0);
					Mat Rz1 = rotation_matrix_degrees(Vec(0,0,1),120.0);
					Mat Rz2 = rotation_matrix_degrees(Vec(0,0,1),240.0);
					for(Size ir = 1; ir <= base.n_residue(); ++ir) {
						if(ir==ibpy) continue;
						if(!psym.residue(ir).is_protein()) continue;
						for(Size ia = 1; ia <= psym.residue(ir).nheavyatoms(); ia++) {
						    Vec const x0 = Rz0*psym.xyz(AtomID(ia,ir));
							Vec const x1 = Rz1*psym.xyz(AtomID(ia,ir));
							Vec const x2 = Rz2*psym.xyz(AtomID(ia,ir));
							for(Size ja = 7; ja <= psym.residue(ibpy).nheavyatoms(); ja++) {
							  Vec Xbpy = psym.xyz(AtomID(ja,ibpy));
								if( x0.distance_squared(Xbpy) < 25.0 ) bpy_mono_atom_nbrs++;
								if( x0.distance_squared(Xbpy) < 25.0 ) bpy_tri_atom_nbrs++;
								if( x1.distance_squared(Xbpy) < 25.0 ) bpy_tri_atom_nbrs++;
								if( x2.distance_squared(Xbpy) < 25.0 ) bpy_tri_atom_nbrs++;
							}
						}
					}

					// continue;




					// repack(psym,ibpy,dimersub);
					// //psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_RPK.pdb");
					// Real s0 = sfsym->score(psym);
					// for(Size ir = 1; ir <= base.n_residue(); ++ir) {
					// for(Size ia = 1; ia <= psym.residue_type(ir).natoms(); ++ia) {
					// 	core::id::AtomID const aid(core::id::AtomID(ia,ir));
					// 	if("ICS"==symtag) psym.set_xyz( aid, psym.xyz(aid) + 20*Vec(-0.178411, 0.309017, 0.934172 ) );
					// 	if("OCT"==symtag) psym.set_xyz( aid, psym.xyz(aid) + 20*Vec( 0.408248, 0.408248, 0.816497 ) );
					// 	if("TET"==symtag) psym.set_xyz( aid, psym.xyz(aid) + 20*Vec( 0.816497, 0.000000, 0.577350 ) );
					// }
					// }
					// repack(psym,ibpy,dimersub);
					// //psym.dump_pdb(option[OptionKeys::out::file::o]()+"/"+outfname+"_DDG_RPK.pdb");
					// Real s1 = sfsym->score(psym);
					//
					// TR << "HIT " << outfname << " " <<	ang << " " << ibpy << " " << bch1 << " " << jaxs << " ddG: " << s0-s1 << " sc: " << sc << " " << int_area <<std::endl;
					//
					// core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
					// ss_out->fill_struct(psym_bare,outfname);
					// ss_out->add_energy("rms",rms);
					// ss_out->add_energy("nres",base.n_residue());
					// ss_out->add_energy("batr",batr);
					// ss_out->add_energy("ddg",s0-s1);
					// ss_out->add_energy("sc",sc);
					// ss_out->add_energy("int_area",int_area);
					// sfd.write_silent_struct( *ss_out, option[OptionKeys::out::file::o]()+"/"+option[basic::options::OptionKeys::out::file::silent ]() );

					int nmut2 = 0; for(Size i = 1; i <= base.n_residue(); ++i) if(base.residue(i).name3()!=psym.residue(i).name3()) nmut2++;
					int nala2 = 0; for(Size i = 1; i <= base.n_residue(); ++i) if(base.residue(i).name3()!="ALA" && "ALA"==psym.residue(i).name3()) nala2++;
					sfsym->score(psym);
					new_sc(psym,intra_subs,int_area,sc);

					cout << "DEBUG " << (angerr*isct.distance(c3f1))*0.01745418/option[bpytoi::max_sym_error]() << std::endl;

					cout << "DOCK_HIT "
						 << option[OptionKeys::out::file::o]()+"/"+outfname << " "
						 << F(5,3,angerr) << " "
						 << F(5,3,drms) << " "
						 << F(5,3,bpymdis) << " "
						 << F(5,3,bpy_dun) << " "
						 << I(2,neighbor_count(base,ibpy)) << " "
						 << I(3,bpy_mono_bb_atom_nbrs) << " "
						 << I(3,bpy_tri_bb_atom_nbrs ) << " "
						 << I(3,bpy_mono_atom_nbrs ) << " "
						 << I(3,bpy_tri_atom_nbrs ) << " "
						 << F(7,1,int_area) << " "
						 << F(5,3,sc) << " "
						 << F(9,3,sfsym->score(psym)) << " "
						 << I(2,nmut) << " "
						 << I(2,nmut-nala) << " "
						 << I(2,nmut2) << " "
						 << I(2,nmut2-nala2) << " "
						 << endl;
				}
			}
		}
	}
	/*

		string fname = infile+"_"+lead_zero_string_of(irsd,3)+"_"+lead_zero_string_of((Size)(chi1+180.0),3)+"_"+lead_zero_string_of((Size)(chi2+180.0),3)+"_"+lead_zero_string_of(jrsd,3)+"_"+lead_zero_string_of((Size)(totrot),3)+"_"+lead_zero_string_of(iaxs,3)+".pdb";

		sf->score(psym);
		if(psym.energies().total_energies()[core::scoring::rg] > 14.0) {
		TR << "DANGER DANGER DANGER!!!" << std::endl;
		continue;
		}
		psym.dump_pdb(fname);

		core::pose::replace_pose_residue_copying_existing_coordinates(psym,jrsd,rs->name_map("ALA"));
		sf->score(psym);
		core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
		ss_out->fill_struct(psym,fname);
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

int main(int argc, char **argv){
	try {
		register_options();
		devel::init(argc,argv);
		run();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}


















