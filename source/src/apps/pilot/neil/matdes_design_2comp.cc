// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


#include <basic/database/open.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <core/scoring/packing/HolesParams.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <protocols/docking/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <basic/MetricValue.hh>

static thread_local basic::Tracer TR( "2comp_design" );

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;
using core::scoring::ScoreFunctionOP;
using core::pose::Pose;
using std::cout;
using std::endl;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef vector1<Size> Sizes;

OPT_1GRP_KEY( RealVector   , matdes2c, comp2rot )
OPT_1GRP_KEY( RealVector   , matdes2c, comp2disp )
OPT_1GRP_KEY( RealVector   , matdes2c, comp1rot )
OPT_1GRP_KEY( RealVector   , matdes2c, comp1disp )
OPT_1GRP_KEY( IntegerVector, matdes2c, intra_subs1 )
OPT_1GRP_KEY( IntegerVector, matdes2c, intra_subs2 )
OPT_1GRP_KEY( FileVector, matdes2c, I2 )
OPT_1GRP_KEY( FileVector, matdes2c, I3 )
OPT_1GRP_KEY( FileVector, matdes2c, I5 )
OPT_1GRP_KEY( FileVector, matdes2c, O2 )
OPT_1GRP_KEY( FileVector, matdes2c, O3 )
OPT_1GRP_KEY( FileVector, matdes2c, O4 )
OPT_1GRP_KEY( FileVector, matdes2c, T2 )
OPT_1GRP_KEY( FileVector, matdes2c, T3 )
OPT_1GRP_KEY( Boolean   , matdes2c, mutate_intra_bb )
OPT_1GRP_KEY( Boolean   , matdes2c, dump_symmetric )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( matdes2c::comp2rot    , "", -1 );
	NEW_OPT( matdes2c::comp2disp   , "", -1 );
	NEW_OPT( matdes2c::comp1rot    , "", -1 );
	NEW_OPT( matdes2c::comp1disp   , "", -1 );
	NEW_OPT( matdes2c::intra_subs1 , "", -1 );
	NEW_OPT( matdes2c::intra_subs2 , "", -1 );
	NEW_OPT( matdes2c::I2      , "file(s) for icos 2fold comp1"     , ""     );
	NEW_OPT( matdes2c::I3      , "file(s) for icos 3fold comp2"    , ""     );
	NEW_OPT( matdes2c::I5      , "file(s) for icos 5fold pentamer"  , ""     );
	NEW_OPT( matdes2c::O2      , "file(s) for octa 2fold comp1"     , ""     );
	NEW_OPT( matdes2c::O3      , "file(s) for octa 3fold comp2"    , ""     );
	NEW_OPT( matdes2c::O4      , "file(s) for octa 4fold tetramer"  , ""     );
	NEW_OPT( matdes2c::T2      , "file(s) for tetr 2fold comp1"     , ""     );
	NEW_OPT( matdes2c::T3      , "file(s) for tetr 3fold comp2"    , ""     );
	NEW_OPT( matdes2c::mutate_intra_bb      , ""    , false );
	NEW_OPT( matdes2c::dump_symmetric      , ""    , false );
}


core::kinematics::Stub getxform(core::conformation::Residue const & move_resi, core::conformation::Residue const & fixd_resi) {
	core::kinematics::Stub s;
	s.M = alignVectorSets(move_resi.xyz(1)-move_resi.xyz(2),move_resi.xyz(3)-move_resi.xyz(2),fixd_resi.xyz(1)-fixd_resi.xyz(2),fixd_resi.xyz(3)-fixd_resi.xyz(2));
	s.v = fixd_resi.xyz(2)-s.M*move_resi.xyz(2);
	return s;
}

void trans_pose( Pose & pose, Vec const & trans, Size start=1, Size end=0 ) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + (Vec)trans );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Size start=1, Size end=0 ) {
	if(0==end) end = pose.n_residue();
	for(Size ir = start; ir <= end; ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}
void rot_pose( Pose & pose, Mat const & rot, Vec const & cen, Size start=1, Size end=0 ) {
	trans_pose(pose,-cen,start,end);
	rot_pose(pose,rot,start,end);
	trans_pose(pose,cen,start,end);
}
void rot_pose( Pose & pose, Vec const & axis, double const & ang, Size start=1, Size end=0 ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),start,end);
}
void rot_pose( Pose & pose, Vec const & axis, double const & ang, Vec const & cen, Size start=1, Size end=0 ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen,start,end);
}
void alignaxis(core::pose::Pose & pose, Vec newaxis, Vec oldaxis, Vec cen = Vec(0,0,0) ) {
	newaxis.normalize();
	oldaxis.normalize();
	Vec axis = newaxis.cross(oldaxis).normalized();
	Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
	rot_pose(pose,axis,ang,cen);
}


inline Real sqr(Real x) { return x*x; }

void print_movemap(core::kinematics::MoveMap const & movemap) {
	using namespace core::id;
	using namespace core::kinematics;
	TR << "movemap " << std::endl;
	for(std::map< TorsionType, bool >::const_iterator i = movemap.torsion_type_begin(); i != movemap.torsion_type_end(); ++i) {
		TR << "TorsionType " << i->first << " " << i->second << std::endl;
	}
	for(std::map< std::pair< Size, TorsionType >, bool >::const_iterator i = movemap.movemap_torsion_id_begin(); i != movemap.movemap_torsion_id_end(); ++i) {
		TR << "MoveMapTorsionID (" << i->first.first << "," << i->first.second << ") " << i->second << std::endl;
	}
	for(std::map< TorsionID, bool >::const_iterator i = movemap.torsion_id_begin(); i != movemap.torsion_id_end(); ++i) {
		TR << "TorsionID " << i->first << " " << i->second << std::endl;
	}
	for(std::map< DOF_Type, bool >::const_iterator i = movemap.dof_type_begin(); i != movemap.dof_type_end(); ++i) {
		TR << "DOF_Type " << i->first << " " << i->second << std::endl;
	}
	for(std::map< DOF_ID, bool >::const_iterator i = movemap.dof_id_begin(); i != movemap.dof_id_end(); ++i) {
		TR << "DOF_ID " << i->first << " " << i->second << std::endl;
	}
	for(std::map< JumpID, bool >::const_iterator i = movemap.jump_id_begin(); i != movemap.jump_id_end(); ++i) {
		TR << "JumpID " << i->first << " " << i->second << std::endl;
	}

}


void design(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool hphobic_only) {

	using namespace core;
	using namespace pack;
	using namespace task;
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace chemical;

	// Set allowed AAs.
	vector1<bool> allowed_aas(20, true);
//  vector1<bool> allowed_aas(20, false); // Basically dock, repack, and score the input sequence
	allowed_aas[aa_cys] = false;
	allowed_aas[aa_gly] = false;
	allowed_aas[aa_pro] = false;
	// Set used to design O333, for use with fa_elec9000

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	// Set which residues can be designed
	for(Size i=1; i<=pose.n_residue(); i++) {
		if(!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if(pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
			// Don't mess with Pros or Glys at the interfaces
			task->nonconst_residue_task(i).prevent_repacking();
		} else if(find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
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

	make_symmetric_PackerTask_by_truncation(pose, task);
	// Get rid of Nobu's rotamers of death
	core::pack::task::operation::TaskOperationCOP limit_rots = new protocols::toolbox::task_operations::LimitAromaChi2Operation();
	limit_rots->apply(pose, *task);
	// Actually perform design.
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

}


void
design_using_resfile(Pose & pose, ScoreFunctionOP sf, std::string resfile, utility::vector1<Size> & design_pos) {

  using namespace core;
  using namespace pack;
  using namespace task;
  using namespace conformation;
  using namespace conformation::symmetry;
  using namespace scoring;
  using namespace chemical;

  // Get the symmetry info and make the packer task
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
  PackerTaskOP task(TaskFactory::create_packer_task(pose));

	// Modify the packer task according to the resfile
	core::pack::task::parse_resfile(pose,*task, resfile);


	// If design_pos is populated then we are doing design
	if(design_pos.size() > 0){

		for(Size i=1; i<=pose.n_residue(); i++) {
			if(!sym_info->bb_is_independent(i)) {
				task->nonconst_residue_task(i).prevent_repacking();
			} else if(pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
				// Don't mess with Pros or Glys at the interfaces
				task->nonconst_residue_task(i).prevent_repacking();
			} else if(find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
				task->nonconst_residue_task(i).prevent_repacking();
			} else {
				task->nonconst_residue_task(i).initialize_from_command_line();
			}
		}


	} else {


	// Get from the packer task the mutalyze positions
		Size nres_monomer = sym_info->num_independent_residues();
		for(Size i=1; i<=nres_monomer; ++i) {
			if(! task->residue_task(i).has_behavior("AUTO")) {
				design_pos.push_back(i);
				task->nonconst_residue_task(i).initialize_from_command_line();
			}
		}
	}
  make_symmetric_PackerTask_by_truncation(pose, task);
	// Get rid of Nobu's rotamers of death
	core::pack::task::operation::TaskOperationCOP limit_rots = new protocols::toolbox::task_operations::LimitAromaChi2Operation();
	limit_rots->apply(pose, *task);
  // Actually perform design
  protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);

}


void repack(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos) {

	using namespace core;
	using namespace pack;
	using namespace task;
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace chemical;

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	// Set which residues can be repacked
	for(Size i=1; i<=pose.n_residue(); i++) {
		if(!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if(find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else {
			vector1<bool> allowed_aas(20, false);
			allowed_aas[pose.residue(i).aa()] = true;
			task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
			task->nonconst_residue_task(i).initialize_from_command_line();
		}
	}

	// Actually repack.
	make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

}

void minimize(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool move_bb, bool move_sc, bool move_rb) {

	// Initialize a MoveMap
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(move_rb);
	movemap->set_bb(false);
	movemap->set_chi(false);

	// Set allowable move types at interface positions
	// Currently, only sc moves allowed
	for(utility::vector1<Size>::iterator i = design_pos.begin(); i != design_pos.end(); i++) {
		movemap->set_bb (*i, move_bb);
		movemap->set_chi(*i, move_sc);
	}

	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

int which_subsub(int i,Pose const & p1, Pose const & p2) {
	int N = p1.n_residue()+p2.n_residue();
	i = (i-1)%N+1;
	return i > p1.n_residue() ? 2 : 1;
}

utility::vector1<Real> sidechain_sasa(Pose const & pose, Real probe_radius) {
	using core::id::AtomID;
	utility::vector1<Real> rsd_sasa(pose.n_residue(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose,false);
	for(Size i = 1; i <= pose.n_residue(); i++) {
		for(Size j = 1; j <= pose.residue(i).nheavyatoms(); j++) {
			atom_mask[AtomID(j,i)] = true;
		}
	}
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );
	utility::vector1<Real> sc_sasa(pose.n_residue(),0.0);
	for(Size i = 1; i <= pose.n_residue(); i++) {
		// Use CA as the side chain for Glys
		if(pose.residue(i).name3()=="GLY") sc_sasa[i] += atom_sasa[AtomID(2,i)];
		for(Size j = 5; j <= pose.residue(i).nheavyatoms(); j++) {
			sc_sasa[i] += atom_sasa[AtomID(j,i)];
		}
	}
	return sc_sasa;
}

void new_sc(Pose &pose, Sizes isubs1, Sizes isubs2, Pose const & p1, Pose const & p2, Real& int_area, Real& sc) {
	using namespace core;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	core::scoring::sc::ShapeComplementarityCalculator scc; scc.Init();
	Size nres_monomer = symm_info->num_independent_residues();
	std::set<Size> iset,jset;
	for(Size ir=1; ir <= symm_info->num_total_residues_without_pseudo(); ++ir){
		if(std::find(isubs1.begin(),isubs1.end(),symm_info->subunit_index(ir))==isubs1.end()) continue; // sub must be 123
		if(which_subsub(ir,p1,p2)!=1) continue;                                                         // subsub must be 1
		for(Size jr=1; jr <= symm_info->num_total_residues_without_pseudo(); ++jr){
			if(std::find(isubs2.begin(),isubs2.end(),symm_info->subunit_index(jr))==isubs2.end()) continue; // sub 14
			if(which_subsub(jr,p1,p2)!=2) continue;                                                         // subsub 2
			//if(pose.residue(ir).nbr_atom_xyz().distance_squared(pose.residue(jr).nbr_atom_xyz()) >
			// sqr(pose.residue(ir).nbr_radius()+pose.residue(jr).nbr_radius()+4.0)) continue;
			iset.insert(ir);
			jset.insert(jr);
		}
	}
	TR << "SC res sets: " << iset.size() << " " << jset.size() << std::endl;
	for(std::set<Size>::const_iterator i=iset.begin(); i!=iset.end(); ++i) { cout << "+" << *i; scc.AddResidue(0,pose.residue(*i)); } cout << std::endl;
	for(std::set<Size>::const_iterator i=jset.begin(); i!=jset.end(); ++i) { cout << "+" << *i; scc.AddResidue(1,pose.residue(*i)); } cout << std::endl;
	if(scc.Calc()) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	}						TR << std::endl;
	TR << "SC DONE" << std::endl;
	// pose.dump_pdb("test.pdb");
	// utility_exit_with_message("oairsnt");
}

// Pose must be scored in order for this to work.
Pose get_neighbor_subs(Pose const &pose, Sizes intra_subs1, Sizes intra_subs2, Pose const & p1, Pose const & p2){
	// Figure out which chains touch chain A, and return those chains
	Pose sub_pose;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	Size nres_monomer = symm_info->num_independent_residues();
	for(Size i=1; i<=nres_monomer; ++i) {
		if(pose.residue(i).is_lower_terminus()) sub_pose.append_residue_by_jump(pose.residue(i),1);
		else                                    sub_pose.append_residue_by_bond(pose.residue(i));
	}
	for(Size i=1; i<=symm_info->subunits(); ++i) {
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for(Size ir=1; ir<=nres_monomer; ir++) {
			Sizes const & intra_subs(which_subsub(ir,p1,p2)==1?intra_subs1:intra_subs2);
			if(std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end()) continue;
			if(pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0) {
				contact = true;
				break;
			}
		}
		if(contact) {
			for(Size i=1; i<=nres_monomer; ++i) {
				if(pose.residue(i).is_lower_terminus()) sub_pose.append_residue_by_jump(pose.residue(start+i),1);
				else                                    sub_pose.append_residue_by_bond(pose.residue(start+i));
			}
		}
	}
	return sub_pose;
}


// Find residues with buried polar atoms.
Real
get_unsat_polars( Pose const &bound, Pose const &unbound, Size nres_monomer, string fn){

  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_bound = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
  basic::MetricValue< id::AtomID_Map<bool> > bound_Amap;
  unsat_calc_bound->get("atom_bur_unsat", bound_Amap, bound);

  core::pose::metrics::PoseMetricCalculatorOP unsat_calc_unbound = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
  basic::MetricValue< id::AtomID_Map<bool> > unbound_Amap;
  unsat_calc_unbound->get("atom_bur_unsat", unbound_Amap, unbound);

	id::AtomID_Map<bool> bound_am = bound_Amap.value();
	id::AtomID_Map<bool> unbound_am = unbound_Amap.value();
	Size buried_unsat_polars = 0;
	string select_buried_unsat_polars("select buried_unsat_polars, (");
	for (Size ir=1; ir<=nres_monomer; ir++) {
		Size flag = 0;
		for (Size ia=1; ia<=bound.residue(ir).nheavyatoms(); ia++) {
			if (bound_am[id::AtomID(ia,ir)] != unbound_am[id::AtomID(ia,ir)]) {
				buried_unsat_polars++;
				if (flag == 0) {
					TR << "buried unsat polar(s): " << bound.residue(ir).name3() << ir << "\t" << bound.residue(ir).atom_name(ia);
					select_buried_unsat_polars.append("resi " + ObjexxFCL::string_of(ir) + " and name " + bound.residue(ir).atom_name(ia) + "+ ");
					flag = 1;
				} else {
					TR << "," << bound.residue(ir).atom_name(ia);
				}
			}
		}
		if (flag) TR << std::endl;
	}
	select_buried_unsat_polars.erase(select_buried_unsat_polars.end()-2,select_buried_unsat_polars.end());
	TR << select_buried_unsat_polars << ") and " << fn << "_0001";
	TR << std::endl;
	return buried_unsat_polars;

}


Real get_atom_packing_score(Pose const &pose, Sizes intra_subs1, Sizes intra_subs2, Pose const & p1, Pose const & p2, Real cutoff=9.0){
	Pose sub_pose = get_neighbor_subs(pose, intra_subs1,intra_subs2,p1,p2);
	core::scoring::packing::HolesParams hp(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	core::scoring::packing::HolesResult hr(core::scoring::packing::compute_holes_score(sub_pose, hp));
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	Size nres_monomer = symm_info->num_independent_residues();
	Size count = 0; Real if_score = 0;

	Real cutoff2 = cutoff*cutoff;
	for(Size ir=1; ir<=nres_monomer; ir++) {
		for(Size ia = 1; ia<=sub_pose.residue(ir).nheavyatoms(); ia++) {
			bool contact = false;
			for(Size jr=nres_monomer+1; jr<=sub_pose.n_residue(); jr++) {
				for(Size ja = 1; ja<=sub_pose.residue(jr).nheavyatoms(); ja++) {
					if(sub_pose.residue(ir).xyz(ia).distance_squared(sub_pose.residue(jr).xyz(ja)) <= cutoff2)  {
						contact = true;
						break; // ja
					}
				} // ja
				if(contact == true) break;
			} // jr
			if(contact == true) {
				count++;
				if_score += hr.atom_scores[AtomID(ia, ir)];
			}
		} // ia
	} // ir
	return if_score /(Real)count;
}


Real average_degree(Pose const &pose, vector1<Size> mutalyze_pos, Sizes intra_subs1, Sizes intra_subs2,
	                  Pose const & p1, Pose const & p2, Real distance_threshold=10.0){
	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	Size nres_monomer = sym_info->num_independent_residues();
	Size count_neighbors=0;
	for(Size i=1; i <= mutalyze_pos.size(); ++i) {
		Size ires=mutalyze_pos[i];
		core::conformation::Residue const resi( pose.conformation().residue( ires ) );
		Size resi_neighbors=0;
		Sizes const & intra_subs(which_subsub(ires,p1,p2)==1?intra_subs1:intra_subs2);
		for(Size jres = 1; jres <= sym_info->num_total_residues_without_pseudo(); ++jres) {
			if(std::find(intra_subs.begin(),intra_subs.end(),sym_info->subunit_index(jres))==intra_subs.end()) continue;
			core::conformation::Residue const resj( pose.residue( jres ) );
			Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if( distance <= distance_threshold ){	++count_neighbors; ++resi_neighbors; }
		}
		TR << "avg_deg of " << resi.name3() << ires << " = " << resi_neighbors << std::endl;
	}
	return( (Real) count_neighbors / mutalyze_pos.size() );
}

void *dostuff(void*) {
	using namespace core;
	using namespace basic;
	using namespace options;
	using namespace OptionKeys;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace scoring;
	using namespace utility;
	using basic::options::option;

	chemical::ResidueTypeSetCAP resi_set = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::io::silent::SilentFileData sfd;

	vector1<string> compkind;
	vector1<vector1<string> > compfiles;
	if( option[matdes2c::I5].user() ) { compkind.push_back("I5"); compfiles.push_back(option[matdes2c::I5]()); }
	if( option[matdes2c::I3].user() ) { compkind.push_back("I3"); compfiles.push_back(option[matdes2c::I3]()); }
	if( option[matdes2c::I2].user() ) { compkind.push_back("I2"); compfiles.push_back(option[matdes2c::I2]()); }
	if( option[matdes2c::O4].user() ) { compkind.push_back("O4"); compfiles.push_back(option[matdes2c::O4]()); }
	if( option[matdes2c::O3].user() ) { compkind.push_back("O3"); compfiles.push_back(option[matdes2c::O3]()); }
	if( option[matdes2c::O2].user() ) { compkind.push_back("O2"); compfiles.push_back(option[matdes2c::O2]()); }
	if( option[matdes2c::T3].user() ) { compkind.push_back("T3"); compfiles.push_back(option[matdes2c::T3]()); }
	if( option[matdes2c::T2].user() ) { compkind.push_back("T2"); compfiles.push_back(option[matdes2c::T2]()); }
	if(compkind.size()!=2) utility_exit_with_message("must specify two components!");
	if(compkind[1][0]!=compkind[2][0]) utility_exit_with_message("components must be of same sym icos, tetra or octa");
	string cmp1type = compkind[1];
	string cmp2type = compkind[2];
	std::map<string,int> comp_nangle;
	comp_nangle["I2"] = 180; comp_nangle["O2"] = 180; comp_nangle["T2"] = 180;
	comp_nangle["I3"] = 120; comp_nangle["O3"] = 120; comp_nangle["T3"] = 120;
	comp_nangle["I5"] =  72; comp_nangle["O4"] =  90;
	std::map<string,Vec> axismap;
	axismap["I2"] = Vec(-0.166666666666666,-0.28867500000000000, 0.87267800000000 ).normalized();
	axismap["I3"] = Vec( 0.000000000000000, 0.00000000000000000, 1.00000000000000 ).normalized();
	axismap["I5"] = Vec(-0.607226000000000, 0.00000000000000000, 0.79452900000000 ).normalized();
	axismap["O2"] = Vec( 1.000000000000000, 1.00000000000000000, 0.00000000000000 ).normalized();
	axismap["O3"] = Vec( 1.000000000000000, 1.00000000000000000, 1.00000000000000 ).normalized();
	axismap["O4"] = Vec( 1.000000000000000, 0.00000000000000000, 0.00000000000000 ).normalized();
	axismap["T2"] = Vec( 0.816496579408716, 0.00000000000000000, 0.57735027133783 ).normalized();
	axismap["T3"] = Vec( 0.000000000000000, 0.00000000000000000, 1.00000000000000 ).normalized();
	Vec cmp1axs = axismap[cmp1type];
	Vec cmp2axs = axismap[cmp2type];
	Real cmp1nangle = comp_nangle[cmp1type];
	Real cmp2nangle = comp_nangle[cmp2type];

	option[OptionKeys::symmetry::symmetry_definition]("input/"+compkind[1].substr(0,1)+".sym");
	// Create a score function object, turn fa_elec off in the monomer
	ScoreFunctionOP sf = get_score_function();
	// core::scoring::methods::EnergyMethodOptions eo = sf->energy_method_options();
	// eo.exclude_monomer_fa_elec(true);
	// sf->set_energy_method_options(eo);

	Real const contact_dist = option[matdes::design::contact_dist]();
	Real const contact_dist_sq = contact_dist * contact_dist;

	// Define which subs are part of the oligomeric building block. The interfaces between these
	// subunits should not be designed.
	Sizes intra_subs1 = option[matdes2c::intra_subs1]();
	Sizes intra_subs2 = option[matdes2c::intra_subs2]();
	// std::cout << intra_subs1.size() << " " << intra_subs2.size() << std::endl;
	// for(int i = 1; i <= intra_subs1.size(); ++i) { std::cout << "intra_subs1 " << intra_subs1[i] << std::endl; }
	// for(int i = 1; i <= intra_subs2.size(); ++i) { std::cout << "intra_subs2 " << intra_subs2[i] << std::endl; }

	utility::vector1<std::string> files1 = compfiles[1];
	utility::vector1<std::string> files2 = compfiles[2];
	utility::vector1<Real> cmp2rots = option[matdes2c::comp2rot]();
	utility::vector1<Real> cmp1rots = option[matdes2c::comp1rot]();
	utility::vector1<Real> cmp2disps = option[matdes2c::comp2disp]();
	utility::vector1<Real> cmp1disps = option[matdes2c::comp1disp]();
	if(cmp2rots.size()!=cmp1rots.size()||cmp2rots.size()!=cmp2rots.size()) utility_exit_with_message("comp2rot, comp1rot & orirots must be same size");
	if(cmp2rots.size()!=cmp1disps.size()||cmp2rots.size()!=cmp2disps.size()) utility_exit_with_message("comp2rot, comp1rot & orirots must be same size");
	if(files1.size() > 1 && files1.size() != cmp2rots.size() ) utility_exit_with_message("num files1 must be 1 or same num as rots");
	if(files2.size() > 1 && files2.size() != cmp2rots.size() ) utility_exit_with_message("num files2 must be 1 or same num as rots");

	for(Size iconfig = 1; iconfig <= cmp2rots.size(); ++iconfig) {
		std::string file1 = files1[files1.size()==1?1:iconfig];
		std::string file2 = files2[files2.size()==1?1:iconfig];

			// Read in pose
		Pose pose,mono,p1,p2;
		utility::vector1<Real> sc_sasa1,sc_sasa2,sc_sasa;
		{
			import_pose::pose_from_pdb(p1, file1, resi_set);
			import_pose::pose_from_pdb(p2, file2, resi_set);
			sc_sasa1 = sidechain_sasa(p1,2.2); // 110728 -- WAS 2.5 -- Changed for xtal design
			sc_sasa2 = sidechain_sasa(p2,2.2);
			if(sc_sasa1.size() != p1.n_residue()) utility_exit_with_message("SASA BAD");
			if(sc_sasa2.size() != p2.n_residue()) utility_exit_with_message("SASA BAD");
			sc_sasa = sc_sasa1;
			sc_sasa.insert(sc_sasa.end(),sc_sasa2.begin(),sc_sasa2.end());
			rot_pose(p1,Vec(0,0,1),cmp1rots[iconfig]);
			rot_pose(p2,Vec(0,0,1),cmp2rots[iconfig]);
			trans_pose(p1,Vec(0,0,cmp1disps[iconfig]));
			trans_pose(p2,Vec(0,0,cmp2disps[iconfig]));
			alignaxis(p1,cmp1axs,Vec(0,0,1),Vec(0,0,0));
			alignaxis(p2,cmp2axs,Vec(0,0,1),Vec(0,0,0));
			Pose mono;
			for(Size i = 1; i <= p1.n_residue(); ++i) {
				if(p1.residue(i).is_lower_terminus()) mono.append_residue_by_jump(p1.residue(i),1);
				else                                  mono.append_residue_by_bond(p1.residue(i));
			}
			for(Size i = 1; i <= p2.n_residue(); ++i) {
				if(p2.residue(i).is_lower_terminus()) mono.append_residue_by_jump(p2.residue(i),1);
				else                                  mono.append_residue_by_bond(p2.residue(i));
			}
			if(sc_sasa.size() != mono.n_residue()) utility_exit_with_message("SASA BAD");
			pose = mono;
			// mono.dump_pdb("test.pdb");
			// for(int i = 1; i <= sc_sasa.size(); ++i) {
			// 	cout << i << " " << sc_sasa[i] << endl;
			// }
			// utility_exit_with_message("DBG SASA");
		}

		// Handle all of the symmetry stuff
		core::pose::symmetry::make_symmetric_pose(pose);
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		Pose original_pose = pose;

		// std::map<Size,SymDof> dofs = sym_info->get_dofs();
		// int sym_jump = 0;
		// for(std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++) {
		// 	Size jump_num = i->first;
		// 	if(sym_jump == 0) {
		// 		sym_jump = jump_num;
		// 	} else {
		// 		utility_exit_with_message("Can only handle one subunit!");
		// 	}
		// }
		// if(sym_jump == 0) {
		// 	utility_exit_with_message("No jump defined!");
		// }


		// Set up the finer rigid body grid that we sample around each starting configuration
		Real grid_size_angle   = option[matdes::design::grid_size_angle  ]();
		Real grid_size_radius  = option[matdes::design::grid_size_radius ]();
		Size grid_nsamp_angle  = option[matdes::design::grid_nsamp_angle ]();
		Size grid_nsamp_radius = option[matdes::design::grid_nsamp_radius]();
		Real grid_start_angle  = -grid_size_angle / 2.0;
		Real grid_start_radius = 0.0;
		Real grid_incr_angle   = grid_size_angle  / (grid_nsamp_angle -1);
		Real grid_incr_radius  = grid_size_radius / (grid_nsamp_radius-1);
		if( grid_nsamp_angle  == 1 ) { grid_start_angle = 0.0; grid_incr_angle = 0.0; }
		if( grid_nsamp_radius == 1 ) { grid_start_radius = 0.0; grid_incr_radius = 0.0; }

		// Get the favor_native_residue weight from the command line
		Real fav_nat_bonus = 0.0;
		if(option[matdes::design::fav_nat_bonus]()) {
			fav_nat_bonus = option[matdes::design::fav_nat_bonus]();
		}

		for(Size iangle1 = 1; iangle1 <= grid_nsamp_angle; iangle1++) {
			Mat rot1 = numeric::rotation_matrix_degrees(cmp1axs,(grid_start_angle + (iangle1-1)*grid_incr_angle));
			for(Size iradius1 = 1; iradius1 <= grid_nsamp_radius; iradius1++) {
				Vec trans1 = (grid_start_radius + (iradius1-1)*grid_incr_radius) * cmp1axs;
				for(Size iangle2 = 1; iangle2 <= grid_nsamp_angle; iangle2++) {
					Mat rot2 = numeric::rotation_matrix_degrees(cmp2axs,(grid_start_angle + (iangle2-1)*grid_incr_angle));
					for(Size iradius2 = 1; iradius2 <= grid_nsamp_radius; iradius2++) {
						Vec trans2 = (grid_start_radius + (iradius2-1)*grid_incr_radius) * cmp2axs;

						std::cout << iconfig << " " << iangle1 << " " << iradius1 << " " << iangle1 << " " << iradius2 << std::endl;
						Pose pose_for_design = pose;
						rot_pose  (pose_for_design,  rot1,               1,p1.n_residue());
						trans_pose(pose_for_design,trans1,               1,p1.n_residue());
						rot_pose  (pose_for_design,  rot2,p1.n_residue()+1,p1.n_residue()+p2.n_residue());
						trans_pose(pose_for_design,trans2,p1.n_residue()+1,p1.n_residue()+p2.n_residue());

						// pose_for_design.dump_pdb("test_grid_"+ObjexxFCL::string_of(iconfig)+"_"+ObjexxFCL::string_of(iangle1)+"_"+ObjexxFCL::string_of(iradius1)+"_"+ObjexxFCL::string_of(iangle2)+"_"+ObjexxFCL::string_of(iradius2)+".pdb");
						// utility_exit_with_message("aorisht");


						// Find out which positions are near the inter-subunit interfaces
						// These will be further screened below, then passed to design()
						vector1<bool> indy_resis = sym_info->independent_residues();
						Sizes interface_pos;
						for(Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
							if(sym_info->subunit_index(ir) != 1) continue;
							std::string atom_i = (pose_for_design.residue(ir).name3() == "GLY") ? "CA" : "CB";
							for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
								Sizes const & isubs( which_subsub(ir,p1,p2)==1?intra_subs1:intra_subs2);
								if(which_subsub(jr,p1,p2)==which_subsub(ir,p1,p2)&&find(isubs.begin(),isubs.end(),sym_info->subunit_index(jr))!=isubs.end()) continue;
								// if both dimer or trimer, and jr in primary sub
								// std::cout << sym_info->subunit_index(ir) << " " << sym_info->subunit_index(jr) << std::endl;
								std::string atom_j = (pose_for_design.residue(jr).name3() == "GLY") ? "CA" : "CB";
								if(pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
									interface_pos.push_back(ir);
									break;
								}
							}
						}

						// Here we filter the residues that we are selecting for design
						// to get rid of those that make intra-building block interactions
						Sizes nontrimer_pos;
						Real all_atom_contact_thresh2 = 5.0*5.0;
						sf->score(pose_for_design);
						for(Size iip=1; iip<=interface_pos.size(); iip++) {
							Size ir = interface_pos[iip];
							Sizes const & intra_subs(which_subsub(ir,p1,p2)==1?intra_subs1:intra_subs2);
							bool contact = false;
							for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
								if(which_subsub(ir,p1,p2)!=which_subsub(jr,p1,p2)) continue;
								if(find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) == intra_subs.end()) continue;
								if(sym_info->subunit_index(jr) == 1) continue;
								for(Size ia = 1; ia<=pose_for_design.residue(ir).nheavyatoms(); ia++) {
									for(Size ja = 1; ja<=pose_for_design.residue(jr).nheavyatoms(); ja++) {
										if(pose_for_design.residue(ir).xyz(ia).distance_squared(pose_for_design.residue(jr).xyz(ja)) <= all_atom_contact_thresh2)	{
											// However, if the residue in question is clashing badly (usually with a
											// residue from another building block), it needs to be designed.
											core::scoring::EnergyMap em1 = pose_for_design.energies().residue_total_energies(ir);
											Real resi_fa_rep = em1[core::scoring::fa_rep];
											if(resi_fa_rep < 3.0) { contact = true; break; }
											// contact = true; break;
										}
									}
									if(contact == true) break;
								}
								if(contact == true) break;
							}
							if(!contact || option[matdes2c::mutate_intra_bb]()) nontrimer_pos.push_back(ir);
						}

						// cout << "show spheres, name ca and resi ";
						// for(int i = 1; i <= nontrimer_pos.size(); ++i) cout << (i==1?"":"+") << nontrimer_pos[i]; cout << endl;

						// Finally, filter the positions for design based on surface accessibility.
						// We don't want to design positions that are near to the interface but pointed
						// in toward the core of the monomer (e.g., positions on the insides of helices).
						// At the same time, create ResidueTypeConstraints favoring the native residue
						// at each design position if a fav_nat_bonus option is passed.
						Sizes design_pos;
						utility::vector1<core::scoring::constraints::ConstraintOP> favor_native_constraints;
						favor_native_constraints.clear();
						for(Size i = 1; i <= nontrimer_pos.size(); ++i) {
							if( sc_sasa[nontrimer_pos[i]] > 0.0 ) {
								design_pos.push_back(nontrimer_pos[i]);
								if(fav_nat_bonus != 0.0) {
									core::scoring::constraints::ConstraintOP resconstraint = new core::scoring::constraints::ResidueTypeConstraint(pose_for_design, nontrimer_pos[i], fav_nat_bonus);
									favor_native_constraints.push_back(resconstraint);
									resconstraint->show(TR);
									TR << std::endl;
								}
							}
						}
						if(fav_nat_bonus != 0.0) {
							pose_for_design.add_constraints(favor_native_constraints);
						}

						// cout << "show spheres, name ca and resi ";
						// for(int i = 1; i <= design_pos.size(); ++i) cout << (i==1?"":"+") << design_pos[i]; cout << endl;
						// utility_exit_with_message("aorisn");


						// Design
						//design(pose_for_design, sf, design_pos, true);
						utility::vector1<utility::file::FileName> resfile = option[basic::options::OptionKeys::packing::resfile]();
						design_using_resfile(pose_for_design, sf, resfile[1], design_pos);

						// 120206: Get min_rb commandline to implement rigid body minimization
      			bool min_rb = option[matdes::mutalyze::min_rb]();

						// Repack and minimize using scorefxn
						ScoreFunctionOP scorefxn = get_score_function();
						repack(pose_for_design, scorefxn, design_pos);
						minimize(pose_for_design, scorefxn, design_pos, false, true, min_rb);
						scorefxn->score(pose_for_design);

						// Build a filename for the output PDB
						std::string tag = string_of(numeric::random::uniform()).substr(2,4);
						std::ostringstream r_string;
						//r_string << std::fixed << std::setprecision(1) << "0";
						std::string fn = string_of(option[matdes::prefix]()+"_"+option[matdes::pdbID]())+"_"+ObjexxFCL::string_of(iangle1)+"_"+ObjexxFCL::string_of(iradius1)+"_"+ObjexxFCL::string_of(iangle2)+"_"+ObjexxFCL::string_of(iradius2)+"_"+tag;

						// Write the pdb file of the design
						utility::io::ozstream out( option[out::file::o]() + "/" + fn + ".pdb.gz");
						if(option[matdes2c::dump_symmetric]){
							pose_for_design.dump_pdb(out);
						} else {
							Pose pose_out;
							core::pose::symmetry::extract_asymmetric_unit(pose_for_design, pose_out);
							pose_out.dump_pdb(out);
						}
						core::io::pdb::extract_scores(pose_for_design,out);
						out.close();

						// Spit these positions out for visual debugging
						TR << "select interface_pos, " << fn << " and resi ";
						for(Size index=1; index<=interface_pos.size(); index++) {
							TR << interface_pos[index] << "+";
						}
						TR << std::endl;
						// Spit these positions out for visual debugging
						TR << "select nontrimer_pos, " << fn << " and resi ";
						for(Size index=1; index<=nontrimer_pos.size(); index++) {
							TR << nontrimer_pos[index] << "+";
						}
						TR << std::endl;
						// Spit out the final design positions for visual debugging
						TR << "select design_pos, " << fn << " and resi ";
						for(Size index=1; index<=design_pos.size(); index++) {
							TR << design_pos[index] << "+";
						}
						TR << std::endl;

						// Calculate the AverageDegree of the designed positions
						Real avg_deg = average_degree(pose_for_design, design_pos, intra_subs1, intra_subs2, p1, p2);

						// Calculate the surface area and surface complementarity for the interface
						Real int_area = 0; Real sc = 0;
						new_sc(pose_for_design, intra_subs1, intra_subs2, p1, p2, int_area, sc);

						// Get the packing score
						Real packing = get_atom_packing_score(pose_for_design, intra_subs1, intra_subs2, p1, p2, 9.0);

						// Calculate the ddG of the monomer in the assembled and unassembled states
						// protocols::simple_moves::ddG ddG_mover = protocols::simple_moves::ddG(scorefxn, 1, true);
						// ddG_mover.calculate(pose_for_design);
						Real ddG;// = ddG_mover.sum_ddG();
						//must do ddG manually, moving each BB separately

							Pose boundpose(pose_for_design);
							Pose unboundpose(pose_for_design);
							trans_pose(unboundpose,cmp1axs*1000.0,               1,p1.n_residue()               );
							trans_pose(unboundpose,cmp2axs*1000.0,p1.n_residue()+1,p1.n_residue()+p2.n_residue());
							repack(unboundpose, scorefxn, design_pos);
							minimize(unboundpose, scorefxn, design_pos, false, true, false);
							Real ubounde = scorefxn->score(unboundpose);
							Real bounde = scorefxn->score(boundpose);
							ddG = bounde - ubounde;


						// Calculate the number of hbonding groups buried by the interface, also prints them to TR
						Size nres_monomer = sym_info->num_independent_residues();
						Real buried_unsat_polars = get_unsat_polars(boundpose, unboundpose, nres_monomer, fn);


						// Calculate per-residue energies for interface residues
						Real interface_energy = 0;
						core::scoring::EnergyMap em;
						Real avg_interface_energy = 0;
						Size mutations = 0;
						for(Size index=1; index<=design_pos.size(); index++) {
							interface_energy += pose_for_design.energies().residue_total_energy(design_pos[index]);
							em += pose_for_design.energies().residue_total_energies(design_pos[index]);
							// Also, while we're looping, count the number of mutations from the input protein
							if (pose_for_design.residue(design_pos[index]).name3() != original_pose.residue(design_pos[index]).name3()) {
							mutations++;
							}
						}
						avg_interface_energy = interface_energy / design_pos.size();
						// Multiply those energies by the weights
						em *= sf->weights();

						// Create a scorefile struct, add custom metrics to it
						core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
						ss_out->fill_struct(pose_for_design,fn);
						ss_out->add_energy("ddG", ddG);
						ss_out->add_energy("air_energy", avg_interface_energy);
						ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / design_pos.size());
						ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / design_pos.size());
						ss_out->add_energy("air_fa_dun", em[core::scoring::fa_dun] / design_pos.size());
						ss_out->add_energy("des_pos", design_pos.size());
						ss_out->add_energy("mutations", mutations);
						ss_out->add_energy("unsat_pols", buried_unsat_polars);
						ss_out->add_energy("packing", packing);
						ss_out->add_energy("avg_deg", avg_deg);
						ss_out->add_energy("int_area", int_area);
						ss_out->add_energy("sc", sc);

						// Write the scorefile
						sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );

					} // iradius2
				} // iangle2
			} // iradius1
		} // iangle1


	} // iconfig

	return NULL;

}

int
main (int argc, char *argv[])
{
	try{
		register_options();
		devel::init(argc,argv);

		void* (*func)(void*) = &dostuff;

		func(NULL);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


