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
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
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
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
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
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/viewer/viewers.hh>
#include <sstream>
#include <string>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
// Added 100922
// Added 101102
// Added 101103
// AUTO-REMOVED #include <basic/prof.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <pstream.h>
// AUTO-REMOVED #include <time.h>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
//Auto Headers

static basic::Tracer TR("design_symm");

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;
using core::scoring::ScoreFunctionOP;
using core::pose::Pose;
typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Mat;
typedef vector1<Size> Sizes;


void
print_movemap(core::kinematics::MoveMap const & movemap) {
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


void
design(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool hphobic_only) {

	using namespace core;
	using namespace pack;
	using namespace task;
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace chemical;

	// Set allowed aas depending on hphobic_only arg
  // This seems kind of lame to me.
  vector1<bool> allowed_aas(20, true);
  allowed_aas[aa_cys] = false;
  allowed_aas[aa_gly] = false;
	allowed_aas[aa_pro] = false;
  if(hphobic_only == true) {
    allowed_aas[aa_ala] = true;
    //allowed_aas[aa_asp] = true; // FIX LATER
    allowed_aas[aa_asp] = false; // FIX LATER
    //allowed_aas[aa_glu] = true; // FIX LATER
    allowed_aas[aa_glu] = false; // FIX LATER
    allowed_aas[aa_phe] = true;
    allowed_aas[aa_his] = false;
    allowed_aas[aa_ile] = true;
    //allowed_aas[aa_lys] = true; // FIX LATER
    allowed_aas[aa_lys] = false; // FIX LATER
    allowed_aas[aa_leu] = true;
    allowed_aas[aa_met] = true;
    allowed_aas[aa_asn] = false;
    allowed_aas[aa_pro] = false;
    allowed_aas[aa_gln] = false;
    allowed_aas[aa_arg] = false;
    allowed_aas[aa_ser] = false;
    allowed_aas[aa_thr] = false;
    allowed_aas[aa_val] = true;
    allowed_aas[aa_trp] = true;
    allowed_aas[aa_tyr] = true;
  }

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	// Set which residues can be designed
	for (Size i=1; i<=pose.n_residue(); i++) {
		if (!sym_info->bb_is_independent(i)) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if (pose.residue(i).name3() == "PRO" || pose.residue(i).name3() == "GLY") {
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

  // Actually perform design.
	make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

}

void
repack(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos) {

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
  for (Size i=1; i<=pose.n_residue(); i++) {
    if (!sym_info->bb_is_independent(i)) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else if (find(design_pos.begin(), design_pos.end(), i) == design_pos.end()) {
      task->nonconst_residue_task(i).prevent_repacking();
    } else {
      //bool temp = allowed_aas[pose.residue(i).aa()];
      //allowed_aas[pose.residue(i).aa()] = true;
			vector1<bool> allowed_aas(20, false);
      allowed_aas[pose.residue(i).aa()] = true;
      task->nonconst_residue_task(i).restrict_absent_canonical_aas(allowed_aas);
      //task->nonconst_residue_task(i).or_include_current(true);
      task->nonconst_residue_task(i).initialize_from_command_line();
      //allowed_aas[pose.residue(i).aa()] = temp;
    }
  }

  // Actually repack.
  make_symmetric_PackerTask_by_truncation(pose, task);
  protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
  packer->apply(pose);

}

void
minimize(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool move_bb, bool move_sc, bool move_rb) {

	// Initialize a MoveMap
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(move_rb);
	movemap->set_bb(false);
	movemap->set_chi(false);

	// Set allowable move types at interface positions
	// Currently, only sc moves allowed
	for (utility::vector1<Size>::iterator i = design_pos.begin(); i != design_pos.end(); i++) {
		movemap->set_bb (*i, move_bb);
		movemap->set_chi(*i, move_sc);
	}

	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	// print_movemap( *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius) {
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

// Pose needs to be scored before this will work.
void
new_sc(Pose &pose, utility::vector1<Size> intra_subs, Real& int_area, Real& sc) {
//new_sc(Pose &pose, utility::vector1<Size> intra_subs, Real& int_area, Real& sc, Real& q_int_area, Real& q_sc) {

	using namespace core;

	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	core::scoring::sc::ShapeComplementarityCalculator scc;
	//core::scoring::sc::ShapeComplementarityCalculator scq;
	//scq.settings.density = 5; // like the "quick" option
	scc.Init();
	//scq.Init();

	// Figure out which chains touch chain A, and add the residues from those chains
	// into the sc surface objects
	Size nres_monomer = symm_info->num_independent_residues();
	for (Size i=1; i<=nres_monomer; ++i) {
		scc.AddResidue(0, pose.residue(i));
		//scq.AddResidue(0, pose.residue(i));
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
				//scq.AddResidue(1, pose.residue(ir+start));
			}
		}
	}
	if (scc.Calc()) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	}
	/*if (scq.Calc()) {
		q_sc = scq.GetResults().sc;
		q_int_area = scq.GetResults().surface[2].trimmedArea;
	}
	*/
}

// Pose must be scored in order for this to work.
Pose
get_neighbor_subs (Pose const &pose, vector1<Size> intra_subs)
{

  // Figure out which chains touch chain A, and return those chains
	Pose sub_pose;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  Size nres_monomer = symm_info->num_independent_residues();
	sub_pose.append_residue_by_jump(pose.residue(1),1);
  for (Size i=2; i<=nres_monomer; ++i) {
		sub_pose.append_residue_by_bond(pose.residue(i));
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
			sub_pose.append_residue_by_jump(pose.residue(start+1),sub_pose.n_residue());
      for (Size ir=2; ir<=nres_monomer; ir++) {
				sub_pose.append_residue_by_bond(pose.residue(ir+start));
      }
    }
  }

	TR << sub_pose.fold_tree() << std::endl;

	return sub_pose;

}

Real
get_atom_packing_score (Pose const &pose, vector1<Size> intra_subs, Real cutoff=9.0)
{

	Pose sub_pose = get_neighbor_subs(pose, intra_subs);
	core::scoring::packing::HolesParams hp(basic::database::full_name("scoring/rosettaholes/decoy15.params"));
	core::scoring::packing::HolesResult hr(core::scoring::packing::compute_holes_score(sub_pose, hp));
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
  Size nres_monomer = symm_info->num_independent_residues();
	Size count = 0; Real if_score = 0;

	Real cutoff2 = cutoff*cutoff;
	for (Size ir=1; ir<=nres_monomer; ir++) {
   	for (Size ia = 1; ia<=sub_pose.residue(ir).nheavyatoms(); ia++) {
			bool contact = false;
			for (Size jr=nres_monomer+1; jr<=sub_pose.n_residue(); jr++) {
				for (Size ja = 1; ja<=sub_pose.residue(jr).nheavyatoms(); ja++) {
					if (sub_pose.residue(ir).xyz(ia).distance_squared(sub_pose.residue(jr).xyz(ja)) <= cutoff2)  {
						contact = true;
						break; // ja
					}
				} // ja
				if (contact == true) break;
			} // jr
			if (contact == true) {
				count++;
				if_score += hr.atom_scores[AtomID(ia, ir)];
			}
		} // ia
	} // ir

	return if_score / (Real)count;

}

void
*dostuff(void*) {
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

	// Get a random number tag for the tmp dir for sc calculations
	std::string sctag = string_of(numeric::random::uniform()).substr(2,4);

	// Create a score function object, get the fa_rep weight from it (default = 0.44)
	ScoreFunctionOP sf = get_score_function();
	core::scoring::methods::EnergyMethodOptions eo = sf->energy_method_options();
	eo.exclude_monomer_fa_elec(true);
	sf->set_energy_method_options(eo);
	Real orig_rep = sf->get_weight(fa_rep);

	utility::vector1<std::string> files = option[in::file::s]();
	for(Size ifile = 1; ifile <= files.size(); ++ifile) {
		std::string file = files[ifile];

		// Read in pose
		Pose pose;
		import_pose::pose_from_pdb(pose, file, resi_set);
		Pose mono = pose;
		utility::vector1<Real> sc_sasa = sidechain_sasa(mono,2.5);
		core::id::AtomID_Map<Real> bfac;
		core::pose::initialize_atomid_map(bfac,mono);
		for(Size i = 1; i <= mono.n_residue(); i++) {
			for(Size j = 1; j <= bfac.n_atom(i); j++) {
				bfac[AtomID(j,i)] = sc_sasa[i];
			}
		}


		// Parse the input filename so that the output filenames can be constructed
	  	std::vector<std::string> path_fn_vector = string_split(string_of(file), '/');
	  	std::vector<std::string> fn_vector = string_split(path_fn_vector[(path_fn_vector.size()-1)], '_');
	  	Real input_trans = (Real) atoi(fn_vector[3].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.
	  	Real input_angle = (Real) atoi(fn_vector[4].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.
		//TR << "Symm: " << fn_vector[0] << " bb: " << fn_vector[1] << " pdb: " << fn_vector[2] << " radial_disp: " << fn_vector[3].c_str() << " angle:  " << fn_vector[4]  << std::endl;
		//TR << input_trans << " " << input_angle << std::endl;

		// Handle all of the symmetry stuff
		core::pose::symmetry::make_symmetric_pose(pose);
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		std::map<Size,SymDof> dofs = sym_info->get_dofs();
	 	int sym_jump = 0;
	 	for(std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++) {
	   	Size jump_num = i->first;
	   	if (sym_jump == 0) {
		 		sym_jump = jump_num;
	   	} else {
	   		utility_exit_with_message("Can only handle one subunit!");
	   	}
	 	}
	 	if (sym_jump == 0) {
	   	utility_exit_with_message("No jump defined!");
	 	}
		Vec start_trans = pose.jump(sym_jump).get_translation();
		Mat start_rot   = pose.jump(sym_jump).get_rotation();

		// Instantiate the viewer if desired
		if (option[ basic::options::OptionKeys::parser::view ]()) protocols::viewer::add_conformation_viewer(pose.conformation(),"symm_design",1000,1000);

		// Define which subs are part of the oligomeric building block. The interfaces between these
		// subunits should not be designed.
		Sizes intra_subs;
		if (!option[in::olig_design::num_subs_building_block]()) {
			utility_exit_with_message("ERROR: You have not set the required option -in::olig_design::num_subs_building_block");
		} else {
			for(int intrasub=1; intrasub<=option[in::olig_design::num_subs_building_block](); intrasub++) {
				intra_subs.push_back(intrasub);
			}
		}

		// Set up the finer rigid body grid that we sample around each starting configuration
		Real grid_size_angle   = option[in::olig_design::grid_size_angle  ]();
		Real grid_size_radius  = option[in::olig_design::grid_size_radius ]();
		Size grid_nsamp_angle  = option[in::olig_design::grid_nsamp_angle ]();
		Size grid_nsamp_radius = option[in::olig_design::grid_nsamp_radius]();
		Real grid_start_angle   = -grid_size_angle / 2.0;
		Real grid_start_radius  = -grid_size_radius / 2.0;
		Real grid_incr_angle    = grid_size_angle  / (grid_nsamp_angle -1);
		Real grid_incr_radius   = grid_size_radius / (grid_nsamp_radius-1);
		if( grid_nsamp_angle == 1 ) {
			grid_start_angle = 0.0;
			grid_incr_angle = 0.0;
		}
		if( grid_nsamp_radius == 1 ) {
			grid_start_radius = 0.0;
			grid_incr_radius = 0.0;
		}
		if( grid_nsamp_radius == 1 && grid_nsamp_angle == 1 ) {
			pose.dump_pdb(option[out::file::o]() + "/native.pdb");
		}

		Pose const start_pose = pose;

		// Get the favor_native_residue weight from the command line
		Real fav_nat_bonus = 0.0;
		if (option[in::olig_design::fav_nat_bonus]()) {
      fav_nat_bonus = option[in::olig_design::fav_nat_bonus]();
    }

		// Iterate over the rigid body grid, performing design & minimization at each starting point
		for(Size iangle = 1; iangle <= grid_nsamp_angle; iangle++) {
			Real delta_ang = (grid_start_angle + (iangle-1)*grid_incr_angle);
			Mat rot = numeric::x_rotation_matrix_degrees(delta_ang) * start_rot;
			for(Size iradius = 1; iradius <= grid_nsamp_radius; iradius++) {
				Vec trans = start_trans + (grid_start_radius + (iradius-1)*grid_incr_radius) * Vec(1,0,0);
				// Move the pose to the current point on the rigid body grid
				core::kinematics::Jump j = start_pose.jump(sym_jump);
				j.set_translation(trans);
				j.set_rotation(rot);
				Pose pose_for_design = start_pose;
				pose_for_design.set_jump(sym_jump,j);

				// Main design/minimize loop
				// Currently, fa_rep is increased with each iteration
				Size n_iter = 1; // THIS SHOULD PROBABLY BE A COMMAND LINE ARGUMENT
				Size iter = 1;
//				for(Size iter=1; iter<=n_iter; iter++) {
					Real new_rep = pow(0.5, n_iter-iter) * orig_rep;
					sf->set_weight(fa_rep, new_rep);

					// Find out which positions are near the inter-subunit interfaces
					// These will be further screened below, then passed to design() and minimize()
					Real const contact_dist = option[in::olig_design::contact_dist]();
					Real const contact_dist_sq = contact_dist * contact_dist;
					vector1<bool> indy_resis = sym_info->independent_residues();
					Sizes interface_pos;
					for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
						if (!indy_resis[ir]) continue;
						std::string atom_i = "";
						if (pose_for_design.residue(ir).name3() == "GLY") {
							atom_i = "CA";
						} else {
							atom_i = "CB";
						}
						for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
							if (find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) != intra_subs.end()) continue;
							std::string atom_j = "";
            	if (pose_for_design.residue(jr).name3() == "GLY") {
              	atom_j = "CA";
	            } else {
  	            atom_j = "CB";
    	        }
							if (pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
							  interface_pos.push_back(ir);
							  break;
							}
						}
					}

					// Here we filter the residues that we are selecting for design
					// to get rid of those that make inter-building block interactions
					Sizes nontrimer_pos;
					Real all_atom_contact_thresh2 = 5.0*5.0;
					sf->score(pose_for_design);
					for (Size index=1; index<=interface_pos.size(); index++) {
						Size ir = interface_pos[index];
						bool contact = false;
						for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
							if (find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) == intra_subs.end()) continue;
							if (sym_info->subunit_index(jr) == 1) continue;
							for (Size ia = 1; ia<=pose_for_design.residue(ir).nheavyatoms(); ia++) {
								for (Size ja = 1; ja<=pose_for_design.residue(jr).nheavyatoms(); ja++) {
									if (pose_for_design.residue(ir).xyz(ia).distance_squared(pose_for_design.residue(jr).xyz(ja)) <= all_atom_contact_thresh2)	{
										// However, if the residue in question is clashing badly (usually with a
										// residue from another building block), it needs to be designed.
										core::scoring::EnergyMap em1 = pose_for_design.energies().residue_total_energies(ir);
										Real resi_fa_rep = em1[core::scoring::fa_rep];
										if (resi_fa_rep > 3.0) {
											TR << "CLASH: Resis " << ir << " and " << jr << std::endl;
											break;
										} else {
											contact = true;
											break;
										}
									}
								}
								if (contact == true) break;
							}
							if (contact == true) break;
						}
						if (!contact) nontrimer_pos.push_back(ir);
					}

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
							if (fav_nat_bonus != 0.0) {
	              core::scoring::constraints::ConstraintOP resconstraint = new core::scoring::constraints::ResidueTypeConstraint(pose_for_design, nontrimer_pos[i], fav_nat_bonus);
	              favor_native_constraints.push_back(resconstraint);
								resconstraint->show(TR);
								TR << std::endl;
							}
            }
          }
					if (fav_nat_bonus != 0.0) {
	        	pose_for_design.add_constraints(favor_native_constraints);
					}

					// Design
					//design(pose_for_design, sf, design_pos, true);

//				} // iter (fa_rep)

				std::string tag = string_of(numeric::random::uniform()).substr(2,4);
				std::string fn = string_of(fn_vector[0])+"_"+string_of(fn_vector[1])+"_"+string_of(fn_vector[2])+"_"+string_of(input_trans+trans.x())+"_"+string_of(input_angle+delta_ang)+"_"+tag+"_final.pdb.gz";

					// Spit these positions out for visual debugging
					TR << "select interface_pos, " << fn << " and resi ";
					for (Size index=1; index<=interface_pos.size(); index++) {
						TR << interface_pos[index] << "+";
					}
					TR << std::endl;
					// Spit these positions out for visual debugging
					TR << "select nontrimer_pos, " << fn << " and resi ";
					for (Size index=1; index<=nontrimer_pos.size(); index++) {
						TR << nontrimer_pos[index] << "+";
					}
					TR << std::endl;
					// Spit out the final design positions for visual debugging
					TR << "select design_pos, " << fn << " and resi ";
					for (Size index=1; index<=design_pos.size(); index++) {
						TR << design_pos[index] << "+";
					}
					TR << std::endl;

				// Repack using scorefxn
				ScoreFunctionOP scorefxn = get_score_function();
				repack(pose_for_design, scorefxn, design_pos);
				//scorefxn->score(pose_for_design);

				// Calculate per-residue energies for interface residues
				Real interface_energy = 0;
				core::scoring::EnergyMap em;
				Real avg_interface_energy = 0;
        for (Size index=1; index<=design_pos.size(); index++) {
					interface_energy += pose_for_design.energies().residue_total_energy(design_pos[index]);
					em += pose_for_design.energies().residue_total_energies(design_pos[index]);
        }
				avg_interface_energy = interface_energy / design_pos.size();
				// Multiply those energies by the weights
				em *= sf->weights();

				// Create a scorefile struct, add custom metrics to it
				core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
				ss_out->fill_struct(pose_for_design,fn);
				ss_out->add_energy("air_energy", avg_interface_energy);
				//ss_out->add_energy("air_fa_elec", em[core::scoring::fa_elec] / design_pos.size());
				ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / design_pos.size());
				ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / design_pos.size());
				ss_out->add_energy("air_fa_dun", em[core::scoring::fa_dun] / design_pos.size());
				ss_out->add_energy("des_pos", design_pos.size());

				// Calculate the surface area and surface complementarity for the interface
				Real int_area = 0; Real sc = 0; Real q_int_area = 0; Real q_sc = 0;
				//new_sc(pose_for_design, intra_subs, int_area, sc, q_int_area, q_sc);
				new_sc(pose_for_design, intra_subs, int_area, sc);

				// Get the packing score
				Real packing = get_atom_packing_score(pose_for_design, intra_subs, 9.0);

				// Add these to the scorefile.
				ss_out->add_energy("packing", packing);
				//ss_out->add_energy("q_int_area", q_int_area);
				//ss_out->add_energy("q_sc", q_sc);
				ss_out->add_energy("int_area", int_area);
				ss_out->add_energy("sc", sc);

				// Write the scorefile
				sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );

				// Write the pdb file of the design
				utility::io::ozstream out( option[out::file::o]() + "/" + fn );
				pose_for_design.dump_pdb(out);
				core::io::pdb::extract_scores(pose_for_design,out);
				out.close();

			} // iradius
		} // iangle
	} // ifile

	return NULL;

}

int
main (int argc, char *argv[])
{
	try{
	devel::init(argc,argv);

	void* (*func)(void*) = &dostuff;

	if (basic::options::option[ basic::options::OptionKeys::parser::view ]()) {
		protocols::viewer::viewer_main( func );
	} else {
		func(NULL);
	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



