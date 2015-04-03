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
// Added 100922
// Added 101102
// Added 101103
//Auto Headers

static thread_local basic::Tracer TR( "design_symm" );

using std::string;
using ObjexxFCL::string_of;
using namespace core;
using utility::vector1;
using core::id::AtomID;
using namespace core::scoring::packing;
using namespace ObjexxFCL::format;
using core::scoring::ScoreFunctionOP;
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
design(core::pose::Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool hphobic_only) {

	using namespace core;
	using namespace pack;
	using namespace task;
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace scoring;
	using namespace chemical;

/*
	// Set allowed aas depending on hphobic_only arg
  vector1<bool> allowed_aas(20, true);
  allowed_aas[aa_cys] = false;
  allowed_aas[aa_gly] = false;
	allowed_aas[aa_pro] = false;
  if(hphobic_only == true) {
    allowed_aas[aa_ala] = true;
    allowed_aas[aa_asp] = true; // FIX LATER
    allowed_aas[aa_glu] = true; // FIX LATER
    allowed_aas[aa_phe] = true;
    allowed_aas[aa_his] = false;
    allowed_aas[aa_ile] = true;
    allowed_aas[aa_lys] = true; // FIX LATER
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
*/

	// Get the symmetry info and make the packer task
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	PackerTaskOP task( TaskFactory::create_packer_task( pose ));

	// Set which residues can be designed
	for (Size i=1; i<=pose.n_residue(); i++) {
		task->nonconst_residue_task(i).prevent_repacking();
/*
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
*/
	}

  // Actually perform design.
	make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::simple_moves::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

}

utility::vector1<Real> sidechain_sasa(core::pose::Pose const & pose, Real probe_radius) {
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

void
get_sc(core::pose::Pose &pose, Size nums_subs_bb, Real& int_area, Real& sc, std::string sctag) {

	using namespace core;
	using namespace utility::io;
	using namespace core::io::pdb;

	// Set up environment variables and directories for calculating sc
	// locally on syd nodes.
	char * cwd = getenv("PWD");
	std::string blah1 = "/scratch/USERS/neil/sc_tmp_"+sctag+"/";
	char *tmpdir = (char*)blah1.c_str();
	std::string blah2 = "mkdir -p /scratch/USERS/neil/sc_tmp_"+sctag+"/";
	char *mkdir = (char*)blah2.c_str();
	setenv("PATH", "/work/neilpking/bin/ccp4-6.0.2/bin:/bin", 1);
	setenv("CINCL", "/work/neilpking/bin/ccp4-6.0.2/include", 1);
	setenv("CCP4_SCR", "/scratch/USERS/neil/ccp4_scratch", 1);
	setenv("CLIBD", "/work/neilpking/bin/ccp4-6.0.2/lib/data", 1);
	system("mkdir -p /scratch/USERS/neil/ccp4_scratch/");
	//system("mkdir -p /scratch/USERS/neil/sc_tmp/");
	//chdir("/scratch/USERS/neil/sc_tmp/");
	system(mkdir);
	chdir(tmpdir);

	// Dump the temporary pdb (with scores) that calc_sc_inline.php will use.
	std::string tmppdbfn = "tmpfull.pdb";
	utility::io::ozstream out( tmppdbfn );
  pose.dump_pdb(out);
  core::io::pdb::extract_scores(pose,out);
  out.close();

	// Run sc, get the int_area and sc values out of the logfile.
	std::string php_string = "/usr/bin/php /work/neilpking/bin/calc_sc_inline.php tmpfull.pdb 3 "+blah1;
	char *php_char = (char*)php_string.c_str();
	system(php_char);
	utility::io::izstream scfile("sc.log");
	std::string line;
	while( getline(scfile,line) ) {
		std::istringstream line_stream(line);
		line_stream >> int_area >> sc;
		if ( line_stream.fail() ) {
			std::cout << "Error reading two floats from " << "sc.log" << " line " << line << std::endl;
		}
	}
	scfile.close();

	// Move back to where you were.
	chdir(cwd);

	// Clean up.
	//system("rm -rf /scratch/USERS/neil/ccp4_scratch/");
	//system("rm -rf /scratch/USERS/neil/sc_tmp/");

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
			Real input_trans = 0;
			Real input_angle = 0;
//	  	Real input_trans = (Real) atoi(fn_vector[3].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.
//	  	Real input_angle = (Real) atoi(fn_vector[4].c_str()); // THIS REALLY NEEDS TO BE FIXED TO HANDLE NON-INT VALUES.

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

		// Create a score function object, get the fa_rep weight from it (default = 0.44)
		ScoreFunctionOP sf = get_score_function();
		core::scoring::methods::EnergyMethodOptions eo = sf->energy_method_options();
		eo.exclude_monomer_fa_elec(true);
		sf->set_energy_method_options(eo);
		Real orig_rep = sf->get_weight(fa_rep);

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

					Sizes nontrimer_pos;
					// Here we filter the residues that we are selecting for design
					// to get rid of those that make inter-building block interactions
					Real all_atom_contact_thresh2 = 5.0*5.0;
					for (Size index=1; index<=interface_pos.size(); index++) {
						Size ir = interface_pos[index];
						bool contact = false;
						for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
							if (find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) == intra_subs.end()) continue;
							if (sym_info->subunit_index(jr) == 1) continue;
							for (Size ia = 1; ia<=pose_for_design.residue(ir).nheavyatoms(); ia++) {
								for (Size ja = 1; ja<=pose_for_design.residue(jr).nheavyatoms(); ja++) {
									if (pose_for_design.residue(ir).xyz(ia).distance_squared(pose_for_design.residue(jr).xyz(ja)) <= all_atom_contact_thresh2)	{
										contact = true;
										break;
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
					// in toward the core of the monomer (e.g., positions on the insides of helices)
					Sizes design_pos;
          // Spit out the final design positions for visual debugging
          for (Size index=1; index<=design_pos.size(); index++) {
            TR << design_pos[index] << "+";
          }
          TR << std::endl;
					for(Size i = 1; i <= nontrimer_pos.size(); ++i) {
						if( sc_sasa[nontrimer_pos[i]] > 0.0 ) {
							design_pos.push_back(nontrimer_pos[i]);
						}
					}

					// Design
					design(pose_for_design, sf, design_pos, true);

				std::string tag = string_of(numeric::random::uniform()).substr(2,4);
				std::string fn = string_of(fn_vector[0])+"_"+string_of(fn_vector[1])+"_"+string_of(fn_vector[2])+"_"+string_of(input_trans+trans.x())+"_"+string_of(input_angle+delta_ang)+"_"+tag+"_native.pdb.gz";
//				std::string fn = string_of(fn_vector[0])+"_"+string_of(fn_vector[1])+"_"+string_of(fn_vector[2])+"_native.pdb.gz";

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

				// Calculate per-residue energies for interface residues
				Real interface_energy = 0;
				core::scoring::EnergyMap em;
				Real avg_interface_energy = 0;
        for (Size index=1; index<=design_pos.size(); index++) {
					interface_energy += pose_for_design.energies().residue_total_energy(design_pos[index]);
					em += pose_for_design.energies().residue_total_energies(design_pos[index]);
        }
				avg_interface_energy = interface_energy / design_pos.size();
				//em = em / design_pos.size();

				// Create a scorefile struct, add custom metrics to it
				core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
				ss_out->fill_struct(pose_for_design,fn);
				ss_out->add_energy("air_energy", avg_interface_energy);
				ss_out->add_energy("air_fa_elec", em[core::scoring::fa_elec] / design_pos.size());
				ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / design_pos.size());
				ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / design_pos.size());
				ss_out->add_energy("air_fa_dun", em[core::scoring::fa_dun] / design_pos.size());
				ss_out->add_energy("des_pos", design_pos.size());

				// Calculate the surface area and surface complementarity for the interface
				Real int_area = 0;
				Real sc = 0;
				get_sc(pose_for_design, int_area, sc, sctag);

				// Add these to the scorefile.
				ss_out->add_energy("int_area", int_area);
				ss_out->add_energy("sc", sc);

				// Write the scorefile
				sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );

				// Write the pdb file of the design
				utility::io::ozstream out( option[out::file::o]() + "/" + fn );
				pose_for_design.dump_pdb(out);
				//core::scoring::packstat::output_packstat_pdb(pose_for_design, out); // Currently too slow for many structures
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


