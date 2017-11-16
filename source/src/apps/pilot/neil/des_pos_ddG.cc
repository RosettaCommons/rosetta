// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/options/keys/holes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
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
#include <core/io/pdb/pdb_writer.hh>
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
#include <protocols/protein_interface_design/movers/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>
//Auto Headers

static basic::Tracer TR( "des_pos_ddG" );

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
design(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> revert_pos, utility::vector1<string> revert_ids) {

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

	// Set which residues can be designed
	utility::vector1<bool> allowed_aas(20, false);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( find(revert_pos.begin(), revert_pos.end(), i) == revert_pos.end() ) {
			task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	for ( Size ipos = 1; ipos <= revert_pos.size(); ipos++ ) {
		string aa_name = revert_ids[ipos];
		allowed_aas[aa_from_name(aa_name)] = true;
		task->nonconst_residue_task(revert_pos[ipos]).restrict_absent_canonical_aas(allowed_aas);
		task->nonconst_residue_task(revert_pos[ipos]).initialize_from_command_line();
		allowed_aas[aa_from_name(aa_name)] = false;
	}

	// Actually perform design
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
	for ( Size i=1; i<=pose.size(); i++ ) {
		if ( !sym_info->bb_is_independent(i) ) {
			task->nonconst_residue_task(i).prevent_repacking();
		} else if ( find(design_pos.begin(), design_pos.end(), i) == design_pos.end() ) {
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

void
minimize(Pose & pose, ScoreFunctionOP sf, utility::vector1<Size> design_pos, bool move_bb, bool move_sc, bool move_rb) {

	// Initialize a MoveMap
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
	movemap->set_jump(move_rb);
	movemap->set_bb(false);
	movemap->set_chi(false);

	// Set allowable move types at interface positions
	// Currently, only sc moves allowed
	for ( utility::vector1<Size>::iterator i = design_pos.begin(); i != design_pos.end(); i++ ) {
		movemap->set_bb (*i, move_bb);
		movemap->set_chi(*i, move_sc);
	}

	// Make MoveMap symmetric, apply it to minimize the pose
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );
	// print_movemap( *movemap );
	protocols::simple_moves::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
	m.apply(pose);
}

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius) {
	using core::id::AtomID;
	utility::vector1<Real> rsd_sasa(pose.size(),0.0);
	core::id::AtomID_Map<Real> atom_sasa;
	core::id::AtomID_Map<bool> atom_mask;
	core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
	core::pose::initialize_atomid_map(atom_mask,pose,false);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j <= pose.residue(i).nheavyatoms(); j++ ) {
			atom_mask[AtomID(j,i)] = true;
		}
	}
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );
	utility::vector1<Real> sc_sasa(pose.size(),0.0);
	for ( Size i = 1; i <= pose.size(); i++ ) {
		// Use CA as the side chain for Glys
		if ( pose.residue(i).name3()=="GLY" ) sc_sasa[i] += atom_sasa[AtomID(2,i)];
		for ( Size j = 5; j <= pose.residue(i).nheavyatoms(); j++ ) {
			sc_sasa[i] += atom_sasa[AtomID(j,i)];
		}
	}
	return sc_sasa;
}

// Pose needs to be scored before this will work.
void
new_sc(Pose &pose, utility::vector1<Size> intra_subs, Real& int_area, Real& sc) {

	using namespace core;

	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	core::scoring::sc::ShapeComplementarityCalculator scc;
	scc.Init();

	// Figure out which chains touch chain A, and add the residues from those chains
	// into the sc surface objects
	Size nres_monomer = symm_info->num_independent_residues();
	for ( Size i=1; i<=nres_monomer; ++i ) {
		scc.AddResidue(0, pose.residue(i));
	}
	for ( Size i=1; i<=symm_info->subunits(); ++i ) {
		if ( std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end() ) continue;
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for ( Size ir=1; ir<=nres_monomer; ir++ ) {
			if ( pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0 ) {
				contact = true;
				break;
			}
		}
		if ( contact ) {
			for ( Size ir=1; ir<=nres_monomer; ir++ ) {
				scc.AddResidue(1, pose.residue(ir+start));
			}
		}
	}
	if ( scc.Calc() ) {
		sc = scc.GetResults().sc;
		int_area = scc.GetResults().surface[2].trimmedArea;
	}
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
	for ( Size i=2; i<=nres_monomer; ++i ) {
		sub_pose.append_residue_by_bond(pose.residue(i));
	}
	for ( Size i=1; i<=symm_info->subunits(); ++i ) {
		if ( std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end() ) continue;
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for ( Size ir=1; ir<=nres_monomer; ir++ ) {
			if ( pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0 ) {
				contact = true;
				break;
			}
		}
		if ( contact ) {
			sub_pose.append_residue_by_jump(pose.residue(start+1),sub_pose.size());
			for ( Size ir=2; ir<=nres_monomer; ir++ ) {
				sub_pose.append_residue_by_bond(pose.residue(ir+start));
			}
		}
	}

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
	for ( Size ir=1; ir<=nres_monomer; ir++ ) {
		for ( Size ia = 1; ia<=sub_pose.residue(ir).nheavyatoms(); ia++ ) {
			bool contact = false;
			for ( Size jr=nres_monomer+1; jr<=sub_pose.size(); jr++ ) {
				for ( Size ja = 1; ja<=sub_pose.residue(jr).nheavyatoms(); ja++ ) {
					if ( sub_pose.residue(ir).xyz(ia).distance_squared(sub_pose.residue(jr).xyz(ja)) <= cutoff2 )  {
						contact = true;
						break; // ja
					}
				} // ja
				if ( contact == true ) break;
			} // jr
			if ( contact == true ) {
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

	// Iterate through files
	utility::vector1<std::string> files = option[in::file::s]();
	for ( Size ifile = 1; ifile <= files.size(); ++ifile ) {
		std::string file = files[ifile];

		// Read in pose
		Pose pose;
		import_pose::pose_from_file(pose, file, resi_set, core::import_pose::PDB_file);

		// Handle all of the symmetry stuff
		core::pose::symmetry::make_symmetric_pose(pose);
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		std::map<Size,SymDof> dofs = sym_info->get_dofs();
		int sym_jump = 0;
		for ( std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++ ) {
			Size jump_num = i->first;
			if ( sym_jump == 0 ) {
				sym_jump = jump_num;
			} else {
				utility_exit_with_message("Can only handle one subunit!");
			}
		}
		if ( sym_jump == 0 ) {
			utility_exit_with_message("No jump defined!");
		}

		// Define which subs are part of the oligomeric building block.
		Sizes intra_subs;
		if ( !option[matdes::design::num_subs_building_block].user() ) {
			utility_exit_with_message("ERROR: You have not set the required option -matdes::design::num_subs_building_block");
		} else {
			for ( int intrasub=1; intrasub<=option[matdes::design::num_subs_building_block](); intrasub++ ) {
				intra_subs.push_back(intrasub);
			}
		}

		// Make a scorefunction
		ScoreFunctionOP scorefxn = get_score_function();

		// Calculate the ddG of the input structure in the assembled and unassembled states
		protocols::protein_interface_design::movers::ddG ddG_mover = protocols::protein_interface_design::movers::ddG(scorefxn, 1, true);
		ddG_mover.calculate(pose);
		Real ddG = ddG_mover.sum_ddG();
		TR << files[ifile] << "designed ddG = " << ddG << std::endl;
		ddG_mover.report_ddG(TR);

		// Get the design positions and identities from the command line
		Sizes revert_pos = option[matdes::design::revert_pos]();
		utility::vector1<std::string> revert_ids = option[matdes::design::revert_ids]();

		// Design
		design(pose, scorefxn, revert_pos, revert_ids);

		// Repack and minimize using scorefxn
		repack(pose, scorefxn, revert_pos);
		minimize(pose, scorefxn, revert_pos, false, true, false);
		scorefxn->score(pose);

		// Write the pdb file of the design
		std::vector<std::string> path_fn_vector = string_split(string_of(file), '/');
		std::string fn = path_fn_vector[(path_fn_vector.size()-1)] + ".des_pos_stats.pdb";
		utility::io::ozstream out( option[out::file::o]() + "/" + fn );
		pose.dump_pdb(out);
		core::io::pdb::extract_scores(pose,out);
		out.close();

		// Calculate the change in SASA upon complex formation
		Real bound_sasa = core::scoring::calc_total_sasa(pose, 1.4);
		core::kinematics::Jump j = pose.jump(sym_jump);
		j.set_translation(Vec(1000,0,0));
		Pose unbound_pose = pose;
		unbound_pose.set_jump(sym_jump,j);
		Real unbound_sasa = core::scoring::calc_total_sasa(unbound_pose, 1.4);
		Real buried_sasa = (unbound_sasa-bound_sasa)/24;

		// Calculate the Boltzmann probability for the rotamer at each designed position
		for ( Size ipos = 1; ipos <= revert_pos.size(); ++ipos ) {
			protocols::simple_filters::RotamerBoltzmannWeight rbc = protocols::simple_filters::RotamerBoltzmannWeight();
			rbc.scorefxn(scorefxn);
			Real rot_boltz = rbc.compute_Boltzmann_weight(unbound_pose, revert_pos[ipos]);
			std::cout << path_fn_vector[(path_fn_vector.size()-1)] << " " << revert_ids[ipos] << revert_pos[ipos] << " has a rot_boltz of: " << rot_boltz << std::endl;
		}

		// Calculate the surface area and surface complementarity for the interface
		Real int_area = 0; Real sc = 0;
		new_sc(pose, intra_subs, int_area, sc);

		// Get the packing score
		Real packing = get_atom_packing_score(pose, intra_subs, 9.0);

		// Calculate per-residue energies for interface residues
		Real interface_energy = 0;
		core::scoring::EnergyMap em;
		Real avg_interface_energy = 0;
		for ( Size index=1; index<=revert_pos.size(); index++ ) {
			interface_energy += pose.energies().residue_total_energy(revert_pos[index]);
			em += pose.energies().residue_total_energies(revert_pos[index]);
		}
		avg_interface_energy = interface_energy / revert_pos.size();
		// Multiply those energies by the weights
		em *= scorefxn->weights();

		// Calculate the ddG of the monomer in the assembled and unassembled states
		protocols::protein_interface_design::movers::ddG ddG_mover = protocols::protein_interface_design::movers::ddG(scorefxn, 1, true);
		ddG_mover.calculate(pose);
		Real ddG = ddG_mover.sum_ddG();
		TR << files[ifile] << "reverted ddG = " << ddG << std::endl;
		ddG_mover.report_ddG(TR);

		// Create a scorefile struct, add custom metrics to it
		core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
		ss_out->fill_struct(pose,fn);
		ss_out->add_energy("ddG", ddG);
		ss_out->add_energy("air_energy", avg_interface_energy);
		ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / revert_pos.size());
		ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / revert_pos.size());
		ss_out->add_energy("air_fa_dun", em[core::scoring::fa_dun] / revert_pos.size());
		ss_out->add_energy("des_pos", revert_pos.size());
		ss_out->add_energy("packing", packing);
		ss_out->add_energy("sasa_int_area", buried_sasa);
		ss_out->add_energy("sc_int_area", int_area);
		ss_out->add_energy("sc", sc);

		// Write the scorefile
		sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );

		// Loop through the design positions and mutate each non-glycine residue to alanine,
		// calculate the ddG, and subtract from the designed (reverted) ddG to get a measure
		// of the contribution of each residue to the interface
		Sizes pos;
		utility::vector1<std::string> id;
		id.push_back("ALA");
		for ( Size ipos = 1; ipos <= revert_pos.size(); ++ipos ) {
			Pose pose_for_ala_scan = pose;
			pos.clear();
			pos.push_back(revert_pos[ipos]);
			if ( (revert_ids[ipos] == "GLY") || (revert_ids[ipos] == "PRO") ) {
				continue;
			}
			// Design
			design(pose_for_ala_scan, scorefxn, pos, id);

			// Calculate the ddG of the monomer in the assembled and unassembled states
			protocols::protein_interface_design::movers::ddG ddG_mover = protocols::protein_interface_design::movers::ddG(scorefxn, 1, true);
			ddG_mover.calculate(pose_for_ala_scan);
			Real ddG = ddG_mover.sum_ddG();
			TR << files[ifile] << " ddG for mutation " << revert_ids[ipos] << revert_pos[ipos] << "ALA = " << ddG << std::endl;
			ddG_mover.report_ddG(TR);
		}

	} // ifile

	return NULL;

}

int
main (int argc, char *argv[])
{
	try{
		devel::init(argc,argv);

		void* (*func)(void*) = &dostuff;

		if ( basic::options::option[ basic::options::OptionKeys::parser::view ]() ) {
			protocols::viewer::viewer_main( func );
		} else {
			func(NULL);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


