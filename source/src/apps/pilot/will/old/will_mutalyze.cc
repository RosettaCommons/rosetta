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
#include <basic/options/keys/packing.OptionKeys.gen.hh>
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
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
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
#include <protocols/pose_metric_calculators/RotamerBoltzCalculator.hh>
#include <protocols/calc_taskop_filters/RotamerBoltzmannWeight.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <basic/MetricValue.hh>
//Auto Headers

static basic::Tracer TR( "matdes::mutalyze" );

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
	protocols::moves::MoverOP packer = new protocols::minimization_packing::symmetry::SymPackRotamersMover(sf, task);
	packer->apply(pose);

}

void
design_using_resfile(Pose & pose, ScoreFunctionOP sf, std::string resfile, utility::vector1<Size> & mutalyze_pos) {

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
	parse_resfile(pose,*task, resfile);

	// Get from the packer task the mutalyze positions
	Size nres_monomer = sym_info->num_independent_residues();
	for ( Size i=1; i<=nres_monomer; ++i ) {
		if ( ! task->residue_task(i).has_behavior("AUTO") ) {
			mutalyze_pos.push_back(i);
			task->nonconst_residue_task(i).initialize_from_command_line();
		}
	}

	// Actually perform design
	make_symmetric_PackerTask_by_truncation(pose, task);
	protocols::moves::MoverOP packer = new protocols::minimization_packing::symmetry::SymPackRotamersMover(sf, task);
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
	protocols::moves::MoverOP packer = new protocols::minimization_packing::symmetry::SymPackRotamersMover(sf, task);
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
	protocols::minimization_packing::symmetry::SymMinMover m( movemap, sf, "lbfgs_armijo_nonmonotone", 1e-5, true, false, false );
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

Real
average_degree (Pose const &pose, vector1<Size> mutalyze_pos, Size intra_subs, Real distance_threshold=10.0)
{

	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	Size nres_monomer = sym_info->num_independent_residues();
	Size count_neighbors( 0 );

	for ( Size i = 1; i <= mutalyze_pos.size(); ++i ) {
		Size ires = mutalyze_pos[i];
		core::conformation::Residue const resi( pose.conformation().residue( ires ) );
		Size resi_neighbors( 0 );
		for ( Size jres = 1; jres <= (nres_monomer*intra_subs); ++jres ) {
			core::conformation::Residue const resj( pose.residue( jres ) );
			Real const distance( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) );
			if ( distance <= distance_threshold ) {
				++count_neighbors;
				++resi_neighbors;
			}
		}
		TR << "avg_deg of " << resi.name3() << ires << " = " << resi_neighbors << std::endl;
	}

	return( (Real) count_neighbors / mutalyze_pos.size() );

}

Real
get_unsat_polars( Pose const &bound, Pose const &unbound, Size nres_monomer) {

	core::pose::metrics::PoseMetricCalculatorOP unsat_calc_bound = new protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
	basic::MetricValue< id::AtomID_Map<bool> > bound_Amap;
	unsat_calc_bound->get("atom_bur_unsat", bound_Amap, bound);

	core::pose::metrics::PoseMetricCalculatorOP unsat_calc_unbound = new protocols::simple_pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("default", "default");
	basic::MetricValue< id::AtomID_Map<bool> > unbound_Amap;
	unsat_calc_unbound->get("atom_bur_unsat", unbound_Amap, unbound);

	id::AtomID_Map<bool> bound_am = bound_Amap.value();
	id::AtomID_Map<bool> unbound_am = unbound_Amap.value();
	Size buried_unsat_polars = 0;
	for ( Size ir=1; ir<=nres_monomer; ir++ ) {
		Size flag = 0;
		for ( Size ia=1; ia<=bound.residue(ir).nheavyatoms(); ia++ ) {
			if ( bound_am[id::AtomID(ia,ir)] != unbound_am[id::AtomID(ia,ir)] ) {
				buried_unsat_polars++;
				if ( flag == 0 ) {
					TR << "buried unsat polar(s): " << bound.residue(ir).name3() << ir << "\t" << bound.residue(ir).atom_name(ia);
					flag = 1;
				} else {
					TR << "," << bound.residue(ir).atom_name(ia);
				}
			}
		}
		if ( flag ) TR << std::endl;
	}

	return buried_unsat_polars;

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
		// Make a copy so we can dump original resi names at end
		Pose original_pose = pose;

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
		if ( !option[matdes::num_subs_building_block].user() ) {
			utility_exit_with_message("ERROR: You have not set the required option -matdes::num_subs_building_block");
		} else {
			for ( int intrasub=1; intrasub<=option[matdes::num_subs_building_block](); intrasub++ ) {
				intra_subs.push_back(intrasub);
			}
		}

		// Get the total number of subunits
		Size total_subs = 1;
		if ( !option[matdes::num_subs_total].user() ) {
			utility_exit_with_message("ERROR: You have not set the required option -matdes::num_subs_total");
		} else {
			total_subs = option[matdes::num_subs_total]();
		}

		// Make a scorefunction
		ScoreFunctionOP scorefxn = get_score_function();

		// Move the pose to the correct docked configuration before designing
		utility::vector1<Real> radial_disps = option[matdes::radial_disp]();
		utility::vector1<Real> angles = option[matdes::angle]();
		Mat init_rot = pose.jump(sym_jump).get_rotation();
		for ( Size iconfig=1; iconfig<=radial_disps.size(); iconfig++ ) {
			core::kinematics::Jump j = pose.jump(sym_jump);
			j.set_translation(Vec(radial_disps[iconfig],0,0));
			j.set_rotation(Mat(numeric::x_rotation_matrix_degrees(angles[iconfig]) * init_rot));
			pose.set_jump(sym_jump,j); // COMMENT OUT FOR WILLS PENTAMERS

			// Design and get mutalyze_pos and ids
			utility::vector1<Size> mutalyze_pos;
			utility::vector1<utility::file::FileName> resfile = option[basic::options::OptionKeys::packing::resfile]();
			design_using_resfile(pose, scorefxn, resfile[1], mutalyze_pos);
			utility::vector1<std::string> mutalyze_ids;
			for ( Size ipos = 1; ipos <= mutalyze_pos.size(); ++ipos ) {
				mutalyze_ids.push_back(pose.residue(mutalyze_pos[ipos]).name3());
			}

			// Repack and minimize using scorefxn
			repack(pose, scorefxn, mutalyze_pos);
			bool min_rb = option[matdes::mutalyze::min_rb]();
			minimize(pose, scorefxn, mutalyze_pos, false, true, min_rb);
			scorefxn->score(pose);

			// Write the pdb file of the design
			std::ostringstream r_string;
			r_string << std::fixed << std::setprecision(1) << radial_disps[iconfig];
			std::string fn = option[matdes::prefix]() + option[matdes::pdbID]() + "_" +r_string.str()+ "_" + string_of(angles[iconfig]) + "_" + option[matdes::tag] + ".mutalyze.pdb";
			utility::io::ozstream out( option[out::file::o]() + "/" + fn );
			pose.dump_pdb(out);
			core::io::pdb::extract_scores(pose,out);
			out.close();

			// Calculate the AverageDegree of the designed positions
			Real avg_deg = average_degree(pose, mutalyze_pos, intra_subs.size());

			// Calculate the change in SASA upon complex formation
			Real bound_sasa = core::scoring::calc_total_sasa(pose, 1.4);
			j.set_translation(Vec(1000,0,0));
			// NK 110907 for wills pentamers     j.set_translation(Vec(1000,1000,1000));
			Pose unbound_pose = pose;
			unbound_pose.set_jump(sym_jump,j);
			Real unbound_sasa = core::scoring::calc_total_sasa(unbound_pose, 1.4);
			Real buried_sasa = (unbound_sasa-bound_sasa)/total_subs;

			// Calculate the number of hbonding groups buried by the interface, also prints them to TR
			Size nres_monomer = sym_info->num_independent_residues();
			Real buried_unsat_polars = get_unsat_polars(pose, unbound_pose, nres_monomer);

			// Calculate the Boltzmann probability for the rotamer at each designed position
			if ( option[matdes::mutalyze::calc_rot_boltz]() == 1 ) {
				for ( Size ipos = 1; ipos <= mutalyze_pos.size(); ++ipos ) {
					protocols::calc_taskop_filters::RotamerBoltzmannWeight rbc = protocols::calc_taskop_filters::RotamerBoltzmannWeight();
					rbc.scorefxn(scorefxn);
					Real rot_boltz = rbc.compute_Boltzmann_weight(unbound_pose, mutalyze_pos[ipos]);
					std::cout << fn << " " << mutalyze_ids[ipos] << mutalyze_pos[ipos] << " has a rot_boltz of: " << rot_boltz << std::endl;
				}
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
			for ( Size index=1; index<=mutalyze_pos.size(); index++ ) {
				interface_energy += pose.energies().residue_total_energy(mutalyze_pos[index]);
				em += pose.energies().residue_total_energies(mutalyze_pos[index]);
			}
			avg_interface_energy = interface_energy / mutalyze_pos.size();
			// Multiply those energies by the weights
			em *= scorefxn->weights();

			// Calculate the ddG of the monomer in the assembled and unassembled states
			protocols::protein_interface_design::movers::ddG ddG_mover2 = protocols::protein_interface_design::movers::ddG(scorefxn, 1, true);
			ddG_mover2.calculate(pose);
			Real ddG2 = ddG_mover2.sum_ddG();
			TR << files[ifile] << " mutalyzed ddG = " << ddG2 << std::endl;
			ddG_mover2.report_ddG(TR);

			// Create a scorefile struct, add custom metrics to it
			core::io::silent::SilentStructOP ss_out( new core::io::silent::ScoreFileSilentStruct );
			ss_out->fill_struct(pose,fn);
			ss_out->add_energy("ddG", ddG2);
			ss_out->add_energy("air_energy", avg_interface_energy);
			ss_out->add_energy("air_fa_atr", em[core::scoring::fa_atr] / mutalyze_pos.size());
			ss_out->add_energy("air_fa_rep", em[core::scoring::fa_rep] / mutalyze_pos.size());
			ss_out->add_energy("air_fa_dun", em[core::scoring::fa_dun] / mutalyze_pos.size());
			ss_out->add_energy("unsat_pols", buried_unsat_polars);
			ss_out->add_energy("des_pos", mutalyze_pos.size());
			ss_out->add_energy("packing", packing);
			ss_out->add_energy("avg_deg", avg_deg);
			ss_out->add_energy("sasa_int_area", buried_sasa);
			ss_out->add_energy("sc_int_area", int_area);
			ss_out->add_energy("sc", sc);

			// Write the scorefile
			sfd.write_silent_struct( *ss_out, option[out::file::o]() + "/" + option[ out::file::silent ]() );

			// Loop through the design positions and mutate each non-Gly/Pro residue to alanine, then
			// calculate the ddG to get a measure of the contribution of each residue to the interface
			Sizes pos;
			utility::vector1<std::string> id;
			for ( Size ipos = 1; ipos <= mutalyze_pos.size(); ++ipos ) {
				Pose pose_for_ala_scan = pose;
				pos.clear();
				id.clear();
				pos.push_back(mutalyze_pos[ipos]);
				if ( (mutalyze_ids[ipos] == "GLY") || (mutalyze_ids[ipos] == "PRO") ) {
					id.push_back(string_of(mutalyze_ids[ipos]));
				} else {
					id.push_back("ALA");
				}

				// Design
				design(pose_for_ala_scan, scorefxn, pos, id);

				// Calculate the ddG of the monomer in the assembled and unassembled states
				protocols::protein_interface_design::movers::ddG ddG_mover3 = protocols::protein_interface_design::movers::ddG(scorefxn, 1, true);
				ddG_mover3.calculate(pose_for_ala_scan);
				Real ddG3 = ddG_mover3.sum_ddG();
				TR << files[ifile] << " ddG for mutation " << mutalyze_ids[ipos] << mutalyze_pos[ipos] << original_pose.residue(mutalyze_pos[ipos]).name3() << " = " << ddG3 << std::endl;
				ddG_mover3.report_ddG(TR);
			}
		} // for iconfig in radial_disps

	} // ifile

	return NULL;

}

int
main (int argc, char *argv[])
{

	try {


		devel::init(argc,argv);

		void* (*func)(void*) = &dostuff;

		func(NULL);


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


