// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    ligand_transform_with_pcs.cc
/// @brief   find a ligand position that minimizes experimental PCS data through a grid search
///          and that can be used as starting position for ligand docking
/// @details last Modified: 06/02/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Project headers
#include <protocols/nmr/pcs/PCSLigandTransformMover.hh>
#include <core/scoring/nmr/pcs/PCSData.hh>

// Core headers
#include <devel/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

#include <protocols/ligand_docking/LigandArea.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/HighResDocker.hh>
#include <protocols/ligand_docking/FinalMinimizer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/graph/Graph.hh>
#include <ObjexxFCL/format.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <map>

static basic::Tracer TR( "apps.public.nmr.ligand_transform_with_pcs" );

/// @brief Calculate the unweighted center-of-mass of this residue
core::Vector
ligand_centroid(core::conformation::Residue const & residue) {
	utility::vector1<core::Vector> coords;
	coords.reserve(residue.natoms());
	for ( core::Size i(1); i <= residue.natoms(); ++i ) {
		coords.push_back(residue.xyz(i));
	}
	return numeric::center_of_mass(coords);
}

OPT_KEY(Real, trans_stepsize)
OPT_KEY(Real, rot_stepsize)
OPT_KEY(Boolean, opt_transform)
OPT_KEY(Boolean, final_dock)
OPT_KEY(Real, ligand_radius_damping)

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols;
		using namespace protocols::moves;
		using namespace protocols::ligand_docking;
		using namespace protocols::nmr::pcs;
		using namespace core::scoring;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::rotamers;
		using namespace ObjexxFCL::format;

		// Register options
		OPT( in::file::s );
		OPT( in::file::extra_res_fa );
		OPT( nmr::pcs::input_file );
		OPT( out::nstruct );
		NEW_OPT(trans_stepsize, "Translational stepsize of grid search", 5.0);
		NEW_OPT(rot_stepsize, "Rotational stepsize of grid search", 20.0);
		NEW_OPT(opt_transform, "Further optimize initial ligand position and orientation after PCS grid search", false);
		NEW_OPT(final_dock, "Perform local docking after initial pose has been found with PCS grid search", false);
		NEW_OPT(ligand_radius_damping, "Damp the ligand radius that is used for neighbor clash check by this value", 0.0);

		// initialize core and get input
		devel::init(argc, argv);
		utility::vector1<std::string> filenames = option[ in::file::s ]();
		std::string pcs_datafile = option[ OptionKeys::nmr::pcs::input_file ]();
		core::Size nstruct = option[ out::nstruct ]();
		core::Real tstep = option[ trans_stepsize ];
		core::Real rstep = option[ rot_stepsize ];
		bool opt = option[ opt_transform ];
		bool do_docking = option[ final_dock ];
		core::Real damping = option[ ligand_radius_damping ];
		core::Real mc_cycles(5);

		TR << " * * * Starting ligand translational and rotational grid search with PCS * * * " << std::endl;

		// Setup pose, PCSData and scorefunction
		core::pose::PoseOP pose;
		if ( !(filenames.size() > 0) ) {
			utility_exit_with_message( "No input PDB file provided." );
		} else {
			pose = core::import_pose::pose_from_file(filenames[1]);
		}
		core::scoring::nmr::pcs::PCSDataOP pcs_data;
		if ( pcs_datafile.empty() ) {
			utility_exit_with_message( "No PCS input file provided." );
		} else {
			pcs_data = core::scoring::nmr::pcs::PCSDataOP(new core::scoring::nmr::pcs::PCSData(pcs_datafile, *pose));
		}
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		sfxn->set_weight(nmr_pcs,1.0);

		// Create transform movers
		PCSLigandTransformMoverOP pltm( new PCSLigandTransformMover(pcs_data,sfxn));
		pltm->set_trans_step(tstep);
		pltm->set_rot_step(rstep);
		opt ? pltm->do_optimized_transform() : pltm->undo_optimized_transform();
		pltm->set_scorefunction(sfxn);
		pltm->set_resolution_damping(damping);

		// Create docking mover
		core::Size lig_resid = core::pose::get_resnums_for_chain(*pose, 'X')[1];
		LigandAreaOP sc( new LigandArea );
		sc->chain_ = 'X';
		sc->cutoff_ = 6.0;
		sc->Calpha_restraints_ = 0.0;
		sc->minimize_ligand_ = 10.0;
		sc->tether_ligand_ = 0.0;
		sc->high_res_angstroms_ = 0.1;
		sc->high_res_degrees_ = 2.8648;
		sc->add_nbr_radius_ = true;
		sc->all_atom_mode_ = true;
		utility::vector1<LigandAreaOP> sc_areas(1,sc);

		LigandAreaOP bb( new LigandArea );
		bb->chain_ = 'X';
		bb->cutoff_ = 7.0;
		bb->Calpha_restraints_ = 0.3;
		bb->minimize_ligand_ = 10.0;
		bb->tether_ligand_ = 0.0;
		bb->high_res_angstroms_ = 0.1;
		bb->high_res_degrees_ = 2.8648;
		bb->add_nbr_radius_ = false;
		bb->all_atom_mode_ = true;
		utility::vector1<LigandAreaOP> bb_areas(1,bb);

		InterfaceBuilderOP sc_if_builder( new InterfaceBuilder(sc_areas) );
		InterfaceBuilderOP bb_if_builder_dock( new InterfaceBuilder() );
		InterfaceBuilderOP bb_if_builder_min( new InterfaceBuilder(bb_areas, 3) );
		MoveMapBuilderOP mm_builder_dock( new MoveMapBuilder(sc_if_builder, bb_if_builder_dock, false));
		MoveMapBuilderOP mm_builder_min( new MoveMapBuilder(sc_if_builder, bb_if_builder_min, false));
		HighResDockerOP hrd( new HighResDocker(mc_cycles, 1, sfxn, mm_builder_dock) );
		FinalMinimizerOP fmin( new FinalMinimizer(sfxn, mm_builder_min) );

		// Create ligand rotamer vector
		PackerTaskOP task = TaskFactory::create_packer_task(*pose);
		utility::vector1<bool> residues_to_pack(pose->total_residue(), false);
		residues_to_pack[lig_resid]=true;
		task->restrict_to_repacking();
		utility::graph::GraphOP packer_graph = create_packer_graph(*pose,*sfxn,task);
		SingleResidueRotamerLibraryFactory const & rotlib_fac( *SingleResidueRotamerLibraryFactory::get_instance() );
		SingleResidueRotamerLibraryCOP srsd_rotlib(rotlib_fac.get(pose->residue(lig_resid).type()));
		utility::vector1<core::conformation::ResidueOP> lig_rot_all;
		if ( srsd_rotlib ) {
			SingleLigandRotamerLibraryCOP slig_rotlib(utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const > ( srsd_rotlib ));
			if ( slig_rotlib ) {
				bool buried(false);
				utility::vector1< utility::vector1< core::Real > > extra_chi_steps;
				slig_rotlib->fill_rotamer_vector(*pose, *sfxn, *task, packer_graph, pose->residue(lig_resid).type_ptr(), pose->residue(lig_resid), extra_chi_steps, buried, lig_rot_all);
			} else {
				utility_exit_with_message("Failed to create single ligand rotamer library for residue " + pose->residue(lig_resid).type().name3());
			}
		} else {
			utility_exit_with_message("Failed to create single residue rotamer library for residue " + pose->residue(lig_resid).type().name3());
		}

		std::string s(pose->residue(lig_resid).type().name3()+"_centroids_pcs_dock.pdb");
		std::fstream ligand_centers_out;
		ligand_centers_out.open(s.c_str(), std::ios::out);

		for ( core::Size n(1); n <= nstruct; ++n ) {
			// pick a random rotamer
			core::Size i(numeric::random::random_range(1,lig_rot_all.size()));
			pose->replace_residue(lig_resid, *(lig_rot_all[i]), false /*orient_backbone*/);

			// Run PCS grid search
			pltm->apply(*pose);

			// Run optional local docking
			if ( do_docking ) {
				hrd->apply(*pose);
				fmin->apply(*pose);
			}

			// Output PDB
			std::string pdbfile(pose->residue(lig_resid).type().name3()+"_pcs_dock_" + utility::to_string(n) + ".pdb");
			pose->dump_pdb(pdbfile);
			core::Vector centroid(ligand_centroid(pose->residue(lig_resid)));
			if ( ligand_centers_out.is_open() ) {
				ligand_centers_out << "HETATM" << RJ(5,n) << "  O   HOH A" << RJ(4,n) << "    " << F(8,3,centroid.x())
					<< F(8,3,centroid.y()) << F(8,3,centroid.z()) << "  1.00 20.00" << RJ(12,"O") << std::endl;
			} else {
				TR.Warning << "Unable to write ligand center to file." << std::endl;
			}
		}
		if ( ligand_centers_out.is_open() ) { ligand_centers_out << "END" << std::endl; }
		ligand_centers_out.close();

		TR << " * * * DONE * * * " << std::endl;

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
