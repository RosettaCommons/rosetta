// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Ragul Gowthaman

// Protocol Headers
#include <devel/init.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <protocols/pockets/DarcElectrostatics.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <core/import_pose/import_pose.hh>

//reqd minimization headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/simple_pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/simple_pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <utility/file/file_sys_util.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

using namespace std;
using namespace core;
using namespace basic::options;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

//reqd minimization namespace
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace conformation;
using namespace protocols::simple_moves;
using namespace protocols::rigid;

OPT_KEY( String, protein )
OPT_KEY( StringVector, ligand )
OPT_KEY( String, ray_file )
OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Real, steric_weight )
OPT_KEY( Real, missing_point_weight )
OPT_KEY( Real, extra_point_weight )
OPT_KEY( Real, origin_cutoff )
OPT_KEY( Boolean, print_output_complex )
OPT_KEY( Boolean, search_conformers )
OPT_KEY( Boolean, darc_score_only )
OPT_KEY( Boolean, use_ligand_filename )
OPT_KEY( String, ligand_pdb_list)
//reqd minimization options
OPT_KEY( Boolean, minimize_output_complex )
OPT_KEY( Real, cst_force_constant )
//reqd delphi options
OPT_KEY( Real, delphi_grid_midpoint_x )
OPT_KEY( Real, delphi_grid_midpoint_y )
OPT_KEY( Real, delphi_grid_midpoint_z )
OPT_KEY( Integer, delphi_grid_size )
OPT_KEY( Real, delphi_grid_spacing )


static basic::Tracer TR( "apps.pilot.ragul_get_delphi_phimap.main" );

int main( int argc, char * argv [] ) {

	try {

		NEW_OPT( protein, "protein file name", "protein.pdb" );
		NEW_OPT( ligand, "input ligand(s)", "" );
		NEW_OPT( ray_file, "input eggshell(ray) triplet file name", "" );
		NEW_OPT( num_runs, "no. of runs for PSO", 100 );
		NEW_OPT( num_particles, "no. of particles for PSO", 100 );
		NEW_OPT( origin_cutoff, "value for setting minimum and maximum origin cut off", 7.0 );
		NEW_OPT( missing_point_weight, "missing point weight", 21.6 );
		NEW_OPT( extra_point_weight, "extra point weight", 9.5 );
		NEW_OPT( steric_weight, "steric weight for PSO", 1.4 );
		NEW_OPT( print_output_complex, "print DARC output ligand model with protein as PDB file", true );
		NEW_OPT( minimize_output_complex, "minimize the best scoring DARC output model", false );
		NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );
		NEW_OPT( search_conformers, "optimize conformer during docking", true );
		NEW_OPT( darc_score_only, "DARC score only (no pso docking)", false );
		NEW_OPT( use_ligand_filename, "append ligand file name to output files, instead of ligand code", false );
		NEW_OPT( ligand_pdb_list, "List of pdb files of the ligand structures", "" );
		NEW_OPT( delphi_grid_midpoint_x, "delphi_grid_midpoint_x",0.0 );
		NEW_OPT( delphi_grid_midpoint_y, "delphi_grid_midpoint_y",0.0 );
		NEW_OPT( delphi_grid_midpoint_z, "delphi_grid_midpoint_z",0.0 );
		NEW_OPT( delphi_grid_size, "delphi_grid_size", 10 );
		NEW_OPT( delphi_grid_spacing, "delphi_grid_spacing", 1.0 );

		using utility::file::FileName;
		using utility::file::file_exists;

		devel::init(argc, argv);

		std::string const input_protein = option[ protein ];
		std::string const input_ligand = "temp";//change this
		std::string const input_eggshell_triplet = option[ ray_file ];

		int particle_size = option[ num_particles ];
		int run_size = option[ num_runs ];
		core::Real const steric_wt = option[ steric_weight ];
		core::Real const missing_pt_wt = option[ missing_point_weight ];
		core::Real const extra_pt_wt = option[ extra_point_weight ];
		core::Real const origin_space = option[ origin_cutoff ];
		Size const esp_grid_size = option[ delphi_grid_size ];
		core::Real const esp_grid_spacing = option[ delphi_grid_spacing ];
		core::Real const esp_grid_midpoint_x = option[ delphi_grid_midpoint_x ];
		core::Real const esp_grid_midpoint_y = option[ delphi_grid_midpoint_y ];
		core::Real const esp_grid_midpoint_z = option[ delphi_grid_midpoint_z ];

		protocols::pockets::NonPlaidFingerprint npf;
		pose::Pose protein_pose;
		core::import_pose::pose_from_file( protein_pose, input_protein , core::import_pose::PDB_file);
		pose::Pose bound_pose = protein_pose;

		numeric::xyzVector<core::Real> protein_com(0.);
		core::Size total_atoms(0);
		for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
			core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
			if ( curr_rsd.is_protein() ) {
				for ( Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i ) {
					protein_com.x() += curr_rsd.atom(i).xyz()(1);
					protein_com.y() += curr_rsd.atom(i).xyz()(2);
					protein_com.z() += curr_rsd.atom(i).xyz()(3);
					total_atoms++;
				}
			}
		}
		protein_com /= total_atoms;
		std::cout<<"protein_CoM : "<<protein_com.x()<<" "<<protein_com.y()<<" "<<protein_com.z()<<std::endl;

		core::Real minx(999.), miny(999.), minz(999.), maxx(-999.), maxy(-999.), maxz(-999.);
		for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
			core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
			if ( curr_rsd.is_protein() ) {
				for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
					if ( curr_rsd.atom(i).xyz()(1) > maxx ) { maxx = curr_rsd.atom(i).xyz()(1);}
					if ( curr_rsd.atom(i).xyz()(1) < minx ) { minx = curr_rsd.atom(i).xyz()(1);}
					if ( curr_rsd.atom(i).xyz()(2) > maxy ) { maxy = curr_rsd.atom(i).xyz()(2);}
					if ( curr_rsd.atom(i).xyz()(2) < miny ) { miny = curr_rsd.atom(i).xyz()(2);}
					if ( curr_rsd.atom(i).xyz()(3) > maxz ) { maxz = curr_rsd.atom(i).xyz()(3);}
					if ( curr_rsd.atom(i).xyz()(3) < minz ) { minz = curr_rsd.atom(i).xyz()(3);}
				}
			}
		}

		core::Real x_from_grd_cen, y_from_grd_cen, z_from_grd_cen, x_width, y_width, z_width;
		core::Real const cen_x = (maxx + minx)/2;
		core::Real const cen_y = (maxy + miny)/2;
		core::Real const cen_z = (maxz + minz)/2;
		x_width = std::abs(maxx - minx);
		y_width = std::abs(maxy - miny);
		z_width = std::abs(maxz - minz);
		std::cout<<"max, min x : "<<maxx<<" "<<minx<<std::endl;
		std::cout<<"max, min y : "<<maxy<<" "<<miny<<std::endl;
		std::cout<<"max, min z : "<<maxz<<" "<<minz<<std::endl;
		std::cout<<"geometric_enter : "<<cen_x<<" "<<cen_y<<" "<<cen_z<<std::endl;
		std::cout<<"width : "<<x_width<<" "<<y_width<<" "<<z_width<<std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
