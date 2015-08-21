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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/magnesium/minimize_util.hh>
#include <protocols/magnesium/MgOrbitalFrameFinder.hh>
#include <protocols/magnesium/MgWaterHydrogenPacker.hh>
#include <protocols/magnesium/MgHydrater.hh>
#include <protocols/magnesium/MgMonteCarlo.hh>
#include <protocols/magnesium/MgScanner.hh>
#include <protocols/viewer/viewers.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

//silly using/typedef
// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;

static thread_local basic::Tracer TR( "mg_modeler" );

///////////////////////////////////////////////////////////////////////////////
//
// Model Mg(2+), including waters for hexahydrates -- output PDB stats
//  to help define a potential. Now use virtual atoms to help track the 'ligand field'
//  which favors octahedral coordination of lone pair donors.
//
// 1. Orient Mg(2+) ligand-field ('orbital') frames
// 2. Pack hydrogens in existing waters that ligate Mg(2+).
// 3. Pack waters around existing Mg(2+).
// 4. Sample Mg(2+) & water position by monte carlo.
// 5. Dock Mg(2+)
//
//  Rhiju, April 2015
//
///////////////////////////////////////////////////////////////////////////////

void
mg_modeler_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::magnesium;

	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	vector1< std::string > const input_pdb_files( option[ in::file::s ]() );
	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" ) );
	vector1< Size > input_pdb_mg_res = option[ magnesium::mg_res ]();
	Pose pose;

	for ( Size q = 1; q <= input_pdb_files.size(); q++ ) {
		std::string const pdb_file = input_pdb_files[ q ];
		import_pose::pose_from_pdb( pose, *rsd_set,  file_path + '/' + pdb_file );
		PoseCOP reference_pose( pose.clone() ); // can update to actual user-supplied native later.
		std::cout << "Doing input file  ==> " << pdb_file << std::endl;

		if ( q == 1 ) protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 500, 500, false );

		// mg_res usually is alist of all mg(2+) in one pose. But if there are multiple poses, assume a single res specified
		// for each -- makes benchmark code easy.
		vector1< Size > pdb_mg_res( input_pdb_mg_res );
		if ( input_pdb_files.size() > 0 && pdb_mg_res.size() > 0 ) {
			runtime_assert( input_pdb_files.size() == input_pdb_mg_res.size() );
			pdb_mg_res = utility::tools::make_vector1( input_pdb_mg_res[ q ] );
		}

		std::string tag;
		if ( option[ magnesium::fixup ]() ) {
			fixup_magnesiums( pose );
			TR << TR.Blue << "Just did basic fixup -- mg frame orientation & water-bound mg repacking." << TR.Reset << std::endl;
			TR << TR.Blue << "If you want to build mg2+-bound waters, run this app with -hydrate flag." << TR.Reset << std::endl;
			tag = "mg_fixup";
		} else if ( option[ magnesium::pack_water_hydrogens ]() ) {
			MgOrbitalFrameFinder mg_orbital_frame_finder;
			mg_orbital_frame_finder.apply( pose );
			remove_waters_except_mg_bound( pose, get_mg_water_pairs( pose, pdb_to_pose( pose, pdb_mg_res ) ) ); // note that numbering can change again
			MgWaterHydrogenPacker mg_water_hydrogen_packer( pdb_to_pose( pose,  pdb_mg_res ) );
			mg_water_hydrogen_packer.set_use_fast_heuristic( !option[ magnesium::scored_hydrogen_sampling ]() );
			mg_water_hydrogen_packer.apply( pose );
			tag = "pack_water_hydrogens";
		} else if ( option[ magnesium::hydrate ]() ) {
			remove_mg_bound_waters( pose, pdb_to_pose( pose,  pdb_mg_res ), option[ magnesium::leave_other_waters ] ); // remove any mg-bound waters
			MgHydrater mg_hydrater( pdb_to_pose( pose, pdb_mg_res ) );
			mg_hydrater.set_use_fast_frame_heuristic( !option[ magnesium::all_hydration_frames ]() );
			mg_hydrater.set_verbose( true );
			mg_hydrater.apply( pose );
			tag = "hydrate";
			if ( option[ magnesium::minimize ]() ) {
				update_mg_hoh_fold_tree( pose );
				minimize_magnesium_and_hydration_shell( pose, pdb_to_pose( pose, pdb_mg_res ), get_mg_scorefxn(), option[ magnesium::minimize_mg_coord_constraint_distance ]() );
			}
			// produce some diagnostic output.
			get_hydration_stats( pose, *reference_pose, pdb_mg_res,
				option[ out::file::o ].user() ? option[ out::file::o ]() : "default.out" );
		} else if ( option[ magnesium::monte_carlo ]() ) {
			MgMonteCarlo mg_monte_carlo;
			mg_monte_carlo.set_cycles( option[ magnesium::montecarlo::cycles ]() );
			mg_monte_carlo.set_temperature( option[ magnesium::montecarlo::temperature ]() );
			mg_monte_carlo.set_add_delete_frequency( option[ magnesium::montecarlo::add_delete_frequency ]() );
			mg_monte_carlo.set_output_pdb( option[ magnesium::montecarlo::dump ]() );
			mg_monte_carlo.apply( pose );
			tag = "monte_carlo";
		} else {
			strip_out_magnesiums( pose );
			add_single_magnesium( pose );

			// actually could encapsulate mg_positions & mg_energies as extra 'energy' values in pose. Would be way more compact.
			MgScanner mg_scanner;
			mg_scanner.set_score_cut( option[ magnesium::score_cut ] );
			ScoreFunctionOP scorefxn = ( option[ score::weights ].user() ) ? get_score_function() : get_mg_scorefxn();
			mg_scanner.set_scorefxn( scorefxn );
			mg_scanner.set_native_pose( reference_pose );
			if ( option[ in::file::input_res ].user() ) utility_exit_with_message( "Use -pose_ligand_res instead of -input_res." );
			vector1< Size > input_scan_res = option[ magnesium::ligand_res ].user() ? pdb_to_pose( pose, option[ magnesium::ligand_res ].resnum_and_chain() ) :
				vector1<Size>( option[ magnesium::pose_ligand_res ]() );
			mg_scanner.set_input_scan_res( input_scan_res );
			mg_scanner.set_silent_file( option[ out::file::silent ]() );
			mg_scanner.set_hydrate( !option[ magnesium::lores_scan ]() );
			mg_scanner.set_minimize_during_scoring( option[ magnesium::minimize_during_scoring ]() );
			mg_scanner.set_minimize( option[ magnesium::minimize ]() );
			mg_scanner.set_tether_to_closest_res( option[ magnesium::tether_to_closest_res ]() );
			mg_scanner.set_xyz_step( option[ magnesium::xyz_step ]() );
			mg_scanner.set_minimize_mg_coord_constraint_distance( option[ magnesium::minimize_mg_coord_constraint_distance ]() );
			mg_scanner.set_integration_test( option[ magnesium::integration_test ]() );
			if ( option[ out::file::o ].user() ) mg_scanner.set_output_pdb( option[ out::file::o]() );
			mg_scanner.set_score_cut_PDB( option[ magnesium::score_cut_PDB ] );

			mg_scanner.apply( pose );

			tag = "scan";
		}

		std::string outpath = option[ out::path::pdb ]();
		std::string const outfile = outpath + utility::replace_in( pdb_file, ".pdb", "." + tag + ".pdb" );
		std::cout << "Outputting: " << outfile << std::endl;

		pose.dump_pdb( outfile );
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	clock_t const my_main_time_start( clock() );
	mg_modeler_test();
	protocols::viewer::clear_conformation_viewers();
	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		option.add_relevant( out::file::silent );
		option.add_relevant( magnesium::scan );
		option.add_relevant( magnesium::mg_res );
		option.add_relevant( magnesium::minimize_during_scoring );
		option.add_relevant( magnesium::ligand_res );
		option.add_relevant( magnesium::pose_ligand_res );
		option.add_relevant( magnesium::lores_scan );
		option.add_relevant( magnesium::xyz_step );
		option.add_relevant( magnesium::score_cut );
		option.add_relevant( magnesium::score_cut_PDB );
		option.add_relevant( magnesium::integration_test );
		option.add_relevant( magnesium::tether_to_closest_res );
		option.add_relevant( magnesium::fixup );
		option.add_relevant( magnesium::pack_water_hydrogens );
		option.add_relevant( magnesium::hydrate );
		option.add_relevant( magnesium::monte_carlo );
		option.add_relevant( magnesium::scored_hydrogen_sampling );
		option.add_relevant( magnesium::all_hydration_frames );
		option.add_relevant( magnesium::leave_other_waters );
		option.add_relevant( magnesium::minimize );
		option.add_relevant( magnesium::minimize_mg_coord_constraint_distance );
		option.add_relevant( magnesium::montecarlo::temperature );
		option.add_relevant( magnesium::montecarlo::cycles );
		option.add_relevant( magnesium::montecarlo::dump );
		option.add_relevant( magnesium::montecarlo::add_delete_frequency );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////

		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
