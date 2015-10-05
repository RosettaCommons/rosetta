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
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/database/open.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/align/util.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/EtableEnergyCreator.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/BaseEtableEnergy.hh>
#include <core/scoring/etable/BaseEtableEnergy.tmpl.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/atom_pair_energy_inline.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

OPT_KEY( String,  params_file )
OPT_KEY( StringVector, original_input )
OPT_KEY( Boolean, virtualize_free )
OPT_KEY( Boolean, color_by_score )
OPT_KEY( Boolean, soft_rep )

// Move this out to protocols/toolbox/ before checkin.
//
// Could pretty easily get this to read in an actual scorefunction, rather than some
// hard wired combination of fa_atr/fa_rep. Problem is that most scorefunctions
// do not provide functions that compute energies at atom level. Even ones that do
// (geometric_solvation, etable energies, hbond) do not have functions with stereotyped
// input/output. Can fill in by hand, one by one.
//
void
do_color_by_score( core::pose::Pose & pose ) {

	using namespace core::scoring;
	using namespace core::scoring::methods;
	using namespace core::scoring::etable;
	using namespace core::scoring::etable::count_pair;
	using namespace core::conformation;
	using namespace core::id;

	ScoreFunction scorefxn;
	EnergyMethodOptions options( scorefxn.energy_method_options() );
	options.etable_options().no_lk_polar_desolvation = true;
	if ( option[ soft_rep ]() ) options.etable_type( "FA_STANDARD_SOFT" );
	scorefxn.set_energy_method_options( options );
	scorefxn.set_weight( fa_atr, 0.21 );
	scorefxn.set_weight( fa_rep, 0.20 );
	scorefxn.set_weight( fa_sol, 0.25 );
	//scorefxn.set_weight( fa_rep, 1.0);

	core::scoring::etable::EtableCOP etable(core::scoring::ScoringManager::get_instance()->etable( options ) );
	core::scoring::etable::AnalyticEtableEvaluator eval( *etable );
	eval.set_scoretypes( fa_atr, fa_rep, fa_sol );

	EnergyMap emap, emap_total;

	//////////////////////////////////////////////////////////////////
	// stolen from core/scoring/etable/atom_pair_energy_inline.hh
	//////////////////////////////////////////////////////////////////
	DistanceSquared dsq;
	Real weight;
	Size path_dist;
	typedef utility::vector1< Size > const & vect;

	for ( Size m = 1; m <= pose.total_residue(); m++ ) {

		Residue const & res1 = pose.residue( m );

		// get hydrogen interaction cutoff
		Real const Hydrogen_interaction_cutoff2
			( eval.hydrogen_interaction_cutoff2() );

		for ( Size i = 1; i <= res1.natoms(); i++ ) {
			pose.pdb_info()->temperature( m, i, 0.0 );
		}

		Size const res1_start( 1 ), res1_end( res1.nheavyatoms() );
		vect r1hbegin( res1.attached_H_begin() );
		vect r1hend(   res1.attached_H_end()   );

		// Atom pairs
		for ( int i = res1_start, i_end = res1_end; i <= i_end; ++i ) {
			Atom const & atom1( res1.atom(i) );
			//get virtual information
			bool atom1_virtual(res1.atom_type(i).is_virtual());

			emap.zero();

			for ( Size n = 1; n <= pose.total_residue(); n++ ) {

				if ( m == n ) continue; // later could be smart about holding info in fa_intra terms of emap
				Residue const & res2 = pose.residue( n );

				// costly, but I'm having problems with CountPairFactory output.
				// std::cout << m << " " << n << std::endl;
				CountPairGeneric count_pair( res1, res2 );
				count_pair.set_crossover( 4 );
				Size const res2_start( 1 ), res2_end( res2.nheavyatoms() );
				vect r2hbegin( res2.attached_H_begin() );
				vect r2hend(   res2.attached_H_end()   );

				for ( int j=res2_start, j_end = res2_end; j <= j_end; ++j ) {

					Atom const & atom2( res2.atom(j) );

					//     if ( ! count_pair( i, j, weight, path_dist ) ) continue;
					//     eval.atom_pair_energy( atom1, atom2, weight, emap, dsq );

					//     if ( dsq < 16.0 ) std::cout << dsq << " " << weight << " " << emap[ fa_atr ] << std::endl;
					//check if virtual
					bool atom2_virtual(res2.atom_type(j).is_virtual());
					weight = 1.0;
					path_dist = 0;
					if ( atom1_virtual || atom2_virtual ) {
						// NOOP! etable_energy.virtual_atom_pair_energy(emap);
					} else {
						if ( count_pair( i, j, weight, path_dist ) ) {
							eval.atom_pair_energy( atom1, atom2, weight, emap, dsq );
						} else {
							dsq = atom1.xyz().distance_squared( atom2.xyz() );
						}
						if ( dsq < Hydrogen_interaction_cutoff2 ) {
							residue_fast_pair_energy_attached_H(
								res1, i, res2, j,
								r1hbegin[ i ], r1hend[ i ],
								r2hbegin[ j ], r2hend[ j ],
								count_pair, eval , emap);

						}
					} // virtual
				} // j
			} // n

			emap *= 0.5; /* double counting */
			emap_total += emap;

			Real score = emap.dot(  scorefxn.weights() );
			pose.pdb_info()->temperature( m, i, score );


		} // i
	} // m


	( scorefxn )( pose );
	for ( Size n = 1; n <= n_score_types; n++ ) {
		ScoreType st( static_cast< ScoreType >( n ) );
		if ( scorefxn.has_nonzero_weight( st ) ) {
			std::cout << st << ":   conventional score function " << pose.energies().total_energies()[ st ] << "   atomwise " << emap_total[ st ] << std::endl;
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
rna_score_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring::rna::data;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose::full_model_info;
	using namespace protocols::stepwise::modeler;
	using namespace protocols::stepwise::setup;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD /*RNA*/ );

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				) );
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	// native pose setup
	pose::PoseOP native_pose;
	bool native_exists( false );
	if ( option[ in::file::native ].user() ) {
		native_pose = get_pdb_with_full_model_info( option[ in::file::native ](), rsd_set );
		native_exists = true;
	}

	// score function setup
	core::scoring::ScoreFunctionOP scorefxn;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn = core::scoring::get_score_function();
	} else {
		scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}

	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileData silent_file_data;

	FullModelInfoOP my_model;
	utility::vector1< pose::PoseOP > other_poses;

	if ( option[ original_input ].user() ) {
		utility::vector1< std::string > const & original_files = option[ original_input ]();
		utility::vector1< pose::PoseOP > original_poses;

		for ( Size n = 1; n <= original_files.size(); n++ ) {
			original_poses.push_back( get_pdb_and_cleanup( original_files[ n ], rsd_set ) );
		}
		if ( option[ full_model::other_poses ].user() ) get_other_poses( original_poses, option[ full_model::other_poses ](), rsd_set );

		//FullModelInfo (minimal object needed for add/delete)
		fill_full_model_info_from_command_line( original_poses );
		my_model = const_full_model_info( *original_poses[ 1 ] ).clone_info();

	} else {
		// other poses -- for scoring collections of poses connected by (virtual) loops, using full_model_info.
		if ( option[ full_model::other_poses ].user() ) get_other_poses( other_poses, option[ full_model::other_poses ](), rsd_set );
	}

	// if trying to compute stem RMSD
	protocols::farna::RNA_StructureParameters parameters;
	core::io::rna::RNA_DataReader rna_data_reader( option[OptionKeys::rna::data_file ]() );
	RNA_ChemicalMappingEnergyOP rna_chemical_mapping_energy;
	pose::Pose pose,start_pose;

	Size i( 0 );

	if ( native_exists ) ( *scorefxn)( *native_pose );

	while ( input->has_another_pose() ) {

		input->fill_pose( pose, *rsd_set );
		i++;

		if ( option[ virtualize_free ]() ) protocols::stepwise::modeler::rna::virtualize_free_rna_moieties( pose ); // purely for testing.

		if ( !option[ in::file::silent ].user() ) cleanup( pose );

		if ( !full_model_info_defined( pose ) || option[ in::file::fasta ].user() ) {
			if ( ! option[ original_input ].user() ) {
				fill_full_model_info_from_command_line( pose, other_poses ); // only does something if -in:file:fasta specified.
			} else {
				// allows for scoring of PDB along with 'other_poses' supplied from command line. Was used to test loop_close score term.
				utility::vector1< Size > resnum;
				core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
				if ( pdb_info ) {
					//std::cout << std::endl << "PDB Info available for this pose..." << std::endl << std::endl;
					for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( pdb_info->number( n ) );
				} else {
					for ( Size n = 1; n <= pose.total_residue(); n++ ) resnum.push_back( n );
				}
				my_model->set_res_list( resnum );
				my_model->set_other_pose_list( other_poses );
				set_full_model_info( pose, my_model );
			}
		}

		if ( option[params_file].user() ) {
			parameters.initialize(pose, option[params_file],
				basic::database::full_name("sampling/rna/1jj2_RNA_jump_library.dat"),
				false /*ignore_secstruct*/ );
			parameters.setup_base_pair_constraints( pose );
		}

		// do it
		if ( ! option[ score::just_calc_rmsd]() && !rna_data_reader.has_reactivities() ) {
			(*scorefxn)( pose );
		}

		// tag
		std::string tag = tag_from_pose( pose );
		BinarySilentStruct s( pose, tag );

		if ( native_exists ) {
			//Real const rmsd      = all_atom_rmsd( *native_pose, pose );
			Real const rmsd = protocols::stepwise::modeler::align::superimpose_with_stepwise_aligner( pose, *native_pose, option[ OptionKeys::stepwise::superimpose_over_all ]() );
			std::cout << "All atom rmsd over moving residues: " << tag << " " << rmsd << std::endl;
			s.add_energy( "new_rms", rmsd );

			// Stem RMSD
			if ( option[params_file].user() ) {
				std::list< Size > stem_residues( parameters.get_stem_residues( pose ) );
				if ( !stem_residues.empty()/*size() > 0*/ ) {
					Real const rmsd_stems = all_atom_rmsd( *native_pose, pose, stem_residues );
					s.add_energy( "rms_stem", rmsd_stems );
					std::cout << "Stems rmsd: " << rmsd_stems << std::endl;
				}
			}
		}

		// for data_file, don't actually re-score, just compute rna_chem_map score for now.
		if ( rna_data_reader.has_reactivities() ) {
			if ( rna_chemical_mapping_energy == 0 ) rna_chemical_mapping_energy = RNA_ChemicalMappingEnergyOP( new RNA_ChemicalMappingEnergy );
			rna_data_reader.fill_rna_data_info( pose );
			pose.update_residue_neighbors();
			s.add_energy(  "rna_chem_map",       rna_chemical_mapping_energy->calculate_energy( pose, false /*use_low_res*/ ) );
			s.add_energy(  "rna_chem_map_lores", rna_chemical_mapping_energy->calculate_energy( pose , true /*use_low_res*/ ) );
		}

		std::cout << "Outputting " << tag << " to silent file: " << silent_file << std::endl;
		silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );

		if ( option[ score::just_calc_rmsd ]() && tag.find( ".pdb" ) != std::string::npos ) {
			std::string out_pdb_file = utility::replace_in( tag, ".pdb", ".sup.pdb" );
			std::cout << "Creating: " << out_pdb_file << std::endl;
			pose.dump_pdb( out_pdb_file );
		}
		if ( option[ color_by_score ]() ) {
			do_color_by_score( pose );
			pose.dump_pdb( "COLOR_BY_SCORE.pdb" );
		}

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_score_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> " << std::endl;
		std::cout              << "              " << argv[0] << "  -in:file:silent <silent file> " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		utility::vector1< Size > blank_size_vector;
		utility::vector1< std::string > blank_string_vector;
		option.add_relevant( score::weights );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );
		option.add_relevant( in::file::tags );
		option.add_relevant( in::file::native );
		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::input_res );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( score::weights );
		option.add_relevant( score::just_calc_rmsd );
		option.add_relevant( OptionKeys::rna::data_file );
		NEW_OPT( original_input, "If you want to rescore the poses using the original FullModelInfo from a SWM run, input those original PDBs here", blank_string_vector );
		NEW_OPT( virtualize_free, "virtualize no-contact bases (and attached no-contact sugars/phosphates)", false );
		NEW_OPT( params_file, "Input file for pairings", "" );
		NEW_OPT( color_by_score, "color PDB by score (currently handles fa_atr & fa_rep)", false );
		NEW_OPT( soft_rep, "use soft_rep params for color_by_score", false ); // how about for actual scoring?

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
