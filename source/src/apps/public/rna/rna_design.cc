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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <protocols/rna/denovo/movers/RNA_Minimizer.hh>
#include <protocols/rna/denovo/options/RNA_MinimizerOptions.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/rna/denovo/util.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/excn/Exceptions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, disable_o2prime_rotamers )
OPT_KEY( Boolean, disable_include_current )
OPT_KEY( Boolean, sample_chi )
OPT_KEY( Boolean, ss_ds_ts_assign )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, all_RNA )
OPT_KEY( Boolean, rotamer_trials )
OPT_KEY( Boolean, final_minimize )
OPT_KEY( Boolean, multiround )
OPT_KEY( Real, ligand_distance )



///////////////////////////////////////////////////////////////////////////////
void
rna_sequence_recovery_metrics( pose::Pose const & reference_pose, utility::vector1< pose::PoseOP > const & pose_list, std::string const & sequence_recovery_file )
{

	// make a local copy
	pose::Pose pose( reference_pose );

	// Get information on ss, ds, ts residues in native pose.
	Size const nres = pose.size();
	FArray1D_int struct_type( nres, -1 );
	protocols::rna::denovo::check_base_pair( pose, struct_type );

	FArray1D_float recovery( nres, 0.0 );
	for ( Size n = 1; n <= pose_list.size(); n++ ) {
		for ( Size i = 1; i <= nres; i++ ) {
			// For noncanonicals
			if ( (pose_list[n])->residue(i).name3() == pose.residue(i).name3() ) {
				//if ( (pose_list[n])->residue(i).aa() == pose.residue(i).aa() ) {
				recovery( i ) += 1.0;
			}
		}
	}

	recovery /= pose_list.size();

	Size num_ss( 0 ), num_ds( 0 ), num_ts( 0 );
	Real frac_ss( 0.0 ), frac_ds( 0.0 ), frac_ts( 0.0 ), frac_overall( 0.0 );
	for ( Size i = 1; i <= nres; i++ ) {
		frac_overall += recovery(i);
		switch ( struct_type(i) ) {
		case 0 :
			num_ss++;
			frac_ss += recovery(i);
			break;
		case 1 :
			num_ds++;
			frac_ds += recovery(i);
			break;
		case 2 :
			num_ts++;
			frac_ts += recovery(i);
			break;
		}
	}

	if ( num_ss > 0.0 ) frac_ss /= num_ss;
	if ( num_ds > 0.0 ) frac_ds /= num_ds;
	if ( num_ts > 0.0 ) frac_ts /= num_ts;
	if ( nres > 0.0 ) frac_overall /= nres;

	std::map <Size, char > struct_symbol;
	struct_symbol[ 0 ] = 'S';
	struct_symbol[ 1 ] = 'D';
	struct_symbol[ 2 ] = 'T';

	utility::io::ozstream out( sequence_recovery_file );

	for ( Size i = 1; i <= nres; i++ ) {
		out << pose.residue(i).name1() << I(3,i) << " " << F(8,3,recovery(i)) << " " << struct_symbol[ struct_type(i) ] << std::endl;
	}

	out << std::endl;
	out << "SINGLE_STRANDED " << I(3,num_ss) << " " << F(9,4,frac_ss) << std::endl;
	out << "DOUBLE_STRANDED " << I(3,num_ds) << " " << F(9,4,frac_ds) << std::endl;
	out << "TERTIARY_STRUCT " << I(3,num_ts) << " " << F(9,4,frac_ts) << std::endl;
	out << "OVERALL         " << I(3,nres) << " " << F(9,4,frac_overall) << std::endl;

	out.close();

	std::cout << "Wrote stats to: " << sequence_recovery_file << std::endl;

}

bool residues_too_distant( conformation::Residue const & r1, conformation::Residue const & r2, Real const dist ) {
	return r1.xyz( r1.nbr_atom() ).distance_squared( r2.xyz( r2.nbr_atom() ) ) > dist*dist;
}

///////////////////////////////////////////////////////////////////////////////
void
rna_design_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace protocols::toolbox;
	using namespace protocols::rna::denovo;
	using namespace protocols::rna::denovo::movers;
	using namespace protocols::rna::denovo::options;


	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	core::import_pose::pose_from_file( pose, *rsd_set, pdb_file , core::import_pose::PDB_file);
	protocols::rna::denovo::ensure_phosphate_nomenclature_matches_mini( pose );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	RNA_MinimizerOptionsOP min_options( new RNA_MinimizerOptions );
	min_options->initialize_from_command_line();
	RNA_Minimizer rna_minimizer( min_options );

	dump_pdb( pose, "start.pdb");
	pose::Pose save_pose( pose );

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	utility::vector1< std::string > names;
	if ( option[ all_RNA ] ) {
		//auto const RNA_rsd_types = ResidueTypeFinder( *rsd_set ).base_property( RNA ).get_all_possible_residue_types();
		//for ( auto const & type : RNA_rsd_types ) { names.emplace_back( type->name3() ); }
		// Explicitly permit names of particular interest.
		names.emplace_back( "OMA" );
		names.emplace_back( "OMC" );
		names.emplace_back( "OMG" );
		names.emplace_back( "OMU" );
		names.emplace_back( "  I" );
		names.emplace_back( "H2U" );
		names.emplace_back( "PSU" );
		names.emplace_back( "1AP" ); // diaminopurine by another name.
		// Newly added for Dharmacon
		names.emplace_back( "OMI" );
		names.emplace_back( "PUR" );
		names.emplace_back( "2AP" );
		names.emplace_back( "8OG" );
		names.emplace_back( "6MG" );
		names.emplace_back( "5FC" );
		names.emplace_back( "5FU" );
		names.emplace_back( "5BU" );
		names.emplace_back( "5IU" );
		names.emplace_back( "NPU" );
	}
	//utility::vector1< core::chemical::AA > empty_na_vector;

	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		pack::task::parse_resfile(pose, *task);
	} else {
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			// If residue is far from ligand, skip
			if ( residues_too_distant( pose.residue( ii ), pose.residue( pose.size() ), option[ ligand_distance ]() ) ) continue;

			if ( !pose.residue_type( ii ).is_RNA() ) continue;
			//task->nonconst_residue_task( ii ).restrict_absent_nas( empty_na_vector );
			task->nonconst_residue_task( ii ).allow_aa( na_rad );
			task->nonconst_residue_task( ii ).allow_aa( na_ura );
			task->nonconst_residue_task( ii ).allow_aa( na_rgu );
			task->nonconst_residue_task( ii ).allow_aa( na_rcy );
			if ( option[ all_RNA ] ) {
				for ( auto const & name : names ) {
					task->nonconst_residue_task( ii ).allow_noncanonical_aa(name);
				}
			}
			assert( task->design_residue(ii) );
		}
	}

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( !pose.residue_type( ii ).is_RNA() ) continue;
		if ( residues_too_distant( pose.residue( ii ), pose.residue( pose.size() ), option[ ligand_distance ]() ) ) continue;

		//Hmmm, extras.
		//task->nonconst_residue_task( ii ).and_extrachi_cutoff( 0 );

		if ( option[ disable_o2prime_rotamers ]() ) task->nonconst_residue_task(ii).sample_proton_chi( false );

		//if ( option[ sample_chi ]() ) task->nonconst_residue_task(ii).nonconst_rna_task().set_sample_rna_chi( true );
	}

	ScoreFunctionOP scorefxn = get_score_function();

	// scorefxn->energy_method_options().exclude_DNA_DNA( exclude_DNA_DNA );
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );

	scorefxn->show( std::cout, pose );

	pose.dump_pdb( "start.pdb" );

	Size const nstruct = option[ out::nstruct ];
	utility::vector1< std::pair< Real, std::string > > results;
	utility::vector1< pose::PoseOP > pose_list;
	Size pos( pdb_file.find( ".pdb" ) );

	if ( option[ multiround ] ) {
		// TODO: also fill results vector
		core::pose::Pose start_pose = pose;
		for ( Size ii = 1; ii <= nstruct; ++ii ) {
			protocols::simple_moves::PackRotamersMover prm( scorefxn, task );
			prm.apply( pose );
			if ( option[ final_minimize ] ) {
				AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
				atom_level_domain_map->set( true );
				rna_minimizer.set_atom_level_domain_map( atom_level_domain_map );
				rna_minimizer.apply( pose );
			}
			pose_list.push_back( core::pose::PoseOP(&pose) );
			std::string outfile( pdb_file );
			std::stringstream suffix;
			suffix << ii << ".pdb";
			std::string name = outfile;
			name.replace( pos, 4, suffix.str() );
			if ( option[ dump ] ) pose.dump_pdb( name );
		}
	} else if ( option[ rotamer_trials ]() ) {
		// TODO: also fill results vector
		// AMW: this scheme allows the graphics viewer to follow each pose.
		core::pose::Pose start_pose = pose;
		for ( Size ii = 1; ii <= nstruct; ++ii ) {
			pose = start_pose;
			//core::pose::PoseOP poseop( new core::pose::Pose( pose ) );
			pack::rotamer_trials( pose, *scorefxn, task );

			if ( option[ final_minimize ] ) {
				AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( pose ) );
				atom_level_domain_map->set( true );
				rna_minimizer.set_atom_level_domain_map( atom_level_domain_map );
				rna_minimizer.apply( pose );
			}

			pose_list.push_back( core::pose::PoseOP(&pose) );
			std::string outfile( pdb_file );
			std::stringstream suffix;
			suffix << ii << ".pdb";
			std::string name = outfile;
			name.replace( pos, 4, suffix.str() );
			if ( option[ dump ] ) pose.dump_pdb( name );
		}
	} else {
		pack::pack_rotamers_loop( pose, *scorefxn, task, nstruct, results, pose_list);
		if ( option[ final_minimize ] ) {
			for ( auto const & poseop : pose_list ) {
				AtomLevelDomainMapOP atom_level_domain_map( new AtomLevelDomainMap( *poseop ) );
				atom_level_domain_map->set( true );
				rna_minimizer.set_atom_level_domain_map( atom_level_domain_map );
				rna_minimizer.apply( *poseop );
			}
		}
		std::string outfile( pdb_file );
		outfile.replace( pos, 4, ".pack.txt" );
		protocols::rna::denovo::export_packer_results( results, pose_list, scorefxn, outfile, option[ dump ] );
	}

	std::string sequence_recovery_file( pdb_file );
	sequence_recovery_file.replace( pos, 4, ".sequence_recovery.txt" );
	rna_sequence_recovery_metrics( save_pose, pose_list, sequence_recovery_file );

	// std::string const out_file_tag = "S_"+lead_zero_string_of( n, 4 );
	// dump_pdb( pose, out_file_tag + ".pdb" );

	pose = *( pose_list[1] );
	scorefxn->show( std::cout,pose );


}

///////////////////////////////////////////////////////////////
void
ss_ds_ts_assign_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::PoseOP pose_op( new pose::Pose );
	utility::vector1 < std::string > pdb_files  = option[ in::file::s ]();

	for ( Size n = 1; n <= pdb_files.size(); n++ )  {
		std::string const & pdb_file = pdb_files[ n ] ;
		core::import_pose::pose_from_file( *pose_op, *rsd_set, pdb_file , core::import_pose::PDB_file);
		protocols::rna::denovo::ensure_phosphate_nomenclature_matches_mini( *pose_op );

		std::string sequence_recovery_file( pdb_file );
		Size pos( pdb_file.find( ".pdb" ) );
		sequence_recovery_file.replace( pos, 4, ".ss_ds_ts.txt" );

		utility::vector1< pose::PoseOP > pose_list;
		pose_list.push_back( pose_op );
		rna_sequence_recovery_metrics( *pose_op, pose_list, sequence_recovery_file );
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ ss_ds_ts_assign ] ) {
		ss_ds_ts_assign_test();
	} else {
		rna_design_test();
		protocols::viewer::clear_conformation_viewers();
	}
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> [ -resfile <resfile>] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		//Uh, options? MOVE THESE TO OPTIONS NAMESPACE INSIDE CORE/OPTIONS.
		NEW_OPT( disable_o2prime_rotamers, "In designing, don't sample 2'-OH",false);
		NEW_OPT( disable_include_current, "In designing, don't include current",false);
		NEW_OPT( sample_chi,  "In designing RNA, chi torsion sample", false);
		NEW_OPT( ss_ds_ts_assign, "Figure out assignment of residues to single-stranded, double-stranded, tertiary contact categories",false);
		NEW_OPT( dump, "Dump pdb", false );
		NEW_OPT( all_RNA, "Use all RNA", false );
		NEW_OPT( rotamer_trials, "Do rotamer trials instead", false );
		NEW_OPT( final_minimize, "Do a final minimization (compare to minimized starting structures!", false );
		NEW_OPT( multiround, "Do a multiround protocol instead", false );
		NEW_OPT( ligand_distance, "Distance from ligand to design", 10.0 );


		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);


		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
