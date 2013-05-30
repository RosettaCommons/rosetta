// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief app for collection of protein-DNA interaction motifs
/// @author sthyme
///

#include <devel/init.hh>

// Project Headers (protocols)
#include <protocols/dna/util.hh>
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifLibrary.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
static basic::Tracer TR("apps.pilot.dna_motifs_collector");

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <map>
#include <string>

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


///////////////////////////////////////////////////////////////////////////////

void
output_single_motif(
	core::pose::Pose & src_pose,
	std::string & pdb_prefix,
	core::Size prot_pos,
	utility::vector1< core::Size > &  contacts,
	utility::io::ozstream & motif_output_file
)
{

	bool keep_motif_xtal_location = basic::options::option[ basic::options::OptionKeys::motifs::keep_motif_xtal_location ](); // Default false
	core::Size prot_pos_pdb( src_pose.pdb_info()->number( prot_pos ) );
	char prot_pos_chain( src_pose.pdb_info()->chain( prot_pos ) );
	std::string delimiter( "_" );

	for( core::Size ic(1) ; ic <= contacts.size() ; ++ic ) {
		protocols::motifs::Motif motif( src_pose.residue( prot_pos ), src_pose.residue( contacts[ic] ) );
		core::conformation::ResidueOP refdnares = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
		core::conformation::Residue protres( src_pose.residue( prot_pos ) );
		core::conformation::Residue dnares( src_pose.residue( contacts[ic] ) );

		core::Size dna_pos_pdb( src_pose.pdb_info()->number( contacts[ic] ) );
		char dna_pos_chain( src_pose.pdb_info()->chain( contacts[ic] ) );
		std::string motif_name = src_pose.residue_type( prot_pos ).name3() + ObjexxFCL::string_of( prot_pos_pdb ) + ObjexxFCL::string_of( prot_pos_chain) + delimiter +  src_pose.residue_type( contacts[ic] ).name1() + ObjexxFCL::string_of( dna_pos_pdb ) + ObjexxFCL::string_of( dna_pos_chain ) + delimiter + pdb_prefix;
		motif.store_remark( motif_name );

		if( ! keep_motif_xtal_location ) {
			motif.place_residue( *refdnares, protres );
			motif.place_residue( protres, dnares );
		}
		motif_output_file << motif;
		if( basic::options::option[ basic::options::OptionKeys::motifs::motif_output_directory ].user() ) {
			std::string output_path( basic::options::option[ basic::options::OptionKeys::motifs::motif_output_directory ]() ); // No default
			std::string extension( ".pdb" );
			std::string motif_file_name( output_path + src_pose.residue_type( prot_pos ).name3() + ObjexxFCL::string_of( prot_pos_pdb ) + ObjexxFCL::string_of( prot_pos_chain ) + delimiter );
			core::pose::Pose pose;
			pose.append_residue_by_jump( protres, 1 );
			core::Size dna_pos = contacts[ ic ];
			pose.append_residue_by_jump( dnares, 1 );
			core::Size dna_pos_pdb( src_pose.pdb_info()->number( dna_pos ) );
			char dna_pos_chain( src_pose.pdb_info()->chain( dna_pos ) );
			std::string motif_file_name_final = motif_file_name + src_pose.residue_type( dna_pos ).name1() + ObjexxFCL::string_of( dna_pos_pdb ) + ObjexxFCL::string_of( dna_pos_chain ) + delimiter + pdb_prefix + extension;
			TR << "Writing " << motif_file_name_final << std::endl;
			core::io::pdb::dump_pdb( pose, motif_file_name_final );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

void
place_waters_and_minimize(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction & scorefxn,
	std::string & pdb_prefix
)
{
	bool place_adduct_waters = basic::options::option[ basic::options::OptionKeys::motifs::place_adduct_waters ](); // Default true
	bool preminimize_motif_pdbs = basic::options::option[ basic::options::OptionKeys::motifs::preminimize_motif_pdbs ](); // Default false
	bool preminimize_motif_pdbs_sconly = basic::options::option[ basic::options::OptionKeys::motifs::preminimize_motif_pdbs_sconly ](); // Default false
	core::Energy score_orig = scorefxn( pose );
	TR << "Starting pose score " << score_orig << std::endl;

	// The packing is strictly to allow waters to be placed at canonical/adduct positions in the pose
	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	core::kinematics::MoveMap mm;
  mm.set_bb( false );
  mm.set_chi( false );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i) {
		if ( pose.residue(i).is_protein() ) {
			task->nonconst_residue_task(i).prevent_repacking();
			// Only allowing protein minimization currently
			if( preminimize_motif_pdbs ) {
				mm.set_chi( i, true );
				mm.set_bb( i, true );
			}
			if( preminimize_motif_pdbs_sconly) {
				mm.set_chi( i, true );
			}
		}
	}

	if( place_adduct_waters ) {
		task->set_bump_check( true );
		core::pack::pack_rotamers( pose, scorefxn, task);
		core::Energy end_score = scorefxn( pose );
		TR << "Score after addition of waters " << end_score << std::endl;
	}

	if( preminimize_motif_pdbs || preminimize_motif_pdbs_sconly ) {

		core::Real min_tolerance = basic::options::option[ basic::options::OptionKeys::run::min_tolerance ]();
		std::string min_type("dfpmin_armijo_nonmonotone");
		if ( basic::options::option[ basic::options::OptionKeys::run::min_type ].user() ) min_type = basic::options::option[ basic::options::OptionKeys::run::min_type ]();
		core::optimization::MinimizerOptions options( min_type, min_tolerance, true, false );

		core::optimization::AtomTreeMinimizer minimizer;
		minimizer.run( pose, mm, scorefxn, options );
		std::string minimized_pdb = "post_minimization_" + pdb_prefix + ".pdb";
		if ( preminimize_motif_pdbs_sconly) {
			minimized_pdb = "post_minimization_sconly_" + pdb_prefix + ".pdb";
		}
		core::io::pdb::dump_pdb( pose, minimized_pdb );

		core::Energy min_score = scorefxn( pose );
		TR << "Score after minimization " << min_score << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////

core::Real
get_packing_score(
	core::pose::Pose & pose,
	core::Size pos1,
	core::Size pos2,
	core::scoring::ScoreFunction & scorefxn
)
{
	core::scoring::EnergyMap pack_map;
	scorefxn.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, pack_map );
	return pack_map[ core::scoring::fa_atr ];
}

core::Real
get_hbond_score(
	core::pose::Pose & pose,
	core::Size pos1,
	core::Size pos2,
	core::scoring::ScoreFunction & scorefxn
)
{
	core::scoring::EnergyMap hbond_map;
	scorefxn.eval_cd_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, hbond_map );
	return hbond_map[ core::scoring::hbond_sc ];
}

core::Real
get_water_hbond_score(
	core::pose::Pose & pose,
	core::Size pos1,
	core::Size pos2,
	core::scoring::ScoreFunction & scorefxn
)
{
	core::scoring::EnergyMap water_hbond_map;
	scorefxn.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, water_hbond_map );
	return water_hbond_map[ core::scoring::h2o_hbond ];
}

///////////////////////////////////////////////////////////////////////////////

void
process_for_motifs(
	core::pose::Pose & pose,
	std::string & pdb_prefix,
	protocols::motifs::MotifLibrary & motifs,
	utility::io::ozstream & motif_output_file
)
{

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );

	core::Size nres( pose.total_residue() );

	// Potentially place waters and/or minimize
	place_waters_and_minimize( pose, *scorefxn, pdb_prefix );

	// Loop over positions, skipping non-amino acid
	for( core::Size prot_pos = 1 ; prot_pos <= nres ; ++prot_pos ) {
		core::chemical::ResidueType const & res_type( pose.residue_type( prot_pos ) );
		if( ! res_type.is_protein() ) continue;

		// Map will automatically sort the "contacts" with the lowest total_score at the front of map
		std::map< core::Real, core::Size > contacts;

		// Loop over positions, skipping non-DNA
		for( core::Size dna_pos = 1 ; dna_pos <= nres ; ++dna_pos ) {
			core::chemical::ResidueType const & res_type( pose.residue_type( dna_pos ) );
			if( ! res_type.is_DNA() ) continue;

			// Get all the cutoffs that residue interactions have to pass to be counted as a motif
			core::Real pack_score = get_packing_score( pose, dna_pos, prot_pos, *scorefxn );
			core::Real pack_score_cutoff = basic::options::option[ basic::options::OptionKeys::motifs::pack_score_cutoff ](); // Default -0.5
			core::Real hb_score = get_hbond_score( pose, dna_pos, prot_pos, *scorefxn );
			core::Real hb_score_cutoff = basic::options::option[ basic::options::OptionKeys::motifs::hb_score_cutoff ](); // Default -0.3
			core::Real water_score = get_water_hbond_score( pose, dna_pos, prot_pos, *scorefxn );
			core::Real water_score_cutoff = basic::options::option[ basic::options::OptionKeys::motifs::water_score_cutoff ](); // Default -0.3

			core::Real total_score = pack_score + hb_score + water_score;
			if( pack_score < pack_score_cutoff || hb_score < hb_score_cutoff || water_score < water_score_cutoff ) {
				contacts[total_score] = dna_pos;
				core::Size prot_pos_pdb( pose.pdb_info()->number( prot_pos ) );
				char prot_pos_chain( pose.pdb_info()->chain( prot_pos ) );
				core::Size dna_pos_pdb( pose.pdb_info()->number( dna_pos ) );
				char dna_pos_chain( pose.pdb_info()->chain( dna_pos ) );
				TR << "Energies between " << prot_pos_pdb << prot_pos_chain << " and " << dna_pos_pdb << dna_pos_chain <<" are total: " << total_score << " pack: " << pack_score << " hbond: " << hb_score << " water: " << water_score << std::endl;
			}
		} // End loop over DNA residues

		core::Size contacts_size( contacts.size() );
		bool eliminate_weak_motifs = basic::options::option[ basic::options::OptionKeys::motifs::eliminate_weak_motifs ](); // Default true
		if( contacts_size != 0 ) {
			utility::vector1< core::Size > contacts_strength_filter;
			core::Size contacts_size( contacts.size() );
			for( std::map< core::Real, core::Size >::const_iterator it( contacts.begin() ),
					end( contacts.end() ); it != end; ++it ) {
				if( contacts_size == 1 ) {
					contacts_strength_filter.push_back( it->second );
					break;
				}
				if( eliminate_weak_motifs) {
					core::Real first( (it)->first );
					core::Real firstb( (++it)->first );
					core::Real divided( first / firstb );
					if( divided > 2.0 ) {
						contacts_strength_filter.push_back( (--it)->second );
						break;
					} else {
						contacts_strength_filter.push_back( it->second );
						contacts_strength_filter.push_back( (--it)->second );
						break;
					}
				}
				TR << "Adding all motif contacts that pass the filters because elimination of weak contacts is turned off. Warning, this is a lot of motifs." << std::endl;
				contacts_strength_filter.push_back( it->second );
			} // End loop over the inital map of potential contacts, getting rid of weaker contacts by default

			utility::vector1< core::Size > final_contacts;

			core::Size prot_pos_pdb( pose.pdb_info()->number( prot_pos ) );
			char prot_pos_chain( pose.pdb_info()->chain( prot_pos ) );
			std::string delimiter( "_" );
			core::Real duplicate_motif_cutoff = basic::options::option[ basic::options::OptionKeys::motifs::duplicate_motif_cutoff ](); // Default 0.2
			for( core::Size ic=1 ; ic <= contacts_strength_filter.size() ; ++ic ) {
				if( pose.residue( prot_pos ).name3() == "GLY" ) continue; // GLY cannot currently be a motif
				protocols::motifs::Motif motif( pose.residue( prot_pos ), pose.residue( contacts_strength_filter[ic] ) );

				// Making canonical Rosetta residues for placing as motifs and determining if motif already exists in library
				// You need to use canonical residues in order to do an RMSD test, unless your test was only over the 6 motif atoms (and you used place_atoms instead of place_residue)
				// The only reason to change this code would be if you really wanted a speed increase in the motif finding
				core::conformation::ResidueOP protres = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( motif.restype_name1() ) );
				core::conformation::ResidueOP dnares = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
				core::conformation::ResidueOP dnares2 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motif.restype_name2()) ) );
				motif.place_residue( *protres, *dnares );
				bool broke( false );

				core::Size dna_pos = contacts_strength_filter[ ic ];
				core::Size dna_pos_pdb( pose.pdb_info()->number( dna_pos ) );
				char dna_pos_chain( pose.pdb_info()->chain( dna_pos ) );
				std::string motif_name = pose.residue_type( prot_pos ).name3() + ObjexxFCL::string_of( prot_pos_pdb ) + ObjexxFCL::string_of( prot_pos_chain) + delimiter +  pose.residue_type( dna_pos ).name1() + ObjexxFCL::string_of( dna_pos_pdb ) + ObjexxFCL::string_of( dna_pos_chain ) + delimiter + pdb_prefix;
				motif.store_remark( motif_name );

				for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motifs.begin(), end_itr = motifs.end();
						motifcop_itr != end_itr; ++motifcop_itr ) {
					protocols::motifs::MotifCOP motifcop( *motifcop_itr );
					if( motifcop->restype_name1() != motif.restype_name1() ) continue;
					if( motifcop->restype_name2() != motif.restype_name2() ) continue;
					motifcop->place_residue( *protres, *dnares2 );
					core::Real rmsdtest = core::scoring::automorphic_rmsd( *dnares, *dnares2, false );
					if( rmsdtest < duplicate_motif_cutoff ) {
						TR << "Skipping motif " << motif.remark() << " because it matches motif " << motifcop->remark() << " already in motif library, with an RMSD = " << rmsdtest << std::endl;
						broke = true;
						break;
					}
				}
				if( !broke ) {
					TR << "Adding motif " << motif.remark() << std::endl;
					final_contacts.push_back( contacts_strength_filter[ic] );
					motifs.add_to_library( motif );
				}
			}  // End loop that adds motifs to the MotifLibrary and skips the motif if a very similar motif has already been found
			if( final_contacts.size() >= 1 ) {
				output_single_motif( pose, pdb_prefix, prot_pos, final_contacts, motif_output_file );
			}
		} // if there were actually contacts found for the particular protein residue
	} // End loop over protein residues

}

///////////////////////////////////////////////////////////////////////////////

void
process_file_list()
{

	utility::vector1< std::string > pdb_files( basic::options::start_files() );
	protocols::motifs::MotifLibrary motifs;
	protocols::motifs::MotifLibrary motif_library_previous( protocols::motifs::get_MotifLibrary_user() );

	std::string MotifLibraryFileName = "MotifLibrary.motifs";
	if ( basic::options::option[ basic::options::OptionKeys::motifs::output_file ].user() ) {
		MotifLibraryFileName = basic::options::option[ basic::options::OptionKeys::motifs::output_file ]();
	}
	utility::io::ozstream motif_output_file( MotifLibraryFileName );

	if ( motif_library_previous.nmotifs() > 0 ) {
		TR << "Adding motifs from a previously made library" << std::endl;
		for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motif_library_previous.begin(), end_itr = motif_library_previous.end();
				motifcop_itr != end_itr; ++motifcop_itr ) {
				protocols::motifs::MotifCOP motifcop( *motifcop_itr );
				motifs.add_to_library( *motifcop );
				motif_output_file << *motifcop;
		}
	}

	for ( utility::vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
	      pdb_file != pdb_files.end(); ++pdb_file ) {
		std::string pdb_name( *pdb_file );

		TR << "Working on file: " << pdb_name << std::endl;

		std::string pdb_prefix( utility::string_split( utility::string_split( pdb_name, '/' ).back(), '.' ).front() );

		core::pose::PoseOP pose = new core::pose::Pose;
		core::import_pose::pose_from_pdb( *pose, pdb_name );

		process_for_motifs( *pose, pdb_prefix, motifs, motif_output_file );
	}

}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{

	try {

	devel::init( argc, argv );
	process_file_list();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
