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
/// @brief Application to extract motifs from a set of pdb files

// libRosetta headers

#include <core/types.hh>

//#include <core/chemical/AtomTypeSet.hh>
//#include <core/chemical/MMAtomTypeSet.hh>

//#include <core/id/AtomID.hh>

//#include <core/chemical/AA.hh>
//#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
//#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/ResidueSelector.hh>
//#include <core/chemical/ChemicalManager.hh>
//#include <core/conformation/ResidueFactory.hh>
//#include <core/chemical/residue_io.hh>
//#include <core/chemical/VariantType.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

//#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreType.hh>
//#include <core/scoring/Ramachandran.hh>
//#include <core/pack/dunbrack/RotamerLibrary.hh>
//#include <core/scoring/hbonds/HBondSet.hh>
//#include <core/scoring/hbonds/hbonds.hh>
//#include <core/scoring/etable/count_pair/CountPairFunction.hh>
//#include <core/scoring/LREnergyContainer.hh>
//#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
//#include <core/scoring/TenANeighborGraph.hh>
//#include <core/scoring/methods/EnergyMethodOptions.hh>


//#include <core/pack/rotamer_trials.hh>
//#include <core/pack/pack_rotamers.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/rotamer_set/RotamerSet.hh>
//#include <core/pack/rotamer_set/RotamerSetFactory.hh>
//#include <core/pack/packer_neighbors.hh>

//#include <core/pack/rotamer_set/OptEData.hh>

//#include <core/graph/Graph.hh>

//#include <core/kinematics/FoldTree.hh>
//#include <core/kinematics/MoveMap.hh>
//#include <core/id/AtomID_Map.hh>
//#include <core/id/AtomID_Map.Pose.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>

//#include <core/scoring/mm/MMTorsionLibrary.hh>
//#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>

//#include <core/optimization/types.hh>
//#include <core/optimization/Multifunc.hh>
//#include <core/optimization/AtomTreeMinimizer.hh>
//#include <core/optimization/Minimizer.hh>
//#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

//#include <basic/database/open.hh>

//#include <utility/vector1.hh>
//#include <utility/exit.hh>
//#include <utility/pointer/owning_ptr.hh>
//#include <utility/pointer/ReferenceCount.hh>

//#include <numeric/xyzVector.hh>
//#include <numeric/random/random.hh>

//#include <ObjexxFCL/string.functions.hh>

//#include <protocols/dna/util.hh>
//#include <protocols/dna/classes.hh>

#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/SingleMotif.hh>
#include <protocols/motifs/MotifLibrary.hh>
//#include <protocols/motifs/motif_utils.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
//#include <iostream>
//#include <string>
//#include <algorithm>

//silly using/typedef


//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

//using namespace basic;
using namespace core;
using namespace pose;
using namespace chemical;
using namespace scoring;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;


// These are for aligning all interactions to a common coordinate frame.
// They don't need to be realistic

Vector const atom1_vector(   1.500,  0.000,  0.000 );
Vector const atom2_vector(   0.000,  0.000,  0.000 );
Vector const atom3_vector(   0.000,  1.500,  0.000 );

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "extract_motifs" );

void
fetch_atom_names(
	core::chemical::AA const & aa,
	std::string & oa1,
	std::string & oa2,
	std::string & oa3
)
{
	switch( aa ) {
	case aa_ala:
		oa1 = "1HB";
		oa2 = "CB";
		oa3 = "2HB";
		break;
	case aa_cys:
		oa1 = "SG";
		oa2 = "CB";
		oa3 = "CA";
		break;
	case aa_asp:
		oa1 = "OD1";
		oa2 = "CG";
		oa3 = "OD2";
		break;
	case aa_glu:
		oa1 = "OE1";
		oa2 = "CD";
		oa3 = "OE2";
		break;
	case aa_phe:
		oa1 = "CG";
		oa2 = "CD1";
		oa3 = "CZ";
		break;
	case aa_gly:
		oa1 = "HA1";
		oa2 = "CA";
		oa3 = "HA2";
		break;
	case aa_his:
		oa1 = "ND1";
		oa2 = "CE1";
		oa3 = "NE2";
		break;
	case aa_ile:
		oa1 = "CD1";
		oa2 = "CG1";
		oa3 = "CB";
		break;
	case aa_lys:
		oa1 = "NZ";
		oa2 = "CE";
		oa3 = "CD";
		break;
	case aa_leu:
		oa1 = "CD1";
		oa2 = "CG";
		oa3 = "CD2";
		break;
	case aa_met:
		oa1 = "CG";
		oa2 = "SD";
		oa3 = "CE";
		break;
	case aa_asn:
		oa1 = "OD1";
		oa2 = "CG";
		oa3 = "ND2";
		break;
	case aa_gln:
		oa1 = "OE1";
		oa2 = "CD";
		oa3 = "NE2";
		break;
	case aa_arg:
		oa1 = "NH1";
		oa2 = "CZ";
		oa3 = "NH2";
		break;
	case aa_ser:
		oa1 = "HG";
		oa2 = "OG";
		oa3 = "CB";
		break;
	case aa_thr:
		oa1 = "HG1";
		oa2 = "OG1";
		oa3 = "CB";
		break;
	case aa_val:
		oa1 = "CG1";
		oa2 = "CB";
		oa3 = "CG2";
		break;
	case aa_trp:
		oa1 = "CG";
		oa2 = "CE2";
		oa3 = "CZ3";
		break;
	case aa_tyr:
		oa1 = "CG";
		oa2 = "CD1";
		oa3 = "CZ";
		break;
	default:
		TR << "Did not find residue type for " << name_from_aa( aa ) << std::endl;
		oa1 = "X";
		oa2 = "X";
		oa3 = "X";
	}

	return;
}

void
output_single_motif(
	Pose & src_pose,
	AA const & target_aa,
	std::string & pdb_name,
	int prot_pos,
	std::vector< Size > &  contacts
)
{

	// Make a temporary pose
	Pose pose;

	// These need to be assigned be looking at the target amino acid
	std::string atom1_name;
	std::string atom2_name;
	std::string atom3_name;
	fetch_atom_names( target_aa, atom1_name, atom2_name, atom3_name );

	pose.append_residue_by_jump( src_pose.residue( prot_pos ), 1 );
	for( Size ic = 0 ; ic < contacts.size() ; ++ic ) {
		Size dna_pos = contacts[ ic ];
		pose.append_residue_by_jump( src_pose.residue( dna_pos ), 1 );
	}
	pose.conformation().insert_chain_ending( 1 );

	// Make a second pose _just for orientation_

	Pose centered_pose;
	centered_pose = pose;

	// Fix key residues
	centered_pose.set_xyz( id::AtomID( centered_pose.residue( 1 ).atom_index( atom1_name ), 1 ), atom1_vector );
	centered_pose.set_xyz( id::AtomID( centered_pose.residue( 1 ).atom_index( atom2_name ), 1 ), atom2_vector );
	centered_pose.set_xyz( id::AtomID( centered_pose.residue( 1 ).atom_index( atom3_name ), 1 ), atom3_vector );

	// Move the pose onto this fixed point

  core::kinematics::Stub end_stub(
    centered_pose.residue( 1 ).atom( atom2_name ).xyz(),
    centered_pose.residue( 1 ).atom( atom1_name ).xyz(),
    centered_pose.residue( 1 ).atom( atom2_name ).xyz(),
    centered_pose.residue( 1 ).atom( atom3_name ).xyz()
  );

  core::kinematics::Stub start_stub(
    pose.residue( 1 ).atom( atom2_name ).xyz(),
    pose.residue( 1 ).atom( atom1_name ).xyz(),
    pose.residue( 1 ).atom( atom2_name ).xyz(),
    pose.residue( 1 ).atom( atom3_name ).xyz()
  );

	for( Size ires = 1 ; ires <= pose.total_residue() ; ++ires ) {
		for( Size iatom = 1 ; iatom <= pose.residue(ires).natoms() ; ++iatom ) {
			pose.set_xyz( id::AtomID( iatom, ires ),
					end_stub.local2global( start_stub.global2local( pose.residue(ires).xyz( iatom ) ) ) );
		}
	}

	std::string const output_path( option[ out::path::pdb ]() );
	std::string delimiter( "_" );
	std::string extension( ".pdb" );
	std::string motif_file_name( output_path + src_pose.residue_type( prot_pos ).name3() +
			string_of( src_pose.pdb_info()->number( prot_pos ) ) +
			string_of( src_pose.pdb_info()->chain( prot_pos ) ) + delimiter );

	for( Size ic = 0 ; ic < contacts.size() ; ++ic ) {
		Size dna_pos = contacts[ ic ];
		motif_file_name = motif_file_name +
											src_pose.residue_type( dna_pos ).name1() +
											string_of( src_pose.pdb_info()->number( dna_pos ) ) +
											string_of( src_pose.pdb_info()->chain( dna_pos ) ) + delimiter;
	}

	//motif_file_name = motif_file_name + pdb_name + extension;
	motif_file_name = motif_file_name + pdb_name;

	TR << "Writing " << motif_file_name << std::endl;
	io::pdb::dump_pdb( pose, motif_file_name );

}

void
motif_distances
(
	protocols::motifs::Motif const & m1,
	protocols::motifs::Motif const & m2,
	core::Real & dist_diff,
	core::Real & angl_diff

)
{

	// Handle when motifs have different residue pairs
	if( m1.restype_name1() != m2.restype_name1() ||
			m1.restype_name2() != m2.restype_name2() ) {
		dist_diff = angl_diff = 9999.0;
		return;
	}

	core::kinematics::jump_distance( m1.forward_jump(), m2.forward_jump(), dist_diff, angl_diff );

	return;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

Real
get_packing_score(
	Pose & pose,
	Size pos1,
	Size pos2,
	ScoreFunction & sf
)
{

	EnergyMap pack_map;

	sf.eval_ci_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, pack_map );

	return pack_map[ fa_atr ] + pack_map[ fa_rep ];
}

Real
get_hbond_score(
	Pose & pose,
	Size pos1,
	Size pos2,
	ScoreFunction & sf
)
{

	EnergyMap hbond_map;

	sf.eval_cd_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, hbond_map );

	return hbond_map[ hbond_sc ];
}

Real
get_elec_score(
	Pose & pose,
	Size pos1,
	Size pos2,
	ScoreFunction & sf
)
{

	EnergyMap elec_map;

	sf.eval_cd_2b_sc_sc( pose.residue( pos1 ), pose.residue( pos2 ), pose, elec_map );

	return elec_map[ fa_elec ];

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

void
process_for_motifs(
	Pose & pose,
	std::string & pdb_name,
	chemical::AA const target_aa,
	protocols::motifs::MotifLibrary & motif_lib
)
{
	int nres( pose.total_residue() );

	Real dist_threshold( option[ motifs::duplicate_dist_cutoff ]  );
	Real angl_threshold( option[ motifs::duplicate_angle_cutoff ] );

	Real pack_min_threshold( option[ motifs::pack_min_threshold ] );
	Real pack_max_threshold( option[ motifs::pack_max_threshold ] );
	Real hbond_min_threshold( option[ motifs::hbond_min_threshold ] );
	Real hbond_max_threshold( option[ motifs::hbond_max_threshold ] );
	Real elec_min_threshold( option[ motifs::elec_min_threshold ] );
	Real elec_max_threshold( option[ motifs::elec_max_threshold ] );

	// Get a simple scorefunction to screen for thresholds
	ScoreFunction scorefxn;
	// These are talaris2013 values.  Since the thresholds can be
	// chosen arbitrarily, the only thing that is truly hard-wired
	// here is the ratio between fa_atr and fa_rep, which are
	// combined in get_packing_score above.
  scorefxn.set_weight( fa_atr, 0.8 );
  scorefxn.set_weight( fa_rep, 0.44 );
  scorefxn.set_weight( hbond_sc, 1.1 );
  scorefxn.set_weight( fa_elec, 0.7 );
	scorefxn.setup_for_scoring( pose );

	// These need to be assigned be looking at the target amino acid
	std::string atom1_name;
	std::string atom2_name;
	std::string atom3_name;
	fetch_atom_names( target_aa, atom1_name, atom2_name, atom3_name );

	// Loop over positions, skipping things that aren't the target amino acid
	for( int aaoi_pos = 1 ; aaoi_pos <= nres ; ++aaoi_pos ) {
		ResidueType const & res_type( pose.residue_type( aaoi_pos ) );

		if( res_type.aa() != target_aa ) continue;

		std::vector< Size > contacts;

		Real total_pack_score( 0.0 );
		Real total_hb_score( 0.0 );
		Real check_hb_score( 0.0 );

		// Loop over all other positions
		for( int other_pos = 1 ; other_pos <= nres ; ++other_pos ) {
			ResidueType const & other_type( pose.residue_type( other_pos ) );
			if( other_pos == aaoi_pos ) continue;

			// This prevents against getting the same interaction twice
			if( ( other_type.aa() == target_aa ) && ( other_pos < aaoi_pos ) ) continue;

			Real pack_score = get_packing_score( pose, other_pos, aaoi_pos, scorefxn );
			Real hb_score = get_hbond_score( pose, other_pos, aaoi_pos, scorefxn );
			Real elec_score = get_elec_score( pose, other_pos, aaoi_pos, scorefxn );

			int num_contacts( 0 );

			// Loop over all amino acid atoms
			for( Size this_atom = 1 ; this_atom <= res_type.nheavyatoms() ; ++this_atom ) {
				// Loop over heavy atoms
				for( Size other_atom = other_type.first_sidechain_atom() ;
						other_atom <= other_type.nheavyatoms() ; ++other_atom ) {

					Real const dis2( pose.residue( other_pos ).xyz( other_atom ).distance_squared(
																	pose.residue( aaoi_pos ).xyz( this_atom ) ) );
					if( dis2 < (3.9*3.9) ) {
						++num_contacts;
					}

				}
			} // End loop over atoms

			check_hb_score += hb_score;

//			if( num_contacts > 0 )
//			if( pack_score <= -0.5 && hb_score > -0.2 )
//			if( hb_score <= -1.00 ) {
//			if( pack_score <= -1.00 ) {

			if( ( pack_score > pack_min_threshold && pack_score < pack_max_threshold) &&
					( hb_score > hbond_min_threshold && hb_score < hbond_max_threshold) &&
					( elec_score > elec_min_threshold && elec_score < elec_max_threshold) ) {

				TR << "Energies between " << aaoi_pos << " and " << other_pos << " are pack: " << pack_score << " hbond: " << hb_score << std::endl;

				total_pack_score += pack_score;
				total_hb_score += hb_score;

				// Store the motif info
				std::string other_atom1;
				std::string other_atom2;
				std::string other_atom3;
				fetch_atom_names( pose.aa( other_pos ), other_atom1, other_atom2, other_atom3 );
				if( other_atom1 != "X" ) {
					contacts.push_back( other_pos );
					protocols::motifs::SingleMotif new_motif( pose, aaoi_pos, atom1_name, atom2_name, atom3_name, other_pos, other_atom1, other_atom2, other_atom3 );
					std::string motif_origin( pose.residue_type( aaoi_pos ).name3() + " " + string_of( pose.pdb_info()->number( aaoi_pos ) ) + string_of( pose.pdb_info()->chain( aaoi_pos ) ) + " and " + pose.residue_type( other_pos ).name3() + " " + string_of( pose.pdb_info()->number( other_pos ) ) + string_of( pose.pdb_info()->chain( other_pos ) ) + " pdb code:  " + pdb_name );
					new_motif.store_remark( motif_origin );
					TR << new_motif.remark() << std::endl;

					protocols::motifs::Motif work_motif( new_motif );

					bool unique_motif( true );

					// Check against other motifs
					core::Size check_i( 1 );
					for( protocols::motifs::MotifCOPs::const_iterator this_motif = motif_lib.begin() ; this_motif != motif_lib.end() ; ++this_motif ) {
						core::Real dist_diff( 0.0 );
						core::Real angl_diff( 0.0 );
						motif_distances( work_motif, **this_motif, dist_diff, angl_diff );
//						TR << "Motif dist from current motif " << check_i << " is translation " << dist_diff << " and angle " << angl_diff << std::endl;
						check_i++;
						if( dist_diff < dist_threshold && angl_diff < angl_threshold ) {
							unique_motif = false;
							break;
						}
					}

					if( unique_motif ) {
						motif_lib.add_to_library( new_motif );
					}
				}
			}

		} // End loop over amino acids

		if( contacts.size() > 0 ) {
			TR << "Outputing motif with " << contacts.size() << " contacts " << std::endl;
			TR << "Energies are pack: " << total_pack_score << " hbond: " << total_hb_score << std::endl;
			TR << "Check HB Energy is: " << check_hb_score << std::endl;
			output_single_motif( pose, target_aa, pdb_name, aaoi_pos, contacts );
		}


	} // End loop over residues

	return;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	using namespace pose;
	using namespace conformation;
	using namespace chemical;
	using namespace import_pose;

	//using namespace core;
	devel::init( argc, argv );

	// Get the name of the file that contains the pdb files for extraction
	std::string const list_of_pdb_files( option[ in::file::l ]()[1] );
	std::ifstream file_data( list_of_pdb_files.c_str() );

	// find the amino acid for which interaction motifs are requested
	std::string target_3code( option[ motifs::target_aa ]() );
	chemical::AA target_aa( aa_from_name( uppercase( target_3code ) ) );

	std::string pdb_code;

	protocols::motifs::MotifLibrary motif_lib;

	std::string const output_path( option[ out::path::pdb ]() );
	std::string motif_path_and_file( output_path + "motif_library" );
	std::ofstream motif_ostream( motif_path_and_file.c_str() );

	file_data >> pdb_code;
	while( !file_data.eof() ) {

		TR << "Working on pdb " << pdb_code << std::endl;
		std::string const pdb_file_path( option[ in::path::pdb ]( 1 ) );
		std::string pdb_file( pdb_file_path + pdb_code );
		TR << "PDB file " << pdb_file << std::endl;

		Pose pose;
		pose_from_pdb( pose, pdb_file );

		process_for_motifs( pose, pdb_code, target_aa, motif_lib );

		file_data >> pdb_code;
	}

	motif_ostream << motif_lib;
	motif_ostream.close();

	return 0;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception" << e.msg() << std::endl;
		return -1;
	}
}

///////////////////////////////////////////////////////////////////////////////
