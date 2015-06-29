// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/motifs/LigandMotifSearch.cc
/// @brief Ligand motif searching protocol
/// @author sthyme (sthyme@gmail.com); mdsmith (mdwsmith@u.washington.edu)

// Unit Headers
#include <protocols/motifs/LigandMotifSearch.hh>

// Package Headers
#include <protocols/motifs/BuildPosition.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifHit.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/dna/util.hh>
#include <protocols/simple_moves/MinMover.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh> //Need this to prevent the compiling error: invalid use of incomplete type 'const struct core::chemical::AtomType
#include <core/chemical/AtomTypeSet.hh> //Need this to get AtomType integers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/chemical/VariantType.hh>

#include <protocols/toolbox/rotamer_set_operations/SpecialRotamerRotSetOps.hh>
#include <core/pose/util.hh>

// Utility Headers
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyzVector.hh>

// C++ Headers
#include <iostream>

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace motifs {

static thread_local basic::Tracer ms_tr( "protocols.motifs.LigandMotifSearch", basic::t_info );

LigandMotifSearch::~LigandMotifSearch()
{}

LigandMotifSearch::LigandMotifSearch()
	: motif_library_(),
		target_positions_(),
		build_positionOPs_(/* 0 */),
		target_conformers_map_(),
		// Option flags/parameters: default to command line options
		ztest_cutoff_1_(	basic::options::option[ basic::options::OptionKeys::motifs::z1 ]() ),
		ztest_cutoff_2_(	basic::options::option[ basic::options::OptionKeys::motifs::z2 ]() ),
		rmsd_cutoff_1_(	basic::options::option[ basic::options::OptionKeys::motifs::r1 ]() ),
		rmsd_cutoff_2_(	basic::options::option[ basic::options::OptionKeys::motifs::r2 ]() ),
		dtest_cutoff_(	basic::options::option[ basic::options::OptionKeys::motifs::dtest ]() ),
		rot_level_(	basic::options::option[ basic::options::OptionKeys::motifs::rotlevel ]() ),
		minimize_(	basic::options::option[ basic::options::OptionKeys::motifs::minimize ]() )
{
	init_options();
}

LigandMotifSearch::LigandMotifSearch( LigandMotifSearch const & src ) :
	utility::pointer::ReferenceCount( src )
{
	(*this) = src;
}

LigandMotifSearch const &
LigandMotifSearch::operator = ( LigandMotifSearch const & src )
{
	if( this != &src ) {
		motif_library_ = src.motif_library_ ;
		target_positions_ = src.target_positions_ ;
		build_positionOPs_ = src.build_positionOPs_ ;
		target_conformers_map_ = src.target_conformers_map_ ;
		ztest_cutoff_1_ = src.ztest_cutoff_1_ ;
		ztest_cutoff_2_ = src.ztest_cutoff_2_ ;
		rmsd_cutoff_1_ = src.rmsd_cutoff_1_ ;
		rmsd_cutoff_2_ = src.rmsd_cutoff_2_ ;
		dtest_cutoff_ = src.dtest_cutoff_ ;
		rot_level_ = src.rot_level_ ;
		minimize_ = src.minimize_;
		bpdata_ = src.bpdata_;
		bpdata_filename_ = src.bpdata_filename_;
		output_ = src.output_;
		output_filename_ = src.output_filename_;
		data_ = src.data_;
		data_filename_ = src.data_filename_;
		quick_and_dirty_ = src.quick_and_dirty_;
		dump_motifs_ = src.dump_motifs_;
		clear_bprots_ = src.clear_bprots_;
		rots2add_ = src.rots2add_;
	}
	return *this;
}

void
LigandMotifSearch::run(
	Pose const & pose,
	utility::vector1< Size > & input_BPs
)
{
	for (  utility::vector1< Size  >::const_iterator position( input_BPs.begin() );
				 position != input_BPs.end(); ++position ) {
		ms_tr << "In run, Design position: " << *position << std::endl;
	}
	initialize( pose, input_BPs );
	incorporate_motifs( pose );
}

void
LigandMotifSearch::run(
	Pose const & pose,
	core::Real & ligand_motif_sphere
)
{
	using utility::file::file_exists;

	utility::io::ozstream motif_output_file;

	utility::vector1< Size > target_positions_sphere ;
	target_positions_sphere = get_sphere_aa( pose, ligand_motif_sphere ) ;
	initialize( pose, target_positions_sphere );

//Move functions back down to after open_append if we have problems

	if( output_ ) {
		if ( file_exists( output_filename_  ) ) {
			Pose & pose2 = const_cast< Pose & >( pose );
			//Size const ligand_marker = 0;
			ms_tr << "Motif search has already run, we will get rotamers from output file" << std::endl;
			//We already have the output file, let's just take the rotamers from there rather than searching again.
			for( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
					 ir != end_ir; ++ir ) {
				load_build_position_data( **ir, output_filename_, pose2, LIGAND );
				if( clear_bprots_ ) {
					// I could avoid adding them in the first place, but that loading function is more complicated than just clearing and I might want them later
					(*ir)->clear_rots();
				}
			}
			if( ! motif_library_.empty() ) {
				ms_tr << "Loading BPData, but also loaded a MotifLibrary of all motifs.  Probably not a problem for ligand motifs (we will always hit this condition)" << std::endl;
			}

		} //If the output file already exists
		else { //If the output file doesn't already exist, let's find some motifs!
			ms_tr << "Starting motif search" << std::endl;
			motif_output_file.open_append( output_filename_ );


			incorporate_motifs( pose );
		}

	} else {
		//User gave no output file, this is a problem!
		ms_tr << "You didn't give a -motifs:output filename, but we need that for storage! Please specify a filename." << std::endl;

	}
	//This is run for Enzdes, we have no input_BPs
}


//run as task operation
void
LigandMotifSearch::run(
	Pose const & pose,
	PackerTask & task
)
{
	using utility::file::file_exists;

	utility::io::ozstream motif_output_file;

	utility::vector1< core::Size > target_positions;

	for(core::Size i=1; i<=pose.total_residue(); ++i) {
		if (pose.residue(i).is_protein() && task.being_designed(i) ) {
			target_positions.push_back(i);
		}
	}

	initialize( pose, target_positions );

	if( output_ ) {
		if ( file_exists( output_filename_  ) ) {
			Pose & pose2 = const_cast< Pose & >( pose );
			//Size const ligand_marker = 0;
			ms_tr << "Motif search has already run, we will get rotamers from output file" << std::endl;
			//We already have the output file, let's just take the rotamers from there rather than searching again.
			for( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
					 ir != end_ir; ++ir ) {
				load_build_position_data( **ir, output_filename_, pose2, LIGAND );
				if( clear_bprots_ ) {
					// I could avoid adding them in the first place, but that loading function is more complicated than just clearing and I might want them later
					(*ir)->clear_rots();
				}
			}
			if( ! motif_library_.empty() ) {
				ms_tr << "Loading BPData, but also loaded a MotifLibrary of all motifs.  Probably not a problem for ligand motifs (we will always hit this condition)" << std::endl;
			}

		} //If the output file already exists
		else { //If the output file doesn't already exist, let's find some motifs!
			ms_tr << "Starting motif search" << std::endl;
			motif_output_file.open_append( output_filename_ );


			incorporate_motifs( pose );
		}
		//OK, we have either loaded motifs from file or we have found them.  Now we will add rotamers to task.


		std::set< core::Size > src_pos;


		for(core::Size i=1; i<=pose.total_residue(); ++i) {
			if (pose.residue(i).is_protein() && task.being_designed(i) ) {

				std::set< std::string > name3set;
				core::pack::rotamer_set::Rotamers motif_rotamers( bp_rotamers( i ) );
				src_pos.insert( i );
				core::pack::rotamer_set::Rotamers variant_rotamers;

				if ( ! motif_rotamers.empty() ) {
					ms_tr << "yes, we have motif rotamers for this position!" << std::endl;
					for( core::Size rot(1); rot <= motif_rotamers.size(); ++rot ) {
						core::conformation::ResidueOP variant_rot( core::pose::add_variant_type_to_residue( *(motif_rotamers[rot]), core::chemical::SPECIAL_ROT, pose ) ); //This is where we're crashing (original 546)
						variant_rotamers.push_back( variant_rot );

						name3set.insert( (motif_rotamers[rot])->name3() );
					}
					//rotamer_map[i] = variant_rotamers;
					protocols::toolbox::rotamer_set_operations::SpecialRotamerRSOOP ms_rsoop( new protocols::toolbox::rotamer_set_operations::SpecialRotamerRSO( i ) );
					ms_rsoop->set_new_rots( variant_rotamers );
					//		taskfactory->push_back( new AppendRotamerSet( ms_rsoop ) );
					//task.nonconst_residue_task( lig_it->first ).append_rotamerset_operation( rb_rotsetop );
					task.nonconst_residue_task( i ).append_rotamerset_operation(  ms_rsoop  );
				} else {
					ms_tr << "no motif rotamers for this position" << std::endl;
				} //end else
			} //end if position is protein and designed
		} //end protein position loop
///////////////////////////////// END ENZDES TASK OPERATION ///////////////////////////////////////////////

	} else {
		//User gave no output file, this is a problem!
		ms_tr << "You didn't give a -motifs:output filename, but we need that for storage! Please specify a filename." << std::endl;

	}
}

void
LigandMotifSearch::initialize(
	Pose const & pose
)
{
//This is initialize for solo app (not from Enzdes)
	// Obtain all necessary user input
	if( motif_library_.empty() ) {
		MotifLibrary motifs( get_LigandMotifLibrary_user() );
		MotifCOPs motifcops = motifs.library();
		motif_library_ = motifcops;
	} // if it's not empty that means that the app must have filled the motif_library_
	core::conformation::ResidueOPs conformers( get_targetconformers_user() );
	target_conformers_map_ = setup_conformer_map( conformers );

	DnaDesignDefOPs build_position_defs;
	position_vector_setup( pose );
}
void
LigandMotifSearch::initialize(
	Pose const & pose,
	utility::vector1< Size > & input_BPs
)
{
//This is initialize for solo app (not from Enzdes)
	// Obtain all necessary user input
	for (  utility::vector1< Size  >::const_iterator position( input_BPs.begin() );
				 position != input_BPs.end(); ++position ) {
		ms_tr << "In init, Design position: " << *position << std::endl;
	}
	if( motif_library_.empty() ) {
		MotifLibrary motifs( get_LigandMotifLibrary_user() );
		MotifCOPs motifcops = motifs.library();
		motif_library_ = motifcops;
	} // if it's not empty that means that the app must have filled the motif_library_
	core::conformation::ResidueOPs conformers( get_targetconformers_user() );
	target_conformers_map_ = setup_conformer_map( conformers );

	DnaDesignDefOPs build_position_defs;
	position_vector_setup( pose );

	if ( ! input_BPs.empty() ) {
		for( Size i(1); i <= input_BPs.size(); ++i ) {
			BuildPosition_from_Size( pose, input_BPs[i] );
		}
	}	else {
		ms_tr << "Didn't give any BPs, so I'm gonna quit." << std::endl;
	}
}

void
LigandMotifSearch::incorporate_motifs(
	Pose const & pose
)
{
	core::pose::Pose posecopy( pose );
	core::pose::Pose posecopy2( pose );

	// Setup output files for motifs and rotamers
	utility::io::ozstream motif_output_file;
	utility::io::ozstream data_output_file;
	utility::io::ozstream qd_output_file;
	if( output_ ) {
		motif_output_file.open_append( output_filename_ );
	}
	if( data_ ) {
		data_output_file.open_append( data_filename_ );
	}
	if( quick_and_dirty_ ) {
		std::string motif_filename( pose.pdb_info()->name() );
		motif_filename.erase( motif_filename.end()-4, motif_filename.end() );
		std::string qd_output_filename = motif_filename + ".qd_motifs";
		qd_output_file.open_append( qd_output_filename );
	}

	//Here we're going to get triplets from the ligand in the protein which is being searched
	using namespace core;
	using namespace ObjexxFCL;
	using namespace pose;
	using namespace chemical;
	using namespace scoring;
	using namespace optimization;

	using utility::vector1;
	int ligand_resi_number = 0;
	utility::vector1< utility::vector1< Size  > > search_lig_triplets;
//	utility::vector1< utility::vector1< Size > > motif_indices_list;
	utility::vector1< utility::vector1< utility::vector1< Size > > > motif_indices_list; //Now, instead of having a vector of triplet vectors containing a size for the atom number, change size to vector of size to contain atom number and atomtype size
	//motif_indices_list[all triplets][atom i/j/k][1 is atom number, 2 is AtomType integer]
	utility::vector1< utility::vector1< utility::vector1< Size > > > all_motif_indices;

	int nres( pose.total_residue() );
	for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
		ResidueType const & lig_type( pose.residue_type( lig_pos ) );

		if(  !lig_type.is_ligand() ) continue;
		ms_tr << "in ligand splitter block, found my ligand, lig_pos is " << lig_pos << std::endl;
		ligand_resi_number = lig_pos;
// This is to make a ligres object once we find our ligand
		core::conformation::ResidueOP ligres( new core::conformation::Residue (pose.residue( lig_pos ) ) );

// This is to make an atomtypeset to get atomtype integers
					core::chemical::AtomTypeSetCOP atset = core::chemical::ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );

		for(core::Size atom_i = 1; atom_i <= ligres->natoms(); ++atom_i) {

			//std::cout << "in atom iterate block, atom num is" << atom_i <<   std::endl;
			std::string const atom_name = ligres->atom_name(atom_i);
			if( ligres->atom_is_hydrogen(atom_i) ) { continue; }
			//std::cout << "atom name is " << atom_name <<  std::endl;

			// This is a for loop to iterate over each atom's connected atoms:
			core::conformation::Residue::AtomIndices atom_i_connects(  ligres->bonded_neighbor( atom_i ) );
			for(core::Size atom_j = 1; atom_j <= atom_i_connects.size(); ++ atom_j ) {
				if( ligres->atom_is_hydrogen(atom_i_connects[atom_j]) ) { continue; }
				// std::cout << "ATOM j: " << atom_i_connects[atom_j] << " Name: " << ligres.atom_name(atom_i_connects[atom_j]) << std::endl;

				// This is the next for loop to find connects for the second atom, giving us the final atom number (atom k)
				core::conformation::Residue::AtomIndices atom_j_connects(  ligres->bonded_neighbor( atom_i_connects[atom_j] ) );
				for(core::Size atom_k = 1; atom_k <= atom_j_connects.size(); ++ atom_k ) {
					if( ligres->atom_is_hydrogen(atom_j_connects[atom_k]) ) { continue; }
					//std::cout << "ATOM k: " << atom_j_connects[atom_k] << " Name: " << ligres.atom_name(atom_j_connects[atom_k]) << std::endl;

					chemical::AtomType atom_i_type(ligres->atom_type(atom_i));
					std::string atom_i_name = atom_i_type.atom_type_name();
					//Size atom_i_int = atset->atom_type_index(atom_i_name);
				// std::cout << "ATOM j: " << atom_i << " Name: " << atom_i_name << " Int: " << atom_i_int << std::endl;


					//std::cout << "Connected triplet is: " << atom_i << ", type is " << atom_i_name  << ", ";
					//std::cout << atom_i_connects[atom_j] << ", type is " << ligres.atom_type(atom_i_connects[atom_j]).atom_type_name() << ", " ;
					//std::cout << atom_j_connects[atom_k] << ", type is " << ligres.atom_type(atom_j_connects[atom_k]).atom_type_name() << " " << std::endl;
					if ( atom_i != atom_j_connects[atom_k] ) {

						//make the 3 atom vector
						utility::vector1< utility::vector1< Size > > cur_motif_indices;

						utility::vector1< Size > atom_i_vector;
						atom_i_vector.push_back( atom_i );
						atom_i_vector.push_back( atset->atom_type_index( ligres->atom_type(atom_i).atom_type_name() ) );

						utility::vector1< Size > atom_j_vector;
						atom_j_vector.push_back( atom_i_connects[atom_j] );
						atom_j_vector.push_back( atset->atom_type_index( ligres->atom_type(atom_i_connects[atom_j]).atom_type_name() ) );

						utility::vector1< Size > atom_k_vector;
						atom_k_vector.push_back( atom_j_connects[atom_k] );
						atom_k_vector.push_back( atset->atom_type_index( ligres->atom_type(atom_j_connects[atom_k]).atom_type_name() ) );

						cur_motif_indices.push_back( atom_i_vector);
						cur_motif_indices.push_back( atom_j_vector);
						cur_motif_indices.push_back( atom_k_vector);

						//check if current index is a duplicate
						//don't need to do that here, but I forget what code I can delete.  I'm just going to change bool assignment.

						bool current_is_duplicate ( false );
						if ( !motif_indices_list.empty() ) {
							for(core::Size cur_motif_check = 1; cur_motif_check <= motif_indices_list.size(); ++ cur_motif_check) {
								// run a test to see if vectors are identical
								utility::vector1< utility::vector1< Size > > cur_mainlist_parent ( motif_indices_list[cur_motif_check] );
								utility::vector1< Size > cur_mainlist;
								utility::vector1< Size > cur_index;
								for(core::Size slice_parent = 1; slice_parent <= 3; ++ slice_parent) {
								cur_mainlist.push_back( cur_mainlist_parent[slice_parent][1] ); }
								for(core::Size slice_cur = 1; slice_cur <= 3; ++ slice_cur) {
								cur_index.push_back( cur_motif_indices[slice_cur][1] ); }
								utility::vector1< Size > sort_from_curindex (cur_index) ;
								utility::vector1< Size > sort_from_mainlist (cur_mainlist) ;
								std::sort(  sort_from_mainlist.begin(),  sort_from_mainlist.end() );
								std::sort( sort_from_curindex.begin(), sort_from_curindex.end() );
								if ( sort_from_mainlist[1] == sort_from_curindex[1] && sort_from_mainlist[2] == sort_from_curindex[2] && sort_from_mainlist[3] == sort_from_curindex[3] ) {
									current_is_duplicate = false;
								}
							}
						}
						//Making master list to see what's getting pruned with my check
						all_motif_indices.push_back( cur_motif_indices );
						//if current index isn't a duplicate, add it to the vector
						if ( !current_is_duplicate ) {
							motif_indices_list.push_back( cur_motif_indices );
						}
					}
				}
			}
		}
		std::cout << "Total 3 atoms in unpruned indices list is: " << all_motif_indices.size()  << std::endl;
		for( int prot_pos = 1 ; prot_pos <= nres ; ++prot_pos ) {
			ResidueType const & prot_type( pose.residue_type( prot_pos ) );
			if(  !prot_type.is_protein() ) continue;
			if( pose.residue( prot_pos ).name3() == "GLY" ) continue;

			// map will automatically sort the "contacts" with the lowest total_score at the front of map
			std::map< Real, Size > contacts;
			std::map< Real, Size  > distance_sorter;

			// Loop over positions, skipping non-ligand
			for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
				ResidueType const & lig_type( pose.residue_type( lig_pos ) );

				if(  !lig_type.is_ligand() ) continue;
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////////////
////////////// end of triplet search for current ligand //////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//Now we want to prune the motif library as well as the triplet list.
//Few reasons: first, test to see if we have motifs for the current ligand.  This is not a given.
//Also pruning will speed up considerably.  Don't look for motifs from orphan triplets, don't search through unrepresented motifs.

	MotifCOPs motif_library;
	// If the BP has motifs in it at this point then they came from an input file from the cmd line
	motif_library = motif_library_;

//Make new motif library
	MotifCOPs pruned_motif_library;

//Make new triplet list
// utility::vector1< utility::vector1< Size > > pruned_indices_list;
// std::map< Real, Size > contacts;

//For each motif in library

	for( protocols::motifs::MotifCOPs::iterator motifcop_itr = motif_library.begin(), end_itr = motif_library.end();
			 motifcop_itr != end_itr; ++motifcop_itr ) {

		protocols::motifs::MotifCOP motifcop( *motifcop_itr );
		int motif_atom1_int(motifcop->res2_atom1_int());
		int motif_atom2_int(motifcop->res2_atom2_int());
		int motif_atom3_int(motifcop->res2_atom3_int());
// ms_tr << "Motif atom 1 type: " << motif_atom1_name << ";  Motif atom 2 type: " << motif_atom2_name <<  ";  Motif atom 3 type: " << motif_atom3_name <<  std::endl;
		//For each triplet in ligand
		for(core::Size cur_lig_trip = 1; cur_lig_trip <= motif_indices_list.size(); ++ cur_lig_trip) {
			int ligand_atom1_int(motif_indices_list[cur_lig_trip][1][2]);
			int ligand_atom2_int(motif_indices_list[cur_lig_trip][2][2]);
			int ligand_atom3_int(motif_indices_list[cur_lig_trip][3][2]);
		//	 ms_tr << "Motif types: " <<  motif_atom1_int << " " <<  motif_atom2_int << " " <<  motif_atom3_int << ", Ligand types: " << ligand_atom1_int << " " << ligand_atom2_int << " " << ligand_atom3_int << " " << std::endl;
			//Check to see if match
//	 ms_tr << "Checking if match..."  <<  std::endl;
			if (
				motif_atom1_int	== ligand_atom1_int && motif_atom2_int == ligand_atom2_int  && motif_atom3_int == ligand_atom3_int ) {
//also OK if it's reversed
 // ms_tr << "It's a match!"  <<  std::endl;
				//also add motif to new motif library
				pruned_motif_library.push_back(motifcop);
				break;
			}

			//If match, add motif to new library, and also add triplet to new triplet list
		} // end ligand triplet iteration
	} // end library iteration
	ms_tr << "Unpruned motifs: "  << motif_library.size() << std::endl;
	Size motif_library_size = pruned_motif_library.size();
	Size motif_percent_chunk;
	double double_chunk = std::floor( (double) motif_library_size / 10 ) ;  // REQUIRED FOR WINDOWS
	motif_percent_chunk = Size ( double_chunk );
	//Size next_motif_percent = motif_percent_chunk;
	ms_tr << "Pruned motifs: "  << motif_library_size << std::endl;

///////////////////////////////////////////////////////////////////////////
/////////////////// done with pruning, ready to search! ///////////////////
///////////////////////////////////////////////////////////////////////////

	// for every protein backbone position (motif build position)
	for( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
			 ir != end_ir; ++ir ) {
		Size motif_counter=0;
		Size motif_percent=0;
		Size next_motif_percent = motif_percent_chunk;

		Size  trip_atom_1;
		Size  trip_atom_2;
		Size  trip_atom_3;
		// Map of all of the very best residues for each amino acid type, to make sure I don't add 2,000 Args and only 2 Tyrs
		std::map< std::string, std::map< Real, MotifHitOP > > best_mhits_all;

		// If we have rotamers coming in from files and they weren't cleared in initialization, then the search won't happen on this BuildPosition
		if ( ! ((*ir)->best_rotamers()).empty() ) continue;
		ms_tr << "WORKING ON PROTEIN POSITION " << (*ir)->seqpos() << std::endl;
		ms_tr << "NATIVE AA IS " <<  pose.residue( (*ir)->seqpos()  ).name3()  << std::endl;

//Setting up rmsd list
		std::map< Real, std::string > rmsd_list;


		MotifCOPs bp_best_motifs( (*ir)->best_motifs() );
		core::pack::rotamer_set::Rotamers bp_best_rotamers( (*ir)->best_rotamers() );
		// Need to clear the best_rotamers and best_motifs before I collect new ones
		(*ir)->clear_data();

		Size seqpos( (*ir)->seqpos() );
		std::stringstream firstline;
		firstline << "POSITION " << seqpos;
		if( output_ ) {
			motif_output_file << firstline.str() << "\n";
		}
		if( data_ ) {
			data_output_file << firstline.str() << "\n";
		}
		if( quick_and_dirty_ ) {
			qd_output_file << firstline.str() << "\n";
		}


		std::map< core::Size, core::pack::rotamer_set::RotamerSetOP > rotamer_sets;
		if( bp_best_rotamers.empty() ) {
			for( Size i(1); i <= core::chemical::num_canonical_aas; ++i ) {
				utility::vector1< bool > aa_info( core::chemical::num_canonical_aas, false );
				aa_info[i] = true;
				bool bump_yes( true );
				//need to make this an option
				core::pack::rotamer_set::RotamerSetOP rotset = build_rotamers_lite( posecopy, seqpos, aa_info, rot_level_, bump_yes );
				rotamer_sets[i] = rotset;
			}
		} else {
			Size bp_rots( bp_best_rotamers.size() );
			for( Size i(1); i <= core::chemical::num_canonical_aas; ++i ) {
				core::pack::rotamer_set::RotamerSetFactory rsf;
				core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( posecopy.residue((*ir)->seqpos()) );
				for( Size r(1); r <= bp_rots; ++r ) {
					if( bp_best_rotamers[r]->name3() == core::chemical::name_from_aa(core::chemical::AA(i)) ) {
						rotset->add_rotamer( *((bp_best_rotamers)[r]) );
					}
				}
				rotamer_sets[i] = rotset;
			}
		}

		// For every motif in the motif_library_
		for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = pruned_motif_library.begin(), end_itr = pruned_motif_library.end();
				 motifcop_itr != end_itr; ++motifcop_itr ) {
			++motif_counter;
			if ( motif_counter == next_motif_percent ) {
				next_motif_percent = motif_percent_chunk + next_motif_percent ;
				motif_percent = 10 + motif_percent;
				ms_tr << "On motif number: " << motif_counter << " and we are " << motif_percent << "% through library" << std::endl;


			}
			// WARNING: everything in this code assumes that residue 1 in the motif is the protein position and residue 2 is the dna
			bool passed_quick_and_dirty(false);
			protocols::motifs::MotifCOP motifcop( *motifcop_itr );
			//ms_tr << "WORKING ON MOTIF: " << motifcop->remark() << std::endl;

			// The BuildPosition may at some point have the ability to restrict ahead of time
			// currently I am not sure how to get restrictions from resfiles and so on . . .
			// it may not be an issue because there are ways other than resfiles, such as the defs for protein positions
			// and the main reason I would want this functionality is for homology models and theoretically the input pose
			// should already have all of the position types set as wild-type and be relaxed and minimized
			std::set< std::string > allowedtypes( (*ir)->allowed_types() );

			// Tighter cutoffs for residues with 3 and 4 chi angles
			Real rmsd_cutoff_1 = rmsd_cutoff_1_;
			Real rmsd_cutoff_2 = rmsd_cutoff_2_;
			Real dtest_cutoff = dtest_cutoff_;
			if( (motifcop->restype_name1() == "ARG") && ( ! quick_and_dirty_ ) ) {
				rmsd_cutoff_1 = (rmsd_cutoff_1 / 2.0 );
				rmsd_cutoff_2 = (rmsd_cutoff_2 / 2.0 );
				dtest_cutoff = (dtest_cutoff / 2.0 );
			}
			if( ( (motifcop->restype_name1() == "MET") || (motifcop->restype_name1() == "LYS") || (motifcop->restype_name1() == "GLU") ||  (motifcop->restype_name1() == "GLN") ) && ( ! quick_and_dirty_ ) ) {
				rmsd_cutoff_1 = (rmsd_cutoff_1 / 1.33 );
				rmsd_cutoff_2 = (rmsd_cutoff_2 / 1.33 );
				dtest_cutoff = (dtest_cutoff / 1.33 );
			}
			//utility::vector1< std::string > atoms;
			/*
						// Atoms required for the parallel base test
						// maybe these vectors should be stored as a part of this class?
						// for non-DNA motif searches I won't be using the z-test at all
						if( protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "GUA" || protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "ADE" ) {
							atoms.push_back("C5");
							atoms.push_back("C6");
							atoms.push_back("N3");
							atoms.push_back("C2");
							atoms.push_back("N1");
							atoms.push_back("C4");
						} else if( protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "CYT" || protocols::dna::dna_full_name3( motifcop->restype_name2() ) == "THY" ) {
							atoms.push_back("C5");
							atoms.push_back("C4");
							atoms.push_back("N1");
							atoms.push_back("C2");
							atoms.push_back("N3");
							atoms.push_back("C6");
						} else {
							ms_tr << "Residue you are planning to do parallel base test with is not a DNA base!" << std::endl;
						}
			*/ //commented out parallel base test--doesn't work for ligands
			bool automorphism(false);
			bool passed_automorphism(false);
			if(
				motifcop->restype_name1() == "ASP" ||
				motifcop->restype_name1() == "GLU" ||
				motifcop->restype_name1() == "PHE" ||
				motifcop->restype_name1() == "LEU" ||
				motifcop->restype_name1() == "ARG" ||
				motifcop->restype_name1() == "TYR" ||
				motifcop->restype_name1() == "VAL"
			) {
				automorphism = true;
			}

			// Related to only using motifs that include an allowed amino acid at the BuildPosition
			// At this point the allowedtypes is always going to be empty
			bool allowed(false);
			if( allowedtypes.empty() ) {
				allowed = true;
			}
			for( std::set< std::string >::const_iterator ir2(allowedtypes.begin() ), end_ir2( allowedtypes.end() );
					 ir2 != end_ir2; ++ir2 ) {
				if( (*ir2) == motifcop->restype_name1() ) {
					allowed = true;
				}
			}
			if( ! allowed ) continue;

			// Create a standard dna base of type in motif
			// For ligand motifs we will need to create a residue of the type we are trying to target
			// Since there won't be a second residue in the motif
			// restype_name2() will need to be set to something in the motif?
			std::string basetype( motifcop->restype_name2() );
			// If there are conformers
			core::conformation::ResidueOPs DNAResidueOPs( target_conformers_map_[basetype] );
			//core::conformation::ResidueOPs LigandResidueOPs( target_conformers_map_[basetype]);
			bool noconformers( false );
			if( DNAResidueOPs.empty() ) {
				noconformers = true;
			}
			core::conformation::ResidueOP ligres = core::conformation::ResidueFactory::create_residue( posecopy.residue( ligand_resi_number ).type() );
			// Get the rotamer set with type matching motif at the motif build position
			core::pack::rotamer_set::RotamerSetOP rotset = rotamer_sets[ core::chemical::aa_from_name(motifcop->restype_name1()) ];

			Real final(100);
			std::pair< core::conformation::ResidueOP, core::conformation::ResidueOP > bestpair;
			bool b_bestpair( false );
			Size rs1( rotset->num_rotamers() );

			Real rmsdtest_ir2(100.0);

			// This atom type is only C1' because that atom is common
			// to all DNA bases and it is basically in the plane of the base
			// Right now we're looking at a motif -- we don't have triplet yet.  Don't place because maybe nothing in common (need to do test once we're in triplet loop).

			// Important in case of the very strange situation where more than one
			// target_position is close enough to pass the tests
			bool tftest(false);
			utility::vector1< utility::vector1< Size > > bestpos;
			Real test(1000);
// Used to iterate over target_positions, now we will iterate over triplets in current ligand
			Sizes target_positions( (*ir)->target_positions() );


	//utility::vector1< utility::vector1< utility::vector1< Size > > > motif_indices_list; //Now, instead of having a vector of triplet vectors containing a size for the atom number, change size to vector of size to contain atom number and atomtype size
	//motif_indices_list[all triplets][atom i/j/k][1 is atom number, 2 is AtomType integer]
			for( utility::vector1< utility::vector1< utility::vector1< Size > > >::const_iterator curtrip( motif_indices_list.begin() ), end( motif_indices_list.end() ); //prof II
					 curtrip != end; ++curtrip ) {

				//before: bpos now: curtrip
				//bpos was DNA build position
				//		std::set< std::string > allowed_types( target_positions_[*bpos] ); nonsensical now
				// Should have a test to ensure that it has the atom types I'll be using?
				// At this point, I know what my triplet is, and I know what my motif is--let's se if they match.  Otherwise, we'll continue to the next triplet.
				utility::vector1<  utility::vector1< Size > > deref_trip( *curtrip );
				core::conformation::ResidueOP check_ligand = ligres;
				int motif_atom1_int(motifcop->res2_atom1_int());
				int motif_atom2_int(motifcop->res2_atom2_int());
				int motif_atom3_int(motifcop->res2_atom3_int());
				//For each triplet in ligand
				int ligand_atom1_int( deref_trip[1][2]);
				int ligand_atom2_int( deref_trip[2][2]);
				int ligand_atom3_int( deref_trip[3][2]);
				//Check to see if match
				if (
					motif_atom1_int  == ligand_atom1_int && motif_atom2_int == ligand_atom2_int  && motif_atom3_int == ligand_atom3_int )
					//not OK if it's reversed--we'll find it because we don't prune ligand triplets
				{
					//std::cout << "Motif atom 1 type: " << motif_atom1_name << ";  Motif atom 2 type: " << motif_atom2_name <<  ";  Motif atom 3 type: " << motif_atom3_name <<  std::endl;
					//std::cout << "Lig atom 1 type: " << ligand_atom1_name << ";  Lig atom 2 type: " << ligand_atom2_name <<  ";  Lig atom 3 type: " << ligand_atom3_name <<  std::endl;
					// ms_tr << "It's a match! (INSIDE LOOPS NOW)"  <<  std::endl;
				} else {
					continue;
				}


				//Instead of parallel base test, let's check the similarity of the motifs.  A1-A2-A3, B1-B2-B3.  Compare D(A1,A2)~=D(B1,B2); D(A2,A3)~=D(B2,B3); Angle(A1-A2-A3)~=Angle(B1-B2-B3)
				//These tests will uniquely specify the three atoms on either side, and so they will show us if the motif and triplet are a match.
				//Want 0.3 A agreement between distances
				//Want 30 degree agreement between angles
				//Can tighten later if needed
				//Turns out this is a real bitch to code, not going to do it until later -Matt
				/*     numeric::xyzVector motif_vector();
						numeric::xyzVector ligand_vector();
						//For each atom in motif triplet
				//	Vector const res2_C  = res2.xyz(  "C" );
				  //   Vector const motif_atom1_xyz( motifcop->forward_jump() );

					//For each atom in ligand triplet
				    // Vector const ligand_atom1_xyz( check_ligand.atom_type(deref_trip[1]).atom_type_name() );

							ligand_vector += check_ligand.xyz( deref_trip[1]) ;
							ligand_vector += check_ligand.xyz( deref_trip[2]) ;
							ligand_vector += check_ligand.xyz( deref_trip[3]) ;

						  Real ligand_distance_AB( check_ligand.xyz( deref_trip[1]  ).distance(  check_ligand.xyz( deref_trip[2]  )  ) );
						  Real ligand_distance_BC( check_ligand.xyz( deref_trip[2]  ).distance(  check_ligand.xyz( deref_trip[3]  )  ) );
				//						src/numeric/xyzTriple.hh:       angle_of( xyzTriple const & a, xyzTriple const & b )

				*/

				//We are looping over the rotamer set here.  Matt put the rotamer loop here so that we wouldn't be checking motif/triplet matching for every rotamer--that would take forever!
				for( Size ir2(1); ir2 <= rs1; ++ir2 ) { //rotamer loop

					//Used to be before this loop, now put it here because we didn't know which atom to place yet
					core::conformation::Atom atm( check_ligand->atom( deref_trip[2][1] ) );
					core::conformation::Atom auto_atm( check_ligand->atom( deref_trip[2][1] ) );

					Size  trip_atom_1(deref_trip[1][1]);
					Size  trip_atom_2(deref_trip[2][1]);
					Size  trip_atom_3(deref_trip[3][1]);
					//utility::vector1< Size > atoms = new utility::vector1< Size >( {trip_atom_1, trip_atom_2, trip_atom_3} );
					//Size myatoms[] = {trip_atom_1, trip_atom_2, trip_atom_3};
					utility::vector1< Size > atoms;
					atoms.push_back(trip_atom_1);
					atoms.push_back(trip_atom_2);
					atoms.push_back(trip_atom_3);

//					motifcop->place_atom( *(rotset->nonconst_rotamer(ir2)), *check_ligand, atm, trip_atom_1, trip_atom_2, trip_atom_3, atm_str ); //prof II
					motifcop->place_atom( *(rotset->nonconst_rotamer(ir2)), *check_ligand, atm, trip_atom_1, trip_atom_2, trip_atom_3, trip_atom_2 );
					if( automorphism ) {
						motifcop->place_atom( *(rotset->nonconst_rotamer(ir2)), *check_ligand, auto_atm, trip_atom_1, trip_atom_2, trip_atom_3, trip_atom_2, false ); //the false just flips the placement
					}


					Real dtest1( atm.xyz().distance( posecopy.residue( ligand_resi_number  ).xyz( deref_trip[2][1] ) ) );
// For dtest1, we are getting the distance from the motif base which has already been placed to the DNA base in the build.  This is very similar to what we're going
					Real dtest1_auto(100);
					if( automorphism ) {
						dtest1_auto = ( auto_atm.xyz().distance( posecopy.residue( ligand_resi_number ).xyz( deref_trip[2][1] ) ) );


						//ms_tr << "RMSD between ligand resi (rosetta #) " <<  ligand_resi_number << " and motif ligand = " << dtest1 << " dtestauto is " << dtest1_auto << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
					}

					if( ! automorphism ) {
						if( dtest1 > dtest_cutoff ) continue;
					} else {
						if( dtest1 > dtest_cutoff && dtest1_auto > dtest_cutoff ) continue;
						if( dtest1 < dtest_cutoff && dtest1_auto < dtest_cutoff ) {
							if( dtest1_auto < dtest1 ) {
								passed_automorphism = true;
							} else {
								passed_automorphism = false;
							}
						} else if( dtest1_auto < dtest_cutoff ) {
							passed_automorphism = true;
						} else {
							passed_automorphism = false;
						}
					} //All of these tests are good, we'll keep for ligand motifs

					//					ms_tr << "Passed first round tests on 712: RMSD between ligand resi (rosetta #) " <<  ligand_resi_number << " and motif ligand = " << dtest1 << " dtestauto is " << dtest1_auto << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;

					core::conformation::ResidueOP posebase( new core::conformation::Residue( posecopy.residue( ligand_resi_number ) ) );
					if( passed_automorphism ) {
//  rmsd_list[dtest1_auto] = motifcop->restype_name1() ; //(ADD TO THE MAP)
						motifcop->place_atoms( *(rotset->nonconst_rotamer(ir2)), *posebase, atoms, trip_atom_1, trip_atom_2, trip_atom_3, false );
					} else {
//  rmsd_list[dtest1] = motifcop->restype_name1() ; //(ADD TO THE MAP)
						motifcop->place_atoms( *(rotset->nonconst_rotamer(ir2)), *posebase, atoms, trip_atom_1, trip_atom_2, trip_atom_3 );
					}

					if( quick_and_dirty_ ) {
						qd_output_file << *motifcop;
						passed_quick_and_dirty = true;
						break;
					}
					Real rmsdtest = atom_specific_rms( *posebase, posecopy.residue(ligand_resi_number), atoms );
					if( rmsdtest > rmsd_cutoff_1 ) continue;
					//			ms_tr << "at 779, RMSD cutoff is " << rmsd_cutoff_1  << " and rmsd is " <<  rmsdtest << std::endl;
					//ms_tr << "Passed 728! " << std::endl;

					rmsdtest_ir2 = rmsdtest;
					if( rmsdtest < test ) {
						test = rmsdtest;
						bestpos = deref_trip;
						posecopy2 = posecopy;
						tftest = true;
					} // if( rmsdtest < test )
					if( passed_quick_and_dirty ) break;

					if( passed_quick_and_dirty ) break;
					if( ! tftest ) continue;
					//ms_tr << "Passed 741! " << std::endl;

					MotifHitOP motifhit( new MotifHit( *motifcop, ligand_resi_number, passed_automorphism ) );

					// If there are no conformers will use the base from the pose, so rmsd can be 0 theoretically
					if( noconformers ) {
						core::conformation::ResidueOP posebase2( new core::conformation::Residue( posecopy2.residue( ligand_resi_number ) ) );
						if( passed_automorphism ) {
							motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), *posebase2,trip_atom_1, trip_atom_2, trip_atom_3 , false );
						} else {
							motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), *posebase2,trip_atom_1, trip_atom_2, trip_atom_3  );
						}
						Real rmsdtest2 = rmsdtest;
						//Real rmsdtest2 = core::scoring::automorphic_rmsd( *posebase2, posecopy2.residue(ligand_resi_number), false );
						Real finaltest = ( ( rmsdtest2)  );
						//std::cout << "FINAL: " << finaltest << std::endl;
						//std::cout << "cutoff: " << rmsd_cutoff_2 << " ir2 " << rmsdtest_ir2 << std::endl;
						//NOTE: Do I want to keep any statistics about percentages of passing certain cutoffs??
						//ms_tr << "Passed 755! " << std::endl;
						if( (rmsdtest_ir2 < rmsd_cutoff_2)  ) {
			//ms_tr << "adding to map, RMSD cutoff is " << rmsd_cutoff_2 << " and rmsdtest is " << rmsdtest_ir2 << std::endl;
							rmsd_list[rmsdtest_ir2] = motifcop->restype_name1() ; //(ADD TO THE MAP)
							//ms_tr << "Passed 757! RMSD between DNA resi (rosetta #), no conformers " <<  ligand_resi_number << " and motif DNA = " << rmsdtest_ir2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
							if( data_ ) {
								data_output_file << "Passed 759! RMSD between DNA resi (rosetta #), no conformers " << ligand_resi_number << " and motif DNA = " << rmsdtest_ir2  << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
							}
							motifhit->final_test( finaltest );
							motifhit->build_rotamer( *(rotset->nonconst_rotamer(ir2) ) );
							motifhit->target_conformer( *posebase2 );
							best_mhits_all[motifcop->restype_name1()][finaltest] = motifhit->clone();
							if ( finaltest < final ) {
								b_bestpair = true;
								bestpair = std::make_pair( (rotset->nonconst_rotamer(ir2))->clone(), posebase2->clone() );
								final = finaltest;
							}
						} // if passed the second round of tests
					} else {
						for( core::conformation::ResidueOPs::const_iterator resop( DNAResidueOPs.begin() ), end( DNAResidueOPs.end() );
								 resop != end; ++resop ) {
							if( passed_automorphism ) {
								motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), **resop, trip_atom_1, trip_atom_2, trip_atom_3 , false );
							} else {
								motifcop->place_residue(*(rotset->nonconst_rotamer(ir2)), **resop,trip_atom_1, trip_atom_2, trip_atom_3  );
							}
							Real rmsdtest2 = rmsdtest;
							//Real rmsdtest2 = core::scoring::automorphic_rmsd( **resop, posecopy2.residue(ligand_resi_number), false );
							if( rmsdtest2 > rmsd_cutoff_2 ) continue;
							Real finaltest = ( ( rmsdtest2 * 100 )  );
							Real finaltestc = ( ( rmsdtest_ir2 * 100 )  );
							//NOTE: Do I want to keep any statistics about percentages of passing certain cutoffs??
							//ms_tr << "Passed 784! " << std::endl;
							if( (rmsdtest2 < rmsd_cutoff_2)  ) {
								rmsd_list[rmsdtest2] = motifcop->restype_name1() ; //(ADD TO THE MAP)
								//		ms_tr << "at 843, adding to map, RMSD cutoff is " << rmsd_cutoff_2 << " and rmsdtest is " << rmsdtest_ir2 << std::endl;
								//ms_tr << "Passed 781! RMSD between DNA resi (rosetta #) " << ligand_resi_number << " and motif DNA = " << rmsdtest2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
								if( data_ ) {
									data_output_file << "Passed 783! RMSD between DNA resi (rosetta #) " << ligand_resi_number << " and motif DNA = " << rmsdtest2 << " and combined score = " << finaltest << " for residue type " << motifcop->restype_name1() << ", rotamer # " << ir2 << ", motif named " << motifcop->remark() << std::endl;
								}
								motifhit->final_test( finaltest );
								motifhit->build_rotamer( *(rotset->nonconst_rotamer(ir2) ) );
								motifhit->target_conformer( **resop );
								best_mhits_all[motifcop->restype_name1()][finaltestc] = motifhit->clone();
								if ( finaltest < final ) {
									b_bestpair = true;
									bestpair = std::make_pair( (rotset->nonconst_rotamer(ir2))->clone(), (*resop)->clone() );
									final = finaltest;
								}
							} // if conformer passed the second round of tests
						} // loop over conformers
					} // if not noconfomers
				}  // loop over first rotamer set
			} // loop over triplets in current ligand

			// Dump the best rotamer and conformer pair for each motif, mainly for debugging purposes
			if( dump_motifs_ ) {
				if( b_bestpair ) {
					core::pose::Pose pose_dump2;
					pose_dump2.append_residue_by_jump( *(bestpair.second), 1);
					pose_dump2.append_residue_by_jump( *(bestpair.first), 1);
					std::stringstream pose2_name_full;
					if( passed_automorphism ) {
						pose2_name_full << "Test_auto_" << motifcop->restype_name2()[0] << "_" << (*ir)->seqpos() << motifcop->restype_name1() << "_" << motifcop->remark() << ".pdb";
					} else {
						pose2_name_full << "Test_" << motifcop->restype_name2()[0] << "_" << (*ir)->seqpos() << motifcop->restype_name1()  << "_" << motifcop->remark() << ".pdb";
					}
					core::io::pdb::dump_pdb( pose_dump2, pose2_name_full.str() );
				}
			}
		}
		if( ! best_mhits_all.empty() ) {
			core::pose::Pose pose_dump( pose );
			for( std::map< std::string, std::map< Real, MotifHitOP > >::const_iterator bh( best_mhits_all.begin() ),
					 end( best_mhits_all.end() ); bh != end; ++bh ) {
				Size hits = 0;
				for( std::map< Real, MotifHitOP >::const_iterator bh2( (bh->second).begin() ),
						 end2( (bh->second).end() ); bh2 != end2; ++bh2 ) {
					MotifHitOP motifhitop( bh2->second );
					if( ! minimize_ ) {
						(*ir)->keep_rotamer( *(motifhitop->build_rotamer()) );
						++hits;
						if( hits > rots2add_ ) break;
					} else {
						using namespace core::scoring;
						// Need a copy of the pose to generate the constraints
						// Maybe there's a way to rewrite the constraint making code to avoid this?
						// No, because the minimize itself needs the residues to be a part of the pose?
						// Or maybe instead of copying the pose you could save the residue you are replacing and replace it back?
						//core::pose::Pose pose_dump( pose );
						core::conformation::ResidueOP build_rotamer( new core::conformation::Residue( *(motifhitop->build_rotamer()) ) );
						pose_dump.replace_residue( (*ir)->seqpos(), *build_rotamer, false );
						// The residue is probably already placed . . .
						//	if( passed_automorphism ) {
						//		motifcop->place_residue( motifhitop->build_rotamer(), motifhitop->target_conformer(), false );
						//	} else {
						//		motifcop->place_residue( motifhitop->build_rotamer(), motifhitop->target_conformer() );
						//		}
						if ( protocols::dna::dna_full_name3( pose_dump.residue(motifhitop->vbpos()).name3() ) != protocols::dna::dna_full_name3( (motifhitop->target_conformer())->name3() ) ) {
							make_base_pair_mutation( pose_dump, motifhitop->vbpos(), core::chemical::aa_from_name( protocols::dna::dna_full_name3( motifhitop->target_conformer()->name3() ) ) );
						}
						if( motifhitop->passed_automorphism() ) {
							(motifhitop->motifcop())->place_residue( pose_dump.residue( motifhitop->vbpos() ), *build_rotamer, trip_atom_1, trip_atom_2, trip_atom_3 , false );
						} else {
							(motifhitop->motifcop())->place_residue( pose_dump.residue( motifhitop->vbpos() ), *build_rotamer, trip_atom_1, trip_atom_2, trip_atom_3  );
							//Place_residue wants (fixed, mobile, ...)
						}
						/*core::pose::Pose pose_dump2( pose );
						pose_dump2.replace_residue( (*ir)->seqpos(), *build_rotamer, false );
						std::stringstream pose2_name_full;
						std::stringstream pose_name_full;
						std::stringstream pose3_name_full;
						pose2_name_full << "AfterReplace_" << bh2->first << ".pdb";
						pose_name_full << "BeforeReplace_" << bh2->first << ".pdb";
						pose3_name_full << "AfterMinBeforeReplace_" << bh2->first << ".pdb";
						core::io::pdb::dump_pdb( pose_dump2, pose2_name_full.str() );
						core::io::pdb::dump_pdb( pose_dump, pose_name_full.str() );*/

						constraints::ConstraintSetOP sc_cst_set( new constraints::ConstraintSet() );
						add_motif_sc_constraints( sc_cst_set, pose_dump, (*ir)->seqpos(), *build_rotamer, motifhitop->motifcop(), false );
						//add_motif_sc_constraints( sc_cst_set, pose_dump2, (*ir)->seqpos(), *build_rotamer, motifhitop->motifcop(), false );
						ScoreFunctionOP score_fxn( get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
						methods::EnergyMethodOptions options( score_fxn->energy_method_options() );
						score_fxn->set_energy_method_options( options );
						score_fxn->set_weight( coordinate_constraint, 10.0 );
						pose_dump.constraint_set( sc_cst_set );
						//pose_dump2.constraint_set( sc_cst_set );
						//core::Real pre_sc_constraint_check( pose_dump.energies().total_energies()[ coordinate_constraint ] );
						//core::Real pre_sc_constraint_check2( pose_dump2.energies().total_energies()[ coordinate_constraint ] );
						//ms_tr << "Before sidechain refinement constraints score is " << pre_sc_constraint_check << std::endl;
						//ms_tr << "2Before sidechain refinement constraints score is " << pre_sc_constraint_check2 << std::endl;
						/*if( data_ ) {
							data_output_file << "Before sidechain refinement constraints score is " << pre_sc_constraint_check << std::endl;
						}*/
						core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
						movemap->set_chi( (*ir)->seqpos(), true );
						protocols::simple_moves::MinMoverOP minmover( new protocols::simple_moves::MinMover( movemap, score_fxn, "dfpmin_armijo_nonmonotone_atol", 0.000001, true ) );
						minmover->apply( pose_dump );
						//core::io::pdb::dump_pdb( pose_dump, pose3_name_full.str() );
						core::Real sc_constraint_check( pose_dump.energies().total_energies()[ coordinate_constraint ] );
						ms_tr << "After sidechain refinement constraints score is " << sc_constraint_check << std::endl;
						/*if( data_ ) {
							data_output_file << "After sidechain refinement constraints score is " << sc_constraint_check << std::endl;
						}*/
						/*if( sc_constraint_check < 10.0 && sc_constraint_check < pre_sc_constraint_check ) {
							(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
							++hits;
							if( output_ ) {
								motif_output_file << *(motifhitop->motifcop());
								motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
							}
						} else if( pre_sc_constraint_check < 10.0 && sc_constraint_check > pre_sc_constraint_check ) {
							pose_dump.replace_residue( (*ir)->seqpos(), *(motifhitop->build_rotamer()), true );
							(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
							++hits;
							if( output_ ) {
								motif_output_file << *(motifhitop->motifcop());
								motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
							}
						}*/

					}		//Here is end of minimize code block, get rid of it eventually

					pose_dump.replace_residue( (*ir)->seqpos(), *(motifhitop->build_rotamer()), true );
					(*ir)->keep_rotamer( (pose_dump.residue((*ir)->seqpos())) );
					if( output_ ) {

						motif_output_file << *(motifhitop->motifcop());
						motif_output_file << "RESIDUE " << (pose_dump.residue((*ir)->seqpos()));
					}
					if( hits > rots2add_ ) break;
				}
			}
		}


//print out placed motifs by RMSD

		ms_tr << "RMSD, amino acid"  << std::endl;
		for( std::map< Real,  std::string >::const_iterator rmsd_it( rmsd_list.begin() ),  end( rmsd_list.end() ); rmsd_it != end; ++rmsd_it ) {
			ms_tr << rmsd_it->first << ", " << rmsd_it->second << std::endl;

		}
	}
	if( output_ ) {
		motif_output_file.close();
	}
	if( data_ ) {
		data_output_file.close();
	}
	if( quick_and_dirty_ ) {
		qd_output_file.close();
	}
}

core::pack::rotamer_set::Rotamers
LigandMotifSearch::get_rotamers()
{
	core::pack::rotamer_set::Rotamers best_rotamers;
	for ( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
				ir != end_ir; ++ir ) {
		if ( ! ((*ir)->best_rotamers()).empty() ) {
			best_rotamers = (*ir)->best_rotamers();
		}
	}
	return best_rotamers;
}
core::pack::rotamer_set::Rotamers
LigandMotifSearch::bp_rotamers(
	Size const seqpos
)
{
	core::pack::rotamer_set::Rotamers best_rotamers;
	for ( BuildPositionOPs::const_iterator ir( build_positionOPs_.begin() ), end_ir( build_positionOPs_.end() );
				ir != end_ir; ++ir ) {
		if( (*ir)->seqpos() != seqpos ) continue;
		if ( ! ((*ir)->best_rotamers()).empty() ) {
			best_rotamers = (*ir)->best_rotamers();
		} else {
			ms_tr << "There were no rotamers to be included for position " << seqpos << std::endl;
		}
	}
	return best_rotamers;
}

// Maybe this belongs with Motif.cc or MotifLibrary.cc, or both, but not here
bool
LigandMotifSearch::protein_dna_motif()
{
	using namespace core::chemical;
	// This function only works for a motif that is made up of only two residues
	bool protein_dna( false );
	// motif_library_ has to be filled with the input MotifLibrary before you can call this fxn
	if ( !	motif_library_.empty() ) {
		// Check to see if motif has a protein component and a dna component
		//THIS IS RIDICULOUS, DON'T MAKE RESIDUES, MAKE RESIDUE TYPE FROM NAME3
		core::conformation::ResidueOP res1 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3( motif_library_[1]->restype_name1() ) ) );
		core::conformation::ResidueOP res2 = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map( protocols::dna::dna_full_name3( motif_library_[1]->restype_name2() ) ) );
		if	( ( res1->is_protein() && res2->is_DNA() ) || ( res2->is_protein() && res1->is_DNA() ) ) {
			protein_dna = true;
		}
	} else {
		ms_tr << "MotifLibrary has not been initialized yet, cannot yet identify the type of motifs being used, assuming protein-DNA with possible disastrous consequences." << std::endl;
		protein_dna = true;
	}
	return protein_dna;
}

// Will make this more general!  Should keep the chains separate
// Maybe make a map of vectors, keeping track of the different chains
void
LigandMotifSearch::position_vector_setup(
	Pose const & pose
)
{
	for ( Size i(1), end( pose.total_residue() ); i <= end; ++i ) {
		if ( pose.residue_type(i).is_protein() ) {
			protein_positions_.push_back(i);
		}
		if ( pose.residue_type(i).is_DNA() ) {
			dna_positions_.push_back(i);
		}
	}
}

void
LigandMotifSearch::identify_motif_build_positions(
	Pose const & pose,
	utility::vector1< Size > & build_positions
)
{
	if ( protein_dna_motif() ) {
		protein_DNA_motif_build_positions_JA( pose, build_positions, dna_positions_ ); //could easily change this line to take one of Phil's DNA interface fxns that also fills a vector
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}


utility::vector1< core::Size >
LigandMotifSearch::get_sphere_aa(
	Pose const & pose,
	core::Real cut1
)
{
//This is where Matt will put enzdes sphere finding function to make build_position list
	using namespace core;
	using namespace ObjexxFCL;
	using namespace pose;
	using namespace chemical;
	using namespace scoring;
	using namespace optimization;
	int nres( pose.total_residue() );
	// std::cout << "In get_sphere_aa, about to find ligand " << std::endl;
	std::set< core::Size > interface_target_res;
	for( int lig_pos = 1 ; lig_pos <= nres ; ++lig_pos ) {
		ResidueType const & lig_type( pose.residue_type( lig_pos ) );
		// std::cout << "In get_sphere_aa, made res type, aa number is " << lig_pos << std::endl;

		if(  lig_type.is_ligand() ) {
			ms_tr << "in get_sphere_aa, found my ligand, lig_pos is " << lig_pos << std::endl;
			interface_target_res.insert( lig_pos) ;
		}
	}
	utility::vector1< core::Size > sphere_resi;

	{

		core::Real cut2 = cut1 + 1;
		core::Real cut3 = cut2 + 1;
		core::Real cut4 = cut3 + 1;
		core::Real cut1_sq = cut1 * cut1;
		core::Real cut2_sq = cut2 * cut2;
		core::Real cut3_sq = cut3 * cut3;
		core::Real cut4_sq = cut4 * cut4;

		for( std::set< core::Size >::const_iterator targ_it( interface_target_res.begin()),targ_end(interface_target_res.end());
				 targ_it != targ_end; ++targ_it ) {

			// on protein side, have to do distance check
			core::conformation::Residue const & targ_rsd = pose.residue( *targ_it );
			core::Size targ_res_atom_start = 1;
			if( targ_rsd.is_protein() ) {
				targ_res_atom_start = targ_rsd.first_sidechain_atom();
			}

			for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {
				core::conformation::Residue const & prot_rsd = pose.residue(i);
				for(core::Size k = targ_res_atom_start, k_end = targ_rsd.nheavyatoms(); k <= k_end; ++k) {
					core::Vector prot_cb, prot_ca;
					if( prot_rsd.has("CB") ) prot_cb = prot_rsd.xyz("CB");
					if( prot_rsd.has("CA") ) prot_ca = prot_rsd.xyz("CA"); // GLY
					core::Real ca_dist2 = targ_rsd.xyz(k).distance_squared( prot_ca );
					if( ca_dist2 <= cut4_sq ) {
						if( ca_dist2 <= cut3_sq ) {
							if( ca_dist2 <= cut2_sq ) {
								if( ca_dist2 <= cut1_sq) {
									sphere_resi.push_back(i);
									break;
								} // cut1
								else if( prot_rsd.has("CB") ) {
									core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
									//                tr.Info << "cb_dist2 is " << cb_dist2 << "; ";
									if( cb_dist2 < ca_dist2 ) {
										sphere_resi.push_back(i);
										break;
									}
								}  // end of non-gly residues
								else if ( prot_rsd.has("2HA") ) {   //glycine doesn't have a CB, so use 2HA to get position where CB would be
									// use the name "cb" to describe the 2HA atom; design if 2HA < CA
									prot_cb = prot_rsd.xyz("2HA");
									core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
									if( cb_dist2 < ca_dist2 ) {   // 2HA is closer than CA
										sphere_resi.push_back(i);
										break;
									}
								}  // end of gly residues
								else {  // Exception handling case for residue without CB or 2HA
									ms_tr << "Weird residue without CB or 2HA. Watch out! Residue: " << i << std::endl;
									break;
								} // end of exception catching for neither CB nor 2HA
							} //cut2
						} //cut3
					} //cut4
				} //loop over target res atoms
			} //loop over protein residues
		} //loop over target residues

		ms_tr << "Finished sphere check ok" << std::endl;

	} //find_design_interface

	return sphere_resi;
} // get_sphere_aa function

/*
void
LigandMotifSearch::fill_bp_allowed_types(
	Pose const & pose,
	Size const seqpos,
	std::set< std::string > & allowed_types
)
{
	// This doesn't work for some reason . . .
	// Not quite sure how to do it, so that it knows about resfile restrictions
	core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory;
	tf->push_back( new core::pack::task::ReadResfileOperation );//new protocols::dna::RestrictDesignToProteinDNAInterface );
	core::pack::task::PackerTaskOP task( tf->create_task_and_apply_taskoperations( pose ) );
	for( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter
			allowed_iter = task->residue_task( seqpos ).allowed_residue_types_begin(),
			allowed_end = task->residue_task( seqpos ).allowed_residue_types_end();
			allowed_iter != allowed_end; ++allowed_iter ) {
		std::cout << ' ' << (*allowed_iter)->name3();
	}
	std::cout << std::endl;
}
*/

void
LigandMotifSearch::identify_motif_BuildPositions(
	Pose const & pose
)
{
	Sizes positions(0);
	if ( protein_dna_motif() ) {
		Sizes target_positions( map2keyvector( target_positions_ ) );
		protein_DNA_motif_build_positions_JA( pose, positions, target_positions );
		for ( Sizes::const_iterator pos( positions.begin() ), end( positions.end() );
					pos != end; ++pos ) {
			Size seqpos(*pos);
			Sizes short_target_positions( shorten_target_list( pose, seqpos, target_positions ) );
			std::set< std::string > allowed_types; // this vector will remain empty in this sitatuation since there is no input Def to limit the types of amino acids allowed
			//fill_bp_allowed_types( pose, seqpos, allowed_types );
			BuildPositionOP build_position( new BuildPosition( seqpos, short_target_positions, allowed_types ) );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

void
LigandMotifSearch::BuildPosition_from_Size(
	Pose const & pose,
	Size const input_BP
)
{
	Sizes target_positions( map2keyvector( target_positions_ ) );
	Sizes short_target_positions( shorten_target_list( pose, input_BP, target_positions ) );
	std::set< std::string > allowed_types; // this set will remain empty in this sitatuation since there is no input Def to limit the types of amino acids allowed
	// MAKE ANOTHER FUNCTION THAT FILLS ALLOWED_TYPES VIA CHECKING WHAT IS ALLOWED
	//fill_bp_allowed_types( pose, input_BP, allowed_types );
	BuildPositionOP build_position( new BuildPosition( input_BP, short_target_positions, allowed_types ) );
	build_positionOPs_.push_back( build_position );
}

void
LigandMotifSearch::defs2BuildPositions(
	Pose const & pose,
	DnaDesignDefOPs const & defs
)
{
	if ( protein_dna_motif() ) {
		Sizes full_tl( map2keyvector( target_positions_ ) );
		std::map< Size, std::set< std::string > > mappositions( defs2map( pose, defs ) );
		for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
					end( mappositions.end() ); it != end; ++it ) {
			BuildPositionOP build_position( new BuildPosition( it->first, full_tl, it->second ) );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

void
LigandMotifSearch::defs2BuildPositions_findts(
	Pose const & pose,
	DnaDesignDefOPs const & defs
)
{
	if ( protein_dna_motif() ) {
		Sizes full_tl( map2keyvector( target_positions_ ) );
		std::map< Size, std::set< std::string > > mappositions( defs2map( pose, defs ) );
		for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
					end( mappositions.end() ); it != end; ++it ) {
			Size test(it->first);
			Sizes short_tl( shorten_target_list( pose, test, full_tl ) );
			BuildPositionOP build_position( new BuildPosition( it->first, short_tl, it->second ) );
			build_positionOPs_.push_back( build_position );
		}
	} else {
		ms_tr << "ERROR! These motifs are not protein-DNA, need to add another else if statement that allows for a different type of interface finding function." << std::endl;
	}
}

utility::vector1< core::Size >
LigandMotifSearch::map2keyvector(
	std::map< Size, std::set< std::string > > mappositions
)
{
	Sizes positions(0);
	for ( std::map<Size, std::set< std::string > >::const_iterator it( mappositions.begin() ),
				end( mappositions.end() ); it != end; ++it ) {
		positions.push_back( it->first );
	}
	return positions;
}

utility::vector1< core::Size >
LigandMotifSearch::shorten_target_list(
	Pose const & pose,
	Size const bp,
	Sizes & full_tl
)
{
	Sizes short_tl(0);
	Sizes bps(0);
	bps.push_back( bp );
	protein_DNA_motif_target_positions_JA( pose, bps, full_tl, short_tl );
	return short_tl;
}

void
LigandMotifSearch::protein_DNA_motif_build_positions_JA(
	Pose const & pose,
	Sizes & build_positions,
	Sizes & target_positions
)
{
	protocols::dna::DnaInterfaceFinderOP interface( new protocols::dna::DnaInterfaceFinder( 10*10, 3.9*3.9, 6., true ) );  //how does the z_axis cutoff work here that JA uses, I looked up once, but need ot read again
	if ( ! target_positions.empty() ) { // again, this won't work well if there are multiple target positions in vector
		interface->determine_protein_interface( pose, protein_positions_, target_positions ); // unless target_positions_ is empty - will actually deal with that later, will fill with all DNA in initialize if no protein pos or dna pos are given
		protocols::dna::DnaNeighbors protein_neighbors = interface->protein_neighbors();
		for ( protocols::dna::DnaNeighbors::const_iterator itr( protein_neighbors.begin() ),
					end( protein_neighbors.end() ); itr != end; ++itr ) {
			if ( (*itr).second.contact() ) {
				build_positions.push_back( itr->first );
				ms_tr << "Positions being targeted for motif design " << itr->first << std::endl;
			}
		}
//		ms_tr << "Attempting to identify build positions when there are no target positions." << std::endl;
	}
}

void
LigandMotifSearch::protein_DNA_motif_target_positions_JA(
	Pose const & pose,
	Sizes & build_positions,
	Sizes & target_positions,
	Sizes & short_tl
)
{
	// JA used 3.7, I picked 3.9 to access a position I knew was important
	protocols::dna::DnaInterfaceFinderOP interface( new protocols::dna::DnaInterfaceFinder( 10*10, 3.9*3.9, 10., true ) );  //how does the z_axis cutoff work here that JA uses, I looked up once, but need ot read again
	if ( ! build_positions.empty() ) {
		interface->determine_dna_interface( pose, build_positions, target_positions );
		protocols::dna::DnaNeighbors dna_neighbors = interface->dna_neighbors();
		for ( protocols::dna::DnaNeighbors::const_iterator itr( dna_neighbors.begin() ),
					end( dna_neighbors.end() ); itr != end; ++itr ) {
			if ( (*itr).second.contact() ) {
				short_tl.push_back( itr->first );
				ms_tr << "Positions (DNA) being targeted for motif design " << itr->first << std::endl;
			}
		}
//		ms_tr << "Attempting to identify build positions when there are no target positions." << std::endl;
	}
}

void
LigandMotifSearch::override_option_input(
	Real const & r1,
	Real const & z1,
	Real const & r2,
	Real const & z2,
	Real const & d1,
	Size const & rlevel
)
{
	ztest_cutoff_1_ = z1;
	ztest_cutoff_2_ = z2;
	rmsd_cutoff_1_ = r1;
	rmsd_cutoff_2_ = r2;
	dtest_cutoff_ = d1;
	rot_level_ = rlevel;
}

void
LigandMotifSearch::reset_option_input()
{
	ztest_cutoff_1_ = basic::options::option[ basic::options::OptionKeys::motifs::z1 ]();
	ztest_cutoff_2_ =	basic::options::option[ basic::options::OptionKeys::motifs::z2 ]();
	rmsd_cutoff_1_ =	basic::options::option[ basic::options::OptionKeys::motifs::r1 ]();
	rmsd_cutoff_2_ =	basic::options::option[ basic::options::OptionKeys::motifs::r2 ]();
	dtest_cutoff_ =	basic::options::option[ basic::options::OptionKeys::motifs::dtest ]();
	rot_level_ =	basic::options::option[ basic::options::OptionKeys::motifs::rotlevel ]();
}

void
LigandMotifSearch::set_motif_library(
	MotifLibrary & motiflibrary
)
{
	for( protocols::motifs::MotifCOPs::const_iterator motifcop_itr = motiflibrary.begin(), end_itr = motiflibrary.end();
			 motifcop_itr != end_itr; ++motifcop_itr ) {
		protocols::motifs::MotifCOP motifcop( *motifcop_itr );
		motif_library_.push_back( motifcop );
	}
}

// should probably be a private fxn, so should other ones . . . need to organize code
void
LigandMotifSearch::init_options()
{
	if( basic::options::option[ basic::options::OptionKeys::motifs::BPData ].user() ) {
		bpdata_ = true;
		bpdata_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::BPData ]();
	} else {
		bpdata_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::output_file ].user() ) {
		output_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::output_file ]();
		output_ = true;
	} else {
		output_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::data_file ].user() ) {
		data_filename_ = basic::options::option[ basic::options::OptionKeys::motifs::data_file ]();
		data_ = true;
	} else {
		data_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::quick_and_dirty ]).user() ) {
		quick_and_dirty_ = true;
	} else {
		quick_and_dirty_ = false;
	}
	if( (basic::options::option[ basic::options::OptionKeys::motifs::dump_motifs ]).user() ) {
		dump_motifs_ = true;
	} else {
		dump_motifs_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::clear_bprots ].user() ) {
		clear_bprots_ = true;
	} else {
		clear_bprots_ = false;
	}
	if( basic::options::option[ basic::options::OptionKeys::motifs::rots2add ].user() ) {
		rots2add_ = basic::options::option[ basic::options::OptionKeys::motifs::rots2add ]();
	} else {
		rots2add_ = 100;
	}
}

} // motifs
} // protocols
