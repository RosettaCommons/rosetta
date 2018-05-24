// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/requirements/LigandBindingResPlacer.hh
/// @brief Using the file format set in zinc_statistic_generator, checks for residues that could coordinate a metal.
/// @author guffysl (guffy@email.unc.edu)

#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
#include <protocols/sewing/data_storage/SmartAssembly.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/data_storage/LigandSegment.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/kinematics/Stub.hh>

#include <basic/Tracer.hh>

#include <numeric/HomogeneousTransform.hh>
#include <numeric/constants.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/io/izstream.hh>


static basic::Tracer TR( "protocols.sewing.hashing.LigandBindingResPlacer" );

namespace protocols {
namespace sewing {
namespace hashing {



///@brief Default constructor
LigandBindingResPlacer::LigandBindingResPlacer(){
	atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
}

///@brief Constructor with arguments
LigandBindingResPlacer::LigandBindingResPlacer( utility::pointer::shared_ptr< std::map< char, LigandHashMap > > coords, data_storage::LigandResidueOP ligand_residue ){
	ligand_ = ligand_residue;
	ligand_local_coords_ = coords;
	//Calculate contacts_to_make_
	core::Size desired_contacts = 0;
	for ( std::pair< core::Size, hashing::IdealContact > ideal: ligand_->get_ideal_contacts() ) {
		desired_contacts += ideal.second.max_contacts;
	}
	//core::Size current_contacts = ligand_->get_current_contacts().size();
	core::Size current_contacts = 0;
	for ( data_storage::LigandContactOP contact: ligand_->get_current_contacts() ) {
		if ( contact->segment_id != 0 || (contact->segment_id == 0 && contact->residue_number != 0 ) ) { //Second part refers to partner contacts
			++current_contacts;
		}
	}
	contacts_to_make_ = desired_contacts - current_contacts;
	atom_types_ = core::chemical::ChemicalManager::get_instance()->atom_type_set( "fa_standard" );
}




///@brief Parses the ligand file format
//TODO: Make this a static function that returns the data. Store the data in LigandBindingAssemblyMover (a map from LigandID to this data) and just pass a COP for it to the LBRP
void
LigandBindingResPlacer::parse_ligand_files( utility::vector1< data_storage::LigandDescription > & ligands ){
	///@details
	//The form it should be parsed into should fill the following criteria:
	//1) Organized by secondary structure
	//2) Organized by residue to avoid unnecessary checks
	//3) Ligand coordinates hashed for fast comparison using numeric::hash_value( xyzVector ) for first step; stored as vector for second step
	//4) We won't need chi angles until the second step

	///First format: Map of (secstruct, vector< LigandCoordInfo> )
	///Second format: Map of secstruct to a std::set of zinc HashKeys (generated the same basic way we do in Hasher)


	for ( data_storage::LigandDescription & ligand: ligands ) {


		utility::pointer::shared_ptr< std::map< char, LigandHashMap > > output( new std::map< char, LigandHashMap > );

		//Iterate over the files
		for ( std::string file_name: ligand.ligand_coord_files ) {

			//Open the file
			utility::io::izstream ligand_file( file_name, std::ios_base::in );
			if ( !ligand_file.good() ) {
				utility_exit_with_message( "File " + file_name + " invalid!" );
			}
			std::string line;
			//Get the DSSP char for that file (should be the first line of the file )
			ligand_file.getline( line );
			utility::trim( line, "\t\n ");
			char dssp = line.at( 0 );

			//Check that the ligand res name (second line of the file ) matches the amino acid name of ligand_. If not, exit with an error.
			ligand_file.getline( line );
			utility::trim( line, "\t\n ");
			TR.Debug << "File content: " << line << " " << dssp << std::endl;

			//Skip the third line (column labels)
			ligand_file.getline( line );

			//Iterate over the lines
			while ( ligand_file.getline( line ) ) {
				//Separate each line into fields "Coord_res_name\tCoord_atom_number\tLigand_atom_number\tLocal_x\tLocal_y\tLocal_z\tChi1\tChi2\tChi3\tChi4"
				utility::vector1< std::string > fields = utility::string_split_multi_delim( line, " \t" );
				numeric::xyzVector< core::Real > local_coords( utility::string2Real( fields.at( 4 ) ), utility::string2Real( fields.at( 5 ) ), utility::string2Real( fields.at( 6 ) ) );
				utility::fixedsizearray1< core::Real, 4 > chis;
				chis[ 1 ] = utility::string2Real( fields.at( 7 ) );
				chis[ 2 ] = utility::string2Real( fields.at( 8 ) );
				chis[ 3 ] = utility::string2Real( fields.at( 9 ) );
				chis[ 4 ] = utility::string2Real( fields.at( 10 ) );
				//Create a LigandCoordInfo from the fields
				LigandCoordInfo coord;
				coord.local_coords = local_coords;
				coord.chis = chis;
				coord.coord_res_name = fields.at( 1 );
				coord.coord_atom_index = utility::string2Size( fields.at( 2 ) );
				coord.ligand_atom_index = utility::string2Size( fields.at( 3 ) );

				//Convert the coordinates into a HashKey
				HashKey key = generate_ligand_key( local_coords );

				//Populate the entry in ligand_local_coords_ for that line
				(*output)[ dssp ][ key ].insert( coord );

			} //End iterate over lines
			//Close file
			ligand_file.close();
		} //End iterate over files
		ligand.ligand_coords = output;
	}//end iterate over ligands
}



core::conformation::ResidueOP
LigandBindingResPlacer::make_aligned_residue( data_storage::SmartSewingResidueOP base_res, std::string type ){
	//Make the (mobile) residue
	core::chemical::ResidueTypeSetCOP res_type_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	TR << "Generating residue of type " << type << std::endl;
	core::conformation::ResidueOP mobile_res( new core::conformation::Residue( res_type_set->name_mapOP( type), true ) );

	//Make the homogeneous transform for the mobile residue
	utility::vector1< core::conformation::Atom > & mobile_atoms = mobile_res->atoms();
	numeric::HomogeneousTransform< core::Real > mobile_ht( mobile_atoms[ 3 ].xyz(), mobile_atoms[ 1 ].xyz(), mobile_atoms[ 2 ].xyz() );

	//Make the homogeneous transform for the stationary residue
	utility::vector1< core::conformation::Atom > & stationary_atoms = base_res->get_atom_vector();
	numeric::HomogeneousTransform< core::Real > stationary_ht( stationary_atoms[ 3 ].xyz(), stationary_atoms[ 1 ].xyz(), stationary_atoms[ 2 ].xyz() );

	//Make the finalized HTs that we will apply to go from one to the other
	numeric::HomogeneousTransform< core::Real > inverse_mobile_ht = mobile_ht.inverse();
	numeric::HomogeneousTransform< core::Real > mobile_to_stationary_ht = stationary_ht * inverse_mobile_ht;

	//Apply the HTs

	for ( core::conformation::Atom & current_atom: mobile_atoms ) {
		numeric::xyzVector< core::Real > old_coords = current_atom.xyz();
		numeric::xyzVector< core::Real > new_coords = mobile_to_stationary_ht * old_coords;
		current_atom.xyz( new_coords );
	}
	return mobile_res;
}



///@brief Checks all residue stubs in assembly vs metal coordinates to determine which residues are potentially metal-coordinating. Return a set of SmartSewingResidueOP.
//std::set< data_storage::SmartSewingResidueOP >
std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> >
LigandBindingResPlacer::identify_possible_binders( data_storage::SmartAssemblyCOP assembly ){
	///@details Iterate over all segments
	// << "Begin identify_possible_binders" << std::endl;
	std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> > coord_per_residue;

	//We don't want to test the whole assembly. Instead, we will look at the new segments starting at the appropriate terminus and going until we reach the chimaeric segment's basis residue.


	data_storage::SmartSegmentOP current_segment;
	if ( assembly->get_last_change_was_n_terminal() ) {
		current_segment = assembly->get_n_terminal_segment();
	} else {
		current_segment = assembly->get_c_terminal_segment();
	}
	numeric::xyzVector< core::Real > origin;
	origin.x( 0);
	origin.y( 0);
	origin.z( 0);
	bool done = false;
	core::Size stopping_point = 1;
	while ( !done && current_segment != nullptr ) { //End when we either pass the chimaera or get to the end
		///Iterate over all residues in the segment
		if ( current_segment->is_chimaeric() ) {
			done = true; //We'll be done after this segment
			//Figure out what the breaking point will be
			//If this is a *double* chimaera, things will be a little weird, but only if the parent receiving the ligand residue is the chimaeric parent
			if ( assembly->get_last_change_was_n_terminal() ) {
				//Simple--we'll go until we reach the second basis residue (still in same order)
				stopping_point = current_segment->get_basis_pair().second.resnum();
			} else {
				//Less simple--go BACKWARDS until we reach the first basis residue OR just skip all residues <= that one
				stopping_point = current_segment->get_basis_pair().first.resnum();
			}
		}
		//We'll actually want to leave room for the window width! I've been seeing lots of sites with histidines at the very end of the helix--we don't want that.
		//if we're close to the n terminus, we don't want to start until we're at residue # window_width + 1
		//We also don't want to get too close to the C terminus--leave at least window_width residues after

		if ( current_segment->get_length() <= assembly->get_window_width() ) {
			if ( assembly->get_last_change_was_n_terminal() ) {
				current_segment = current_segment->get_c_terminal_neighbor();
			} else {
				current_segment = current_segment->get_n_terminal_neighbor();
			}
			continue;
		}
		for ( core::Size i = assembly->get_window_width() + 1; i <= ( current_segment->get_length() - assembly->get_window_width() ); ++i ) {
			if ( done && assembly->get_last_change_was_n_terminal() && i >= stopping_point ) {
				break; //We're done
			}
			if ( done && !assembly->get_last_change_was_n_terminal() && i <= stopping_point ) {
				continue; // Move on to the next residue
			}
			//Skip SS elements we aren't considering (e.g. loops)
			if ( ligand_local_coords_->count( current_segment->get_dssp_code() ) == 0 ) {
				continue;
			}
			///If it's a required residue, skip it
			std::set< core::Size > vital = current_segment->get_vital_residues();
			if ( vital.find( i ) != vital.end() ) {
				continue;
			}
			///Otherwise, convert it to a stub
			data_storage::SmartSewingResidueOP res = current_segment->get_residue( i );
			//TR << "Checking segment " << current_segment->get_segment_id() << " DSSP " << current_segment->get_dssp_code() << " residue " << i << std::endl;
			core::kinematics::Stub res_stub( res->get_atom( 2 ).xyz(), res->get_atom( 1 ).xyz(), res->get_atom( 2 ).xyz(), res->get_atom( 3 ).xyz() );
			//We'll need to check all the atoms against the LigandHashMap--iterate over atoms
			core::Size count = 1;
			for ( core::conformation::Atom const & ligat: ligand_->get_atom_vector() ) {
				if ( ( *atom_types_.lock() )[ ligat.type() ].atom_type_name() == "VIRT" ) {
					continue;
				}
				//     if( ligat.type() == 147 ) continue; //These are virtual atoms
				///Transform the ligand's coordinates into the stub's local coordinates
				numeric::xyzVector< core::Real > transformed_coords = res_stub.global2local( ligat.xyz() );
				//TR << "transformed zinc coordinates: " << transformed_coords.x() << " " << transformed_coords.y() << " " << transformed_coords.z() << std::endl;
				//TEMP
				//This will slow things down some but is being used for debugging
				core::Real distance = transformed_coords.distance( origin );
				if ( distance >= 8 ) continue;
				//END TEMP
				//We'll actually test 27 different hash keys ( 3 in each direction)
				HashKey main_key = generate_ligand_key( transformed_coords );
				HashKey current_key;
				for ( int delta_x = -1; delta_x <= 1; ++delta_x ) {
					current_key[ 1 ] = main_key[ 1 ] + delta_x;
					for ( int delta_y = -1; delta_y <= 1; ++delta_y ) {
						current_key[ 2 ] = main_key[ 2 ] + delta_y;
						for ( int delta_z = -1; delta_z <= 1; ++delta_z ) {
							current_key[ 3 ] = main_key[ 3 ] + delta_z;
							if ( ligand_local_coords_->at( current_segment->get_dssp_code() ).find( current_key ) == ligand_local_coords_->at( current_segment->get_dssp_code() ).end()  ) {
								continue;
							}
							///Compare the ligand atom name to see if it's a real match
							//std::map< char, LigandHashMap > ligand_local_coords_;
							//typedef boost::unordered_map< HashKey, std::set< LigandCoordInfo >, coord_hash, coord_equal_to > LigandHashMap;
							//TR << "Found a match!" << std::endl;
							TR.Debug << "Match in segment " << current_segment->get_segment_id() << " residue " << i << std::endl;
							///If it's there, add this residue to the map with the corresponding entries from LigandCoordInfo
							for ( LigandCoordInfo coord: ligand_local_coords_->at( current_segment->get_dssp_code() ).at( current_key ) ) {
								if ( coord.ligand_atom_index != count ) {
									//TR << "Atom names do not match!" << std::endl;
									continue;
								}
								//TR << "Adding match to coord_per_residue" << std::endl;
								( coord_per_residue[ std::make_pair( current_segment->get_segment_id(), i  )] ).insert( coord );
							}//end iterate over coord info
						}//end deltaz
					}//end deltay
				}//end deltax
				++count;
			}//End iterate over ligand atoms
		}//End iterate over residues
		if ( assembly->get_last_change_was_n_terminal() ) {
			current_segment = current_segment->get_c_terminal_neighbor();
		} else {
			current_segment = current_segment->get_n_terminal_neighbor();
		}
	}//End iterate over segments
	return coord_per_residue;
}

///@brief Finds the best residue/rotamer to coordinate metal (Main function for class)
std::pair< data_storage::SmartSegmentOP, bool >
LigandBindingResPlacer::choose_best_metal_coordinator( data_storage::SmartAssemblyOP assembly ){
	//TR << "Begin choose_best_metal_coordinator" << std::endl;
	TR.Debug << "Trying to make " << contacts_to_make_ << " contacts" << std::endl;
	///@details First identify possible metal binding residues
	if ( contacts_to_make_ == 0 ) {
		return std::make_pair( nullptr, false );
	}
	possible_coordinations_.clear();
	//This is a vector of pairs instead of a map so we can sort it based on scores
	utility::vector1< std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > > > res_scores; //This is ugly
	//The first pair is (segment ID, resnum ) so we know which residue to change
	//The second pair is (LigandCoordInfo, score)--LigandCoordInfo tells us the residue and rotamer to use, and the score is used for sorting
	possible_coordinations_ = identify_possible_binders( assembly );
	//WHY ARE THESE NOT RUNNING?????
	TR.Debug << "Done finding possible coordinations. Size: " << possible_coordinations_.size();
	if ( possible_coordinations_.size() == 0 ) {
		//TR << "No contacts found!" << std::endl;
		return std::make_pair( nullptr, false );
	}
	///Iterate over all residues in the map
	for ( std::pair< std::pair< core::Size, core::Size>, std::set< LigandCoordInfo > > res_coord: possible_coordinations_ ) {
		//TR << "Scoring match!" << std::endl;
		///Find best coordinating rotamer for each residue per restype
		res_scores.push_back(  std::make_pair( res_coord.first, best_rotamer_for_residue( res_coord.first, res_coord.second, assembly ) ) );
	}
	if ( res_scores.size() == 0 ) {
		return std::make_pair( nullptr, false );
	}
	//Sort the vector by score
	//std::sort( res_scores.begin(), res_scores.end(), boost::bind( &std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > >::second::second, _1) < boost::bind(  &std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > >::second::second, _2 ) );
	//TR << "Sorting vector by score" << std::endl;
	//see if this lambda function will work
	std::sort( res_scores.begin(), res_scores.end(), []( std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > > const &left, std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > > const &right ){
		return left.second.second < right.second.second;
	});

	//Check that we actually found one that works
	if ( res_scores.at( 1 ).second.second >= geometry_score_weight_ ) {
		TR << "None of the matches passed the score threshold!" << std::endl;
		TR << "Top failed score: " << res_scores.at( 1 ).second.second << std::endl;
		if ( res_scores.size() > 1 ) {
			TR << "Next best: " << res_scores.at( 2 ).second.second << std::endl;
		}
		return std::make_pair( nullptr, false );
	}
	TR << "Found contact! Making mutation." << std::endl;
	//At this point we should just be able to make the mutation at the top of the list, add that info to the ligand contacts and to the segment (as an unowned contact)
	data_storage::SmartSewingResidueOP res_to_change = assembly->get_segment( res_scores.at( 1 ).first.first )->get_residue( res_scores.at( 1 ).first.second  );
	core::conformation::ResidueOP reference_res = make_aligned_residue( res_to_change, res_scores.at( 1 ).second.first.coord_res_name );
	for ( core::Size i = 1; i <= reference_res->nchi(); ++i ) {
		reference_res->set_chi( i, res_scores.at( 1 ).second.first.chis[ i ] );
	}

	//Now we'll need to change the SmartSewingResidue
	//Set the atom vector
	res_to_change->set_atom_vector( reference_res->atoms() );
	//Set the amino acid type
	res_to_change->set_amino_acid_type( reference_res->type().base_name() );
	res_to_change->set_full_type_name( reference_res->type().name() );
	//Set the chis??
	res_to_change->set_chis( reference_res->chi() );

	//That segment may not/ probably will not be a ligand segment, though, so if it is not then we will need to replace it in the assembly
	data_storage::SmartSegmentOP chosen_seg = assembly->get_segment( res_scores.at( 1 ).first.first );
	if ( std::dynamic_pointer_cast< data_storage::LigandSegment >( chosen_seg ) != nullptr ) {
		//We won't need to replace the segment--just add the new contact and ligand
		data_storage::LigandSegmentOP chosen_ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment >( chosen_seg );
		chosen_ligseg->set_is_vital( true );
		chosen_ligseg->attach_ligand( ligand_, false );
		chosen_ligseg->add_ligand_contact( res_scores.at( 1 ).first.second );
		//Note that it's still entirely possible that this ligseg is a chimaera, and so we'll need to modify its parent(s) as well
		if ( chosen_seg->is_chimaeric() ) {
			bool change_n_terminal_parent = assembly->get_last_change_was_n_terminal();
			bool first_basis_is_n_terminal;
			//In this case, we know the last change, so it's easy to know which basis is which
			if ( change_n_terminal_parent ) {
				first_basis_is_n_terminal = false;
			} else {
				first_basis_is_n_terminal = true;
			}
			replace_chimaera_ligand_parent( assembly, chosen_seg, chosen_ligseg, change_n_terminal_parent, first_basis_is_n_terminal, res_scores.at( 1 ).first.second, res_scores );
		}
		//Just in case
		assembly->local_segments()[ chosen_seg->get_segment_id() ] = chosen_seg;
	} else {
		std::set< core::Size > contact_res;
		contact_res.insert( res_scores.at( 1 ).first.second );
		data_storage::LigandSegmentOP replacement_ligseg( new data_storage::LigandSegment( *chosen_seg, true, contact_res, ligand_->get_ligand_id() ) );

		//NOTE!!!! Think about giving these replacements new segIDs and adding them to pdb_segments if problems persist
		//Doesn't look like it should be an issue, but stay alert
		replacement_ligseg->set_is_in_Assembly( true );

		//Now to replace the old segment
		//We can directly modify local_segments, but we'll do that last
		//First replace linkages from the old segment's neighbors
		if ( chosen_seg->get_n_terminal_neighbor() != nullptr ) {
			data_storage::SmartSegment::link_to( chosen_seg->get_n_terminal_neighbor(), replacement_ligseg );
		}
		if ( chosen_seg->get_c_terminal_neighbor() != nullptr ) {
			data_storage::SmartSegment::link_to( replacement_ligseg, chosen_seg->get_c_terminal_neighbor() );
		}
		//If this segment is a chimaera, we'll need to rewrite history just a little bit and turn one of its parents into a ligand segment
		//The parent that will need to become a ligand segment is the one that was just added
		if ( chosen_seg->is_chimaeric() ) {
			//TODO
			bool change_n_terminal_parent = assembly->get_last_change_was_n_terminal();
			bool first_basis_is_n_terminal;
			//In this case, we know the last change, so it's easy to know which basis is which
			if ( change_n_terminal_parent ) {
				first_basis_is_n_terminal = false;
			} else {
				first_basis_is_n_terminal = true;
			}
			//This will also recursively handle multichimaerae
			core::Size chimaeric_contact_resnum = res_scores.at( 1 ).first.second;
			replace_chimaera_ligand_parent( assembly, chosen_seg, replacement_ligseg, change_n_terminal_parent, first_basis_is_n_terminal, chimaeric_contact_resnum, res_scores );
		}
		//Replace in local_segments
		assembly->local_segments()[ chosen_seg->get_segment_id() ] = data_storage::SmartSegmentOP( replacement_ligseg );
		delete chosen_seg.get();
		chosen_seg = data_storage::SmartSegmentOP( replacement_ligseg );
	}
	//Add the new LigandContact to the ligand
	data_storage::LigandContactOP contact( new data_storage::LigandContact( res_scores.at( 1 ).first.first, res_scores.at( 1 ).first.second, res_scores.at( 1 ).second.first.coord_atom_index, res_scores.at( 1 ).second.first.ligand_atom_index ) );
	best_geometry_score_ = res_scores.at( 1 ).second.second;
	ligand_->add_contact( contact );
	last_contact_added_ = contact;
	--contacts_to_make_;
	//delete chosen_seg;
	return std::make_pair( chosen_seg, true );
}

HashKey
LigandBindingResPlacer::generate_ligand_key( numeric::xyzVector< core::Real > coords){
	HashKey key;
	/*
	key[1] = boost::math::iround( coords.x()*4 );
	key[2] = boost::math::iround( coords.y()*4 );
	key[3] = boost::math::iround( coords.z()*4 );
	*/
	key[1] = std::round( coords.x()*4 );
	key[2] = std::round( coords.y()*4 );
	key[3] = std::round( coords.z()*4 );
	return key;
}

///@brief Finds the best residue/rotamer to coordinate metal for a given residue and returns for each residue type a pair of (LigandCoordInfo, distance) where distance is the distance of the metal from its ideal position
std::pair< LigandCoordInfo, core::Real >
LigandBindingResPlacer::best_rotamer_for_residue( std::pair< core::Size, core::Size > res, std::set< LigandCoordInfo > res_coord, data_storage::SmartAssemblyOP assembly ){
	//At this point we've already created ligand_local_coords. We'll need to also know the DSSP of this residue, so maybe we should just take a SmartSegment and residue number?
	///@details For the given residue, do the following:
	///Iterate over LigandCoordInfo and for each one set the chi angles
	//TR << "Begin best_rotamer_for_residue" << std::endl;
	//Initialize return value
	LigandCoordInfo dummy;
	std::pair< LigandCoordInfo, core::Real > best = std::make_pair( dummy, 10e12 );
	core::chemical::ResidueTypeSetCOP res_type_set =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	//Get ligand contacts
	utility::vector1< data_storage::LigandContactOP > lig_contacts = ligand_->get_current_contacts();

	core::conformation::ResidueOP test_residue;
	std::string current_type = "";

	for ( LigandCoordInfo coord: res_coord ) {
		//If the residue type hasn't changed, no need to make a new residue
		//If it has, make a new one and orient the backbone coordinates
		if ( coord.coord_res_name != current_type ) {
			///Make a residue of the appropriate type and transform it onto res
			test_residue = make_aligned_residue( assembly->get_segment( res.first )->get_residue( res.second ), coord.coord_res_name );
			current_type = coord.coord_res_name;
		}
		for ( core::Size i = 1; i <= test_residue->nchi(); ++i ) {
			test_residue->set_chi( i, coord.chis[ i ] );
		}
		if ( ligand_->get_ideal_contacts().count( coord.ligand_atom_index ) == 0 ) {
			continue;
		}

		IdealContact goal_values = ligand_->get_ideal_contacts().at( coord.ligand_atom_index );

		//Number of contacts
		core::Size num_contacts = 1;
		utility::vector1< data_storage::LigandContactOP > other_atom_contacts;
		for ( data_storage::LigandContactOP ligcon: lig_contacts ) {
			if ( ligcon->ligand_atom == coord.ligand_atom_index ) {
				++ num_contacts;
				other_atom_contacts.push_back( ligcon );
			}
		}
		if ( num_contacts > goal_values.max_contacts ) {
			continue; //This contact would exceed the maximum number of contacts for this atom, so skip it
		}

		//Now check for clashes in this rotamer
		//TEMPORARILY COMMENTED OUT FOR DEBUGGING

		bool clash = false;
		//Iterate over atoms in test residue
		//for( core::conformation::Atom test_at: test_residue->atoms() ){
		for ( core::Size test_i = 1; test_i <= test_residue->atoms().size(); ++test_i ) {
			if ( clash ) {
				break;
			}
			core::conformation::Atom test_at = test_residue->atoms()[ test_i ];

			if ( ( *atom_types_.lock() )[ test_at.type() ].atom_type_name() == "VIRT" ) {
				continue;
			}
			//Get the LJ radius
			core::Real test_at_LJ = ( *atom_types_.lock() )[ test_at.type() ].lj_radius();
			numeric::xyzVector< core::Real > test_at_xyz = test_at.xyz();
			//Iterate over other residues that are in contact with the ligand
			std::string test_at_element = ( *atom_types_.lock() )[ test_at.type() ].element();

			for ( data_storage::LigandContactOP ligcon: other_atom_contacts ) {
				if ( clash ) {
					break;
				}

				//If this atom is in the ligand, we'll skip it since those are tested later
				if ( ligcon->segment_id == 0 && ligcon->residue_number == 0 ) {
					continue;
				} else if ( ligcon->segment_id == 0 ) {
					//If this is a partner contact, it will need to be treated specially
					//Partner contact--contacts are Rosetta residues, not SSRes
					core::conformation::Residue const & contact_res = assembly->get_partner()->residue( ligcon->residue_number );
					//Iterate over atoms
					for ( core::Size contact_res_i = 1; contact_res_i <= contact_res.natoms(); ++contact_res_i ) {
						//Get atom
						core::conformation::Atom contact_res_at = contact_res.atom( contact_res_i );
						if ( ( *atom_types_.lock() )[ contact_res_at.type() ].atom_type_name() == "VIRT" ) {
							continue;
						}
						core::Real contact_res_at_LJ = ( *atom_types_.lock() )[ contact_res_at.type() ].lj_radius();
						numeric::xyzVector< core::Real > contact_res_at_xyz = contact_res_at.xyz();
						core::Real max_ok_distance = 0.8 * ( test_at_LJ + contact_res_at_LJ );
						std::string contact_res_at_element = ( *atom_types_.lock() )[ contact_res_at.type() ].element();
						if ( test_at_xyz.distance( contact_res_at_xyz ) < max_ok_distance ) {
							TR << "Possible clash detected with coordinating res: " << std::endl;
							TR << "Test residue = segment " << res.first << " residue " << res.second << " atom " << test_residue->atom_name( test_i ) << " " << test_at.xyz().x() << " " <<  test_at.xyz().y() << " " <<  test_at.xyz().z() << std::endl;
							TR << "Clashes with segment " << ligcon->segment_id << " residue " << ligcon->residue_number << " atom " << contact_res.type().atom_name( contact_res_i ) << " " << contact_res_at.xyz().x() << " " << contact_res_at.xyz().y() << " " << contact_res_at.xyz().z() << " " << std::endl;
							clash = true;
							break;
						}
					}//end iterate over atoms
				} else { //end if
					//Retrieve the SSResidue for this contact
					data_storage::SmartSewingResidueOP ssres = assembly->get_segment( ligcon->segment_id )->get_residue( ligcon->residue_number );
					//Iterate over atoms in ssres
					//for( core::conformation::Atom ssres_at: ssres->get_const_atom_vector() ){
					for ( core::Size ssres_i = 1; ssres_i <= ssres->get_const_atom_vector().size(); ++ssres_i ) {
						core::conformation::Atom ssres_at = ssres->get_const_atom_vector()[ ssres_i ];
						if ( ( *atom_types_.lock() )[ ssres_at.type() ].atom_type_name() == "VIRT" ) {
							continue;
						}
						//Get the LJ radius
						core::Real ssres_at_LJ = ( *atom_types_.lock() )[ ssres_at.type() ].lj_radius();
						numeric::xyzVector< core::Real > ssres_at_xyz = ssres_at.xyz();
						core::Real max_ok_distance = 0.8 * (test_at_LJ + ssres_at_LJ );
						//TEMP////
						std::string ssres_at_element = ( *atom_types_.lock() )[ ssres_at.type() ].element();
						if ( test_at_xyz.distance( ssres_at_xyz ) < max_ok_distance ) {
							TR << "Possible clash detected with coordinating res: " << std::endl;
							TR << "Test residue = segment " << res.first << " residue " << res.second << " atom " << test_residue->atom_name( test_i ) << " " << test_at.xyz().x() << " " <<  test_at.xyz().y() << " " <<  test_at.xyz().z() << std::endl;
							TR << "Clashes with segment " << ligcon->segment_id << " residue " << ligcon->residue_number << " atom " << res_type_set->name_map( ssres->get_full_type_name() ).atom_name( ssres_i ) << " " << ssres_at.xyz().x() << " " << ssres_at.xyz().y() << " " << ssres_at.xyz().z() << " " << std::endl;
							clash = true;
							break;
						}//end if clash
					}///end iterate over atoms
					//end else
				}
			}//end iterate over residues

			//We should also check for clashes of non-binding atoms with the ligand

			if ( test_residue->atom( coord.coord_atom_index ).xyz() != test_at.xyz() ) {
				//So if this atom is not the coordinating atom
				//Iterate over ligand atoms
				for ( core::Size ligand_i = 1; ligand_i <= ligand_->get_atom_vector().size(); ++ligand_i ) {
					core::conformation::Atom ligand_at = ligand_->get_atom_vector()[ ligand_i ];
					//for( core::conformation::Atom ligand_at: ligand_->get_atom_vector() ){
					if ( ( *atom_types_.lock() )[ ligand_at.type() ].atom_type_name() == "VIRT" ) {
						continue;
					}
					std::string ligat_element = ( *atom_types_.lock() )[ ligand_at.type() ].element();
					numeric::xyzVector< core::Real > ligat_xyz = ligand_at.xyz();
					core::Real ligat_LJ = ( *atom_types_.lock() )[ ligand_at.type() ].lj_radius();
					core::Real max_ok_distance = 0.8 * ( test_at_LJ + ligat_LJ );

					if ( test_at_xyz.distance( ligat_xyz ) < max_ok_distance ) {
						TR << "Possible clash detected with ligand: " << std::endl;
						TR << "Test residue = segment " << res.first << " residue " << res.second << " atom " << test_residue->atom_name( test_i ) << " " << test_at.xyz().x() << " " <<  test_at.xyz().y() << " " <<  test_at.xyz().z() << std::endl;
						std::string ligand_atom_name = res_type_set->name_map( ligand_->get_full_type_name() ).atom_name( ligand_i );
						TR << "Clashes with ligand atom " << ligand_atom_name << " " << ligand_at.xyz().x() << " " << ligand_at.xyz().y() << " "  << ligand_at.xyz().z() << " " << std::endl;
						clash = true;
						break;
					}
				}
			}
		}//end iterate over test_at
		if ( clash ) {
			TR << "Found possible clash; skipping rotamer" << std::endl;
			continue;
		}



		//BEGIN SCORING

		core::Size coord_base = test_residue->atom_base( coord.coord_atom_index );
		core::Size coord_base_base = test_residue->atom_base( coord_base );

		core::Real score = 0;
		//Measure the distance between the ligand atom and the coordinating atom
		core::Real distance = test_residue->atom( coord.coord_atom_index ).xyz().distance( ligand_->get_atom( coord.ligand_atom_index ).xyz() );
		//Compare the distance to the optimal distance for this contact
		TR << "Comparing distance " << distance << " to goal value " << goal_values.preferred_distance << std::endl;
		score += ( distance - goal_values.preferred_distance ) * ( distance - goal_values.preferred_distance );

		//Angles between contacts
		//For these to be weighted properly vs the distance, these should be measured in radians, which is what this function returns anyway
		for ( data_storage::LigandContactOP ligcon: other_atom_contacts ) {
			//Measure angle between test_residue->atom( coord.coord_atom_index ).xyz(), ligand_->get_atom( coord.ligand_atom_index ).xyz(), and assembly.get_segment( ligcon->segment_id )->get_residue( ligcon->residue_number )->get_atom( ligcon->residue_atom ).xyz() )
			if ( ligcon->segment_id != 0 && ligcon->residue_number > assembly->get_segment( ligcon->segment_id )->get_length() ) {
				utility_exit_with_message( "Invalid access to ligand contacts!" );
			}
			numeric::xyzVector< core::Real > test_atom_coords = test_residue->atom( coord.coord_atom_index ).xyz();
			numeric::xyzVector< core::Real > ligand_atom_coords = ligand_->get_atom( coord.ligand_atom_index ).xyz();
			numeric::xyzVector< core::Real > other_atom_coords;
			if ( ligcon->segment_id != 0 ) {
				other_atom_coords = assembly->get_segment( ligcon->segment_id )->get_residue( ligcon->residue_number )->get_atom( ligcon->residue_atom ).xyz();
			} else if ( ligcon->residue_number == 0 ) { //Contact is internal to ligand
				other_atom_coords = ligand_->get_atom( ligcon->residue_atom ).xyz();
			} else { //Contact is in partner
				other_atom_coords = assembly->get_partner()->residue( ligcon->residue_number ).atom( ligcon->residue_atom ).xyz();
			}

			core::Real angle = numeric::angle_of( test_atom_coords, ligand_atom_coords, other_atom_coords );

			TR << "Comparing angle " << angle << " to ideal angle " << goal_values.preferred_contact_angle << std::endl;
			score += ( angle - goal_values.preferred_contact_angle ) * ( angle - goal_values.preferred_contact_angle ) * 10; //Upweight these substantially



			//Measure dihedral 1 about the ligand atom (test base - test - ligand atom - other contact atom )--provide in ideal contacts

			core::Real dihedral_1 = numeric::dihedral_radians( test_residue->atom( coord_base ).xyz(), test_atom_coords, ligand_atom_coords, other_atom_coords );
			TR << "Comparing dihedral " << dihedral_1 << " to ideal dihedral " << goal_values.preferred_dihedral_1 << std::endl;
			core::Real val_to_compare = std::min( std::fmod( dihedral_1, goal_values.preferred_dihedral_1 ), std::fabs( std::fmod( dihedral_1, goal_values.preferred_dihedral_1 ) - goal_values.preferred_dihedral_1 ) );
			score +=  val_to_compare * val_to_compare * 5;


			//Measure dihedral 2 about the ligand atom (test - ligand atom - other contact atom - other contact base )--provide in ideal contacts
			//This will be less straightforward since other contact isn't on a real residue
			numeric::xyzVector< core::Real > other_base_coords;

			if ( ligcon->segment_id != 0 && ligcon->residue_number != 0 ) {
				core::chemical::ResidueType const & other_type = res_type_set->name_map( assembly->get_segment( ligcon->segment_id )->get_residue( ligcon->residue_number )->get_full_type_name() );
				core::Size other_base = other_type.atom_base( ligcon->residue_atom );
				other_base_coords = assembly->get_segment( ligcon->segment_id )->get_residue( ligcon->residue_number )->get_atom( other_base ).xyz();
			} else if ( ligcon->residue_number != 0 ) { //Contact is in partner
				core::chemical::ResidueType const & other_type = assembly->get_partner()->residue( ligcon->residue_number ).type();
				core::Size other_base = other_type.atom_base( ligcon->residue_atom );
				other_base_coords = assembly->get_partner()->residue( ligcon->residue_number ).atom( other_base ).xyz();
			} else {
				core::chemical::ResidueType const & other_type = res_type_set->name_map( ligand_->get_amino_acid_type() );
				core::Size other_base = other_type.atom_base( ligcon->residue_atom );
				other_atom_coords = ligand_->get_atom( other_base ).xyz();
			}

			core::Real dihedral_2 = numeric::dihedral_radians( test_atom_coords, ligand_atom_coords, other_atom_coords, other_base_coords );
			TR << "Comparing dihedral " << dihedral_2 << " to ideal dihedral " << goal_values.preferred_dihedral_2 << std::endl;
			core::Real val2_to_compare = std::min( std::fmod( dihedral_2, goal_values.preferred_dihedral_2 ), std::fabs( std::fmod( dihedral_2, goal_values.preferred_dihedral_2 ) - goal_values.preferred_dihedral_2 ) );
			score += val2_to_compare * val2_to_compare * 5;

		}

		//Now we'll do something similar for both the his dihedral and the angle about the histidine. These won't be upweighted as much (we'll do 5x each)
		//ISSUE: Ideal angle and dihedral are dependent on the ligand residue. Therefore, we can't really score them directly.
		//Do residues hold any information about this? Atoms have properties--if they are aromatic or if this & base both sp2, ideal dihedral will be %180=0, else %30=0 (actual ideal value is rotamer dependent) . If they are SP2, ideal angle will be +-120, else +-109.5.

		//Determine whether coord atom and coord base are sp2/aromatic
		utility::vector1< std::string > coord_properties = ( *atom_types_.lock() )[ test_residue->atom( coord.coord_atom_index ).type() ].get_all_properties();
		utility::vector1< std::string > coord_base_properties = ( *atom_types_.lock() )[ test_residue->atom( coord_base ).type() ].get_all_properties();

		bool coord_is_sp2 = std::find( coord_properties.begin(), coord_properties.end(), "SP2_HYBRID" ) != coord_properties.end();
		bool coord_is_aromatic = std::find( coord_properties.begin(), coord_properties.end(), "AROMATIC" ) != coord_properties.end();
		bool coord_base_is_sp2 = std::find( coord_base_properties.begin(), coord_base_properties.end(), "SP2_HYBRID" ) != coord_base_properties.end();

		core::Real ideal_angle_with_residue;
		core::Real ideal_dihedral;
		if ( coord_is_aromatic || coord_is_sp2 ) {
			ideal_angle_with_residue = 120.0;
		} else {
			ideal_angle_with_residue = 109.5;
		}
		if ( coord_is_aromatic || ( coord_is_sp2 && coord_base_is_sp2 ) ) {
			ideal_dihedral = 180.0;
		} else {
			//We're actually just going to take modulus with these values, so I'll put 30 to indicate any %30
			ideal_dihedral = 30.0;
		}

		//Angle with residue
		core::Real angle_with_residue = numeric::angle_degrees( ligand_->get_atom( coord.ligand_atom_index ).xyz(), test_residue->atom( coord.coord_atom_index ).xyz(), test_residue->atom( coord_base ).xyz() );
		core::Real angleval_to_compare = std::min( std::fmod( angle_with_residue, ideal_angle_with_residue ), std::fabs( std::fmod( angle_with_residue, ideal_angle_with_residue ) - ideal_angle_with_residue ) );
		score += angleval_to_compare * angleval_to_compare * 5 * numeric::constants::r::degrees_to_radians * numeric::constants::r::degrees_to_radians;

		//Dihedral
		core::Real dihedral = numeric::dihedral_degrees( ligand_->get_atom( coord.ligand_atom_index ).xyz(), test_residue->atom( coord.coord_atom_index ).xyz(), test_residue->atom( coord_base).xyz(), test_residue->atom( coord_base_base ).xyz() );
		core::Real dival_to_compare = std::min( std::fmod( dihedral, ideal_dihedral ), std::fabs( std::fmod( dihedral, ideal_dihedral ) - ideal_dihedral ) );
		score += dival_to_compare * dival_to_compare * numeric::constants::r::degrees_to_radians * numeric::constants::r::degrees_to_radians;

		//If the score is better than the current best score, replace it
		TR << "Final score for match: " << score << std::endl;
		if ( score < best.second ) {
			best = std::make_pair( coord, score );
		}
	}//End iterate over LigandCoordInfo
	return best;
}

void
LigandBindingResPlacer::replace_chimaera_ligand_parent( data_storage::SmartAssemblyOP assembly, data_storage::SmartSegmentOP chosen_seg, data_storage::LigandSegmentOP replacement_ligseg, bool change_n_terminal_parent, bool first_basis_is_n_terminal, core::Size chimaeric_contact_resnum, utility::vector1< std::pair< std::pair< core::Size, core::Size >, std::pair< LigandCoordInfo, core::Real > > > res_scores ){
	//Figure out which parent should be replaced
	//If the last change was N terminal, then the chimaera's N-terminal segment will be the new one

	data_storage::SmartSegmentOP bad_parent;
	if ( change_n_terminal_parent ) {
		bad_parent = chosen_seg->get_n_terminal_parent();
	} else {
		bad_parent = chosen_seg->get_c_terminal_parent();
	}
	//Create the parent ligseg
	//The challenge here is figuring out what the contact residue number will be
	//Luckily I already wrote the code to figure it out in Assembly in reset_chimaera_contacts
	//If the residue came from the N-terminal segment, the residue number will be unchanged
	core::Size new_contact_resnum;
	if ( change_n_terminal_parent ) {
		new_contact_resnum = chimaeric_contact_resnum;
	} else {
		//We need to know the basis pair of chosen_seg to figure this out
		if ( first_basis_is_n_terminal ) {
			//We know the nterm basis will be 1 and cterm will be 2 because we added to the c terminus of the assembly
			new_contact_resnum = chosen_seg->get_basis_pair().second.resnum() + chimaeric_contact_resnum - chosen_seg->get_basis_pair().first.resnum();
		} else {
			new_contact_resnum = chosen_seg->get_basis_pair().first.resnum() + chimaeric_contact_resnum - chosen_seg->get_basis_pair().second.resnum();
		}
	}
	std::set< core::Size > new_contact_res;
	new_contact_res.insert( new_contact_resnum );
	data_storage::LigandSegmentOP replacement_parent( new data_storage::LigandSegment( *bad_parent, true, new_contact_res, ligand_->get_ligand_id() ) );
	replacement_parent->set_is_in_Assembly( true );
	//Change  replacement_ligseg's parent pointer
	if ( change_n_terminal_parent ) {
		replacement_ligseg->set_n_terminal_parent( data_storage::SmartSegmentOP( replacement_parent ) );
	} else {
		replacement_ligseg->set_c_terminal_parent( data_storage::SmartSegmentOP( replacement_parent ) );
	}

	//Now we have to replace the actual residue in the parent

	data_storage::SmartSewingResidueOP res_to_change = replacement_parent->get_residue( new_contact_resnum  );
	core::conformation::ResidueOP reference_res = make_aligned_residue( res_to_change, res_scores.at( 1 ).second.first.coord_res_name );
	for ( core::Size i = 1; i <= reference_res->nchi(); ++i ) {
		reference_res->set_chi( i, res_scores.at( 1 ).second.first.chis[ i ] );
	}

	//Now we'll need to change the SmartSewingResidue
	//Set the atom vector
	res_to_change->set_atom_vector( reference_res->atoms() );
	//Set the amino acid type
	res_to_change->set_amino_acid_type( reference_res->type().base_name() );
	res_to_change->set_full_type_name( reference_res->type().name() );
	//Set the chis??
	res_to_change->set_chis( reference_res->chi() );


	//Replace the parent in local_segments
	assembly->local_segments()[ bad_parent->get_segment_id() ] = data_storage::SmartSegmentOP( replacement_parent );
	if ( bad_parent->is_chimaeric() ) {
		//1) Figure out which of the parent's parents was N terminal

		bool first_parent_basis_is_n_terminal;
		if ( bad_parent->get_basis_pair().first.segment_id() == bad_parent->get_n_terminal_parent()->get_segment_id() ) {
			first_parent_basis_is_n_terminal = true;
		} else {
			first_parent_basis_is_n_terminal = false;
		}
		//2) Stick it on whichever side of the basis residue it falls on
		bool parent_has_bad_n_terminal_parent;
		if ( first_parent_basis_is_n_terminal ) {
			if ( new_contact_resnum <= bad_parent->get_basis_pair().first.resnum() ) {
				parent_has_bad_n_terminal_parent = true;
			} else {
				parent_has_bad_n_terminal_parent = false;
			}
		} else if ( new_contact_resnum <= bad_parent->get_basis_pair().second.resnum() ) {
			parent_has_bad_n_terminal_parent = true;
		} else {
			parent_has_bad_n_terminal_parent = false;
		}

		//Call this function on bad_parent
		replace_chimaera_ligand_parent( assembly, bad_parent, replacement_parent, parent_has_bad_n_terminal_parent, first_parent_basis_is_n_terminal, new_contact_resnum, res_scores);
		//replace_chimaera_ligand_parent( assembly, bad_parent, replacement_parent, parent_has_bad_n_terminal_parent, new_contact_resnum );
		//I think this should be safe here
		delete bad_parent.get();
		//*bad_parent = NULL;
	}
	//Otherwise we're done
}



void
LigandBindingResPlacer::revert_last_added_contact( data_storage::SmartAssemblyOP assembly ){
	//Our stored last_contact_added_ knows which segment/residue it is on
	//Things to do in this function:
	//1) Delete our stored contact from ligand_'s contacts
	if ( std::find( ligand_->get_current_contacts().begin(), ligand_->get_current_contacts().end(), last_contact_added_ ) != ligand_->get_current_contacts().end() ) {
		ligand_->get_nonconst_current_contacts().erase( std::find( ligand_->get_nonconst_current_contacts().begin(), ligand_->get_nonconst_current_contacts().end(), last_contact_added_ ) );
	}
	//2) Increment contacts_remaining_
	++contacts_to_make_;
	//3) Remove stored ligand contact and vital residue from the appropriate segment and its parent(s) (TODO)
	data_storage::SmartSegmentOP segment_with_contact = assembly->get_segment( last_contact_added_->segment_id );
	core::Size current_resnum = last_contact_added_->residue_number;

	//Delete current_resnum as a contact
	data_storage::LigandSegmentOP ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment >( segment_with_contact );
	if ( ligseg->get_nonconst_ligand_contact_indices().find( current_resnum ) != ligseg->get_nonconst_ligand_contact_indices().end() ) {
		ligseg->get_nonconst_ligand_contact_indices().erase( ligseg->get_nonconst_ligand_contact_indices().find( current_resnum ) );
	}
	//Delete current_resnum as a vital residue
	if ( segment_with_contact->get_vital_residues().find( current_resnum ) != segment_with_contact->get_vital_residues().end() ) {
		segment_with_contact->get_vital_residues().erase( segment_with_contact->get_vital_residues().find( current_resnum ) );
	}

	//Delete the ligand ID as a ligand residue unless there are other contacts
	if ( ligseg->get_ligand_contact_indices().size() == 0 ) {
		ligseg->get_nonconst_ligand_residues().clear();
	}

	//Maps the parent segment to the residue number in that parent that forms a contact
	std::map< data_storage::SmartSegmentOP, core::Size > parent_contacts;
	data_storage::SmartSegmentOP current_seg = segment_with_contact;
	while ( current_seg->is_chimaeric() ) {
		//Figure out which order the basis pair is in
		bool first_basis_is_n_terminal = ( current_seg->get_basis_pair().first.segment_id() == current_seg->get_n_terminal_parent()->get_segment_id()  );
		//Figure out if the contact is N-terminal or C-terminal of the basis residue
		bool change_n_terminal_parent;
		if ( first_basis_is_n_terminal && current_resnum <= current_seg->get_basis_pair().first.resnum() ) {
			change_n_terminal_parent = true;
		} else if ( first_basis_is_n_terminal ) {
			change_n_terminal_parent = false;
		} else if ( current_resnum <= current_seg->get_basis_pair().second.resnum() ) {
			change_n_terminal_parent = true;
		} else {
			change_n_terminal_parent = false;
		}
		//Set the new values for current_resnum and current_seg
		if ( !change_n_terminal_parent ) { //If we're changing the n terminal parent, current_resnum will stay the same
			if ( first_basis_is_n_terminal ) {
				current_resnum = current_seg->get_basis_pair().second.resnum() + current_resnum - current_seg->get_basis_pair().first.resnum();
			} else {
				current_resnum = current_seg->get_basis_pair().first.resnum() + current_resnum - current_seg->get_basis_pair().second.resnum();
			}
		}
		if ( change_n_terminal_parent ) {
			current_seg = current_seg->get_n_terminal_parent();
		} else {
			current_seg = current_seg->get_c_terminal_parent();
		}

		//Add them to the map
		parent_contacts[ current_seg ] = current_resnum;
		//Delete current_resnum as a contact
		data_storage::LigandSegmentOP parent_ligseg = std::dynamic_pointer_cast< data_storage::LigandSegment >( current_seg );
		if ( parent_ligseg->get_nonconst_ligand_contact_indices().find( current_resnum ) != parent_ligseg->get_nonconst_ligand_contact_indices().end() ) {
			parent_ligseg->get_nonconst_ligand_contact_indices().erase( parent_ligseg->get_nonconst_ligand_contact_indices().find( current_resnum ) );
		}
		//Delete current_resnum as a vital residue
		if ( current_seg->get_vital_residues().find( current_resnum ) != current_seg->get_vital_residues().end() ) {
			current_seg->get_vital_residues().erase( current_seg->get_vital_residues().find( current_resnum ) );
		}
		//Delete the ligand ID as a ligand residue unless there are other contacts
		if ( parent_ligseg->get_ligand_contact_indices().size() == 0 ) {
			parent_ligseg->get_nonconst_ligand_residues().clear();
		}
	}
	//4) Make the appropriate segment and parent(s) no longer vital
	//CAVEATS
	//If this segment was already vital before, then we want it to still be vital
	//If it has any vital residues and/or any ligand contacts, then we want it to still be vital
	//If either adjacent segment (loop) is vital, it needs to remain vital
	if ( segment_with_contact->get_vital_residues().size() == 0 ) {
		//Maybe this isn't the most intuitive way to write this, but I don't feel like changing it now
		if ( !( ( segment_with_contact->get_n_terminal_neighbor() && segment_with_contact->get_n_terminal_neighbor()->is_vital() ) || ( segment_with_contact->get_c_terminal_neighbor() && segment_with_contact->get_c_terminal_neighbor()->is_vital() ) ) ) {
			segment_with_contact->set_is_vital( false );
		}
	}
	for ( std::pair< data_storage::SmartSegmentOP, core::Size > parent: parent_contacts ) {
		if ( parent.first->get_vital_residues().size() != 0 ) {
			continue;
		}
		if ( ( parent.first->get_n_terminal_neighbor() && parent.first->get_n_terminal_neighbor()->is_vital() ) || ( parent.first->get_c_terminal_neighbor() && parent.first->get_c_terminal_neighbor()->is_vital() ) ) {
			continue;
		}
		//Check all of the same caveats as above
		parent.first->set_is_vital( false );
	}
}



//Getters
data_storage::LigandResidueOP
LigandBindingResPlacer::get_ligand() const{
	return ligand_;
}
/*
std::map< core::Size, IdealContact >
LigandBindingResPlacer::get_ideal_contacts() const{
return ideal_contacts_;
}
*/
std::map< std::pair< core::Size, core::Size >, std::set< LigandCoordInfo> > const &
LigandBindingResPlacer::get_possible_coordinations() const{
	return possible_coordinations_;
}

utility::pointer::shared_ptr< std::map< char, LigandHashMap > const >
LigandBindingResPlacer::get_ligand_local_coords() const{
	return ligand_local_coords_;
}

void
LigandBindingResPlacer::set_ligand_local_coords( utility::pointer::shared_ptr< std::map< char, LigandHashMap > > local ){
	ligand_local_coords_ = local;
}

core::Size
LigandBindingResPlacer::contacts_remaining() const{
	return contacts_to_make_;
}

core::Real
LigandBindingResPlacer::get_geometry_score_weight() const{
	return geometry_score_weight_;
}

core::Real
LigandBindingResPlacer::get_best_geometry_score() const{
	return best_geometry_score_;
}
//Setters

void
LigandBindingResPlacer::set_ligand( data_storage::LigandResidueOP lig ){
	ligand_ = lig;
	//Calculate contacts_to_make_
	core::Size desired_contacts = 0;
	for ( std::pair< core::Size, hashing::IdealContact > ideal: ligand_->get_ideal_contacts() ) {
		desired_contacts += ideal.second.max_contacts;
	}
	core::Size current_contacts = ligand_->get_current_contacts().size();
	contacts_to_make_ = desired_contacts - current_contacts;
}
void
LigandBindingResPlacer::set_geometry_score_weight( core::Real gsw ){
	geometry_score_weight_ = gsw;
}

void
LigandBindingResPlacer::set_best_geometry_score( core::Real best ){
	best_geometry_score_ = best;
}
//XML-related stuff

void
LigandBindingResPlacer::parse_ideal_contacts( utility::tag::TagCOP ligands_tag, utility::vector1< data_storage::LigandDescription > & ligands ) {
	core::Size counter = 1;
	//Iterate over Ligand tags
	for ( utility::tag::TagCOP ligand_tag: ligands_tag->getTags() ) {
		data_storage::LigandDescription & current_ligand = ligands[ counter ];


		for ( utility::tag::TagCOP coord_tag: ligand_tag->getTags() ) {
			if ( coord_tag->getName() != "Coordination" ) {
				continue;
			}
			std::string file_names_string = coord_tag->getOption< std::string >( "coordination_files" );
			current_ligand.ligand_coord_files = utility::string_split_simple( file_names_string, ',' );
			current_ligand.geometry_score_threshold = coord_tag->getOption< core::Real >( "geometry_score_threshold", 1 );

			for ( utility::tag::TagCOP ideal_contacts_tag: coord_tag->getTags() ) {
				//Each IdealContacts subtag has attributes for ligand atom name, ideal contact distance, ideal angle between contacts, and max # of contacts
				IdealContact current_atom;
				current_atom.preferred_distance = ideal_contacts_tag->getOption< core::Real >( "distance" );
				//Angle will be provided in degrees; we should convert it to radians here since that's how it's processed
				current_atom.preferred_contact_angle = ideal_contacts_tag->getOption< core::Real >( "angle", 109.5 ) * numeric::constants::r::degrees_to_radians;
				current_atom.max_contacts = ideal_contacts_tag->getOption< core::Size >( "max_coordinating_atoms", 4 );
				//Ligand atom is provided as a name but must be converted to a number
				//We won't be able to actually make the conversion until later on when we have the pose and construct the LigandResidue objects
				std::string ligand_atom_name = ideal_contacts_tag->getOption< std::string >( "ligand_atom_name" );
				current_atom.ligand_atom = 1;
				current_atom.preferred_dihedral_1 = ideal_contacts_tag->getOption< core::Real >( "dihedral_1", 30 ) * numeric::constants::r::degrees_to_radians;
				current_atom.preferred_dihedral_2 = ideal_contacts_tag->getOption< core::Real >( "dihedral_2", 30 ) * numeric::constants::r::degrees_to_radians;
				current_ligand.ideal_contacts_str[ ligand_atom_name ] = current_atom;
			}
		}
		++counter;
	}
}


void
LigandBindingResPlacer::define_ideal_contacts_subelement( utility::tag::XMLSchemaSimpleSubelementList & subs, utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;

	//Define the attributes for the IdealContacts subelement
	AttributeList ideal_contact_attributes;
	ideal_contact_attributes
		+ XMLSchemaAttribute::required_attribute( "distance", xsct_real, "Ideal distance between ligand and contact atom" )
		+ XMLSchemaAttribute::attribute_w_default( "angle", xsct_real, "Ideal angle between this atom's contacts",  "109.5" )
		+ XMLSchemaAttribute::attribute_w_default( "dihedral_1", xsct_real, "Ideal dihedral angle: contact_base - contact - ligand_atom - other_contact", "30" )
		+ XMLSchemaAttribute::attribute_w_default( "dihedral_2", xsct_real, "Ideal dihedral angle: contact - ligand_atom - other_contact - other_base", "30" )
		+ XMLSchemaAttribute::required_attribute( "max_coordinating_atoms", xsct_non_negative_integer, "Maximum number of contacts that this atom can form. Note that IdealContacts tags do not need to be defined for atoms with no contacts." )
		+ XMLSchemaAttribute::required_attribute( "ligand_atom_name", xs_string, "Rosetta name of the ligand atom to which this tag applies" );

	//Define IdealContacts subelement list
	XMLSchemaSimpleSubelementList coord_subs;
	coord_subs
		.add_simple_subelement( "IdealContacts", ideal_contact_attributes, "Describes the ideal coordination environment for an atom within the ligand" );

	AttributeList coord_attributes;
	coord_attributes
		+ XMLSchemaAttribute( "coordination_files", xs_string, "Comma-separated list of coordination file names for this ligand" )
		+ XMLSchemaAttribute::attribute_w_default( "geometry_score_threshold", xsct_real, "Maximum score geometry score to allow when forming a contact", "1" );

	//Define the complex type for a Coordination subelement
	XMLSchemaComplexTypeGenerator coord_ct_gen;
	coord_ct_gen.complex_type_naming_func( & ligand_subtag_ct_namer )
		.element_name( "Coordination" )
		.add_attributes( coord_attributes )
		.set_subelements_repeatable( coord_subs )
		.description( "Contains subtags defining ideal coordination environments for atoms in the ligand" )
		.write_complex_type_to_schema( xsd );

	//Add the Coordination subelement to subs as an already defined subelement
	subs.add_already_defined_subelement("Coordination", & ligand_subtag_ct_namer );
}

std::string
LigandBindingResPlacer::ligand_subtag_ct_namer( std::string tag ){
	return "ligand_" + tag + "_complex_type";
}
} //namespace hashing
} //namespace sewing
} //namespace protocols

