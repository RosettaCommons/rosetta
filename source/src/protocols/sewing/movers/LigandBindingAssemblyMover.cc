// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/LigandBindingAssemblyMover.hh
/// @brief an AssemblyMover for adding to existing poses
/// @author frankdt (frankdt@email.unc.edu)


// Unit headers
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/movers/LigandBindingAssemblyMover.hh>
#include <protocols/sewing/movers/LigandBindingAssemblyMoverCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorerCreator.hh>
#include <protocols/sewing/scoring/AssemblyScorerFactory.hh>
#include <protocols/sewing/requirements/AssemblyRequirementCreator.hh>
#include <protocols/sewing/requirements/AssemblyRequirementFactory.hh>
#include <protocols/sewing/requirements/LigandAssemblyRequirement.hh>
#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/hashing/EdgeMapGenerator.hh>
#include <protocols/sewing/hashing/AlignmentFileGeneratorMover.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/hashing/LigandBindingResPlacer.hh>
#include <protocols/sewing/data_storage/LigandResidue.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/moves/mover_schemas.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
//#include <shared_ptr>

static basic::Tracer TR( "protocols.sewing.movers.LigandBindingAssemblyMover" );

namespace protocols {
namespace sewing {
namespace movers {


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
LigandBindingAssemblyMover::LigandBindingAssemblyMover():
	AppendAssemblyMover()
{
	//By default, no ligand conformers, alignment atoms, or ideal contacts
}

/// @brief Copy constructor (not needed unless you need deep copies)
LigandBindingAssemblyMover::LigandBindingAssemblyMover( LigandBindingAssemblyMover const & src ):
	AppendAssemblyMover( src )
{
	set_distance_cutoff( src.get_distance_cutoff() );
	set_segment_distance_cutoff( src.get_segment_distance_cutoff() );
	//Do a deep copy of this map? Should be unnecessary since it's const & on the heap
	set_ligand_requirements( src.get_ligand_requirements() );
	set_non_ligand_requirements( src.get_non_ligand_requirements() );
	set_build_site_only( src.get_build_site_only() );
	expanded_ligands_ = src.expanded_ligands_;
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
LigandBindingAssemblyMover::apply( core::pose::Pose & pose ){


	core::Size starttime = time(NULL);
	set_last_move_status(protocols::moves::MS_SUCCESS);
	data_storage::SmartAssemblyOP assembly = set_up_assembly( pose );
	LigandBindingAssemblyMover::generate_assembly( assembly, pose);
	//If we failed, set to FAIL_DO_NOT_RETRY
	if ( assembly == nullptr ) {
		TR << "Failed to generate an Assembly" << std::endl;
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		return;
	} else if ( get_last_move_status() == protocols::moves::FAIL_DO_NOT_RETRY ) {
		TR << "Failed to find ligand contacts" << std::endl;
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		return;
	}

	core::Size endttime = time(NULL);

	TR << "Assembly successfully generated in " << endttime - starttime << " seconds" << std::endl;

	AssemblyMover::print_statistics( assembly );

	//Now convert the assembly to a pose
	pose = assembly->to_pose("fa_standard");

	add_motif_scorers_to_score_file( pose, assembly );

	//At some point in here we'll need to output the side chain coordinates
	//Again, no idea what this function will look like
	assembly->dump_side_chains();


}


data_storage::SmartAssemblyOP
LigandBindingAssemblyMover::set_up_assembly(core::pose::Pose & pose){
	data_storage::SmartAssemblyOP assembly;
	//Annoyingly, for the references to work I have to copy these vectors here
	utility::vector1< data_storage::LigandDescription > ligands = get_ligands();
	std::map< core::Size, data_storage::LigandResidueCOP > partner_ligands;
	//utility::vector1< data_storage::LigandDescription > expanded_ligands;
	std::string required_resnums = get_required_resnums();
	bool ssfd = get_segments_from_dssp();
	core::select::residue_selector::ResidueSelectorCOP required_selector = get_required_selector();
	if ( get_hashed() ) {
		data_storage::HashedSmartAssemblyOP hashed_assembly( new data_storage::HashedSmartAssembly( get_segment_vector() ) );
		hashed_assembly->set_basis_map_generator( get_basis_map_generator() );
		assembly = hashed_assembly;

		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector( pose, get_partner_pdb(), get_basis_map_generator(), ligands, partner_ligands, expanded_ligands_, required_resnums, required_selector, ssfd );
		this->set_ligands( ligands );
		assembly->pdb_segments( get_basis_map_generator()->pdb_segments() ); //This also adds them to the local segment vector
		assembly->set_partner( get_partner_pdb() );
		assembly->set_partner_ligands( partner_ligands );

		for ( data_storage::LigandDescription & ligdes: expanded_ligands_ ) {
			assembly->load_initial_conformers( ligdes );
		}

	} else {
		assembly = data_storage::SmartAssemblyOP( new data_storage::SmartAssembly( get_segment_vector(), get_window_width() ) );
		std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs;
		pdbsegs.clear();
		hashing::AlignmentFileGeneratorMover::add_pose_segments_to_segment_vector( pose, get_partner_pdb(), assembly->get_segment_vector(), pdbsegs, get_pose_segment_starts_string(), get_pose_segment_ends_string(), get_pose_segment_dssp(), ligands, partner_ligands, expanded_ligands_, required_resnums, required_selector, get_strict_dssp_changes() );
		this->set_ligands( ligands );
		TR.Debug << "PDB segments added" << std::endl;
		assembly->pdb_segments( pdbsegs );
		assembly->set_partner( get_partner_pdb() );
		assembly->set_partner_ligands( partner_ligands );

		for ( data_storage::LigandDescription & ligdes: expanded_ligands_ ) {
			assembly->load_initial_conformers( ligdes );
		}
		//Our starting segment will be one of the segments from our input pose (either the first or the last if it's not hashed
	}
	//We want to wait until now to do this so we won't lose our required residues
	data_storage::SmartSegmentOP starting_segment = AppendAssemblyMover::find_starting_segment( assembly );
	assembly->set_starting_segment( starting_segment, assembly->get_start_node_vital_segments() );

	return assembly;
}


core::Real
LigandBindingAssemblyMover::get_min_distance( data_storage::SmartSegmentOP seg, data_storage::LigandResidueOP ligand){
	//Should we compare to the "center of mass" of the ligand? Or to all ligand atoms?
	utility::vector1< numeric::xyzVector< core::Real > > ligand_coordinates;
	core::Real min_distance = 100;
	for ( core::conformation::Atom atom: ligand->get_atom_vector() ) {
		ligand_coordinates.push_back( atom.xyz() );
	}
	numeric::xyzVector< core::Real > ligand_center = numeric::center_of_mass( ligand_coordinates );
	//Now iterate over the residues in the segment
	for ( data_storage::SmartSewingResidueOP res: seg->get_residue_vector() ) {
		numeric::xyzVector< core::Real > calpha = res->get_atom( 3 ).xyz();
		core::Real distance = ligand_center.distance( calpha );
		if ( distance < min_distance ) {
			min_distance = distance;
		}
	}
	return min_distance;
}

void
LigandBindingAssemblyMover::generate_assembly(data_storage::SmartAssemblyOP assembly, core::pose::Pose & pose ){

	//PHASE 1: Add ligand contacts
	//This will involve some amount of add/delete/switch moves (with LigandBindingResPlacer after each) plus conformer sampling with some probability

	//It will run until all ligand contacts are formed OR until it reaches a maximum # cycles (and then fails)
	//This stage will also have its own requirement set with a LengthRequirement, ClashRequirement, LigandClashRequirement (combine with ClashRequirement??),DsspSpecificLengthRequirement, LigandContactRequirement

	//Have a map of Ligand ID to LigandBindingResPlacer
	//Add functions to SmartAssembly to get ligand from ID (like we have with getting segments from SegmentID

	//Step 0: Initialize variables
	core::Real score = 1000;
	core::Real old_score = 1000;
	binder_finders_.clear();
	core::Real add_prob = 0.5;
	bool revertable_move = true;
	core::Size adds_to_n_term = 0;
	core::Size adds_to_c_term = 0;
	bool add_to_n_terminus = true;
	bool contacts_done = false;
	data_storage::SmartAssemblyOP best_so_far = assembly;
	core::Real best_score_so_far = 1000;
	std::map< core::Size, core::Size > remaining_contacts_to_make;
	//TODO--weight_it is not guaranteed to be correct

	for ( data_storage::LigandDescription & ligand_des: this->get_ligands() ) {
		if ( ligand_des.ideal_contacts_num.size() == 0 ) {
			continue;
		}
		binder_finders_[ ligand_des.ligand_id ] = hashing::LigandBindingResPlacerOP( new hashing::LigandBindingResPlacer( ligand_des.ligand_coords, assembly->get_local_ligands().at( ligand_des.ligand_id ) ) );
		binder_finders_[ ligand_des.ligand_id ]->set_geometry_score_weight( ligand_des.geometry_score_threshold );
	}

	//Now, in a loop over the maximum number of binder finding cycles:
	core::Size cycle = 1;
	TR.Debug << "Will attempt to find ligand contacts for " << binding_cycles_ << " cycles" << std::endl;
	while ( cycle <= binding_cycles_ ) {
		//for( core::Size cycle = 1; cycle <= binding_cycles_; ++cycle ){
		old_score = score;
		TR.Debug << "Binding Cycle: " << cycle << std::endl;
		AssemblyMover::print_statistics( assembly );

		//Step 1: Make move (add/delete/switch)
		//Step 1A: Decide add probability and which terminus we can/should add to
		//Switches and deletes could occur anywhere
		if ( !assembly->can_delete() ) {
			TR.Debug << "N terminus " << assembly->get_n_terminal_segment()->get_segment_id() << "returns " << assembly->get_n_terminal_segment()->is_vital() << std::endl;
			TR.Debug << "C terminus " << assembly->get_c_terminal_segment()->get_segment_id() << "returns " << assembly->get_c_terminal_segment()->is_vital() << std::endl;
			TR.Debug << "All termini vital. Adding." << std::endl;
			add_prob = 1.0;
		} else if ( adds_to_n_term >= segment_distance_cutoff_ && adds_to_c_term >= segment_distance_cutoff_ ) {
			TR.Debug << "Assembly at max size for finding ligand contacts. Cannot add." << std::endl;
			add_prob = 0.0;
		} else if ( adds_to_n_term >= segment_distance_cutoff_ ) {
			add_to_n_terminus = false;
			add_prob = AssemblyMover::get_add_probability();
		} else if ( adds_to_c_term >= segment_distance_cutoff_ ) {
			add_to_n_terminus = true;
			add_prob = AssemblyMover::get_add_probability();
		} else {
			//Randomly decide which terminus we should add to
			add_to_n_terminus = assembly->get_modifiable_terminus( 'A' );
			add_prob = AssemblyMover::get_add_probability();
		}

		TR.Debug << "Add probability: " << add_prob << std::endl;
		//Step 1b: Decide which type of move to make
		core::Real random_action;
		random_action = numeric::random::rg().uniform();
		TR.Debug << "random action " << random_action << std::endl;
		//Step 1c: Make our move
		if ( random_action < add_prob ) {
			TR.Debug << "Adding a segment" << std::endl;
			if ( !(assembly->add_segment( add_to_n_terminus ) ) ) { // add happens here
				TR.Debug << "Move failed!" << std::endl;
				revertable_move = false;
				//cycle--;
			} else {
				++cycle;
				revertable_move = true;
				if ( add_to_n_terminus ) {
					++adds_to_n_term;
					TR.Debug << adds_to_n_term << " adds to n term" << std::endl;
				} else { //end added to n term
					++adds_to_c_term;
					TR.Debug << adds_to_c_term << " adds to c term" << std::endl;
				} //end added to c term
			}//end move succeeded
		} else if ( random_action < ( add_prob + AssemblyMover::get_delete_probability() ) ) { //end added
			TR.Debug << "Deleting a segment" << std::endl;
			bool n_terminal_delete =  assembly->get_modifiable_terminus('D');
			if ( !( assembly->delete_segment( n_terminal_delete ) ) ) { //delete happens here
				TR.Debug << "Move failed!" << std::endl;
				revertable_move = false;
				//cycle--;
			} else { //end if failed
				++cycle;
				revertable_move = true;
				if ( n_terminal_delete ) {
					--adds_to_n_term;
					TR.Debug << adds_to_n_term << " adds to n term" << std::endl;
				} else { //end if N term
					--adds_to_c_term;
					TR.Debug << adds_to_c_term << " adds to c term" << std::endl;
				} //end if c term
			} //end succeeded
		} else { //end delete
			TR.Debug << "Switching a segment" << std::endl;
			if ( !( assembly->switch_segment( assembly->get_modifiable_terminus('S') ) ) ) {
				TR.Debug << "Move failed!" << std::endl;
				revertable_move = false;
				//cycle--;
			} else { //end if failed
				++cycle;
				revertable_move = true;
			} //end if succeeded
		}
		//If the move (whatever it was) failed, try again
		if ( !revertable_move ) {
			TR.Debug << "This cycle didn't count!" << std::endl;
			continue;
		}
		//Step 2: Apply non-ligand requirements
		char last_move = assembly->get_last_change();
		bool revert = false;
		for ( requirements::AssemblyRequirementOP requirement: non_ligand_requirements_ ) {
			if ( !requirement->test( assembly ).first ) {
				TR.Debug << "Requirement " << requirement->get_name() << " failed!" << std::endl;
				revert = true;
				break;
			}
		}
		//Step 3: score assembly with non-ligand scorers (get boolean, will revert if false)
		if ( !revert ) {
			score = AssemblyMover::score( assembly );
			TR.Debug << "Old assembly score: " << old_score << std::endl;
			TR.Debug << "New assembly score: " << score << std::endl;
			//40 and -40 are just there as safeguards
			if ( score > old_score && numeric::random::rg().uniform() > std::exp( ( old_score - score )/ AssemblyMover::get_start_temperature() ) ) {
				revert = true;
			}
		}
		//Step 4: If the last move was an add or switch that will not be reverted:
		bool added_contact = false;
		utility::vector1< hashing::LigandBindingResPlacerOP > binder_finders_that_added_contacts;
		if ( !revert && ( last_move == 'A' || last_move == 'S'  ) ) {
			TR.Debug << "Unreverted add/switch. Checking for ligand contacts." << std::endl;
			//Step a: Check all the segments that we added
			//Check all non-loop segments
			data_storage::SmartSegmentOP added_segment;
			if ( assembly->get_last_change_was_n_terminal() ) {
				added_segment = assembly->get_n_terminal_segment();
			} else {
				added_segment = assembly->get_c_terminal_segment();
			}
			//Go until you reach a chimaera, but check chimaeric segments too
			bool done = false;
			while ( !done ) {
				TR.Debug << "Checking segment " << added_segment->get_segment_id() << std::endl;
				if ( added_segment->is_chimaeric() ) {
					done = true;
				}
				if ( added_segment->get_dssp_code() != 'L' ) {
					//Step i: Check minimum distance from Calpha in new segment (if applicable) to ligand(s) (false if not)
					core::Real min_distance = 100;
					for ( std::pair< core::Size, hashing::LigandBindingResPlacerOP > bf: binder_finders_ ) {
						TR.Debug << "Checking ligand " << bf.first << std::endl;
						core::Real distance = get_min_distance( added_segment, assembly->get_local_ligands()[ bf.first ] );
						if ( distance < min_distance ) {
							//min_ligand_id = ligand.first;
							min_distance = distance;
						}
						TR.Debug << "Minimum distance: " << min_distance << std::endl;
						TR.Debug << "Distance cutoff: " << distance_cutoff_ << std::endl;
						if ( min_distance > distance_cutoff_ ) {
							done = true;
							revert = true;
							break;
						}
						//Step ii: Try to build a contact from the new non-loop segments to ligand(s) in range
						//if( min_distance <= distance_cutoff_ ){
						TR.Debug << "Should we check for new contacts for ligand " << bf.first << "?" << std::endl;
						core::Size segID = added_segment->get_segment_id();
						std::pair< data_storage::SmartSegmentOP, bool > cbm_results = bf.second->choose_best_metal_coordinator( assembly );
						bool new_contact = cbm_results.second;
						if ( new_contact ) {
							if ( segID == cbm_results.first->get_segment_id() ) { //added_segment has been replaced
								added_segment = cbm_results.first;
							}
							binder_finders_that_added_contacts.push_back( bf.second );
							TR.Debug << "Added a new contact to ligand " << bf.first << std::endl;
						}
						added_contact = ( added_contact || new_contact );
					}



				}
				if ( assembly->get_last_change_was_n_terminal() ) {
					added_segment = added_segment->get_c_terminal_neighbor();
				} else {
					added_segment = added_segment->get_n_terminal_neighbor();
				}
			}
			//Step b: Decide whether to revert based on whether we added a contact
			//Might not need to do anything since we can always switch/delete either way
			if ( !added_contact ) {
				TR.Debug << "Did not add ligand contact" << std::endl;
				if ( score < best_score_so_far ) {
					best_so_far = assembly;
					best_score_so_far = score;
				}
			}
			if ( added_contact ) {
				best_so_far = assembly;
				best_score_so_far = score;
			}
		} else if ( !revert ) {
			if ( score < best_score_so_far ) {
				best_so_far = assembly;
				best_score_so_far = score;
			}
		}
		//Step 5: Apply ligand requirements
		for ( requirements::AssemblyRequirementOP requirement: ligand_requirements_ ) {
			if ( !requirement->test( assembly ).first ) {
				TR.Debug << "Requirement " << requirement->get_name() << " failed!" << std::endl;
				revert = true;
				break;
			}
		}
		//Step 5b: Revert addition of contacts if necessary
		if ( added_contact && revert ) {
			//Which binder finder was it from?
			for ( hashing::LigandBindingResPlacerOP bf: binder_finders_that_added_contacts ) {
				bf->revert_last_added_contact( assembly );
			}
		}
		//Step 6: Revert if necessary:
		if ( revert ) {
			TR.Debug << "Reverting change." << std::endl;
			print_statistics( assembly );
			if ( last_move == 'A' ) {
				if ( assembly->get_last_change_was_n_terminal() ) {
					--adds_to_n_term;
					TR.Debug << adds_to_n_term << " adds to n term" << std::endl;
				} else {
					--adds_to_c_term;
					TR.Debug << adds_to_c_term << " adds to c term" << std::endl;
				}
			} else if ( last_move == 'D' ) {
				if ( assembly->get_last_change_was_n_terminal() ) {
					++adds_to_n_term;
					TR.Debug << adds_to_n_term << " adds to n term" << std::endl;
				} else {
					++adds_to_c_term;
					TR.Debug << adds_to_c_term << " adds to c term" << std::endl;
				}
			}
			//If it was a switch, we don't need to change any numbering
			assembly->revert();
			score = old_score;
		}
		contacts_done = true;
		for ( std::pair< core::Size, hashing::LigandBindingResPlacerOP > finder: binder_finders_ ) {
			if ( finder.second->contacts_remaining() != 0 ) {
				contacts_done = false;
			}
		}
		if ( contacts_done ) {
			//We can move onto the next phase of assembly
			break;
		}
	}
	//assembly = best_so_far;
	//PHASE 1.5: Figure out if the ligand contacts were satisfied properly. If they were not, we should not go on to the next step and should instead FAIL_RETRY
	//Figure out if we're done
	if ( !contacts_done ) {
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		assembly = nullptr;
		return;
	}
	//PHASE 2: More or less the same as AssemblyMover's generate_assembly (with conformer switch moves)
	if ( !build_site_only_ ) {
		AssemblyMover::generate_assembly( assembly , pose);
	}
}


/// @brief Show the contents of the Mover
std::string
LigandBindingAssemblyMover::class_name(){
	return "LigandBindingAssemblyMover";
}

void
LigandBindingAssemblyMover::show( std::ostream & output ) const{
	protocols::moves::Mover::show( output );
}

/// @brief Get the name of the Mover
std::string
LigandBindingAssemblyMover::get_name() const{
	return LigandBindingAssemblyMover::class_name();
}

//Getters
void
LigandBindingAssemblyMover::add_motif_scorers_to_score_file( core::pose::Pose & pose, data_storage::SmartAssemblyOP & assembly ){
	AssemblyMover::add_motif_scorers_to_score_file( pose, assembly );
	for (  std::pair< core::Size, hashing::LigandBindingResPlacerOP > bf:  binder_finders_ ) {
		std::string scorer_name = "ligand_" + bf.second->get_ligand()->get_amino_acid_type() + "_" + utility::to_string( bf.first );
		core::pose::setPoseExtraScore( pose, scorer_name, bf.second->get_best_geometry_score() );
	}
}

core::Real
LigandBindingAssemblyMover::get_distance_cutoff() const{
	return distance_cutoff_;
}

core::Size
LigandBindingAssemblyMover::get_segment_distance_cutoff() const{
	return segment_distance_cutoff_;
}

core::Size
LigandBindingAssemblyMover::get_binding_cycles() const{
	return binding_cycles_;
}
utility::vector1< requirements::AssemblyRequirementOP >
LigandBindingAssemblyMover::get_ligand_requirements() const{
	return ligand_requirements_;
}

utility::vector1< requirements::AssemblyRequirementOP >
LigandBindingAssemblyMover::get_non_ligand_requirements() const{
	return non_ligand_requirements_;
}

bool
LigandBindingAssemblyMover::get_build_site_only() const{
	return build_site_only_;
}



//Setters

void
LigandBindingAssemblyMover::set_distance_cutoff( core::Real cutoff ){
	distance_cutoff_ = cutoff;
}

void
LigandBindingAssemblyMover::set_segment_distance_cutoff( core::Size cutoff ){
	segment_distance_cutoff_ = cutoff;
}

void
LigandBindingAssemblyMover::set_binding_cycles( core::Size cycles ){
	binding_cycles_ = cycles;
}

void
LigandBindingAssemblyMover::set_ligand_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements ){
	ligand_requirements_ = requirements;
}
void
LigandBindingAssemblyMover::set_non_ligand_requirements( utility::vector1 <requirements::AssemblyRequirementOP > requirements ){
	non_ligand_requirements_ = requirements;
}

void
LigandBindingAssemblyMover::set_build_site_only( bool build_only ){
	build_site_only_ = build_only;
}

///////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
LigandBindingAssemblyMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	AppendAssemblyMover::parse_my_tag( tag, data, filters, movers, pose );

	//Parse the extra attributes
	set_distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 10.0 ) );
	set_segment_distance_cutoff( tag->getOption< core::Size >( "segment_distance_cutoff", 1 ) );
	binding_cycles_ = tag->getOption< core::Size >( "binding_cycles", 1000 );
	build_site_only_ = tag->getOption< bool >( "build_site_only", false );
	TR << "Number of binding cycles: " << binding_cycles_ << std::endl;
	//Call parse_ideal_contacts on the tag
	if ( tag->hasTag( "Ligands" ) ) {
		utility::vector1< data_storage::LigandDescription > ligands = get_ligands();
		//ligands = this->get_ligands();
		utility::tag::TagCOP ligands_tag = tag->getTag( "Ligands" );
		hashing::LigandBindingResPlacer::parse_ideal_contacts( ligands_tag, ligands );
		hashing::LigandBindingResPlacer::parse_ligand_files( ligands );
		this->set_ligands( ligands );
	}
	//get our ligand requirements and non-ligand requirements organized
	//Requirements were set in the parent's parse_my_tag
	for ( requirements::AssemblyRequirementOP req: AssemblyMover::get_requirements() ) {
		if ( std::dynamic_pointer_cast< requirements::LigandAssemblyRequirementOP >( req ) != nullptr ) {
			ligand_requirements_.push_back( req );
		} else {
			non_ligand_requirements_.push_back( req );
		}
	}
}

//LigandBindingAssemblyMover & operator=( LigandBindingAssemblyMover const & src );

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LigandBindingAssemblyMover::fresh_instance() const{
	return protocols::moves::MoverOP( new LigandBindingAssemblyMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LigandBindingAssemblyMover::clone() const{
	return protocols::moves::MoverOP( new LigandBindingAssemblyMover( *this ) );
}

void
LigandBindingAssemblyMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	scoring::AssemblyScorerFactory::get_instance()->define_assembly_scorer_subtag( xsd );
	requirements::AssemblyRequirementFactory::get_instance()->define_assembly_requirement_subtag( xsd );
	//The ligands subtag will be somewhat different--we should probably do it on a ligand-by-ligand basis (???)
	//Solution: Provide a boolean to the function telling it whether or not to include the ideal_contacts stuff for each ligand
	XMLSchemaSimpleSubelementList mover_subelements;
	mover_subelements
		.add_already_defined_subelement( "AssemblyScorers", & AssemblyMover::assembly_mover_subtag_ct_namer )
		.add_already_defined_subelement( "AssemblyRequirements", & AssemblyMover::assembly_mover_subtag_ct_namer );
	hashing::AlignmentFileGeneratorMover::append_ligands_subelement( mover_subelements, xsd, true );

	AttributeList attributes;
	AppendAssemblyMover::attributes_for_append_assembly_mover( attributes );
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "distance_cutoff", xsct_real, "Cutoff for distance in Angstroms between segment and ligand to check for contacts", "10.0" )
		+ XMLSchemaAttribute::attribute_w_default( "segment_distance_cutoff", xsct_non_negative_integer, "Maximum number of segments away from the starting structure to look for contacts", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "binding_cycles", xsct_non_negative_integer, "How many cycles to check for ligand contacts before giving up", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "build_site_only", xsct_rosetta_bool, "Should we stop after finding all of the desired contacts?", "false" );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( class_name() )
		.add_attributes( attributes )
		.description( "AssemblyMover designed to add contacts to specified ligand atoms and then build an assembly around the ligand." )
		.add_optional_name_attribute( "Name to identify this mover" )
		.set_subelements_single_appearance_optional( mover_subelements )
		.write_complex_type_to_schema( xsd ); //We won't have regular score functions/task operations

}



////////////////////////////////////////////////////////////////////////////////
/// Creator ///
///////////////

/////////////// Creator ///////////////

protocols::moves::MoverOP
LigandBindingAssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new LigandBindingAssemblyMover );
}

//Why do we have three functions in the creator that do exactly the same thing??

std::string
LigandBindingAssemblyMoverCreator::keyname() const
{
	return LigandBindingAssemblyMover::class_name();
}
std::string
LigandBindingAssemblyMoverCreator::class_name()
{
	return LigandBindingAssemblyMover::class_name();
}
std::string
LigandBindingAssemblyMoverCreator::mover_name()
{
	return LigandBindingAssemblyMover::class_name();
}

void
LigandBindingAssemblyMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LigandBindingAssemblyMover::provide_xml_schema( xsd );
}



std::ostream &
operator<<( std::ostream & os, LigandBindingAssemblyMover const & mover ){
	mover.show( os );
	return os;
}

} //protocols
} //sewing
} //movers

