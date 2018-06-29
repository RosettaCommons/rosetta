// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/guffysl/metal_site_finder.cc
/// @brief Application to find metal binding sites in PDBs.
/// @author Sharon Guffy


// Headers
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/conformation/AtomGraph.hh>
#include <core/conformation/AtomGraphData.hh>
#include <core/conformation/Residue.hh>
//#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>


#include <string>
//tracers
using basic::Error;
using basic::Warning;
//using basic::T;
static basic::Tracer TR("apps.pilot.guffysl.metal_site_finder");

//NOTE: Structures must be cleaned to remove all non-CA ligands before running this protocol.



//Define local options
namespace local {
basic::options::StringOptionKey const metal_type( "metal_type" );
basic::options::RealOptionKey const metal_ligand_distance_cutoff( "metal_ligand_distance_cutoff" );
basic::options::IntegerOptionKey const min_total_coord_atoms ( "min_total_coord_atoms" );
basic::options::IntegerOptionKey const min_local_coord_atoms ( "min_local_coord_atoms" );
basic::options::IntegerOptionKey const min_coord_res ( "min_coord_res" );
basic::options::IntegerOptionKey const max_res_span ( "max_res_span" );
basic::options::BooleanOptionKey const require_bordering_ss ( "require_bordering_ss" );
basic::options::IntegerOptionKey const max_loop_size ( "max_loop_size" );
basic::options::BooleanOptionKey const helices_only ( "helices_only" );
basic::options::BooleanOptionKey const strict_ss_changes( "strict_ss_changes" );
basic::options::BooleanOptionKey const append_nonlocal_residues( "append_nonlocal_residues" );
basic::options::StringOptionKey const allowed_elements( "allowed_elements" );
}

//Header file information
class MetalSiteFinderMover: public protocols::moves::Mover{
public:
	//Default constructor
	MetalSiteFinderMover();
	//Copy constructor
	MetalSiteFinderMover (MetalSiteFinderMover const & other);
	//Destructor
	virtual ~MetalSiteFinderMover();
	//Set private data
	//Get private data
	//Virtual methods
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual void apply(core::pose::Pose & pose);
	virtual std::string get_name() const;
	//For RosettaScripts
	//Any other methods
private:
	//private data
};

//Owning pointers
typedef utility::pointer::shared_ptr< MetalSiteFinderMover > MetalSiteFinderMoverOP;
typedef utility::pointer::shared_ptr< MetalSiteFinderMover const > MetalSiteFinderMoverCOP;


//Define methods
MetalSiteFinderMover::MetalSiteFinderMover(){
}

MetalSiteFinderMover::MetalSiteFinderMover(MetalSiteFinderMover const & other):
	protocols::moves::Mover(other)
{
}
MetalSiteFinderMover::~MetalSiteFinderMover(){}

//Virtual methods
protocols::moves::MoverOP MetalSiteFinderMover::clone() const{
	return protocols::moves::MoverOP( new MetalSiteFinderMover(*this) );
}
protocols::moves::MoverOP MetalSiteFinderMover::fresh_instance() const{
	return protocols::moves::MoverOP( new MetalSiteFinderMover);
}
std::string MetalSiteFinderMover::get_name() const{
	return "MetalSiteFinderMover";
}

void MetalSiteFinderMover::apply(core::pose::Pose & pose){
	TR << "Begin structure" << std::endl;
	//Get DSSP for pose
	core::scoring::dssp::Dssp pose_dssp_object( pose );


	//Create a dummy pose that will be used for pdb dumps
	core::pose::PoseOP dummy_pose;

	//Set max_res_span and min_local_coord_atoms and min_coord_res from options
	core::Size max_res_span = basic::options::option[ local::max_res_span ].value();
	core::Size min_total_coord_atoms = basic::options::option[ local::min_total_coord_atoms ].value();
	core::Size min_local_coord_atoms = basic::options::option[ local::min_local_coord_atoms ].value();
	core::Size min_coord_res = basic::options::option[ local::min_coord_res ].value();
	core::Size max_loop_size = basic::options::option[ local::max_loop_size ].value(); // ClangSA: Value stored to 'max_loop_size' during its initialization is never read
	std::string allowed_elements = basic::options::option[ local::allowed_elements ].value();
	utility::vector1< std::string > allowed_elements_list = utility::string_split_simple( allowed_elements, ',' );
	bool strict = basic::options::option[ local::strict_ss_changes ].value(); // ClangSA: Value stored to 'strict' during its initialization is never read
	core::Real distance_cutoff = basic::options::option[ local::metal_ligand_distance_cutoff ].value(); //All known calcium ligands are 1.5 to 3.5 angstroms from the ion
	//Declare the AtomGraphOP pose_atom_graph

	core::conformation::AtomGraphOP pose_atom_graph( new core::conformation::AtomGraph );


	//Fill the AtomGraph

	core::conformation::atom_graph_from_conformation( pose.conformation(), *pose_atom_graph );
	if ( !pose_atom_graph->get_vertex( 1 ).data().residue_id() ) {
		TR << "AtomGraph did NOT properly initialize." << std::endl;
	} else {
		TR << "AtomGraph loaded data" << std::endl;
	}
	//The AtomGraph is an UpperEdgeGraph

	///TODO BEGIN  ================================================================

	//Find neighboring atoms within cutoff using find_neighbors
	core::conformation::find_neighbors< core::conformation::AtomGraphVertexData, core::conformation::AtomGraphEdgeData >( pose_atom_graph, distance_cutoff );
	//Now there are edges b/w atoms separated by no more distance_cutoff
	core::Size num_atoms = pose_atom_graph->num_vertices();

	//Find all metal residues
	//For each metal
	core::Size ca_site_number = 0;
	for ( core::Size ca_resnum = 1; ca_resnum <= pose.total_residue(); ++ca_resnum ) {
		//TR << "Checking residue " << ca_resnum << " " << pose.residue( ca_resnum ).name() << std::endl;
		if ( pose.residue( ca_resnum ).name() == basic::options::option[ local::metal_type ].value() ) {
			TR << "Found metal" << std::endl;
			//This residue should only contain one atom (the metal atom)
			ca_site_number += 1;
			core::Size calcium_index = 0;
			for ( core::Size i = 1; i <= num_atoms; ++i ) {
				//Find the index for this metal
				//TR << "Current residue: " <<  pose_atom_graph->get_vertex( i ).data().residue_id() << std::endl;
				if ( pose_atom_graph->get_vertex( i ).data().residue_id() == ca_resnum ) {
					calcium_index = i;
					TR << "Found metal atom at " << calcium_index << std::endl;
					break;//There should only be 1 atom in this residue
				}
			}
			runtime_assert( calcium_index != 0);
			//Print the total number of metal ligands
			TR << "Ligands for metal at position " << ca_resnum << ": " << pose_atom_graph->get_vertex( calcium_index ).num_neighbors() << std::endl;
			//Find all ligand atoms
			utility::vector1< core::id::AtomID > ligand_atoms;

			//This will use the AtomGraph. What is the position of our metal atom in this graph?
			//The graph knows atom names & residue positions ( atom_name() and residue_id() )

			//First, get all the edges that connect FROM another atom TO this metal (can't easily do w/graph setup)
			utility::graph::UEVertex< core::conformation::AtomGraphVertexData, core::conformation::AtomGraphEdgeData > current_lower_atom;
			for ( core::Size lower_atom_num = 1; lower_atom_num < calcium_index; ++lower_atom_num ) {
				current_lower_atom = pose_atom_graph->get_vertex( lower_atom_num );

				if ( current_lower_atom.edge_exists( calcium_index ) ) {
					ligand_atoms.push_back( core::id::AtomID(  pose.residue( current_lower_atom.data().residue_id() ).atom_index( current_lower_atom.data().atom_name() ) , current_lower_atom.data().residue_id() ) );
				}

			}
			//Now get all the edges that connect FROM this metal TO another atom
			for ( core::conformation::AtomGraph::UpperEdgeListConstIter iter = pose_atom_graph->get_vertex( calcium_index ).const_upper_edge_list_begin(); iter != pose_atom_graph->get_vertex( calcium_index ).const_upper_edge_list_end(); ++iter ) {
				ligand_atoms.push_back( core::id::AtomID(  pose.residue( pose_atom_graph->get_vertex( iter->upper_vertex() ).data().residue_id() ).atom_index( pose_atom_graph->get_vertex( iter->upper_vertex()).data().atom_name() ) , pose_atom_graph->get_vertex( iter->upper_vertex() ).data().residue_id() ) );
			}




			///TODO END ================================================================



			runtime_assert(ligand_atoms.size() ==  pose_atom_graph->get_vertex(calcium_index ).num_neighbors() );
			//If the metal has no ligands, move on to the next one
			if ( ligand_atoms.size() == 0 ) {
				continue;
			}
			//Toss out anything that isn't a metal coordinating atom
			utility::vector1< core::id::AtomID > ligand_oxygens;
			ligand_oxygens.clear();
			//get_metalbinding_atoms modifies a vector of core::Size to return indices of all metal binding atoms
			utility::vector1< core::Size > possible_ligand_atom_indices;

			for ( core::Size ligand_no = 1; ligand_no <= ligand_atoms.size(); ++ligand_no ) {
				possible_ligand_atom_indices.clear();
				pose.residue( ligand_atoms[ ligand_no ].rsd() ).type().get_metal_binding_atoms( possible_ligand_atom_indices );
				if ( std::find( possible_ligand_atom_indices.begin(), possible_ligand_atom_indices.end(), ligand_atoms[ ligand_no ].atomno() ) != possible_ligand_atom_indices.end() ) {
					//Remove this entry from ligand_atoms
					ligand_oxygens.push_back( ligand_atoms[ ligand_no ] );
				}
			}
			TR << "Number of coordinating atoms for metal at position " << ca_resnum << ": " << ligand_oxygens.size() << std::endl;
			//A residue may appear more than once if it is bidentate
			//This vector should be sorted--first by residue & then by atom. Comparison operators for AtomID work this way.
			//If there are no coordinating oxygens, move on
			if ( ligand_oxygens.size() == 0 ) {
				continue;
			}
			/*
			//Sort ligand_oxygens
			utility::vector1< core::Size >  index_order( ligand_oxygens.size() );
			utility::arg_least_several(ligand_oxygens, index_order ); //Should return indices in order least to greatest
			//Make temporary copy for sorting purposes
			utility::vector1< core::id::AtomID > temp_copy;
			for (core::Size ii = 1; ii <= index_order.size(); ++ii ){
			temp_copy.push_back( ligand_oxygens[ index_order[ ii ] ] );
			}*/

			//Make temporary copy for sorting purposes
			utility::vector1< core::id::AtomID > temp_copy;
			core::Size min_index;
			//core::Size previous_min_index = 0;
			while ( ligand_oxygens.size() > 0 ) {
				min_index = utility::arg_min( ligand_oxygens );
				temp_copy.push_back( ligand_oxygens[ min_index ] );
				ligand_oxygens.erase( ligand_oxygens.begin() + min_index - 1 );
			}
			TR << "Confirm number of ligand atoms: " << temp_copy.size() << std::endl;
			ligand_oxygens = temp_copy;




			//Check if it has the required # of coordinating residues
			core::Size num_coord_res = 0;
			core::Size last_checked_res = 0;
			for ( core::Size ligand_atom_no = 1; ligand_atom_no <= ligand_oxygens.size(); ++ligand_atom_no ) {
				if ( ligand_oxygens[ ligand_atom_no ].rsd() != last_checked_res ) {
					++num_coord_res;
					last_checked_res = ligand_oxygens[ ligand_atom_no ].rsd();
				}
			}
			if ( num_coord_res < min_coord_res ) { //If the binding site contains too few residues, move on
				TR << "The binding site contains too few residues!" << std::endl;
				continue;
			}




			//Make another vector of the same length.
			utility::vector1< core::Size > num_coord_oxygens;



			for ( core::Size ox_num = 1; ox_num <= ligand_oxygens.size(); ++ox_num ) {
				core::Size current_resnum = ligand_oxygens[ ox_num ].rsd();
				core::Size current_local_ligands = 1;
				//Loop through all downstream residues
				for ( core::Size next_ox_num = ox_num + 1; next_ox_num <= ligand_oxygens.size(); ++next_ox_num ) {
					//Count how many downstream coordinating residues are within max_res_span
					if ( next_ox_num > ligand_oxygens.size() ) {
						TR << "Faulty loop for next_ox_num" << std::endl;
					}
					if ( ligand_oxygens[ next_ox_num ].rsd() < (current_resnum + max_res_span ) ) {
						++current_local_ligands;
					}
				}
				num_coord_oxygens.push_back( current_local_ligands );
			}
			if ( num_coord_oxygens.size() < min_total_coord_atoms ) {
				TR << "Not enough nearby atoms are ligand atoms." << std::endl;
				continue;
			}
			//Make sure we have enough ligands
			if ( utility::max( num_coord_oxygens ) < min_local_coord_atoms ) {
				TR << "Not enough local ligand atoms." << std::endl;
			}
			{
				continue; //Move on to the next metal
			}
			//Find which residue has the most local ligands
			core::Size max_index = utility::arg_max( num_coord_oxygens );
			//Span = [residue in coord_res with same index, coord_res[index + # - 1]
			core::Size span_start = ligand_oxygens[ max_index ].rsd();
			core::Size span_end = ligand_oxygens[ max_index + num_coord_oxygens[ max_index ] - 1 ].rsd();
			TR << "Span start " << span_start << " Span end " << span_end << std::endl;

			//Make note of whether the span includes all of the coordinating residues
			bool local_ca_site = ( ligand_oxygens.size() == num_coord_oxygens[ max_index ] );
			utility::vector1<core::id::AtomID> nonlocal_contacts;
			if ( !local_ca_site ) {
				for ( core::Size ligand_oxy_index = 1; ligand_oxy_index < max_index; ++ligand_oxy_index ) {
					nonlocal_contacts.push_back( ligand_oxygens[ ligand_oxy_index ] );
				}
				for ( core::Size ligand_oxy_index = max_index + utility::max( num_coord_oxygens ); ligand_oxy_index <= ligand_oxygens.size() ; ++ligand_oxy_index ) {
					nonlocal_contacts.push_back( ligand_oxygens[ ligand_oxy_index ] );
				}
				TR << "Found : " << nonlocal_contacts.size() << " nonlocal contacts: ";
				for ( core::Size nonlocal_contact_index = 1; nonlocal_contact_index <= nonlocal_contacts.size(); ++nonlocal_contact_index ) {
					TR << nonlocal_contacts[nonlocal_contact_index].rsd() << " ";
				}
				TR << std::endl;
			}


			//Get DSSP secondary structure for span residue 1
			std::string secstruct = pose_dssp_object.get_dssp_reduced_IG_as_L_secstruct();
			TR << pose.pdb_info()->name() << " Secondary structure: " << secstruct << std::endl;
			//If = H or E
			//NOTE: Should require a 2-residue consistent change in DSSP to constitute beginning of loop in SSE--avoid 1-residue loops/strands
			core::Size sse_begin = 0;
			core::Size sse_end = 0;
			core::Size site_end = 0;
			core::Size next_sse_begin = 0;
			core::Size chain_start = pose.conformation().chain_begin(pose.residue(span_start).chain());
			core::Size chain_stop = pose.conformation().chain_end(pose.residue(span_start).chain());

			if ( secstruct[ span_start - 1] != 'L' ) { //The first coordinating residue is already in an SSE, so must still find sse_end
				//Find first residue in that SSE


				if ( strict ) {
					for ( core::Size ss_res = span_start; ss_res > chain_start + 1; --ss_res ) {
						//Check if ss_res - 1 AND ss_res - 2 are a different SSE (i.e. ss_res is the first residue in the same SSE as span_start)
						if ( secstruct[ span_start - 1 ] != secstruct[ ss_res - 2 ] && secstruct[ ss_res - 3 ] !=  secstruct[ span_start - 1 ] ) {
							//Store that resnum as sse_begin
							sse_begin = ss_res;
							break;
						} else if ( ss_res == chain_start + 2 && secstruct[ ss_res - 3 ] != secstruct[ span_start - 1] ) {

							sse_begin = ss_res - 1;

						} else if ( ss_res == chain_start + 2 ) {
							sse_begin = chain_start;
							TR << "Could not locate beginning of current SSE. Setting sse_begin to first protein residue." << std::endl;
						}
					}
				} else { //Now we're just checking if ss_res - 1 is a different SSE than ss_res
					for ( core::Size ss_res = span_start; ss_res > chain_start; --ss_res ) {

						if ( secstruct[ span_start - 1] != secstruct[ ss_res - 2 ] ) {
							sse_begin = ss_res;
							break;
						} else if ( ss_res == chain_start ) {
							sse_begin = chain_start;
							TR << "Could not locate beginning of current SSE. Setting sse_begin to first protein residue." << std::endl;
						}
					}
				}
				//Find the end of that SSE (for loop size)






				//Detecting end of SSE that contains the first residue

				if ( strict ) {
					//In the strict case, both ss_res and ss_res - 1 (two residues after the possible SSE) need to be different from the SSE of span_start

					for ( core::Size ss_res = span_start + 2; ss_res <= chain_stop; ++ss_res ) {
						if ( secstruct[ span_start - 1] != secstruct[ ss_res - 2 ] && secstruct[ span_start - 1] != secstruct[ ss_res - 1] ) {

							sse_end = ss_res - 2;
							break;
						} else if ( ss_res == chain_stop ) { //We have reached the end of the pose, and at least 1 of the last 2 residues still matches our SSE
							sse_end = chain_stop;
						}
					}
				} else {
					for ( core::Size ss_res = span_start + 1; ss_res <= chain_stop; ++ss_res ) {
						//Now we only really need to check that ss_res is different from span_start
						if ( secstruct[ span_start - 1 ] != secstruct[ ss_res - 1 ] ) {
							sse_end = ss_res - 1;
							break;
						} else if ( ss_res == chain_stop ) {
							sse_end = chain_stop;
						}
					}
				}




			} else { //If span_start is on a loop


				//NOTE: In this case, even when "strict" is not enabled, I still require at least two residues of the same non-loop SSE type to be considered an SSE



				for ( core::Size ss_res = span_start; ss_res > chain_start + 1; --ss_res ) {
					//Loop backwards to the previous SSE
					//Detecting end of SSE (i.e. first non-loop residue) before the first Ca-coordinating residue
					if ( secstruct[ ss_res - 2] != 'L' && secstruct[ ss_res - 3] == secstruct[ ss_res - 2 ] ) {
						sse_end = ss_res - 1;
						break;
					} else if ( ss_res == chain_start + 2 && secstruct[ chain_start] != 'L' && secstruct[chain_start] == secstruct[chain_start - 1] ) { //If the first 2 residues constitute an SSE
						sse_end = chain_start + 1;
						break;
					}
				}

				//If we reached the beginning of the pose without finding an SSE, just set beginning of motif to the first non-CA residue
				core::Size current_res_num = chain_start;
				while ( sse_end == 0 && sse_begin == 0 && current_res_num <= span_start ) {
					if ( pose.residue( current_res_num ).name() != "CA" ) {
						sse_begin = current_res_num;
						TR << "Could not locate previous SSE. Setting sse_begin to first protein residue." << std::endl;
					}
					++current_res_num;
				}
				if ( sse_begin == 0 && sse_end == 0 ) {
					TR << "Could not locate previous SSE. Setting sse_begin to span_start." << std::endl;
					sse_begin = span_start;
				}



				//If span_start is in a loop *and* we found an SSE before it, we need to find where that SSE begins
				if ( sse_end != 0 ) {
					//Else find the first residue in the previous SSE (just before the DSSP changes)


					if ( strict ) {
						for ( core::Size prev_ss_res = sse_end; prev_ss_res > chain_start + 1; --prev_ss_res ) {
							if ( secstruct[ sse_end - 1] != secstruct[ prev_ss_res - 2] && secstruct[ sse_end - 1] != secstruct[ prev_ss_res - 3] ) {
								sse_begin = prev_ss_res;
								break;
							} else if ( prev_ss_res == chain_start + 2 && secstruct[ prev_ss_res - 2 ] == secstruct[ sse_end - 1] && secstruct[ prev_ss_res - 3] != secstruct[ sse_end - 1] ) {
								//In the special case where only the first residue in the pose is a mismatch, we'll start the SSE at residue 2
								sse_begin = prev_ss_res - 1;
								break;
							} else if ( prev_ss_res == chain_start + 2 ) {
								sse_begin = chain_start;
							}
						}
					} else {
						for ( core::Size prev_ss_res = sse_end; prev_ss_res > 1; --prev_ss_res ) {
							if ( secstruct[ sse_end - 1] != secstruct[ prev_ss_res - 2] ) {
								sse_begin = prev_ss_res;
								break;
							}
						}
					}
					//If we didn't find sse_begin, set it to 1 (Ligands will never be in an SSE)
					if ( sse_begin == 0 ) {
						sse_begin = chain_start;
					}
				}
			}
			//At this point, the only case where sse_end should be 0 is if we didn't find any SSE before this span
			//In that case, set the value to sse_begin - 1, since sse_begin will just be the beginning of the loop (may still be zero, but there may be metals)
			if ( sse_end == 0 ) {
				sse_end = sse_begin - 1;
			}

			//Now we need to find the end of the stretch that we will take

			//If span_end is on an SSE, go to the end of that SSE and find the beginning of that SSE for loop length







			//Maybe have some flag that tells us whether sse_end is greater than span_end--if it is, this SSE contains the entire site, and we need to figure out whether to move forward or backwards to take the next one


			if ( sse_end > span_end ) { //This has to mean that span_end is on the same sse as span_start
				////In this case, we can just output this helix!!!
				site_end = sse_end;
				/*
				//Step 1) Find which distance is smaller, span_start to sse_begin or sse_end to span_end
				if ( (span_start - sse_begin ) >= (span_end - sse_end ) ) {
				//Step 2) If the first distance is smaller, we actually want to move backwards. Set site_end to sse_end, next_sse_begin to sse_begin, sse_begin and sse_end to zero, and find the loop and SSE before them.

				site_end = sse_end;
				next_sse_begin = sse_begin;
				sse_begin = 0;
				sse_end = 0;

				//Starting at one residue before next_sse_begin and going back to the beginning of the pose, loop through and look for 2 consecutive residues of the same secondary structure
				for ( core::Size prev_loop_res = next_sse_begin - 2; prev_loop_res >= chain_start; ++prev_loop_res ) {
				if ( secstruct[ prev_loop_res - 1] != 'L' && secstruct[ prev_loop_res - 1] == secstruct[ prev_loop_res ] ) {
				sse_end = prev_loop_res + 1;
				}
				//This time, if we reach the beginning of the pose without finding another SSE, we just want to move on, so no special case here
				}
				//Now, if we *did* find the end of a previous SSE, we want to now find where that SSE began.
				if ( strict && sse_end > 0 ) {
				for ( core::Size prev_ss_res = sse_end; prev_ss_res > chain_start + 1; --prev_ss_res ) {
				if ( secstruct[ sse_end - 1] != secstruct[ prev_ss_res - 2] && secstruct[ sse_end - 1 ] != secstruct[ prev_ss_res - 3] ) {
				sse_begin = prev_ss_res;
				break;
				} else if ( prev_ss_res == chain_start + 2 && secstruct[ prev_ss_res - 2 ] == secstruct[ sse_end - 1 ] && secstruct[ prev_ss_res - 3] != secstruct[ sse_end - 1 ] ) {
				//In the special case where only the first residue in the pose is a mismatch, we'll start the SSE at residue 2
				sse_begin = prev_ss_res - 1;
				break;
				} else if ( prev_ss_res == chain_start + 2 ) {
				sse_begin = chain_start;
				}
				}
				} else if ( sse_end > 0 ) {
				for ( core::Size prev_ss_res = sse_end; prev_ss_res > chain_start; --prev_ss_res ) {
				if ( secstruct[ sse_end - 1] != secstruct[ prev_ss_res - 2] ) {
				sse_begin = prev_ss_res;
				break;
				}
				}
				}

				//2b) If we don't find an SSE before this one, then we'll find one afterwards instead. Set sse_begin to next_sse_begin, sse_end to site_end, and set next_sse_end and site_end back to zero.
				if ( sse_end == 0 && sse_begin == 0 ) {
				sse_end = site_end;
				sse_begin = next_sse_begin;
				site_end = 0;
				next_sse_begin = 0;
				}
				}



				//Step 3) If next_sse_begin and site_end are still zero, then we need to find an SSE after this one.
				if ( next_sse_begin == 0 && site_end == 0 ) {
				//We know that sse_end is greater than span_end
				//Starting at sse_end + 2
				for ( core::Size next_loop_res = sse_end + 2; next_loop_res <= chain_stop ;++next_loop_res ) {
				if ( secstruct[ next_loop_res - 1] != 'L' && secstruct[ next_loop_res - 1] == secstruct[ next_loop_res - 2 ] ) {
				next_sse_begin = next_loop_res - 1;
				}
				}

				//Now we need to find where this last sse ends if we found one.
				if ( next_sse_begin > 0 && strict ) {
				for ( core::Size next_ss_res = next_sse_begin + 2; next_ss_res <= chain_stop; ++next_ss_res ) {
				//Normal case
				if ( secstruct[ next_ss_res - 1] != secstruct[ next_sse_begin - 1] && secstruct[ next_ss_res - 2] != secstruct[ next_sse_begin - 1] ) {
				site_end = next_ss_res - 2;
				} else if ( next_ss_res == chain_stop && secstruct[ next_ss_res - 2] == secstruct[ next_sse_begin - 1 ] && secstruct[ next_ss_res - 1 ] != secstruct[ next_sse_begin - 1] ) {
				//Special case where only the last residue in the pose doesn't match
				site_end = next_ss_res - 1;
				}

				}
				} else if ( next_sse_begin > 0 ) {
				for ( core::Size next_ss_res = next_sse_begin + 1; next_ss_res <= chain_stop; ++next_ss_res ) {
				if ( secstruct[ next_ss_res - 1 ] != secstruct[ next_sse_begin - 1 ] ) { //If the ss of this residue does not match our SSE, then the previous residue was the last residue in the SSE
				site_end = next_ss_res - 1;
				}
				}
				}
				//If we didn't find one, just set site_end to the end of the pose.
				if ( site_end == 0 ) {
				site_end = chain_stop;
				}
				}*/
			} else if ( secstruct[ span_end  - 1 ] != 'L' ) { //In this case, span_end is in a different SSE from span_begin

				//Find last residue in that SSE
				if ( strict ) {
					for ( core::Size ss_res = span_end + 2; ss_res <= chain_stop; ++ss_res ) {
						if ( secstruct[ span_end  - 1 ] != secstruct[ ss_res - 2 ] && secstruct[ span_end - 1 ] != secstruct[ ss_res - 1 ] ) {
							site_end = ss_res - 2;
							break;
						} else if ( ss_res == chain_stop  && secstruct[ ss_res - 1 ] != secstruct[ span_end - 1 ] ) {
							site_end = chain_stop - 1;
						} else if ( ss_res == chain_stop ) {
							site_end = chain_stop;
						}
					}
				} else {
					for ( core::Size ss_res = span_end + 1; ss_res <= chain_stop; ++ss_res ) {
						if ( secstruct[ span_end - 1 ] != secstruct[ ss_res - 1 ] ) {
							site_end = ss_res - 1;
							break;
						} else if ( ss_res == chain_stop ) {
							site_end = chain_stop;
						}
					}
				}
				//Now find the first residue in that SSE
				if ( strict ) {
					for ( core::Size ss_res = span_end; ss_res > chain_start + 1; --ss_res ) {
						//Check if ss_res - 1 is a different SSE (i.e. ss_res is the first residue in the same SSE as span_start)
						if ( secstruct[ span_end - 1] != secstruct[ ss_res - 2 ] && secstruct[ ss_res - 3 ] !=  secstruct[ span_end - 1] ) {
							//Store that resnum as sse_begin
							next_sse_begin = ss_res;
							break;
						} else if ( ss_res == chain_start + 2 && secstruct[ ss_res - 2 ] == secstruct[ sse_end - 1 ] && secstruct[ ss_res - 3] != secstruct[ sse_end - 1 ] ) {
							//In the special case where only the first residue in the pose is a mismatch, we'll start the SSE at residue 2
							next_sse_begin = ss_res - 1;
							break;
						} else if ( ss_res == chain_start + 2 ) {
							next_sse_begin = chain_start;
						}
					}
				} else {
					for ( core::Size ss_res = span_end; ss_res > chain_start; --ss_res ) {
						if ( secstruct[ span_end - 1 ] != secstruct[ ss_res - 2 ] ) {
							next_sse_begin = ss_res;
							break;
						} else if ( ss_res == chain_start + 1 ) {
							next_sse_begin = chain_start;
						}
					}
				}
			} else { //If the last residue is on a loop
				next_sse_begin = 0;
				//To count as an SSE, we must have at least 2 consecutive residues of the same SS type
				for ( core::Size ss_res = span_end + 2; ss_res <= chain_stop; ++ss_res ) {
					if ( secstruct[ ss_res - 2 ] != 'L' && secstruct[ ss_res - 1 ] == secstruct[ ss_res - 2 ] ) { //If we reach something that is no longer a loop
						next_sse_begin = ss_res - 1;
						break;
					}
				}
				//If we didn't find an SSE after this residue, then set site_end to the last non-CA residue
				core::Size current_res_num = chain_stop;
				while ( next_sse_begin == 0 && site_end == 0 && current_res_num > span_end ) {
					if ( pose.residue( current_res_num ).name() != basic::options::option[ local::metal_type ].value() ) {
						site_end = current_res_num;
					}
					--current_res_num;
				}
				//If we got back to span_end and still didn't find one, set it to span_end
				if ( next_sse_begin == 0 && site_end == 0 ) {
					site_end = span_end;
				}
				//In cases where we DID find another SSE after this loop
				if ( next_sse_begin != 0 ) {
					//Find last residue in this SSE
					if ( strict ) {
						for ( core::Size next_ss_res = next_sse_begin + 1; next_ss_res <= chain_stop - 1; ++next_ss_res ) {
							if ( secstruct[ next_ss_res - 1 ] != secstruct[ next_sse_begin - 1 ] && secstruct[ next_ss_res ] != secstruct[ next_sse_begin - 1 ] ) { //If we have moved on to a new SSE
								site_end = next_ss_res - 1;
								break;
							} else if ( next_ss_res == chain_stop - 1 && secstruct[ next_ss_res  - 1 ] == secstruct[ next_sse_begin - 1 ] ) {
								//Special case where we reach the end and *only* the last residue in the pose doesn't match
								site_end = next_ss_res;
							} else if ( next_ss_res == chain_stop ) {
								site_end = chain_stop;
							}
						}
					} else {
						for ( core::Size next_ss_res = next_sse_begin + 1; next_ss_res <= chain_stop; ++next_ss_res ) {
							if ( secstruct[ next_ss_res - 1 ] != secstruct[ next_sse_begin - 1 ] ) {
								site_end = next_ss_res - 1;
								break;
							} else if ( next_ss_res == chain_stop ) {
								site_end = chain_stop;
							}
						}
					}
					//If we didn't find it, set to the last residue in the pose
					if ( site_end == 0 ) {
						site_end = chain_stop;
					}
				}
			}
			//At this point, next_sse_begin should only be 0 if there is no sse after this span, so set it to site_end + 1
			if ( next_sse_begin == 0 ) {
				next_sse_begin = site_end + 1;
			}
			//Now we know where both the span (residues which must be conserved) and the site (all residues to take) begin and end.

			//Dump the corresponding stretch of residues to a PDB with the original pdb name + _ca_ + (local_ or split_ ) +  site number
			std::string output_name = pose.pdb_info()->name();
			core::Size name_end = output_name.rfind(".pdb");//Returns index of "."
			output_name = output_name.substr(0, name_end ) + "_ca_";
			//Convert site number to string
			std::stringstream ss;
			ss << ca_site_number;
			std::string siteno = ss.str();
			if ( local_ca_site ) {
				output_name += "local_" + siteno + ".pdb";
			} else {
				output_name += "split_" + siteno + ".pdb";
			}
			//Filter on loop length here
			//Check that the loop is short enough
			core::Size loop_size;
			if ( sse_end >= next_sse_begin ) {
				loop_size = 0;
			} else {
				loop_size = next_sse_begin - sse_end - 1;
			}
			if ( loop_size > max_loop_size ) {
				TR << output_name << " rejected for loop size greater than maximum. Size: " << loop_size <<  std::endl;
				continue; //Move on to the next ca site
			}
			///ADD CHECK HERE FOR HELICES
			//Possibilities we want to allow output for:
			//No SS required
			//SS required but no helices required and both ends are not loops
			//Helices required and both ends are helices
			if ( !( basic::options::option[ local::require_bordering_ss ].value() ) || ( !( basic::options::option[ local::helices_only ].value() ) && secstruct[ sse_begin - 1 ] != 'L' && secstruct[ site_end - 1 ] != 'L' && next_sse_begin < site_end) || ( secstruct[ sse_begin - 1 ] == 'H' && secstruct[ site_end - 1 ] == 'H' && next_sse_begin < site_end) ) { //Check to make sure there is bordering SS if required



				//All the local contacts must be of the allowed type
				//If empty string (default), ignore
				//Otherwise only allow contacts of the specified element type
				core::Size good_local_ligands = 0;
				//bool site_failed = false;

				for ( core::Size ligand_no = 1; ligand_no <=ligand_oxygens.size(); ++ligand_no ) {
					std::string element = pose.residue( ligand_oxygens[ ligand_no ].rsd() ).atom_type_set()[ pose.residue( ligand_oxygens[ ligand_no ].rsd() ).atom_type_index( ligand_oxygens[ ligand_no ].atomno() ) ].element();
					//If this element is NOT allowed, skip it
					if ( allowed_elements != "" && std::find( allowed_elements_list.begin(), allowed_elements_list.end(), element ) == allowed_elements_list.end() ) {
						continue;
					}
					++good_local_ligands;
				}
				if ( good_local_ligands <= min_local_coord_atoms ) {
					TR << "Not enough local ligands!" << std::endl;
					continue;
				}

				//TR << output_name << " Begin " << sse_begin << " End " << site_end << " First coordinating " << span_start << " Last coordinating " << span_end << " Number coordinating " <<  utility::max( num_coord_oxygens ) << std::endl;
				TR << output_name << "Site in pose begins at " << sse_begin << " and ends at " << site_end << std::endl;
				TR << output_name << " First coordinating " << span_start - sse_begin + 1 << " Last coordinating " << span_end - sse_begin + 1 << " Number coordinating " <<  utility::max( num_coord_oxygens ) << std::endl;
				TR << output_name << " Secondary structure " << secstruct.substr(sse_begin - 1, site_end - sse_begin + 1) << std::endl;
				//Recall: sse_end gives the end of the previous SSE, and next_sse_begin gives the beginning of the next SSE
				TR << output_name << " Loop length " << next_sse_begin - sse_end << std::endl;
				//Now print all the coordinating residue numbers by the new residue numbering system
				TR << output_name << " Coordinating residues ";
				core::Size residue_number = 0;
				for ( core::Size res_id = 1; res_id <=ligand_oxygens.size(); ++res_id ) {
					if ( residue_number == ligand_oxygens[ res_id ].rsd() ) {
						continue;
					} else {
						residue_number = ligand_oxygens[ res_id ].rsd();
						TR << residue_number << " ";
					}
				}
				TR << std::endl;
				dummy_pose = core::pose::PoseOP(new core::pose::Pose(pose, sse_begin, site_end) );
				//core::pose::add_lower_terminus_type_to_pose_residue( *dummy_pose, 1 );
				//TR << "Added lower terminus to site" << std::endl;
				//core::pose::add_upper_terminus_type_to_pose_residue( *dummy_pose, dummy_pose->total_residue() );
				//TR << "Added upper terminus to site" << std::endl;
				//We also want to include the metal in this though
				core::conformation::ResidueOP ca_residue_copy = pose.residue(ca_resnum).clone();
				dummy_pose->append_residue_by_jump(*ca_residue_copy, (span_start - sse_begin + 1) );
				if ( basic::options::option[ local::append_nonlocal_residues ].value() ) {
					std::set<core::Size> nonlocal_residues;
					for ( core::Size nonlocal_cont_id = 1; nonlocal_cont_id <= nonlocal_contacts.size(); ++nonlocal_cont_id ) {
						nonlocal_residues.insert(nonlocal_contacts[nonlocal_cont_id].rsd());
					}

					for ( std::set<core::Size>::iterator nonlocal_itr = nonlocal_residues.begin(); nonlocal_itr != nonlocal_residues.end(); ++nonlocal_itr ) {
						core::conformation::ResidueOP nonlocal_residue_copy = pose.residue(*nonlocal_itr).clone();
						if ( nonlocal_itr == nonlocal_residues.begin() ) {
							dummy_pose->append_residue_by_jump(*nonlocal_residue_copy, (span_start - sse_begin + 1), "", "", true );
						} else {
							dummy_pose->append_residue_by_jump(*nonlocal_residue_copy, dummy_pose->total_residue());
						}
					}

					/*core::pose::add_lower_terminus_type_to_pose_residue( *dummy_pose, dummy_pose->total_residue() );
					TR << "Added lower terminus to nonlocal residue " << *nonlocal_itr << std::endl;
					core::pose::add_upper_terminus_type_to_pose_residue( *dummy_pose, dummy_pose->total_residue() );
					TR << "Added upper terminus to nonlocal residue " << *nonlocal_itr << std::endl;*/
				}
				dummy_pose->dump_pdb(output_name);
				//We also want to include the metal in this though
				/*core::conformation::ResidueOP ca_residue_copy = pose.residue(ca_resnum).clone();
				dummy_pose->append_residue_by_jump(*ca_residue_copy, (span_start - sse_begin + 1) );
				dummy_pose->dump_pdb(output_name);*/
			}
		}
		//Repeat for each metal site
	}
}
int main( int argc, char* argv[] ){
	try{
		basic::options::option.add( local::metal_type, "Residue name for the metal you want the app to search for (e.g. \"CA\" for calcium" );
		basic::options::option.add( local::metal_ligand_distance_cutoff, "Maximum distance between metal atom and metal-coordinating atom for a contact to be detected" ).def( 3.5 );
		basic::options::option.add( local::min_total_coord_atoms, "Minimum number of coordinating oxygen atoms to constitute a calcium binding site").def(5);
		basic::options::option.add( local::min_local_coord_atoms, "Minimum number of local coordinating oxygen atoms to constitute a calcium binding site" ).def(3);
		basic::options::option.add( local::min_coord_res, "Minimum number of residues involved in coordination" ).def(3);
		basic::options::option.add( local::max_res_span, "Maximum length for a calcium binding site sequence" ).def(20);
		basic::options::option.add( local::require_bordering_ss, "Requires sites to begin and end with an SSE" ).def( true );
		basic::options::option.add( local::max_loop_size, "Maximum separation between 2 secondary structure elements" ).def(12);
		basic::options::option.add( local::helices_only, "Bordering secondary structure elements must be helices. Only takes effect if require_bordering_ss is also enabled" ).def( false );
		basic::options::option.add( local::strict_ss_changes, "Detects change in secondary structure only after 2 consistently changed residues" ).def( false );
		basic::options::option.add( local::append_nonlocal_residues, "Output nonlocal contacts as a separate chain in the same PDB file" ).def( false );
		basic::options::option.add( local::allowed_elements, "Restrict metal binding atoms to the specified elements. By default, all potential metal binding atoms are considered." ).def( "" );
		//basic::options::option.add( local::mc_temperature, "Temperature to use in MonteCarlo object" ).def(0.6);
		devel::init(argc, argv);
		MetalSiteFinderMoverOP ca_finder (new MetalSiteFinderMover );
		//Apply using JD2
		protocols::jd2::JobDistributor::get_instance()->go(ca_finder);
		return 0;
	}
catch ( utility::excn::Exception const &e ){
	std::cout << "Caught exception " << e.msg() << std::endl;
	return -1;
}
}
