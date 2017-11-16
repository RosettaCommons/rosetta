// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/task/operation/InterfaceVectorDefinition.cc
/// @brief  Calculates the residues at an interface between two protein chains or jump.
/// The calculation is done in the following manner.  First the point graph
/// is used to find all residues within some big cutoff of residues on the other chain.
/// For these residues near the interface, two metrics are used to decide if they are actually
/// possible interface residues.  The first metric is to iterate through all the side chain
/// atoms in the residue of interest and check to see if their distance is less than the nearby
/// atom cutoff, if so then they are an interface residue.  If a residue does not pass that
/// check, then two vectors are drawn, a CA-CB vector and a vector from CB to a CB atom on the
/// neighboring chain.  The dot product between these two vectors is then found and if the angle
/// between them is less than some cutoff then they are classified as interface.
/// @author Ben Stranges (stranges@unc.edu)


// Unit headers
#include <core/select/util/interface_vector_calculate.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/kinematics/FoldTree.hh>
#include <ObjexxFCL/FArray1D.hh>

// Utility headers
#include <utility/string_util.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>
#include <cmath>
#include <set>

// option key includes


#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <utility/vector1.hh>

namespace core {
namespace select {
namespace util {

// typedefs
typedef std::pair< std::set<core::Size>,std::set<core::Size> > InterfacePair;
typedef numeric::HomogeneousTransform< core::Real > HTReal;

//static basic::Tracer TR("core.select.util.interface_vector_calculate");


//forward declarations of funtions that do the work
/// @brief looks at the big set and figures out what is actually pointing towards the interface
void find_interface_pointing_residues_from_neighbs(
	core::pose::Pose const & pose, InterfacePair const & interface_pair,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff,
	utility::vector1_bool & interface_residues
);

/// @brief find nearby atoms to other in interface
bool any_atoms_within_cutoff(core::conformation::Residue & res1,
	core::conformation::Residue & res2,
	core::Real cutoff);

/// @brief neighbors to look for vectors within (big set here)
InterfacePair find_neighbors_within_CB_cutoff( core::pose::Pose const & pose, core::Real big_cutoff,
	core::Size chain1, core::Size chain2 );

/// @brief find neighbors to look for vectors within using a big cutoff for CBs
InterfacePair find_jump_partners_within_CB_cutoff( core::pose::Pose const & pose, core::Real big_cutoff,
	int jump_num );

/// @brief the Cbeta vector(s) from on rsd to another
numeric::xyzVector<core::Real> cbeta_vector( core::conformation::Residue & res);
/// @brief the action coordinate for each residue
numeric::xyzVector<core::Real> select_coord_for_residue(core::conformation::Residue & res);
/// @brief out if res1 and res2 are pointing at eachother
bool res1_pointed_at_res2( core::conformation::Residue & res1,
	core::conformation::Residue & res2,
	core::Real angle_cutoff /*degrees*/,
	core::Real dist_cutoff);
void fill_in_chain_terminii( core::pose::Pose const & pose, core::Size chain1, core::Size chain2 );


/// @details minimal chain number definition
utility::vector1_bool
calc_interface_vector( core::pose::Pose const & pose, core::Size const chain1_number, core::Size const chain2_number ){
	// set some logical defaults and run the full calc function
	core::Real CB_dist_cutoff(10.0);
	core::Real nearby_atom_cutoff(5.5);
	core::Real vector_angle_cutoff(75.0);
	core::Real vector_dist_cutoff(9.0);

	return calc_interface_vector( pose, chain1_number, chain2_number, CB_dist_cutoff,
		nearby_atom_cutoff, vector_angle_cutoff, vector_dist_cutoff );
}

/// @details full runner that takes all of the inputs for chains
utility::vector1_bool
calc_interface_vector(
	core::pose::Pose  const & pose,
	core::Size const chain1_number, core::Size const chain2_number,
	core::Real const CB_dist_cutoff, core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff, core::Real const vector_dist_cutoff
){
	//set all residues in pose to false
	utility::vector1_bool at_interface(pose.size(), false);
	//do stuff
	//first find all the neighbors within some Cbeta cutoff distance from eachother
	InterfacePair CB_pairs_list;
	CB_pairs_list = find_neighbors_within_CB_cutoff(pose,  CB_dist_cutoff, chain1_number, chain2_number );
	//second, prune this set down to what matters
	find_interface_pointing_residues_from_neighbs( pose, CB_pairs_list, nearby_atom_cutoff,
		vector_angle_cutoff, vector_dist_cutoff, at_interface );
	// //debugging
	// for(core::Size ii = 1; ii<= at_interface.size(); ++ii)
	//  std::cout << "Residue number: " << ii << " at_interface value: " << at_interface[ii] << std::endl;

	return at_interface;
}

/// @details full runner that takes the jump
utility::vector1_bool
calc_interface_vector(
	core::pose::Pose const & pose,
	int const interface_jump,
	core::Real const CB_dist_cutoff,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff
){
	//set all residues in pose to false
	utility::vector1_bool at_interface(pose.size(), false);
	//do stuff
	//first find all the neighbors within some Cbeta cutoff distance from eachother
	InterfacePair CB_pairs_list;
	CB_pairs_list = find_jump_partners_within_CB_cutoff( pose, CB_dist_cutoff, interface_jump );
	//second, prune this set down to what matters
	find_interface_pointing_residues_from_neighbs( pose, CB_pairs_list, nearby_atom_cutoff,
		vector_angle_cutoff, vector_dist_cutoff, at_interface );
	// //debugging
	// for(core::Size ii = 1; ii<= at_interface.size(); ++ii)
	//  std::cout << "Residue number: " << ii << " at_interface value: " << at_interface[ii] << std::endl;

	return at_interface;
}

/// @details minimal jump runner
utility::vector1_bool
calc_interface_vector( core::pose::Pose const & pose, int const interface_jump ){
	// set some logical defaults and run the full calc function
	core::Real const CB_dist_cutoff(10.0);
	core::Real const nearby_atom_cutoff(5.5);
	core::Real const vector_angle_cutoff(75.0);
	core::Real const vector_dist_cutoff(9.0);
	return calc_interface_vector( pose, interface_jump, CB_dist_cutoff,
		nearby_atom_cutoff, vector_angle_cutoff, vector_dist_cutoff );
}

/// @details calc_interacting_vector does the same thing except does not need interface separation
/// I
utility::vector1_bool calc_interacting_vector(
	core::pose::Pose const & pose,
	std::set< core::Size > const & part1res,
	std::set< core::Size > const & part2res,
	core::Real const CB_dist_cutoff,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff )
{
	//set all residues in pose to false
	utility::vector1_bool at_interface(pose.size(), false);
	//first find all the neighbors within some Cbeta cutoff distance from eachother
	InterfacePair CB_pairs_list;

	//setup
	std::set<core::Size> side1_within_cutoff, side2_within_cutoff;
	//use neighbor atoms from point graph
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg);
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, CB_dist_cutoff );

	// for all nodes in chain1 == for all residues in chain 1
	// all this is setup by verify_chain_setup in InterfaceDefinitionBase
	for ( std::set<Size>::const_iterator side1_it = part1res.begin(); side1_it != part1res.end(); ++side1_it ) {
		for ( conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( *side1_it ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex( *side1_it ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			// get node on other edge of that node == 2nd residue index
			Size const edge_res = edge_iter->upper_vertex();
			// if that node(residue) is in the second set of residues
			if ( part2res.count( edge_res ) ) {
				side1_within_cutoff.insert( *side1_it ); // add partner1 residue
				side2_within_cutoff.insert( edge_res ); // add partner2 residue
			} else continue;
		} // end - for all edges of node
	} // end - for all nodes in chain1

	//return the pair of interface side sets
	CB_pairs_list = std::make_pair( side1_within_cutoff, side2_within_cutoff );

	//Now that we have the pairs list we can get the residues that are in proximity
	find_interface_pointing_residues_from_neighbs( pose, CB_pairs_list, nearby_atom_cutoff,
		vector_angle_cutoff, vector_dist_cutoff, at_interface );
	// //debugging
	// for(core::Size ii = 1; ii<= at_interface.size(); ++ii)
	//  std::cout << "Residue number: " << ii << " at_interface value: " << at_interface[ii] << std::endl;

	return at_interface;

}

/// @details does the real work, looks at the big set and figures out what is actually pointing towards the interface
///sets a vector bool value to true if a residue is at the interface
void find_interface_pointing_residues_from_neighbs(
	core::pose::Pose const & pose,
	InterfacePair const & interface_pairs,
	core::Real const nearby_atom_cutoff,
	core::Real const vector_angle_cutoff,
	core::Real const vector_dist_cutoff,
	utility::vector1_bool & interface_residues ){

	using namespace utility;
	using namespace core;
	//itterate over pairs
	for ( std::set<Size>::const_iterator chain1_it = interface_pairs.first.begin();  chain1_it!=interface_pairs.first.end(); ++chain1_it ) {
		conformation::Residue ch1residue(pose.residue(*chain1_it));
		for (  std::set<Size>::const_iterator chain2_it = interface_pairs.second.begin();  chain2_it!=interface_pairs.second.end(); ++chain2_it ) {
			conformation::Residue ch2residue(pose.residue(*chain2_it));
			//check to see if both are already in the interface
			if ( interface_residues[*chain1_it] && interface_residues[*chain2_it] ) {
				continue;
			}
			//use a stricter distance cutoff to define any very close to other chain
			bool r1_near_r2 (any_atoms_within_cutoff( ch1residue, ch2residue, nearby_atom_cutoff ) );
			bool r2_near_r1 (any_atoms_within_cutoff( ch2residue, ch1residue, nearby_atom_cutoff ) );
			//are residue1 sidechain atoms near residue2 atoms?
			if ( r1_near_r2 ) {
				interface_residues[*chain1_it] = true;
				// #ifndef NDEBUG
				//      std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain1_it) << " is included because Residue: "
				//           << pose.pdb_info()->pose2pdb(*chain2_it) <<" meets the nearby cutoff."<<std::endl;
				// #endif
			}
			//are residue2 sidechain atoms near residue1 atoms?
			if ( r2_near_r1 ) {
				interface_residues[*chain2_it] = true;
				// #ifndef NDEBUG
				//      std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain2_it) << " is included because Residue: "
				//           << pose.pdb_info()->pose2pdb(*chain1_it) <<" meets the nearby cutoff."<<std::endl;
				// #endif
			}
			//check to see if both are NOW at the interface
			if ( interface_residues[*chain1_it] && interface_residues[*chain2_it] ) {
				continue;
			}
			//check to see if the vectors are in the same direction
			//chain1 res pointed at chain 2 res?
			if ( res1_pointed_at_res2( ch1residue, ch2residue, vector_angle_cutoff, vector_dist_cutoff  ) ) {
				interface_residues[*chain1_it] = true;
				// #ifndef NDEBUG
				//      std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain1_it) << " is included because Residue: "
				//           << pose.pdb_info()->pose2pdb(*chain2_it) <<" is pointed at it."<<std::endl;
				// #endif
			}
			//chain 2 residue pointed at chain1 residue?
			if ( res1_pointed_at_res2( ch2residue, ch1residue, vector_angle_cutoff, vector_dist_cutoff ) ) {
				interface_residues[*chain2_it] = true;
				// #ifndef NDEBUG
				//      std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain2_it) << " is included because Residue: "
				//           << pose.pdb_info()->pose2pdb(*chain1_it) <<" is pointed at it."<<std::endl;
				// #endif
			}
			//space here for Cgamma defintion if wanted


		}//end itterate over chain 2
	}//end itterate over chain 1
	// //debugging
	// for(core::Size ii = 1; ii<= interface_residues.size(); ++ii)
	//  std::cout << "End finding: Residue number: " << ii << " at_interface value: " << interface_residues[ii] << std::endl;

}//end find_interface_pointing_residues_from_neighbs


/// @details looks are residue 1 and sees if any of the side chain atoms are within the cutoff distance to residue 2
bool any_atoms_within_cutoff(core::conformation::Residue & res1,
	core::conformation::Residue & res2,
	core::Real cutoff){
	using namespace core;
	bool within_cutoff(false);
	Real cutoff_squared( cutoff * cutoff );
	//only look at side chain atoms  that are <cutoff from ANY atom in the other chain
	//res1 sidechain atoms vs all res2 atoms
	//if there is a CB atom itterate through all the side chain atms in res 1
	if ( res1.type().has("CB") ) {
		//itterate over res1 sidechain atoms
		for ( conformation::Atoms::const_iterator res1atm = res1.sidechainAtoms_begin(); res1atm != res1.heavyAtoms_end(); ++res1atm ) {
			//itterate over res2 all atoms
			for ( conformation::Atoms::const_iterator res2atm = res2.atom_begin(); res2atm != res2.heavyAtoms_end(); ++res2atm ) {
				if ( (*res1atm).xyz().distance_squared( (*res2atm).xyz() )  < cutoff_squared ) {
					within_cutoff =  true;
					break;
				}
			}
			if ( within_cutoff ) {
				break;
			}
		}
	} else { //end if residue has CB atom
		//else if there is no CB atom in that residue make one up and check it against all res2 atms
		numeric::xyzVector<core::Real> pretend_CB(select_coord_for_residue( res1 ) );
		for ( conformation::Atoms::const_iterator res2atm = res2.atom_begin(); res2atm != res2.heavyAtoms_end(); ++res2atm ) {
			if (  pretend_CB.distance_squared( (*res2atm).xyz() )  < cutoff_squared ) {
				within_cutoff =  true;
				break;
			}
		}//end loop over res2
	} //end if no CB
	return within_cutoff;
} //end any_atoms_within_cutoff


/// @details find based on chains neighbors to look for vectors within using a big cutoff for CBs
InterfacePair
find_neighbors_within_CB_cutoff( core::pose::Pose const & pose,
	core::Real big_cutoff,
	core::Size chain1 , core::Size chain2 ){
	//setup chain begin and end
	core::Size ch1_begin_num = pose.conformation().chain_begin( chain1 );
	core::Size ch1_end_num = pose.conformation().chain_end( chain1 );
	core::Size ch2_begin_num = pose.conformation().chain_begin( chain2 );
	core::Size ch2_end_num = pose.conformation().chain_end( chain2);

	//setup
	std::set<core::Size> side1_within_cutoff, side2_within_cutoff;

	//use neighbor atoms from point graph
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg);
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, big_cutoff );

	// for all nodes in chain1 == for all residues in chain 1
	// all this is setup by verify_chain_setup in InterfaceDefinitionBase
	for ( Size partner1_res = ch1_begin_num; partner1_res <= ch1_end_num; ++partner1_res ) {
		for ( conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( partner1_res ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex( partner1_res ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			// get node on other edge of that node == 2nd residue index
			Size const partner2_res = edge_iter->upper_vertex();
			// if that node(residue) is in chain 2
			if ( ( partner2_res >= ch2_begin_num ) && (partner2_res <= ch2_end_num ) ) {
				side1_within_cutoff.insert( partner1_res ); // add partner1 residue
				side2_within_cutoff.insert( partner2_res ); // add partner2 residue
			} else continue;
		} // end - for all edges of node
	} // end - for all nodes in chain1

	//return the pair of interface side sets
	return std::make_pair( side1_within_cutoff, side2_within_cutoff );
} // end find_neighbors_within_CB_cutoff

/// @details find neighbors to look for vectors within using a big cutoff for CBs
InterfacePair
find_jump_partners_within_CB_cutoff( core::pose::Pose const & pose, core::Real big_cutoff, int jump_num ){
	std::set<core::Size> side1_within_cutoff, side2_within_cutoff;

	//use neighbor atoms from point graph
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg);
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, big_cutoff );

	/// create a dummy array to initialize all the members of the partner
	/// and is_interface array to false
	/// make foldtree copy
	core::kinematics::FoldTree foldtree_copy( pose.fold_tree() );
	ObjexxFCL::FArray1D_bool partners_array ( pose.size(), false );
	//set one side of jump to true
	foldtree_copy.partition_by_jump( jump_num, partners_array );

	//itterate through all residues
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		//check edges to that residue
		for ( conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( ii ).upper_edge_list_begin(),
				edge_end_iter = pg->get_vertex( ii ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			// get node on other edge of that node == 2nd residue index
			Size const pointgraph_res = edge_iter->upper_vertex();
			//see if they are in the same residue group or not
			//need to subtract 1 because FArray1D indexes from 0
			if ( partners_array[ii - 1] == partners_array[ pointgraph_res -1 ] ) continue;
			//figure out what side of jump this is on.
			if ( partners_array[ ii - 1 ] ) {
				side1_within_cutoff.insert( ii ); // add partner1 residue
				side2_within_cutoff.insert( pointgraph_res ); // add partner2 residue
			} else {
				side1_within_cutoff.insert( pointgraph_res ); // add partner1 residue
				side2_within_cutoff.insert( ii ); // add partner2 residue
			}
		}//end itterate through neighbors point graph
	}//end itterate through all residues

	//return the pair of interface side sets
	return std::make_pair( side1_within_cutoff, side2_within_cutoff );
} // end find_jump_partners_within_CB_cutoff

/// @details find the Cbeta vector(s) from one residue to another, returns the normalized vector needed
numeric::xyzVector<core::Real>
cbeta_vector( core::conformation::Residue & res){
	//std::string const atom_to_use( "CA" );
	if ( ! res.has("CA") ) {
		return numeric::xyzVector< core::Real> (0.0);
	}
	numeric::xyzVector< core::Real > CA_position( res.atom("CA").xyz() );
	numeric::xyzVector< core::Real > CB_position( select_coord_for_residue( res ) );
	//subtract CB position from CA position to get the right vector, then .normalize()
	numeric::xyzVector< core::Real > cbvector( CB_position - CA_position );
	return cbvector.normalize();
}

/// @details selects the action position for a given residue
/// Generally CB for everything but gly, and an imaginary CB position for gly.
numeric::xyzVector<core::Real>
select_coord_for_residue(core::conformation::Residue & res){
	using namespace numeric;
	using namespace core;
	//if there is a CB then use it
	if ( res.type().has("CB") ) {
		return res.atom("CB").xyz();
		//otherwise estimate where one would be.
	} else if ( res.has("CA") && res.has("C") && res.has("N") ) {
		//locations of other bb atoms
		xyzVector< core::Real > CA_xyz ( res.atom("CA").xyz() );
		xyzVector< core::Real >  C_xyz ( res.atom("C").xyz() );
		xyzVector< core::Real >  N_xyz ( res.atom("N").xyz() );
		//   //figure out where CB would be
		//   xyzVector< Real > v1( midpoint( C_xyz, N_xyz)  );
		//   xyzVector< Real > v2( 2*CA_xzy - v1 ); //reflection of v1 about CA position
		//   Real d1( magnitdue( v1 - C_xyz) );  //distance from midpoint (v1 or v2) to an atom
		//   xyzVector< Real > dir_CB( cross(N_xyz - CA_xyz, C_xyz - CA_xyz).normalize() ); //direction of CB from V2
		//   xyzVector< Real > CB_xyz( v2 + dir_cb * d1);
		//   return CB_xyz;

		//another way to try, this gets the transform from the ideal residue type frame to the existing location
		xyzVector< core::Real > halfpoint_input = midpoint( C_xyz, N_xyz );
		HTReal input_frame( N_xyz, halfpoint_input, CA_xyz );
		//now find CB position in ideal space
		//lets use alanine as the ideal here
		chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		chemical::ResidueType const & restype= rts->name_map("ALA");//define ala
		xyzVector< core::Real > idealN  =  (  restype.atom( restype.atom_index("N" ) ).ideal_xyz() );
		xyzVector< core::Real > idealCA =  (  restype.atom( restype.atom_index("CA") ).ideal_xyz() );
		xyzVector< core::Real > idealC  =  (  restype.atom( restype.atom_index("C")  ).ideal_xyz() );
		xyzVector< core::Real > idealCB  = (  restype.atom( restype.atom_index("CB") ).ideal_xyz() );
		//now use the HT to map from ideal space to the current residue position
		xyzVector< core::Real > ideal_halfpoint = 0.5 * ( idealN + idealC );
		HTReal ideal_frame( idealN, ideal_halfpoint, idealCA );
		//do not normalize because the distance is important
		xyzVector< core::Real > CB_ideal_local = ideal_frame.to_local_coordinate( idealCB  );
		xyzVector< core::Real > CB_xyz( input_frame * CB_ideal_local );
		return CB_xyz;
	} else {
		return res.xyz( res.nbr_atom() ); // neighbor atom for non protein
	}//end if no CB

} //end select_coord_for_residue

/// @details figures out if res1 and res2 are pointing at eachother
/// @details the angle is the max angle between the two residues, dist_cutoff is how far the coords are from eachother
bool res1_pointed_at_res2( core::conformation::Residue & res1,
	core::conformation::Residue & res2,
	core::Real angle_cutoff /*degrees*/,
	core::Real dist_cutoff ){
	using namespace numeric;
	bool is_pointed_at(false);
	core::Real dist_squared(dist_cutoff * dist_cutoff);
	//get the vectors for the residues in question
	xyzVector< core::Real > res1_vector ( cbeta_vector( res1) );
	//find CB to other action coordinate residue
	xyzVector< core::Real > base_position = select_coord_for_residue ( res1 );
	xyzVector< core::Real > dest_position = select_coord_for_residue ( res2 );
	//see if the base and destination are close enough to be considered
	if ( base_position.distance_squared( dest_position ) <= dist_squared ) {
		//find the vector between residues, then calculate the dot product
		xyzVector< core::Real > base_to_dest = (dest_position - base_position).normalize();
		core::Real r1_dot_r2 = dot_product (res1_vector, base_to_dest);
		core::Real costheta = std::cos( conversions::to_radians( angle_cutoff ) );
		// is projection larger than cos(theta)?
		if ( r1_dot_r2 > costheta ) {
			is_pointed_at = true;
			// #ifndef NDEBUG
			//    std::cout<< "Residue meets angle cutoff: "<< r1_dot_r2 << " cutoff( costheata=" << costheta <<") and has distance "
			//     << numeric::distance_squared( base_position, dest_position ) <<" cutoff( "<< dist_squared<<")."<<std::endl;
			// #endif
		}
	}
	return is_pointed_at;
}//end res1_pointed_at_res2

}//end namespace util
}//end namespace select
}//end namespace core

