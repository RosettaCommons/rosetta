// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/InterfaceVectorDefinitionCalculator.cc
/// @brief  Calculates the residues at an interface between two protein chains.  The calculation is done in the following manner.  First the point graph is used to find all residues within some big cutoff of residues on the other chain.  For these residues near the interface, two metrics are used to decide if they are actually possible interface residues.  The first metric is to itterate through all the side chain atoms in the residue of interest and check to see if their distance is less than the nearby atom cutoff, if so then they are an interface residue.  If a residue does not pass that check, then two vectors are drawn, a CA-CB vector and a vector from CB to a CB atom on the neighboring chain.  The dot product between these two vectors is then found and if the angle between them is less than some cutoff then they are classified as interface.
/// @author Ben Stranges (stranges@unc.edu)

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/InterfaceVectorDefinitionCalculator.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>


// Utility headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
//#include <utility/exit.hh>
#include <utility/string_util.hh>
//#include <basic/options/option.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/HomogeneousTransform.hh>
#include <cmath>
//#include <core/pose/PDBInfo.hh> //debugging only

//#include <cassert>

// option key includes

#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

//Auto Headers
//#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

//constructors
//chain number definition
InterfaceVectorDefinitionCalculator::InterfaceVectorDefinitionCalculator( core::Size const chain1_number, core::Size const chain2_number ) :
	InterfaceDefinitionCalculator(chain1_number, chain2_number),
	CB_dist_cutoff_(10.0), //distance for big CB cutoff
	nearby_atom_cutoff_(5.5), // used for finding atoms that are close to other chain
	vector_angle_cutoff_(75.0), // used for cutoff for res1 CB to res2 CB angle cutoff
	vector_dist_cutoff_(9.0) //used for check distance between CB atoms for vector calc
{}
//full constructor that takes all of the inputs
InterfaceVectorDefinitionCalculator::InterfaceVectorDefinitionCalculator( core::Size const chain1_number,
																																					core::Size const chain2_number,
																																					core::Real CB_dist_cutoff,
																																					core::Real nearby_atom_cutoff,
																																					core::Real vector_angle_cutoff,
																																					core::Real vector_dist_cutoff ) :
	InterfaceDefinitionCalculator(chain1_number, chain2_number),
	CB_dist_cutoff_(CB_dist_cutoff), //distance for big CB cutoff
	nearby_atom_cutoff_(nearby_atom_cutoff), // used for finding atoms that are close to other chain
	vector_angle_cutoff_(vector_angle_cutoff), // used for cutoff for res1 CB to res2 CB angle cutoff
	vector_dist_cutoff_(vector_dist_cutoff)  //used for check distance between CB atoms for vector calc
{}
//chain character definition
InterfaceVectorDefinitionCalculator::InterfaceVectorDefinitionCalculator( char const chain1_letter,
																																					char const chain2_letter ) :
	InterfaceDefinitionCalculator(chain1_letter, chain2_letter),
	CB_dist_cutoff_(10.0), //distance for big CB cutoff
	nearby_atom_cutoff_(5.5), // used for finding atoms that are close to other chain
	vector_angle_cutoff_(75.0), // used for cutoff for res1 CB to res2 CB angle cutoff
	vector_dist_cutoff_(9.0)
{}
//full constructor that takes all of the inputs and uses chain characters
InterfaceVectorDefinitionCalculator::InterfaceVectorDefinitionCalculator( char const chain1_letter,
																																					char const chain2_letter,
																																					core::Real CB_dist_cutoff,
																																					core::Real nearby_atom_cutoff,
																																					core::Real vector_angle_cutoff,
																																					core::Real vector_dist_cutoff ) :
	InterfaceDefinitionCalculator(chain1_letter, chain2_letter),
	CB_dist_cutoff_(CB_dist_cutoff), //distance for big CB cutoff
	nearby_atom_cutoff_(nearby_atom_cutoff), // used for finding atoms that are close to other chain
	vector_angle_cutoff_(vector_angle_cutoff), // used for cutoff for res1 CB to res2 CB angle cutoff
	vector_dist_cutoff_(vector_dist_cutoff)
{}
// clone needs to return *this because clone is used to call this functionality
core::pose::metrics::PoseMetricCalculatorOP InterfaceVectorDefinitionCalculator::clone() const
{ return new InterfaceVectorDefinitionCalculator( *this ); }

// pose metric lookups
void InterfaceVectorDefinitionCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "interface_residues" ) {
		basic::check_cast( valptr, &interface_residues_, "interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( interface_residues_ );

	} else if ( key == "first_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain1_interface_residues_, "first_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain1_interface_residues_ );

	} else if ( key == "second_chain_interface_residues" ) {
		basic::check_cast( valptr, &chain2_interface_residues_, "second_chain_interface_residues expects to return a std::set< Size >" );
		(static_cast<basic::MetricValue<std::set<Size> > *>(valptr))->set( chain2_interface_residues_ );

	} else if ( key == "num_interface_residues" ) {
		basic::check_cast( valptr, &num_interface_residues_, "num_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_interface_residues_ );

	} else if ( key == "num_first_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain1_interface_residues_, "num_first_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain1_interface_residues_ );

	} else if ( key == "num_second_chain_interface_residues" ) {
		basic::check_cast( valptr, &num_chain2_interface_residues_, "num_second_chain_interface_residues expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( num_chain2_interface_residues_ );

	} else if ( key == "first_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch1_begin_num_, "first_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_begin_num_ );

	} else if ( key == "first_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch1_end_num_, "first_chain_last_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch1_end_num_ );

	} else if ( key == "second_chain_first_resnum" ) {
		basic::check_cast( valptr, &ch2_begin_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_begin_num_ );

	} else if ( key == "second_chain_last_resnum" ) {
		basic::check_cast( valptr, &ch2_end_num_, "second_chain_first_resnum expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( ch2_end_num_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string InterfaceVectorDefinitionCalculator::print( std::string const & key ) const {

	if ( key == "interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "first_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "second_chain_interface_residues" ) {
		basic::Error() << "No output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "num_interface_residues" ) {
		return utility::to_string( num_interface_residues_ );
	} else if ( key == "num_first_chain_interface_residues" ) {
		return utility::to_string( num_chain1_interface_residues_ );
	} else if ( key == "num_second_chain_interface_residues" ) {
		return utility::to_string( num_chain2_interface_residues_ );
	} else if ( key == "first_chain_first_resnum" ) {
		return utility::to_string( ch1_begin_num_ );
	} else if ( key == "first_chain_last_resnum" ) {
		return utility::to_string( ch1_end_num_ );
	} else if ( key == "second_chain_first_resnum" ) {
		return utility::to_string( ch2_begin_num_ );
	} else if ( key == "second_chain_last_resnum" ) {
		return utility::to_string( ch2_end_num_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}

//PRIVATE functions to do the work
///@details main function that goes through all of the important calculations
/// Takes input of possible interface residues and prunes them down using the nearby atom cutoff and the vector cutoff
void InterfaceVectorDefinitionCalculator::find_interface_pointing_residues_from_neighbs(core::pose::Pose pose, InterfacePair interface_pairs ){
	using namespace utility;
	using namespace core;
	//itterate over pairs
	for( std::set<Size>::const_iterator chain1_it = interface_pairs.first.begin();  chain1_it!=interface_pairs.first.end(); ++chain1_it ){
		conformation::Residue ch1residue(pose.residue(*chain1_it));
		for(  std::set<Size>::const_iterator chain2_it = interface_pairs.second.begin();  chain2_it!=interface_pairs.second.end(); ++chain2_it ){
			conformation::Residue ch2residue(pose.residue(*chain2_it));
			//check to see if both are already in the interface set
			if( chain1_interface_residues_.count(*chain1_it) != 0 && chain2_interface_residues_.count(*chain2_it) != 0 )
				continue;
			//Real cutoff(5.0);//debugging
			//use a stricter distance cutoff to define any very close to other chain
			bool r1_near_r2 (any_atoms_within_cutoff( ch1residue, ch2residue, nearby_atom_cutoff_ ) );
			bool r2_near_r1 (any_atoms_within_cutoff( ch2residue, ch1residue, nearby_atom_cutoff_ ) );
			//are residue1 sidechain atoms near residue2 atoms?
			if( r1_near_r2 ){
				chain1_interface_residues_.insert(*chain1_it);
				interface_residues_.insert(*chain1_it);
// #ifndef NDEBUG
// 				std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain1_it) << " is included because Residue: "
// 								 <<	pose.pdb_info()->pose2pdb(*chain2_it) <<" meets the nearby cutoff."<<std::endl;
// #endif
			}
			//are residue2 sidechain atoms near residue1 atoms?
			if( r2_near_r1 ){
				chain2_interface_residues_.insert(*chain2_it);
				interface_residues_.insert(*chain2_it);
// #ifndef NDEBUG
// 				std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain2_it) << " is included because Residue: "
// 								 <<	pose.pdb_info()->pose2pdb(*chain1_it) <<" meets the nearby cutoff."<<std::endl;
// #endif
// 				continue;
			}
			//see if both are in the interface set NOW
			if( chain1_interface_residues_.count(*chain1_it) != 0 && chain2_interface_residues_.count(*chain2_it) != 0 )
				continue;
			//chain1 res pointed at chain 2 res?
			if( res1_pointed_at_res2( ch1residue, ch2residue, vector_angle_cutoff_, vector_dist_cutoff_  ) ){
				chain1_interface_residues_.insert(*chain1_it);
				interface_residues_.insert(*chain1_it);
// #ifndef NDEBUG
// 				std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain1_it) << " is included because Residue: "
// 								 <<	pose.pdb_info()->pose2pdb(*chain2_it) <<" is pointed at it."<<std::endl;
// #endif
			}
			//chain 2 residue pointed at chain1 residue?
			if( res1_pointed_at_res2( ch2residue, ch1residue, vector_angle_cutoff_, vector_dist_cutoff_ ) ){
				chain2_interface_residues_.insert(*chain2_it);
				interface_residues_.insert(*chain2_it);
// #ifndef NDEBUG
// 				std::cout<< "Residue:  "<<pose.pdb_info()->pose2pdb(*chain2_it) << " is included because Residue: "
// 								 <<	pose.pdb_info()->pose2pdb(*chain1_it) <<" is pointed at it."<<std::endl;
// #endif
				continue;
			}
			//space here for Cgamma defintion

		} //end residues in chain 2
	} //end residues in chain 1
} //end  find_interface_pointing_residues_from_neighbs


///@details  looks are residue 1 and sees if any of the side chain atoms are within the cutoff distance to residue 2
bool InterfaceVectorDefinitionCalculator::any_atoms_within_cutoff(core::conformation::Residue & res1,
																																	core::conformation::Residue & res2,
																																	core::Real & cutoff){
	using namespace core;
	bool within_cutoff(false);
	Real cutoff_squared( cutoff * cutoff );
	//only look at side chain atoms  that are <cutoff from ANY atom in the other chain
	//res1 sidechain atoms vs all res2 atoms
	//if there is a CB atom itterate through all the side chain atms in res 1
	if ( res1.type().has("CB") ){
		//itterate over res1 sidechain atoms
		for( conformation::Atoms::const_iterator res1atm = res1.sidechainAtoms_begin(); res1atm <= res1.heavyAtoms_end(); ++res1atm ){
			//itterate over res2 all atoms
			for( conformation::Atoms::const_iterator res2atm = res2.atom_begin(); res2atm != res2.heavyAtoms_end(); ++res2atm ){
				if ( (*res1atm).xyz().distance_squared( (*res2atm).xyz() )  < cutoff_squared){
					within_cutoff =  true;
					break;
				}
			}
			if(within_cutoff)
				break;
		}
	} //end if residue has CB atom
	//else if there is no CB atom in that residue make one up and check it against all res2 atms
	else{
		numeric::xyzVector<core::Real> pretend_CB(select_coord_for_residue( res1 ) );
		for( conformation::Atoms::const_iterator res2atm = res2.atom_begin(); res2atm != res2.heavyAtoms_end(); ++res2atm ){
			if (  pretend_CB.distance_squared( (*res2atm).xyz() )  < cutoff_squared){
				within_cutoff =  true;
				break;
			}
		}//end loop over res2
	} //end if no CB
	return within_cutoff;
}


///@details find neighbors to look for vectors within using a big cutoff for CBs
InterfaceVectorDefinitionCalculator::InterfacePair
InterfaceVectorDefinitionCalculator::find_neighbors_within_CB_cutoff( core::pose::Pose pose, core::Real big_cutoff ){
	std::set<core::Size> side1_within_cutoff, side2_within_cutoff;

 	//use neighbor atoms from point graph
	conformation::PointGraphOP pg( new conformation::PointGraph );
	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg);
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, big_cutoff );

	// for all nodes in chain1 == for all residues in chain 1
	// all this is setup by verify_chain_setup in InterfaceDefinitionBase
	utility::vector1< Size > chain1_interface, chain2_interface;
	for ( Size partner1_res = ch1_begin_num_; partner1_res <= ch1_end_num_; ++partner1_res ) {
		for ( conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( partner1_res ).upper_edge_list_begin(),
						edge_end_iter = pg->get_vertex( partner1_res ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {
			// get node on other edge of that node == 2nd residue index
			Size const partner2_res = edge_iter->upper_vertex();
			// if that node(residue) is in chain 2
			if ( ( partner2_res >= ch2_begin_num_ ) && (partner2_res <= ch2_end_num_ ) ) {
				side1_within_cutoff.insert( partner1_res );	// add partner1 residue
				side2_within_cutoff.insert( partner2_res );	// add partner2 residue
			}
			else continue;
		} // end - for all edges of node
	}	// end - for all nodes in chain1

	//return the pair of interface side sets
	return std::make_pair( side1_within_cutoff, side2_within_cutoff );
} // end find_neighbors_within_CB_cutoff


///@details find the Cbeta vector(s) from one residue to another, returns the normalized vector needed
numeric::xyzVector<core::Real>
InterfaceVectorDefinitionCalculator::cbeta_vector( core::conformation::Residue & res){
	//std::string const atom_to_use( "CA" );
	numeric::xyzVector< core::Real > CA_position( res.atom("CA").xyz() );
	numeric::xyzVector< core::Real > CB_position( select_coord_for_residue( res ) );
	//subtract CB position from CA position to get the right vector, then .normalize()
	numeric::xyzVector< core::Real > cbvector( CB_position - CA_position );
	return cbvector.normalize();
}


///@details selects the action position for a given residue
/// Generally CB for everything but gly, and an imaginary CB position for gly.
numeric::xyzVector<core::Real>
InterfaceVectorDefinitionCalculator::select_coord_for_residue(core::conformation::Residue & res){
	using namespace numeric;
	using namespace core;
	//if there is a CB then use it
	if ( res.type().has("CB") )
		return res.atom("CB").xyz();
	//otherwise estimate where one would be.
	else{
 		//locations of other bb atoms
 		xyzVector< core::Real > CA_xyz ( res.atom("CA").xyz() );
 		xyzVector< core::Real >  C_xyz ( res.atom("C").xyz() );
 		xyzVector< core::Real >  N_xyz ( res.atom("N").xyz() );
// 		//figure out where CB would be
// 		xyzVector< Real > v1( midpoint( C_xyz, N_xyz)  );
// 		xyzVector< Real > v2( 2*CA_xzy - v1 ); //reflection of v1 about CA position
// 		Real d1( magnitdue( v1 - C_xyz) );  //distance from midpoint (v1 or v2) to an atom
// 		xyzVector< Real > dir_CB( cross(N_xyz - CA_xyz, C_xyz - CA_xyz).normalize() ); //direction of CB from V2
// 		xyzVector< Real > CB_xyz( v2 + dir_cb * d1);
// 		return CB_xyz;

		//another way to try, this gets the transform from the ideal residue type frame to the existing location
		xyzVector< core::Real > halfpoint_input = midpoint( C_xyz, N_xyz );
		HTReal input_frame( N_xyz, halfpoint_input, CA_xyz );
		//now find CB position in ideal space
		//lets use alanine as the ideal here
		chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		chemical::ResidueType restype= rts->name_map("ALA");//define ala
		xyzVector< core::Real > idealN  =  (  restype.xyz( restype.atom_index("N" ) ) );
		xyzVector< core::Real > idealCA =  (  restype.xyz( restype.atom_index("CA") ) );
		xyzVector< core::Real > idealC  =  (  restype.xyz( restype.atom_index("C")  ) );
		xyzVector< core::Real > idealCB  = (  restype.xyz( restype.atom_index("CB") ) );
		//now use the HT to map from ideal space to the current residue position
		xyzVector< core::Real > ideal_halfpoint = 0.5 * ( idealN + idealC );
		HTReal ideal_frame( idealN, ideal_halfpoint, idealCA );
		//do not normalize because the distance is important
		xyzVector< core::Real > CB_ideal_local = ideal_frame.to_local_coordinate( idealCB  );
		xyzVector< core::Real > CB_xyz( input_frame * CB_ideal_local );
		return CB_xyz;
	} //end if no CB

} //end select_coord_for_residue

///@details figures out if res1 and res2 are pointing at eachother
///@details the angle is the max angle between the two residues, dist_cutoff is how far the coords are from eachother
bool InterfaceVectorDefinitionCalculator::res1_pointed_at_res2( core::conformation::Residue & res1,
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
	if ( base_position.distance_squared( dest_position ) <= dist_squared){
		//find the vector between residues, then calculate the dot product
		xyzVector< core::Real > base_to_dest = (dest_position - base_position).normalize();
		core::Real r1_dot_r2 = dot_product (res1_vector, base_to_dest);
		core::Real costheta = std::cos( conversions::to_radians( angle_cutoff ) );
		// is projection larger than cos(theta)?
		if (r1_dot_r2 > costheta){
			is_pointed_at = true;
// #ifndef NDEBUG
// 			std::cout<< "Residue meets angle cutoff: "<< r1_dot_r2 << " cutoff( costheata=" << costheta <<") and has distance "
// 				<< numeric::distance_squared( base_position, dest_position ) <<" cutoff( "<< dist_squared<<")."<<std::endl;
// #endif
		}
	}
	return	is_pointed_at;

}//end res1_pointed_at_res2


void InterfaceVectorDefinitionCalculator::recompute( Pose const & pose ) {
	//clear all the old stuff
	interface_residues_.clear();
	chain1_interface_residues_.clear();
	chain2_interface_residues_.clear();
	//Base class setup
	verify_chain_setup( pose );

	//first find all the neighbors within some Cbeta cutoff distance from eachother
	InterfacePair CB_pairs_list = find_neighbors_within_CB_cutoff(pose,  CB_dist_cutoff_ );
	//second, prune this set down to what matters
	find_interface_pointing_residues_from_neighbs( pose, CB_pairs_list );

	num_interface_residues_ = interface_residues_.size();
	num_chain1_interface_residues_ = chain1_interface_residues_.size();
	num_chain2_interface_residues_ = chain2_interface_residues_.size();

}

//set from options function
void set_from_options(){
	//use setters
}


} // PoseMetricCalculators
} // toolbox
} // protocols
