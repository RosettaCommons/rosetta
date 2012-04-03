// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNode.cc
/// @brief  .cc file for inverse rotamer tree node
/// @author Florian Richter, flosopher@gmail.com, mar 2012


//unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

//project headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>

#include <basic/Tracer.hh>

//#include <core/io/pdb/pose_io.hh> //debug only include
//#include <fstream> //debug only include
//#include <utility/string_util.hh> //debug only include

//c++ headers

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static basic::Tracer tr( "protocols.toolbox.match_enzdes_util.InvrotTreeNode" );

InvrotTreeNode::InvrotTreeNode( InvrotTreeNodeBaseCAP parent )
	: InvrotTreeNodeBase( parent ),
		geom_cst_(0) //location_in_parent_node_(0)
{
	invrots_and_next_nodes_.clear();
}

InvrotTreeNode::~InvrotTreeNode() {}


/// @details
/// assumptions: this node represents the upstream_res
/// of this geom_cst (invrot_geomcst), and the given
/// target_residue corresponds to the downstream_res of
/// invrot_geomcst
bool
InvrotTreeNode::initialize_from_enzcst_io(
	core::conformation::Residue const & target_residue,
	EnzConstraintIOCOP enzcst_io,
	Size invrot_geomcst,
	Size target_geomcst
)
{
	invrots_and_next_nodes_.clear(); //safety
	geom_cst_ = invrot_geomcst;

	//1. create all invrots against target residue
	//mcfi(1) is hardcoded because matcher can't
	//deal with upstream-upstream matching when the
	//target residue has explicit ambiguity
	//we have to put in an early catch for this somewhere..
	std::list< core::conformation::ResidueCOP > all_invrots( enzcst_io->mcfi_list( geom_cst_)->inverse_rotamers_against_residue( enzcst_io->mcfi_list( geom_cst_ )->mcfi( 1 )->downstream_res(), &target_residue ) );

	//2. do clash check against all parent residues
	//if all get removed, this means this is a dead end
	//and we return false
	Size num_rots_total( all_invrots.size() );
	this->remove_invrots_clashing_with_parent_res( all_invrots );
	if( all_invrots.size() == 0 ){
		return false;
	}
	Size num_invrots_clashing( num_rots_total - all_invrots.size() );
	tr << "When initializing a node for geomcst " << geom_cst_ << ", " << num_invrots_clashing << " of a total of " << num_rots_total << " inverse rotamers were found to clash with something." << std::endl;

	//3. now we have to figure out if the upstream_res of
	//this geom_cst (i.e. what this node represents) is the
	//downstream_res in a later mcfi
	//i.e. do we need to build inverse rotamers against
	//the inverse rotamers of this node?
	utility::vector1< Size > dependent_mcfi;
	for( Size i = geom_cst_ +1; i <= enzcst_io->num_mcfi_lists(); ++i){
		std::pair< Size, Size> const & target_res( enzcst_io->target_downstream_res()[i] );
		if( (target_res.first == geom_cst_) && (target_res.second == 2 ) ) dependent_mcfi.push_back( i );
	}

	//if there are no dependent mcfi, we're done
	//simply put the inverse rots and an empty vector
	//into the invrots_and_next_nodes_ vector
	if( dependent_mcfi.size() == 0 ){
		invrots_and_next_nodes_.push_back( std::pair< std::list<core::conformation::ResidueCOP >, utility::vector1< InvrotTreeNodeOP > > () );
		invrots_and_next_nodes_[1].first = all_invrots;
		return true;
	}

	//4. if we get here, that means we need to build
	//inverse_rotamer_trees against the inverse rotamers
	//in this node
	//it's gettin' brutishly complicated, because we have to figure
	//out the redundancy groups of this nodes invrots', i.e. against
	//how many of the invrots here do we need to build separate
	//dependent inverse rotamer trees...
	//ok, let's go
	//4.1 first invrot is always non redundant
	invrots_and_next_nodes_.push_back( std::pair< std::list<core::conformation::ResidueCOP >, utility::vector1< InvrotTreeNodeOP > > () );
	invrots_and_next_nodes_[1].first.push_back( *(all_invrots.begin()) );

	//4.2 now loop through remaining invrots and
	//check for redundancy to other invrots with respect to
	//positioning of dependent inverse rotamer trees
	std::list<core::conformation::ResidueCOP>::iterator invrot_it1 = all_invrots.begin();
	invrot_it1++;
	while( invrot_it1 != all_invrots.end() ){
		bool this_invrot_redundant(true);
		for( core::Size i = 1; i <= invrots_and_next_nodes_.size(); ++i ){

			//4.3 really complicated code goes here
			//not written yet, for now we just push back into first list, i.e. everything redundant
			invrots_and_next_nodes_[1].first.push_back( *invrot_it1 );

		}//loop over all so-far non redundant rots

		//if this rotamer is non redundant, we throw it into a new list
		if( !this_invrot_redundant ){
			invrots_and_next_nodes_.push_back( std::pair< std::list<core::conformation::ResidueCOP >, utility::vector1< InvrotTreeNodeOP > > () );
			invrots_and_next_nodes_[invrots_and_next_nodes_.size() ].first.push_back( *invrot_it1 );
		}
		invrot_it1++;
	} //while loop over invrot_it1

	//4.4 redundancy determined, almost done, now we need to
	//initialize the daughter nodes for each set of non-redundant
	//invrots. note that if any of these nodes fail to initialize,
	//that renders the corresponding set of invrots a dead end
	for( utility::vector1< invrots_node_ptrs_pair >::iterator pair_it( invrots_and_next_nodes_.begin() ); pair_it != invrots_and_next_nodes_.end(); /*increment happening in loop*/ ){

		core::conformation::Residue const & this_target( **(pair_it->first.begin()) );
		bool all_initialization_successful(true);
		for( Size j = 1; j <= dependent_mcfi.size(); ++j ){
			InvrotTreeNodeOP child = new InvrotTreeNode( this );
			pair_it->second.push_back( child );
			if( ! child->initialize_from_enzcst_io( this_target, enzcst_io, dependent_mcfi[j],geom_cst_ ) ){
				//this_node_ptr_pair_dead_end = true;
				pair_it = invrots_and_next_nodes_.erase( pair_it ); //note: erasing from vector, not ideal, but the vectors should usually be fairly small and the initialization shenanigans are only called once
				all_initialization_successful = false;
				break;
			}
		} // loop over dependent mcfi
		if( all_initialization_successful ) ++pair_it;
	}//loop over all non-redundant invrot groups
	//i guess we have to set the locations in this node in the parent nodes
	for( Size i =1; i <= invrots_and_next_nodes_.size(); ++i ){
		for( utility::vector1< InvrotTreeNodeBaseOP >::iterator child_it( invrots_and_next_nodes_[i].second.begin() ); child_it != invrots_and_next_nodes_[i].second.end(); ++child_it ){
			(*child_it)->set_location_in_parent_node( i );
		}
	}

	//4.5 in case there were dead ends only, we return false
	if( invrots_and_next_nodes_.size() == 0 ) return false;

	return true;
}

/// @details the real meat of this thing
core::scoring::constraints::ConstraintCOP
InvrotTreeNode::generate_constraints(
	core::pose::Pose const &,  //pose,  commented out to prevent compiler warning
	AllowedSeqposForGeomCstOP //geomcst_seqpos commented out to prevent compiler warning
) const
{
	//bunch of code, not written yet
	for( core::Size i = 1; i <= invrots_and_next_nodes_.size(); ++i ) {}
	return NULL;  // required for compilation on Windows
}


utility::vector1< std::list< core::conformation::ResidueCOP > >
InvrotTreeNode::all_target_residues( InvrotTreeNodeBaseCAP child_node ) const
{

	utility::vector1< std::list< core::conformation::ResidueCOP > > to_return;

	//1. get the target residues of the parent node
	InvrotTreeNodeBaseCAP parent(this->parent_node() );
	if( parent ) to_return = parent->all_target_residues( this );

	//2. add the target residues from this node corresponding
	//corresponding to the asking child node
	bool child_found(false);
	for( Size i =1; i <= invrots_and_next_nodes_.size(); ++i ){
		for( Size j =1; j <= invrots_and_next_nodes_[i].second.size(); ++j ){
			if( child_node == invrots_and_next_nodes_[i].second[j] ){
				to_return.push_back( invrots_and_next_nodes_[i].first );
				child_found = true;
				break;
			}
		}
		if( child_found ) break;
	}

	return to_return;
}


void
InvrotTreeNode::remove_invrots_clashing_with_parent_res(
	std::list< core::conformation::ResidueCOP > & invrots
) const
{

	InvrotTreeNodeBaseCAP parent(this->parent_node() );

	//in case there's no parent, i.e. in the unit test
	//or potential use of this outside the tree, there
	//are no target to clash with, so we'll return
	if( !parent ) return;

	//1. get all parent rotamers
	utility::vector1< std::list< core::conformation::ResidueCOP > > all_targets ( (*parent).all_target_residues( this ) );

	//safety: if for some reason all_targets is empty,
	//we return
	if( all_targets.size() == 0 ) return;

	//2.clash check between invrots and parent residues
	//all invrots that clash with all rotamers from any
	//of the all_targets list will be removed
	for( std::list <core::conformation::ResidueCOP >::iterator invrot_it( invrots.begin() ); invrot_it != invrots.end(); /*increment in loop*/ ){
		bool cur_invrot_clashes(false);

		for( Size i =1; i <= all_targets.size(); ++i){

			//safety: if for some reason this list is empty, we continue
			if( all_targets[i].size() == 0 ) continue;

			bool this_list_all_clash( true );
			for(  std::list <core::conformation::ResidueCOP >::const_iterator target_rot_it( all_targets[i].begin() ), target_rot_it_end( all_targets[i].end() ); target_rot_it != target_rot_it_end; ++target_rot_it ){

				//the actual clash check
				//hardcoded as hadrsphere here for now, but in the future this could be
				//combined with matcher filters
				bool these_two_clash(false);
				core::Real cutoff_sq = 2.5 * 2.5;
				core::conformation::Residue const & res1(**invrot_it), res2(**target_rot_it);
				for( Size res1at = 1; res1at <= res1.nheavyatoms(); ++res1at ){
					if( res1.type().is_virtual( res1at ) ) continue;
					for( Size res2at = 1; res2at <= res2.nheavyatoms(); ++res2at){
						if( res2.type().is_virtual( res2at ) ) continue;
						if( res1.atom( res1at ).xyz().distance_squared( res2.atom( res2at ).xyz() ) < cutoff_sq ){
							these_two_clash = true;

							/* //debug
							static Size clashrescount(0);
							clashrescount++;
							tr << "ACKACK clash num " << clashrescount << " between res1 atom " << res1.type().atom_name( res1at ) << " and res2 atom " << res2.type().atom_name( res2at ) << " which are sqrt(" << res1.atom( res1at ).xyz().distance_squared( res2.atom( res2at ).xyz() ) << " apart." << std::endl;
							std::ofstream file_out( ("clashno_"+utility::to_string( clashrescount )+".pdb").c_str() );
							Size atomcounter(0);
							core::io::pdb::dump_pdb_residue( res1, atomcounter, file_out );
							core::io::pdb::dump_pdb_residue( res2, atomcounter, file_out );
							file_out.close();
							//debug over */

							break; //breaks loop over res2 atoms
						}
					} //loop over res2 atoms
					if( these_two_clash ) break; //breaks loop over res 1 atoms
				}//loop over res1 atoms
				if( !these_two_clash){
					this_list_all_clash = false;
					break; //breaks loop over rots in this list
				}
			} //for loop over target_rot
			if( this_list_all_clash ){
				cur_invrot_clashes = true;
				break; // breaks loop over target_lists
			}
		} //for loop over target_rot_lists

		if( cur_invrot_clashes ) invrot_it = invrots.erase( invrot_it );
		else ++invrot_it;
	} //loop over all invrots
}


void
InvrotTreeNode::collect_all_inverse_rotamers(
	utility::vector1< InvrotCollectorOP > & invrot_collectors
) const
{

	runtime_assert( invrots_and_next_nodes_.size() != 0 );
	//1. find all invrot collectors that have stuff from this node's parent
	//and location in them
	utility::vector1< Size > empty_parent_collectors;
	utility::vector1< Size > filled_parent_collectors;

	for( Size i =1; i <= invrot_collectors.size(); ++i ){
		std::map< InvrotTreeNodeBaseCOP, Size > const & collector_map( invrot_collectors[i]->owner_nodes_and_locations() );
		//safety
		if( collector_map.find( this ) != collector_map.end() ) utility_exit_with_message("InvrotTreeNode being asked to fill a collector that it apparently already previously filled.");

		std::map< InvrotTreeNodeBaseCOP, Size >::const_iterator map_it( collector_map.find( &(*(this->parent_node())) ) );

		if( map_it != collector_map.end() ){
			if( map_it->second == this->location_in_parent_node() ){
				if( invrot_collectors[i]->invrots()[ geom_cst_].size() == 0 ) empty_parent_collectors.push_back( i );
				else filled_parent_collectors.push_back( i );
			}
		}
	}
	tr << "Collecting inverse rotamers for a node from geomcst " << geom_cst_ << ". " << empty_parent_collectors.size() << " empty parent collectors and " << filled_parent_collectors.size() << " filled parent collectors were found." << std::endl;

	//1.b make sure we found at least one parent
	Size total_parents( empty_parent_collectors.size() + filled_parent_collectors.size() );
	runtime_assert( total_parents != 0);

	//2 for every parent, we need to add rotamers for every definiton
	utility::vector1< Size > collectors_to_fill( empty_parent_collectors );
	for( Size i =1; i <= filled_parent_collectors.size(); ++ i ){
		invrot_collectors.push_back( invrot_collectors[ filled_parent_collectors[i] ]->clone() );
		collectors_to_fill.push_back( invrot_collectors.size() );
	}

	//3. fill 'er up
	for( Size i = 1; i <= collectors_to_fill.size(); ++i ){
		Size overflow_start( invrot_collectors.size());
		for( Size j = 2; j <= invrots_and_next_nodes_.size(); ++j ){
			invrot_collectors.push_back( invrot_collectors[ collectors_to_fill[i] ]->clone() );
		}

		invrot_collectors[ collectors_to_fill[ i ] ]->set_invrots_for_listnum( geom_cst_, invrots_and_next_nodes_[1].first, this, 1 );

		for( Size j =2; j <= invrots_and_next_nodes_.size(); ++j ){
			invrot_collectors[ overflow_start + j - 1 ]->set_invrots_for_listnum( geom_cst_, invrots_and_next_nodes_[j].first, this, j );
		}
	}

	//4. call this shit on all daughter nodes
	for( Size i = 1; i <= invrots_and_next_nodes_.size(); ++i){
		for( utility::vector1< InvrotTreeNodeBaseOP >::const_iterator child_it( invrots_and_next_nodes_[i].second.begin() ); child_it != invrots_and_next_nodes_[i].second.end(); ++child_it ){
			(*child_it)->collect_all_inverse_rotamers( invrot_collectors );
		}
	}
}

}
}
}
