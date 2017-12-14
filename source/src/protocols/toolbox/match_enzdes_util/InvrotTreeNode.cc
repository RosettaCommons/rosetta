// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNode.cc
/// @brief  .cc file for inverse rotamer tree node
/// @author Florian Richter, flosopher@gmail.com, mar 2012


//unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>

//project headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <basic/Tracer.hh>

//  //debug only include
//#include <fstream> //debug only include
//#include <utility/string_util.hh> //debug only include

//c++ headers

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static basic::Tracer tr( "protocols.toolbox.match_enzdes_util.InvrotTreeNode" );

InvrotTreeNode::InvrotTreeNode( InvrotTreeNodeBaseCAP parent )
: InvrotTreeNodeBase( parent ),
	geom_cst_(0), generate_invrot_csts_(true)
{
	invrots_and_next_nodes_.clear();
}

InvrotTreeNode::~InvrotTreeNode() = default;


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
	core::pose::PoseCOP pose
)
{
	invrots_and_next_nodes_.clear(); //safety
	std::list< core::conformation::ResidueCOP > all_invrots;
	all_invrots.clear();
	if ( pose ) {
		core::conformation::ResidueCOP cst_res( cst_residue_in_pose( *pose, invrot_geomcst,  enzcst_io->mcfi_list( invrot_geomcst )->mcfi( 1 )->upstream_res() ) );
		if ( cst_res ) all_invrots.push_back( cst_res );
	}

	//1. create all invrots against target residue
	//mcfi(1) is hardcoded because matcher can't
	//deal with upstream-upstream matching when the
	//target residue has explicit ambiguity
	//we have to put in an early catch for this somewhere..
	if ( all_invrots.empty() ) { //size() == 0 ){
		all_invrots =  enzcst_io->mcfi_list( invrot_geomcst)->inverse_rotamers_against_residue( enzcst_io->mcfi_list( invrot_geomcst )->mcfi( 1 )->downstream_res(), target_residue.get_self_ptr() );

		//2. do clash check against all parent residues
		//if all get removed, this means this is a dead end
		//and we return false
		Size num_rots_total( all_invrots.size() );
		this->remove_invrots_clashing_with_parent_res( all_invrots, enzcst_io->mcfi_list(invrot_geomcst)->mcfi(1)->is_covalent() );
		Size num_invrots_clashing( num_rots_total - all_invrots.size() );
		tr << "When initializing a node for geomcst " <<invrot_geomcst << ", " << num_invrots_clashing << " of a total of " << num_rots_total << " inverse rotamers were found to clash with something." << std::endl;
		if ( all_invrots.empty() ) { //size() == 0 ){
			return false;
		}
	}

	return this->initialize_from_enzcst_io_and_invrots( all_invrots, enzcst_io, invrot_geomcst, pose );
}

bool
InvrotTreeNode::initialize_from_enzcst_io_and_invrots(
	std::list< core::conformation::ResidueCOP > const & all_invrots,
	EnzConstraintIOCOP enzcst_io,
	Size invrot_geomcst,
	core::pose::PoseCOP pose
)
{

	geom_cst_ = invrot_geomcst;
	//3. now we have to figure out if the upstream_res of
	//this geom_cst (i.e. what this node represents) is the
	//downstream_res in a later mcfi
	//i.e. do we need to build inverse rotamers against
	//the inverse rotamers of this node?
	utility::vector1< Size > dependent_mcfi;
	for ( Size i = geom_cst_ +1; i <= enzcst_io->num_mcfi_lists(); ++i ) {
		std::pair< Size, Size> const & target_res( enzcst_io->target_downstream_res()[i] );
		if ( (target_res.first == geom_cst_) && (target_res.second == 2 ) ) dependent_mcfi.push_back( i );
	}

	//if there are no dependent mcfi, we're done
	//simply put the inverse rots and an empty vector
	//into the invrots_and_next_nodes_ vector
	if ( dependent_mcfi.size() == 0 ) {
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
	auto invrot_it1 = all_invrots.begin();
	++invrot_it1; //first one was already taken care of
	while ( invrot_it1 != all_invrots.end() ) {
		bool this_invrot_redundant(false);
		for ( core::Size i = 1; i <= invrots_and_next_nodes_.size(); ++i ) {

			//4.3 really complicated code goes here
			//let's see: we have to see whether the two rotamers
			//are redundant according to each of the dependent mcfi,
			//i.e. we have to ask the downstream EnzCstTemplateRes
			//of the dependent mcfi whether the two rotamers are redundant
			core::conformation::Residue const & cur_rot(**invrot_it1), test_rot( **(invrots_and_next_nodes_[i].first.begin()) );
			bool these_two_redundant(true);
			for ( Size j = 1; j <= dependent_mcfi.size(); ++j ) {

				if ( !enzcst_io->mcfi_list( dependent_mcfi[j])->mcfi(1)->enz_cst_template_res( enzcst_io->mcfi_list( dependent_mcfi[j])->mcfi(1)->downstream_res() )->residue_conformations_redundant( cur_rot, test_rot ) ) {
					these_two_redundant = false;
					break;
				}
			}
			if ( these_two_redundant ) {
				invrots_and_next_nodes_[i].first.push_back( *invrot_it1 );
				this_invrot_redundant = true;
				break;
			}
		}//loop over all so-far non redundant rots

		//if this rotamer is non redundant, we throw it into a new list
		if ( !this_invrot_redundant ) {
			invrots_and_next_nodes_.push_back( std::pair< std::list<core::conformation::ResidueCOP >, utility::vector1< InvrotTreeNodeOP > > () );
			invrots_and_next_nodes_[invrots_and_next_nodes_.size() ].first.push_back( *invrot_it1 );
		}
		++invrot_it1;
	} //while loop over invrot_it1
	tr << "Node initialization for geomcst " << geom_cst_ << ". After redundancy determination, " << invrots_and_next_nodes_.size() << " non-redundant sets of inverse rotamers exist." << std::endl;
	//4.4 redundancy determined, almost done, now we need to
	//initialize the daughter nodes for each set of non-redundant
	//invrots. note that if any of these nodes fail to initialize,
	//that renders the corresponding set of invrots a dead end
	for ( auto pair_it( invrots_and_next_nodes_.begin() ); pair_it != invrots_and_next_nodes_.end(); /*increment happening in loop*/ ) {

		core::conformation::Residue const & this_target( **(pair_it->first.begin()) );
		bool all_initialization_successful(true);
		for ( Size j = 1; j <= dependent_mcfi.size(); ++j ) {
			InvrotTreeNodeOP child( new InvrotTreeNode( get_self_weak_ptr() ) );
			pair_it->second.push_back( child );
			if ( ! child->initialize_from_enzcst_io( this_target, enzcst_io, dependent_mcfi[j], pose ) ) {
				pair_it = invrots_and_next_nodes_.erase( pair_it ); //note: erasing from vector, not ideal, but the vectors should usually be fairly small and the initialization shenanigans are only called once
				all_initialization_successful = false;
				break;
			}
		} // loop over dependent mcfi
		if ( all_initialization_successful ) {
			++pair_it;
		}
	}//loop over all non-redundant invrot groups
	//I guess we have to set the locations in this node in the child nodes

	tr << "Node initialization for geomcst " << geom_cst_ << ". After child node initialization, " << invrots_and_next_nodes_.size() << " non-redundant sets of inverse rotamers exist." << std::endl;
	for ( Size i =1; i <= invrots_and_next_nodes_.size(); ++i ) {
		tr << invrots_and_next_nodes_[i].first.size() << " inverse rotamers for non-redundant set " << i << std::endl;
		for ( auto child_it( invrots_and_next_nodes_[i].second.begin() ); child_it != invrots_and_next_nodes_[i].second.end(); ++child_it ) {
			(*child_it)->set_location_in_parent_node( i );
		}
	}

	//4.5 in case there were dead ends only, we return false
	if ( invrots_and_next_nodes_.size() == 0 ) return false;

	return true;
}

/// @details the real meat of this thing
/// see brief description in .hh file
core::scoring::constraints::ConstraintCOP
InvrotTreeNode::generate_constraints(
	core::pose::Pose const & pose,
	AllowedSeqposForGeomCstCOP geomcst_seqpos
) const
{
	core::id::AtomID fixed_pt( this->get_fixed_pt( pose ) );
	utility::vector1< core::scoring::constraints::ConstraintCOP > constraints_this_node;

	if ( !generate_invrot_csts_ ) runtime_assert( invrots_and_next_nodes_.size() == 1 ); //sanity check, relation exists because in this case we should only have one invrot

	for ( core::Size i = 1; i <= invrots_and_next_nodes_.size(); ++i ) {

		//utility::vector1< core::scoring::constraints::ConstraintCOP > constraints_this_invrot_set;
		utility::vector1< core::scoring::constraints::ConstraintCOP > constraints_this_invrot_node_pair;

		//1a. create the ambiguous constraint for this set of inverse rotamers
		if ( generate_invrot_csts_ ) constraints_this_invrot_node_pair.push_back( constrain_pose_res_to_invrots( invrots_and_next_nodes_[i].first, geomcst_seqpos->seqpos_for_geomcst( geom_cst_), pose, fixed_pt ) );

		//1b.
		for ( auto const & node_it : invrots_and_next_nodes_[i].second ) {
			core::scoring::constraints::ConstraintCOP this_child_cst( node_it->generate_constraints( pose, geomcst_seqpos ) );
			if ( this_child_cst ) constraints_this_invrot_node_pair.push_back( this_child_cst );
		}

		//1c.
		if ( constraints_this_invrot_node_pair.size() == 1 ) constraints_this_node.push_back( constraints_this_invrot_node_pair[1] );
		else if ( constraints_this_invrot_node_pair.size() > 1 ) constraints_this_node.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::MultiConstraint( constraints_this_invrot_node_pair ) ) );

	}//loop over node_pointer_pairs_

	if ( constraints_this_node.size() == 0 ) return nullptr;

	if ( constraints_this_node.size() == 1 ) return constraints_this_node[1];

	return core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint( constraints_this_node ) ) );
}


/// @details approach: get the target residues,
/// then find a residue in the pose that has the
/// same name as the first target residue and a
/// neighbor atom in the same position. this will
/// be the fixed point. exit w error if not found
core::id::AtomID
InvrotTreeNode::get_fixed_pt( core::pose::Pose const & pose ) const
{
	InvrotTreeNodeBaseCOP parent(this->parent_node() );
	if ( !parent ) utility_exit_with_message("the impossible just happened");

	utility::vector1< std::list< core::conformation::ResidueCOP > > parent_res( parent->all_target_residues( get_self_weak_ptr() ) );
	std::string target_name3( (*parent_res[1].begin())->name3() );
	core::Vector const & xyz_to_find( (*parent_res[1].begin())->nbr_atom_xyz() );

	for ( Size i = 1; i<= pose.size(); ++i ) {
		if ( pose.residue_type(i).name3() == target_name3 ) {
			for ( Size j =1; j<= pose.residue(i).atoms().size(); ++j ) {
				if ( xyz_to_find.distance_squared( pose.residue(i).atom(j).xyz() ) < 0.1 ) {
					return core::id::AtomID( j, i );
				}
			}
		}
	}
	//if we get here, that means no atom was found and shit's fucked up somewhere
	utility_exit_with_message("No success when trying to find a fixed pt in the InvrotTree that's also in the pose. Something is profoundly broken.");
	return core::id::AtomID( 0, 0 ); // to pacify compiler
}


utility::vector1< std::list< core::conformation::ResidueCOP > >
InvrotTreeNode::all_target_residues( InvrotTreeNodeBaseCAP child_node ) const
{

	utility::vector1< std::list< core::conformation::ResidueCOP > > to_return;

	//1. get the target residues of the parent node
	InvrotTreeNodeBaseCOP parent = this->parent_node().lock();
	if ( parent ) to_return = parent->all_target_residues( get_self_weak_ptr() );
	//2. add the target residues from this node
	//corresponding to the asking child node
	bool child_found(false);
	for ( Size i =1; i <= invrots_and_next_nodes_.size(); ++i ) {
		for ( Size j =1; j <= invrots_and_next_nodes_[i].second.size(); ++j ) {
			if ( utility::pointer::equal(child_node, invrots_and_next_nodes_[i].second[j]) ) {
				to_return.push_back( invrots_and_next_nodes_[i].first );
				child_found = true;
				break;
			}
		}
		if ( child_found ) break;
	}

	return to_return;
}


/// @details if covalent is true, clashes will not be checked
/// for the last vector that comes down from the parent.
/// kinda crude, could be made better, i.e. only except the actual
/// constrained atoms from clash check
void
InvrotTreeNode::remove_invrots_clashing_with_parent_res(
	std::list< core::conformation::ResidueCOP > & invrots,
	bool covalent
) const
{

	InvrotTreeNodeBaseCOP parent = this->parent_node().lock();

	//in case there's no parent, i.e. in the unit test
	//or potential use of this outside the tree, there
	//are no target to clash with, so we'll return
	if ( !parent ) return;

	//1. get all parent rotamers
	utility::vector1< std::list< core::conformation::ResidueCOP > > all_targets ( parent->all_target_residues( get_self_weak_ptr() ) );

	if ( covalent ) all_targets.pop_back(); //covalent: we ignore the last list, since this represents the immediate parent residues

	//safety: if for some reason all_targets is empty,
	//we return
	if ( all_targets.size() == 0 ) return;

	//2.clash check between invrots and parent residues
	//all invrots that clash with all rotamers from any
	//of the all_targets list will be removed
	for ( auto invrot_it( invrots.begin() ); invrot_it != invrots.end(); /*increment in loop*/ ) {
		bool cur_invrot_clashes(false);

		for ( Size i =1; i <= all_targets.size(); ++i ) {

			//safety: if for some reason this list is empty, we continue
			if ( all_targets[i].size() == 0 ) continue;

			bool this_list_all_clash( true );
			for (  std::list <core::conformation::ResidueCOP >::const_iterator target_rot_it( all_targets[i].begin() ), target_rot_it_end( all_targets[i].end() ); target_rot_it != target_rot_it_end; ++target_rot_it ) {

				//the actual clash check
				//hardcoded as hadrsphere here for now, but in the future this could be
				//combined with matcher filters
				bool these_two_clash(false);
				core::Real cutoff_sq = 2.5 * 2.5;
				core::conformation::Residue const & res1(**invrot_it), res2(**target_rot_it);
				for ( Size res1at = 1; res1at <= res1.nheavyatoms(); ++res1at ) {
					if ( res1.type().is_virtual( res1at ) ) continue;
					for ( Size res2at = 1; res2at <= res2.nheavyatoms(); ++res2at ) {
						if ( res2.type().is_virtual( res2at ) ) continue;
						if ( res1.atom( res1at ).xyz().distance_squared( res2.atom( res2at ).xyz() ) < cutoff_sq ) {
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
					if ( these_two_clash ) break; //breaks loop over res 1 atoms
				}//loop over res1 atoms
				if ( !these_two_clash ) {
					this_list_all_clash = false;
					break; //breaks loop over rots in this list
				}
			} //for loop over target_rot
			if ( this_list_all_clash ) {
				cur_invrot_clashes = true;
				break; // breaks loop over target_lists
			}
		} //for loop over target_rot_lists

		if ( cur_invrot_clashes ) invrot_it = invrots.erase( invrot_it );
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

	for ( Size i =1; i <= invrot_collectors.size(); ++i ) {
		std::map< InvrotTreeNodeBaseCOP, Size > const & collector_map( invrot_collectors[i]->owner_nodes_and_locations() );
		//safety
		if ( collector_map.find( get_self_ptr() ) != collector_map.end() ) utility_exit_with_message("InvrotTreeNode being asked to fill a collector that it apparently already previously filled.");

		auto map_it( collector_map.find( this->parent_node().lock() ) );

		if ( map_it != collector_map.end() ) {
			if ( map_it->second == this->location_in_parent_node() ) {
				if ( invrot_collectors[i]->invrots()[ geom_cst_].size() == 0 ) empty_parent_collectors.push_back( i );
				else filled_parent_collectors.push_back( i );
			}
		}
	}
	tr << "Collecting inverse rotamers for a node from geomcst " << geom_cst_ << ". " << empty_parent_collectors.size() << " empty parent collectors and " << filled_parent_collectors.size() << " filled parent collectors were found." << std::endl;

	//1.b make sure we found at least one parent
	Size total_parents( empty_parent_collectors.size() + filled_parent_collectors.size() );
	runtime_assert( total_parents != 0);

	//2a for every parent, we need to add rotamers for every definiton
	utility::vector1< Size > collectors_to_fill( empty_parent_collectors );
	for ( Size i =1; i <= filled_parent_collectors.size(); ++ i ) {
		invrot_collectors.push_back( invrot_collectors[ filled_parent_collectors[i] ]->clone() );
		collectors_to_fill.push_back( invrot_collectors.size() );
	}

	//2b. fill 'er up
	for ( Size i = 1; i <= collectors_to_fill.size(); ++i ) {
		Size overflow_start( invrot_collectors.size());
		for ( Size j = 2; j <= invrots_and_next_nodes_.size(); ++j ) {
			invrot_collectors.push_back( invrot_collectors[ collectors_to_fill[i] ]->clone() );
		}

		invrot_collectors[ collectors_to_fill[ i ] ]->set_invrots_for_listnum( geom_cst_, invrots_and_next_nodes_[1].first, get_self_ptr(), 1 );

		for ( Size j =2; j <= invrots_and_next_nodes_.size(); ++j ) {
			invrot_collectors[ overflow_start + j - 1 ]->set_invrots_for_listnum( geom_cst_, invrots_and_next_nodes_[j].first, get_self_ptr(), j );
		}
	}

	//3. call this shit on all daughter nodes
	for ( Size i = 1; i <= invrots_and_next_nodes_.size(); ++i ) {
		for ( auto const & child_it : invrots_and_next_nodes_[i].second ) {
			child_it->collect_all_inverse_rotamers( invrot_collectors );
		}
	}
}

}
}
}
