// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTarget.hh
/// @brief  .hh file for inverse rotamer target
/// @author Florian Richter, flosopher@gmail.com, mar 2012


//unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/MultiConstraint.hh>

#include <basic/Tracer.hh>

//utility headers
#include <utility/string_util.hh>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static thread_local basic::Tracer tr( "protocols.toolbox.match_enzdes_util.InvrotTarget" );

InvrotTarget::InvrotTarget()
	: InvrotTreeNodeBase( InvrotTreeNodeBaseCAP() )
{
	representative_target_res_for_geom_cst_.clear();
	all_target_res_.clear();
	next_nodes_.clear();
}

InvrotTarget::~InvrotTarget(){}

core::conformation::ResidueCOP
InvrotTarget::target_res_for_geom_cst( core::Size geom_cst ) const
{
	return representative_target_res_for_geom_cst_[ geom_cst ];
}

std::list< core::conformation::ResidueCOP >
InvrotTarget::all_target_res() const
{
	return all_target_res_;
}

core::scoring::constraints::ConstraintCOP
InvrotTarget::generate_constraints(
    core::pose::Pose const & pose,
    AllowedSeqposForGeomCstCOP geomcst_seqpos
) const
{
	if( next_nodes_.size() == 0 ) utility_exit_with_message("generate constraints function called on InvrotTraget that's not pointing to any next_nodes, something's wrong somewhere.");

	utility::vector1< core::scoring::constraints::ConstraintCOP > node_constraints;

	for( Size i =1; i <= next_nodes_.size(); ++i ){
		core::scoring::constraints::ConstraintCOP this_child_cst( (next_nodes_[i]->generate_constraints( pose, geomcst_seqpos ) ) );
		if( this_child_cst) node_constraints.push_back( this_child_cst );
	}
	if( node_constraints.size() == 0 ) return NULL;

	return core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::MultiConstraint( node_constraints ) );
}


/// @details the meat...
bool
InvrotTarget::initialize_tree_nodes_from_enzcst_io(
	EnzConstraintIOCOP enzcst_io )
{
	next_nodes_.clear();
	representative_target_res_for_geom_cst_.clear();
	this->generate_representative_target_res_for_geom_cst( enzcst_io->num_mcfi_lists() );
	if( representative_target_res_for_geom_cst_.size() != enzcst_io->num_mcfi_lists() ){
		utility_exit_with_message("child class failed to properly implement 'generate_representative_target_res_for_geom_cst' method.");
	}

	utility::vector1< Size > mcfi_to_build;
	for( Size i =1; i <= enzcst_io->num_mcfi_lists(); ++i){
		//make sure this geom cst is concerned with the downstream residue / object
		std::pair< Size, Size> const & target_res( enzcst_io->target_downstream_res()[i] );
		if( (target_res.first == 1 ) && (target_res.second == 1 ) ){
			mcfi_to_build.push_back( i );
			next_nodes_.push_back( InvrotTreeNodeOP( new InvrotTreeNode( get_self_weak_ptr() ) ) );
		}
	}

	for( Size i =1; i <= next_nodes_.size(); ++i ){
		if( !next_nodes_[i]->initialize_from_enzcst_io( *representative_target_res_for_geom_cst_[ mcfi_to_build[i] ], enzcst_io, mcfi_to_build[i]) ){
			tr << "InvrotTarget could not initialize because node initialization for geom_cst " << mcfi_to_build[i] << " failed." << std::endl;
			return false;
		}
	}
	return true;
}

utility::vector1< std::list< core::conformation::ResidueCOP > >
InvrotTarget:: all_target_residues( InvrotTreeNodeBaseCAP /*child_node*/ ) const
{
	utility::vector1< std::list< core::conformation::ResidueCOP > > to_return;
	to_return.push_back( all_target_res_ );
	return to_return;
}


void
InvrotTarget::collect_all_inverse_rotamers(
	utility::vector1< InvrotCollectorOP > & invrot_collectors
) const
{
	//1. make space
	Size num_residue_lists( representative_target_res_for_geom_cst_.size() + 1 ); //+1 bc we're also counting the ligand now
	Size input_size( invrot_collectors.size() );
	invrot_collectors.push_back(  utility::pointer::shared_ptr<class protocols::toolbox::match_enzdes_util::InvrotCollector>( new InvrotCollector( num_residue_lists ) ) );

	//2. put target res into 0th element
	invrot_collectors[ invrot_collectors.size() ]->set_invrots_for_listnum( 0, all_target_res_, get_self_ptr(), 1 );

	//3. collect daughter node invrots
	for( Size i = 1; i <= next_nodes_.size(); ++i ){
		next_nodes_[i]->collect_all_inverse_rotamers( invrot_collectors );
	}

	//4. some kind of sanity check to make sure this shit actually worked maybe?
	Size output_size( invrot_collectors.size() );
	for( Size i = input_size + 1; i <= output_size; ++i ){

		if( invrot_collectors[ i ]->invrots().size() != num_residue_lists ){
			utility_exit_with_message("Tree definition "+utility::to_string( i )+" does not contain the necessary "+utility::to_string( num_residue_lists )+" residue lists.");
		}

		for( Size j = 0; j < num_residue_lists; ++j){
			if( invrot_collectors[i]->invrots()[j].size() == 0 ) utility_exit_with_message("Tree definition "+utility::to_string( i )+" does not contain rotamers for residue list "+utility::to_string( j )+".");
		}
		//other sanity checks necessary?
	}
}


void
InvrotTarget::set_all_target_res( std::list< core::conformation::ResidueCOP > const & all_target_res )
{
	all_target_res_ = all_target_res;
}

void
InvrotTarget::set_representative_target_res_for_geom_cst(
	utility::vector1< core::conformation::ResidueCOP > const & representative_res )
{
	representative_target_res_for_geom_cst_ = representative_res;
}


SingleResidueInvrotTarget::SingleResidueInvrotTarget(
	utility::vector1< core::conformation::ResidueCOP > const & single_res )
{
	//this might seem a bit clumsy, but at some point we have to change from vector to list
	std::list< core::conformation::ResidueCOP > all_targets;
	for( utility::vector1< core::conformation::ResidueCOP >::const_iterator res_it( single_res.begin() ); res_it != single_res.end(); ++res_it ) all_targets.push_back( *res_it );
	this->set_all_target_res( all_targets );
}

SingleResidueInvrotTarget::~SingleResidueInvrotTarget(){}

void
SingleResidueInvrotTarget::generate_representative_target_res_for_geom_cst( Size const num_geom_cst )
{
	utility::vector1< core::conformation::ResidueCOP > representative_res;
	for( Size i =1; i <= num_geom_cst; ++i ) representative_res.push_back( (*this->all_target_res().begin()) );
	this->set_representative_target_res_for_geom_cst( representative_res );
}

}
}
}
