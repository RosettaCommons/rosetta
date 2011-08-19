// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/SetAtomTree.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/SetAtomTree.hh>
#include <protocols/protein_interface_design/movers/SetAtomTreeCreator.hh>
#include <protocols/docking/util.hh>

// Package headers
#include <protocols/protein_interface_design/util.hh>

// Project headers
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <protocols/rosetta_scripts/util.hh>
//Auto Headers
#include <core/chemical/AtomType.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.SetAtomTree" );

std::string
SetAtomTreeCreator::keyname() const
{
	return SetAtomTreeCreator::mover_name();
}

protocols::moves::MoverOP
SetAtomTreeCreator::create_mover() const {
	return new SetAtomTree;
}

std::string
SetAtomTreeCreator::mover_name()
{
	return "AtomTree";
}

//initializing
core::Size anchor_num = 0;
std::string connect_from("");
std::string connect_to("");

SetAtomTree::SetAtomTree() :
	protocols::moves::Mover( SetAtomTreeCreator::mover_name() ),
	fold_tree_( NULL )
{}

SetAtomTree::SetAtomTree( core::kinematics::FoldTreeOP ft ) :
	protocols::moves::Mover( SetAtomTreeCreator::mover_name() )
{
	fold_tree_ = ft;
}
SetAtomTree::~SetAtomTree() {}

protocols::moves::MoverOP
SetAtomTree::clone() const {
	return (protocols::moves::MoverOP( new SetAtomTree( *this ) ) );
}


void
SetAtomTree::fold_tree( core::kinematics::FoldTreeOP ft ) { fold_tree_ = ft; }

core::kinematics::FoldTreeOP
SetAtomTree::fold_tree() const { return( fold_tree_ ); }


void
SetAtomTree::parse_my_tag( TagPtr const tag, DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	if( tag->getOption< bool >( "docking_ft", 0 ) ){//just generate a docking default foldtree
		core::pose::Pose nonconst_pose( pose );
		core::SSize const jump( tag->getOption< core::SSize >( "jump", 1 ));
		utility::vector1< core::SSize > jumps;
		jumps.clear();
		jumps.push_back( jump );
		std::string const partners( "_" );
		protocols::docking::setup_foldtree( nonconst_pose, partners, jumps );
		fold_tree_ = new core::kinematics::FoldTree( nonconst_pose.fold_tree() );
		TR<<"Setting up docking foldtree over jump "<<jump<<'\n';
		TR<<*fold_tree_;
		return;
	}
	core::Size const resnum( protocols::rosetta_scripts::get_resnum( tag, pose ) );
	core::conformation::Residue const res_central( pose.residue( resnum ) );;

	if( tag->hasOption( "connect_to" ) ){
		connect_to = tag->getOption<string>( "connect_to" );
		TR<<"USER DEFINED connect to atom : "<<connect_to<<std::endl;
	}
	else {
    connect_to = optimal_connection_point( res_central.name3() );
    TR<<"connect_to not defined by user. Defaulting to "<<connect_to<<std::endl;
	}

	if ( tag->hasOption( "anchor_res" ) ) {
			std::string anchor_string( tag->getOption<std::string>("anchor_res" ) );
			core::pose::PDBPoseMap const pose_map( pose.pdb_info()->pdb2pose() );
    	char const chain( anchor_string[ anchor_string.length() - 1 ] );
			std::stringstream ss( anchor_string.substr( 0, anchor_string.length() - 1 ) );
      core::Size number;
      ss >> number;
      anchor_num  = pose_map.find( chain, number );
		  TR<<"anchor_num::: " << anchor_num << "and pdbnum:::" << resnum <<std::endl;
			core::conformation::Residue const anchor_res( pose.residue( anchor_num ) );

			if( tag->hasOption( "connect_from" ) ) {
			 connect_from = tag->getOption<string>( "connect_from" );
			 TR<<"USER-DEFINED connect from atom : "<<connect_from<<std::endl;
			}
			else {
			 connect_from = optimal_connection_point( anchor_res.name3() );
			 TR<<"DEFAULTING connection from : "<<connect_from<<std::endl;
			}
		}	//end anchor reside

	core::Size const host_chain( tag->getOption<core::Size>( "host_chain", 2 ) );

	create_atom_tree( pose, host_chain, resnum, anchor_num , connect_to, connect_from );
	TR<<"AtomTree saving the following atom tree:\n"<<*fold_tree_<<std::endl;
	TR<<"resnum: "<<resnum<<" anchor: "<< anchor_num<<std::endl;
}//end parse my tag

core::kinematics::FoldTreeOP
SetAtomTree::create_atom_tree( core::pose::Pose const & pose, core::Size const host_chain, core::Size const resnum, core::Size const anchor_num_in, std::string connect_to_in/*=""*/, std::string connect_from_in/*=""*/ )
{
	fold_tree_ = new core::kinematics::FoldTree;

	core::Size const begin( pose.conformation().chain_begin( host_chain == 1 ? 2 : 1 ) );
	core::Size const end( pose.conformation().chain_end( host_chain == 1 ? 2 : 1 ) );
	core::Real min_dist(10000);
	core::conformation::Residue const res_central( pose.residue( resnum ) );
	core::Size anchor_num = anchor_num_in;
	std::string connect_to( connect_to_in );
	std::string connect_from( connect_from_in );

	if( connect_to == "" )
		connect_to = optimal_connection_point( res_central.name3() );

	core::Size nearest_res( 0 );

	if ( anchor_num == 0 ) {
		for( core::Size res = begin+1; res <= end-1; ++res ) {
			core::conformation::Residue const res2( pose.residue(res) );
			core::Real const distance( res_central.xyz( res_central.nbr_atom() ).distance( res2.xyz( res2.nbr_atom() ) ) );
			if( distance<=min_dist ) {
				min_dist = distance;
				nearest_res = res;
				TR<<"anchor defaults to: "<<nearest_res<<std::endl;
			}
		}
	}//end no anchor specified, hence finds the closest one
	else
		 nearest_res = anchor_num;

	core::Size const rb_jump( 1 );
	core::Size const jump_pos1( host_chain == 1 ? resnum : nearest_res );
	core::Size const jump_pos2( host_chain == 1 ? nearest_res : resnum );
	fold_tree_->clear();
	fold_tree_->add_edge( jump_pos1, jump_pos2, rb_jump );
	fold_tree_->add_edge( 1, jump_pos1, kinematics::Edge::PEPTIDE );
	fold_tree_->add_edge( jump_pos1, pose.conformation().chain_end( 1 ), kinematics::Edge::PEPTIDE );
	fold_tree_->add_edge( pose.conformation().chain_begin( 2 ), jump_pos2, kinematics::Edge::PEPTIDE );
	fold_tree_->add_edge( jump_pos2, pose.total_residue(), kinematics::Edge::PEPTIDE );
	TR<<"CONNECT_FROM: "<<connect_from<<"and  CONNECT TO: " <<connect_to<<std::endl;

	if ( connect_from == "" ) connect_from = optimal_connection_point( pose.residue( jump_pos1 ).name3() );
  fold_tree_->set_jump_atoms( rb_jump, connect_from , connect_to );
	fold_tree_->reorder( 1 );
	return( fold_tree_ );
}

void
SetAtomTree::apply( core::pose::Pose & pose )
{
	runtime_assert( fold_tree_ );
	TR<<"Previous fold tree: "<<pose.fold_tree()<<'\n';
	pose.fold_tree( *fold_tree_ );
	TR<<"New fold tree: "<<pose.fold_tree()<<std::endl;
}

std::string
SetAtomTree::get_name() const {
	return SetAtomTreeCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

