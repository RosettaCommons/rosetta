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
#include <utility/vector0.hh>
#include <utility/vector1.hh>


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
	docking_ft_( false ),
	simple_ft_( false ),
	jump_( 1 ),
	resnum_( "" ),
	connect_to_( "" ),
	anchor_res_( "" ),
	connect_from_( "" ),
	host_chain_( 2 )
{}

SetAtomTree::~SetAtomTree() {}

protocols::moves::MoverOP
SetAtomTree::clone() const {
	return (protocols::moves::MoverOP( new SetAtomTree( *this ) ) );
}

void
SetAtomTree::parse_my_tag( TagPtr const tag, DataMap &, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	docking_ft_ = tag->getOption< bool >("docking_ft", 0 );
	simple_ft( tag->getOption< bool >( "simple_ft", 0 ) );
	jump_ = tag->getOption< core::Size >( "jump", 1);
	if( docking_ft_ ) return;
	/// resnum & pdb_num are now equivalent
	if( tag->hasOption( "resnum" ) )
		resnum_ = tag->getOption< std::string > ("resnum" );
	else if( tag->hasOption( "pdb_num" ) ){
		resnum_ = tag->getOption< std::string > ("pdb_num" );
		connect_to_ = tag->getOption< std::string >( "connect_to","" );
	}
	if( tag->hasOption( "anchor_res" ) ){
		anchor_res_ = tag->getOption< std::string > ( "anchor_res" );
		if( tag->hasOption( "connect_from" ) )
			connect_from_ = tag->getOption< std::string >( "connect_from" );
	}
	host_chain_ = tag->getOption< core::Size >( "host_chain", 2);
	TR<<"resnum: "<<resnum_<<" anchor: "<< anchor_res_<<std::endl;
}//end parse my tag

core::kinematics::FoldTreeOP
SetAtomTree::create_atom_tree( core::pose::Pose const & pose, core::Size const host_chain, core::Size const resnum, core::Size const anchor_num_in, std::string connect_to_in/*=""*/, std::string connect_from_in/*=""*/ )
{
	core::kinematics::FoldTreeOP fold_tree = new core::kinematics::FoldTree;

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
	if ( connect_from == "" ) connect_from = optimal_connection_point( pose.residue( jump_pos1 ).name3() );
	fold_tree->clear();
	fold_tree->add_edge( jump_pos1, jump_pos2, rb_jump );
	fold_tree->add_edge( 1, jump_pos1, kinematics::Edge::PEPTIDE );
	fold_tree->add_edge( jump_pos1, pose.conformation().chain_end( 1 ), kinematics::Edge::PEPTIDE );
	fold_tree->add_edge( pose.conformation().chain_begin( 2 ), jump_pos2, kinematics::Edge::PEPTIDE );
	fold_tree->add_edge( jump_pos2, pose.total_residue(), kinematics::Edge::PEPTIDE );
	TR<<"CONNECT_FROM: "<<connect_from<<"and  CONNECT TO: " <<connect_to<<std::endl;

  fold_tree->set_jump_atoms( rb_jump, connect_from , connect_to );
	fold_tree->reorder( 1 );
	return( fold_tree );
}

void
SetAtomTree::apply( core::pose::Pose & pose )
{
	if( docking_ft_ ){
		utility::vector1< core::SSize > jumps;
		jumps.clear();
		jumps.push_back( jump_ );
		std::string const partners( "_" );
		protocols::docking::setup_foldtree( pose, partners, jumps );
		TR<<"Setting up docking foldtree over jump "<<jump_<<'\n';
		TR<<"Docking foldtree: "<<pose.fold_tree()<<std::endl;
		return;
	}
	if( simple_ft() ){
		core::kinematics::FoldTree new_ft;
		new_ft.clear();
		core::conformation::Conformation const conf( pose.conformation() );
		core::Size jump( 1 );
		for( core::Size chain = 1; chain <= conf.num_chains(); ++chain ){
			new_ft.add_edge( conf.chain_begin( chain ), conf.chain_end( chain ), -1 );
			if( chain > 1 ){
				new_ft.add_edge( conf.chain_begin( chain - 1 ), conf.chain_begin( chain ), jump );
				jump++;
			}
		}
		new_ft.reorder( 1 );
		TR<<new_ft<<std::endl;
		pose.fold_tree( new_ft );
		TR<<"Simple tree: "<<pose.fold_tree()<<std::endl;
		return;
	}

	core::Size const resnum( protocols::rosetta_scripts::parse_resnum( resnum_, pose ) );
	core::conformation::Residue const res_central( pose.residue( resnum ) );;

  std::string connect_to( connect_to_ );
	if( connect_to_ == "" ){
    connect_to = optimal_connection_point( res_central.name3() );
    TR<<"connect_to not defined by user. Defaulting to "<<connect_to<<std::endl;
	}

	core::Size anchor_num( 0 );
	std::string connect_from( connect_from_ );
	if ( anchor_res_ != "" ) {
			core::pose::PDBPoseMap const pose_map( pose.pdb_info()->pdb2pose() );
    	char const chain( anchor_res_[ anchor_res_.length() - 1 ] );
			std::stringstream ss( anchor_res_.substr( 0, anchor_res_.length() - 1 ) );
      core::Size number;
      ss >> number;
      anchor_num  = pose_map.find( chain, number );
		  TR<<"anchor_num::: " << anchor_num << "and pdbnum:::" << resnum_ <<std::endl;

			core::conformation::Residue const res_anchor( pose.residue( anchor_num ) );;
			if( connect_from_ == "" )
				connect_from = optimal_connection_point( res_anchor.name3() );
		}	//end anchor reside


	TR<<"Previous fold tree: "<<pose.fold_tree()<<'\n';
	pose.fold_tree( *create_atom_tree( pose, host_chain_, resnum, anchor_num , connect_to, connect_from ) );
	TR<<"New fold tree: "<<pose.fold_tree()<<std::endl;
}

std::string
SetAtomTree::get_name() const {
	return SetAtomTreeCreator::mover_name();
}

} //movers
} //protein_interface_design
} //protocols

