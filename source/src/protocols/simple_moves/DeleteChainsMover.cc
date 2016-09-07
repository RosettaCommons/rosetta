// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteChainsMover.cc
/// @brief

// Unit headers
#include <protocols/simple_moves/DeleteChainsMover.hh>
#include <protocols/simple_moves/DeleteChainsMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.DeleteChainsMover" );

namespace protocols {
namespace simple_moves {

DeleteChainsMover::DeleteChainsMover()
: moves::Mover("DeleteChainsMover")
{
	set_defaults();
}

DeleteChainsMover::DeleteChainsMover( std::string const & chains, core::pose::Pose const & pose)
: moves::Mover("DeleteChainsMover")
{
	set_chains( chains, pose );
	set_defaults();
}

DeleteChainsMover::DeleteChainsMover( utility::vector1< core::Size > const & chains )
: moves::Mover("DeleteChainsMover")
{
	set_chains( chains );
	set_defaults();
}

void
DeleteChainsMover::set_defaults(){
	set_detect_bonds( true );
	set_detect_pseudobonds( true );
}

void
DeleteChainsMover::set_chains( std::string const & chains, core::pose::Pose const & pose ){
	chains_.clear();
	for ( char chain : chains ) {
		core::Size chain_id = core::pose::get_chain_id_from_chain( chain, pose );
		chains_.push_back(chain_id);
	}
}

void
DeleteChainsMover::set_chains( utility::vector1< core::Size > const & chains ){
	chains_ = chains;
}

utility::vector1< core::Size >
DeleteChainsMover::chains() const {
	return chains_;
}

void
DeleteChainsMover::set_detect_bonds(bool detect_bonds) {
	detect_bonds_ = detect_bonds;
}

bool
DeleteChainsMover::detect_bonds() const {
	return detect_bonds_;
}

void
DeleteChainsMover::set_detect_pseudobonds(bool detect_pseudobonds) {
	detect_pseudobonds_ = detect_pseudobonds;
}

bool
DeleteChainsMover::detect_pseudobonds() const {
	return detect_pseudobonds_;
}

void
DeleteChainsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	set_defaults();
	if ( ! tag->hasOption( "chains") ) {
		utility_exit_with_message("Must pass chains tag to DeleteChainsMover...");
	}

	set_chains( tag->getOption< std::string >( "chains"), pose );

	set_detect_bonds( tag->getOption< bool >( "detect_bonds", detect_bonds_) );
	set_detect_pseudobonds( tag->getOption< bool >( "detect_pseudobonds", detect_pseudobonds_ ) );
}

void
DeleteChainsMover::apply( Pose & pose )
{
	if ( chains_.size() == 0 ) {
		utility_exit_with_message("DeleteChainsMover requires chains to be set...");
	}

	for ( core::Size i = 1; i <= chains_.size(); ++i ) {

		pose.conformation().delete_residue_range_slow( pose.conformation().chain_begin( chains_ [ i ] ), pose.conformation().chain_end( chains_ [ i ] ) );
	}

	pose.pdb_info()->obsolete( false );

	if ( pose.is_fullatom() ) {

		//Same order of detection as import pdb.
		if ( detect_bonds_ ) {
			pose.conformation().detect_bonds();
		}
		pose.conformation().detect_disulfides();
		if ( detect_pseudobonds_ ) {
			pose.conformation().detect_pseudobonds();
		}
	}

	core::pose::set_reasonable_fold_tree( pose );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
DeleteChainsMoverCreator::keyname() const
{
	return DeleteChainsMoverCreator::mover_name();
}

protocols::moves::MoverOP
DeleteChainsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeleteChainsMover );
}

std::string
DeleteChainsMoverCreator::mover_name()
{
	return "DeleteChainsMover";
}


std::string
DeleteChainsMover::get_name() const {
	return DeleteChainsMoverCreator::mover_name();
}

moves::MoverOP
DeleteChainsMover::clone() const
{
	return moves::MoverOP( new DeleteChainsMover( *this ) );
}

moves::MoverOP
DeleteChainsMover::fresh_instance() const
{
	return moves::MoverOP( new DeleteChainsMover );
}

} // simple_moves
} // protocols
