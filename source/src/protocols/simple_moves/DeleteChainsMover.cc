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
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/simple_moves/DeleteChainsMover.hh>
#include <protocols/simple_moves/DeleteChainsMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.simple_moves.DeleteChainsMover" );

namespace protocols {
namespace simple_moves {

DeleteChainsMover::DeleteChainsMover()
: moves::Mover("DeleteChainsMover")
{
	set_defaults();
}

DeleteChainsMover::DeleteChainsMover( utility::vector1< std::string > const & chains )
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

/// @brief Provide the citation.
void
DeleteChainsMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		"DeleteChainsMover", basic::citation_manager::CitedModuleType::Mover
		)
	);
	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1371/journal.pcbi.1006112" ) );
	citations.add( cc );
}

void
DeleteChainsMover::set_chains( utility::vector1< std::string > const & chains ){
	chains_ = chains;
}

utility::vector1< std::string > const &
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
	basic::datacache::DataMap &
) {
	set_defaults();
	if ( ! tag->hasOption( "chains") ) {
		utility_exit_with_message("Must pass chains tag to DeleteChainsMover...");
	}

	std::string chain_string = tag->getOption< std::string >( "chains");
	// Because of the input format, we're somewhat stuck with single letter chains
	utility::vector1< std::string > chain_vector;
	for ( char c: chain_string ) {
		chain_vector.push_back( std::string{c} );
	}
	set_chains( chain_vector );

	set_detect_bonds( tag->getOption< bool >( "detect_bonds", detect_bonds_) );
	set_detect_pseudobonds( tag->getOption< bool >( "detect_pseudobonds", detect_pseudobonds_ ) );
}

void
DeleteChainsMover::apply( Pose & pose )
{
	if ( chains_.size() == 0 ) {
		utility_exit_with_message("DeleteChainsMover requires chains to be set...");
	}

	for ( std::string const & chain : chains_ ) {
		core::Size chain_id = core::pose::get_chain_id_from_chain( chain, pose );

		pose.conformation().delete_residue_range_slow( pose.conformation().chain_begin( chain_id ), pose.conformation().chain_end( chain_id ) );

		//JAB - put this in inner loop to prevent segfaults.  Delete each chain one at a time.
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


}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





moves::MoverOP
DeleteChainsMover::clone() const
{
	return utility::pointer::make_shared< DeleteChainsMover >( *this );
}

moves::MoverOP
DeleteChainsMover::fresh_instance() const
{
	return utility::pointer::make_shared< DeleteChainsMover >();
}

std::string DeleteChainsMover::get_name() const {
	return mover_name();
}

std::string DeleteChainsMover::mover_name() {
	return "DeleteChainsMover";
}

void DeleteChainsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "chains", xs_string, "delete these chains.  The format is unseparated chain letters, like 'AHLR'") //XRW TODO, I guess this could have a restriction, but it's not clear that the PDB rules for what you can put in the chain column are really that valid here
		+ XMLSchemaAttribute( "detect_bonds", xsct_rosetta_bool, "detect and delete broken bonds afterwards")
		+ XMLSchemaAttribute( "detect_pseudobonds", xsct_rosetta_bool, "detect and delete broken pseudobonds afterwards");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"remove chains from a pose", attlist );
}

std::string DeleteChainsMoverCreator::keyname() const {
	return DeleteChainsMover::mover_name();
}

protocols::moves::MoverOP
DeleteChainsMoverCreator::create_mover() const {
	return utility::pointer::make_shared< DeleteChainsMover >();
}

void DeleteChainsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeleteChainsMover::provide_xml_schema( xsd );
}


} // simple_moves
} // protocols
