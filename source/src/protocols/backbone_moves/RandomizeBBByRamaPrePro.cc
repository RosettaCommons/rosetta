// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/RandomizeBBByRamaPrePro.cc
/// @brief A simple mover to randomize a backbone, or a portion of a backbone, biased by the rama_prepro score of each residue.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Unit headers
#include <protocols/backbone_moves/RandomizeBBByRamaPrePro.hh>
#include <protocols/backbone_moves/RandomizeBBByRamaPreProCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/id/TorsionID.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.backbone_moves.RandomizeBBByRamaPrePro" );

namespace protocols {
namespace backbone_moves {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
RandomizeBBByRamaPrePro::RandomizeBBByRamaPrePro():
	protocols::moves::Mover( RandomizeBBByRamaPrePro::mover_name() ),
	selector_(nullptr)
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
RandomizeBBByRamaPrePro::RandomizeBBByRamaPrePro( RandomizeBBByRamaPrePro const & src ):
	protocols::moves::Mover( src ),
	selector_( src.selector_==nullptr ? nullptr : src.selector_->clone() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
RandomizeBBByRamaPrePro::~RandomizeBBByRamaPrePro(){}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RandomizeBBByRamaPrePro::clone() const
{
	return utility::pointer::make_shared< RandomizeBBByRamaPrePro >( *this );
}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
RandomizeBBByRamaPrePro::apply( core::pose::Pose& pose ){

	using namespace core::chemical;

	//Select residues:
	core::select::residue_selector::ResidueSubset selection;
	if ( selector_ != nullptr ) {
		selection = selector_->apply(pose);
		debug_assert( selection.size() == pose.total_residue() );
	}

	//Get RamaPrePro energy:
	core::scoring::RamaPrePro const &rama( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );

	//Loop through and randomize each residue.
	for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
		if ( selector_ == nullptr || selection[i] ) {
			utility::vector1 < core::Real > tors_vals;
			core::Size const nextres_index( pose.residue(i).connected_residue_at_upper() );
			core::chemical::ResidueTypeCOP next_res(
				nextres_index != 0 ?
				pose.residue_type_ptr(nextres_index) :
				pose.residue_type_set_for_pose( FULL_ATOM_t )->name_mapOP("ALA")
			);
			rama.random_mainchain_torsions( pose.conformation(), pose.residue_type_ptr(i), next_res, tors_vals );
			for ( core::Size j(1), jmax(tors_vals.size()); j<=jmax; ++j ) {
				pose.set_torsion( core::id::TorsionID(i, core::id::BB, j), tors_vals[j] );
			}
		}
	}
	pose.update_residue_neighbors();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
RandomizeBBByRamaPrePro::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
	if ( selector_ != nullptr ) {
		output << "ResidueSelector of type " << selector_->get_name() << " provided." << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
RandomizeBBByRamaPrePro::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap
) {
	set_residue_selector( core::select::residue_selector::parse_residue_selector( tag, datamap ) );
}

void RandomizeBBByRamaPrePro::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "An optional residue selector.  If provided, the mover will only randomize the subset of residues that are selected.  If not provided, then the mover is applied to the whole pose." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "The RandomizeBBByRamaPrePro mover randomizes a backbone, or selected residues in a backbone, biased by the Ramachandran preferences of each amino acid, using the RamaPrePro energy term (which uses different scoring tables for positions before proline/peptoids than it does for other positions).", attlist );
}

/// @brief Set the residue selector that this mover will use.
/// @details The selector is cloned.
void
RandomizeBBByRamaPrePro::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP const & selector_in
) {
	runtime_assert( selector_in != nullptr );
	selector_ = selector_in->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
RandomizeBBByRamaPrePro::fresh_instance() const
{
	return protocols::moves::MoverOP( new RandomizeBBByRamaPrePro );
}

std::string RandomizeBBByRamaPrePro::get_name() const {
	return mover_name();
}

std::string RandomizeBBByRamaPrePro::mover_name() {
	return "RandomizeBBByRamaPrePro";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
RandomizeBBByRamaPreProCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RandomizeBBByRamaPrePro );
}

std::string
RandomizeBBByRamaPreProCreator::keyname() const
{
	return RandomizeBBByRamaPrePro::mover_name();
}

void RandomizeBBByRamaPreProCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RandomizeBBByRamaPrePro::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::ostream &
operator<<( std::ostream & os, RandomizeBBByRamaPrePro const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //backbone_moves
