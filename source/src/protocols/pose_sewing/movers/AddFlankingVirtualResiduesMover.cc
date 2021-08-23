// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.cc
/// @brief adds virtual residues to either side of a given Pose
/// @author frankdt (frankdt@email.unc.edu)
/// @author Jared Adolf-Bryfogle

// Unit headers
#include <protocols/pose_sewing/movers/AddFlankingVirtualResiduesMover.hh>
#include <protocols/pose_sewing/movers/AddFlankingVirtualResiduesMoverCreator.hh>
#include <protocols/pose_sewing/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/scoring/methods/util.hh>
#include <core/scoring/dssp/Dssp.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.pose_sewing.movers.AddFlankingVirtualResiduesMover" );

namespace protocols {
namespace pose_sewing {
namespace movers {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AddFlankingVirtualResiduesMover::AddFlankingVirtualResiduesMover():
	protocols::moves::Mover( AddFlankingVirtualResiduesMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
AddFlankingVirtualResiduesMover::AddFlankingVirtualResiduesMover( AddFlankingVirtualResiduesMover const & src ):
	protocols::moves::Mover( src )
{
	N_term_length_ = src.get_N_term_length();
	C_term_length_ = src.get_C_term_length();
	vital_selector_ = src.get_vital_selector();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
AddFlankingVirtualResiduesMover::~AddFlankingVirtualResiduesMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
AddFlankingVirtualResiduesMover::apply( core::pose::Pose& pose){

	add_flanking_virtual_residues(pose);

}
void
AddFlankingVirtualResiduesMover::add_flanking_virtual_residues( core::pose::Pose& pose){

	core::pose::Pose original_pose = pose;
	if ( remove_pre_pose_ ) {
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			pose.pdb_info()->clear_reslabel(i, "PRE_POSE_SEWING_START");
			pose.pdb_info()->clear_reslabel(i, "PRE_POSE_SEWING_END");
		}
	}
	//Make refpose for use with labels.
	std::string refpose_name = "starting_model_refpose";
	pose.reference_pose_from_current(refpose_name, true);
	if ( chain_to_modify_ > pose.num_chains() ) {
		utility_exit_with_message("AddFlankingVirtResidues: Chain to modify is greater than the number of chains in the pose!");
	}
	if ( N_term_length_ == 0 and C_term_length_ == 0 ) {
		utility_exit_with_message("AddFlankingVirtResidues: Please select N and/or C term residues to add virts for!");
	}
	core::chemical::ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose();
	core::chemical::ResidueTypeCOP res_type = rsd_set->get_representative_type_name1( 'A');
	if ( N_term_length_ ) {
		core::Size start = pose.chain_begin(chain_to_modify_);
		for ( core::Size i = 1; i <= N_term_length_; ++i ) {
			core::Size current_resnum = start ;
			core::pose::remove_lower_terminus_type_from_pose_residue(pose, current_resnum);
			core::conformation::ResidueOP res = core::conformation::ResidueFactory::create_residue(*res_type);
			pose.prepend_polymer_residue_before_seqpos(*res, current_resnum, true);
			pose.pdb_info()->chain(current_resnum,pose.pdb_info()->chain(current_resnum+1));
			pose.pdb_info()->number(current_resnum,pose.pdb_info()->number(current_resnum+1)-1);
			pose.real_to_virtual(current_resnum );
			if ( i == N_term_length_ ) {
				pose.pdb_info()->add_reslabel(current_resnum , "PRE_POSE_SEWING_START");
			}
		}
	}
	if ( C_term_length_ ) {
		core::Size start = pose.chain_end(chain_to_modify_);
		for ( core::Size i = 1; i <= C_term_length_; ++i ) {
			core::Size current_resnum = start + i -1;
			core::pose::remove_upper_terminus_type_from_pose_residue(pose, current_resnum);
			core::conformation::ResidueOP res = core::conformation::ResidueFactory::create_residue(*res_type);
			pose.append_polymer_residue_after_seqpos(*res, current_resnum, true);
			pose.pdb_info()->chain(current_resnum+1,pose.pdb_info()->chain(current_resnum));
			pose.pdb_info()->number(current_resnum+1,pose.pdb_info()->number(current_resnum)+1);
			pose.real_to_virtual(current_resnum+1);
			if ( i == C_term_length_ ) {
				pose.pdb_info()->add_reslabel(current_resnum + 1, "PRE_POSE_SEWING_END");
			}
		}
	}
	assert(pose.size() == original_pose.size() + C_term_length_ + N_term_length_);
	if ( vital_selector_ ) {
		core::select::residue_selector::ResidueSubset subset = vital_selector_->apply(original_pose);
		for ( core::Size i = 1; i <= subset.size(); ++i ) {
			if ( subset[i] ) {
				core::Size new_resnum = pose.corresponding_residue_in_current(i, refpose_name);
				pose.pdb_info()->add_reslabel(new_resnum, "VITAL_RESIDUE");
			}
		}
	}
}

core::Size
AddFlankingVirtualResiduesMover::get_N_term_length() const {
	return N_term_length_;
}

core::Size
AddFlankingVirtualResiduesMover::get_C_term_length() const {
	return C_term_length_;
}

core::select::residue_selector::ResidueSelectorCOP
AddFlankingVirtualResiduesMover::get_vital_selector() const {
	return vital_selector_;
}

core::Size
AddFlankingVirtualResiduesMover::get_chain_to_modify() const {
	return chain_to_modify_;
}

bool
AddFlankingVirtualResiduesMover::get_remove_pre_pose() const {
	return remove_pre_pose_;
}

void
AddFlankingVirtualResiduesMover::set_N_term_length(core::Size N_length){
	N_term_length_ = N_length;
}

void
AddFlankingVirtualResiduesMover::set_C_term_length(core::Size C_length){
	C_term_length_ = C_length;
}


void
AddFlankingVirtualResiduesMover::set_vital_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	vital_selector_ = selector;
}

void
AddFlankingVirtualResiduesMover::set_chain_to_modify(core::Size chain){
	chain_to_modify_ = chain;
}

void
AddFlankingVirtualResiduesMover::set_remove_pre_pose(bool in){
	remove_pre_pose_ = in;
}
////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
AddFlankingVirtualResiduesMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AddFlankingVirtualResiduesMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap
)
{
	if ( tag->hasOption("N_term_residues") ) {
		N_term_length_ = tag->getOption< core::Size >("N_term_residues");
	}
	if ( tag->hasOption("C_term_residues") ) {
		C_term_length_ = tag->getOption< core::Size >("C_term_residues");
	}
	if ( tag->hasOption("chain_to_modify") ) {
		chain_to_modify_ = tag->getOption< core::Size >("chain_to_modify");
	}
	if ( tag->hasOption( "vital_selector" ) ) {
		vital_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "vital_selector" ), datamap );
	}
	remove_pre_pose_ = tag->getOption("remove_pre_pose", remove_pre_pose_);


}
void AddFlankingVirtualResiduesMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	attlist + XMLSchemaAttribute::attribute_w_default( "N_term_residues", xsct_non_negative_integer,  "How many residues should we add to the N terminus?", "0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "C_term_residues", xsct_non_negative_integer,  "How many residues should we add to the C terminus?", "0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "chain_to_modify", xsct_non_negative_integer,  "Number of the chain to be prepared", "1" );
	attlist + XMLSchemaAttribute::attribute_w_default( "remove_pre_pose", xsct_non_negative_integer,  "Remove any current pre_pose res labels", "0" );

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "vital_selector",  "vital residue selector. These can't move or mutate." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add virtual residues (alanines) to either termini of a given chain. Used by SewAnythingAddMover to guide design to either end.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AddFlankingVirtualResiduesMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AddFlankingVirtualResiduesMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AddFlankingVirtualResiduesMover::clone() const
{
	return protocols::moves::MoverOP( new AddFlankingVirtualResiduesMover( *this ) );
}

std::string AddFlankingVirtualResiduesMover::get_name() const {
	return mover_name();
}

std::string AddFlankingVirtualResiduesMover::mover_name() {
	return "AddFlankingVirtualResiduesMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
AddFlankingVirtualResiduesMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddFlankingVirtualResiduesMover );
}

std::string
AddFlankingVirtualResiduesMoverCreator::keyname() const
{
	return AddFlankingVirtualResiduesMover::mover_name();
}

void AddFlankingVirtualResiduesMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddFlankingVirtualResiduesMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, AddFlankingVirtualResiduesMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //pose_sewing
} //movers
