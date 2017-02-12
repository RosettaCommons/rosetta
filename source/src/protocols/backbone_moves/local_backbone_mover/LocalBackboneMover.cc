// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/LocalBackboneMover.cc
/// @brief LocalBackboneMover moves a stretch of backbone locally.
/// @author xingjiepan (xingjiepan@gmail.com)

// Unit headers
#include <protocols/backbone_moves/local_backbone_mover/LocalBackboneMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/LocalBackboneMoverCreator.hh>
#include <protocols/backbone_moves/local_backbone_mover/GapCloser.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.hh>
#include <protocols/backbone_moves/local_backbone_mover/LocalBackboneMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.hh>

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.hh>

// Utility headers
#include <utility/exit.hh>

// Core headers
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// STD
#include <string>


static THREAD_LOCAL basic::Tracer TR( "protocols.backbone_moves.local_backbone_mover.LocalBackboneMover" );

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
LocalBackboneMover::LocalBackboneMover():
	protocols::moves::Mover( LocalBackboneMover::mover_name() )
{
	gap_closer_ = GapCloserOP(new GapCloser);
	set_max_trial_num(10);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
LocalBackboneMover::~LocalBackboneMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////

/// @brief Apply the mover
void
LocalBackboneMover::apply( core::pose::Pose& pose){
	for(Size i=1; i<= max_trial_num_; ++i){
	
		set_free_peptide(pose, pivot1_, pivot2_);

		// Move the free peptide
		for(free_peptide_movers::FreePeptideMoverOP f_mover : free_peptide_movers_){
			f_mover->apply(*free_peptide_);
		}

		// Close gaps
		gap_closer_->solve_gaps(*free_peptide_);

		// Apply changes
		if(gap_closer_->gap_solved()){
			gap_closer_->apply_closure(pose, *free_peptide_);
			return;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
LocalBackboneMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
LocalBackboneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	using namespace free_peptide_movers;

	// Set the pivots
	
	pivot1_ = tag->getOption<Size>("pivot1");
	pivot2_ = tag->getOption<Size>("pivot2");

	// Set the free peptide movers

	clear_free_peptide_movers();

	std::string config = tag->getOption<std::string>("move_type", "translate");
	if("translate" == config){
		add_free_peptide_mover(FreePeptideMoverOP( new TranslationFreePeptideMover(0.5) ));
	}else if("rotate" == config){
		add_free_peptide_mover(FreePeptideMoverOP( new LongAxisRotationFreePeptideMover(1, true) ));
	}else if("shear" == config){
		add_free_peptide_mover(FreePeptideMoverOP( new ShearFreePeptideMover(pivot1_ + 1, 120, true) ));
	}else if("circular_permute" == config){
		add_free_peptide_mover(FreePeptideMoverOP( new CircularPermuteFreePeptideMover(1) ));
	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
LocalBackboneMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new LocalBackboneMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
LocalBackboneMover::clone() const
{
	return protocols::moves::MoverOP( new LocalBackboneMover( *this ) );
}

std::string LocalBackboneMover::get_name() const {
	return mover_name();
}

std::string LocalBackboneMover::mover_name() {
	return "LocalBackboneMover";
}

void LocalBackboneMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default(
		"move_type", xs_string,
		"Set the type of backbone movement. Available settings are 'translate', 'rotate', 'shear' and 'circular_permute'. 'translate' is the default setting.",
		"translate");

	attlist + XMLSchemaAttribute::required_attribute(
		"pivot1", xsct_non_negative_integer,
		"First pivot residue.");

	attlist + XMLSchemaAttribute::required_attribute(
		"pivot2", xsct_non_negative_integer,
		"Second pivot residue.");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Move the backbone between two selected pivots. Conformation changes are restricted two reisidues between pivot1 - 1 and pivot2 + 1.", attlist );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
LocalBackboneMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new LocalBackboneMover );
}

std::string
LocalBackboneMoverCreator::keyname() const
{
	return LocalBackboneMover::mover_name();
}

void LocalBackboneMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LocalBackboneMover::provide_xml_schema( xsd );
}

void LocalBackboneMover::set_free_peptide(core::pose::Pose &pose, Size pivot1, Size pivot2){
	if(pivot1 < 3){
		utility_exit_with_message("Residue pivot1 - 2 must exist!");	
	}
	if(pivot2 > pose.size() - 2){
		utility_exit_with_message("Residue pivot2 + 2 must exist!");	
	}
	if(pivot1 + 3 > pivot2){
		utility_exit_with_message("Residue there must be at least 2 residues between pivot 1 and pivot2!");	
	}

	free_peptide_ = FreePeptideOP(new FreePeptide(pose, pivot1, pivot2));
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


	std::ostream &
	operator<<( std::ostream & os, LocalBackboneMover const & mover )
	{
		mover.show(os);
		return os;
	}

} //protocols
} //backbone_moves
} //local_backbone_mover
