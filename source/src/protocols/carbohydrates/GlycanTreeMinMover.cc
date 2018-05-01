// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeMinMover.cc
/// @brief A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/carbohydrates/GlycanTreeMinMover.hh>
#include <protocols/carbohydrates/GlycanTreeMinMoverCreator.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/minimization_packing/MinMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/kinematics/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/symmetry/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.carbohydrates.GlycanTreeMinMover" );

namespace protocols {
namespace carbohydrates {

using namespace core::kinematics;
using namespace core::scoring;
using namespace protocols::minimization_packing;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycanTreeMinMover::GlycanTreeMinMover():
	protocols::moves::Mover("GlycanTreeMinMover")
{

}

GlycanTreeMinMover::GlycanTreeMinMover(
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool min_rings /* False */,
	bool min_bb /* False */,
	bool min_chi /* False */):
	protocols::moves::Mover("GlycanTreeMinMover"),
	min_rings_(min_rings),
	min_bb_(min_bb),
	min_chi_(min_chi)

{
	selector_ = selector->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanTreeMinMover::GlycanTreeMinMover( GlycanTreeMinMover const & src ):
	protocols::moves::Mover( src ),
	min_rings_(src.min_rings_),
	min_bb_(src.min_bb_),
	min_chi_(src.min_chi_)

{
	if ( src.selector_ ) {
		selector_ = src.selector_->clone();
	}

	if ( src.min_mover_ ) {
		set_minmover( src.min_mover_ );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanTreeMinMover::~GlycanTreeMinMover()= default;


///@brief Minimize Rings? Default False
void
GlycanTreeMinMover::set_min_rings( bool min_rings ){
	min_rings_ = min_rings;
}

///@brief Minimize Chi? Default True
void
GlycanTreeMinMover::set_min_chi( bool min_chi ){
	min_chi_ = min_chi;
}

///@brief Minimize BB? Default True
void
GlycanTreeMinMover::set_min_bb( bool min_bb ){
	min_bb_ = min_bb;
}

///@brief Set a pre-configured MinMover for this class.
///  Will OVERRIDE movemap settings.
void
GlycanTreeMinMover::set_minmover( protocols::minimization_packing::MinMoverCOP min_mover){
	min_mover_ = protocols::minimization_packing::MinMoverOP( new protocols::minimization_packing::MinMover( *min_mover));
}



////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
GlycanTreeMinMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanTreeMinMover::clone() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover( *this ) );
}

/// @brief Get the name of the Mover
std::string
GlycanTreeMinMover::get_name() const
{
	return GlycanTreeMinMover::class_name();
}

std::string
GlycanTreeMinMover::class_name()
{
	return "GlycanTreeMinMover";
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycanTreeMinMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GlycanTreeMinMover );
}

std::string
GlycanTreeMinMoverCreator::keyname() const
{
	return GlycanTreeMinMover::class_name();
}

void
GlycanTreeMinMover::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	runtime_assert( selector != nullptr );
	selector_ = selector->clone();

}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////


/// @brief Apply the mover
void
GlycanTreeMinMover::apply( core::pose::Pose& pose){

	using namespace core::select::residue_selector;
	using namespace core::scoring;

	GlycanResidueSelectorOP branch_selector = GlycanResidueSelectorOP( new GlycanResidueSelector() );
	AndResidueSelectorOP and_selector = AndResidueSelectorOP( new AndResidueSelector());


	branch_selector->set_include_root( true ); //Since we wont really have the ASN root, and we want the possibility to sample the whole glycan.

	if ( ! selector_ ) {
		selector_ = GlycanResidueSelectorOP( new GlycanResidueSelector());
	}

	utility::vector1< core::Size > mm_residues = selection_positions( selector_->apply(pose));

	//We need to set reasonable defaults here!




	//Randomly select a residue from the movemap.
	core::Size index = numeric::random::rg().random_range( 1, mm_residues.size() );
	core::Size resnum = mm_residues[ index ];

	TR << "Minimizing from carbohydrate root: " << resnum << std::endl;


	//Get downstream Tree.
	branch_selector->set_select_from_branch_residue( resnum );

	and_selector->add_residue_selector( selector_ );
	and_selector->add_residue_selector( branch_selector );


	utility::vector1< bool > subset = and_selector->apply(pose);
	core::Size min_residue_n = count_selected( subset );

	MoveMapOP mm = core::pose::carbohydrates::create_glycan_movemap_from_residue_selector(
		pose,
		and_selector,
		min_chi_,
		min_rings_,
		min_bb_
	);

	if ( ! min_mover_ ) {
		//Make Symmetric or Non-symmetric versions of the MinMover.
		ScoreFunctionCOP scorefxn = core::scoring::get_score_function();

		min_mover_ = MinMoverOP( new MinMover( mm->clone(), scorefxn, "dfpmin_armijo_nonmonotone", 0.01, true /* use_nblist*/ ) );
	} else {

		min_mover_->set_movemap( mm );
	}

	//Minimize
	TR << "Minimizing " << min_residue_n << " residues " << std::endl;
	min_mover_->apply(pose);
}

void
GlycanTreeMinMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	using namespace core::select::residue_selector;

	min_rings_ = tag->getOption< bool >("min_rings", min_rings_);
	min_bb_ = tag->getOption< bool >("min_bb", min_bb_);
	min_chi_ = tag->getOption< bool >("min_chi", min_chi_);


	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	}
}

void
GlycanTreeMinMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("min_rings", xsct_rosetta_bool, "Minimize Ring Torsions?", "false")
		+ XMLSchemaAttribute::attribute_w_default("min_bb", xsct_rosetta_bool, "Minimize the Backbone?", "true")

		+ XMLSchemaAttribute::attribute_w_default("min_chi", xsct_rosetta_bool, "Minimize the Chi angles?", "true");

	rosetta_scripts::attributes_for_parse_residue_selector( attlist );

	XMLSchemaSimpleSubelementList subelements;
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, "GlycanTreeMinMover",

		" A class that selects the downstream branch from residues in a movemap/selector, and minimizes those residues if on in the primary glycan movemap. Multiple Applies randomly select a different residue in the movemap/selector", attlist, subelements );

}

void GlycanTreeMinMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanTreeMinMover::provide_xml_schema( xsd );
}




} //protocols
} //carbohydrates

