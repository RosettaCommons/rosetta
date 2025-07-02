// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon , adapted from the ResfileReader code
/// by Steven Lewis (smlewi@gmail.com) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/FinalMinimizerCreator.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/InterfaceBuilder.fwd.hh>
#include <protocols/ligand_docking/MinimizeBackbone.hh>

//Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/CitationManager.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.ligand_options.FinalMinimizer" );




FinalMinimizer::FinalMinimizer():
	Mover("FinalMinimizer"),
	score_fxn_(/* NULL */),
	movemap_builder_(/* NULL */),
	remove_bb_constraints_(false)
{}

FinalMinimizer::FinalMinimizer(
	core::scoring::ScoreFunctionOP score_fxn,
	MoveMapBuilderOP movemap_builder,
	bool remove_constraints
):
	Mover("FinalMinimizer"),
	score_fxn_(score_fxn),
	movemap_builder_(movemap_builder),
	remove_bb_constraints_( remove_constraints )
{}

FinalMinimizer::FinalMinimizer(FinalMinimizer const & /*that*/) = default;

FinalMinimizer::~FinalMinimizer() = default;

protocols::moves::MoverOP FinalMinimizer::clone() const {
	return utility::pointer::make_shared< FinalMinimizer >( *this );
}

protocols::moves::MoverOP FinalMinimizer::fresh_instance() const {
	return utility::pointer::make_shared< FinalMinimizer >();
}


//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
FinalMinimizer::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->getName() != "FinalMinimizer" ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'FinalMinimizer' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get_ptr<core::scoring::ScoreFunction>( "scorefxns", scorefxn_name);

	/// MoveMapBuilder///
	if ( ! tag->hasOption("movemap_builder") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'FinalMinimizer' requires 'movemap_builder' tag");
	std::string movemap_builder_name= tag->getOption<std::string>("movemap_builder");
	movemap_builder_= datamap.get_ptr<protocols::ligand_docking::MoveMapBuilder>( "movemap_builders", movemap_builder_name);

	remove_bb_constraints_ = tag->getOption<bool>("remove_constraints", false);
	cartesian_ = tag->getOption<bool>("cartesian", false);

}

void
FinalMinimizer::apply( core::pose::Pose & pose ){
	debug_assert(movemap_builder_);

	TR << "Energy prior to minimizing: " << (*score_fxn_)(pose) << std::endl;
	if ( movemap_builder_->minimize_backbone() ) {
		TR.Debug << "Setting up FoldTree for backbone minimization" << std::endl;
		TR.Debug << "FoldTree Before: " << pose.fold_tree() << std::endl;
		core::kinematics::FoldTree fold_tree_copy( pose.fold_tree() );
		InterfaceBuilderOP bb_interface_builder= movemap_builder_->get_bb_interface_builder();
		MinimizeBackbone backbone_foldtree_setup(bb_interface_builder);
		backbone_foldtree_setup.apply(pose);
		TR.Debug << "FoldTree Reordered: " << pose.fold_tree() << std::endl;

		protocols::minimization_packing::MinMoverOP const dfpMinTightTol = get_final_min_mover(pose, true);
		dfpMinTightTol->min_options()->nblist_auto_update(true);
		dfpMinTightTol->apply(pose);
		pose.fold_tree(fold_tree_copy);
		backbone_foldtree_setup.remove_cutpoints(pose); //We redid the foldtree - we should redo the cutpoints.
		if ( remove_bb_constraints_ ) {
			backbone_foldtree_setup.remove_constraints(pose);
		}
		TR.Debug << "FoldTree Restored: " << pose.fold_tree() << std::endl;
	} else {
		protocols::minimization_packing::MinMoverOP const dfpMinTightTol = get_final_min_mover(pose, false);
		dfpMinTightTol->min_options()->nblist_auto_update(true);
		dfpMinTightTol->apply(pose);
	}
	TR << "Energy after minimizing: " << (*score_fxn_)(pose) << std::endl;

}

protocols::minimization_packing::MinMoverOP const
FinalMinimizer::get_final_min_mover(core::pose::Pose const & pose, bool backbone) const{
	std::string min_type= "lbfgs_armijo_nonmonotone_atol";
	core::Real tolerance= 0.02;
	bool use_nb_list=true;
	core::kinematics::MoveMapOP movemap= movemap_builder_->build(pose);
	core::scoring::ScoreFunctionCOP scorefxn( score_fxn_ );
	if ( backbone && scorefxn->get_weight( core::scoring::chainbreak ) == 0 ) {
		TR.Warning << "FinalMinimizer with backbone minimization used without chainbreak scoreterm - setting chainbreak to 1" << std::endl;
		core::scoring::ScoreFunctionOP new_scorefxn( scorefxn->clone() );
		new_scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
		scorefxn = new_scorefxn;
	}
	if ( cartesian_ && scorefxn->get_weight( core::scoring::cart_bonded ) == 0 ) {
		min_type = "lbfgs_armijo_nonmonotone"; // Recommended type for cartesian minimization.
		TR.Warning << "FinalMinimizer with Cartesian min used without cart_bonded scoreterm - setting cart_bonded to 1" << std::endl;
		core::scoring::ScoreFunctionOP new_scorefxn( scorefxn->clone() );
		new_scorefxn->set_weight( core::scoring::cart_bonded, 1.0 );
		new_scorefxn->set_weight( core::scoring::pro_close, 0.0 );
		scorefxn = new_scorefxn;
	}
	movemap->show(TR.Debug, pose.size());
	TR.Debug << std::endl;
	auto min_mover = utility::pointer::make_shared< protocols::minimization_packing::MinMover >(movemap, scorefxn, min_type, tolerance, use_nb_list);
	min_mover->cartesian(cartesian_);
	return min_mover;
}

std::string FinalMinimizer::get_name() const {
	return mover_name();
}

std::string FinalMinimizer::mover_name() {
	return "FinalMinimizer";
}

void FinalMinimizer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("scorefxn", xs_string, "Used scorefunction. Required.")
		+ XMLSchemaAttribute::required_attribute("movemap_builder", xs_string, "Name of a previously defined MoveMapBuilder. Required.")
		+ XMLSchemaAttribute::attribute_w_default("remove_constraints", xsct_rosetta_bool, "Remove any added constraints after minimization.", "false" )
		+ XMLSchemaAttribute::attribute_w_default("cartesian", xsct_rosetta_bool, "Use Cartesian minimization.", "false" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Do a gradient based minimization of the final docked pose.", attlist );
}

std::string FinalMinimizerCreator::keyname() const {
	return FinalMinimizer::mover_name();
}

protocols::moves::MoverOP
FinalMinimizerCreator::create_mover() const {
	return utility::pointer::make_shared< FinalMinimizer >();
}

void FinalMinimizerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FinalMinimizer::provide_xml_schema( xsd );
}

/// @brief Provide the citation.
void
FinalMinimizer::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationCollectionOP cc(
		utility::pointer::make_shared< basic::citation_manager::CitationCollection >(
		"FinalMinimizer", basic::citation_manager::CitedModuleType::Mover
		)
	);
	cc->add_citation( basic::citation_manager::CitationManager::get_instance()->get_citation_by_doi( "10.1007/978-1-61779-465-0_10" ) );

	citations.add( cc );
	//TODO ADD CITATION FOR MOVEMAP BUILDER AND SCOREFUNCTION!!!
}

} //namespace ligand_docking
} //namespace protocols
