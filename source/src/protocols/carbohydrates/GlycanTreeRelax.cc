// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeRelax.cc
/// @brief A protocol for optimizing glycan trees using GlycanRelax from the base of the tree out to the leaves.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Raemisch (raemisch@scripps.edu)

// Unit headers
#include <protocols/carbohydrates/GlycanTreeRelax.hh>
#include <protocols/carbohydrates/GlycanTreeRelaxCreator.hh>
#include <protocols/carbohydrates/GlycanRelaxMover.hh>
#include <protocols/carbohydrates/util.hh>

#include <protocols/simple_moves/ConvertRealToVirtualMover.hh>
#include <protocols/simple_moves/ConvertVirtualToRealMover.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/GlycanLayerSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/kinematics/MoveMap.hh>

// Protocl headers
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.carbohydrates.GlycanTreeRelax" );

namespace protocols {
namespace carbohydrates {
using namespace protocols::simple_moves;
using namespace core::scoring;
using namespace core::pose;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycanTreeRelax::GlycanTreeRelax():
	protocols::moves::Mover( GlycanTreeRelax::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanTreeRelax::GlycanTreeRelax( GlycanTreeRelax const & src ):
	protocols::moves::Mover( src ),
	layer_size_(src.layer_size_),
	window_size_(src.window_size_),
	rounds_(src.rounds_),
	completed_quenches_(src.completed_quenches_),
	trees_to_model_(src.trees_to_model_),
	glycan_relax_rounds_(src.glycan_relax_rounds_),
	refine_(src.refine_),
	quench_mode_(src.quench_mode_),
	final_min_pack_min_(src.final_min_pack_min_)

{
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanTreeRelax::~GlycanTreeRelax(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycanTreeRelax::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
GlycanTreeRelax::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	layer_size_ = tag->getOption< core::Size >( "layer_size", layer_size_);
	window_size_ = tag->getOption< core::Size >("window_size", window_size_);
	rounds_ = tag->getOption< core::Size >("rounds", rounds_);

	refine_ = tag->getOption< bool >("refine", refine_);
	quench_mode_ = tag->getOption< bool >("quench_mode", quench_mode_);
	final_min_pack_min_ = tag->getOption< bool >("final_min_pack_min", final_min_pack_min_);
	glycan_relax_rounds_ = tag->getOption< core::Size >("glycan_relax_rounds", glycan_relax_rounds_);

	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	}
	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
GlycanTreeRelax::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanTreeRelax );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanTreeRelax::clone() const
{
	return protocols::moves::MoverOP( new GlycanTreeRelax( *this ) );
}

std::string GlycanTreeRelax::get_name() const {
	return mover_name();
}

std::string GlycanTreeRelax::mover_name() {
	return "GlycanTreeRelax";
}

void GlycanTreeRelax::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
	attlist + XMLSchemaAttribute("layer_size", xsct_non_negative_integer,
		"Brief: Set the layer size we will be using.  A layer is a set of glycan residues that we will be optimizing.\n"
		"  We work our way through the layers, while the rest of the residues are virtual (not scored).\n"
		" \n"
		"Details: \n"
		" \n"
		"  The distance that make up a layer.  If we have a distance of 2,\n"
		"  we first model all glycans that are equal to or less than 2 residue distance to the root.\n"
		"  We then slide this layer up.  So we take all residues that have a distance between 3 and 1, and so on.\n")


		+ XMLSchemaAttribute("window_size", xsct_non_negative_integer,
		"\tBrief: Set the window size.  This is the overlap of the layers during modeling. \n"
		" \n"
		"  Details: \n"
		"  \n"
		"  The layers are slid down throught the tree of the glycan.  The window size represents the overlap in the layers.\n"
		"  A window size of 1, means that the last residue (or residues of layer 1) from the last modeling effort, will be used again as \n"
		"  part of the next layer.  A window size of 0, means that no residues will be re-modeled. \n"
		"  Typically, we would want at least a window size of 1.\n")

		+ XMLSchemaAttribute("rounds", xsct_non_negative_integer, "Set the number of rounds.  We start with the forward direction, from the "
		"root/start out to the glycan tips.  Next round we go backward, then forward.  The number of rounds corresponds to how many back "
		"and forth directions we go")

		+ XMLSchemaAttribute::attribute_w_default("refine", xsct_rosetta_bool, "Do a refinement instead of a denovo model", "false")
		+ XMLSchemaAttribute::attribute_w_default("quench_mode", xsct_rosetta_bool, "Do quench mode for each glycan tree?", "false")
		+ XMLSchemaAttribute::attribute_w_default("final_min_pack_min", xsct_rosetta_bool, "Do a final set of cycles of min/pack", "true")
		+ XMLSchemaAttribute::attribute_w_default("glycan_relax_rounds", xsct_non_negative_integer, "Round Number for the internal GlycanRelax.  Default is the default of GlycanRelax.", "25");


	std::string documentation =
		"Brief: A protocol for optimizing glycan trees using GlycanRelax from the base of the tree out to the leaves.\n"
		"Details: Works by making all other residues virtual except the ones it is working on (current Layer).\n"
		"A virtual residue is not scored.\n"
		"It will start at the first glycan residues, and then move out to the edges.\n"
		"\n"
		"GENERAL ALGORITHM\n"
		"\n"
		"We start at the roots, and make all other glycan residues virtual.\n"
		"We first model towards the leaves and this is considered the forward direction.\n"
		"GlycanRelax is used for the actual modeling, we only model a layer at a time, until we reach the tips.\n"
		"If more than one round is set, the protocol will move backwards on the next round, from the leafs to the roots.\n"
		"A third round will involve relaxation again in the forward direction.\n"
		"So we go forward, back, forward, etc. for how ever many rounds you set.\n"
		"\n"
		"QUECHING\n"
		"\n"
		"By default, we model all glycans simultaneously. First, all glycan roots (the start of the tree), and slowly unvirtualize \n"
		"all glycan residues, while only modeling each layer. \n"
		"Alternatively, we can choose a particular glycan tree, run the algorithm, and then choose another glycan tree randomly until all \n"
		"glycan trees have been optimized. \n"
		"Here, we call this quenching. \n"
		"\n"
		"GLYCAN LAYERS \n"
		"\n"
		"Draw a tree on a paper.  We start with the beginning N residues, and work our way out towards the leaves. \n"
		"Layers are defined by the glycan residue distance to the rooot.  This enables branching residues to be considered the same \n"
		"layer conceptually and computationally, and allows them to be modeled together. \n"
		"\n"
		"--LAYER SIZE-- \n"
		"\n"
		"The distance that make up a layer.  If we have a distance of 2, \n"
		"we first model all glycans that are equal to or less than 2 residue distance to the root. \n"
		"We then slide this layer up.  So we take all residues that have a distance between 3 and 1, and so on. \n"
		"\n"
		"--WINDOW SIZE-- \n"
		"\n"
		"The layers are slid down throught the tree of the glycan.  The window size represents the overlap in the layers. \n"
		"A window size of 1, means that the last residue (or residues of layer 1) from the last modeling effort, will be used again as \n"
		"part of the next layer.  A window size of 0, means that no residues will be re-modeled. \n"
		"Typically, we would want at least a window size of 1. \n";

	rosetta_scripts::attributes_for_parse_residue_selector( attlist );
	rosetta_scripts::attributes_for_get_score_function_name( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), documentation, attlist );
}


void
GlycanTreeRelax::set_layer_size( const core::Size layer_size ){
	layer_size_ = layer_size;
}

void
GlycanTreeRelax::set_window_size( core::Size const window_size ){
	window_size_ = window_size;
}

void
GlycanTreeRelax::set_refine(const bool refine){
	refine_ = refine;
}

void
GlycanTreeRelax::set_final_min_pack_min( const bool minpackmin ){
	final_min_pack_min_ = minpackmin;
}

void
GlycanTreeRelax::set_quench_mode( bool quench_mode ){
	quench_mode_ = quench_mode;
}

void
GlycanTreeRelax::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

void
GlycanTreeRelax::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector->clone();
}

void
GlycanTreeRelax::set_rounds(const core::Size rounds){
	rounds_ = rounds;
}

///@brief Override Glycan Relax rounds.
void
GlycanTreeRelax::set_glycan_relax_rounds( core::Size glycan_relax_rounds){
	glycan_relax_rounds_ = glycan_relax_rounds;
}

bool
GlycanTreeRelax::is_quenched() const {
	if ( ! quench_mode_ && completed_quenches_ == 1 ) {
		return true;
	} else if ( ! quench_mode_ && completed_quenches_ == 0 ) {
		return false;
	} else if ( quench_mode_ && completed_quenches_ == trees_to_model_ ) {
		return true;
	} else {
		return false;
	}

}

/// @brief Apply the mover
void
GlycanTreeRelax::apply( core::pose::Pose & pose){
	using namespace core::select::residue_selector;
	using namespace core::kinematics;
	using namespace core::pack::task;

	debug_assert( layer_size_ != window_size_ );

	//Setup everything we need.

	// Setup Selectors
	GlycanResidueSelectorOP glycan_selector = GlycanResidueSelectorOP( new GlycanResidueSelector() );
	GlycanLayerSelectorOP modeling_layer_selector = GlycanLayerSelectorOP( new GlycanLayerSelector() );
	GlycanLayerSelectorOP virtual_layer_selector = GlycanLayerSelectorOP( new GlycanLayerSelector() );
	ReturnResidueSubsetSelectorOP store_subset = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector() );
	ReturnResidueSubsetSelectorOP store_to_de_virt_subset =  ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector() );
	AndResidueSelectorOP and_selector = AndResidueSelectorOP( new AndResidueSelector() );
	AndResidueSelectorOP and_selector_virt = AndResidueSelectorOP( new AndResidueSelector() );




	ResidueSubset starting_subset;
	//ResidueSubset current_subset;

	if ( selector_ && ! quench_mode_ ) {
		starting_subset = selector_->apply( pose );

	} else {
		starting_subset = glycan_selector->apply( pose );
	}

	ResidueSubset min_subset = starting_subset;

	if ( ! scorefxn_ ) {
		scorefxn_ = get_score_function();
	}

	core::Real starting_score = scorefxn_->score( pose );
	TR << "Starting Score: " << starting_score << std::endl;


	// Setup GlycanRelax
	GlycanRelaxMover glycan_relax = GlycanRelaxMover();
	if ( scorefxn_ ) {
		glycan_relax.set_scorefunction( scorefxn_ );
	}
	glycan_relax.set_refine( refine_ );

	if ( glycan_relax_rounds_ != 0 ) {
		glycan_relax.set_rounds(glycan_relax_rounds_);
	}

	ConvertRealToVirtualMover real_to_virt = ConvertRealToVirtualMover();
	ConvertVirtualToRealMover virt_to_real = ConvertVirtualToRealMover();

	QuenchDirection current_direction = forward;

	utility::vector1< core::Size > all_trees = pose.glycan_tree_set()->get_start_points();
	utility::vector1< core::Size > tree_subset;

	//Setup for quench mode
	if ( quench_mode_ ) {
		if ( selector_ ) {
			utility::vector1< bool > tree_subset = selector_->apply( pose );
			for ( core::Size tree_start : all_trees ) {
				if ( tree_subset[ tree_start ] ) {
					tree_subset.push_back( tree_start );
				}
			}
		} else {
			tree_subset = all_trees;
		}
		trees_to_model_ = tree_subset.size();
	}



	while ( ! is_quenched() ) {

		if ( quench_mode_ ) {
			//Need to set total glycan trees

			//set current subset down to a single tree randomly.  Remove one of the trees.
			TR << "Quench mode: " << tree_subset.size() << " left to model " << std::endl;
			core::Size tree_index = numeric::random::rg().random_range( 1, tree_subset.size() );
			core::Size tree_start = all_trees[ tree_index ];

			glycan_selector->set_include_root(true);
			glycan_selector->set_select_from_branch_residue( tree_start );
			starting_subset = glycan_selector->apply( pose );

			//Remove the tree from the list to model.
			tree_subset.pop( tree_start );

		}
		for ( core::Size i = 1; i <= rounds_; ++i ) {
			TR << "starting round " << i << std::endl;

			if ( i % 2 == 0 ) {
				current_direction = backward;
			} else {
				current_direction = forward;
			}
			if ( i > 1 ) {

				///Don't randomize if we already have something.  Just try to make it better.
				glycan_relax.set_refine( true );
			}

			switch ( current_direction ) {
			case forward :
				{
				TR << "Going in the forward direction " << std::endl;

				bool first_modeling = true;
				core::Size current_start = 0;
				core::Size current_end = layer_size_ - 1;

				core::Size max_end = pose.glycan_tree_set()->get_largest_glycan_tree_layer();
				TR << "Modeling up to max end: " << max_end << std::endl;
				ResidueSubset pre_virtualized( pose.total_residue(), false );
				ResidueSubset needed_virtualization( pose.total_residue(), false );



				while ( current_start <= max_end ) {


					modeling_layer_selector->set_layer( current_start, current_end );
					//ResidueSubset modeling_layer = modeling_layer_selector->apply( pose );

					virtual_layer_selector->set_layer_as_greater_than_or_equal_to( current_end+1 );
					//ResidueSubset virtual_layer = virtual_layer_selector->apply( pose );

					//Combine modeling layer with whatver the current_subset is.
					store_subset->set_residue_subset( starting_subset );
					and_selector->clear();
					and_selector_virt->clear();

					and_selector->add_residue_selector( store_subset );
					and_selector->add_residue_selector( modeling_layer_selector );

					and_selector_virt->add_residue_selector( store_subset );
					and_selector_virt->add_residue_selector( virtual_layer_selector );

					needed_virtualization = and_selector_virt->apply( pose );

					///De-virtualize only those residues that were virtualized before and shouldn't still be virtualized.
					/// This is to make the visualization of the protocol clean.
					if ( ! first_modeling ) {
						TR << "De-Virtualizing foliage" << std::endl;
						ResidueSubset to_de_virtualize( pose.total_residue(), false );
						for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
							if ( pre_virtualized[ i ] == true && needed_virtualization[ i ] == false ) {
								to_de_virtualize[ i ] = true;
							}
						}
						store_to_de_virt_subset->set_residue_subset( to_de_virtualize );
						virt_to_real.set_residue_selector( store_to_de_virt_subset );
						virt_to_real.apply( pose );
					}

					TR << "Virtualizing foliage " << std::endl;
					pre_virtualized = and_selector_virt->apply( pose );
					real_to_virt.set_residue_selector( and_selector_virt );
					real_to_virt.apply( pose );

					TR << "Running GlycanRelax on layer: " << current_start << " " << current_end << std::endl;
					glycan_relax.set_selector( and_selector );
					glycan_relax.apply( pose );

					// Set the next layer.
					current_start = current_start + layer_size_ - window_size_;
					current_end = current_end + layer_size_ - window_size_;
					
					//Here, we just model the rest of the glycans at the end.
					if (current_end > max_end) {
						current_end = max_end;
					}

					first_modeling = false;

				}


				//Make sure everything is de-virtualized.
				store_subset->set_residue_subset( starting_subset );
				virt_to_real.set_residue_selector( store_subset );
				virt_to_real.apply( pose );
				break;

			}
			case backward :
				{
				TR << "Going in the backward direction " << std::endl;

				core::Size max_end = pose.glycan_tree_set()->get_largest_glycan_tree_layer();
				core::SSize current_end = max_end;
				core::SSize current_start = max_end - layer_size_ + 1;


				while ( current_start >= 0 ) {

					modeling_layer_selector->set_layer( current_start, current_end );
					//ResidueSubset modeling_layer = modeling_layer_selector->apply( pose );

					TR << "Running GlycanRelax on layer: " << current_start << " " << current_end << std::endl;

					//Combine modeling layer with whatver the current_subset is.
					store_subset->set_residue_subset( starting_subset );
					and_selector->clear();

					and_selector->add_residue_selector( store_subset );
					and_selector->add_residue_selector( modeling_layer_selector );

					glycan_relax.set_selector( and_selector );
					glycan_relax.apply( pose );

					current_end = current_end - layer_size_ + window_size_;
					current_start = current_start - layer_size_ + window_size_;
				}
				break;
			}

			}
			TR << "round " << i << " " << scorefxn_->score( pose ) << std::endl;

		} //End for loop

		completed_quenches_+=1;


	}

	TR << "Final Score: " << scorefxn_->score( pose ) << std::endl;
	if ( final_min_pack_min_ ) {

		PackRotamersMover packer = PackRotamersMover();
		MinMover min_mover = MinMover();

		//For now, we create a movemap from a selector manually.  Refactor this to Andrew's new code!
		MoveMapOP mm = MoveMapOP( new MoveMap());
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( ! min_subset[ resnum ] ) continue;
			mm->set_bb( resnum, true);
			mm->set_chi( resnum, true);
		}

		//Configure the scorefunction:
		packer.score_function( scorefxn_ );
		min_mover.score_function( scorefxn_ );


		//Configure the TaskFactory.
		TaskFactoryOP tf = get_all_glycans_and_neighbor_res_task_factory( min_subset );
		packer.task_factory( tf );

		TR << "Running final min-pack." << std::endl;
		core::Real prefinal = scorefxn_->score( pose );
		///Run min/pack/min/pack/min round.
		min_mover.apply( pose );
		packer.apply( pose );
		min_mover.apply( pose );
		packer.apply( pose);
		min_mover.apply( pose );

		TR << "Start- : " << starting_score << std::endl;
		TR << "Pre  - : " << prefinal << std::endl;
		TR << "Post - : " << scorefxn_->score( pose ) << std::endl;

	}

}




/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycanTreeRelaxCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GlycanTreeRelax );
}

std::string
GlycanTreeRelaxCreator::keyname() const
{
	return GlycanTreeRelax::mover_name();
}

void GlycanTreeRelaxCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanTreeRelax::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, GlycanTreeRelax const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //carbohydrates
