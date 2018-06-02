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
#include <protocols/moves/MonteCarlo.hh>


// Core headers
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/GlycanLayerSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>

// Protocl headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.carbohydrates.GlycanTreeRelax" );

namespace protocols {
namespace carbohydrates {

using namespace protocols::simple_moves;
using namespace protocols::minimization_packing;
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
	final_min_pack_min_(src.final_min_pack_min_),
	cartmin_(src.cartmin_),
	min_rings_(src.min_rings_)

{
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanTreeRelax::~GlycanTreeRelax()= default;

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

	min_rings_ = tag->getOption< core::Real >("min_rings", min_rings_);
	cartmin_ = tag->getOption< bool >("cartmin", cartmin_);
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
		"  Details: \n"
		"  The distance that make up a layer.  If we have a distance of 2,\n"
		"  we first model all glycans that are equal to or less than 2 residue distance to the root.\n"
		"  We then slide this layer up.  So we take all residues that have a distance between 3 and 1, and so on.\n")


		+ XMLSchemaAttribute("window_size", xsct_non_negative_integer,
		"Brief: Set the window size.  This is the overlap of the layers during modeling. \n"
		"  \n"
		"  Details: \n"
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
		+ XMLSchemaAttribute::attribute_w_default("min_rings", xsct_rosetta_bool, "Minimize Carbohydrate Rings during minimization. ", "false")
		+ XMLSchemaAttribute::attribute_w_default("cartmin", xsct_rosetta_bool, "Use Cartesian Minimization instead of Dihedral Minimization during packing steps.", "false")
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
		"Typically, we would want at least a window size of 1. \n"
		"--Residue Selection-- \n"
		"\n"
		" You do not need a ResidueSelector passed in.  It will select all glycan residues automatically.\n"
		" However, if you do, you must only pass in glycan residues.  See the GlycanResidueSelector and the GlycanLayerSelector for a very easy way to select specific glycan trees and residues. ";

	rosetta_scripts::attributes_for_parse_residue_selector( attlist, "Residue Selector containing only glycan residues.  This is not needed, as this class will automatically select ALL glycan residues in the pose to model.  See the GlycanResidueSelector and the GlycanLayerSelector for control glycan selection.  Note that the ASN is not technically a glycan.  Since dihedral angles are defined for a sugar from the upper to lower residue, the dihedral angles between the first glycan and the ASN are defined by the first glycan. " );

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
	using namespace protocols::moves;

	debug_assert( layer_size_ != window_size_ );

	//Setup everything we need.

	// Setup Selectors
	GlycanResidueSelectorOP glycan_selector = GlycanResidueSelectorOP( new GlycanResidueSelector() );
	GlycanLayerSelectorOP modeling_layer_selector = GlycanLayerSelectorOP( new GlycanLayerSelector() );
	GlycanLayerSelectorOP virtual_layer_selector =  GlycanLayerSelectorOP( new GlycanLayerSelector() );
	ReturnResidueSubsetSelectorOP store_subset = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector() );
	ReturnResidueSubsetSelectorOP store_subset2 =ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector() );
	ReturnResidueSubsetSelectorOP store_to_de_virt_subset =  ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector() );
	AndResidueSelectorOP and_selector =      AndResidueSelectorOP( new AndResidueSelector() );
	AndResidueSelectorOP and_selector_virt = AndResidueSelectorOP( new AndResidueSelector() );


	Pose starting_pose = pose;


	ResidueSubset starting_subset;
	//ResidueSubset current_subset;

	if ( selector_ ) {
		starting_subset = selector_->apply( pose );

	} else {
		starting_subset = glycan_selector->apply( pose );
	}

	ResidueSubset min_subset = starting_subset;

	if ( ! scorefxn_ ) {
		scorefxn_ = get_score_function();
	}
	scorefxn_->score( pose );


	//Check that we don't have FUNKY things in our ResidueSelector.
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( starting_subset[ i ] ) {
			if ( ! pose.residue_type( i ).is_carbohydrate() ) {
				utility_exit_with_message(" GlycanTreeRelax: Residue "+utility::to_string(i)+" set in ResidueSelector, but not a carbohdyrate residue!");
			}
		}
	}

	core::Size total_glycan_residues = count_selected(starting_subset);
	if ( total_glycan_residues == 0 ) {
		utility_exit_with_message("GlycanTreeRelax: No glycan residues to model.  Cannot continue.");
	}

	// Setup GlycanRelax
	GlycanRelaxMover glycan_relax = GlycanRelaxMover();
	if ( scorefxn_ ) {
		glycan_relax.set_scorefunction( scorefxn_ );
	}
	glycan_relax.set_refine( true );

	if ( ! refine_ ) {
		glycan_relax.randomize_glycan_torsions(pose, starting_subset);

	}

	//Only override cmd-line settings of GlycanRelax if it is set here.  Otherwise, cmd-line controls do not work.
	if ( glycan_relax_rounds_ != 0 ) {
		glycan_relax.set_rounds(glycan_relax_rounds_);
	}

	glycan_relax.use_cartmin( cartmin_ );
	glycan_relax.set_min_rings( min_rings_ );

	ConvertRealToVirtualMover real_to_virt = ConvertRealToVirtualMover();
	ConvertVirtualToRealMover virt_to_real = ConvertVirtualToRealMover();

	QuenchDirection current_direction = forward;

	utility::vector1< core::Size > all_trees = pose.glycan_tree_set()->get_start_points();
	utility::vector1< core::Size > tree_subset; //Quench-mode trees.

	//JAB - Check to make sure we are using all glycans.
	for ( core::Size i = 1; i <= starting_subset.size(); ++i ) {
		if ( starting_subset[ i ] && !pose.residue( i).is_carbohydrate() ) {
			utility_exit_with_message("Non-Carbohydrate residue in selection. Please see the GlycanResidueSelector for custom input.  Not setting a selector will automatically model ALL glycan residues in the pose. Residue: "+utility::to_string( i )+" : "+pose.pdb_info()->pose2pdb( i ));
		}

	}
	core::Size max_end = pose.glycan_tree_set()->get_largest_glycan_tree_layer( starting_subset );
	core::Size min_layer = pose.glycan_tree_set()->get_smallest_glycan_tree_layer( starting_subset );

	TR << "Smallest glycan layer: " << min_layer << std::endl;
	TR << "Largest  glycan layer: " << max_end << std::endl;

	//Setup for quench mode
	if ( selector_ ) {
		for ( core::Size tree_start : all_trees ) {

			if ( starting_subset[ tree_start] ) {
				tree_subset.push_back( tree_start );
				continue;
			} else {
				utility::vector1< core::Size > tree_residues = pose.glycan_tree_set()->get_tree( tree_start )->get_residues();

				for ( core::Size & resnum : tree_residues ) {

					if ( starting_subset[ resnum ] ) {
						tree_subset.push_back( tree_start );
						break;
					}
				}
				// For Tree Residues
			}

		} //For All Trees
	} else {
		tree_subset = all_trees;
	}

	trees_to_model_ = tree_subset.size();

	TR << "Ntrees to model " << utility::to_string(trees_to_model_) << std::endl;

	core::Real starting_score = scorefxn_->score( pose );
	TR << "Starting Score: " << starting_score << std::endl;

	ResidueSubset mask = starting_subset; //All the glycan residues we will be modeling

	MonteCarloOP mc = MonteCarloOP( new MonteCarlo(pose, *scorefxn_, 1.0) );
	mc->set_last_accepted_pose(pose);

	completed_quenches_ = 0;
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
			store_subset->set_residue_subset( starting_subset );
			store_subset2->set_residue_subset( mask );
			and_selector->clear();

			and_selector->add_residue_selector( store_subset );
			and_selector->add_residue_selector( store_subset2 );

			starting_subset = and_selector->apply( pose ); //Combine the particular Tree and the overall glycan Mask.

			and_selector->clear();


			//Remove the tree from the list to model.
			tree_subset.pop( tree_start );

			min_layer = pose.glycan_tree_set()->get_smallest_glycan_tree_layer( starting_subset );
			max_end = pose.glycan_tree_set()->get_largest_glycan_tree_layer( starting_subset );

		}
		for ( core::Size i = 1; i <= rounds_; ++i ) {
			TR << "starting round " << i << std::endl;

			if ( i % 2 == 0 ) {
				current_direction = backward;
			} else {
				current_direction = forward;
			}

			switch ( current_direction ) {
			case forward :
				{
				TR << "Going in the forward direction " << std::endl;

				bool first_modeling = true;
				core::Size current_start = min_layer;
				core::Size current_end = layer_size_ - 1;

				TR << "Modeling up to max end: " << max_end << std::endl;

				if ( max_end + 1 <= layer_size_ ) {
					if ( ! quench_mode_ ) {
						TR.Error << "Maximum number of layers smaller than the layer size.  Either decrease the layer_size or use GlycanRelax instead.  Layer-based sampling does not make sense here." << std::endl;
						set_last_move_status(FAIL_DO_NOT_RETRY);
						pose = starting_pose;
						return;
					}
				}
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
						TR << "De-Virtualizing current foliage layer" << std::endl;
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

					TR << "Virtualizing new foliage layer" << std::endl;
					pre_virtualized = and_selector_virt->apply( pose );
					real_to_virt.set_residue_selector( and_selector_virt );
					real_to_virt.apply( pose );

					TR << "Running GlycanRelax on layer [ start -> end (including) ]: " << current_start << " " << current_end << std::endl;

					core::Size n_in_layer = count_selected( and_selector->apply( pose ) );

					if ( n_in_layer > 0 ) {
						glycan_relax.set_selector( and_selector );
						glycan_relax.apply( pose );
					} else {
						TR << "No residues in layer, continueing." << std::endl;
					}

					// Set the next layer.
					current_start = current_start + layer_size_ - window_size_;
					current_end = current_end + layer_size_ - window_size_;

					//Here, we just model the rest of the glycans at the end.
					if ( current_end > max_end ) {
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



					//If there are no residues in the layer, we continue.  With glycan masking, this is actually possible.
					// Not recommended by any means, but definitely posssible..
					core::Size n_in_layer = count_selected( and_selector->apply( pose ) );

					if ( n_in_layer > 0 ) {
						glycan_relax.set_selector( and_selector );
						glycan_relax.apply( pose );
					} else {
						TR << "No residues in layer, continueing." << std::endl;
					}

					current_end = current_end - layer_size_ + window_size_;
					current_start = current_start - layer_size_ + window_size_;
				}
				break;
			}

			} //End Switch, end of Round.

			TR << "round " << i << " " << scorefxn_->score( pose ) << std::endl;
			bool accepted = mc->boltzmann( pose );
			TR << "accepted " << accepted << std::endl << std::endl;
			TR << "round " << i << " " << scorefxn_->score( pose ) << " final " << std::endl << std::endl;

		} //End for loop

		completed_quenches_+=1;

	}

	TR << "Final Score: " << scorefxn_->score( pose ) << std::endl;
	if ( final_min_pack_min_ ) {

		PackRotamersMover packer = PackRotamersMover();
		MinMover min_mover = MinMover();

		//For now, we create a movemap from a selector manually.  Refactor this to Andrew's new code!
		ReturnResidueSubsetSelectorOP return_subset_min = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector(min_subset));

		MoveMapOP mm = core::pose::carbohydrates::create_glycan_movemap_from_residue_selector(pose,return_subset_min,true /*chi*/,min_rings_,true /*bb*/, cartmin_);

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

		mc->boltzmann(pose);

	}
	mc->recover_low(pose);

	TR << "Final- : " << scorefxn_->score( pose ) << std::endl;

	if ( rounds_ > 1 ) {
		mc->show_counters();
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
