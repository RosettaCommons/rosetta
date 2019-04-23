// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanTreeModeler.cc
/// @brief A protocol for optimizing glycan trees using the GlycanSampler from the base of the tree out to the leaves.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Raemisch (raemisch@scripps.edu)

// Unit headers
#include <protocols/carbohydrates/GlycanTreeModeler.hh>
#include <protocols/carbohydrates/GlycanTreeModelerCreator.hh>
#include <protocols/carbohydrates/GlycanSampler.hh>
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
#include <core/select/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanTree.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>

// Protocl headers
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>


// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.carbohydrates.GlycanTreeModeler" );

namespace protocols {
namespace carbohydrates {

using namespace basic::options;
using namespace protocols::simple_moves;
using namespace protocols::minimization_packing;
using namespace core::scoring;
using namespace core::pose;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycanTreeModeler::GlycanTreeModeler():
	protocols::moves::Mover( GlycanTreeModeler::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycanTreeModeler::GlycanTreeModeler( GlycanTreeModeler const & src ):
	protocols::moves::Mover( src ),
	layer_size_(src.layer_size_),
	window_size_(src.window_size_),
	rounds_(src.rounds_),
	completed_quenches_(src.completed_quenches_),
	trees_to_model_(src.trees_to_model_),
	glycan_sampler_rounds_(src.glycan_sampler_rounds_),
	refine_(src.refine_),
	quench_mode_(src.quench_mode_),
	final_min_pack_min_(src.final_min_pack_min_),
	cartmin_(src.cartmin_),
	min_rings_(src.min_rings_),
	idealize_(src.idealize_),
	use_conformer_populations_(src.use_conformer_populations_),
	force_virts_for_refine_(src.force_virts_for_refine_),
	hybrid_protocol_(src.hybrid_protocol_),
	use_gaussian_sampling_( src.use_gaussian_sampling_),
	use_shear_( src.use_shear_)

{
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.selector_ ) selector_ = src.selector_->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycanTreeModeler::~GlycanTreeModeler()= default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycanTreeModeler::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
GlycanTreeModeler::parse_my_tag(
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
	glycan_sampler_rounds_ = tag->getOption< core::Size >("glycan_sampler_rounds", glycan_sampler_rounds_);

	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	}
	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

	min_rings_ = tag->getOption< core::Real >("min_rings", min_rings_);
	cartmin_ = tag->getOption< bool >("cartmin", cartmin_);
	idealize_ = tag->getOption< bool >("idealize", idealize_);

	set_force_virts_for_refinement(tag->getOption< bool >("force_virts_for_refinement", force_virts_for_refine_));
	set_use_conformer_probabilities( tag->getOption< bool >( "use_conformer_probs", use_conformer_populations_));
	set_hybrid_protocol( tag->getOption< bool >("hybrid_protocol", hybrid_protocol_));
	set_use_gaussian_sampling( tag->getOption< bool >( "use_gaussian_sampling", use_gaussian_sampling_));
	use_shear_ = tag->getOption< bool >("shear", use_shear_);
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
GlycanTreeModeler::fresh_instance() const
{
	return protocols::moves::MoverOP( new GlycanTreeModeler );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycanTreeModeler::clone() const
{
	return protocols::moves::MoverOP( new GlycanTreeModeler( *this ) );
}

std::string GlycanTreeModeler::get_name() const {
	return mover_name();
}

std::string GlycanTreeModeler::mover_name() {
	return "GlycanTreeModeler";
}

void GlycanTreeModeler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
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

		+ XMLSchemaAttribute::attribute_w_default("refine", xsct_rosetta_bool, "Do a refinement instead of a denovo model.  This means NO randomization at the start of the protocol and modeling is done in the context of the full glycan tree instead of using virtual (non-scored) residues during tree growth", "false")

		+ XMLSchemaAttribute("rounds", xsct_non_negative_integer, "Set the number of rounds.  We start with the forward direction, from the "
		"root/start out to the glycan tips.  Next round we go backward, then forward.  The number of rounds corresponds to how many back "
		"and forth directions we go")

		+ XMLSchemaAttribute::attribute_w_default("glycan_sampler_rounds", xsct_non_negative_integer, "Round Number for the internal the GlycanSampler.  Default is the default of the GlycanSampler.", "25")

		+ XMLSchemaAttribute::attribute_w_default("use_conformer_probs", xsct_rosetta_bool, "Use conformer probabilities instead of doing uniform sampling", "false")
		+ XMLSchemaAttribute::attribute_w_default("use_gaussian_sampling", xsct_rosetta_bool, "Set whether to build conformer torsions using a gaussian of the angle or through uniform sampling up to 1 SD (default)", "false")

		+ XMLSchemaAttribute::attribute_w_default("force_virts_for_refinement", xsct_rosetta_bool, "Refinement now models layers in the context of the full glycan tree.\n Turn this option on to use non-scored residues (virtuals) at the ends of the tree that are not being modeled.\n  This is how the de-novo protocol works.", "false")

		+ XMLSchemaAttribute::attribute_w_default("idealize", xsct_rosetta_bool, "Attempt to idealize the bond lengths/angles of glycan residues being modeled", "false")

		+ XMLSchemaAttribute::attribute_w_default("quench_mode", xsct_rosetta_bool, "Do quench mode for each glycan tree?", "false")
		+ XMLSchemaAttribute::attribute_w_default("final_min_pack_min", xsct_rosetta_bool, "Do a final set of cycles of min/pack", "true")
		+ XMLSchemaAttribute::attribute_w_default("min_rings", xsct_rosetta_bool, "Minimize Carbohydrate Rings during minimization. ", "false")
		+ XMLSchemaAttribute::attribute_w_default("cartmin", xsct_rosetta_bool, "Use Cartesian Minimization instead of Dihedral Minimization during packing steps.", "false")


		+ XMLSchemaAttribute::attribute_w_default("hybrid_protocol", xsct_rosetta_bool ,"Set to use an experimental protocol where we build out the glycans, but re-model the full tree during the building process" , "false")
		+ XMLSchemaAttribute::attribute_w_default("shear", xsct_rosetta_bool, "Use the Shear Mover that is now compatible with glycans at an a probability of 10 percent", "false");

	std::string documentation =
		"Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)\n"
		"Brief: A protocol for optimizing glycan trees using the GlycanSampler from the base of the tree out to the leaves.\n"
		"Details: Works by making all other residues virtual except the ones it is working on (current Layer).\n"
		"A virtual residue is not scored.\n"
		"It will start at the first glycan residues, and then move out to the edges.\n"
		"\n"
		"GENERAL ALGORITHM\n"
		"\n"
		"We start at the roots, and make all other glycan residues virtual.\n"
		"We first model towards the leaves and this is considered the forward direction.\n"
		"the GlycanSampler is used for the actual modeling, we only model a layer at a time, until we reach the tips.\n"
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
GlycanTreeModeler::set_layer_size( const core::Size layer_size ){
	layer_size_ = layer_size;
}

void
GlycanTreeModeler::set_window_size( core::Size const window_size ){
	window_size_ = window_size;
}

void
GlycanTreeModeler::set_refine(const bool refine){
	refine_ = refine;
}

void
GlycanTreeModeler::set_final_min_pack_min( const bool minpackmin ){
	final_min_pack_min_ = minpackmin;
}

void
GlycanTreeModeler::set_quench_mode( bool quench_mode ){
	quench_mode_ = quench_mode;
}

void
GlycanTreeModeler::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn->clone();
}

void
GlycanTreeModeler::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector->clone();
}

void
GlycanTreeModeler::set_rounds(const core::Size rounds){
	rounds_ = rounds;
}

void
GlycanTreeModeler::set_use_gaussian_sampling(bool gaussian_sampling){
	use_gaussian_sampling_ = gaussian_sampling;
}

///@brief Override Glycan Relax rounds.
void
GlycanTreeModeler::set_glycan_sampler_rounds( core::Size glycan_sampler_rounds){
	glycan_sampler_rounds_ = glycan_sampler_rounds;
}

bool
GlycanTreeModeler::is_quenched() const {
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

void
GlycanTreeModeler::setup_score_function() {

	//Create Scorefunction if needed.
	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}

	//If Cartmin is enabled, and scorefxn_ needs the terms, enable them!
	if ( cartmin_ ) {
		core::scoring::ScoreFunctionOP local= scorefxn_->clone();
		setup_cartmin(local);
		scorefxn_ = local;
	}
}

void
GlycanTreeModeler::setup_cartmin(core::scoring::ScoreFunctionOP scorefxn) const {

	scorefxn->set_weight_if_zero(core::scoring::cart_bonded, .5);
	scorefxn->set_weight(core::scoring::pro_close, 0);

}

void
GlycanTreeModeler::set_use_conformer_probabilities( bool use_conformer_stats ){
	use_conformer_populations_ = use_conformer_stats;
}

void
GlycanTreeModeler::set_force_virts_for_refinement(bool force_virts) {
	force_virts_for_refine_ = force_virts;
}

void
GlycanTreeModeler::set_hybrid_protocol(bool hybrid_protocol){
	hybrid_protocol_ = hybrid_protocol;
}

void
GlycanTreeModeler::set_use_shear(bool use_shear){
	use_shear_ = use_shear;
}

/// @brief Apply the mover
void
GlycanTreeModeler::apply( core::pose::Pose & pose){
	using namespace core::select::residue_selector;
	using namespace core::kinematics;
	using namespace core::pack::task;
	using namespace protocols::moves;
	using namespace protocols::minimization_packing;

	debug_assert( layer_size_ != window_size_ );

	//Setup everything we need.

	// Setup Selectors
	GlycanResidueSelectorOP glycan_selector = GlycanResidueSelectorOP( new GlycanResidueSelector() );
	GlycanLayerSelectorOP modeling_layer_selector = GlycanLayerSelectorOP( new GlycanLayerSelector() );
	GlycanLayerSelectorOP full_window_selector = GlycanLayerSelectorOP( new GlycanLayerSelector());
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

	if ( core::pose::symmetry::is_symmetric( pose )  )  {
		TR << "De-Symmetrizing selector" << std::endl;
		starting_subset = core::select::get_master_subunit_selection(pose, starting_subset);
	}

	ResidueSubset min_subset = starting_subset;

	setup_score_function();
	scorefxn_->score( pose );


	//Check that we don't have FUNKY things in our ResidueSelector.
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( starting_subset[ i ] ) {
			if ( ! pose.residue_type( i ).is_carbohydrate() ) {
				utility_exit_with_message(" GlycanTreeModeler: Residue "+utility::to_string(i)+" set in ResidueSelector, but not a carbohdyrate residue!");
			}
		}
	}

	core::Size total_glycan_residues = count_selected(starting_subset);
	if ( total_glycan_residues == 0 ) {
		utility_exit_with_message("GlycanTreeModeler: No glycan residues to model.  Cannot continue.");
	}

	// Setup the GlycanSampler
	GlycanSampler glycan_sampler = GlycanSampler();
	if ( scorefxn_ ) {
		glycan_sampler.set_scorefunction( scorefxn_ );
	}

	//Control our own refinement here as we may have multiple apply calls to Glycan Sampler.
	// Without this, we would randomize EVERY time GlycaSampler's apply function is called.
	glycan_sampler.set_randomize_first( false );
	if ( ! refine_ ) {
		glycan_sampler.randomize_glycan_torsions(pose, starting_subset);

	}

	//Only override cmd-line settings of the GlycanSampler if it is set here.  Otherwise, cmd-line controls do not work.
	if ( glycan_sampler_rounds_ != 0 ) {
		glycan_sampler.set_rounds(glycan_sampler_rounds_);
	}

	//Setup Glycan Sampler Options
	glycan_sampler.set_use_shear(use_shear_);
	glycan_sampler.set_refine( refine_ );
	glycan_sampler.use_cartmin( cartmin_ );
	glycan_sampler.set_min_rings( min_rings_ );
	glycan_sampler.set_population_based_conformer_sampling(use_conformer_populations_);
	glycan_sampler.set_use_gaussian_sampling( use_gaussian_sampling_);

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

	if ( idealize_ ) {
		TR << "Idealizing full glycan trees " << std::endl;
		glycan_sampler.idealize_glycan_residues(pose, tree_subset);
	}

	core::Real starting_score = scorefxn_->score( pose );
	TR << "Starting Score: " << starting_score << std::endl;

	ResidueSubset mask = starting_subset; //All the glycan residues we will be modeling

	MonteCarloOP mc = MonteCarloOP( new MonteCarlo(pose, *scorefxn_, 1.0) );
	mc->set_last_accepted_pose(pose);

	core::Size default_sampler_rounds = glycan_sampler.get_glycan_sampler_rounds();

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
						TR.Error << "Maximum number of layers smaller than the layer size.  Either decrease the layer_size or use the GlycanSampler instead.  Layer-based sampling does not make sense here." << std::endl;
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
						for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
							if ( pre_virtualized[ ii ] && ! needed_virtualization[ ii ] ) {
								to_de_virtualize[ ii ] = true;
							}
						}
						store_to_de_virt_subset->set_residue_subset( to_de_virtualize );

						if ( (! refine_)  && (! force_virts_for_refine_) ) {
							virt_to_real.set_residue_selector( store_to_de_virt_subset );
							virt_to_real.apply( pose );
						}
					}

					if ( (! refine_)  && (! force_virts_for_refine_) ) {
						TR << "Virtualizing new foliage layer" << std::endl;
						pre_virtualized = and_selector_virt->apply( pose );
						real_to_virt.set_residue_selector( and_selector_virt );
						real_to_virt.apply( pose );
					}

					TR << "Running the GlycanSampler on layer [ start -> end (including) ]: " << current_start << " " << current_end << std::endl;

					core::Size n_in_layer = count_selected( and_selector->apply( pose ) );

					if ( n_in_layer > 0 ) {

						if ( hybrid_protocol_ && ! first_modeling ) {
							glycan_sampler.set_rounds(default_sampler_rounds/2);
							glycan_sampler.set_selector( and_selector );
							glycan_sampler.apply( pose );

							TR << "Running full window protocol sampling" << std::endl;

							//Second, we do the layer plus everything else

							full_window_selector->set_layer( 0, current_end);
							and_selector->clear();
							and_selector->add_residue_selector(store_subset);

							//This is so we get exactly equivalent amount of random sampling from both normal and the full window
							// protocol - so that we may actually benchmark it properly.
							core::Size total_rounds = n_in_layer * (default_sampler_rounds/2);
							glycan_sampler.force_total_rounds(total_rounds);

							store_subset->set_residue_subset( starting_subset );
							and_selector->clear();
							and_selector->add_residue_selector(store_subset);
							and_selector->add_residue_selector( full_window_selector );
							glycan_sampler.set_selector( and_selector );
							glycan_sampler.apply(pose);
						} else {
							TR << "Applying normal protocol" << std::endl;
							glycan_sampler.set_selector( and_selector );
							glycan_sampler.apply( pose );
						}



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
				if ( (! refine_)  && (! force_virts_for_refine_) )  {
					store_subset->set_residue_subset( starting_subset );
					virt_to_real.set_residue_selector( store_subset );
					virt_to_real.apply( pose );
				}
				break;

			}
			case backward :
				{
				TR << "Going in the backward direction " << std::endl;

				core::Size max_end2 = pose.glycan_tree_set()->get_largest_glycan_tree_layer();
				core::SSize current_end = max_end2;
				core::SSize current_start = max_end2 - layer_size_ + 1;


				while ( current_start >= 0 ) {

					modeling_layer_selector->set_layer( current_start, current_end );
					//ResidueSubset modeling_layer = modeling_layer_selector->apply( pose );

					TR << "Running the GlycanSampler on layer: " << current_start << " " << current_end << std::endl;

					//Combine modeling layer with whatver the current_subset is.
					store_subset->set_residue_subset( starting_subset );
					and_selector->clear();

					and_selector->add_residue_selector( store_subset );
					and_selector->add_residue_selector( modeling_layer_selector );



					//If there are no residues in the layer, we continue.  With glycan masking, this is actually possible.
					// Not recommended by any means, but definitely posssible..
					core::Size n_in_layer = count_selected( and_selector->apply( pose ) );

					if ( n_in_layer > 0 ) {
						glycan_sampler.set_selector( and_selector );
						glycan_sampler.apply( pose );
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


		//For now, we create a movemap from a selector manually.  Refactor this to Andrew's new code!
		ReturnResidueSubsetSelectorOP return_subset_min = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector(min_subset));
		core::Size n_glycan_residues = count_selected(min_subset);

		MoveMapOP mm = core::pose::carbohydrates::create_glycan_movemap_from_residue_selector(pose,return_subset_min,true /*chi*/,min_rings_,true /*bb*/, cartmin_);

		MinMover min_mover = MinMover(mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.01, true);

		if ( (! option [OptionKeys::run::nblist_autoupdate].user()) && (! option [OptionKeys::run::nblist_autoupdate]() ) ) {
			min_mover.min_options()->nblist_auto_update( true );
		}

		//Configure the scorefunction:
		packer.score_function( scorefxn_ );

		//Configure the TaskFactory.
		TaskFactoryOP tf = get_all_glycans_and_neighbor_res_task_factory( min_subset );
		packer.task_factory( tf );

		//Configure the ShearMover.
		ShearMover shear = ShearMover();
		shear.set_residue_selector(return_subset_min);

		TR << "Running final min-pack." << std::endl;
		core::Real prefinal = scorefxn_->score( pose );
		///Run min/pack/min/pack/min round.

		run_shear_min_pack(min_mover, packer, shear, *mc, n_glycan_residues, pose, use_shear_);
		run_shear_min_pack(min_mover, packer, shear, *mc, n_glycan_residues, pose, use_shear_);

		min_mover.apply(pose);

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
GlycanTreeModelerCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GlycanTreeModeler );
}

std::string
GlycanTreeModelerCreator::keyname() const
{
	return GlycanTreeModeler::mover_name();
}

void GlycanTreeModelerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanTreeModeler::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, GlycanTreeModeler const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //carbohydrates
