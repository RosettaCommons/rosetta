// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanSampler.cc
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)

#include <protocols/carbohydrates/GlycanSampler.hh>
#include <protocols/carbohydrates/GlycanSamplerCreator.hh>
#include <protocols/carbohydrates/LinkageConformerMover.hh>
#include <protocols/carbohydrates/GlycanTreeMinMover.hh>
#include <protocols/carbohydrates/util.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/symmetry/util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/id/types.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/RandomGlycanFoliageSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/GlycanLayerSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/SymmetricalResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/simple_metrics/metrics/TimingProfileMetric.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>
#include <protocols/simple_moves/bb_sampler/SmallBBSampler.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.carbohydrates.GlycanSampler" );


namespace protocols {
namespace carbohydrates {

using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::bb_sampler;
using namespace protocols::simple_moves::rational_mc;
using namespace protocols::minimization_packing;

using namespace core::pack::task;
using namespace basic::options;
using namespace core::select::residue_selector;
using namespace core::kinematics;
using namespace core::pose::carbohydrates;

GlycanSampler::GlycanSampler():
	protocols::moves::Mover( "GlycanSampler" )
{
	set_defaults();
}

GlycanSampler::GlycanSampler(
	core::select::residue_selector::ResidueSelectorCOP selector,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size rounds):
	protocols::moves::Mover("GlycanSampler")

{
	scorefxn_ = scorefxn->clone();
	selector_ = selector->clone();
	set_defaults();
	set_rounds(rounds);
}

GlycanSampler::~GlycanSampler()= default;


void
GlycanSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	using namespace core::select::residue_selector;

	kt_ = tag->getOption< core::Real >("kt", kt_);
	rounds_ = tag->getOption< core::Size >("rounds", rounds_);

	final_min_ = tag->getOption< bool >("final_min", final_min_);

	pymol_movie_ = tag->getOption< bool >("pymol_movie", pymol_movie_);

	//Scorefunction.
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, datamap );
	}

	if ( tag->hasOption("task_operations") ) {
		TaskFactoryOP tf( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
		set_taskfactory( tf );
	}

	pack_distance_ = tag->getOption< core::Real >("pack_distance", pack_distance_);

	cartmin_ = tag->getOption< bool >("cartmin", cartmin_);

	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	}

	tree_based_min_pack_ = tag->getOption< bool >("tree_based_min_pack", tree_based_min_pack_);
	population_based_conformer_sampling_ = tag->getOption< bool >("use_conformer_probs", population_based_conformer_sampling_);
	use_gaussian_sampling_ = tag->getOption< bool >("use_gaussian_sampling", use_gaussian_sampling_);
	min_rings_ = tag->getOption< core::Real >("min_rings", min_rings_);
	use_shear_ = tag->getOption< bool >("shear", use_shear_);
	randomize_first_ = tag->getOption< bool >("randomize_torsions", randomize_first_);
	inner_ncycles_ = tag->getOption< bool >("inner_bb_cycles", inner_ncycles_);
	match_sampling_of_modeler_ = tag->getOption< bool >("match_sampling_of_modeler", match_sampling_of_modeler_);

}

void GlycanSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute("kt", xsct_real, "Temperature for metropolis criterion")
		+ XMLSchemaAttribute("rounds", xsct_non_negative_integer, "Number of relax rounds to perform.  Will be multiplied by the number of glycan residues.")
		+ XMLSchemaAttribute::attribute_w_default("final_min", xsct_rosetta_bool, "Perform a final minimization", "true")
		+ XMLSchemaAttribute::attribute_w_default("pymol_movie", xsct_rosetta_bool, "Output a PyMOL movie of the run", "false")

		+ XMLSchemaAttribute::attribute_w_default("refine", xsct_rosetta_bool, "Do not start with a random glycan conformation.", "false")

		+ XMLSchemaAttribute("pack_distance", xsct_real, "Neighbor distance for packing")
		+ XMLSchemaAttribute::attribute_w_default("cartmin", xsct_rosetta_bool, "Use Cartesian Minimization instead of Dihedral Minimization during packing steps.", "false")
		+ XMLSchemaAttribute::attribute_w_default("tree_based_min_pack", xsct_rosetta_bool, "Use Tree-based minimization and packing instead of minimizing and packing ALL residues each time we min.  Significantly impacts runtime.  If you are seeing crappy structures for a few sugars, turn this off.  This is default-on to decrease runtime for a large number of glycans.", "true")
		+ XMLSchemaAttribute::attribute_w_default("use_conformer_probs", xsct_rosetta_bool, "Use the populations of the conformers as probabilities during our linkage conformer sampling.  This makes it harder to overcome energy barriers with more-rare conformers!", "false")
		+ XMLSchemaAttribute::attribute_w_default("use_gaussian_sampling", xsct_rosetta_bool, "Set whether to build conformer torsions using a gaussian of the angle or through uniform sampling up to 1 SD (default)", "false")
		+ XMLSchemaAttribute::attribute_w_default("min_rings", xsct_rosetta_bool, "Minimize Carbohydrate Rings during minimization. ", "false")
		+ XMLSchemaAttribute::attribute_w_default("shear", xsct_rosetta_bool, "Use the Shear Mover that is now compatible with glycans at an a probability of 10 percent", "false")
		+ XMLSchemaAttribute::attribute_w_default("randomize_torsions", xsct_rosetta_bool, "If NOT doing refinement, control whether we randomize torsions at the start, which helps to achieve low energy structures.", "true")
		+ XMLSchemaAttribute::attribute_w_default("match_sampling_of_modeler", xsct_rosetta_bool, "Option that matches the amount of sampling in a default GlycanTreeModeler run.  Used for benchmarking.", "0")
		+ XMLSchemaAttribute::attribute_w_default("inner_bb_cycles", xsct_non_negative_integer, "Inner cycles for each time we call BB sampling either through small/medium/large moves or through the SugarBB Sampler.  This is multiplied by the number of glycans.  Scales poorly with GlycanModeler.  If 0 (default), we run the protocol normally", "0");
	//Append for MoveMap, scorefxn, and task_operation tags.
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_get_score_function_name( attlist );
	rosetta_scripts::attributes_for_parse_residue_selector( attlist, "Residue Selector containing only glycan residues.  This is not needed, as this class will automatically select ALL glycan residues in the pose to model.  See the GlycanResidueSelector and the GlycanLayerSelector for control glycan selection.  Note that the ASN is not technically a glycan.  Since dihedral angles are defined for a sugar from the upper to lower residue, the dihedral angles between the first glycan and the ASN are defined by the first glycan. " );

	XMLSchemaSimpleSubelementList subelements;
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(),

		"Authors: Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason Labonte (WLabonte@jhu.edu)\n"
		"Main mover for Glycan Relax, which optimizes glycans in a pose. "
		"Each round optimizes either one residue for BB sampling, linkage, or a [part of a ] branch for minimization and packing."
		"Minimization and packing work by default by selecting a random glycan residue from any set movemap or all of them "
		"and then selecting the rest of the downstream branch.  Those residues are then minimized or packed.  Packing includes a neighbor packing shell. "
		"Currently uses a random sampler with a set of weights to each mover for sampling.", attlist, subelements );
}

void
GlycanSampler::set_defaults(){
	test_ = false;

	final_min_ = true;

	set_cmd_line_defaults();

	total_glycan_residues_ = 0;
	pack_distance_ = 6.0;
	cartmin_ = false;

}


void
GlycanSampler::set_cmd_line_defaults(){

	rounds_ = option [OptionKeys::carbohydrates::glycan_sampler::glycan_sampler_rounds]();
	test_ = option [OptionKeys::carbohydrates::glycan_sampler::glycan_sampler_test]();
	final_min_ = option[ OptionKeys::carbohydrates::glycan_sampler::final_min_glycans]();
	pymol_movie_ = option[ OptionKeys::carbohydrates::glycan_sampler::glycan_sampler_movie]();
	kt_ = option[ OptionKeys::carbohydrates::glycan_sampler::glycan_sampler_kt]();
	cartmin_ = option[ OptionKeys::carbohydrates::glycan_sampler::cartmin]();
	refine_ = option[ OptionKeys::carbohydrates::glycan_sampler::glycan_sampler_refine]();
	tree_based_min_pack_ = option[ OptionKeys::carbohydrates::glycan_sampler::tree_based_min_pack]();
	population_based_conformer_sampling_ = option[ OptionKeys::carbohydrates::glycan_sampler::population_based_conformer_sampling]();
	use_gaussian_sampling_ = option[ OptionKeys::carbohydrates::glycan_sampler::use_gaussian_sampling]();
	min_rings_ = option[ OptionKeys::carbohydrates::glycan_sampler::min_rings]();

}

GlycanSampler::GlycanSampler( GlycanSampler const & src ):
	protocols::moves::Mover(src),
	rounds_( src.rounds_ ),
	kt_( src.kt_ ),
	accept_log_( src.accept_log_ ),
	test_( src.test_ ),
	final_min_( src.final_min_ ),
	refine_( src.refine_ ),
	total_glycan_residues_( src.total_glycan_residues_ ),
	pymol_movie_( src.pymol_movie_ ),
	parsed_positions_( src.parsed_positions_),
	pack_distance_( src.pack_distance_ ),
	cartmin_( src.cartmin_ ),
	tree_based_min_pack_( src.tree_based_min_pack_ ),
	population_based_conformer_sampling_(src.population_based_conformer_sampling_),
	use_gaussian_sampling_(src.use_gaussian_sampling_),
	min_rings_(src.min_rings_),
	forced_total_rounds_(src.forced_total_rounds_),
	use_shear_(src.use_shear_),
	randomize_first_(src.randomize_first_),
	inner_ncycles_( src.inner_ncycles_),
	match_sampling_of_modeler_( src.match_sampling_of_modeler_),
	final_residue_subset_( src.final_residue_subset_)
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();

	if ( src.min_mover_ ) {
		min_mover_ = utility::pointer::make_shared< protocols::minimization_packing::MinMover >( *src.min_mover_);
	}
	if ( src.packer_ ) {
		packer_ = utility::pointer::make_shared< PackRotamersMover >( * src.packer_);
	}
	if ( src.shear_ ) {
		shear_ = utility::pointer::make_shared< ShearMover>( *src.shear_ );
	}
	if ( src.linkage_mover_ ) linkage_mover_ = utility::pointer::make_shared< LinkageConformerMover >( * src.linkage_mover_);
	if ( src.weighted_random_mover_ ) weighted_random_mover_ = utility::pointer::make_shared< moves::RandomMover >( *src.weighted_random_mover_);
	if ( src.tf_ ) tf_ = src.tf_->clone();
	if ( src.mc_ ) mc_ = src.mc_->clone();



}



protocols::moves::MoverOP
GlycanSampler::clone() const{
	return utility::pointer::make_shared< GlycanSampler >( *this );
}

moves::MoverOP
GlycanSampler::fresh_instance() const
{
	return utility::pointer::make_shared< GlycanSampler >();
}

// XRW TEMP std::string
// XRW TEMP GlycanSampler::get_name() const {
// XRW TEMP  return "GlycanSampler";
// XRW TEMP }

void
GlycanSampler::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, GlycanSampler const &mover)
{
	mover.show(os);
	return os;
}

void
GlycanSampler::set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector ){
	selector_ = selector->clone();
}

void
GlycanSampler::set_taskfactory(core::pack::task::TaskFactoryCOP tf){
	tf_ = tf;
}

void
GlycanSampler::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

void
GlycanSampler::set_inner_bb_cycles( core::Size inner_ncycles ){
	inner_ncycles_ = inner_ncycles;
}

void
GlycanSampler::setup_default_task_factory(utility::vector1< bool > const & glycan_positions, core::pose::Pose const & pose ){
	using namespace core::pack::task::operation;
	using namespace core::select::residue_selector;

	TaskFactoryOP tf = utility::pointer::make_shared< TaskFactory >();




	//Select all the NON-glycan and non-neighbors and then turn them off.
	if ( tree_based_min_pack_ ) {

		TR << "Setting up Tree-based packing." << std::endl;

		tf->push_back(utility::pointer::make_shared< InitializeFromCommandline >());
		RandomGlycanFoliageSelectorOP select_random_foliage= utility::pointer::make_shared< RandomGlycanFoliageSelector >(glycan_positions);
		PreventRepackingRLTOP prevent_repacking = utility::pointer::make_shared< PreventRepackingRLT >();

		NeighborhoodResidueSelectorOP neighbor_selector = utility::pointer::make_shared< NeighborhoodResidueSelector >(select_random_foliage, pack_distance_, true /* include focus */);

		OperateOnResidueSubsetOP subset_op = utility::pointer::make_shared< OperateOnResidueSubset >( prevent_repacking, neighbor_selector, true /* flip */);


		//Prevent repacking of virtual residues.  Not sure if they are actually repacked, so we make sure they are not here.
		ResidueSubset virtual_residues( pose.total_residue(), false);
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue( i ).is_virtual_residue() ) {
				virtual_residues[ i ] = true;
			}
		}

		ReturnResidueSubsetSelectorOP store_subset = utility::pointer::make_shared< ReturnResidueSubsetSelector >();
		store_subset->set_residue_subset( virtual_residues );
		OperateOnResidueSubsetOP operate_on_virt_subset = utility::pointer::make_shared< OperateOnResidueSubset >( prevent_repacking, store_subset );

		tf->push_back( subset_op );
		tf->push_back( utility::pointer::make_shared< RestrictToRepacking >());

	} else {
		tf = get_all_glycans_and_neighbor_res_task_factory( glycan_positions );
	}

	set_taskfactory(tf);
}

core::Size
GlycanSampler::get_glycan_sampler_rounds(){
	return rounds_;
}

void
GlycanSampler::force_total_rounds(core::Size total_rounds){
	forced_total_rounds_ = total_rounds;
}

void
GlycanSampler::set_scorefunction(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

void
GlycanSampler::set_kt(core::Real kt){
	kt_ = kt;
}

void
GlycanSampler::set_rounds(core::Size rounds){
	rounds_ = rounds;
}

void
GlycanSampler::use_cartmin( bool use_cartmin ){
	cartmin_ = use_cartmin;
}

void
GlycanSampler::set_refine( bool refine ){
	refine_ = refine;
}

///@brief set to minimize ring torsions during minimzation.  Default false.
void
GlycanSampler::set_min_rings( bool min_rings){
	min_rings_ = min_rings;
}

void
GlycanSampler::set_population_based_conformer_sampling(bool pop_based_sampling){
	population_based_conformer_sampling_ = pop_based_sampling;
}

void
GlycanSampler::set_use_gaussian_sampling(bool gaussian_sampling){
	use_gaussian_sampling_ = gaussian_sampling;
}

void
GlycanSampler::set_use_shear(bool use_shear){
	use_shear_ = use_shear;
}

void
GlycanSampler::set_randomize_first(bool randomize_first){
	randomize_first_ = randomize_first;
}

void
GlycanSampler::setup_score_function() {

	//Create Scorefunction if needed.
	if ( ! scorefxn_ ) {
		scorefxn_ = core::scoring::get_score_function();
	}

	//If Cartmin is enabled, and scorefxn_ needs the terms, enable them!
	if ( cartmin_ ) {
		core::scoring::ScoreFunctionOP local_scorefxn = scorefxn_->clone();
		setup_cartmin(local_scorefxn);
		scorefxn_ = local_scorefxn;
	}
}


void
GlycanSampler::setup_cartmin(core::scoring::ScoreFunctionOP scorefxn) const {

	scorefxn->set_weight_if_zero(core::scoring::cart_bonded, .5);
	scorefxn->set_weight(core::scoring::pro_close, 0);

}

void
GlycanSampler::randomize_glycan_torsions(core::pose::Pose &pose, utility::vector1<bool> const & subset) const{

	TR << "Randomizing glycan torsions " << std::endl;
	SmallBBSampler random_sampler = SmallBBSampler(  360.0 ); // +/- 180 degrees
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! subset[ i ] ) continue;
		//TR << "Randomizing " << i << std::endl;
		core::Size n_dihedrals = get_n_glycosidic_torsions_in_res( pose, i );
		core::Size parent_res = pose.glycan_tree_set()->get_parent( i ) ;
		for ( core::Size torsion_id = 1; torsion_id <= n_dihedrals; ++torsion_id ) {

			if ( parent_res != 0 ) {
				//TR << "Attempting " << i << " torsion " << torsion_id << std::endl;
				random_sampler.set_torsion_type( static_cast< core::id::MainchainTorsionType >( torsion_id ) );
				random_sampler.set_torsion_to_pose( pose, i );
				//TR << "Randomized " << i << " torsion " << torsion_id << std::endl;
			}
		}
	}
	//pose.dump_pdb("randomized.pdb");
	//TR << "Done Randomizing " << std::endl;

}

void
GlycanSampler::idealize_glycan_residues(core::pose::Pose & , utility::vector1< core::Size > const & ) const {
	TR << "Idealizing doing nothing.  Not currently implemented. " << std::endl;

}

void
GlycanSampler::setup_movers(
	core::pose::Pose & pose,
	utility::vector1< bool > const & dihedral_subset,
	utility::vector1< bool > const & sugar_bb_subset,
	utility::vector1< bool > const & subset)
{
	using namespace utility::pointer;

	ReturnResidueSubsetSelectorOP return_subset_dihedral = utility::pointer::make_shared< ReturnResidueSubsetSelector >( dihedral_subset);
	ReturnResidueSubsetSelectorOP return_subset_sugarbb =  utility::pointer::make_shared< ReturnResidueSubsetSelector >( sugar_bb_subset);
	ReturnResidueSubsetSelectorOP return_subset = utility::pointer::make_shared< ReturnResidueSubsetSelector >(subset);
	MoveMapOP min_mover_movemap = core::pose::carbohydrates::create_glycan_movemap_from_residue_selector(
		pose,
		return_subset,
		true /*chi move */,
		min_rings_,
		cartmin_);

	// Make sure we only allow symmetrical degrees of freedom to move and convert the local_movemap
	// to a local movemap
	if ( core::pose::symmetry::is_symmetric( pose )  )  {
		core::pose::symmetry::make_symmetric_movemap( pose, *min_mover_movemap );
	}

	//Dihedral masks to tell out samplers which torsions each residue has.
	std::map< core::Size, utility::vector1< core::Size >> dihedral_mask;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		dihedral_mask[ i ];
		if ( pose.residue( i ).is_carbohydrate() ) {
			for ( core::Size x = 1; x <= get_n_glycosidic_torsions_in_res(pose, i); ++x ) {
				dihedral_mask[ i ].push_back( x );
			}
		} else {
			dihedral_mask[ i ].push_back( 1 );
			dihedral_mask[ i ].push_back( 2 );
		}
	}

	//Totally Reset the weighted random mover
	weighted_random_mover_ = utility::pointer::make_shared< RandomMover >();

	//Inner cycles for simple bb movers
	core::Size total_inner_cycles = inner_ncycles_ * total_glycan_residues_;

	//////        //////
	// Linkage Movers //
	//////        //////
	linkage_mover_ = utility::pointer::make_shared< LinkageConformerMover >( return_subset_dihedral );
	linkage_mover_->set_use_conformer_population_stats(population_based_conformer_sampling_); //Uniform sampling of conformers.
	linkage_mover_->set_use_gaussian_sampling(use_gaussian_sampling_);

	weighted_random_mover_->add_mover(linkage_mover_, .20);

	//////        //////
	// SugarBB Movers //
	//////        //////
	core::Real sugar_bb_offset = 0;
	if ( count_selected( sugar_bb_subset ) != 0 ) {
		BBDihedralSamplerMoverOP sugar_sampler_mover = utility::pointer::make_shared< BBDihedralSamplerMover >();
		SugarBBSamplerOP phi_sugar_sampler = utility::pointer::make_shared< SugarBBSampler >( core::id::phi_dihedral );
		SugarBBSamplerOP psi_sugar_sampler = utility::pointer::make_shared< SugarBBSampler >( core::id::psi_dihedral );
		SugarBBSamplerOP omega_sugar_sampler = utility::pointer::make_shared< SugarBBSampler >( core::id::omega_dihedral);

		sugar_sampler_mover->add_sampler( phi_sugar_sampler );
		sugar_sampler_mover->add_sampler( psi_sugar_sampler );
		sugar_sampler_mover->add_sampler( omega_sugar_sampler );

		sugar_sampler_mover->set_residue_selector( return_subset_sugarbb );
		sugar_sampler_mover->set_dihedral_mask(dihedral_mask);

		core::Real sugar_bb_sampler_weight;
		if ( tree_based_min_pack_ /* default true */ ) {

			sugar_bb_sampler_weight = .30;
		} else {
			sugar_bb_sampler_weight = .40;
		}
		//If we are doing layer 1 that has no sugar bb data, we add an offset so that the linkage conformer mover
		// does not dominate.  This will be better when we derive sugarbb data from the PDB and make the conformers better.
		TR.Debug << " Adding sugar bb mover. " << std::endl;

		//We downweight sugarbb to do shear movements.
		if ( use_shear_ ) {
			sugar_bb_sampler_weight = sugar_bb_sampler_weight - .10;
		}

		if ( total_inner_cycles != 0 ) {
			RationalMonteCarloOP mc_cycle = make_shared< RationalMonteCarlo >(
				sugar_sampler_mover,
				scorefxn_,
				total_inner_cycles,
				kt_,
				false /*recover_low*/);

			weighted_random_mover_->add_mover(mc_cycle, sugar_bb_sampler_weight);
		} else {
			weighted_random_mover_->add_mover(sugar_sampler_mover, sugar_bb_sampler_weight);
		}
	} else {
		TR << "Increasing SmallMover weights as there is no data for SugarBB sampler to use..." << std::endl;
		sugar_bb_offset = .30;
	}

	//////        //////
	// Shear Movers //
	//////        //////
	shear_ = make_shared< ShearMover >();
	shear_->set_residue_selector(return_subset_dihedral);
	if ( use_shear_ ) {
		if ( total_inner_cycles != 0 ) {
			RationalMonteCarloOP mc_cycle_shear = make_shared< RationalMonteCarlo >(
				shear_,
				scorefxn_,
				total_inner_cycles,
				kt_,
				false /*recover_low*/);

			weighted_random_mover_->add_mover(mc_cycle_shear, .10);
		} else {
			weighted_random_mover_->add_mover(shear_, .10);
		}
	}

	//////        //////
	// SmallBB Movers //
	//////        //////
	BBDihedralSamplerMoverOP glycan_small_mover = utility::pointer::make_shared< BBDihedralSamplerMover >();
	BBDihedralSamplerMoverOP glycan_medium_mover= utility::pointer::make_shared< BBDihedralSamplerMover >();
	BBDihedralSamplerMoverOP glycan_large_mover = utility::pointer::make_shared< BBDihedralSamplerMover >();

	glycan_small_mover->set_residue_selector(  return_subset_dihedral );
	glycan_medium_mover->set_residue_selector(  return_subset_dihedral );
	glycan_large_mover->set_residue_selector(  return_subset_dihedral );

	glycan_small_mover->set_dihedral_mask( dihedral_mask );
	glycan_medium_mover->set_dihedral_mask( dihedral_mask );
	glycan_large_mover->set_dihedral_mask( dihedral_mask );

	core::Size max_glycan_dihedrals = 0;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( dihedral_subset[ i ] ) {
			core::Size n_dihedrals = get_n_glycosidic_torsions_in_res( pose, i );
			if ( max_glycan_dihedrals < n_dihedrals ) max_glycan_dihedrals = n_dihedrals;
		}
	}
	TR.Trace << "Max Dihedrals: " << max_glycan_dihedrals << std::endl;
	for ( core::Size i =1; i <= max_glycan_dihedrals; ++i ) {
		auto dih_type = static_cast< core::id::MainchainTorsionType >( i );

		SmallBBSamplerOP small_sampler = utility::pointer::make_shared< SmallBBSampler >( dih_type, 30 ); // +/- 15 degrees
		SmallBBSamplerOP medium_sampler= utility::pointer::make_shared< SmallBBSampler >( dih_type, 90 ); // +/- 45 degrees
		SmallBBSamplerOP large_sampler = utility::pointer::make_shared< SmallBBSampler >( dih_type, 180); // +/- 90 degrees

		glycan_small_mover->add_sampler( small_sampler );
		glycan_medium_mover->add_sampler( medium_sampler );
		glycan_large_mover->add_sampler( large_sampler );
	}

	///Wrap movers in RationalMonteCarlo objects so that we may call them a few times since they are so fast.


	//
	//BB Sampling:
	//  4:2:1 ratio of sampling Small,Medium,Large
	//
	//  .3 = 4X (small) + 2X (medium) + X(large)

	core::Real total_prob = .3;
	constexpr core::Real  X = .3/(1+2+4);

	core::Real   small_prob = 4*X + sugar_bb_offset/3;
	core::Real  medium_prob = 2*X + sugar_bb_offset/3;
	core::Real   large_prob =   X + sugar_bb_offset/3;

	//Initialize inner MC movers
	if ( total_inner_cycles != 0 ) {
		RationalMonteCarloOP mc_cycle_small = make_shared< RationalMonteCarlo >(
			glycan_small_mover,
			scorefxn_,
			total_inner_cycles,
			kt_,
			false /*recover_low*/);

		RationalMonteCarloOP mc_cycle_medium = make_shared< RationalMonteCarlo >(
			glycan_medium_mover,
			scorefxn_,
			total_inner_cycles,
			kt_,
			false /*recover_low*/);

		RationalMonteCarloOP mc_cycle_large = make_shared< RationalMonteCarlo >(
			glycan_large_mover,
			scorefxn_,
			total_inner_cycles,
			kt_,
			false /*recover_low*/);

		//If we are refining, do not make huge random movements.
		if ( ! refine_ ) {
			weighted_random_mover_->add_mover( mc_cycle_small,  small_prob);
			weighted_random_mover_->add_mover( mc_cycle_medium,  medium_prob);
			weighted_random_mover_->add_mover( mc_cycle_large, large_prob );
		} else {
			weighted_random_mover_->add_mover( mc_cycle_small, total_prob);
		}
	} else {
		if ( ! refine_ ) {
			weighted_random_mover_->add_mover( glycan_small_mover,  small_prob);
			weighted_random_mover_->add_mover( glycan_medium_mover,  medium_prob);
			weighted_random_mover_->add_mover( glycan_large_mover, large_prob );
		} else {
			weighted_random_mover_->add_mover( glycan_small_mover, total_prob);
		}
	}

	//////    //////
	// Min Movers //
	//////    //////
	min_mover_ = utility::pointer::make_shared< MinMover >( min_mover_movemap->clone(), scorefxn_, "dfpmin_armijo_nonmonotone", 0.01, true /* use_nblist*/ );

	if ( (! option [OptionKeys::run::nblist_autoupdate].user()) && (! option [OptionKeys::run::nblist_autoupdate]() ) ) {
		min_mover_->min_options()->nblist_auto_update( true );
	}

	if ( cartmin_ ) {
		min_mover_->min_type("lbfgs_armijo_nonmonotone");
		min_mover_->cartesian(true);
		min_mover_->min_options()->max_iter(200); //Great suggestion from Patrick Conway
	}

	core::Real min_pack_weight;
	if ( tree_based_min_pack_ ) {
		min_pack_weight = .10;

		GlycanTreeMinMoverOP tree_min_mover = utility::pointer::make_shared< GlycanTreeMinMover >();
		tree_min_mover->set_minmover(min_mover_);
		tree_min_mover->set_min_rings(min_rings_);
		tree_min_mover->set_residue_selector(return_subset);

		weighted_random_mover_->add_mover(tree_min_mover, min_pack_weight);
	} else {
		min_pack_weight = .05;
		weighted_random_mover_->add_mover(min_mover_, min_pack_weight);
	}
}

void
GlycanSampler::setup_packer(
	core::pose::Pose & pose,
	utility::vector1< bool > const & full_subset)
{
	if ( ! tf_ ) {
		setup_default_task_factory(full_subset, pose);
	}

	core::pack::task::TaskFactoryOP full_task = get_all_glycans_and_neighbor_res_task_factory(full_subset);

	//Setup pack rotamers mover, add it.
	PackRotamersMoverOP tree_packer;
	packer_     = utility::pointer::make_shared< PackRotamersMover >();
	tree_packer = utility::pointer::make_shared< PackRotamersMover >();
	packer_->score_function(scorefxn_);
	tree_packer->score_function(scorefxn_);

	packer_->task_factory(full_task);
	tree_packer->task_factory(tf_);

	core::Real min_pack_weight = 0;
	if ( tree_based_min_pack_ ) {
		min_pack_weight = .10;
	} else {
		min_pack_weight = .05;
	}
	weighted_random_mover_->add_mover( tree_packer, min_pack_weight );
}



///@brief Initialize all objects.  Called at apply time!
void
GlycanSampler::init_objects(core::pose::Pose & pose ){


	using namespace core::pose::carbohydrates;
	using namespace protocols::moves;
	using namespace protocols::simple_moves;
	using namespace core::kinematics;

	TR << "initializing objects " << std::endl;

	final_residue_subset_.clear();

	setup_score_function();
	scorefxn_->score(pose);

	if ( ! selector_ ) selector_ = utility::pointer::make_shared< GlycanResidueSelector >();

	// Make sure our selector is symmetrical.  We don't extra non-symmetric DOFs moving (BBSampler/SugarBBSampler/Conformer) machinery.
	utility::vector1< bool > subset = selector_->apply(pose );
	if ( core::pose::symmetry::is_symmetric( pose )  )  {
		subset = core::select::get_master_subunit_selection(pose, subset);
	}

	utility::vector1< bool > dihedral_subset = subset; //Start of tree has no dihedrals to move.
	utility::vector1< bool > sugar_bb_subset = subset; //SugarBB not for residues attached to protein.


	//Setup SugarBB subset and filter the subset.
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( subset[ i ] ) {
			if ( ! pose.residue( i ).is_carbohydrate() ) {
				utility_exit_with_message(" GlycanSampler: Residue "+utility::to_string(i)+" set in ResidueSelector, but not a carbohdyrate residue!");
			}
			core::Size parent = pose.glycan_tree_set()->get_parent( i );
			//Residue has no dihedrals.  So turn it off so we don't choose it and then do nothing with it.
			if ( parent == 0 ) {
				dihedral_subset[ i ] = false;
				sugar_bb_subset[ i ] = false;
			} else if ( ! pose.residue( parent).is_carbohydrate() ) {
				//Residue is the first residue of the tree, attached to protein.
				//SugarBB not for residues attached to protein.
				sugar_bb_subset[i] = false;
			}
		}
	}

	if ( ! refine_ && randomize_first_ ) randomize_glycan_torsions( pose, subset);


	total_glycan_residues_ = count_selected(subset);
	final_residue_subset_ = subset;

	TR << "Modeling " << total_glycan_residues_ << " glycan residues" << std::endl;

	if ( total_glycan_residues_ == 0 ) {
		utility_exit_with_message( " No glycan residues in pose.  Cannot continue. ");
	}

	setup_movers( pose, dihedral_subset, sugar_bb_subset, subset);
	setup_packer( pose, subset );
}

void
GlycanSampler::apply( core::pose::Pose& pose ){
	using namespace core::kinematics;
	using namespace core::chemical::carbohydrates;
	using namespace core::scoring;
	using namespace protocols::moves;
	using namespace protocols::minimization_packing;
	using utility::to_string;

	accept_log_.clear();

	init_objects( pose );

	PyMOLMover pmm_accepts = PyMOLMover();
	PyMOLMover pmm_trials  = PyMOLMover();
	pmm_accepts.keep_history( true );
	pmm_trials.keep_history( true );
	pmm_accepts.set_PyMOL_model_name( "accepts_"+ pmm_accepts.get_PyMOL_model_name( pose ));
	pmm_trials.set_PyMOL_model_name(  "trials_" + pmm_trials.get_PyMOL_model_name( pose ));

	TR << "Initialized" << std::endl;
	utility::vector1< core::Size > bb_residues;

	core::Real energy = scorefxn_->score( pose );
	TR << "starting energy: "<< energy << std::endl;

	core::Size total_rounds = 0;
	if ( forced_total_rounds_ == 0 ) {
		total_rounds = total_glycan_residues_ * rounds_;
		TR << "Total Rounds = "<< total_rounds << " ( " << total_glycan_residues_ << " glycans * " << rounds_ << " )"<<std::endl;
	} else if ( match_sampling_of_modeler_ ) {
		//Matches sampling to that of the GlycanTreeSampler with default settings (window_size 2, overlap 1)
		// Used for benchmarking

		//R*nL0 + R*nLm1 + R*(N-(nL0+nL1)*2
		//
		//R = Number of set rounds
		//N = Number of total glycans
		//nL0 = Number of glycan residues in layer 0
		//nLm1 = Number of glycan residues in last Layer (-1 index)

		//Setup
		GlycanLayerSelector layer_selector = GlycanLayerSelector();
		layer_selector.set_layer(0, 0);
		utility::vector1< bool > const layer_0 = layer_selector.apply(pose);

		core::Size max_end_layer = pose.glycan_tree_set()->get_largest_glycan_tree_layer( final_residue_subset_ );
		layer_selector.set_layer(max_end_layer, max_end_layer);
		utility::vector1< bool > const layer_m1  = layer_selector.apply(pose);

		//Variables
		core::Size const R = rounds_;
		core::Size const N = total_glycan_residues_;
		core::Size const nL0  = count_selected( AND_combine( final_residue_subset_, layer_0 ) );
		core::Size const nLm1 = count_selected( AND_combine( final_residue_subset_, layer_m1) );

		//Calculation
		total_rounds = (R * nL0) + (R * nLm1) + (2 * R * (N - (nL0+nLm1) ));
	} else {
		total_rounds = forced_total_rounds_;
		TR << "Custom total rounds: " << total_rounds << std::endl;
	}


	bool accepted = false;
	mc_ = utility::pointer::make_shared< MonteCarlo >(pose, *scorefxn_, kt_);
	mc_->set_last_accepted_pose(pose);
	core::Real energy_pre_move = 0;

	//pose.dump_pdb("Pre-relax-structure.pdb");

	for ( core::Size round = 1; round <= total_rounds; ++round ) {

		if ( round % 10 == 0 ) {
			TR << "Round: "<< round << std::endl;
		}
		energy_pre_move = scorefxn_->score(pose);
		weighted_random_mover_->apply(pose);
		if ( weighted_random_mover_->get_last_move_status() == protocols::moves::MS_SUCCESS ) {
			core::Real energy2 = scorefxn_->score(pose);
			if ( TR.Debug.visible() ) {
				TR.Debug << "energy pre- move: "<< energy_pre_move << std::endl;
				TR.Debug << "energy post move: "<< energy2 << std::endl;
			}
			//TR << "energy pre- move  "<< energy_pre_move << std::endl;
			//TR << "energy post move: "<< energy2 << std::endl;

			if ( pymol_movie_ ) {
				pmm_trials.apply( pose );
			}
			accepted = mc_->boltzmann( pose );

			if ( TR.Debug.visible() ) TR.Debug << "Accepted? " << accepted << std::endl;
			if ( pymol_movie_ && accepted ) {
				pmm_accepts.apply( pose );
			}
			std::string out = to_string( round )+" "+to_string( energy2 )+" "+to_string( accepted );
			accept_log_.push_back( out );

			out = "FINAL "+ to_string( round )+" "+to_string( scorefxn_->score( pose ));
			accept_log_.push_back( out );
		} else {
			TR << "Last mover failed.  Continueing!" << std::endl;
			continue;
		}
	}
	mc_->recover_low( pose );

	TR << "energy start: "<< energy << std::endl;
	energy = scorefxn_->score( pose );
	TR << "energy final: "<< energy << std::endl;

	if ( final_min_ ) {

		run_shear_min_pack(*min_mover_, *packer_, *shear_, *mc_, total_glycan_residues_, pose, use_shear_);
		run_shear_min_pack(*min_mover_, *packer_, *shear_, *mc_, total_glycan_residues_, pose, use_shear_);

		min_mover_->apply( pose );
		energy = scorefxn_->score(pose);
		TR << "energy final post min: " << energy << std::endl;
		mc_->boltzmann( pose );
		mc_->recover_low( pose);
		energy = scorefxn_->score(pose);
		TR << "energy recov post min: " << energy << std::endl;
	} else {
		energy = scorefxn_->score( pose );
		TR << "energy final: "<< energy << std::endl;
	}

	std::string out = "FINAL END "+to_string( scorefxn_->score(pose) );
	accept_log_.push_back(out);

	//Add Accept log to pose, have it written:
	for ( core::Size i = 1; i <= accept_log_.size(); ++i ) {
		core::pose::add_comment(pose, "ACCEPT LOG "+utility::to_string(i), accept_log_[i]);
	}
	mc_->show_counters();
}

/////////////// Creator ///////////////


std::string GlycanSampler::get_name() const {
	return mover_name();
}

std::string GlycanSampler::mover_name() {
	return "GlycanSampler";
}


std::string GlycanSamplerCreator::keyname() const {
	return GlycanSampler::mover_name();
}

protocols::moves::MoverOP
GlycanSamplerCreator::create_mover() const {
	return utility::pointer::make_shared< GlycanSampler >();
}

void GlycanSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycanSampler::provide_xml_schema( xsd );
}


} //protocols
} //carbohydrates
