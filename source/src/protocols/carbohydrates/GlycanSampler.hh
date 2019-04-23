// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanSampler.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_carbohydrates_GlycanSampler_hh
#define INCLUDED_protocols_carbohydrates_GlycanSampler_hh

// Unit headers
#include <protocols/carbohydrates/GlycanSampler.fwd.hh>
#include <protocols/carbohydrates/LinkageConformerMover.hh>
#include <protocols/moves/Mover.hh>


#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/minimization_packing/MinMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>



#include <basic/datacache/DataMap.fwd.hh>


namespace protocols {
namespace carbohydrates {

///@brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// Each round optimizes either one residue for BB sampling, linkage, or multiple for minimization.
/// Currently uses a random sampler with a set of weights to each mover for sampling.
///
/// Weights are currently as follows:
///  .40 Phi/Psi Sugar BB Sampling
///  .20 Linkage Conformer Sampling
///  .30 Small BB Sampling - equal weight to phi, psi, or omega
///    -> .17 +/- 15 degrees
///    -> .086 +/- 45 degrees
///    -> .044 +/- 90 degrees
///  .05 GlycanTreeMinMover
///  .05 PackRotamersMover
///
/// Supports Symmetry
///
class GlycanSampler : public protocols::moves::Mover {

public:

	GlycanSampler();

	//@brief constructor with arguments
	GlycanSampler( core::select::residue_selector::ResidueSelectorCOP selector,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size rounds = 75);

	// copy constructor
	GlycanSampler( GlycanSampler const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~GlycanSampler() override;

	void
	apply( core::pose::Pose & pose ) override;

public:


	///@brief Set the Movemap
	void
	set_residue_selector(core::select::residue_selector::ResidueSelectorCOP selector );

	///@brief Set a ResidueSelector for glycan residue selection, instead of the typical movemap.
	void
	set_selector(core::select::residue_selector::ResidueSelectorCOP selector);

	///@brief Set the TaskFactory to control side-chain packing of surrounding amino acids and the OH groups of the glycans.
	void
	set_taskfactory(core::pack::task::TaskFactoryCOP tf);

	///@brief Set the ScoreFunction
	void
	set_scorefunction( core::scoring::ScoreFunctionCOP scorefxn);

	///@brief Each round applys a random mover to pose.
	/// This setting is multiplied by the number of glycan residues in the movemap for the total number of rounds
	void
	set_rounds( core::Size rounds);

	///@brief Set how our Conformer mover samples.  Default is to do uniform sampling on the conformers instead of using the population as probabilities.
	void
	set_population_based_conformer_sampling(bool pop_based_sampling);

	///@brief Set whether if we are sampling torsions uniformly within an SD for the LinkageConformerMover (false)
	///  or sampling the gaussian (true).
	///  Default false
	void
	set_use_gaussian_sampling(bool use_gaussian_sampling);

	void
	set_kt( core::Real kt);

	void
	set_defaults();

	///@brief Use Cartesian minimization instead of Dihedral minimization.
	/// Default (for now) is dihedral.
	void
	use_cartmin( bool use_cartmin );

	///@brief Set refinement mode instead of denovo modeling mode.
	void
	set_refine( bool refine );

	///@brief set to minimize ring torsions during minimzation.  Default false.
	void
	set_min_rings( bool min_rings);

	///@brief Set the protocol to use the refactored Shear Mover for glycan torsions at 10 % probability.
	/// Default false.
	void
	set_use_shear( bool use_shear);

	///@brief Set the protocol to randomize torsions before beginning.
	/// This actually helps get to lower energy models.
	/// Default True.  If doing refinement, this is automatically turned off.
	void
	set_randomize_first( bool randomize_first );

	///@brief Set the number of inner cycles for BB sampling through small/sugarBB.
	/// This is multiplied by the number of glycan residues
	/// Default 1
	void
	set_inner_bb_cycles( core::Size inner_bb_cycles );

public:
	void
	show( std::ostream & output=std::cout ) const override;


	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;



	//GlycanSampler & operator=( GlycanSampler const & src );

	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	bool
	reinitialize_for_each_job() const override {
		return true;
	}

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	///@brief Randomize all torsions of the subset.  Used to start the protocol.
	void
	randomize_glycan_torsions( core::pose::Pose & pose, utility::vector1< bool > const & subset ) const;

	///@brief Attempt to idealize all residues in of a set of trees. Experimental!
	void
	idealize_glycan_residues( core::pose::Pose & pose, utility::vector1< core::Size > const & tree_subset ) const;


public:

	///@brief Get the number of glycan sampler rounds this class is set to run.
	core::Size
	get_glycan_sampler_rounds();

	///@brief This allows us to force a number of rounds instead of doing rounds*glycan residues
	void
	force_total_rounds( core::Size total_rounds );

private:

	void
	init_objects( core::pose::Pose & pose );

	void
	set_cmd_line_defaults();

	void
	apply_to_res(
		core::pose::Pose & pose,
		core::Size resnum,
		core::kinematics::MoveMapOP mm,
		core::scoring::ScoreFunctionOP score,
		core::Size round);

	void
	setup_default_task_factory(utility::vector1< bool > const & glycan_residues, core::pose::Pose const & pose );

	void
	setup_score_function();

	void
	setup_cartmin(core::scoring::ScoreFunctionOP scorefxn) const;


	//Setup the final WeightedMover from our subsets and movemap.
	void
	setup_movers(
		core::pose::Pose & pose,
		utility::vector1< bool > const & dihedral_subset,
		utility::vector1< bool > const & sugar_bb_subset,
		utility::vector1< bool > const & subset);

	void
	setup_packer(
		core::pose::Pose & pose,
		utility::vector1< bool > const & full_subset );

private:

	core::pack::task::TaskFactoryCOP tf_              = nullptr;

	moves::MonteCarloOP mc_                           = nullptr;
	core::scoring::ScoreFunctionCOP scorefxn_         = nullptr;

	LinkageConformerMoverOP linkage_mover_            = nullptr;
	moves::RandomMoverOP weighted_random_mover_       = nullptr;

	minimization_packing::MinMoverOP min_mover_       = nullptr;
	minimization_packing::PackRotamersMoverOP packer_ = nullptr;
	simple_moves::ShearMoverOP shear_                 = nullptr;

	core::Size rounds_ = 25; // cmdline
	core::Real kt_ = 2.0; // cmdline

	utility::vector1<std::string> accept_log_;

	bool test_ = false; // cmdline
	bool final_min_ = true; // cmdline
	bool refine_ = false; // cmdline

	core::Size total_glycan_residues_ = 0;
	bool pymol_movie_ = false; // cmdline

	utility::vector1< std::string > parsed_positions_;
	core::Real pack_distance_ = 6.0;
	bool cartmin_ = false; // cmdline
	bool tree_based_min_pack_ = true; // cmdline

	core::select::residue_selector::ResidueSelectorCOP selector_;
	bool population_based_conformer_sampling_ = false;
	bool use_gaussian_sampling_ = false;
	bool min_rings_ = false;

	core::Size forced_total_rounds_ = 0;
	bool use_shear_ = false;
	bool randomize_first_ = true;
	core::Size inner_ncycles_ = 0; //For individual bb movements, multiply this by n glycans.
	bool match_sampling_of_modeler_ = false; //For benchmarking
	utility::vector1< bool > final_residue_subset_;

};

std::ostream &operator<< (std::ostream &os, GlycanSampler const &mover);


} //protocols
} //carbohydrates


#endif //protocols/carbohydrates_GlycanSampler_hh







