// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanRelaxMover.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@jhu.edu)


#ifndef INCLUDED_protocols_carbohydrates_GlycanRelaxMover_hh
#define INCLUDED_protocols_carbohydrates_GlycanRelaxMover_hh

// Unit headers
#include <protocols/carbohydrates/GlycanRelaxMover.fwd.hh>
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
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>

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
///  .05 MinMover
///  .05 PackRotamersMover
///
class GlycanRelaxMover : public protocols::moves::Mover {

public:

	GlycanRelaxMover();

	//@brief constructor with arguments
	GlycanRelaxMover( core::kinematics::MoveMapCOP mm,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size rounds = 75);

	// copy constructor
	GlycanRelaxMover( GlycanRelaxMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~GlycanRelaxMover() override;

	void
	apply( core::pose::Pose & pose ) override;

public:


	///@brief Set the Movemap
	void
	set_movemap(core::kinematics::MoveMapCOP movemap);

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

	///@brief Number of SDs to sample within during conformer sampling. Default is 2.0 SDs
	///@details Sample within X standard_deviations of the means when building [non-idealized] conformers
	void
	set_conformer_sampling_sd(core::Real conformer_sampling_sd);

	///@brief Set whether if we are sampling uniform within the set number of standard deviations or by uniform within the SD.
	/// Default True
	void
	set_uniform_sd_sampling(bool uniform_sd_sampling);

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



	//GlycanRelaxMover & operator=( GlycanRelaxMover const & src );

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
	setup_cartmin(core::scoring::ScoreFunctionOP scorefxn) const;

private:

	core::kinematics::MoveMapOP full_movemap_;
	core::kinematics::MoveMapOP glycan_movemap_;

	core::pack::task::TaskFactoryCOP tf_;

	moves::MonteCarloOP mc_;
	core::scoring::ScoreFunctionCOP scorefxn_;

	LinkageConformerMoverOP linkage_mover_;
	moves::RandomMoverOP weighted_random_mover_;

	simple_moves::MinMoverOP min_mover_;
	simple_moves::PackRotamersMoverOP packer_;

	core::Size rounds_ = 25; // cmdline
	core::Real kt_ = 2.0; // cmdline

	utility::vector1<std::string> accept_log_;

	bool test_ = false; // cmdline
	bool final_min_ = true; // cmdline
	bool refine_ = false; // cmdline

	core::Size total_glycan_residues_ = 0;
	bool pymol_movie_ = false; // cmdline

	std::string ref_pose_name_;
	bool use_branches_ = false;

	utility::vector1< std::string > parsed_positions_;
	core::Real pack_distance_ = 6.0;
	bool cartmin_ = false; // cmdline
	bool tree_based_min_pack_ = true; // cmdline

	core::select::residue_selector::ResidueSelectorCOP selector_;  //Residue selector to pass residues to relax.  Currently used for RS only.
	bool population_based_conformer_sampling_ = false;
	core::Real conformer_sd_ = 2.0;
	bool uniform_conformer_sd_sampling_ = true;

};

std::ostream &operator<< (std::ostream &os, GlycanRelaxMover const &mover);


} //protocols
} //carbohydrates


#endif //protocols/carbohydrates_GlycanRelaxMover_hh







