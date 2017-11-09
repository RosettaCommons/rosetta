// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/backrub/BackrubSidechainMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_backrub_BackrubSidechainMover_hh
#define INCLUDED_protocols_backrub_BackrubSidechainMover_hh

// Unit Headers
#include <protocols/backrub/BackrubSidechainMover.fwd.hh>

// Project Headers
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric Headers

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.fwd.hh>
#include <utility/vector0.hh>
#include <numeric/MultiDimensionalHistogram.fwd.hh>

#ifdef PYROSETTA
#include <numeric/MultiDimensionalHistogram.hh>
#endif


namespace protocols {
namespace backrub {

/// @details
class BackrubSidechainMover : public protocols::canonical_sampling::ThermodynamicMover {

public:

	/// @brief
	BackrubSidechainMover();

	/// @brief
	BackrubSidechainMover(
		BackrubSidechainMover const & mover
	);


	~BackrubSidechainMover() override;


	protocols::moves::MoverOP
	clone() const override;


	protocols::moves::MoverOP
	fresh_instance() const override;




	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

	void
	update_segments(
		core::pose::Pose const & pose
	);


	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	) override;


	void
	apply(
		core::pose::Pose & pose
	) override;

	/// @brief callback after the Metropolis criterion is evaluated

	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	) override;


	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	) override;

	//utility::vector1<core::Size> const &
	//pivot_residues() const;

	void
	set_pivot_residues(
		utility::vector1<core::Size> const & pivot_residues
	);

	void
	set_pivot_residue_selector(
		core::select::residue_selector::ResidueSelectorCOP pivot_residues
	);

	core::pack::task::TaskFactoryCOP
	task_factory() const;

	void
	set_task_factory(
		core::pack::task::TaskFactoryOP task_factory
	);

	/// @brief get the probability of uniformly sampling chi angles
	core::Real
	prob_uniform() const;

	/// @brief set the probability of uniformly sampling chi angles
	void
	set_prob_uniform(
		core::Real prob_uniform
	);

	/// @brief get the probability of sampling within the same rotamer
	core::Real
	prob_withinrot() const;

	/// @brief set the probability of sampling within the same rotamer
	void
	set_prob_withinrot(
		core::Real prob_withinrot
	);

	core::Real
	prob_random_pert_current() const;

	void
	set_prob_random_pert_current(
		core::Real prob_pert
	);

	/// @brief get whether detailed balance is preserved

	bool
	preserve_detailed_balance() const override;

	/// @brief set whether detailed balance is preserved

	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	) override;

	/// @brief get whether to exit during initialize_simulation if the mm_bend term isn't turned on
	bool
	require_mm_bend() const;

	/// @brief set whether to exit during initialize_simulation if the mm_bend term isn't turned on
	void
	set_require_mm_bend(
		bool require_mm_bend
	);

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges

	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	) override;

	/// @brief get whether to record mover statistics or not
	bool
	record_statistics() const;

	/// @brief set whether to record mover statistics or not
	void
	set_record_statistics(
		bool record_statistics
	);

	/// @brief get filename for statistics output
	std::string const &
	statistics_filename() const;

	/// @brief set filename for statistics output
	void
	set_statistics_filename(
		std::string const & statistics_filename
	);

	void
	reset_statistics();

	void
	output_statistics(
		std::ostream & out
	);

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
	setup_histograms();

	void
	record_histograms(
		bool accept
	);

	void
	update_type();

	protocols::backrub::BackrubMoverOP backrub_mover_;
	protocols::simple_moves::sidechain_moves::SidechainMoverOP sidechain_mover_;
	utility::vector1<core::Size> valid_segments_;
	core::Size last_valid_segment_index_;
	core::Real last_chi1_pre_;
	core::Real last_chi1_post_;
	bool record_statistics_;
	std::string statistics_filename_;
	utility::vector1<numeric::MultiDimensionalHistogram> proposal_hists_;
	utility::vector1<numeric::MultiDimensionalHistogram> accept_hists_;

}; //end BackrubSidechainMover

} //namespace backrub
} //namespace protocols

#endif // INCLUDED_protocols_backrub_BackrubSidechainMover_HH
