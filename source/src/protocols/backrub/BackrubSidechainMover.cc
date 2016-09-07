// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backrub/BackrubSidechainMover.cc
/// @brief BackrubSidechainMover methods implemented
/// @author

// Unit Headers
#include <protocols/backrub/BackrubSidechainMover.hh>
#include <protocols/backrub/BackrubSidechainMoverCreator.hh>

// Package Headers

// Project Headers
#include <protocols/backrub/BackrubSegment.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <numeric/MultiDimensionalHistogram.fwd.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

//Auto Headers
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/MultiDimensionalHistogram.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/id/TorsionID_Range.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.backrub.BackrubSidechainMover" );

namespace protocols {
namespace backrub {

std::string
protocols::backrub::BackrubSidechainMoverCreator::keyname() const {
	return BackrubSidechainMoverCreator::mover_name();
}

protocols::moves::MoverOP
protocols::backrub::BackrubSidechainMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BackrubSidechainMover );
}

std::string
protocols::backrub::BackrubSidechainMoverCreator::mover_name() {
	return "BackrubSidechain";
}

protocols::backrub::BackrubSidechainMover::BackrubSidechainMover(
) :
	protocols::canonical_sampling::ThermodynamicMover(),
	backrub_mover_(protocols::backrub::BackrubMoverOP( new protocols::backrub::BackrubMover )),
	sidechain_mover_(protocols::simple_moves::sidechain_moves::SidechainMoverOP( new protocols::simple_moves::sidechain_moves::SidechainMover )),
	record_statistics_(false),
	statistics_filename_("brsc_stats.txt")
{
	backrub_mover_->set_min_atoms(7);
	backrub_mover_->set_max_atoms(7);
}

protocols::backrub::BackrubSidechainMover::BackrubSidechainMover(
	BackrubSidechainMover const & mover
) :
	//utility::pointer::ReferenceCount(),
	protocols::canonical_sampling::ThermodynamicMover(mover),
	valid_segments_(mover.valid_segments_),
	last_valid_segment_index_(mover.last_valid_segment_index_),
	last_chi1_pre_(mover.last_chi1_pre_),
	last_chi1_post_(mover.last_chi1_post_),
	record_statistics_(mover.record_statistics_),
	statistics_filename_(mover.statistics_filename_),
	proposal_hists_(mover.proposal_hists_),
	accept_hists_(mover.accept_hists_)
{
	if ( mover.backrub_mover_ ) {
		backrub_mover_ = utility::pointer::dynamic_pointer_cast< protocols::backrub::BackrubMover > ( mover.backrub_mover_->clone() );
		runtime_assert(backrub_mover_ != nullptr);
	}
	if ( mover.sidechain_mover_ ) {
		sidechain_mover_ = utility::pointer::dynamic_pointer_cast< protocols::simple_moves::sidechain_moves::SidechainMover > ( mover.sidechain_mover_->clone() );
		runtime_assert(sidechain_mover_ != nullptr);
	}
}

protocols::backrub::BackrubSidechainMover::~BackrubSidechainMover()= default;

protocols::moves::MoverOP
protocols::backrub::BackrubSidechainMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::backrub::BackrubSidechainMover( *this ) );
}

protocols::moves::MoverOP
protocols::backrub::BackrubSidechainMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new BackrubSidechainMover );
}

std::string
protocols::backrub::BackrubSidechainMover::get_name() const
{
	return "BackrubSidechainMover";
}

void
protocols::backrub::BackrubSidechainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & pose
)
{
	if ( tag->hasOption("pivot_residues") ) {
		set_pivot_residues(core::pose::get_resnum_list(tag, "pivot_residues", pose));
	}

	core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );

	if ( tag->hasOption("task_operations") ) {

		std::string const t_o_val( tag->getOption<std::string>("task_operations") );
		typedef utility::vector1< std::string > StringVec;
		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		 for ( auto const & t_o_key : t_o_keys ) {
			if ( data.has( "task_operations", t_o_key ) ) {
				new_task_factory->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", t_o_key ) );
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
			}
		}

	} else {

		new_task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	}

	set_task_factory(new_task_factory);

	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", prob_uniform() ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", prob_withinrot() ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", prob_random_pert_current() ) );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
	set_require_mm_bend( tag->getOption<bool>( "require_mm_bend", require_mm_bend() ) );
	set_record_statistics( tag->getOption<bool>( "record_statistics", record_statistics() ) );
	set_statistics_filename( tag->getOption<std::string>( "statistics_filename", statistics_filename() ) );

	update_segments(pose);
}

void
protocols::backrub::BackrubSidechainMover::update_segments(
	core::pose::Pose const & pose
)
{
	if ( !(backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree()) ) {
		backrub_mover_->set_input_pose(core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose) ) ));
	}

	backrub_mover_->set_input_pose(core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose) ) ));
	backrub_mover_->clear_segments();
	backrub_mover_->add_mainchain_segments();

	sidechain_mover_->init_task(pose);

	valid_segments_.clear();

	for ( core::Size i = 1; i <= backrub_mover_->num_segments(); ++i ) {
		BackrubSegment const & segment(backrub_mover_->segment(i));
		if ( pose.residue(segment.start_atomid().rsd()).atom_name(segment.start_atomid().atomno()) == " CA " &&
				pose.residue(segment.end_atomid().rsd()).atom_name(segment.end_atomid().atomno()) == " CA " &&
				segment.size() == 7 ) {
			core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;
			if ( sidechain_mover_->residue_packed()[middle_rsd] ) {
				//TR << "Adding segment start:" << segment.start_atomid() << "end:"
				//  << segment.end_atomid() << "size: " << segment.size() << std::endl;
				valid_segments_.push_back(i);
			}
		}
	}

	if ( record_statistics_ ) reset_statistics();
}

void
protocols::backrub::BackrubSidechainMover::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
)
{
	if ( !(valid_segments_.size() && backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree()) ) {
		update_segments(pose);
	}

	// because the segments have already been set up, they shouldn't be set up again in the next call
	backrub_mover_->initialize_simulation(pose, metropolis_hastings_mover,cycle);
	sidechain_mover_->initialize_simulation(pose, metropolis_hastings_mover,cycle);

	if ( record_statistics_ ) reset_statistics();
}

void
protocols::backrub::BackrubSidechainMover::apply(
	core::pose::Pose & pose
)
{
	if ( !(valid_segments_.size() && backrub_mover_->get_input_pose() && backrub_mover_->get_input_pose()->fold_tree() == pose.fold_tree()) ) {
		update_segments(pose);
	}

	last_valid_segment_index_ = numeric::random::rg().random_range(1, valid_segments_.size());
	BackrubSegment const & segment(backrub_mover_->segment(valid_segments_[last_valid_segment_index_]));
	core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;

	last_chi1_pre_ = numeric::conversions::radians(pose.chi(1, middle_rsd));

	sidechain_mover_->next_resnum(middle_rsd);
	sidechain_mover_->apply(pose);

	last_chi1_post_ = numeric::conversions::radians(pose.chi(1, middle_rsd));

	backrub_mover_->set_next_segment_id(valid_segments_[last_valid_segment_index_]);
	backrub_mover_->apply(pose);

	update_type();
}

void
protocols::backrub::BackrubSidechainMover::observe_after_metropolis(
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	if ( record_statistics_ ) record_histograms(metropolis_hastings_mover.monte_carlo()->mc_accepted());
}

void
protocols::backrub::BackrubSidechainMover::finalize_simulation(
	core::pose::Pose & /*pose*/,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	if ( record_statistics_ ) {

		std::ostringstream filename;
		if ( metropolis_hastings_mover.output_name() != "" ) {
			filename << metropolis_hastings_mover.output_name() << "_";
		}
		filename << statistics_filename();

		utility::io::ozstream statistics_stream(filename.str());
		output_statistics(statistics_stream);
		statistics_stream.close();
	}
}

utility::vector1<core::Size> const &
protocols::backrub::BackrubSidechainMover::pivot_residues() const
{
	return backrub_mover_->pivot_residues();
}

void
protocols::backrub::BackrubSidechainMover::set_pivot_residues(
	utility::vector1<core::Size> const & pivot_residues
)
{
	backrub_mover_->set_pivot_residues(pivot_residues);
}

core::pack::task::TaskFactoryCOP
protocols::backrub::BackrubSidechainMover::task_factory() const
{
	return sidechain_mover_->task_factory();
}

void
protocols::backrub::BackrubSidechainMover::set_task_factory(
	core::pack::task::TaskFactoryOP task_factory
)
{
	task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	sidechain_mover_->set_task_factory(task_factory);
}

core::Real
protocols::backrub::BackrubSidechainMover::prob_uniform() const
{
	return sidechain_mover_->prob_uniform();
}

void
protocols::backrub::BackrubSidechainMover::set_prob_uniform(
	core::Real prob_uniform
)
{
	sidechain_mover_->set_prob_uniform(prob_uniform);
}

core::Real
protocols::backrub::BackrubSidechainMover::prob_withinrot() const
{
	return sidechain_mover_->prob_withinrot();
}

void
protocols::backrub::BackrubSidechainMover::set_prob_withinrot(
	core::Real prob_withinrot
)
{
	sidechain_mover_->set_prob_withinrot(prob_withinrot);
}

core::Real
protocols::backrub::BackrubSidechainMover::prob_random_pert_current() const
{
	return sidechain_mover_->prob_random_pert_current();
}

void
protocols::backrub::BackrubSidechainMover::set_prob_random_pert_current(
	core::Real prob_pert
)
{
	sidechain_mover_->set_prob_random_pert_current(prob_pert);
}

bool
protocols::backrub::BackrubSidechainMover::preserve_detailed_balance() const
{
	return backrub_mover_->preserve_detailed_balance() && sidechain_mover_->preserve_detailed_balance();
}

void
protocols::backrub::BackrubSidechainMover::set_preserve_detailed_balance(
	bool preserve_detailed_balance
)
{
	backrub_mover_->set_preserve_detailed_balance(preserve_detailed_balance);
	sidechain_mover_->set_preserve_detailed_balance(preserve_detailed_balance);
}

bool
protocols::backrub::BackrubSidechainMover::require_mm_bend() const
{
	return backrub_mover_->require_mm_bend();
}

void
protocols::backrub::BackrubSidechainMover::set_require_mm_bend(
	bool require_mm_bend
)
{
	backrub_mover_->set_require_mm_bend(require_mm_bend);
}

utility::vector1<core::id::TorsionID_Range>
protocols::backrub::BackrubSidechainMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

bool
protocols::backrub::BackrubSidechainMover::record_statistics() const
{
	return record_statistics_;
}

void
protocols::backrub::BackrubSidechainMover::set_record_statistics(
	bool record_statistics
)
{
	bool const needs_reset(record_statistics_ != record_statistics);
	record_statistics_ = record_statistics;
	if ( needs_reset ) reset_statistics();
}

std::string const &
protocols::backrub::BackrubSidechainMover::statistics_filename() const
{
	return statistics_filename_;
}

void
protocols::backrub::BackrubSidechainMover::set_statistics_filename(
	std::string const & statistics_filename
)
{
	statistics_filename_ = statistics_filename;
}

void
protocols::backrub::BackrubSidechainMover::reset_statistics()
{
	setup_histograms();
}

void
protocols::backrub::BackrubSidechainMover::output_statistics(
	std::ostream & out
)
{
	for ( core::Size i = 1; i <= proposal_hists_.size(); ++i ) {
		out << proposal_hists_[i] << accept_hists_[i];
	}
}

void
protocols::backrub::BackrubSidechainMover::setup_histograms()
{
	if ( !(record_statistics_ && backrub_mover_->get_input_pose()) ) {
		proposal_hists_.resize(0);
		accept_hists_.resize(0);
		return;
	}

	core::pose::Pose const & pose(*backrub_mover_->get_input_pose());

	core::Real max_angle_disp(ceil(numeric::conversions::degrees(backrub_mover_->max_angle_disp_7())));
	core::Size const num_bins(static_cast<core::Size>(2*max_angle_disp));
	numeric::conversions::to_radians(max_angle_disp);

	proposal_hists_.resize(valid_segments_.size());
	accept_hists_.resize(valid_segments_.size());

	for ( core::Size i = 1; i <= valid_segments_.size(); ++i ) {

		core::Size res_num(backrub_mover_->segment(valid_segments_[i]).start_atomid().rsd()+1);

		std::ostringstream res_label_stream;
		res_label_stream << pose.residue(res_num).name3() << " "
			<< pose.pdb_info()->chain(res_num) << " "
			<< pose.pdb_info()->number(res_num);

		std::ostringstream proposal_label_stream;
		proposal_label_stream << res_label_stream.str() << " Proposal";
		proposal_hists_[i].label(proposal_label_stream.str());
		proposal_hists_[i].reset_counts();
		proposal_hists_[i].num_dimensions(3);
		proposal_hists_[i].set_dimension(1, num_bins, -max_angle_disp, max_angle_disp, "backrub_disp");
		proposal_hists_[i].set_dimension(2, 3, 0, numeric::constants::r::pi_2, "chi1_pre");
		proposal_hists_[i].set_dimension(3, 3, 0, numeric::constants::r::pi_2, "chi1_post");

		std::ostringstream accept_label_stream;
		accept_label_stream << res_label_stream.str() << " Accept";
		accept_hists_[i].label(accept_label_stream.str());
		accept_hists_[i].reset_counts();
		accept_hists_[i].num_dimensions(3);
		accept_hists_[i].set_dimension(1, num_bins, -max_angle_disp, max_angle_disp, "backrub_disp");
		accept_hists_[i].set_dimension(2, 3, 0, numeric::constants::r::pi_2, "chi1_pre");
		accept_hists_[i].set_dimension(3, 3, 0, numeric::constants::r::pi_2, "chi1_post");
	}
}

void
protocols::backrub::BackrubSidechainMover::record_histograms(
	bool accepted
)
{
	if ( proposal_hists_.size() == 0 ) setup_histograms();

	//BackrubSegment const & segment(backrub_mover_->segment(backrub_mover_->last_segment_id()));
	//core::Size middle_rsd = (segment.start_atomid().rsd() + segment.end_atomid().rsd())/2;
	//TR << "Changed residue " << middle_rsd << " chi1 from "
	//   << numeric::conversions::degrees(last_chi1_pre_) << " to "
	//   << numeric::conversions::degrees(last_chi1_post_) << " backrub angle "
	//   << numeric::conversions::degrees(backrub_mover_->last_angle()) << " accepted "
	//   << accepted << std::endl;

	utility::vector1<core::Real> values(3);
	values[1] = backrub_mover_->last_angle();
	values[2] = numeric::nonnegative_principal_angle_radians(last_chi1_pre_);
	values[3] = numeric::nonnegative_principal_angle_radians(last_chi1_post_);

	proposal_hists_[last_valid_segment_index_].record(values);
	if ( accepted ) accept_hists_[last_valid_segment_index_].record(values);
}

void
protocols::backrub::BackrubSidechainMover::update_type()
{
	std::stringstream mt;

	char bin_letters[] = {'p', 't', 'm'};

	core::Size chi1_pre_bin(static_cast<core::Size>(floor(numeric::nonnegative_principal_angle_radians(last_chi1_pre_)/numeric::constants::r::pi_2_over_3)));
	if ( chi1_pre_bin == 3 ) chi1_pre_bin = 2;
	core::Size chi1_post_bin(static_cast<core::Size>(floor(numeric::nonnegative_principal_angle_radians(last_chi1_post_)/numeric::constants::r::pi_2_over_3)));
	if ( chi1_post_bin == 3 ) chi1_post_bin = 2;

	mt << "brsc_" << bin_letters[chi1_pre_bin] << bin_letters[chi1_post_bin] << "_"
		<< (sidechain_mover_->last_nchi() ? (sidechain_mover_->last_uniform() ? "unif" : (sidechain_mover_->last_withinrot() ? "withinrot" : "rot")) : "none");

	std::string const new_type(mt.str());
	type(new_type);
}

} //moves
} //protocols

