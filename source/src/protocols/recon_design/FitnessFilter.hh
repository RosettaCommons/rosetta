// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/FitnessFilter.hh
/// @brief Returns the sum of energy of input poses. Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_recon_design_FitnessFilter_hh
#define INCLUDED_protocols_recon_design_FitnessFilter_hh

#include <protocols/recon_design/FitnessFilter.fwd.hh>
#include <protocols/filters/VectorPoseFilter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>

namespace protocols {
namespace recon_design {

/// @brief Returns the sum of energy of input poses. Only accessible through recon application.
class FitnessFilter : public filters::VectorPoseFilter
{
public:
	/// @brief empty constructor
	FitnessFilter();

	~FitnessFilter() override;

	/// @brief Movers derived from VectorPoseFilter must define two apply methods:
	/// apply and apply_mpi, depending on whether the protocol is run in MPI or not
	/// Running in MPI will distribute each pose to a separate node, so the mover will
	/// operate on one pose
	bool apply( core::pose::Pose const & pose ) const override;
	bool apply_mpi( core::pose::Pose const & pose ) const override;

	void report( std::ostream &, core::pose::Pose const & ) const override;

	core::Real report_sm( core::pose::Pose const & ) const override;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	filters::FilterOP clone() const override;
	filters::FilterOP fresh_instance() const override;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	/// @brief Internal functions to get the fitness over all poses
	core::Real calculate_fitness( core::pose::Pose const & ) const;
	core::Real calculate_fitness_mpi( core::pose::Pose const & pose ) const;

	core::Real threshold() const;
	void threshold( core::Real );

	core::scoring::ScoreFunctionOP sfxn() const;
	void sfxn( core::scoring::ScoreFunctionOP );

private:
	core::scoring::ScoreFunctionOP sfxn_;
	bool output_to_scorefile_;
	core::Real threshold_;
};


} //recon_design
} //protocols
#endif
