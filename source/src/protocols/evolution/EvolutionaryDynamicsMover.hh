// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/evolution/EvolutionaryDynamicsMover.hh
/// @brief perform a given mover and sample structures by MonteCarlo
/// @details The "score" evaluation of pose during MC after applying mover is done by
/// ither FilterOP that can do report_sm() or ScoreFunctionOP you gave.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_protocols_evolution_EvolutionaryDynamicsMover_hh
#define INCLUDED_protocols_evolution_EvolutionaryDynamicsMover_hh

// Unit header
#include <protocols/evolution/EvolutionaryDynamicsMover.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/moves/MonteCarloStatus.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered_map.hpp>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.fwd.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <protocols/moves/MoverApplyingMover.hh>
#include <protocols/monte_carlo/GenericMonteCarloMover.hh>

namespace protocols {
namespace evolution {

/// @brief Trigger API definition
typedef boost::function<bool(core::Size,
	core::Size,
	const core::pose::Pose&,
	core::scoring::ScoreFunctionOP)> EvolutionaryDynamicsMoverTrigger;


class EvolutionaryDynamicsMover : public protocols::monte_carlo::GenericMonteCarloMover {
	typedef protocols::monte_carlo::GenericMonteCarloMover Super;
public:
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef protocols::filters::FilterOP FilterOP;
	typedef protocols::moves::MoverOP MoverOP;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::PackerTaskOP PackerTaskOP;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef protocols::rosetta_scripts::ParsedProtocolOP ParsedProtocolOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::Movers_map Movers_map;

	/// @brief default constructor
	EvolutionaryDynamicsMover();

	/// @brief destructor
	~EvolutionaryDynamicsMover() override;

	/// @brief create copy constructor
	MoverOP clone() const override;

	/// @brief create this type of objectt
	MoverOP fresh_instance() const override;

	/// @brief apply EvolutionaryDynamicsMover (Mover)
	void apply( Pose & pose ) override;


	/// @brief core of MC -- evaulates a pose based on the scores/filters + temperatures. random_num is a vector of random numbers between 0 and 1 with size equal to the number of MC criteria
	bool boltzmann( Pose & pose, utility::vector1< core::Real > const & random_nums ) override;



public: // mutators
	void parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		Movers_map const & movers,
		Pose const &
	) override;

	void population_size( core::Real const s ){ population_size_ = s; };
	core::Real population_size() const{ return population_size_; }

	void disable_fitness_evaluation( bool const d ){ disable_fitness_evaluation_ = d; };
	bool disable_fitness_evaluation() const{ return disable_fitness_evaluation_; }

	void n_nucleotide_mut_trials_corrected( core::Real const n ){ n_nucleotide_mut_trials_corrected_ = n; };
	core::Real n_nucleotide_mut_trials_corrected() const{ return n_nucleotide_mut_trials_corrected_; }

	void mutation_rate( core::Real const m ){ mutation_rate_ = m; };
	core::Real mutation_rate() const{ return mutation_rate_; }

	void branch_length( core::Real const m ){ branch_length_ = m; };
	core::Real branch_length() const{ return branch_length_; }

	void total_trials( core::Real const m ){ total_trials_ = m; };
	core::Real total_trials() const{ return total_trials_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	core::Real population_size_; // dlft 10^6
	bool disable_fitness_evaluation_; // dlft false, this is used for benchmarks only
	core::Real n_nucleotide_mut_trials_corrected_;
	core::Real branch_length_;
	core::Real mutation_rate_;
	core::Real total_trials_;
};

} // namespace evolution
} // namespace protocols

#endif
