// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/TrajectoryRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_moves_TrajectoryRecorder_hh
#define INCLUDED_protocols_moves_TrajectoryRecorder_hh


// Project forward headers
#include <protocols/moves/TrajectoryRecorder.fwd.hh>


// Project headers
#include <protocols/moves/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/io/ozstream.hh>


// External library headers


// C++ headers
#include <string>

#include <utility/vector1.hh>



// Operating system headers


// Forward declarations


namespace protocols {
namespace moves {


	/// @brief
class TrajectoryRecorder : public ThermodynamicObserver {
public:
	/// @brief Constructor
	TrajectoryRecorder();

	/// @brief Destructor
	~TrajectoryRecorder();

	/// @brief Copy constructor
	TrajectoryRecorder( TrajectoryRecorder const & );

	/// @brief operator=
	TrajectoryRecorder&
	operator=( TrajectoryRecorder const & );

	virtual
	MoverOP
	clone() const;

	virtual
	MoverOP
	fresh_instance() const;

	virtual
	std::string
	get_name() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	std::string const &
	file_name() const
	{
		return file_name_;
	}

	void
	file_name(
		std::string const & file_name
	)
	{
		file_name_ = file_name;
	}

	core::Size
	stride() const
	{
		return stride_;
	}

	void
	stride(
		core::Size stride
	)
	{
		stride_ = stride;
	}

	void
	reset(
		protocols::moves::MonteCarlo const & mc
	);

	void
	update_after_boltzmann(
		core::pose::Pose const & pose
	);

	void
	update_after_boltzmann(
		protocols::moves::MonteCarlo const & mc
	);

	virtual
	void
	apply(
		core::pose::Pose & pose
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	void
	observe_after_metropolis(
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::moves::MetropolisHastingsMover const & metropolis_hastings_mover
	);

private:

	void
	write_model(
		core::pose::Pose const & pose
	);

	std::string file_name_;
	core::Size stride_;
	core::Size model_count_;
	core::Size step_count_;
	utility::io::ozstream trajectory_stream_;

}; // TrajectoryRecorder


} // namespace moves
} // namespace protocols


#endif // INCLUDED_protocols_moves_TrajectoryRecorder_HH
