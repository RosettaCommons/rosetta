// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/PDBTrajectoryRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_PDBTrajectoryRecorder_hh
#define INCLUDED_protocols_canonical_sampling_PDBTrajectoryRecorder_hh


// Project forward headers
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>


// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/io/ozstream.hh>


// External library headers


// C++ headers
#include <string>


// Operating system headers


// Forward declarations


namespace protocols {
namespace canonical_sampling {


	/// @brief
class PDBTrajectoryRecorder : public protocols::canonical_sampling::TrajectoryRecorder {
	typedef TrajectoryRecorder Parent;
public:
	/// @brief Constructor
	PDBTrajectoryRecorder();

	/// @brief Destructor
	~PDBTrajectoryRecorder();

	/// @brief Copy constructor
	PDBTrajectoryRecorder( PDBTrajectoryRecorder const & );

private:
	//assignment not allowed
	/// @brief operator=
	PDBTrajectoryRecorder&
	operator=( PDBTrajectoryRecorder const & );

public:
	virtual
	protocols::moves::MoverOP
	clone() const;

	virtual
	protocols::moves::MoverOP
	fresh_instance() const;

	virtual
	std::string
	get_name() const;

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual void reset(
		protocols::moves::MonteCarlo const& mc,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	virtual
	void
	finalize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
	);

private:

	void
	write_model(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	utility::io::ozstream trajectory_stream_;

}; // PDBTrajectoryRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_PDBTrajectoryRecorder_HH
