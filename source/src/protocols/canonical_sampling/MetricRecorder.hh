// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/MetricRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_canonical_sampling_MetricRecorder_hh
#define INCLUDED_protocols_canonical_sampling_MetricRecorder_hh


// Project forward headers
#include <protocols/canonical_sampling/MetricRecorder.fwd.hh>


// Project headers
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/io/ozstream.hh>


// External library headers


// C++ headers
#include <ctime>
#include <string>

#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>



// Operating system headers


// Forward declarations


namespace protocols {
namespace canonical_sampling {


	/// @brief
class MetricRecorder : public ThermodynamicObserver
{
	// Friends


public: // Types


private: // Types




public: // Constants


private: // Constants




public: // Creation


	/// @brief Constructor
	MetricRecorder();


	/// @brief Destructor
	~MetricRecorder();


	/// @brief Copy constructor
	MetricRecorder( MetricRecorder const & );


private: // Creation




public: // Methods: assignment


	/// @brief operator=
	MetricRecorder&
	operator=( MetricRecorder const & );


public: // Methods: comparison



public: // Methods

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

	virtual
	bool
	reinitialize_for_each_job() const { return true; };

	std::string const &
	file_name() const;

	void
	file_name(
		std::string const & file_name
	);

	core::Size
	stride() const;

	void
	stride(
		core::Size stride
	);

	bool
	cumulate_jobs() const;

	void
	cumulate_jobs(
		bool cumulate_jobs
	);

	bool
	cumulate_replicas() const;

	void
	cumulate_replicas(
		bool cumulate_replicas
	);

	bool
	prepend_output_name() const;

	void
	prepend_output_name(
		bool prepend_output_name
	);

	void
	add_torsion(
		core::id::TorsionID const & torsion_id,
		std::string name = ""
	);

	void
	add_torsion(
		core::pose::Pose const & pose,
		std::string const & rsd,
		std::string type,
		core::Size torsion,
		std::string name = ""
	);

	void
	reset(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
	);

	void
	update_after_boltzmann(
		core::pose::Pose const & pose,
		protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover = 0
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
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	virtual
	void
	observe_after_metropolis(
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
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
		core::pose::Pose const & pose
	);



public: // Properties




private: // Fields

	std::string file_name_;
	core::Size stride_;
	bool cumulate_jobs_;
	bool cumulate_replicas_;
	bool prepend_output_name_;
	core::Size step_count_;
	utility::io::ozstream recorder_stream_;
	utility::vector1<std::pair<core::id::TorsionID, std::string> > torsion_ids_;
	time_t last_flush_;

}; // MetricRecorder


} // namespace canonical_sampling
} // namespace protocols


#endif // INCLUDED_protocols_canonical_sampling_MetricRecorder_HH
