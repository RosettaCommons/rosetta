// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MetricRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_moves_MetricRecorder_hh
#define INCLUDED_protocols_moves_MetricRecorder_hh


// Project forward headers
#include <protocols/moves/MetricRecorder.fwd.hh>


// Project headers
#include <protocols/moves/ThermodynamicObserver.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/io/ozstream.hh>


// External library headers


// C++ headers
#include <string>

#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>



// Operating system headers


// Forward declarations


namespace protocols {
namespace moves {


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
		core::pose::Pose const & pose
	);

	void
	update_after_boltzmann(
		core::pose::Pose const & pose
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



public: // Properties




private: // Fields

	std::string file_name_;
	core::Size stride_;
	core::Size step_count_;
	utility::io::ozstream recorder_stream_;
	utility::vector1<std::pair<core::id::TorsionID, std::string> > torsion_ids_;

}; // MetricRecorder


} // namespace moves
} // namespace protocols


#endif // INCLUDED_protocols_moves_MetricRecorder_HH
