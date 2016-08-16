// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/mmt_msd/MMTPackingJob.hh
/// @brief  declaration for class MMTPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTPackingJob_HH
#define INCLUDED_devel_mmt_msd_MMTPackingJob_HH

// Unit headers
#include <devel/mmt_msd/MMTPackingJob.fwd.hh>

// core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

namespace devel {
namespace mmt_msd {

class MMTPackingJob : public utility::pointer::ReferenceCount
{
public:
	typedef std::list< std::pair< core::Size, std::string > > required_npds;
	typedef std::pair< core::Size, core::Real > npd_property_and_value;
	typedef std::list< npd_property_and_value > npd_properties;

public:
	MMTPackingJob();
	virtual ~MMTPackingJob();

	void set_pose( core::pose::Pose const & pose );
	void set_sfxn( core::scoring::ScoreFunction const & sfxn );
	void set_packer_task( core::pack::task::PackerTask const & task );
	void set_npd_properties( required_npds const & npd_to_do_list );

	core::pose::Pose const &             get_pose() const;
	core::scoring::ScoreFunction const & get_sfxn() const;
	core::pack::task::PackerTaskCOP      get_task() const;

	bool has_pose() const;
	bool has_sfxn() const;
	bool has_task() const;

	virtual void setup() = 0;
	virtual void optimize() = 0;
	virtual void update_pose( core::pose::Pose & pose ) = 0;

	void go();

	core::Real running_time() const;
	virtual core::Real final_energy() const = 0;
	bool optimization_complete() const;

	core::Size n_npd_properties() const;
	npd_properties::const_iterator npd_properties_begin() const;
	npd_properties::const_iterator npd_properties_end() const;

protected:
	core::pose::Pose & pose();
	core::scoring::ScoreFunction & sfxn();
	core::pack::task::PackerTaskOP task();

	/// @brief Compute the set of npd properties for the Pose that
	/// results when the optimized rotamers are assigned to the pose_.
	/// Performs work only if there are any npd properties for this job
	void compute_npd_properties();

	/// @brief deallocates the pose, sfxn, and task.
	virtual void clean_up();

private:
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP sfxn_;
	core::pack::task::PackerTaskOP task_;
	required_npds npd_to_do_list_;
	npd_properties computed_npd_properties_;

	core::Real running_time_; // in seconds, negative if it has not yet been run.

};

}
}

#endif
