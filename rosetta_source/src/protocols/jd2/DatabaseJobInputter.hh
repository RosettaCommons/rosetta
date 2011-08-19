// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/DatabaseJobInputter.hh
/// @brief  header file for DatabaseJobInputter class
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_jd2_DatabaseJobInputter_hh
#define INCLUDED_protocols_jd2_DatabaseJobInputter_hh

// Unit Headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/DatabaseJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/features/ProteinSilentReport.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <utility/vector1.hh>

//C++ Headers
#include <string>

namespace protocols {
namespace jd2 {

///@details This is the simplest implementation of JobInputter, which
///reads from -s/-l and Database files.
class DatabaseJobInputter : public protocols::jd2::JobInputter
{
public:

	DatabaseJobInputter();

	virtual ~DatabaseJobInputter();

	static void register_options();

	/// @brief load options from option sytem
	void
	load_options_from_option_system();

	/// @brief Set database file name
	void
	set_database_fname(std::string const & database_fname);

	/// @brief Get database file name
	std::string
	get_database_fname() const;

	/// @brief Get score function
	core::scoring::ScoreFunctionOP
	get_scorefunction();

	/// @brief Set score function
	void
	set_scorefunction(core::scoring::ScoreFunctionOP scorefunction );

	/// @brief Set input tags.  If none are specified, use all tags
	void
	set_tags(utility::vector1< std::string > const & tags);

	/// @brief Get input tags
	void
	get_tags(utility::vector1< std::string > & tags);

	/// @brief this function is responsible for filling the pose reference with
	/// the pose indicated by the job.  The Job object (within its InnerJob)
	/// contains a PoseCOP.  This function needs to either fill the pose
	/// reference from the InnerJob or, on first demand of a pose from that
	/// InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill
	/// the reference.  This implementation uses pose_from_pdb
	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	/// @brief this function determines what jobs exist from -in::file::silent and
	/// -in::file::tags.
	virtual void fill_jobs( Jobs & jobs );

	/// @brief Return the type of input source that the DatabaseJobInputter is currently
	///  using.
	/// @return Always <em>DATABASE</em>.
	virtual JobInputterInputSource::Enum input_source() const;

private:
	core::scoring::ScoreFunctionOP scfxn_;
	protocols::features::ProteinSilentReportOP protein_silent_report_;
	std::string database_fname_;
	utility::vector1< std::string > tags_;

}; // DatabaseJobInputter

} // namespace
} // namespace

#endif //include guard
