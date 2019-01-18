// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/job_summaries/EnergyJobSummary.hh
/// @brief A JobSummary that simply stores the energy from a job.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_job_summaries_EnergyJobSummary_HH
#define INCLUDED_protocols_jd3_job_summaries_EnergyJobSummary_HH

// Unit headers
#include <protocols/jd3/job_summaries/EnergyJobSummary.fwd.hh>
#include <protocols/jd3/JobSummary.hh>
#include <core/pose/Pose.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jd3 {
namespace job_summaries {

///@brief A JobSummary that simply stores the energy from a job.
class EnergyJobSummary : public JobSummary
{
public:
	EnergyJobSummary();

	EnergyJobSummary( core::Real energy );

	EnergyJobSummary( core::pose::Pose const & pose );

	virtual ~EnergyJobSummary();

public:

	///@brief Extract the energy fromt the Energies object of the pose
	void
	extract_energy( core::pose::Pose const & pose );

	///@brief Set an energy of the job.
	void
	energy( core::Real setting );

	///@brief Get the energy stored in this JobSummary.
	core::Real
	energy() const;


private:
	core::Real energy_ = 0.0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //protocols
} //jd3
} //job_summaries

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_job_summaries_EnergyJobSummary )
#endif // SERIALIZATION

#endif //protocols_jd3_job_summaries_EnergyJobSummary_HH

