// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/job_summaries/EnergyJobSummary.cc
/// @brief A JobSummary that simply stores the energy from a job.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/jd3/job_summaries/EnergyJobSummary.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>



#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jd3 {
namespace job_summaries {

using namespace protocols::jd3;

EnergyJobSummary::EnergyJobSummary() :
	JobSummary()
{

}

EnergyJobSummary::EnergyJobSummary( core::Real energy ) :
	JobSummary(),
	energy_( energy )
{

}

EnergyJobSummary::EnergyJobSummary( core::pose::Pose const & pose )
{
	extract_energy( pose );
}

EnergyJobSummary::~EnergyJobSummary() = default;

void
EnergyJobSummary::extract_energy(core::pose::Pose const & pose ){
	if ( pose.energies().energies_updated() ) {
		energy_ = pose.energies().total_energy();
	}
}

core::Real
EnergyJobSummary::energy() const
{
	return energy_;
}

void
EnergyJobSummary::energy( core::Real setting )
{
	energy_ = setting;
}

} //protocols
} //jd3
} //job_summaries

#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::job_summaries::EnergyJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( CEREAL_NVP( energy_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::job_summaries::EnergyJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( energy_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::job_summaries::EnergyJobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::job_summaries::EnergyJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_job_summaries_EnergyJobSummary )
#endif // SERIALIZATION
