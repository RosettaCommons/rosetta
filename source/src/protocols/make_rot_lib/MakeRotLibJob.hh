// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/make_rot_lib/MakeRotLibJob.fwd.hh
/// @brief  Header file for MakeRotLibJob class
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

#ifndef INCLUDED_protocols_make_rot_lib_MakeRotLibJob_hh
#define INCLUDED_protocols_make_rot_lib_MakeRotLibJob_hh

// unit headers
#include <protocols/make_rot_lib/MakeRotLibJob.fwd.hh>
#include <protocols/make_rot_lib/MakeRotLibOptionsData.hh>

#include <utility/vector1.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.fwd.hh>

namespace protocols {
namespace make_rot_lib {

class MakeRotLibJob : public protocols::jd2::Job
{
public:

	MakeRotLibJob( jd2::InnerJobOP inner_job, core::Size nstruct_index,
		core::Real omg, utility::vector1< core::Real> const & bbs, utility::vector1< core::Size > const & bb_ids, core::Real eps,
		MakeRotLibOptionsDataOP mrlod, bool semirotameric );

	~MakeRotLibJob() override = default;

	/// @brief acessors
	core::Real get_omg() const { return omg_; }
	core::Real get_bb(core::Size i) const { return bbs_[ i ]; }
	utility::vector1< core::Real > get_bbs() const { return bbs_; }
	utility::vector1< core::Size > get_bb_ids() const { return bb_ids_; }
	core::Real get_eps() const { return eps_; }
	MakeRotLibOptionsDataOP get_options_data() const { return mrlod_; }
	bool get_semirotameric() const { return semirotameric_; }

private:

	/// @brief backbone torsion angle values that are unique to this job
	core::Real omg_;
	utility::vector1< core::Real > bbs_;
	utility::vector1< core::Size > bb_ids_;
	core::Real eps_;

	/// @brief access to the options data info that are not unique to this job
	MakeRotLibOptionsDataOP mrlod_;

	bool semirotameric_;
};

}//make_rot_lib
}//protocols

#endif //INCLUDED_protocols_make_rot_lib_MakeRotLibJob_hh
