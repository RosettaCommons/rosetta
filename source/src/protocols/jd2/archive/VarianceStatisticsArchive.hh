// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @author Oliver Lange
/// Archive class to collect structures such that variances of scores can be computed to determine normalized weights


#ifndef INCLUDED_protocols_jd2_archive_VarianceStatisticsArchive_hh
#define INCLUDED_protocols_jd2_archive_VarianceStatisticsArchive_hh

// Unit Headers
//#include <protocols/abinitio/IterativeAbrelax.fwd.hh>

// Package Headers
#include <protocols/jd2/archive/EvaluatedArchive.hh>
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/jd2/archive/VarianceStatisticsArchive.fwd.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// Third-party Headers

//// C++ headers

#include <string>


namespace protocols {
namespace jd2 {
namespace archive {

class VarianceStatisticsArchive : public EvaluatedArchive {
	typedef EvaluatedArchive Parent;
	typedef utility::vector1< core::io::silent::SilentStructOP > SilentStructVector;

public:
	VarianceStatisticsArchive( std::string name );

	virtual bool add_evaluated_structure(
		core::io::silent::SilentStructOP,
		core::io::silent::SilentStructOP alternative_decoy,
		Batch const&
	);

	virtual void generate_batch() {};
	/// @brief overloaded to make input decoys appear the same as decoys coming from batches
	virtual void init_from_decoy_set( core::io::silent::SilentFileData const& ) {};

	void set_insertion_prob( core::Real setting ) {
		insertion_prob_ = setting;
	}

protected:

private:
	core::Real insertion_prob_;
};


}
}
}

#endif
