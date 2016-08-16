// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/JD2ResourceManagerJobInputterCreator.hh
/// @brief
/// @author Mattthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_jd2_JD2ResourceManagerJobInputterCreator_hh
#define INCLUDED_protocols_jd2_JD2ResourceManagerJobInputterCreator_hh

#include <protocols/jd2/JobInputterCreator.hh>

namespace protocols {
namespace jd2 {

/// @brief the %JD2ResourceManagerJobInputterCreator is responsible for instantiating
/// the JD2ResourceManagerJobInputter for the JobDistributorFactory
class JD2ResourceManagerJobInputterCreator : public protocols::jd2::JobInputterCreator {
public:

	/// @brief Create the JD2ResourceManagerJobInputter
	virtual JobInputterOP create_JobInputter() const;

	/// @brief Return the name of the JobInputter this class instantiates:
	/// "JD2ResourceManagerJobInputter"
	virtual std::string keyname() const;

};

} //jd2
} //protocols

#endif //INCLUDED_protocols_jd2_JD2ResourceManagerJobInputterCreator_hh

