// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/membrane/benchmark/SampleTiltAnglesCreator.hh
/// @brief Calculates the energy at all possible tilt angles (0->180 degrees)
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_SampleTiltAnglesCreator_hh
#define INCLUDED_protocols_membrane_benchmark_SampleTiltAnglesCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace benchmark {

class SampleTiltAnglesCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();

};

} // benchmark
} // membrane
} // protocols

#endif // INCLUDED_protocols/membrane/benchmark_SampleTiltAngles_fwd_hh



