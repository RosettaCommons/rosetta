// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/sample_stream/SampleStream.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sample_stream_SampleStream_HH
#define INCLUDED_protocols_sample_stream_SampleStream_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/sample_stream/SampleStream.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace sample_stream {

	class SampleStream: public utility::pointer::ReferenceCount {

	public:

	//constructor
	SampleStream();

	//destructor
	~SampleStream();
	public:

		virtual bool get_next_sample( core::pose::Pose & ) = 0;

		virtual void reset() = 0;

	private:

	};

} //sample_stream
} //protocols

#endif
