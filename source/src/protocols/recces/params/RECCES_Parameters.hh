// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/params/RECCES_Parameters.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_parameters_RECCES_Parameters_HH
#define INCLUDED_protocols_recces_parameters_RECCES_Parameters_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/recces/params/RECCES_Parameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace recces {
namespace params {

	class RECCES_Parameters: public utility::pointer::ReferenceCount {

	public:

		//constructor
		RECCES_Parameters( core::pose::Pose const & pose );

		//destructor
		~RECCES_Parameters();

	public:

		void set_bp_res( utility::vector1< core::Real > const & setting ){ bp_res_ = setting; }
		utility::vector1< core::Real > bp_res() const { return bp_res_; }

		void set_dangling_res( utility::vector1< core::Real > const & setting ){ dangling_res_ = setting; }
		utility::vector1< core::Real > dangling_res() const { return dangling_res_; }

	private:

		utility::vector1< core::Real > bp_res_;
		utility::vector1< core::Real > dangling_res_;

	};

} //params
} //recces
} //protocols

#endif
