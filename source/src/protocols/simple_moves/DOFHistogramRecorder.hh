// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_moves/DOFHistogramRecorder.hh
///
/// @brief
/// @author


#ifndef INCLUDED_protocols_simple_moves_DOFHistogramRecorder_hh
#define INCLUDED_protocols_simple_moves_DOFHistogramRecorder_hh


// Project forward headers
#include <protocols/simple_moves/DOFHistogramRecorder.fwd.hh>


// Project headers
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


// External library headers


// C++ headers

#include <core/id/DOF_ID.fwd.hh>
#include <core/id/types.hh>
#include <utility/vector0.hh>
#include <numeric/MultiDimensionalHistogram.fwd.hh>
#include <numeric/types.hh>


// Operating system headers


// Forward declarations


namespace protocols {
namespace simple_moves {


/// @brief
class DOFHistogramRecorder
{
	// Friends


public: // Types


private: // Types


public: // Constants


private: // Constants


public: // Creation


	/// @brief Constructor
	DOFHistogramRecorder();


	/// @brief Destructor
	~DOFHistogramRecorder();


	/// @brief Copy constructor - not allowed
	DOFHistogramRecorder( DOFHistogramRecorder const & ) = delete;


private: // Creation


public: // Methods: assignment


	/// @brief operator= - not allowed
	DOFHistogramRecorder&
	operator=( DOFHistogramRecorder const & ) = delete;


public: // Methods: comparison


public: // Methods

	core::Size
	num_bins() const
	{
		return num_bins_;
	}

	void
	num_bins(
		core::Size num_bins
	)
	{
		num_bins_ = num_bins;
	}

	void
	insert_dofs_by_residue(
		core::pose::Pose const & pose,
		utility::vector1<core::id::DOF_ID_Range> dof_ranges
	);

	utility::vector1<utility::vector1<core::id::DOF_ID> > const &
	dofs() const
	{
		return dofs_;
	}

	utility::vector1<numeric::MultiDimensionalHistogram> const &
	histograms() const
	{
		return histograms_;
	}

	void
	update_after_boltzmann(
		core::pose::Pose const & pose
	);

	void
	write_mse_summary(
		std::ostream & os
	) const;


private:


public: // Properties


private: // Fields

	utility::vector1<utility::vector1<core::id::DOF_ID> > dofs_;
	utility::vector1<numeric::MultiDimensionalHistogram> histograms_;
	utility::vector1<utility::vector1<core::Real> > dof_values_;
	numeric::Size num_bins_;

}; // DOFHistogramRecorder


std::ostream & operator << (
	std::ostream & os,
	DOFHistogramRecorder const & dof_recorder
);


utility::vector1<core::Real>
uniform_dof_distribution(
	core::id::DOF_Type dof_type,
	core::Size num_bins,
	core::Real min,
	core::Real max
);


} // namespace moves
} // namespace protocols


#endif // INCLUDED_protocols_simple_moves_DOFHistogramRecorder_HH
