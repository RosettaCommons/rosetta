// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/DOFHistogramRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/simple_moves/DOFHistogramRecorder.hh>

// Other project headers or inline function headers
#include <core/chemical/ResidueType.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>

// External library headers

// C++ headers
#include <iomanip>
#include <sstream>
#include <set>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/MultiDimensionalHistogram.hh>


// Operating system headers

// Forward declarations


namespace protocols {
namespace simple_moves {

DOFHistogramRecorder::DOFHistogramRecorder() :
	num_bins_(10)
{}

DOFHistogramRecorder::~DOFHistogramRecorder()
{}

DOFHistogramRecorder::DOFHistogramRecorder( DOFHistogramRecorder const & /* other */ )
{
	// copy constructor not allowed
	runtime_assert(false);
}

DOFHistogramRecorder&
DOFHistogramRecorder::operator=( DOFHistogramRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return * this;
}

void
DOFHistogramRecorder::insert_dofs_by_residue(
	core::pose::Pose const & pose,
	utility::vector1<core::id::DOF_ID_Range> dof_ranges
)
{
	utility::vector1<std::set<core::id::DOF_ID_Range> > residue_dof_ranges;
	residue_dof_ranges.resize(pose.total_residue());

	for ( core::Size i = 1; i <= dof_ranges.size(); ++i ) {

		core::id::AtomID parent_id( pose.atom_tree().atom(dof_ranges[i].dof_id().atom_id()).parent()->id() );
		residue_dof_ranges[parent_id.rsd()].insert(dof_ranges[i]);
	}

	for ( core::Size i = 1; i <= residue_dof_ranges.size(); ++i ) {

		if ( residue_dof_ranges[i].size() ) {

			dofs_.resize(dofs_.size()+1);
			histograms_.resize(histograms_.size()+1);
			dof_values_.resize(dof_values_.size()+1);

			utility::vector1<core::id::DOF_ID> & dofs(dofs_[dofs_.size()]);
			numeric::MultiDimensionalHistogram & histogram(histograms_[histograms_.size()]);
			utility::vector1<core::Real> & dof_values(dof_values_[dof_values_.size()]);

			std::ostringstream res_label;
			res_label << pose.residue_type(i).name3() << ' ' << i;
			histogram.label(res_label.str());

			utility::vector1<core::id::DOF_ID_Range> r_dof_ranges(residue_dof_ranges[i].begin(), residue_dof_ranges[i].end());

			dofs.resize(r_dof_ranges.size());
			histogram.num_dimensions(r_dof_ranges.size());
			dof_values.resize(r_dof_ranges.size());

			for ( core::Size j = 1; j <= r_dof_ranges.size(); ++j ) {

				dofs[j] = r_dof_ranges[j].dof_id();

				std::ostringstream label;
				label << r_dof_ranges[j].dof_id();
				histogram.set_dimension(j, num_bins_, r_dof_ranges[j].min(), r_dof_ranges[j].max(), label.str());
			}
		}
	}
}

void
DOFHistogramRecorder::update_after_boltzmann(
	core::pose::Pose const & pose
)
{
	for ( core::Size i = 1; i <= dofs_.size(); ++i ) {
		for ( core::Size j = 1; j <= dofs_[i].size(); ++j ) {
			dof_values_[i][j] = pose.dof(dofs_[i][j]);
		}
		histograms_[i].record(dof_values_[i]);
	}
}

void
DOFHistogramRecorder::write_mse_summary(
	std::ostream & os
) const
{
	for ( core::Size i = 1; i <= dofs_.size(); ++i ) {

		numeric::MultiDimensionalHistogram const & histogram(histograms_[i]);

		utility::vector1<utility::vector1<core::Real> > expected_frequencies_all(dofs_[i].size());
		utility::vector1<utility::vector1<core::Real> > expected_frequencies_subset(1);

		for ( core::Size j = 1; j <= dofs_[i].size(); ++j ) {
			expected_frequencies_all[j] = uniform_dof_distribution(
				dofs_[i][j].type(), histogram.num_bins()[j], histogram.start()[j], histogram.end()[j]
			);

			expected_frequencies_subset[1] = expected_frequencies_all[j];
			utility::vector1<core::Size> const dims(1, j);

			core::Real mse = histogram.collapse(dims).mean_squared_error(expected_frequencies_subset);

			os << histogram.dim_labels()[j] << ": " << mse << std::endl;
		}
	}
}

std::ostream & operator << (
	std::ostream & os,
	DOFHistogramRecorder const & dof_recorder)
{
	for ( core::Size i = 1; i <= dof_recorder.histograms().size(); ++i ) {

		os << dof_recorder.histograms()[i];

		// print out collapsed histograms for debugging purposes
		//   utility::vector1<core::Size> dim_to_keep(1);
		//   for (core::Size j = 1; j <= dof_recorder.histograms()[i].num_dimensions(); ++j) {
		//
		//    dim_to_keep[1] = j;
		//    numeric::MultiDimensionalHistogram mdhist(dof_recorder.histograms()[i].collapse(dim_to_keep));
		//    os << mdhist;
		//
		//    core::id::DOF_Type dof_type(dof_recorder.dofs()[i][j].type());
		//    utility::vector1<core::Real> expected(
		//     uniform_dof_distribution(dof_type, mdhist.num_bins()[1], mdhist.start()[1], mdhist.end()[1])
		//    );
		//
		//    for (core::Size i = 1; i < expected.size(); ++i) os << expected[i] << " ";
		//    os << std::endl;
		//   }
	}

	return os;
}

utility::vector1<core::Real>
uniform_dof_distribution(
	core::id::DOF_Type dof_type,
	core::Size num_bins,
	core::Real min,
	core::Real max
)
{
	utility::vector1<core::Real> frequencies;

	if ( dof_type == core::id::PHI ) {

		frequencies.resize(num_bins, 1.0/num_bins);

	} else if ( dof_type == core::id::THETA ) {

		frequencies.resize(num_bins);
		core::Real frequencies_total(0);

		core::Real const delta((max-min)/num_bins);
		for ( core::Size i = 1; i <= num_bins; ++i ) {

			core::Real const interval_start(min+(i-1)*delta);
			core::Real const interval_end(min+i*delta);

			// integral of sin function
			frequencies[i] = cos(interval_start) - cos(interval_end);
			frequencies_total += frequencies[i];
		}

		for ( core::Size i = 1; i <= num_bins; ++i ) {
			frequencies[i] /= frequencies_total;
		}

	} else {
		// nothing else supported yet
		runtime_assert(false);
	}

	return frequencies;
}


} // namespace moves
} // namespace protocols
