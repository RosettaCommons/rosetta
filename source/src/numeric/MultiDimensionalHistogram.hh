// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/MultiDimensionalHistogram.hh
///
/// @brief  a class for accumulating a histogram of one or more numeric variables
/// @author Colin A. Smith <colin.smith@ucsf.edu>


#ifndef INCLUDED_numeric_MultiDimensionalHistogram_hh
#define INCLUDED_numeric_MultiDimensionalHistogram_hh


// numeric forward headers
#include <numeric/MultiDimensionalHistogram.fwd.hh>

// numeric headers
#include <numeric/types.hh>

// External library headers
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// C++ headers
#include <ostream>
#include <set>

namespace numeric {


/// @brief a class for accumulating a histogram of one or more numeric variables
class MultiDimensionalHistogram
{

public: // Creation


	/// @brief Constructor
	MultiDimensionalHistogram():
		total_counts_(0)
	{
		num_dimensions(1);
	}

	MultiDimensionalHistogram(
		numeric::Size num_dims
	):
		total_counts_(0)
	{
		num_dimensions(num_dims);
	}


	/// @brief Destructor
	~MultiDimensionalHistogram() {};


	/// @brief Copy constructor
	MultiDimensionalHistogram( MultiDimensionalHistogram const & mdhist )
	{
		*this = mdhist;
	};


private: // Creation


public: // Methods: assignment


	/// @brief operator=
	MultiDimensionalHistogram&
	operator=( MultiDimensionalHistogram const & rhs)
	{
		num_bins_ = rhs.num_bins_;
		start_ = rhs.start_;
		end_ = rhs.end_;
		label_ = rhs.label_;
		dim_labels_ = rhs.dim_labels_;
		total_counts_ = rhs.total_counts_;
		counts_ = rhs.counts_;

		return *this;
	};


public: // Methods: comparison


public: // Methods
	// Further subsections of methods allowed
	numeric::Size
	num_dimensions() const
	{
		return num_bins_.size();
	}

	void
	num_dimensions(
		numeric::Size num
	)
	{
		assert(num > 0);

		num_bins_.resize(num);
		start_.resize(num);
		end_.resize(num);
		end_.resize(num);
		dim_labels_.resize(num);

		num_bins(1);
	}

	utility::vector1<numeric::Size>
	num_bins() const
	{
		utility::vector1<numeric::Size> num_bins1;

		num_bins1.resize(num_bins_.size());
		for ( numeric::Size i = 0; i < num_bins_.size(); ++i ) {
			num_bins1[i+1] = num_bins_[i];
		}

		return num_bins1;
	}

	void
	num_bins(
		numeric::Size num_bins
	)
	{
		assert(num_bins > 0);

		for ( numeric::Size i = 0; i < num_bins_.size(); ++i ) {
			num_bins_[i] = num_bins;
		}

		resize_counts();
	}

	void
	num_bins(
		utility::vector1<numeric::Size> const & num_bins1
	)
	{
		if ( num_bins1.size() != num_dimensions() ) num_dimensions(num_bins1.size());

		for ( numeric::Size i = 0; i < num_bins_.size(); ++i ) {
			assert(num_bins1[i+1] > 0);
			num_bins_[i] = num_bins1[i+1];
		}

		resize_counts();
	}

	void
	num_bins(
		numeric::Size dim,
		numeric::Size num_bins
	)
	{
		assert(dim >= 1 && dim <= num_dimensions());
		assert(num_bins > 0);

		num_bins_[dim-1] = num_bins;

		resize_counts();
	}

	void
	start(
		numeric::Real start
	)
	{
		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			start_[i] = start;
		}
	}

	utility::vector1<numeric::Real>
	start() const
	{
		utility::vector1<numeric::Real> start1;

		start1.resize(start_.size());
		for ( numeric::Size i = 0; i < start_.size(); ++i ) {
			start1[i+1] = start_[i];
		}

		return start1;
	}

	void
	start(
		utility::vector1<numeric::Real> const & start1
	)
	{
		assert(start1.size() == num_dimensions());

		for ( numeric::Size i = 0; i < start_.size(); ++i ) {
			start_[i] = start1[i+1];
		}
	}

	void
	start(
		numeric::Size dim,
		numeric::Real start
	)
	{
		assert(dim >= 1 && dim <= num_dimensions());

		start_[dim-1] = start;
	}

	utility::vector1<numeric::Real>
	end() const
	{
		utility::vector1<numeric::Real> end1;

		end1.resize(end_.size());
		for ( numeric::Size i = 0; i < end_.size(); ++i ) {
			end1[i+1] = end_[i];
		}

		return end1;
	}

	void
	end(
		numeric::Real end
	)
	{
		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			end_[i] = end;
		}
	}

	void
	end(
		utility::vector1<numeric::Real> const & end1
	)
	{
		assert(end1.size() == num_dimensions());

		for ( numeric::Size i = 0; i < end_.size(); ++i ) {
			end_[i] = end1[i+1];
		}
	}

	void
	end(
		numeric::Size dim,
		numeric::Real end
	)
	{
		assert(dim >= 1 && dim <= num_dimensions());

		end_[dim-1] = end;
	}

	void
	range(
		numeric::Size dim,
		numeric::Real start,
		numeric::Real end
	)
	{
		assert(dim >= 1 && dim <= num_dimensions());
		assert(start <= end);

		start_[dim-1] = start;
		end_[dim-1] = end;
	}

	void
	set_dimension(
		numeric::Size dim,
		numeric::Size num_bins,
		numeric::Real start,
		numeric::Real end,
		std::string label = ""
	)
	{
		assert(dim >= 1 && dim <= num_dimensions());
		assert(start <= end);

		num_bins_[dim-1] = num_bins;
		resize_counts();
		start_[dim-1] = start;
		end_[dim-1] = end;
		dim_labels_[dim-1] = label;
	}

	std::string const &
	label() const
	{
		return label_;
	}

	void
	label(
		std::string const & label
	)
	{
		label_ = label;
	}

	utility::vector1<std::string>
	dim_labels() const
	{
		utility::vector1<std::string> dim_labels1;

		dim_labels1.resize(dim_labels_.size());
		for ( numeric::Size i = 0; i < dim_labels_.size(); ++i ) {
			dim_labels1[i+1] = dim_labels_[i];
		}

		return dim_labels1;
	}

	void
	dim_labels(
		utility::vector1<std::string> const & dim_labels1
	)
	{
		assert(dim_labels1.size() == num_dimensions());

		for ( numeric::Size i = 0; i < dim_labels_.size(); ++i ) {
			dim_labels_[i] = dim_labels1[i+1];
		}
	}

	void
	reset_counts()
	{
		for ( numeric::Size i = 0; i < counts_.size(); ++i ) {
			counts_[i] = 0;
		}

		total_counts_ = 0;
	}

	void
	record(
		utility::vector1<numeric::Real> const & values
	)
	{
		assert(values.size() == num_dimensions());

		++counts_[bin_index(values)];
		++total_counts_;
	}

	void
	record(
		numeric::Real value
	)
	{
		assert(num_dimensions() == 1);

		++counts_[bin_index(0, value)];
		++total_counts_;
	}

	utility::vector0<numeric::Size> const &
	counts() const
	{
		return counts_;
	}

	numeric::Size
	total_counts() const
	{
		return total_counts_;
	}

	MultiDimensionalHistogram
	collapse(
		utility::vector1<numeric::Size> dimensions
	) const
	{
		MultiDimensionalHistogram new_mdhist;

		std::set<numeric::Size> dimensions_set;

		for ( numeric::Size i = 1; i <= dimensions.size(); ++i ) {
			if ( dimensions[i] <= num_dimensions() ) {
				dimensions_set.insert(dimensions[i]);
			}
		}

		new_mdhist.num_dimensions(dimensions_set.size());

		utility::vector0<numeric::Size> dim_to_keep;
		utility::vector0<numeric::Size> dim_to_collapse;

		numeric::Size new_counter = 0;
		for ( numeric::Size i = 0; i < num_bins_.size(); ++i ) {

			if ( dimensions_set.count(i+1) ) {

				new_mdhist.num_bins_[new_counter] = num_bins_[i];
				new_mdhist.start_[new_counter] = start_[i];
				new_mdhist.end_[new_counter] = end_[i];
				new_mdhist.dim_labels_[new_counter] = dim_labels_[i];

				++new_counter;

				dim_to_keep.push_back(i);
			} else {
				dim_to_collapse.push_back(i);
			}
		}

		new_mdhist.resize_counts();

		new_mdhist.label_ = label_;
		new_mdhist.total_counts_ = total_counts_;

		if ( dim_to_collapse.size() ) {

			utility::vector0<numeric::Size> idx_current(num_dimensions(), 0);
			utility::vector0<numeric::Size> idx_collapsed(new_mdhist.num_dimensions(), 0);

			// iterate over all counts of the new histogram
			while ( idx_collapsed.back() < new_mdhist.num_bins_.back()+2 ) {

				// copy idx_collapsed into idx_current
				for ( numeric::Size i = 0; i < idx_collapsed.size(); ++i ) {
					idx_current[dim_to_keep[i]] = idx_collapsed[i];
				}

				// reset counters of indicies of idx_current to be collapsed
				for ( numeric::Size i = 0; i < dim_to_collapse.size(); ++i ) {
					idx_current[dim_to_collapse[i]] = 0;
				}

				// add up all the counts from the dimensions to collapse
				while ( idx_current[dim_to_collapse.back()] < num_bins_[dim_to_collapse.back()]+2 ) {

					new_mdhist.counts_[new_mdhist.bin_index(idx_collapsed)] += counts_[bin_index(idx_current)];

					// increment idx_current
					for ( numeric::Size i = 0; i < dim_to_collapse.size(); ++i ) {
						numeric::Size const dim(dim_to_collapse[i]);
						if ( idx_current[dim] == num_bins_[dim]+1 && i < dim_to_collapse.size()-1 ) {
							idx_current[dim] = 0;
						} else {
							++idx_current[dim];
							break;
						}
					}
				}

				// increment idx_collapsed
				for ( numeric::Size i = 0; i < idx_collapsed.size(); ++i ) {
					numeric::Size const dim(i);
					if ( idx_collapsed[dim] == new_mdhist.num_bins_[dim]+1 && i < idx_collapsed.size()-1 ) {
						idx_collapsed[dim] = 0;
					} else {
						++idx_collapsed[dim];
						break;
					}
				}
			}

			numeric::Size count_total = 0;
			for ( numeric::Size i = 0; i < new_mdhist.counts_.size(); ++i ) {
				count_total += new_mdhist.counts_[i];
			}
			runtime_assert(new_mdhist.total_counts_ == count_total);

		} else {

			new_mdhist.counts_ = counts_;
		}

		return new_mdhist;
	}

	numeric::Real
	mean_squared_error(
		utility::vector1<utility::vector1<numeric::Real> > const & expected_1d_frequencies
	)
	{
		runtime_assert(expected_1d_frequencies.size() == num_dimensions());

		for ( numeric::Size i = 1; i <= expected_1d_frequencies.size(); ++i ) {
			runtime_assert(expected_1d_frequencies[i].size() == num_bins()[i]);
		}

		numeric::Real const total_counts(total_counts_);
		numeric::Real mse(0);

		utility::vector0<numeric::Size> idx(num_dimensions(), 1);

		// iterate over all non-overflow counts of the histogram
		while ( idx.back() <= num_bins_.back() ) {

			numeric::Real const frequency(counts_[bin_index(idx)]/total_counts);

			numeric::Real expected_frequency(expected_1d_frequencies[1][idx[0]]);
			for ( numeric::Size i = 1; i < expected_1d_frequencies.size(); ++i ) {
				expected_frequency *= expected_1d_frequencies[i+1][idx[i]];
			}

			numeric::Real const difference(frequency-expected_frequency);
			mse += difference*difference;

			// increment idx
			for ( numeric::Size i = 0; i < idx.size(); ++i ) {
				numeric::Size const dim(i);
				if ( idx[dim] == num_bins_[dim] && i < idx.size()-1 ) {
					idx[dim] = 1;
				} else {
					++idx[dim];
					break;
				}
			}
		}

		numeric::Size num_non_overflow_bins(1);
		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			num_non_overflow_bins *= num_bins_[i];
		}

		mse /= num_non_overflow_bins;

		return mse;
	}

private: // Methods
	// Further subsections of methods allowed

	void
	resize_counts()
	{
		numeric::Size total_bins = 1;

		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			total_bins *= num_bins_[i]+2;
		}

		if ( total_bins != counts_.size() ) counts_.resize(total_bins);
	}

	numeric::Size
	bin_index(
		numeric::Size dim,
		numeric::Real value
	) const
	{
		if ( start_[dim] == end_[dim] ) {
			if ( value < start_[dim] ) return 0;
			if ( value > end_[dim] ) return num_bins_[dim]+1;
			assert(num_bins_[dim] == 1);
			return 1;
		}

		assert(start_[dim] < end_[dim]);
		numeric::Real index = (value-start_[dim])/(end_[dim]-start_[dim])*num_bins_[dim]+1;

		if ( index < 1 ) return 0;
		if ( index > num_bins_[dim]+1 ) return num_bins_[dim]+1;
		if ( index == num_bins_[dim]+1 ) return num_bins_[dim];

		return static_cast<numeric::Size> (index);
	}

	numeric::Size
	bin_index(
		utility::vector1<numeric::Real> const & values
	) const
	{
		assert(values.size() == num_dimensions());

		numeric::Size index = 0;
		numeric::Size interval = 1;

		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			index += bin_index(i, values[i+1])*interval;
			interval *= num_bins_[i]+2;
		}

		assert(index < counts_.size());

		return index;
	}

	numeric::Size
	bin_index(
		utility::vector0<numeric::Size> const & indices
	) const
	{
		assert(indices.size() == num_dimensions());

		numeric::Size index = 0;
		numeric::Size interval = 1;

		for ( numeric::Size i = 0; i < num_dimensions(); ++i ) {
			index += indices[i]*interval;
			interval *= num_bins_[i]+2;
		}

		return index;
	}

private: // Fields

	utility::vector0<numeric::Size> num_bins_;
	utility::vector0<numeric::Real> start_;
	utility::vector0<numeric::Real> end_;

	std::string label_;
	utility::vector0<std::string> dim_labels_;

	numeric::Size total_counts_;
	utility::vector0<numeric::Size> counts_;

}; // MultiDimensionalHistogram

std::ostream & operator << (
	std::ostream & os,
	MultiDimensionalHistogram const & mdhist
);

} // numeric


#endif // INCLUDED_numeric_MultiDimensionalHistogram_HH
