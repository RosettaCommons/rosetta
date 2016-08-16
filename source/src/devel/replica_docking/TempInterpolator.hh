// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief This class is used to generate a value according to temperateure level

/// @author Zhe Zhang
#ifndef INCLUDED_devel_replica_docking_TempInterpolator_hh
#define INCLUDED_devel_replica_docking_TempInterpolator_hh

#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <numeric/types.hh>

#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace replica_docking {
class TempInterpolatorBase :  public utility::pointer::ReferenceCount {
public:
	virtual core::Real get_value(core::Size temp_level) = 0;
};

class TempFixValue : public TempInterpolatorBase {
public:
	TempFixValue( core::Real the_value );
	virtual ~TempFixValue();
	virtual core::Real get_value( core::Size ) { return value_; }
private:
	core::Real value_;
};

class TempInterpolator : public TempInterpolatorBase {

public:
	TempInterpolator( core::Size n_levels, core::Real start, core::Real end, std::string curve="exponential" );

	TempInterpolator( utility::tag::TagCOP tag, core::Size n_levels );

	TempInterpolator( TempInterpolator const & temp_interpolator );

	virtual ~TempInterpolator();

	//   void set_end( core::Real end ) {
	//     end_ = end;
	//   }

	//   void set_start( core::Real s ) {
	//     start_ = s;
	//   }

	//   core::Real get_end() {
	//     return end_;
	//   }

	//   core::Real get_start() {
	//     return start_;
	//   }

	void interpolate();

	virtual core::Real get_value( core::Size temp_level );

private:
	core::Real start_;
	core::Real end_;
	std::string curve_;
	core::Size n_levels_;
	//  protocols::canonical_sampling::TemperatureControllerOP tempering_;
	utility::vector1< core::Real > interpolated_nums_;
	bool calculated_;
	//  protocols::canonical_sampling::TemperatureControllerOP tempering_;

};

} // namespace of replica_docking
}

#endif
