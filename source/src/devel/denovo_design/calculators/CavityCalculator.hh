// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/calculators/CavityCalculator.hh
/// @brief header file for CavityCalculator class.
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_devel_denovo_design_calculators_CavityCalculator_hh
#define INCLUDED_devel_denovo_design_calculators_CavityCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Utility headers
#include <numeric/xyzVector.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace denovo_design {
namespace calculators {

class CavityCalculator : public core::pose::metrics::StructureDependentCalculator {
public:// constructor/destructor

	/// @brief default constructor
	CavityCalculator();

	/// @brief destructor
	virtual ~CavityCalculator();

public:// virtual constructor
	/// @brief make clone
  core::pose::metrics::PoseMetricCalculatorOP
	clone() const { return new CavityCalculator( *this ); }

public:// mutators

protected:
  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
  virtual std::string print( std::string const & key ) const;
  virtual void recompute( core::pose::Pose const & this_pose );

private: // private member functions

private:// member variables
	core::Real total_volume_;
  utility::vector1< core::scoring::packstat::CavityBallCluster > clusters_;
}; //CavityCalculator


} // ns calculators
} // ns denovo_design
} // ns devel

#endif
