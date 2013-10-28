// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/toolbox/PoseMetricCalculator/FragQualCalculator.hh
/// @brief header file for FragQualCalculator class.
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_FragQualCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_FragQualCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

//// C++ headers
// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

class FragQualCalculator : public core::pose::metrics::StructureDependentCalculator {
public:


	typedef core::pose::metrics::StructureDependentCalculator Super;
	typedef std::string String;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pose::metrics::PoseMetricCalculatorOP PoseMetricCalculatorOP;
  typedef basic::MetricValueBase MetricValueBase;
	typedef core::fragment::FragSet FragSet;
	typedef core::fragment::FragSetOP FragSetOP;

	typedef utility::tag::TagCOP TagCOP;
	typedef protocols::filters::Filters_map Filters_map;
	typedef basic::datacache::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	/// @brief default constructor
	FragQualCalculator();

	/// @brief default constructor
	FragQualCalculator( FragSetOP const & frag, Real const rmsd=1.0, Real const ratio=30.0 );

	/// @brief copy constructor
	FragQualCalculator( FragQualCalculator const & rval );

	/// @brief destructor
	virtual ~FragQualCalculator();


public:// virtual constructor


	/// @brief make clone
  PoseMetricCalculatorOP clone() const { return new FragQualCalculator( *this ); }


public:// mutator


	/// @brief set fragments
	void set_fragset( FragSetOP const & frags );

	/// @brief rmsd cutoff of good fragments
	void rmsd_cutoff( Real const & val );

	/// @brief
	void ratio_cutoff( Real const & val );

	/// @brief
	void set_region( Size const val1, Size const val2 );

	/// @brief
	void begin( Size const begin ) { begin_ = begin; }

	/// @brief
	void end( Size const end ) { end_ = end; }


public:


	void parse_my_tag( TagCOP const tag,
										 basic::datacache::DataMap & data,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & pose );


protected:


  virtual void lookup( String const & key, MetricValueBase * valptr ) const;
  virtual std::string print( String const & key ) const;
  virtual void recompute( Pose const & this_pose );


private:

	/// @brief
	Real rmsd_cutoff_goodfrag_;

	/// @brief
	Real ratio_cutoff_goodfrag_;

	/// @brief
	Real total_goodfrags_;
	/// @brief
	Real coverage_;

	/// @brief
	utility::vector1< Size > goodfrags_;

	/// @brief
	FragSetOP frag_;

	/// @brief
	Size begin_;

	/// @brief
	Size end_;

	/// @brief
	bool verbose_;


}; //FragQualCalculator


} // ns PoseMetricCalculators
} // ns toolbox
} // ns protocols

#endif
