// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/toolbox/PoseMetricCalculator/FragQualCalculator.hh
/// @brief header file for FragQualCalculator class.
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_FragQualCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_FragQualCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/ResidueIndexDescription.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

//// C++ headers
// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
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
	typedef basic::datacache::DataMap DataMap;


public:// constructor/destructor


	/// @brief default constructor
	FragQualCalculator();

	/// @brief default constructor
	FragQualCalculator( FragSetOP frag, Real const rmsd=1.0, Real const ratio=30.0 );

	/// @brief copy constructor
	FragQualCalculator( FragQualCalculator const & rval );

	/// @brief destructor
	~FragQualCalculator() override;


public:// virtual constructor


	/// @brief make clone
	PoseMetricCalculatorOP clone() const override { return utility::pointer::make_shared< FragQualCalculator >( *this ); }


public:// mutator


	/// @brief set fragments
	void set_fragset( FragSetOP const & frags );

	/// @brief rmsd cutoff of good fragments
	void rmsd_cutoff( Real const & val );

	/// @brief
	void ratio_cutoff( Real const & val );

	/// @brief
	void set_region( core::Size const val1, core::Size const val2 );

	/// @brief
	void begin( core::Size const begin );

	/// @brief
	void end( core::Size const end );


public:

	// Note that this isn't used in a general Factory system
	// but is instead called directly from FragQualFilter
	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data
	);


protected:


	void lookup( String const & key, MetricValueBase * valptr ) const override;
	std::string print( String const & key ) const override;
	void recompute( Pose const & this_pose ) override;


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
	utility::vector1< core::Size > goodfrags_;

	/// @brief
	FragSetOP frag_;

	/// @brief
	core::pose::ResidueIndexDescriptionCOP begin_;

	/// @brief
	core::pose::ResidueIndexDescriptionCOP end_;

	/// @brief
	bool verbose_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //FragQualCalculator


} // ns PoseMetricCalculators
} // ns protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_FragQualCalculator )
#endif // SERIALIZATION


#endif
