// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/calculators/CavityCalculator.hh
/// @brief header file for CavityCalculator class.
/// Roughly, fragment quality is number of fragments which are close to a pose in rmsd
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) -- Moved from devel to protocols.


#ifndef INCLUDED_protocols_denovo_design_calculators_CavityCalculator_hh
#define INCLUDED_protocols_denovo_design_calculators_CavityCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/packstat/compute_sasa.hh>

// Utility headers
#include <numeric/xyzVector.fwd.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

//// C++ headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace denovo_design {
namespace calculators {

class CavityCalculator : public core::pose::metrics::StructureDependentCalculator {
public:// constructor/destructor

	/// @brief default constructor
	CavityCalculator();

	/// @brief destructor
	~CavityCalculator() override;

public:// virtual constructor
	/// @brief make clone
	core::pose::metrics::PoseMetricCalculatorOP
	clone() const override { return utility::pointer::make_shared< CavityCalculator >( *this ); }

public:// mutators

protected:
	void lookup( std::string const & key, basic::MetricValueBase * valptr ) const override;
	std::string print( std::string const & key ) const override;
	void recompute( core::pose::Pose const & this_pose ) override;

private: // private member functions

private:// member variables
	core::Real total_volume_;
	utility::vector1< core::scoring::packstat::CavityBallCluster > clusters_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; //CavityCalculator


} // ns calculators
} // ns denovo_design
} // ns protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_denovo_design_calculators_CavityCalculator )
#endif // SERIALIZATION


#endif
