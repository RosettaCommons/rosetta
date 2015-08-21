// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/ClashCountCalculator.cc
/// @brief  ClashCountCalculator class
/// @author Oliver Lange

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/ClashCountCalculator.hh>

//  Project headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoringManager.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/MetricValue.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

static thread_local basic::Tracer tr( "protocols.metrics.ClashCountCalculator" );

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {
ClashCountCalculator::ClashCountCalculator( core::Real clash_threshold ) :
	clash_threshold_( clash_threshold ),
	vdw_scale_factor_( 0.8 ) // hack from rosetta++
{}

void ClashCountCalculator::lookup( std::string const & key, basic::MetricValueBase *valptr ) const {
	if ( key == "total" ) {
		basic::check_cast( valptr, &total_clashes_, "total_clashes expects to return a Real" );
		( static_cast< basic::MetricValue<core::Size> *>(valptr))->set( total_clashes_ );
	} else if ( key == "bb" ) {
		basic::check_cast( valptr, &bb_clashes_, "bb_clashes expects to return a Size" );
		( static_cast< basic::MetricValue<core::Size> *>(valptr))->set( bb_clashes_ );
	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string ClashCountCalculator::print( std::string const & /*key*/) const {
	/*
	if ( key == "total_sasa" ) {
	return utility::to_string( total_sasa_ );
	} else if ( key == "atom_sasa" ) {
	basic::Error() << "id::AtomID_Map< Real > has no output operator, for metric " << key << std::endl;
	utility_exit();
	} else if ( key == "residue_sasa" ) {
	return utility::to_string( residue_sasa_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	*/
	return "";

}


void ClashCountCalculator::recompute( Pose const& pose ) {
	total_clashes_ = 0;
	bb_clashes_ = 0;
	std::string atom_type_set_name;
	if ( pose.is_fullatom() ) {
		return;
		//  atom_type_set_name = chemical::FA_STANDARD;
	} else {
		atom_type_set_name = chemical::CENTROID;
	}
	core::scoring::AtomVDW const& atom_vdw( core::scoring::ScoringManager::get_instance()->get_AtomVDW( atom_type_set_name ));
	total_clashes_ = 0;
	bb_clashes_ = 0;
	for ( Size ipos = 1; ipos <= pose.total_residue(); ipos++ ) {
		for ( Size jpos = ipos+2; jpos <= pose.total_residue(); jpos++ ) {
			conformation::Residue const & rsd1( pose.residue( ipos ) );
			conformation::Residue const & rsd2( pose.residue( jpos ) );
			if ( !( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) ) {
				for ( Size i = 1, i_end = rsd1.natoms(); i <= i_end; ++i ) {
					Vector const & i_xyz( rsd1.xyz(i) );
					Size const i_type( rsd1.atom_type_index(i) );
					utility::vector1< Real > const & i_atom_vdw( atom_vdw( i_type ) );
					for ( Size j = 1, j_end = rsd2.natoms(); j <= j_end; ++j ) {
						Real const bump_dsq( i_atom_vdw[ rsd2.atom_type_index(j) ] );
						Real const clash( bump_dsq - i_xyz.distance_squared( rsd2.xyz(j) ) );
						if ( clash > 0.0 ) {
							core::Real score( clash*clash/bump_dsq*vdw_scale_factor_ );
							if ( score > clash_threshold_ )  {
								using namespace ObjexxFCL::format;
								std::string const TAG( ( i <= 4 && j <= 4 ) ? "BB BUMP: " : "   BUMP: ");
								tr.Info << TAG << I(4,rsd1.seqpos() ) << I(4,rsd2.seqpos() )
									<<' ' << rsd1.atom_name(i) << ' ' << rsd2.atom_name(j) << ' '
									<< ( clash * clash ) / bump_dsq * vdw_scale_factor_
									<< ' ' << i_xyz.distance_squared( rsd2.xyz(j) ) <<  std::endl;

								++total_clashes_;
								if ( i <= 5 && j <= 5 ) ++bb_clashes_;
							}

						} //clash>0.0
					}//j-atoms
				}//i-atoms
			}// rsd is not bonded
		}// jpos
	}//ipos
}


} // PoseMetricCalculators
} // toolbox
} //protocols
