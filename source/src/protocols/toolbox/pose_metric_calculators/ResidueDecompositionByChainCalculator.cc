// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/ResidueDecompositionByChainCalculator.cc
/// @brief  ResidueDecompositionByChainCalculator class
/// @author Colin A. Smith

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/ResidueDecompositionByChainCalculator.hh>

#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/string_util.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {


ResidueDecompositionByChainCalculator::ResidueDecompositionByChainCalculator() :
	ResidueDecompositionCalculator(),
	use_numbers_(false)
{
}

ResidueDecompositionByChainCalculator::ResidueDecompositionByChainCalculator(
	ResidueDecompositionByChainCalculator const & calculator
):
	ResidueDecompositionCalculator(calculator),
	chain_letters_(calculator.chain_letters()),
	chain_numbers_(calculator.chain_numbers()),
	use_numbers_(calculator.use_numbers())
{
}

core::pose::metrics::PoseMetricCalculatorOP
ResidueDecompositionByChainCalculator::clone() const {
	return core::pose::metrics::PoseMetricCalculatorOP( new ResidueDecompositionByChainCalculator(*this) );
}

void
ResidueDecompositionByChainCalculator::recompute(
	core::pose::Pose const & this_pose
)
{
	if ( use_numbers_ ) {
		// create a map from chain number to set number
		std::map<core::Size, core::Size> chain_map;
		if ( chain_numbers_.size() ) {
			set_names_.assign(chain_numbers_.size(), "");
			for ( core::Size i = 1; i <= chain_numbers_.size(); ++i ) {
				for ( std::set<core::Size>::iterator iter = chain_numbers_[i].begin(); iter != chain_numbers_[i].end(); ++iter ) {
					// there shouldn't be any duplicate chains in the the sets
					runtime_assert(chain_map.find(*iter) == chain_map.end());
					chain_map[*iter] = i;
					if ( set_names_[i] != "" ) set_names_[i] += ",";
					set_names_[i] += utility::to_string(*iter);
				}
			}
			residue_decomposition_.resize(chain_numbers_.size());
		} else {
			// create separate sets for each chain
			set_names_.resize(this_pose.conformation().num_chains());
			for ( core::Size i = 1; i <= this_pose.conformation().num_chains(); ++i ) {
				chain_map[i] = i;
				set_names_[i] = utility::to_string(i);
			}
			residue_decomposition_.resize(this_pose.conformation().num_chains());
		}

		residue_set_numbers_.assign(this_pose.total_residue(), 0);
		for ( core::Size i = 1; i <= this_pose.total_residue(); ++i ) {
			std::map<core::Size, core::Size>::iterator iter(chain_map.find(this_pose.chain(i)));
			if ( iter != chain_map.end() ) {
				residue_decomposition_[iter->second].insert(i);
				residue_set_numbers_[i] = iter->second;
			}
		}

	} else {
		runtime_assert(this_pose.pdb_info() != 0);
		// create a map from chain letter to set number
		std::map<char, core::Size> chain_map;
		if ( chain_letters_.size() ) {
			set_names_.assign(chain_letters_.size(), "");
			for ( core::Size i = 1; i <= chain_letters_.size(); ++i ) {
				for ( std::set<char>::iterator iter = chain_letters_[i].begin(); iter != chain_letters_[i].end(); ++iter ) {
					// there shouldn't be any duplicate chains in the the sets
					runtime_assert(chain_map.find(*iter) == chain_map.end());
					chain_map[*iter] = i;
					if ( set_names_[i] != "" ) set_names_[i] += ",";
					set_names_[i] += utility::to_string(*iter);
				}
			}
			residue_decomposition_.resize(chain_letters_.size());
		} else {
			// create separate sets for each chain
			set_names_.clear();
			for ( core::Size i = 1; i <= this_pose.total_residue(); ++i ) {
				if ( chain_map.find(this_pose.pdb_info()->chain(i)) == chain_map.end() ) {
					core::Size const chain_map_size(chain_map.size());
					chain_map[this_pose.pdb_info()->chain(i)] = chain_map_size+1;
					set_names_.push_back(utility::to_string(this_pose.pdb_info()->chain(i)));
				}
			}
			residue_decomposition_.resize(chain_map.size());
		}

		residue_set_numbers_.assign(this_pose.total_residue(), 0);
		for ( core::Size i = 1; i <= this_pose.total_residue(); ++i ) {
			std::map<char, core::Size>::iterator iter(chain_map.find(this_pose.pdb_info()->chain(i)));
			if ( iter != chain_map.end() ) {
				residue_decomposition_[iter->second].insert(i);
				residue_set_numbers_[i] = iter->second;
			}
		}
	}
}


} // PoseMetricCalculators
} // toolbox
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator >( this ) );
	arc( CEREAL_NVP( chain_letters_ ) ); // utility::vector1<std::set<char> >
	arc( CEREAL_NVP( chain_numbers_ ) ); // utility::vector1<std::set<core::Size> >
	arc( CEREAL_NVP( use_numbers_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculator::load( Archive & arc ) {
	arc( cereal::base_class< protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator >( this ) );
	arc( chain_letters_ ); // utility::vector1<std::set<char> >
	arc( chain_numbers_ ); // utility::vector1<std::set<core::Size> >
	arc( use_numbers_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::ResidueDecompositionByChainCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_ResidueDecompositionByChainCalculator )
#endif // SERIALIZATION
