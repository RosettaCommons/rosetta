// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file edge reweighting
/// @brief collection of routines to assign different weights to IG edges
/// @author Florian Richter, floric@u.washington.edu, june 08

// Unit headers
#include <protocols/toolbox/IGEdgeReweighters.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {

using namespace core;

core::Real
IGLigandDesignEdgeUpweighter::get_edge_reweight(
	pose::Pose const & pose,
	pack::task::PackerTask const & task,
	Size res1,
	Size res2
) const
{

	if ( ( pose.residue( res1 ).is_ligand() && task.design_residue( res2 ) )
			||( pose.residue( res2 ).is_ligand() && task.design_residue( res1 ) ) ) {
		return weight_factor_;
	} else return default_weight_;

}

/*
template <class T>
core::Real
ResidueGroupIGEdgeUpweighter::get_edge_reweight(
pose::Pose const & pose,
pack::task::PackerTask const & task,
Size res1,
Size res2
) const
{

if( (group1_.find(res1) != group1_.end()) && (group2_.find(res2) != group2_.end() )
||(group2_.find(res1) != group2_.end()) && (group1_.find(res2) != group1_.end() ) ){
return weight_factor_;
}
else return default_weight_;

}
*/

} //namespace toolbox
} //namespace protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::IGLigandDesignEdgeUpweighter::IGLigandDesignEdgeUpweighter() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::IGLigandDesignEdgeUpweighter::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( CEREAL_NVP( weight_factor_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::IGLigandDesignEdgeUpweighter::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( weight_factor_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::IGLigandDesignEdgeUpweighter );
CEREAL_REGISTER_TYPE( protocols::toolbox::IGLigandDesignEdgeUpweighter )


/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::ResidueGroupIGEdgeUpweighter::ResidueGroupIGEdgeUpweighter() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::ResidueGroupIGEdgeUpweighter::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( CEREAL_NVP( weight_factor_ ) ); // core::Real
	arc( CEREAL_NVP( group1_ ) ); // std::set<core::Size>
	arc( CEREAL_NVP( group2_ ) ); // std::set<core::Size>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::ResidueGroupIGEdgeUpweighter::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( weight_factor_ ); // core::Real
	arc( group1_ ); // std::set<core::Size>
	arc( group2_ ); // std::set<core::Size>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::ResidueGroupIGEdgeUpweighter );
CEREAL_REGISTER_TYPE( protocols::toolbox::ResidueGroupIGEdgeUpweighter )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_IGEdgeReweighters )
#endif // SERIALIZATION
