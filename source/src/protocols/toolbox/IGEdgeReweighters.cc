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
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>


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

IGInterfaceEdgeUpweighter::IGInterfaceEdgeUpweighter( core::Real weight_factor ) :
	weight_factor_ ( weight_factor )
{
	set_skip_loop_chains( "" );
	set_sec_str("");
}

IGInterfaceEdgeUpweighter::IGInterfaceEdgeUpweighter( core::Real weight_factor, std::string const & skip_loop_in_chain, std::string const & sec_str_in ) :
	weight_factor_ ( weight_factor )
{
	set_skip_loop_chains( skip_loop_in_chain );
	set_sec_str( sec_str_in );
}

void
IGInterfaceEdgeUpweighter::set_skip_loop_chains( std::string const & chain_string )
{
	chains_to_ignore_loops_.clear();
	utility::vector1<std::string> chains_vec = utility::string_split( chain_string, ',' );
	for ( std::string const & iichain_string : chains_vec ) {
		if ( iichain_string.size() == 0 ) continue;
		if ( iichain_string.size() != 1 ) {
			utility_exit_with_message("Chain identifier should not be more than one character!!!");
		}
		chains_to_ignore_loops_.insert( iichain_string[0] );
	}
}

core::Real
IGInterfaceEdgeUpweighter::get_edge_reweight(
	pose::Pose const & pose,
	pack::task::PackerTask const &,
	Size res1,
	Size res2
) const
{

	if ( pose.chain(res1) != pose.chain(res2) ) {
		if ( chains_to_ignore_loops_.size() != 0 ) {
			if ( !pose.pdb_info() ) {
				utility_exit_with_message("InterfaceUpweightor received a pose without a valid PDBInfo--chains cannot be selected.");
			}
			if ( sec_str_.size() != pose.size() ) {
				utility_exit_with_message("Secondary struture string doesn't match with the length of the pose!!!! So the pose was changed?");
			}

			if ( ( chains_to_ignore_loops_.count( pose.pdb_info()->chain( res1 ) ) && (sec_str_[res1-1] == 'L') ) ||
					( chains_to_ignore_loops_.count( pose.pdb_info()->chain( res2 ) ) && (sec_str_[res2-1] == 'L') ) ) {
				return default_weight_;
			} else {
				return weight_factor_;
			}
		} else {
			return weight_factor_;
		}
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

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::IGInterfaceEdgeUpweighter::IGInterfaceEdgeUpweighter() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::IGInterfaceEdgeUpweighter::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( CEREAL_NVP( weight_factor_ ) ); // core::Real
	arc( CEREAL_NVP( chains_to_ignore_loops_ ) );
	arc( CEREAL_NVP( sec_str_ ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::IGInterfaceEdgeUpweighter::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::task::IGEdgeReweighter >( this ) );
	arc( weight_factor_ ); // core::Real
	arc( chains_to_ignore_loops_ );
	arc( sec_str_ );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::IGInterfaceEdgeUpweighter );
CEREAL_REGISTER_TYPE( protocols::toolbox::IGInterfaceEdgeUpweighter )


CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_IGEdgeReweighters )
#endif // SERIALIZATION
