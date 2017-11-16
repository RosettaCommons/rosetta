// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file edge reweighting for Interaction Graphs
/// @brief
/// @author Florian Richter, floric@u.washington.edu, june 08

// Unit headers
#include <core/pack/task/IGEdgeReweightContainer.hh>
// Package headers

#include <basic/Tracer.hh>

#include <utility/string_util.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {

static basic::Tracer tr( "core.pack.task.IGEdgeReweightContainer" );

IGEdgeReweighter::~IGEdgeReweighter() {}

IGEdgeReweightContainer::IGEdgeReweightContainer( Size nres )
{
	nres_ = nres;
	edge_reweighters_.clear();
}
IGEdgeReweightContainer::~IGEdgeReweightContainer() {}


Real
IGEdgeReweightContainer::res_res_weight(
	pose::Pose const & pose,
	PackerTask const & task,
	core::Size res1id,
	core::Size res2id
) const {

	debug_assert( res1id <= nres_ );
	debug_assert( res2id <= nres_ );

	Real reweight = 1.0;
	bool firstpass = true;

	for ( utility::vector1< IGEdgeReweighterOP >::const_iterator re_it = edge_reweighters_.begin();
			re_it != edge_reweighters_.end();
			++re_it
			) {

		Real weight_this_upweighter = (*re_it)->get_edge_reweight( pose, task, res1id, res2id );

		if ( firstpass ) {
			reweight = weight_this_upweighter;
			if ( reweight != 1.0 ) firstpass = false;
		} else if ( (reweight != weight_this_upweighter) && (weight_this_upweighter != 1.0)  ) {

			reweight = ( reweight > weight_this_upweighter) ? reweight : weight_this_upweighter ;

			tr.Info << "WARNING WARNING!!! Conflicting IG-reweighting factors specified for residues "+utility::to_string( res1id )+" and "+utility::to_string( res2id )+"are given, will user the larger value ( "+utility::to_string( reweight )+" )." << std::endl;

		}
	}

	return reweight;
}

void
IGEdgeReweightContainer::add_reweighter( IGEdgeReweighterOP reweighter ){

	edge_reweighters_.push_back( reweighter );

}

} //namespace task
} //namespace pack
} //namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::IGEdgeReweighter::save( Archive & ) const {
	// The default weight is set to 1.0 in the constructor and cannot be otherwise set.
	// EXEMPT default_weight_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::IGEdgeReweighter::load( Archive & ) {
	// EXEMPT default_weight_
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::IGEdgeReweighter );
CEREAL_REGISTER_TYPE( core::pack::task::IGEdgeReweighter )


/// @brief Default constructor required by cereal to deserialize this class
core::pack::task::IGEdgeReweightContainer::IGEdgeReweightContainer() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::IGEdgeReweightContainer::save( Archive & arc ) const {
	arc( CEREAL_NVP( edge_reweighters_ ) ); // utility::vector1<IGEdgeReweighterOP>
	arc( CEREAL_NVP( nres_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::IGEdgeReweightContainer::load( Archive & arc ) {
	arc( edge_reweighters_ ); // utility::vector1<IGEdgeReweighterOP>
	arc( nres_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::IGEdgeReweightContainer );
CEREAL_REGISTER_TYPE( core::pack::task::IGEdgeReweightContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_IGEdgeReweightContainer )
#endif // SERIALIZATION

