// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pack/task/rna/RNA_ResidueLevelTask.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pack/task/rna/RNA_ResidueLevelTask.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.pack.task.rna.RNA_ResidueLevelTask" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {
namespace rna {

//Constructor
RNA_ResidueLevelTask::RNA_ResidueLevelTask():
	sample_rna_chi_( false ),
	sample_five_prime_phosphate_( false ),
	sample_three_prime_phosphate_( false ),
	allow_phosphate_virtualization_( false )
{}

//Destructor
RNA_ResidueLevelTask::~RNA_ResidueLevelTask()
{}

} //rna
} //task
} //pack
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::rna::RNA_ResidueLevelTask::save( Archive & arc ) const {
	arc( CEREAL_NVP( sample_rna_chi_ ) ); // _Bool
	arc( CEREAL_NVP( sample_five_prime_phosphate_ ) ); // _Bool
	arc( CEREAL_NVP( sample_three_prime_phosphate_ ) ); // _Bool
	arc( CEREAL_NVP( allow_phosphate_virtualization_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::rna::RNA_ResidueLevelTask::load( Archive & arc ) {
	arc( sample_rna_chi_ ); // _Bool
	arc( sample_five_prime_phosphate_ ); // _Bool
	arc( sample_three_prime_phosphate_ ); // _Bool
	arc( allow_phosphate_virtualization_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::rna::RNA_ResidueLevelTask );
CEREAL_REGISTER_TYPE( core::pack::task::rna::RNA_ResidueLevelTask )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_rna_RNA_ResidueLevelTask )
#endif // SERIALIZATION
