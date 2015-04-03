// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStateEntity.hh
/// @brief
/// @author Colin A. Smith

#ifndef INCLUDED_protocols_multistate_design_MultiStateEntity_hh
#define INCLUDED_protocols_multistate_design_MultiStateEntity_hh

#include <protocols/genetic_algorithm/Entity.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>

#include <protocols/multistate_design/SingleStateEntityData.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

class MultiStateEntity : public protocols::genetic_algorithm::Entity {

public:
	typedef genetic_algorithm::Entity::OP EntityOP;
	typedef genetic_algorithm::EntityElementOP EntityElementOP;

	MultiStateEntity();
	MultiStateEntity( MultiStateEntity const & entity );
	virtual ~MultiStateEntity();

	virtual EntityOP clone() const;

	virtual void show( std::ostream & os ) const;

	virtual void write_checkpoint( std::ostream & os ) const;
	virtual bool read_checkpoint( std::istream & is );

	utility::vector1< SingleStateEntityData > const & single_state_entity_data() const;
	utility::vector1< SingleStateEntityData > & single_state_entity_data();

private:
	utility::vector1< SingleStateEntityData > single_state_entity_data_;
};


} // namespace multistate_design
} // namespace protocols

#endif
