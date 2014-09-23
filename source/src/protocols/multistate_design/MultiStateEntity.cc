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

// Unit headers
#include <protocols/multistate_design/MultiStateEntity.hh>

#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/multistate_design/SingleStateEntityData.hh>

#include <core/types.hh>
#include <basic/MetricValue.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

// C++ Headers
#include <map>

#include <utility/vector1.hh>



namespace protocols {
namespace multistate_design {

MultiStateEntity::MultiStateEntity() : protocols::genetic_algorithm::Entity() {}

MultiStateEntity::MultiStateEntity( MultiStateEntity const & entity ) :
	protocols::genetic_algorithm::Entity(entity),
	single_state_entity_data_(entity.single_state_entity_data())
{}

MultiStateEntity::~MultiStateEntity() {}

MultiStateEntity::EntityOP
MultiStateEntity::clone() const { return MultiStateEntity::EntityOP( new MultiStateEntity(*this) ); }

void MultiStateEntity::show( std::ostream & os ) const
{
	using namespace ObjexxFCL::format;

	os << "MultiStateEntity with traits:";
	genetic_algorithm::EntityElements const & seq( this->traits() );
	for ( genetic_algorithm::EntityElements::const_iterator
			it( seq.begin() ), end( seq.end() );
			it != end; ++it ) {
		os << " " << (*it)->to_string();
	}
	os << " and fitness " << F(6,3,this->fitness());
	for (core::Size i = 1; i <= single_state_entity_data_.size(); ++i) {
		os << '\n' << " SingleState " << i << " with fitness: " << single_state_entity_data_[i].fitness();
	}
}


void
MultiStateEntity::write_checkpoint(
	std::ostream & os
) const
{
	protocols::genetic_algorithm::Entity::write_checkpoint(os);

	os << " states " << single_state_entity_data_.size();
	for (core::Size i = 1; i <= single_state_entity_data_.size(); ++i) {
		os << "\n ";
		single_state_entity_data_[i].write_checkpoint(os);
	}
}

bool
MultiStateEntity::read_checkpoint(
	std::istream & is
)
{
	if (!protocols::genetic_algorithm::Entity::read_checkpoint(is)) return false;

	std::string word;
	if (!(is >> word)) return false;
	if (word != "states") return false;

	core::Size num_states;
	if (!(is >> num_states)) return false;
	single_state_entity_data_.resize(num_states);

	for (core::Size i = 1; i <= single_state_entity_data_.size(); ++i) {
		if (!single_state_entity_data_[i].read_checkpoint(is)) return false;
	}

	return true;
}

utility::vector1< SingleStateEntityData > const &
MultiStateEntity::single_state_entity_data() const { return single_state_entity_data_; }

utility::vector1< SingleStateEntityData > &
MultiStateEntity::single_state_entity_data() { return single_state_entity_data_; }


} // namespace multistate_design
} // namespace protocols

