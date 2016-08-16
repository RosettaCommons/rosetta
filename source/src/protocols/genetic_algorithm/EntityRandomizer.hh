// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file EntityRandomizer.hh
/// @brief controls the alteration of the traits that define Entity
/// @author ashworth

#ifndef INCLUDED_protocols_genetic_algorithm_EntityRandomizer_hh
#define INCLUDED_protocols_genetic_algorithm_EntityRandomizer_hh

#include <protocols/genetic_algorithm/Entity.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


#include <algorithm> // std::copy

namespace protocols {
namespace genetic_algorithm {

////////////////////////////////////////////////////////////////////////////////////////////////////
class EntityRandomizer : public utility::pointer::ReferenceCount {

public:
	typedef utility::pointer::shared_ptr< EntityRandomizer > OP;
	typedef utility::pointer::shared_ptr< EntityRandomizer const > COP;
	typedef Entity::OP EntityOP;
	typedef Entity::COP EntityCOP;

	EntityRandomizer();
	virtual ~EntityRandomizer();
	virtual EntityOP random_entity();
	virtual void mutate( Entity & entity ) = 0;
	virtual void crossover( Entity & entity1, Entity & entity2 );

	virtual core::Size entity_length() const { return entity_length_; }
	virtual void set_mutation_rate( core::Real rate ) { mutation_rate_ = rate; }
	virtual core::Real mutation_rate() const { return mutation_rate_; }
	virtual core::Size library_size() const = 0;
	virtual EntityCOP entity_template() const;
	virtual void set_entity_template(EntityCOP entity);

protected:
	virtual void set_entity_length( core::Size length );

private:
	core::Size entity_length_;
	core::Real mutation_rate_;
	EntityCOP entity_template_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
class DiscreteRandomizer : public EntityRandomizer {
public:

	virtual ~DiscreteRandomizer();
	virtual void add_choice( EntityElementOP const & choice );
	virtual void set_choices( EntityElements const & choices );
	virtual void mutate( Entity & entity );
	virtual core::Size library_size() const;
	virtual EntityElements const & choices() const;
private:
	EntityElements choices_;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief different set of choices at each position in Entity's traits
class PositionSpecificRandomizer : public EntityRandomizer {
public:
	typedef utility::pointer::shared_ptr< PositionSpecificRandomizer > OP;
	typedef utility::pointer::shared_ptr< PositionSpecificRandomizer const > COP;

	virtual ~PositionSpecificRandomizer();
	virtual void append_choices( EntityElements const & choices );
	virtual void mutate( Entity & entity );
	virtual core::Size library_size() const;
	virtual utility::vector1< EntityElements > const & choices() const;
private:
	utility::vector1< EntityElements > choices_;
};


} // namespace genetic_algorithm
} // namespace protocols

#endif
