// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MultiStatePacker.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_multistate_design_MultiStatePacker_hh
#define INCLUDED_protocols_multistate_design_MultiStatePacker_hh

#include <protocols/genetic_algorithm/Entity.hh>
#include <protocols/multistate_design/MultiStatePacker.fwd.hh>
#include <protocols/multistate_design/MultiStateFitnessFunction.hh>
#include <protocols/multistate_design/SingleState.fwd.hh>
#include <protocols/multistate_design/PackingState.fwd.hh>

#include <core/chemical/AA.hh>

#include <core/types.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector1.fwd.hh>

#include <iosfwd>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

// the element to be stored by genetic_algorithm::Entity<T> for multistate optimization of residue sequence
class PosType : public genetic_algorithm::EntityElement {
public:
	typedef genetic_algorithm::EntityElement parent;
	typedef genetic_algorithm::EntityElementOP EntityElementOP;

public:
	PosType();
	virtual ~PosType();
	PosType( core::Size index, core::chemical::AA type );
	PosType( std::string word );

	virtual EntityElementOP clone();
	virtual EntityElementOP fresh_instance();
	virtual Size hash() const;
	virtual bool operator <  ( EntityElement const & rhs ) const;
	virtual bool operator == ( EntityElement const & rhs ) const;
	virtual EntityElement & operator =  ( EntityElement const & rhs );
	virtual std::string to_string() const;
	virtual std::string name() const; // Each entity element must have a distinct name

	//core::Size index() const; -- base class defines this now
	core::chemical::AA type() const;

private:
	//core::Size index_; -- base class defines this now
	core::chemical::AA type_;
};

class PosTypeCreator : public genetic_algorithm::EntityElementCreator {
public:
	typedef genetic_algorithm::EntityElementOP EntityElementOP;
public:
	virtual ~PosTypeCreator();
	virtual std::string widget_name() const;
	virtual EntityElementOP new_entity( std::string const & word );
};

////////////////////////////////////////////////////////////////////////////////////////////////////
class MultiStatePacker : public MultiStateFitnessFunction {

public:
	MultiStatePacker()
	: MultiStateFitnessFunction(),
		num_packs_(1)
	{}

	virtual ~MultiStatePacker(){}

	MultiStatePacker( core::Size num_packs )
	: MultiStateFitnessFunction(),
		num_packs_(num_packs)
	{}

	virtual void single_state_design( bool restrict_to_canonical = true );
	virtual core::Real evaluate(
		protocols::genetic_algorithm::Entity & entity,
		core::Size single_state_num
	);

	virtual void set_num_packs( core::Size num ) { num_packs_ = num; }

private:
	core::Size num_packs_;
};

void
limit_rotamer_set(
	utility::vector0< int > & rot_to_pack,
	PackingState const & state,
	genetic_algorithm::EntityElements const & seq // Each EntityElement must be castable to PosType.
);

void
restrict_to_canonical_aas(
	PackingState const & state,
	utility::vector0< int > & rot_to_pack
);

} // namespace multistate_design
} // namespace protocols

#endif
