// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Entity.hh
/// @brief the unit employed/optimized by GeneticAlgorithm
/// @author ashworth

#ifndef INCLUDED_protocols_genetic_algorithm_Entity_fwd_hh
#define INCLUDED_protocols_genetic_algorithm_Entity_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace protocols {
namespace genetic_algorithm {

class EntityElement;
typedef utility::pointer::shared_ptr< EntityElement > EntityElementOP;
typedef utility::pointer::shared_ptr< EntityElement const > EntityElementCOP;

typedef utility::vector1< EntityElementOP > EntityElements;


class EntityElementCreator;
typedef utility::pointer::shared_ptr< EntityElementCreator > EntityElementCreatorOP;
typedef utility::pointer::shared_ptr< EntityElementCreator const > EntityElementCreatorCOP;

class EntityElementFactory;

class Entity;
typedef utility::pointer::shared_ptr< Entity > EntityOP;
typedef utility::pointer::shared_ptr< Entity const > EntityCOP;

struct Vec1Hash;

} // namespace genetic_algorithm
} // namespace protocols

#endif
