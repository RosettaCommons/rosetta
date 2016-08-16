// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file GeneticAlgorithm.hh
/// @brief genetic algorithm protocols forward declaration
/// @author ashworth, based on template "pseudo"code by Colin Smith

#ifndef INCLUDED_protocols_genetic_algorithm_GeneticAlgorithm_fwd_hh
#define INCLUDED_protocols_genetic_algorithm_GeneticAlgorithm_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace genetic_algorithm {

class GeneticAlgorithmBase;

typedef  utility::pointer::shared_ptr< GeneticAlgorithmBase >  GeneticAlgorithmBaseOP;
typedef  utility::pointer::shared_ptr< GeneticAlgorithmBase const >  GeneticAlgorithmBaseCOP;

class GeneticAlgorithm;

typedef  utility::pointer::shared_ptr< GeneticAlgorithm >  GeneticAlgorithmOP;
typedef  utility::pointer::shared_ptr< GeneticAlgorithm const >  GeneticAlgorithmCOP;

} // namespace genetic_algorithm
} // namespace protocols

#endif
