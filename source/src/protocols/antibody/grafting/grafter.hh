// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/grafter.hh
/// @brief Grafter implementation: take SCS results and create antibody Pose with CDR loops from results
/// @author Sergey Lyskov


#ifndef INCLUDED_protocols_antibody_grafting_grafter_hh
#define INCLUDED_protocols_antibody_grafting_grafter_hh

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/scs_blast.hh>

#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace antibody {
namespace grafting {


/// @brief graft cdr-loops using best scs-results and write results into specified output_prefix
core::pose::PoseOP graft_cdr_loops(AntibodySequence const &A, SCS_ResultSet const &, std::string const & output_prefix, std::string const & suffix, std::string const & database);

// class Grafter
// {
// public:
// 	Grafter() {}

// 	virtual ~Grafter() {}

// 	virtual core::pose::PoseOP graft(AntibodySequence const &A, SCS_ResultSet const &) const = 0;
// };

// class SimpleGrafter : public Grafter
// {
// public:

// 	virtual core::pose::PoseOP graft(AntibodySequence const &A, SCS_ResultSet const &) const override;
// };





} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__


#endif // INCLUDED_protocols_antibody_grafting_grafter_hh
