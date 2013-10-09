// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/FindEndpointsOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_devel_splice_FindEndpointsOperation_hh
#define INCLUDED_devel_splice_FindEndpointsOperation_hh

// Unit Headers
#include <devel/splice/FindEndpointsOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
namespace devel {
namespace splice {

///@details Find the endpoints on a beta barrel
    class FindEndpointsOperation : public protocols::toolbox::task_operations::RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	FindEndpointsOperation();

	virtual ~FindEndpointsOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagPtr );


    core::Size Cterm_offset() const{ return Cterm_offset_; }
    void Cterm_offset( core::Size const c ){ Cterm_offset_ = c; }
    core::Size Nterm_offset() const{ return Nterm_offset_; }
    void Nterm_offset( core::Size const n ){ Nterm_offset_ = n;}
    bool even() const{ return even_; }
    void even( bool const b ){ even_ = b; }
    bool odd() const{ return odd_; }
    void odd( bool const b ){ odd_ = b; }

    core::Real distance_cutoff() const{ return distance_cutoff_; }
    void distance_cutoff( core::Real const d ) { distance_cutoff_ = d;}

    core::Size neighbors() const{ return neighbors_; }
    void neighbors( core::Size const s ) { neighbors_ = s; }

    bool point_inside() const{ return point_inside_; }
    void point_inside( bool const b ){ point_inside_ = b; }

private:
    core::Size Cterm_offset_, Nterm_offset_; //dflt 0,0; whether to count residues N and C-terminally
    bool even_, odd_; //dflt true, true; report on even and odd blades
    core::Real distance_cutoff_; //dflt 18.0A; what distance for counting nterminal neighbours
    core::Size neighbors_; //dflt 6; how many neighbors to require

    bool point_inside_; //dflt true; select residues that point into the barrel
};

} //namespace splice
} //namespace devel

#endif // INCLUDED_protocols_toolbox_TaskOperations_FindEndpointsOperation_HH
