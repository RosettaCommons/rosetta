// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/enzdes/RemoveLigandFilter.hh
//p
/// @brief Check if the ligand's pocket is stable by removing the ligand, relaxing the structure and calculating rms to the starting structure.
/// @author Javier Castellanos (javiercv@ue.edu)

#ifndef INCLUDED_protocols_enzdes_RemoveLigandFilter_hh
#define INCLUDED_protocols_enzdes_RemoveLigandFilter_hh

// unit headers
#include <protocols/enzdes/RemoveLigandFilter.fwd.hh>

// package headers

// Project Headers
#include <core/io/silent/SilentEnergy.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace enzdes {


class RemoveLigandFilter: public protocols::filters::Filter {

public:

typedef protocols::filters::Filter Filter;
typedef protocols::filters::FilterOP FilterOP;
typedef protocols::filters::Filters_map Filters_map;
typedef protocols::moves::MoverOP MoverOP;
typedef core::pose::Pose Pose;

public:
    //default constructor
    RemoveLigandFilter( );
    // value constructor
    RemoveLigandFilter( core::Real threshold );

    RemoveLigandFilter( RemoveLigandFilter const & rval );

    bool apply( Pose const & pose ) const;

    // Undefined, commenting out to fix PyRosetta build  void set_min_mover( MoverOP min_mover );

    // Undefined, commenting out to fix PyRosetta build  core::Real compute( Pose const & pose ) const;

    core::Real report_sm( Pose const & pose ) const;

    void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

    FilterOP clone() const { return new RemoveLigandFilter( *this ); }
    FilterOP fresh_instance() const { return new RemoveLigandFilter; }
private:
    core::Real threshold_;
    MoverOP mover_;
    FilterOP filter_;
};
} // enzdes
} // protocols


#endif /*INCLUDED_protocols_enzdes_RemoveLigandFilter_HH*/
