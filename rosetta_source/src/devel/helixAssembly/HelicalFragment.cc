// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelicalFragment.cc
///
/// @brief

/// @author Tim jacobs

//Unit Headers
#include <devel/helixAssembly/HelicalFragment.hh>

HelicalFragment::HelicalFragment(){}

HelicalFragment::HelicalFragment(core::Size start, core::Size end):
start_(start),
end_(end)
{}

HelicalFragment::~HelicalFragment(){}

core::Size HelicalFragment::get_start() const
{
    return start_;
}

core::Size HelicalFragment::get_end() const
{
    return end_;
}

core::Size HelicalFragment::get_size() const
{
  return end_-start_+1;
}

std::string HelicalFragment::get_pdb_source() const
{
    return pdb_source_;
}

bool HelicalFragment::get_direction() const
{
    return direction_;
}

void HelicalFragment::set_start(core::Size start_)
{
    this->start_ = start_;
}

void HelicalFragment::set_end(core::Size end_)
{
    this->end_ = end_;
}

void HelicalFragment::set_pdb_source(std::string pdb_source_)
{
    this->pdb_source_ = pdb_source_;
}

void HelicalFragment::set_direction(bool direction_)
{
    this->direction_ = direction_;
}
