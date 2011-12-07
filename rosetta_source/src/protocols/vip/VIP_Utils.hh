// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/vip/VIP_Utils.hh
/// @brief



#include "core/scoring/packstat/types.hh"
#include "core/scoring/packstat/SimplePDB_Atom.hh"
#include "core/scoring/packstat/SimplePDB.hh"
#include "core/scoring/packstat/AtomRadiusMap.hh"
#include "core/scoring/packstat/compute_sasa.hh"

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/PackStatMover.hh>

#include <devel/init.hh>
#include "core/types.hh"
#include <basic/options/option.hh>
#include "basic/Tracer.hh"

#include "utility/vector1.hh"
#include "utility/file/FileName.hh"
#include "utility/io/izstream.hh"
#include "utility/io/ozstream.hh"
#include "numeric/random/random.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

namespace protocols {
namespace vip {

using core::Real;

using namespace core::scoring::packstat;

inline std::string base_name(const std::string& str);
std::string get_out_tag(std::string fname);
core::Real output_packstat( core::pose::Pose & );

}}
