// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file HelicalFragment.hh
///
/// @brief Small helper class that stores the start and end of a helix secondary structure

/// @author Tim jacobs

#ifndef INCLUDED_protocols_features_helixAssembly_HELICALFRAGMENT_HH
#define INCLUDED_protocols_features_helixAssembly_HELICALFRAGMENT_HH

#include <protocols/features/helixAssembly/HelicalFragment.fwd.hh>

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>

//Utility
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

//Devel
//#include <devel/helixAssembly/NativeResidue.hh>

//External
//#include <boost/serialization/access.hpp>
//#include <boost/serialization/map.hpp>
//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/string.hpp>

//C++ Headers
#include <string>
#include <vector>
#include <map>

namespace protocols {
namespace features {
namespace helixAssembly {

class HelicalFragment : public utility::pointer::ReferenceCount{

public:

	HelicalFragment(core::Size start, core::Size end);

	~HelicalFragment();

	core::Size seq_start() const;
	core::Size seq_end() const;
	core::Size end() const;
	std::string get_pdb_source() const;
	core::Size start() const;
	core::Size size() const;
	bool reversed() const;
	
	void sasa(core::Real sasa);
	core::Real sasa() const;

private:
	core::Size start_;
	core::Size end_;
	std::string pdb_source_;
	bool direction_;
	core::Real sasa_;
};

}
}
}

#endif
