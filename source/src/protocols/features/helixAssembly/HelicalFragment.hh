// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

//Numeric
#include <numeric/xyzVector.hh>

//C++ Headers
#include <string>
#include <vector>
#include <map>

namespace protocols {
namespace features {
namespace helixAssembly {

class HelicalFragment : public utility::pointer::ReferenceCount{

public:

	HelicalFragment();
	HelicalFragment(core::Size start, core::Size end);

	~HelicalFragment();

	core::Size seq_start() const;
	core::Size seq_end() const;
	core::Size end() const;
	// Undefined, commenting out to fix PyRosetta build.  std::string get_pdb_source() const;
	core::Size start() const;
	core::Size size() const;
	bool reversed() const;

	void principal_component(numeric::xyzVector<core::Real> principal_component);
	numeric::xyzVector<core::Real> principal_component() const;

	void com(numeric::xyzVector<core::Real> com);
	numeric::xyzVector<core::Real> com() const;

	void p0(numeric::xyzVector<core::Real> p0);
	numeric::xyzVector<core::Real> p0() const;

	void p1(numeric::xyzVector<core::Real> p1);
	numeric::xyzVector<core::Real> p1() const;

	void sasa(core::Real sasa);
	core::Real sasa() const;

	/// @brief output operator
	friend std::ostream & operator <<(std::ostream & os, HelicalFragment const & t);

private:
	core::Size start_;
	core::Size end_;
	std::string pdb_source_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool direction_;
	core::Real sasa_;
	numeric::xyzVector<core::Real> com_;
	numeric::xyzVector<core::Real> principal_component_;
	numeric::xyzVector<core::Real> p0_;
	numeric::xyzVector<core::Real> p1_;
};

}
}
}

#endif
