// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/upstream/ProteinSCSampler.hh
/// @brief  Class declarations for base class and simple derived class for SC sampling
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_ProteinSCSampler_hh
#define INCLUDED_protocols_match_upstream_ProteinSCSampler_hh

// Unit headers
#include <protocols/match/upstream/ProteinSCSampler.fwd.hh>

// Package headers
#include <protocols/match/upstream/ScaffoldBuildPoint.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <iterator>


namespace protocols {
namespace match {
namespace upstream {

class ProteinSCSampler : public utility::pointer::ReferenceCount {
public:
	typedef core::pack::dunbrack::DunbrackRotamerSampleData DunbrackRotamerSampleData;
	typedef utility::vector1< DunbrackRotamerSampleData >      DunbrackRotamerSampleDataVector;

public:
	virtual
	DunbrackRotamerSampleDataVector
	samples(
		ScaffoldBuildPoint const & bb_conf,
		core::chemical::ResidueType const & restype
	) const = 0;

};

/// @brief Basic class for sidechain sampling that pulls data from the
/// Dunbrack rotamer library.  The samples that are returned are the basic
/// rotamers and do not include any expansions by the "ex" flags.
class DunbrackSCSampler : public ProteinSCSampler {

public:
	virtual
	DunbrackRotamerSampleDataVector
	samples(
		ScaffoldBuildPoint const & bb_conf,
		core::chemical::ResidueType const & restype
	) const;

};

}
}
}

#endif
