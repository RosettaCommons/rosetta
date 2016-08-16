// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/IdealAbegoGenerator.hh
/// @brief Logic for selection of abego values
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh
#define INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh

// Unit headers
#include <protocols/denovo_design/components/IdealAbegoGenerator.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Logic for selection of abego values
class IdealAbegoGenerator : public utility::pointer::ReferenceCount {
public:
	typedef Segment Motif;
	typedef SegmentOP MotifOP;
	typedef utility::vector1< MotifOP > MotifOPs;
	typedef std::set< core::Size > LengthSet;

public:
	IdealAbegoGenerator();
	virtual ~IdealAbegoGenerator();

	IdealAbegoGeneratorOP
	clone() const;

	/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
	/// @param extend_ss : If true, "extension" of SS elements by adding residues of the same abego type
	///                    is allowed. For example, connecting A-->A, you might get a 4-residue loop
	///                    of type A-ABGA-A. If false, only different abegos are allowed, and you might
	///                    get a 4-residue loop of type A-BGAB-A
	MotifOPs
	generate(
		char const abego1,
		char const abego2,
		LengthSet const & lenset ) const;

private:

};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_IdealAbegoGenerator_hh
