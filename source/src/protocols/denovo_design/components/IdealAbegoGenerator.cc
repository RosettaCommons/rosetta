// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/IdealAbegoGenerator.cc
/// @brief Logic for selection of abego values
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/IdealAbegoGenerator.hh>

// Protocol headers
#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.IdealAbegoGenerator" );

namespace protocols {
namespace denovo_design {
namespace components {

IdealAbegoGenerator::IdealAbegoGenerator():
	utility::pointer::ReferenceCount()
{

}

IdealAbegoGenerator::~IdealAbegoGenerator()
{}

IdealAbegoGeneratorOP
IdealAbegoGenerator::clone() const
{
	return IdealAbegoGeneratorOP( new IdealAbegoGenerator( *this ) );
}

/// @brief Given desired lengths, compute a set of idealized loop motifs via Nobu/Rie/YuRu rules
/// @param extend_ss : If true, "extension" of SS elements by adding residues of the same abego type
///                    is allowed. For example, connecting A-->A, you might get a 4-residue loop
///                    of type A-ABGA-A. If false, only different abegos are allowed, and you might
///                    get a 4-residue loop of type A-BGAB-A
IdealAbegoGenerator::MotifOPs
IdealAbegoGenerator::generate(
	char const abego1,
	char const abego2,
	LengthSet const & lenset ) const
{
	TR << "Generating motifs for abego1=" << abego1 << " abego2=" << abego2
		<< " lengths=" << lenset << std::endl;
	return MotifOPs();
}

} //protocols
} //denovo_design
} //components

