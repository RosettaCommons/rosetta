// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author  rhiju

#ifndef INCLUDED_protocols_rna_RNA_Fragments_HH
#define INCLUDED_protocols_rna_RNA_Fragments_HH

#include <protocols/farna/fragments/RNA_Fragments.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#endif


namespace protocols {
namespace farna {
namespace fragments {

class RNA_Fragments : public utility::pointer::ReferenceCount {
public:

	//Constructor -- needs vall_torsions_file to get started.
	RNA_Fragments();

	virtual ~RNA_Fragments();

public:

	//Probably the only thing that will actually get called publicly:
	virtual void
	apply_random_fragment(
		core::pose::Pose & pose,
		core::Size const position,
		core::Size const size,
		core::Size const type,
		toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) const;


	virtual bool
	is_fullatom();

};

} //fragments
} //farna
} //protocols

#endif
