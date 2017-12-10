// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author  rhiju

#ifndef INCLUDED_protocols_rna_RNA_Fragments_HH
#define INCLUDED_protocols_rna_RNA_Fragments_HH

#include <protocols/rna/denovo/fragments/RNA_Fragments.fwd.hh>
#include <protocols/rna/denovo/fragments/FragmentLibrary.fwd.hh>
#include <protocols/rna/denovo/fragments/TorsionSet.fwd.hh>
#include <protocols/rna/denovo/fragments/RNA_FragmentHomologyExclusion.fwd.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#ifdef WIN32
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#endif

namespace protocols {
namespace rna {
namespace denovo {
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
		RNA_FragmentHomologyExclusionCOP const & homology_exclusion,
		toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
		core::Size const symm_hack_arity ) const;

	virtual FragmentLibraryOP
	get_fragment_library_pointer(
		std::string const & /*RNA_string*/,
		std::string const & /*RNA_secstruct_string*/,
		RNA_FragmentHomologyExclusionCOP const & /*homology_exclusion*/,
		utility::vector1< SYN_ANTI_RESTRICTION > const & /*restriction*/ /*= utility::vector1< SYN_ANTI_RESTRICTION >()*/,
		Size const /*type*/ /* = MATCH_YR */) const {
		/// STUBBED OUT!
		return nullptr;
	}

	virtual void
	insert_fragment(
		core::pose::Pose & ,//pose,
		Size const ,//position,
		TorsionSet const & ,//torsion_set,
		toolbox::AtomLevelDomainMapCOP /*atom_level_domain_map*/ ) const {
		/// STUBBED OUT!
	}

	virtual bool
	is_fullatom();

};

} //fragments
} //denovo
} //rna
} //protocols

#endif
