// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerComb.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeRotamerComb_HH
#define INCLUDED_protocols_rotamer_sampler_copy_dofs_ResidueAlternativeRotamerComb_HH

#include <protocols/rotamer_sampler/RotamerSizedComb.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerComb.fwd.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamer.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#ifdef WIN32
	#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamer.hh>
#endif


/*
using namespace core::conformation;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/

namespace protocols {
namespace rotamer_sampler {
namespace copy_dofs {

	class ResidueAlternativeRotamerComb: public RotamerSizedComb {

	public:

		//constructor
    ResidueAlternativeRotamerComb();

		//destructor
		~ResidueAlternativeRotamerComb();

	public:

		core::conformation::Residue const &
		get_residue_at_origin( Size const seqpos );

		core::conformation::Residue const &
		get_residue_at_origin_with_matching_type( Size const seqpos, core::conformation::Residue const & rsd_in );

		void
		add_residue_alternative_rotamer( ResidueAlternativeRotamerOP const & rotamer );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ResidueAlternativeRotamerComb"; }

		/// @brief Type of class (see enum in RotamerTypes.hh)
		virtual RotamerType type() const { return RESIDUE_ALTERNATIVE_COMB; }


		bool
		has_resnum( Size const seqpos );

		Size
		find_resnum( Size const seqpos );

		Size
		id_for_resnum( Size const seqpos );

		void
		fast_forward_to_next_residue_pair( Size const i, Size const j);

		void
		fast_forward_to_next_residue( Size const i );

	private:

    using RotamerSizedComb::add_rotamer; // make it private.

		std::map< Size, ResidueAlternativeRotamerOP > residue_alternative_rotamer_map_;

		// following is redundant with rotamer_list_ in parent class, but is useful since it retains ResidueAlternativeRotamer type.
		utility::vector1< ResidueAlternativeRotamerOP > residue_alternative_rotamer_list_;

	};

} //copy_dofs
} //rotamer_sampler
} //protocols

#endif
