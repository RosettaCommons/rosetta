// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/ConformationSelectionToRDMol.hh
/// @brief  This class takes a part of a conformation from Rosetta and converts it into a RDKit RWMol object
/// @author Andy Watkins (watkina6@gene.com)


#ifndef INCLUDED_protocols_drug_design_ConformationSelectionToRDMol_hh
#define INCLUDED_protocols_drug_design_ConformationSelectionToRDMol_hh


#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <rdkit/GraphMol/RWMol.h>

namespace protocols {
namespace drug_design {

class ConformationSelectionToRDMol  {
public:

	/// @brief Convert ResidueType to an RDKit RWMol object
	///
	/// If neutralize is true, attempt to reprotonate into a neutral (non-formally charged) molecule.
	/// (Charges due to heavy atoms, like quaternary amines, will still be present.)
	///
	/// If keep_hydro is true, represent hydrogens as physical atoms, rather than "explicit"
	/// annotations (see below).
	///
	/// Note that due to implementation details, "neutralize = true" will not play well with "keep_hydro = true".
	///
	/// @details RDKit is "aware" of three types of hydrogens.
	///
	/// 1) Physical hydrogens: actual atoms in the atom graph and can have coordinates.
	/// 2) "Explicit" hydrogens: don't have existence in the atom graph or coordinates,
	///    but are instead represented by a field on the heavy atom ("This atom has 3 hydrogens attached to it")
	/// 3) "Implicit" hydrogens: aren't annotated anywhere, but instead are implied by
	///    the difference in the expected valence of the heavy atom and the number of
	///    valences which are currently occupied by bonds and/or "explicit" hydrogens.
	///
	/// Most of RDKit is written around the assumption of non-physical hydrogens, the exception
	/// being things like energy minimization which needs coordinates.
	/// On the other hand, "implicit" hydrogens cause issues in certain cases with kekulization.
	/// Because of this, the default conversion here is to remove the physical hydrogens, and
	/// replace them with "explicit" hydrogen annotations.
	///
	/// Much of RDKit (especially the metric calculations) assumes neutral protonation. (e.g.
	/// what you'd see in an aprotic organic solvent), which is why we neutralize the residues
	/// by default.
	ConformationSelectionToRDMol( core::conformation::Conformation const & conf,
		utility::vector1< core::Size > const & res,
		bool neutralize = true,
		bool keep_hydro = false );

	/// @brief Return an RDKit RWMol object which represents the residue type
	::RDKit::RWMOL_SPTR
	Mol();

private:
	ConformationSelectionToRDMol(); // unimplemented - must provide a residue type

	core::conformation::Conformation const & conf_;
	utility::vector1< core::Size > res_;
	// utility::vector1< MutableResidueType > const & res_;
	bool neutralize_;
	bool keep_hydro_;

};

}
}


#endif
