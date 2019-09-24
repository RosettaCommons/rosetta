// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A Chi record object for a MutableResidueType
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_MutableChiRecord_hh
#define INCLUDED_core_chemical_MutableChiRecord_hh


// Unit headers
#include <core/chemical/MutableChiRecord.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>

//// Project headers
#include <core/types.hh>

//// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <string>

namespace core {
namespace chemical {

/// @brief A class containing bundled info about chis
class MutableChiRecord {
public:
	/// @brief The constructor from four atom
	MutableChiRecord(VD atm1, VD atm2, VD atm3, VD atm4);

	/// @brief The constructor from an atom vector
	MutableChiRecord(utility::vector1< VD > const & atm_vec);

public:

	void set_proton_chi( bool setting = true ) { is_proton_chi_ = setting; }

	void set_proton_chi(
		utility::vector1< Real > const & dihedral_samples,
		utility::vector1< Real > const & extra_samples
	);

	void add_chi_rotamer( Real const mean, Real const sdev );

	void set_chi_rotamers( utility::vector1< std::pair< Real, Real > > const & rots ) {
		chi_rotamers_ = rots;
	}

	void clear_chi_rotamers();

	/// @brief Update the internal VDs based on the provide mapping
	void remap_atom_vds( std::map< VD, VD > const & old_to_new );

public:

	utility::vector1< VD > const & chi_atoms() const { return chi_atoms_; }
	bool is_proton_chi() const { return is_proton_chi_; }
	utility::vector1< Real > const & proton_chi_samples() const { return proton_chi_samples_; }
	utility::vector1< Real > const & proton_chi_extra_samples() const { return proton_chi_extra_samples_; }
	utility::vector1< std::pair< Real, Real > > const & chi_rotamers() const { return chi_rotamers_; }

private:

	/// @brief the four atoms to build each chi angle
	utility::vector1<VD> chi_atoms_;

	/// @brief Is this chi a proton chi?
	bool is_proton_chi_ = false;
	/// @brief For a proton chi, the primary samples to diversify the rotamer library with
	utility::vector1< Real > proton_chi_samples_;
	/// @brief For a proton chi, how to handle extra ex_ levels
	utility::vector1< Real > proton_chi_extra_samples_;

	/// @brief Additional non-Dunbrack rotamer bins
	///
	///    pair<Real,Real>  ==>  mean,sdev
	///    for each chi angle i and rotamer j: chi_rotamers_[i][j]
	///
	utility::vector1< std::pair< Real, Real > > chi_rotamers_;

#ifdef    SERIALIZATION
public:
	MutableChiRecord() = default;

	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core


#endif
