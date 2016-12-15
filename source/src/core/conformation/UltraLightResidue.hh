// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/UltraLightResidue.hh
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)
/// @brief an extremely light residue class used for fast direct cartesian space transforms on ligands.

#ifndef INCLUDED_core_conformation_UltraLightResidue_hh
#define INCLUDED_core_conformation_UltraLightResidue_hh


#include <core/conformation/UltraLightResidue.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/types.hh>
#include <numeric/xyzMatrix.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>


namespace core {
namespace conformation {

class UltraLightResidue : public utility::pointer::ReferenceCount {

public:
	UltraLightResidue(ResidueCOP residue);

	UltraLightResidue(UltraLightResidue const & src);

	~UltraLightResidue() override = default;

	/// @brief update conformation with current coords. Slow.
	void update_conformation(Conformation & conformation) const;

	//@brief apply some transformation matrix and translation perturbation
	void transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector );

	void transform(numeric::xyzMatrix<core::Real> const & rotation_matrix, core::Vector const & translation_vector, core::Vector group_center);


	//@brief align the residue on to some other residue
	void align_to_residue(UltraLightResidue const & other_residue);

	//@brief move residue to some center point
	void slide(core::Vector const & translation_vector);

	//@brief index based access to the xyz coordinates
	PointPosition & operator[](core::Size index)
	{
		return coords_[index];
	}

	//@brief index based access to the xyz coordinates
	PointPosition const & operator[](core::Size index) const
	{
		return coords_[index];
	}

	/// @brief return a constant pointer to the base residue
	ResidueCOP residue() const
	{
		return residue_;
	}

	/// @brief return number of atoms in ultralight residue
	core::Size natoms() const
	{
		return coords_.size();
	}

	/// @brief get centerpoint of residue
	PointPosition center() const
	{
		return center_;
	}

	///@brief get const ref to residue coords
	void set_coords(utility::vector1<PointPosition> const & input_coords)
	{
		coords_ = input_coords;
	}

	///@brief get const ref to residue coords

	utility::vector1<PointPosition> const & coords_vector() const
	{
		return coords_;
	}

	///@brief get const ref to residue heavy atom coords
	utility::vector1<PointPosition> const & heavy_coords_vector() const
	{
		return heavy_coords_;
	}

private:
	//you really don't want to build up one of these from scratch
	UltraLightResidue();
private:

	//datatypes used for easy interfacing with conformation::batch_set_xyz and batch_get_xyz
	utility::vector1<id::AtomID> atom_ids_;
	utility::vector1<PointPosition > coords_;
	utility::vector1<PointPosition > heavy_coords_;

	//read only access to the base residue.
	ResidueCOP residue_;
	PointPosition center_;

};

}
}


#endif /* ULTRALIGHTRESIDUE_HH_ */
