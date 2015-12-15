// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_BallInfo.hh
/// @brief  Orientation dependent variant of the LK Solvation using
/// @author Phil Bradley

#ifndef INCLUDED_core_scoring_methods_LK_BallInfo_HH
#define INCLUDED_core_scoring_methods_LK_BallInfo_HH


// Unit headers
// #include <core/scoring/methods/LK_BallInfo.fwd.hh> ??

// // Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace methods {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @details  Stores the internal coordinates of an ideal water position
class WaterBuilder {
public:
	WaterBuilder(
		Vector const & water,
		conformation::Residue const & rsd,
		Size const atom1,
		Size const atom2,
		Size const atom3
	);

	Vector
	build( conformation::Residue const & rsd ) const;

	// fpd get the derivative of water movement w.r.t. base atom movement
	// fpd SLOW!
	void
	derivatives( conformation::Residue const & rsd,
		numeric::xyzMatrix <Real> &dw_da1, numeric::xyzMatrix <Real> &dw_da2, numeric::xyzMatrix <Real> &dw_da3 ) const;

	Size atom1() const { return atom1_; }
	Size atom2() const { return atom2_; }
	Size atom3() const { return atom3_; }

private:
	Size atom1_;
	Size atom2_;
	Size atom3_;
	Vector xyz_local_;
};

typedef utility::vector1< WaterBuilder > WaterBuilders;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @details  Holds the locations of ideal waters attached to the atoms of a Residue
class LKB_ResidueInfo; // fwd
typedef utility::pointer::shared_ptr< LKB_ResidueInfo > LKB_ResidueInfoOP;

class LKB_ResidueInfo : public basic::datacache::CacheableData {
public:
	virtual ~LKB_ResidueInfo();
	typedef utility::vector1< Vector > Vectors;

public:

	LKB_ResidueInfo( conformation::Residue const & rsd );

	LKB_ResidueInfo( LKB_ResidueInfo const & src );

	LKB_ResidueInfo();

	void
	initialize( chemical::ResidueType const & rsd );

	virtual
	basic::datacache::CacheableDataOP
	clone() const;

	void
	build_waters( conformation::Residue const & rsd );

	// fpd const access to the water builders (to identify stub atoms)
	WaterBuilders const &
	get_water_builder( conformation::Residue const & rsd , Size heavyatom ) const;

	utility::vector1< Vectors > const &
	waters() const { return waters_; }

	//fpd deriv of water position w.r.t. atom 1
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > const &
	atom1_derivs() const { return dwater_datom1_; }

	//fpd deriv of water position w.r.t. atom 2
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > const &
	atom2_derivs() const { return dwater_datom2_; }

	//fpd deriv of water position w.r.t. atom 3
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > const &
	atom3_derivs() const { return dwater_datom3_; }

	bool
	has_waters() const { return has_waters_; }

	utility::vector1< utility::vector1< Real > > const &
	atom_weights() const { return atom_weights_; }

	void
	remove_irrelevant_waters(
		Size const atom,
		chemical::ResidueType const & rsd_type,
		utility::vector1< Vector > & waters
	) const;


	/// danger
	static
	void
	reset_arrays_danger_expert_only();

	bool
	matches_residue_type( chemical::ResidueType const & rsd_type ) const;

	chemical::ResidueType const &
	residue_type() const;

	/////////////////////////////////////////////////////////////////////////////
	// STATIC data
private:
	/////////////////////////////////////////////////////////////////////////////

	typedef std::map< chemical::ResidueType const *, utility::vector1< WaterBuilders > > WaterBuilderMap;
	static WaterBuilderMap water_builder_map_;

	typedef std::map< chemical::ResidueType const *, utility::vector1< utility::vector1< Real > > > AtomWeightsMap;
	static AtomWeightsMap atom_weights_map_;

	void
	initialize_residue_type( chemical::ResidueType const & rsd_type ) const;

	void
	setup_atom_weights(
		chemical::ResidueType const & rsd_type,
		utility::vector1< WaterBuilders > const & rsd_water_builders, // for sanity
		utility::vector1< utility::vector1< Real > > & atom_wts
	) const;

private:
	chemical::ResidueTypeCOP rsd_type_;
	utility::vector1< Vectors > waters_;
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > dwater_datom1_;
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > dwater_datom2_;
	utility::vector1< utility::vector1< numeric::xyzMatrix< Real > > > dwater_datom3_;
	utility::vector1< utility::vector1< Real > > atom_weights_;
	bool has_waters_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::pointer::shared_ptr< LKB_ResidueInfo > LKB_ResidueInfoOP;
typedef utility::pointer::shared_ptr< const LKB_ResidueInfo > LKB_ResidueInfoCOP;


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class LKB_ResiduesInfo : public basic::datacache::CacheableData {
public:
	LKB_ResiduesInfo() {}

	LKB_ResiduesInfo( LKB_ResiduesInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const;

	Size
	size() const { return residues_info_.size(); }

	LKB_ResidueInfo const &
	operator[]( Size const index ) const { return *( residues_info_[ index ] ); }

	LKB_ResidueInfo &
	operator[]( Size const index ) { return *( residues_info_[ index ] ); }

	void
	append( LKB_ResidueInfoOP rsd_info ) { residues_info_.push_back( rsd_info ); }

private:
	utility::vector1< LKB_ResidueInfoOP > residues_info_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef LKB_ResiduesInfo LKB_PoseInfo;
typedef LKB_ResiduesInfo LKB_RotamerSetInfo;

typedef utility::pointer::shared_ptr< LKB_PoseInfo > LKB_PoseInfoOP;
typedef utility::pointer::shared_ptr< LKB_RotamerSetInfo > LKB_RotamerSetInfoOP;


}
}
}
#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_methods_LK_BallInfo )
#endif // SERIALIZATION


#endif
