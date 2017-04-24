// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/chemical/RestypeDestructionEvent.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/thread/ReadWriteMutex.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace lkball {

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
typedef utility::vector1< WaterBuilders > WaterBuildersList;

typedef utility::vector1< utility::vector1< Real > >  AtomWeights;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief A singleton class which stores data for LKBall terms.
/// This is a separate singleton class, rather than static data on the LKB_ResidueInfo class
/// so that the ResidueType destruction observer has a stable object to call back to.
class LKBallDatabase : public utility::SingletonBase< LKBallDatabase >
{
public:
	friend class utility::SingletonBase< LKBallDatabase >;

	~LKBallDatabase();

	/// @brief Returns true if the passed rsd_type is in the database
	bool
	has( chemical::ResidueType const & rsd_type ) const;

	void
	initialize_residue_type( chemical::ResidueType const & rsd_type );

	/// @brief Get the water builders for the ResidueType
	/// A vector of one per heavy atom.
	/// @details IMPORTANT LIFETIME NOTE: The return value is only valid so long as the passed residue type is
	WaterBuildersList const &
	get_water_builders( chemical::ResidueType const & rsd_type ) const;

	/// @brief Get the atom_weights for the ResidueType
	/// @details IMPORTANT LIFETIME NOTE: The return value is only valid so long as the passed residue type is
	AtomWeights const &
	get_atom_weights( chemical::ResidueType const & rsd_type ) const;

	/// danger
	void
	reset_arrays_danger_expert_only();

private:

	void
	setup_atom_weights(
		chemical::ResidueType const & rsd_type,
		WaterBuildersList const & rsd_water_builders, // for sanity
		AtomWeights & atom_wts
	);

	void restype_destruction_observer( core::chemical::RestypeDestructionEvent const & event );

private:

	/// @brief private constructor
	LKBallDatabase();
	LKBallDatabase( LKBallDatabase const & ) = delete;
	LKBallDatabase & operator=( LKBallDatabase const & ) = delete;

private:

	// The raw pointers here are registered in the destruction observer for their respective ResidueType
	typedef std::map< chemical::ResidueType const *, WaterBuildersList > WaterBuilderMap;
	typedef std::map< chemical::ResidueType const *, AtomWeights > AtomWeightsMap;

	WaterBuilderMap water_builder_map_;
	AtomWeightsMap atom_weights_map_;

#ifdef MULTI_THREADED
	static utility::thread::ReadWriteMutex lkball_db_mutex_;
#endif

};

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


	bool
	matches_residue_type( chemical::ResidueType const & rsd_type ) const;

	chemical::ResidueType const &
	residue_type() const;

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
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_lkball_LK_BallInfo )
#endif // SERIALIZATION


#endif
