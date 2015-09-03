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
// #include <core/scoring/methods/LK_BallEnergy.hh>

// // Package headers
// #include <core/scoring/methods/LK_BallEnergy.hh>
// #include <core/scoring/ScoringManager.hh>
// #include <core/scoring/NeighborList.hh>
// #include <core/scoring/EnergyGraph.hh>
// #include <core/scoring/etable/Etable.hh>
// #include <core/scoring/etable/count_pair/CountPairFunction.hh>
// #include <core/scoring/etable/count_pair/CountPairFactory.hh>
// #include <core/scoring/etable/count_pair/types.hh>

// // Project headers
#include <core/pose/Pose.fwd.hh>
// #include <core/scoring/ScoreFunction.hh>
// #include <ObjexxFCL/formatted.o.hh>
#include <core/conformation/Residue.hh>
// #include <core/conformation/ResidueFactory.hh>
// // #include <core/io/pdb/pose_io.hh> // HACK
// // #include <fstream> // HACK

// #include <core/scoring/constraints/AngleConstraint.hh>

// #include <core/options/util.hh> // HACK

// #include <core/util/prof.hh>
// #include <core/util/tracer.hh>
#include <basic/datacache/CacheableData.hh>
// #include <core/conformation/residue_datacache.hh>

// #include <numeric/constants.hh>
// #include <numeric/xyz.functions.hh>

// #include <utility/vector1.functions.hh> // HACK
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
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
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

	utility::vector1< Vectors > const &
	waters() const { return waters_; }

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
	matches_residue_type( chemical::ResidueType const & rsd_type ) { return ( rsd_type_ == &(rsd_type) ); }

	chemical::ResidueType const &
	residue_type() const { return *rsd_type_; }

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
	chemical::ResidueType const * rsd_type_; // bad form
	utility::vector1< Vectors > waters_;
	utility::vector1< utility::vector1< Real > > atom_weights_;
	bool has_waters_;

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

};

typedef LKB_ResiduesInfo LKB_PoseInfo;
typedef LKB_ResiduesInfo LKB_RotamerSetInfo;

typedef utility::pointer::shared_ptr< LKB_PoseInfo > LKB_PoseInfoOP;
typedef utility::pointer::shared_ptr< LKB_RotamerSetInfo > LKB_RotamerSetInfoOP;


}
}
}
#endif
