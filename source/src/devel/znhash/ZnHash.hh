// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/znhash/ZnHash.hh
/// @brief  Declaration of zinc-match hash for use in optimizing zinc coordination
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

#ifndef INCLUDED_devel_znhash_ZnHash_HH
#define INCLUDED_devel_znhash_ZnHash_HH

// Unit headers
#include <devel/znhash/ZnHash.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

// Numeric headers
#include <numeric/polynomial.fwd.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzVector.hh>
#include <numeric/geometry/BoundingBox.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>

// Boost headers
#include <unordered_map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace devel {
namespace znhash {

class ZnCoord {
public:
	typedef core::Vector Vector;

public:
	ZnCoord() : index_(0), nhis_( 0 ), zn_and_orbitals_( Vector( 0 ) ) {}
	ZnCoord( ZnCoord const & ) = default;

	~ZnCoord() = default;

	ZnCoord(
		core::Size index,
		core::Size nhis,
		Vector const & zn,
		Vector const & orb1,
		Vector const & orb2,
		Vector const & orb3,
		Vector const & orb4
	) :
		index_( index ),
		nhis_( nhis )
	{
		zn_and_orbitals_[ 1 ] = zn;
		zn_and_orbitals_[ 2 ] = orb1;
		zn_and_orbitals_[ 3 ] = orb2;
		zn_and_orbitals_[ 4 ] = orb3;
		zn_and_orbitals_[ 5 ] = orb4;
	}

	ZnCoord const &
	operator = ( ZnCoord const & rhs )
	{
		if ( this != &rhs ) {
			index_ = rhs.index_;
			nhis_ = rhs.nhis_;
			zn_and_orbitals_ = rhs.zn_and_orbitals_;
		} return *this;
	}

	Vector const & operator [] ( core::Size ind ) const { return zn_and_orbitals_[ ind ]; }
	core::Size index() const { return index_; }
	core::Size nhis() const { return nhis_; }
	Vector const & zn_center() const { return zn_and_orbitals_[1]; }
	Vector const & orbital( core::Size ind ) const { return zn_and_orbitals_[ind+1]; }
	core::Size norbitals() const { return 4; }
	utility::fixedsizearray1< Vector, 5 > const & zn_and_orbitals() const { return zn_and_orbitals_; }

	static core::Size natoms() { return 5; }

private:
	core::Size index_; // for identifying a particular zn
	core::Size nhis_; // how many histadines are coordinating?
	utility::fixedsizearray1< Vector, 5 > zn_and_orbitals_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class ZnHash : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~ZnHash() override;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef numeric::geometry::hashing::bin_index_hasher bin_index_hasher;
	typedef std::unordered_map< boost::uint64_t, utility::vector1< ZnCoord >, bin_index_hasher > ZnCoordinateHash;
	typedef utility::fixedsizearray1< Size, 3 > Size3;
	typedef core::Vector Vector;

public:
	ZnHash();

	/// @brief must be called before build_hash is called
	void set_uniform_bin_width( Real width );
	/// @brief First, add all zinc coordinates to the hash, then invoke build_hash.
	void add_zn_coordinate( ZnCoord const & zn );

	// @brief Build the hash after all ZnCoordinates have been added.  Must be called before query_hash().
	void build_hash();
	boost::uint64_t bin_for_point( Vector const & query_point ) const;
	ZnCoordinateHash::const_iterator query_hash( Vector const & query_point ) const;
	ZnCoordinateHash::const_iterator hash_begin() const;
	ZnCoordinateHash::const_iterator hash_end() const;

	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > const & bb() const { return bb_; }
	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > const & bb_ext() const { return bb_ext_; }

private:
	Real  grid_size_;
	Real  inv_grid_size_;
	Size3 ngrid_cells_;
	Size3 ndim_prods_;
	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > bb_; // bounding box for the zinc coordinates in the reference frame
	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > bb_ext_; // bounding box for any neighbor query

	bool hash_has_been_built_;
	std::list< ZnCoord > zn_coords_; // Coordinates in reference frame
	ZnCoordinateHash zn_hash_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class ZnMatchData
{
public:
	typedef core::Size Size;

public:
	ZnMatchData();
	~ZnMatchData();
	ZnMatchData( ZnMatchData const & src );
	ZnMatchData const & operator = ( ZnMatchData const & src );

	void index( Size setting ) { index_ = setting; }
	void res1( Size setting ) { res1_ = setting; }
	void res2( Size setting ) { res2_ = setting; }
	void zn_and_orbitals( ZnCoord const & setting ) { zn_and_orbitals_ = setting; }
	void match_pdb_file( std::string const & setting ) { match_pdb_file_ = setting; }
	void match_cst_file( std::string const & setting ) { match_cst_file_ = setting; }
	void res1conf( core::conformation::Residue const & r1 );
	void res2conf( core::conformation::Residue const & r2 );
	void znconf( core::conformation::Residue const & zn );

	Size index() const { return index_; }
	Size res1() const { return res1_; }
	Size res2() const { return res2_; }
	ZnCoord const & zn_and_orbitals() const { return zn_and_orbitals_; }
	std::string const & match_pdb_file() const { return match_pdb_file_; }
	std::string const & match_cst_file() const { return match_cst_file_; }
	core::conformation::Residue const & res1conf() const;
	core::conformation::Residue const & res2conf() const;
	core::conformation::Residue const & znconf() const;

private:

	Size index_;
	Size res1_;
	Size res2_;
	ZnCoord zn_and_orbitals_;
	std::string match_pdb_file_;
	std::string match_cst_file_;
	core::conformation::ResidueCOP res1conf_;
	core::conformation::ResidueCOP res2conf_;
	core::conformation::ResidueCOP znconf_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class ZnCoordinationScorer : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef numeric::HomogeneousTransform< Real > HTReal;
	typedef std::pair< Size, Size > ZnIndexPair;
	typedef std::pair< Real, bool > ZnScoreAndFlipState;
	typedef std::pair< ZnIndexPair, ZnScoreAndFlipState > CoordinationData;

public:
	ZnCoordinationScorer();
	~ZnCoordinationScorer() override;
	ZnCoordinationScorer( ZnCoordinationScorer const & );

	/// @brief set an individual atom ID to use for determining the refrence frame for the asymmetric unit
	/// There are three atoms total which are used to define the reference frame.  By default, these are
	/// initialized to atoms 1,2,&3 on reisude 1.
	void set_asymm_atid( Size which_atid, core::id::AtomID atid );

	/// @brief set an individual atom ID to use for determining the reference frame for the symmetric clone.
	void set_symm_atid( Size which_atid, core::id::AtomID atid );

	/// @brief set the residue on the asymmetric unit and use atoms 1,2,&3 on that residue to define the reference
	/// frame.
	void set_asymm_resid( Size resid );

	/// @brief set the residue on the symmetric clone and use atoms 1,2,&3 on that residue to define the reference
	/// frame.
	void set_symm_resid( Size resid );

	/// @brief Set a third residue (not the asymm resid, nor the symm resid) so that this can be treated as a three-body
	/// interaction and therefore scored during "finalize" instead of
	/// with the two-body energies.
	void set_third_resid( Size resid );

	void set_zn_reach( Real reach );
	void set_orbital_dist( Real dist );
	void set_orbital_reach( Real reach );
	void set_zn_well_depth( Real depth );
	/// @brief If true, then only consider matches where three of
	/// the four coordinating residues are histadine
	void require_3H( bool setting ) { require_3H_ = setting; }

	/// @brief This is used to generate the reference frame.  Set the asymm_atids first (if necessary)
	/// set the reference pdb next, finally add the matches.
	void set_reference_pdb( std::string const & start_pdb );
	void set_matcher_constraint_file_name( std::string const & fname );

	void add_match_from_file(
		std::string const & match_file_name
	);

	void add_matches_from_files(
		std::list< std::string > const & fnames
	);

	void add_match_from_istream( std::istream & input_match );

	utility::fixedsizearray1< core::id::AtomID, 3 > const &
	asymm_atids() const {
		return asymm_atids_;
	}

	utility::fixedsizearray1< core::id::AtomID, 3 > const &
	focused_clone_atids() const {
		return focused_clone_atids_;
	}

	core::id::AtomID const &
	atom_one_on_third_residue() const {
		return third_resid_;
	}

	utility::vector1< ZnMatchData > const & zn_matches() const { return zn_matches_; }
	void set_clash_weight( Real weight ) { clash_weight_ = weight; }
	Real clash_weight() const { return clash_weight_; }
	void set_idealize_input_virtual_atoms( bool setting ) { idealize_input_virtual_atoms_ = setting; }

	void finalize_after_all_matches_added();

	bool optimal_coordination_is_reversed( core::pose::Pose const & p );

	Real score( core::pose::Pose const & p ) const;

	Real score(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	) const;

	Real score_one_zn( ZnCoord const & zn ) const;

	ZnScoreAndFlipState
	score_zn_pair( ZnCoord const & zn1, ZnCoord const & zn2 ) const;

	// @brief Compute the homogeneous transform that will take Zn coordinates from the reference frame,
	// translate them onto the the location of the coordinates on the slave clone, and then
	// translate them further into their effective location in the original reference frame given
	// the new location of the master clone.  Pretty simple, really.  Public for testing purposes only.
	HTReal
	query_frame_to_original_frame(
		core::conformation::Residue const & r1,
		core::conformation::Residue const & r2
	) const;

	HTReal query_frame_to_original_frame( core::pose::Pose const & pose ) const;

	ZnCoord original_frame_coordinate_for_match( Size index, core::pose::Pose const & p ) const;

	/// @ brief score and return the indices of the two Zn matches that produce the best score;
	/// If the result is "p", then the indices are ordered so that p.first is the Zn from
	/// the asymmetric unit, and p.second is the Zn from the symmetric clone.
	ZnIndexPair best_match( core::pose::Pose const & pose ) const;

	ZnIndexPair
	best_match(
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	) const;

	/// @ brief score and return the indices of the two Zn matches that produce the best score;
	/// and the score they produce.  If the result is "p", then p.first is the pair of zn indices,
	/// and p.second is the score.  The indices are ordered so that p.first.first is the Zn from
	/// the asymmetric unit, and p.first.second is the Zn from the symmetric clone
	CoordinationData
	score_and_index_for_best_match(
		core::conformation::Residue const & res1, // <-- from the asymmetric unit
		core::conformation::Residue const & res2 // <-- from the symmetric clone
	) const;

	Real
	clash_score(
		Size m1,
		Size m2,
		HTReal const & to_hash_frame_transform,
		utility::vector1< core::Vector > & symmclone_coords_res1,
		utility::vector1< core::Vector > & symmclone_coords_res2
	) const;

	Real
	clash_score_residue_pair(
		core::conformation::Residue const & r1, // <- from asymm unit
		core::conformation::Residue const & r2, // <- from symm clone
		utility::vector1< core::Vector > const & r2coords // <- transformed coordinates
	) const;

	void
	insert_match_onto_pose(
		core::pose::Pose & p, // <-- should be a symmetric pose if this match is to be applied to all clones
		Size match_index,
		Size chain_insertion_id // 1 for the asymm unit, 2 for the symm clone.
	) const;

public:
	/// @brief quick access to the residue defining the coordinate frame on the
	/// asymmetric unit (residue 1)
	core::Size r1() const { return asymm_atids_[1].rsd(); }

	/// @brief quick access to the residue defining the coordinate frame on the
	/// focused clone (residue 2)
	core::Size r2() const { return focused_clone_atids_[1].rsd(); }

private:
	void reset_reach();
	void reset_znx_orbital_coords();

private:

	/// The coordinates of the zinc and its virtual atoms in idealized geometry.
	utility::vector1< core::Vector > znx_ideal_coords_;
	std::string matchcst_file_name_; // string naming the match geometric constraint file used for all matches

	numeric::Polynomial_1dOP ramp_to_zero_poly_;
	Real well_depth_; // multiplied by the polynomial
	Real reach_; // grid size / threshold distance between Zn centers
	Real reach2_; // reach^2

	HTReal reference_frame_; // built using the asymm_atids_ on the match PDB.
	HTReal inv_reference_frame_; // the inverse reference frame
	bool idealize_input_virtual_atoms_;
	utility::vector1< ZnMatchData > zn_matches_;
	Size max_natoms_;
	Real clash_weight_;
	bool require_3H_;


	Real znreach_; // how far does the zn score extend out to?
	Real orbital_dist_; // how far from the Zn center should the orbitals be located?
	Real orbital_reach_; // at what point are the orbitals not within striking distance?
	Real znwelldepth_; // how much should zn proximity count relative to orbital proximity?  <-- scale the whole score by modifying the metalhash_constraint weight.

	Size asymm_chain_;
	Size focsed_clone_chain_;
	utility::fixedsizearray1< core::id::AtomID, 3 > asymm_atids_; // build a coordinate frame on the asymmetric unit with these three atoms
	utility::fixedsizearray1< core::id::AtomID, 3 > focused_clone_atids_; // build a coordinate frame on the focused clone with these three atoms
	core::id::AtomID third_resid_;

	ZnHashOP hash_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class ZnCoordinationConstraint : public core::scoring::constraints::Constraint
{
public:
	ZnCoordinationConstraint( ZnCoordinationScorerCOP zn_score );
	~ZnCoordinationConstraint() override;

	
	core::scoring::constraints::ConstraintOP
	clone() const override;

	bool operator == ( Constraint const & other ) const override;
	bool same_type_as_me( Constraint const & other ) const override;

	/// @brief Returns the number of atoms involved in defining this constraint.
	
	core::Size
	natoms() const override;

	/// @brief Returns the AtomID referred to by index.
	
	core::id::AtomID const &
	atom( Size const index ) const override;

	
	void
	score(
		core::scoring::func::XYZ_Func const & xyz_func,
		core::scoring::EnergyMap const & weights,
		core::scoring::EnergyMap & emap ) const override;

	/// @brief Fill the f1 and f2 vectors, necessary for considering the
	/// derivative this constraint during minimization. (someone please reference
	/// Bill Wedermeyer's paper here, as I'm in an airport and can't fill it in
	/// myself!)
	
	void
	fill_f1_f2(
		core::id::AtomID const & atom,
		core::scoring::func::XYZ_Func const & xyz_func,
		core::Vector & F1,
		core::Vector & F2,
		core::scoring::EnergyMap const & weights
	) const override;

	void show( std::ostream & /*out*/ ) const override {}

	core::scoring::constraints::ConstraintOP
	remap_resid( core::id::SequenceMapping const &/*seqmap*/ ) const override
	{
		return core::scoring::constraints::ConstraintOP( new ZnCoordinationConstraint( zn_score_ ) );
	}


private:
	/// @brief quick access to the residue defining the coordinate frame on the
	/// asymmetric unit (residue 1)
	core::Size r1() const { return zn_score_->asymm_atids()[1].rsd(); }

	/// @brief quick access to the residue defining the coordinate frame on the
	/// focused clone (residue 2)
	core::Size r2() const { return zn_score_->focused_clone_atids()[1].rsd(); }

private:
	ZnCoordinationScorerCOP zn_score_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ZnCoordinationConstraint();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( devel_znhash_ZnHash )
#endif // SERIALIZATION


#endif
