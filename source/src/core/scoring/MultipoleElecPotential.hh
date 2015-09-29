// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/MultipoleElecPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_MultipoleElecPotential_hh
#define INCLUDED_core_scoring_MultipoleElecPotential_hh

#include <core/scoring/MultipoleElecPotential.fwd.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>

#include <basic/datacache/CacheableData.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/DomainMap.fwd.hh>
//#include <core//.fwd.hh>

// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


class MultipoleElecResidueInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue Residue;
	typedef  numeric::xyzMatrix< Real > Matrix;
	typedef  numeric::xyzVector< Real > Vector;

public:
	virtual ~MultipoleElecResidueInfo();

	MultipoleElecResidueInfoOP
	copy_clone() const
	{
		return core::scoring::MultipoleElecResidueInfoOP( new MultipoleElecResidueInfo( *this ) );
	}

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new MultipoleElecResidueInfo( *this ) );
	}

	///
	MultipoleElecResidueInfo(){}

	///
	MultipoleElecResidueInfo( Residue const & rsd )
	{
		initialize( rsd );
	}

	///
	utility::vector1< Size > const &
	type() const
	{
		return type_;
	}

	///
	Size
	type( Size const atm ) const
	{
		return type_[ atm ];
	}

	///
	void
	set_type( Size const atm, Size in_val )
	{
		type_[ atm ] = in_val;
	}

	///
	utility::vector1< id::AtomID > const &
	coord_frame_ref() const
	{
		return coord_frame_ref_;
	}

	///
	id::AtomID
	coord_frame_ref( Size const atm ) const
	{
		return coord_frame_ref_[ atm ];
	}

	///
	void
	set_coord_frame_ref( Size const atm, id::AtomID in_val )
	{
		coord_frame_ref_[ atm ] = in_val;
	}

	///
	Real
	monopole( Size const atm ) const
	{
		return monopole_[ atm ];
	}

	///
	Real
	rKirkwood( Size const atm ) const
	{
		return rKirkwood_[ atm ];
	}

	///
	Real &
	nonconst_rKirkwood( Size const atm )
	{
		return rKirkwood_[ atm ];
	}

	///
	Vector const &
	dipole( Size const atm ) const
	{
		return dipole_[ atm ];
	}

	///
	Vector const &
	induced_dipole( Size const atm ) const
	{
		return induced_dipole_[ atm ];
	}

	///
	Vector const &
	stored_induced_dipole( Size const atm ) const
	{
		return stored_induced_dipole_[ atm ];
	}

	///
	Vector const &
	induced_rf_dipole( Size const atm ) const
	{
		return induced_rf_dipole_[ atm ];
	}

	///
	Vector const &
	stored_induced_rf_dipole( Size const atm ) const
	{
		return stored_induced_rf_dipole_[ atm ];
	}

	///
	Matrix const &
	quadrupole( Size const atm ) const
	{
		return quadrupole_[ atm ];
	}

	///
	Real &
	nonconst_monopole( Size const atm )
	{
		return monopole_[ atm ];
	}

	///
	Vector &
	nonconst_stored_induced_dipole( Size const atm )
	{
		return stored_induced_dipole_[ atm ];
	}

	///
	Vector &
	nonconst_induced_dipole( Size const atm )
	{
		return induced_dipole_[ atm ];
	}

	///
	Vector &
	nonconst_stored_induced_rf_dipole( Size const atm )
	{
		return stored_induced_rf_dipole_[ atm ];
	}

	///
	Vector &
	nonconst_induced_rf_dipole( Size const atm )
	{
		return induced_rf_dipole_[ atm ];
	}

	///
	Vector &
	Efield_fixed( Size const atm )
	{
		return Efield_fixed_[ atm ];
	}

	///
	utility::vector1< Vector > &
	Efield_fixed()
	{
		return Efield_fixed_;
	}

	///
	Vector &
	Efield_induced( Size const atm )
	{
		return Efield_induced_[ atm ];
	}

	///
	utility::vector1< Vector > &
	Efield_induced()
	{
		return Efield_induced_;
	}

	///
	Vector &
	Efield_rf_fixed( Size const atm )
	{
		return Efield_rf_fixed_[ atm ];
	}

	///
	utility::vector1< Vector > &
	Efield_rf_fixed()
	{
		return Efield_rf_fixed_;
	}

	///
	Vector &
	Efield_rf_induced( Size const atm )
	{
		return Efield_rf_induced_[ atm ];
	}

	///
	utility::vector1< Vector > &
	Efield_rf_induced()
	{
		return Efield_rf_induced_;
	}

	///
	Vector &
	nonconst_dipole( Size const atm )
	{
		return dipole_[ atm ];
	}

	///
	Matrix &
	nonconst_quadrupole( Size const atm )
	{
		return quadrupole_[ atm ];
	}

	///
	utility::vector1< id::AtomID > &
	my_group( Size const atm )
	{
		return my_group_[ atm ];
	}

	///
	utility::vector1< id::AtomID > const &
	const_my_group( Size const atm ) const
	{
		return my_group_[ atm ];
	}

	utility::vector1< id::AtomID > &
	my_local_coord_frame( Size const atm )
	{
		return my_local_coord_frame_[ atm ];
	}

	Matrix &
	local_coord_matrix( Size const atm )
	{
		return local_coord_matrix_[ atm ];
	}

	///
	MultipoleParameterOP const &
	mp_param( Size atm1 ) const
	{
		return mp_param_[atm1];
	}

	MultipoleParameterOP &
	nonconst_mp_param( Size atm1 )
	{
		return mp_param_[atm1];
	}

	///
	std::string const &
	rosetta_res_type() const
	{
		return rosetta_res_type_;
	}

	std::string &
	nonconst_rosetta_res_type()
	{
		return rosetta_res_type_;
	}

	///
	void
	initialize( Residue const & rsd );

private:
	utility::vector1< Size > type_;
	utility::vector1< id::AtomID > coord_frame_ref_;
	utility::vector1< Real > monopole_;
	utility::vector1< Real > rKirkwood_;
	utility::vector1< Vector > dipole_;
	utility::vector1< Vector > induced_dipole_;
	utility::vector1< Vector > stored_induced_dipole_;
	utility::vector1< Vector > induced_rf_dipole_;
	utility::vector1< Vector > stored_induced_rf_dipole_;
	utility::vector1< Vector > Efield_fixed_;
	utility::vector1< Vector > Efield_induced_;
	utility::vector1< Vector > Efield_rf_fixed_;
	utility::vector1< Vector > Efield_rf_induced_;
	utility::vector1< Matrix > quadrupole_;
	utility::vector1< Matrix > local_coord_matrix_;
	utility::vector1< utility::vector1< id::AtomID > > my_group_;
	utility::vector1< utility::vector1< id::AtomID > > my_local_coord_frame_;
	utility::vector1< MultipoleParameterOP > mp_param_;
	std::string rosetta_res_type_;

};

typedef utility::pointer::shared_ptr< MultipoleElecResidueInfo > MultipoleElecResidueInfoOP;
typedef utility::pointer::shared_ptr< const MultipoleElecResidueInfo > MultipoleElecResidueInfoCOP;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class MultipoleElecPoseInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;

public:

	MultipoleElecPoseInfo() {};

	MultipoleElecPoseInfo( MultipoleElecPoseInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new MultipoleElecPoseInfo( *this ) );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	MultipoleElecResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	MultipoleElecResidueInfo const &
	residue_info( Size const i ) const
	{
		return *residue_info_[i];
	}

	///
	bool
	being_packed( Size const seqpos ) const
	{
		return being_packed_[ seqpos ];
	}

	///
	void
	set_placeholder( Size const i, ResidueOP rsd, MultipoleElecResidueInfoOP info );


	///
	MultipoleElecResidueInfo const &
	placeholder_info( Size const seqpos ) const
	{
		assert( placeholder_info_[ seqpos ] );
		return *placeholder_info_[ seqpos ];
	}

	///
	Residue const &
	placeholder_residue( Size const seqpos ) const
	{
		assert( placeholder_residue_[ seqpos ] );
		return *placeholder_residue_[ seqpos ];
	}

	///
	void
	initialize( pose::Pose const & pose );


	///
	void
	set_repack_list( utility::vector1< bool > const & repacking_residues );

private:

	// these are allocated in initialize
	utility::vector1< MultipoleElecResidueInfoOP > residue_info_;

	// these may be null pointers
	utility::vector1< ResidueOP > placeholder_residue_;
	utility::vector1< MultipoleElecResidueInfoOP > placeholder_info_;

	// stores info from the packertask when setup_for_packing calls set_repack_list
	utility::vector1< bool > being_packed_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class MultipoleElecRotamerSetInfo : public basic::datacache::CacheableData {


public:
	typedef conformation::Residue        Residue;
	typedef conformation::ResidueOP      ResidueOP;
	typedef conformation::RotamerSetBase RotamerSetBase;

public:

	///
	MultipoleElecRotamerSetInfo( MultipoleElecRotamerSetInfo const & src ):
		CacheableData()
	{
		residue_info_.resize( src.size() );
		for ( Size i=1; i<= src.size(); ++i ) {
			residue_info_[i] = src.residue_info_[i]->copy_clone();
		}
	}

	///
	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new MultipoleElecRotamerSetInfo( *this ) );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	MultipoleElecResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	MultipoleElecResidueInfo const &
	residue_info( Size const i ) const
	{
		return *residue_info_[i];
	}


	///
	MultipoleElecRotamerSetInfo( RotamerSetBase const & rotamer_set )
	{
		initialize( rotamer_set );
	}

	/// dont forget to 0 the born_radii
	void
	initialize( RotamerSetBase const & rotamer_set );


private:
	utility::vector1< MultipoleElecResidueInfoOP > residue_info_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// Enum for coordinate axis definition types

enum MultipoleAxisType {
	none=1,
	z_axis_only,
	three_fold,
	bisector,
	z_then_bisector,
	z_then_x
};


class MultipoleParameter : public utility::pointer::ReferenceCount {

public:
	typedef utility::pointer::shared_ptr< MultipoleParameter > MultipoleParameterOP;
	typedef numeric::xyzVector< Real > Vector;
	typedef numeric::xyzMatrix< Real > Matrix;

	virtual ~MultipoleParameter(){};

	MultipoleParameter( MultipoleAxisType coord_type_in,
		utility::vector1< Size > & atom_type_in,
		Real chirality_sign_in,
		Real monopole_in,
		Vector & dipole_in,
		Matrix & quadrupole_in ) :
		coord_type_( coord_type_in ),
		atom_type_( atom_type_in ),
		chirality_sign_( chirality_sign_in ),
		monopole_( monopole_in ),
		dipole_( dipole_in ),
		quadrupole_( quadrupole_in ) {}

	MultipoleParameterOP
	clone() const
	{
		return MultipoleParameterOP( new MultipoleParameter( *this ) );
	}

	MultipoleAxisType & coord_type() { return coord_type_; }
	utility::vector1< Size > & atom_type() { return atom_type_; }
	Real & chirality_sign() { return chirality_sign_; }
	Real & monopole() { return monopole_; }
	Vector & dipole() { return dipole_; }
	Matrix & quadrupole() { return quadrupole_; }
	Real & polarity() { return polarity_; }
	Real & thole() { return thole_; }
	Real & pdamp() { return pdamp_; }
	utility::vector1< Size > & my_group_members() { return group_members_; }

private:
	MultipoleAxisType coord_type_;

	utility::vector1< Size > atom_type_;

	Real chirality_sign_;

	Real monopole_;
	Vector dipole_;
	Matrix quadrupole_;

	// The following are for polarizable dipoles
	Real polarity_;
	Real thole_;
	Real pdamp_;
	utility::vector1< Size > group_members_;

};




class MultipoleElecPotential : public utility::pointer::ReferenceCount {
public:
	typedef core::conformation::Residue Residue;

public:
	/// ctor
	MultipoleElecPotential():
		// Dielectric constants for protein (Ep) and for solvent (Ew)
		Ep( 1.0 ),
		Ew( 78.3 ),
		bohr( 0.52917721092 ),
		use_polarization( true ),
		use_gen_kirkwood( false ),
		default_variant_( "NONE" )
	{ read_in_amoeba_parameters();
		read_in_multipole_parameters();
	}


	/// read in parameters for amoeba and mappings between
	/// Rosetta residues/atoms and Amoeba types
	void
	read_in_amoeba_parameters();

	/// read in multipole parameters for amoeba types
	void
	read_in_multipole_parameters();


	/// Find the appropriate multipole params and axis atoms
	void
	find_params_and_neighbors(
		core::pose::Pose const & pose,
		MultipoleParameterOP & mp_param,
		MultipoleElecResidueInfo & mp,
		core::conformation::Residue const & rsd,
		Size const j,
		Size const this_type
	) const;

	/// Build the matrix to rotate from the local frame to the
	/// global one.
	void
	build_frame_and_rotate(
		core::pose::Pose const & pose,
		MultipoleParameterOP & mp_param,
		Size orig_atom,
		MultipoleElecResidueInfo & mp,
		core::conformation::Residue const & rsd
	) const;

	/// Look up Amoeba type by resname/atomname/variant name
	core::Size
	amoeba_type_lookup(
		std::string const & atomname,
		std::string const & resname,
		std::string const & variantname
	) const;

	/// called prior to scoring, eg
	void
	align_residue_multipole_axes( core::pose::Pose const & pose, Residue const & rsd, MultipoleElecResidueInfo & mp) const;

	/// called prior to scoring, eg
	void
	align_multipole_axes( pose::Pose & pose ) const;

	/// called prior to scoring, eg
	void
	assign_residue_amoeba_type( Residue const & rsd, MultipoleElecResidueInfo & mp ) const;

	/// called prior to scoring, eg
	void
	assign_all_amoeba_types( pose::Pose & pose ) const;

	/// called prior to scoring, eg
	void
	determine_polarization_groups( pose::Pose & pose ) const;

	/// called prior to scoring, eg
	void
	induce_polarizable_dipoles( pose::Pose & pose ) const;

	/// Copies over induced dipoles to storage
	void
	store_induced_dipoles( pose::Pose & pose ) const;

	/// Copies over induced dipoles to storage
	core::Real
	relax_induced_dipoles( pose::Pose & pose, Real relax ) const;

	/// Get the electric field due to permanent multipoles
	void
	calculate_fixed_fields_for_polarization( pose::Pose & pose ) const;

	/// Zero out the fields due to induced dipoles
	void
	clear_induced_fields( pose::Pose & pose ) const;

	/// Get the electric field due to permanent multipoles
	void
	calculate_induced_fields_for_polarization( pose::Pose & pose ) const;

	/// Get the electric field due to permanent multipoles
	void
	get_polarization_from_fields( pose::Pose & pose ) const;

	/// Get the effective radii for Generalized Kirkwood
	void
	get_effective_radii( pose::Pose & pose ) const;

	/// Get the amoeba info for rotamers
	void
	get_rotamers_multipole_info( pose::Pose const & pose, conformation::RotamerSetBase & ) const;

	/// Get the effective radii for Generalized Kirkwood
	void
	get_rotamers_effective_radii( pose::Pose const & pose, conformation::RotamerSetBase & ) const;

	/// Get the effective radii for Generalized Kirkwood
	void
	get_single_rotamer_effective_radii(
		Residue const & rsd1,
		pose::Pose const & pose,
		MultipoleElecPoseInfoCOP mp_info,
		MultipoleElecResidueInfo & mp1
	) const;

	///
	void
	setup_for_scoring(
		pose::Pose & pose
	) const;

	///
	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const & repacking_residues
	) const;


	///
	void
	update_residue_for_packing(
		pose::Pose & pose,
		Size const seqpos
	) const;

	///
	Real
	get_res_res_elecE(
		Residue const & rsd1,
		MultipoleElecResidueInfo const & mp1,
		Residue const & rsd2,
		MultipoleElecResidueInfo const & mp2
	) const;

	///
	void
	calculate_res_res_fixed_fields_for_polarization(
		Residue const & rsd1,
		MultipoleElecResidueInfo & mp1,
		Residue const & rsd2,
		MultipoleElecResidueInfo & mp2
	) const;

	///
	void
	calculate_res_res_induced_fields_for_polarization(
		Residue const & rsd1,
		MultipoleElecResidueInfo & mp1,
		Residue const & rsd2,
		MultipoleElecResidueInfo & mp2
	) const;

	///
#ifdef NOTDEF
	void
	eval_atom_derivative(
		id::AtomID const & id,
		Real const weight,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		bool const exclude_DNA_DNA,
		Vector & F1,
		Vector & F2
	) const;
#endif

	void
	calculate_and_store_all_derivs(
		pose::Pose const & pose
	) const;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		MultipoleElecResidueInfo const & mp1,
		MultipoleElecResidueInfo const & mp2,
		pose::Pose const & pose, // provides context
		Real const & factor,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

public:
	// Dielectric constants for protein (Ep) and for solvent (Ew)
	Real Ep;
	Real Ew;
	// Conversion from bohrs to angstroms
	Real const bohr;
	bool use_polarization;
	bool use_gen_kirkwood;
private:

	std::map< std::string, Size > type_lookup_;
	std::string const default_variant_;
	std::multimap< Size, MultipoleParameter::MultipoleParameterOP > multipole_info_;
	mutable utility::vector1< utility::vector1< DerivVectorPair > > cached_atom_derivs_;

};





} // scoring
} // core

#endif
