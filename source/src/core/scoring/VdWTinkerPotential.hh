// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/VdWTinkerPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_VdWTinkerPotential_hh
#define INCLUDED_core_scoring_VdWTinkerPotential_hh

#include <core/scoring/VdWTinkerPotential.fwd.hh>

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


class VdWTinkerResidueInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue Residue;
	typedef  numeric::xyzMatrix< Real > Matrix;
	typedef  numeric::xyzVector< Real > Vector;

public:
	virtual ~VdWTinkerResidueInfo();

	VdWTinkerResidueInfoOP
	copy_clone() const
	{
		return core::scoring::VdWTinkerResidueInfoOP( new VdWTinkerResidueInfo( *this ) );
	}

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new VdWTinkerResidueInfo( *this ) );
	}

	///
	VdWTinkerResidueInfo(){}

	///
	VdWTinkerResidueInfo( Residue const & rsd )
	{
		initialize( rsd );
	}

	///
	Size
	vdw_type( Size const atm ) const
	{
		return vdw_type_[ atm ];
	}

	///
	Size &
	nonconst_vdw_type( Size const atm )
	{
		return vdw_type_[ atm ];
	}

	///
	void
	initialize( Residue const & rsd );

private:
	utility::vector1< Size > vdw_type_;
};

typedef utility::pointer::shared_ptr< VdWTinkerResidueInfo > VdWTinkerResidueInfoOP;
typedef utility::pointer::shared_ptr< const VdWTinkerResidueInfo > VdWTinkerResidueInfoCOP;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class VdWTinkerPoseInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;

public:

	VdWTinkerPoseInfo() {};

	VdWTinkerPoseInfo( VdWTinkerPoseInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return basic::datacache::CacheableDataOP( new VdWTinkerPoseInfo( *this ) );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	VdWTinkerResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	VdWTinkerResidueInfo const &
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
	set_placeholder( Size const i, ResidueOP rsd, VdWTinkerResidueInfoOP info );


	///
	VdWTinkerResidueInfo const &
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
	utility::vector1< VdWTinkerResidueInfoOP > residue_info_;

	// these may be null pointers
	utility::vector1< ResidueOP > placeholder_residue_;
	utility::vector1< VdWTinkerResidueInfoOP > placeholder_info_;

	// stores info from the packertask when setup_for_packing calls set_repack_list
	utility::vector1< bool > being_packed_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class VdWTinkerRotamerSetInfo : public basic::datacache::CacheableData {


public:
	typedef conformation::Residue        Residue;
	typedef conformation::ResidueOP      ResidueOP;
	typedef conformation::RotamerSetBase RotamerSetBase;

public:

	///
	VdWTinkerRotamerSetInfo( VdWTinkerRotamerSetInfo const & src ):
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
		return basic::datacache::CacheableDataOP( new VdWTinkerRotamerSetInfo( *this ) );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	VdWTinkerResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	VdWTinkerResidueInfo const &
	residue_info( Size const i ) const
	{
		return *residue_info_[i];
	}


	///
	VdWTinkerRotamerSetInfo( RotamerSetBase const & rotamer_set )
	{
		initialize( rotamer_set );
	}

	/// dont forget to 0 the born_radii
	void
	initialize( RotamerSetBase const & rotamer_set );


private:
	utility::vector1< VdWTinkerResidueInfoOP > residue_info_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



class VdWTinkerPotential : public utility::pointer::ReferenceCount {
public:
	typedef core::conformation::Residue Residue;

public:
	/// ctor
	VdWTinkerPotential()
	{ read_in_amoeba_parameters();
		read_in_vdw_tinker_parameters();
	}


	/// read in parameters for amoeba and mappings between
	/// Rosetta residues/atoms and Amoeba types
	void
	read_in_amoeba_parameters();

	/// read in vdw parameters for amoeba types
	void
	read_in_vdw_tinker_parameters();

	/// Look up Amoeba type by resname/atomname/variant name
	core::Size
	amoeba_type_lookup(
		std::string const & atomname,
		std::string const & resname,
		std::string const & variantname
	) const;

	/// called prior to scoring, eg
	void
	assign_residue_amoeba_type( Residue const & rsd, VdWTinkerResidueInfo & mp ) const;

	/// called prior to scoring, eg
	void
	assign_all_amoeba_types( pose::Pose & pose ) const;

	/// Get the amoeba info for rotamers
	void
	get_rotamers_vdw_info( pose::Pose const & pose, conformation::RotamerSetBase & ) const;

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
	get_res_res_vdw(
		Residue const & rsd1,
		VdWTinkerResidueInfo const & mp1,
		Residue const & rsd2,
		VdWTinkerResidueInfo const & mp2
	) const;

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
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		VdWTinkerResidueInfo const & mp1,
		VdWTinkerResidueInfo const & mp2,
		pose::Pose const & pose, // provides context
		Real const & factor,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

public:
private:

	std::map< std::string, Size > type_lookup_;
	std::string const default_variant_;
	utility::vector1< Real > vdw_radius_;
	utility::vector1< Real > vdw_depth_;
	utility::vector1< Real > vdw_reduce_;

};





} // scoring
} // core

#endif
