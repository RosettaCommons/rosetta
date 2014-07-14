// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/GenBornPotential.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_GenBornPotential_hh
#define INCLUDED_core_scoring_GenBornPotential_hh

#include <core/scoring/GenBornPotential.fwd.hh>

#include <core/types.hh>
// AUTO-REMOVED #include <core/scoring/types.hh>

#include <basic/datacache/CacheableData.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
//#include <core//.fwd.hh>

// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/**

	 This is a reimplementation of Jim Havranek's original rosetta++ Gen Born code.
	 source files: rosetta++/gb_elec*

**/
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class GenBornResidueInfo : public utility::pointer::ReferenceCount {

public:
	typedef conformation::Residue Residue;

public:
	virtual ~GenBornResidueInfo();

	GenBornResidueInfoOP
	clone() const
	{
		return new GenBornResidueInfo( *this );
	}


	///
	GenBornResidueInfo( Residue const & rsd )
	{
		initialize( rsd );
	}

	///
	Size
	size() const
	{
		return atomic_radius_.size();
	}

	///
	Real
	atomic_radius( Size const atm ) const
	{
		return atomic_radius_[ atm ];
	}

	///
	Real &
	atomic_radius( Size const atm )
	{
		return atomic_radius_[ atm ];
	}

	///
	Real
	born_radius( Size const atm ) const
	{
		return born_radius_[ atm ];
	}

	///
	Real &
	born_radius( Size const atm )
	{
		return born_radius_[ atm ];
	}

	///
	Real
	scale_factor( Size const atm ) const
	{
		return scale_factor_[ atm ];
	}

	///
	Real &
	scale_factor( Size const atm )
	{
		return scale_factor_[ atm ];
	}

	///
	void
	initialize( Residue const & rsd );

private:
	utility::vector1< Real > atomic_radius_;
	utility::vector1< Real > born_radius_;
	utility::vector1< Real > scale_factor_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class GenBornPoseInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;

public:

	GenBornPoseInfo() {};

	GenBornPoseInfo( GenBornPoseInfo const & src );

	basic::datacache::CacheableDataOP
	clone() const
	{
		return new GenBornPoseInfo( *this );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	GenBornResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	GenBornResidueInfo const &
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
	set_placeholder( Size const i, ResidueOP rsd, GenBornResidueInfoOP info );


	///
	GenBornResidueInfo const &
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
	utility::vector1< GenBornResidueInfoOP > residue_info_;

	// these may be null pointers
	utility::vector1< ResidueOP > placeholder_residue_;
	utility::vector1< GenBornResidueInfoOP > placeholder_info_;

	// stores info from the packertask when setup_for_packing calls set_repack_list
	utility::vector1< bool > being_packed_;

};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class GenBornRotamerSetInfo : public basic::datacache::CacheableData {


public:
	typedef conformation::Residue        Residue;
	typedef conformation::ResidueOP      ResidueOP;
	typedef conformation::RotamerSetBase RotamerSetBase;

public:

	///
	GenBornRotamerSetInfo( GenBornRotamerSetInfo const & src ):
		CacheableData()
	{
		residue_info_.resize( src.size() );
		for ( Size i=1; i<= src.size(); ++i ) {
			residue_info_[i] = src.residue_info_[i]->clone();
		}
	}

	///
	basic::datacache::CacheableDataOP
	clone() const
	{
		return new GenBornRotamerSetInfo( *this );
	}

	///
	Size
	size() const
	{
		return residue_info_.size();
	}

	///
	GenBornResidueInfo &
	residue_info( Size const i )
	{
		return *residue_info_[i];
	}

	///
	GenBornResidueInfo const &
	residue_info( Size const i ) const
	{
		return *residue_info_[i];
	}


	///
	GenBornRotamerSetInfo( RotamerSetBase const & rotamer_set )
	{
		initialize( rotamer_set );
	}

	/// dont forget to 0 the born_radii
	void
	initialize( RotamerSetBase const & rotamer_set );


private:
	utility::vector1< GenBornResidueInfoOP > residue_info_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class GenBornPotential : public utility::pointer::ReferenceCount {
public:
	typedef conformation::Residue Residue;

public:
	/// ctor
	GenBornPotential():
		// Dielectric constants for protein (Ep) and for solvent (Ew)
		Ep( 4.0 ),
		Ew( 80.0 ),
		// Parameters for determination of Generalized Born radii
		ParamS( 0.09 ),
		ParamD( 1.00 ),
		ParamB( 0.80 ),
		ParamG( 4.850000 ),

		// More parameters for GB burial
		Param_TA ( 0.333333333333333 ),
		Param_TB ( 0.4 ),
		Param_TC ( 0.42857142857142857143 ),
		Param_TD ( 0.444444444444444 ),
		Param_TDD( 0.45454545454545454545 ),
		// unused Param_TE ( 1.333333333333333 ),
		// unused Param_TF ( 2.4 ),
		// unused Param_TG ( 3.42857142857142857143 ),
		// unused Param_TH ( 4.44444444444444444 ),
		// unused Param_THH( 5.45454545454545454545 ),

		dummy_radius( 2.39 ),
		dummy_scale( 0.80 ),
		dummy_distance( 2.44 )
	{}



	/// called prior to scoring, eg
	void
	get_all_born_radii( pose::Pose & pose ) const;

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
	void
	get_rotamers_born_radii(
		pose::Pose const & pose,
		conformation::RotamerSetBase & rotamer_set
	) const;

	///
	Real
	get_res_res_elecE(
		Residue const & rsd1,
		GenBornResidueInfo const & gb1,
		Residue const & rsd2,
		GenBornResidueInfo const & gb2
	) const;

	///
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

private: // these dont all *have* to be private, just not called by anyone else right now

	///
	void
	res_res_burial(
		Residue const & rsd1,
		GenBornResidueInfo & gb1,
		Residue const & rsd2,
		GenBornResidueInfo const & gb2
	) const;

	/// helper
	void
	finalize_born_radii( GenBornResidueInfo & gb_info ) const;

	///
	Real
	gb_shell_intxn(
		Real const qai,
		Real const rai,
		Real const qbi,
		Real const rbi,
		Real const dist
	) const;


	///
	Real
	gb_shell_intxn_deriv(
		Real const qai,
		Real const rai,
		Real const qbi,
		Real const rbi,
		Real const dist
	) const;


	///
	void
	get_single_rotamer_born_radii(
		Residue const & rsd1,
		pose::Pose const & pose,
		GenBornPoseInfo const & gb_info,
		GenBornResidueInfo & gb1
	) const;

	///
	void
	get_template_born_radii(
		pose::Pose const & pose,
		GenBornPoseInfo & gb_info
	) const;

	///
	void
	build_placeholders(
		pose::Pose const & pose,
		GenBornPoseInfo & gb_info
	) const;

private:

	// Dielectric constants for protein (Ep) and for solvent (Ew)
	Real const Ep;
	Real const Ew;

	// Parameters for determination of Generalized Born radii
	Real const ParamS;
	Real const ParamD;
	Real const ParamB;
	Real const ParamG;

	// More parameters for GB burial
	Real const Param_TA;
	Real const Param_TB;
	Real const Param_TC;
	Real const Param_TD;
	Real const Param_TDD;
	// unused Real const Param_TE;
	// unused Real const Param_TF;
	// unused Real const Param_TG;
	// unused Real const Param_TH;
	// unused Real const Param_THH;

	Real const dummy_radius;
	Real const dummy_scale;
	Real const dummy_distance; // also implicitly defined by the gb placeholder params file

};





} // scoring
} // core

#endif
