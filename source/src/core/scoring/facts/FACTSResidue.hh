// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @file:   core/scoring/facts/FACTSResidue.hh
// @brief:  Header-file for 3 class declartions for the FACTS algorithm
//          FACTS: Fast Analytical Continuum Treatment of Solvation by URS HABERTHUR and AMEDEO CAFLISCH
// @author: Hahnbeom Park

#ifndef INCLUDED_core_scoring_facts_FACTSResidue_HH
#define INCLUDED_core_scoring_facts_FACTSResidue_HH

#include <core/scoring/facts/FACTSPotential.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/chemical/RestypeDestructionEvent.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSRsdTypeInfo class provides all the constants and parameters for given      */
/*            residue type                                                                        */
/*                                                                                                */
/**************************************************************************************************/

class FACTSRsdTypeInfo : public utility::pointer::ReferenceCount {

public:
	void create_info( chemical::ResidueType const & rsd );

	// Accessors
	inline Size natoms() const{ return natoms_; }
	inline Real q( Size const atm ) const{ return q_[ atm ]; }
	inline Real a0(Size const atm) const{ return a0_[atm]; }
	inline Real a1(Size const atm) const{ return a1_[atm]; }
	inline Real a2(Size const atm) const{ return a2_[atm]; }
	inline Real a3(Size const atm) const{ return a3_[atm]; }
	inline Real b1(Size const atm) const{ return b1_[atm]; }
	inline Real b2(Size const atm) const{ return b2_[atm]; }
	inline Real c0(Size const atm) const{ return c0_[atm]; }
	inline Real c1(Size const atm) const{ return c1_[atm]; }
	inline Real c2(Size const atm) const{ return c2_[atm]; }
	inline Real c3(Size const atm) const{ return c3_[atm]; }
	inline Real d1(Size const atm) const{ return d1_[atm]; }
	inline Real d2(Size const atm) const{ return d2_[atm]; }

	inline Real alpha( Size const atm ) const { return alpha_[atm]; }
	inline Real COradius2(Size const atm ) const{ return COradius2_[atm]; }
	inline Real volume(Size const atm ) const{ return volume_[atm]; }
	inline bool not_using( Size const atm ) const {return not_using_[atm]; }
	inline bool charged( Size const atm ) const { return charged_[atm]; }

	//inline bool selfpair( Size const atm1, Size atm2 ) const { return selfpair_[atm1][atm2]; }
	inline bool is_chargedH( Size const atm1 ) const { return is_chargedH_[atm1]; }
	inline bool is_freedof( Size const atm1 ) const { return is_freedof_[atm1]; }

	inline Real intra_solv_scale( Size const atm1, Size const atm2 ) const { return intra_solv_scale_[atm1][atm2]; }
	inline Real intra_elec_scale( Size const atm1, Size const atm2 ) const { return intra_elec_scale_[atm1][atm2]; }

private:
	//This function initializes the vector q to charges of each atom
	void initialize_parameters( chemical::ResidueType const & rsd );
	void initialize_intrascale( chemical::ResidueType const & rsd );

private:
	Size natoms_;
	utility::vector1<bool> not_using_;
	utility::vector1<Real> q_; //list of charges of each atom
	utility::vector1<Real> COradius2_; //list of cut off radius (R^sphere_i) for calculating theta for self-energy
	utility::vector1<Real> volume_; //The volume for atoms with native van der waals
	utility::vector1<Real> alpha_;// Gamma used in equation 12 on page 706 of FACTS paper

	//For evaluating Ci for each atom (see equation 6 on page 704 of FACTS paper)
	utility::vector1<Real> b1_;
	utility::vector1<Real> b2_;

	//For evaluating esolvE_i for each atom (see equation 7 on page 704 of FACTS paper)
	utility::vector1<Real> a0_;
	utility::vector1<Real> a1_;
	utility::vector1<Real> a2_;
	utility::vector1<Real> a3_;

	//The constant for evaluating Di for each atom (see equation 10 on page 706 of FACTS paper)
	utility::vector1<Real> d1_;
	utility::vector1<Real> d2_;

	//The constant for evaluating sasa for each atom (see equation 11 on page 706 of FACTS paper)
	utility::vector1<Real> c0_;
	utility::vector1<Real> c1_;
	utility::vector1<Real> c2_;
	utility::vector1<Real> c3_;
	utility::vector1< utility::vector1<Real> > intra_solv_scale_;
	utility::vector1< utility::vector1<Real> > intra_elec_scale_;
	utility::vector1< bool > is_chargedH_;
	utility::vector1< bool > charged_;

	// Free dof atoms
	utility::vector1< bool > is_freedof_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

/// @brief The FACTSRsdTypeMap is a collection of FACTSRsdTypeInfo for a number of residue types
/// (This is a separate class to shield the raw pointer usage here
class FACTSRsdTypeMap {
public:

	FACTSRsdTypeMap();

	// If you enable copy and assignment operators, you need to account for
	// Detaching old/attaching new destruction observer to the ResidueType
	FACTSRsdTypeMap( FACTSRsdTypeMap const & ) = delete;
	FACTSRsdTypeMap & operator=( FACTSRsdTypeMap const & ) = delete;

	~FACTSRsdTypeMap();

	/// @brief get the info for the residue type, creating it if it doesn't exist
	FACTSRsdTypeInfoCOP
	get_type_info( core::chemical::ResidueType const & rsd_type );

	void restype_destruction_observer( core::chemical::RestypeDestructionEvent const & event );

private:

	// The raw pointers here are registered in the destruction observer for their respective ResidueType
	std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP > type_info_map_;

};

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSResidueInfo class provides all the functions, constants and parameters     */
/*            for different atoms, which are required to calculate the solvation free energy of   */
/*                      of a        molecule embedded in a continuum solvent using FACTS method   */
/*                                                                                                */
/**************************************************************************************************/

//This class provides the information that is specific for each atoms of a residue in FACTS algorithm
// (such as van der waals radii,  Ai_i, Bi_i, lookup_on_radius_, esolvE_i, sasa, ...)
class FACTSResidueInfo : public utility::pointer::ReferenceCount {

public:
	typedef conformation::Residue Residue;

public:
	//The function clone creates an exact copy of the FACTSResidueInfo object.
	//It returns a pointer that points to the new object
	FACTSResidueInfoOP clone() const{
		return FACTSResidueInfoOP( new FACTSResidueInfo( *this ) );
	}

	//The constructor with a residue as its input parameter.
	FACTSResidueInfo( ) { }
	FACTSResidueInfo( Residue const & rsd, FACTSRsdTypeInfoCOP restypeinfo, bool const is_rotamer = false ) {
		initialize( rsd, restypeinfo, is_rotamer );
	}

	inline FACTSRsdTypeInfoCOP restypeinfo() const {
		return restypeinfo_;
	}

	//This function initializes all the private and public variables of the class FACTSResidueInfo
	void initialize( Residue const & rsd, FACTSRsdTypeInfoCOP restypeinfo, bool const is_rotamer = false );

	void refresh_energy_cache( Size const nres );

	// logic for passing calculations for unchanged part
	void store_xyz( Residue const & rsd );
	void set_changed( bool const value ){ changed_ = value; }
	bool changed() const { return changed_; }
	void set_enumeration_shell( bool const value ){ enumeration_shell_ = value; }
	bool enumeration_shell() const { return enumeration_shell_; }

	inline Size natoms() const { return natoms_; }
	inline Real Ai( Size const atm ) const { return Ai_[atm]; }
	inline Real Bi( Size const atm ) const { return Bi_[atm]; }
	inline Real Ci( Size const atm ) const { return Ci_[atm]; }
	inline Real Di( Size const atm ) const { return Di_[atm]; }
	inline Real Ei( Size const atm ) const { return Ei_[atm]; }

	inline utility::vector1< Real > const &Ai() const { return Ai_; }
	inline utility::vector1< Vector > const &nmtr() const { return nmtr_; }
	inline utility::vector1< Real > const &dnmtr() const { return dnmtr_; }

	inline Real esolvE( Size const atm ) const { return esolvE_[atm]; }
	inline Vector nmtr( Size const atm ) const { return nmtr_[atm]; }
	inline Real dnmtr( Size const atm ) const { return dnmtr_[atm]; }
	inline Real sasa( Size const atm ) const { return sasa_[atm]; }
	inline Real BR( Size const atm ) const { return BR_[atm]; }
	inline Real E_elec( Size const res ) const { return E_elec_[res]; }
	inline Real E_solv( Size const res ) const { return E_solv_[res]; }
	inline Real E_solv_self( Size const res ) const { return E_solv_self_[res]; }
	inline Real E_solv_pair( Size const res ) const { return E_solv_pair_[res]; }

	inline Real dG_dCi( Size const atm ) const { return dG_dCi_[atm]; }
	inline Real dSA_dDi( Size const atm ) const { return dSA_dDi_[atm]; }
	inline Real dsolv_dBR( Size const atm ) const { return dsolv_dBR_[atm]; }
	inline Real dB_dBnmtr( Size const atm ) const { return dB_dBnmtr_[atm]; }
	inline Real dB_dBdnmtr( Size const atm ) const { return dB_dBdnmtr_[atm]; }
	inline Real dBR_dG( Size const atm ) const { return dBR_dG_[atm]; }
	inline Vector elecF2( Size const atm ) const { return elecF2_[atm]; }
	inline Vector solvF2d( Size const atm ) const { return solvF2d_[atm]; }
	inline Vector solvF2BR( Size const atm ) const { return solvF2BR_[atm]; }
	inline Vector sasaF2( Size const atm ) const { return sasaF2_[atm]; }
	inline utility::vector1<Vector> const &xyz() const { return xyz_; }

	inline bool flag_for_calculation( Size const atm ) const { return flag_for_calculation_[atm]; }

public: // NO DATA SHOULD BE PUBLIC!
	//const static Size MAXNEIGH = 300;
	//const static Real max_selfdcut2 = 100.0;

	Size natoms_;
	utility::vector1<Real> esolvE_; // DeltaGi (equation 7 on page 704 of FACTS paper)
	utility::vector1<Real> sasa_; // atomic SASA (equation 11 on page 706 of FACTS paper)
	utility::vector1<Real> Ai_; // Ai (equation 3 on page 704 of FACTS paper)
	utility::vector1<Real> Bi_; // Bi (equation 4 on page 704 of FACTS paper)
	utility::vector1<Real> Ci_; // Ci (equation 6 on page 704 of FACTS paper)
	utility::vector1<Real> Di_; // Di (equation 10 on page 706 of FACTS paper)
	utility::vector1<Real> Ei_; // Ei
	utility::vector1<bool> flag_for_calculation_;// this variable is used in res_res_burial and

	utility::vector1<Vector> nmtr_; // nmtr of Bi (equation 4 on page 704 of FACTS paper)
	utility::vector1<Real> dnmtr_; // dnmtr of Bi (equation 4 on page 704 of FACTS paper)
	utility::vector1<Real> BR_;    // effective Born Radius
	utility::vector1<Real> E_elec_;
	utility::vector1<Real> E_solv_;
	utility::vector1<Real> E_solv_self_;
	utility::vector1<Real> E_solv_pair_;

	utility::vector1<Vector> xyz_; // Store coordinate info
	bool changed_;
	bool enumeration_shell_;

	// Derivatives
	utility::vector1<Real> dG_dCi_;
	utility::vector1<Real> dSA_dDi_;
	utility::vector1<Real> dsolv_dBR_;
	utility::vector1<Real> dB_dBnmtr_;
	utility::vector1<Real> dB_dBdnmtr_;
	utility::vector1<Real> dBR_dG_;

	utility::vector1<Vector> elecF2_;
	utility::vector1<Vector> solvF2d_;
	utility::vector1<Vector> solvF2BR_;
	utility::vector1<Vector> sasaF2_;

private:
	FACTSRsdTypeInfoCOP restypeinfo_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

////////////////////////////////////////////////////////////////////////////////////////////////////
class FACTSRotamerSetInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;
	typedef conformation::RotamerSetBase RotamerSet;

public:
	///
	FACTSRotamerSetInfo( FACTSRotamerSetInfo const & src ): CacheableData() {
		residue_info_.resize( src.size() );
		for ( Size i=1; i<= src.size(); ++i ) {
			residue_info_[i] = src.residue_info_[i]->clone();
		}
	}

	basic::datacache::CacheableDataOP clone() const
	{ return basic::datacache::CacheableDataOP( new FACTSRotamerSetInfo( *this ) ); }

	Size size() const
	{ return residue_info_.size(); }

	FACTSResidueInfo & residue_info( Size const i )
	{ return *residue_info_[i]; }

	FACTSResidueInfo const & residue_info( Size const i ) const
	{ return *residue_info_[i]; }

	FACTSRotamerSetInfo( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap )
	{ initialize( rotamer_set, rsdtypemap ); }

	/// dont forget to 0 the born_radii
	void initialize( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap );

private:
	utility::vector1< FACTSResidueInfoOP > residue_info_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	FACTSRotamerSetInfo();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_facts_FACTSResidue )
#endif // SERIALIZATION


#endif
