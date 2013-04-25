// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPotential.hh
// @brief:  Header-file for 3 class declartions for the FACTS algorithm
//          FACTS: Fast Analytical Continuum Treatment of Solvation by URS HABERTHUR and AMEDEO CAFLISCH
// @author: Hahnbeom Park

#ifndef INCLUDED_core_scoring_facts_FACTSPotential_HH
#define INCLUDED_core_scoring_facts_FACTSPotential_HH

// Unit headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <iostream>

using namespace std;

namespace core {
namespace scoring {

// Forward declaration
class FACTSRsdTypeInfo;

// Structure for storing Self-Neighbor(not GB-Neighbor) info
struct SelfNeighInfo {
	Size nneigh;
	utility::vector1<Size> resID;
	utility::vector1<Size> atmID;
};

typedef utility::pointer::owning_ptr< FACTSRsdTypeInfo > FACTSRsdTypeInfoOP;
typedef utility::pointer::owning_ptr< FACTSRsdTypeInfo const > FACTSRsdTypeInfoCOP;
typedef std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP > FACTSRsdTypeMap;

/****************************************************************************************************/
/**                                                                                                **/
/**    @brief: The FACTSRsdTypeInfo class provides all the constants and parameters for given      **/
/**            residue type                                                                        **/
/**                                                                                                **/
/****************************************************************************************************/
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

	inline bool selfpair( Size const atm1, Size atm2 ) const { return selfpair_[atm1][atm2]; }
	inline bool is_chargedH( Size const atm1 ) const { return is_chargedH_[atm1]; }

private:
	//This function initializes the vector q to charges of each atom
	void initialize_parameters( chemical::ResidueType const & rsd );
	void initialize_selfpair( chemical::ResidueType const & rsd );

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
	utility::vector1< utility::vector1<bool> > selfpair_;
	utility::vector1< bool > is_chargedH_;
};

//This class provides the information that is specific for each atoms of a residue in FACTS algorithm
// (such as van der waals radii,  Ai_i, Bi_i, lookup_on_radius_, esolvE_i, sasa, ...)
/****************************************************************************************************/
/**                                                                                                **/
/**    @brief: The FACTSResidueInfo class provides all the functions, constants and parameters     **/
/**            for different atoms, which are required to calculate the solvation free energy of   **/
/**                      of a        molecule embedded in a continuum solvent using FACTS method   **/
/**                                                                                                **/
/****************************************************************************************************/
class FACTSResidueInfo : public utility::pointer::ReferenceCount {

public:
	typedef conformation::Residue Residue;

public:
	//The function clone creates an exact copy of the FACTSResidueInfo object.
	//It returns a pointer that points to the new object
	FACTSResidueInfoOP clone() const{
		return new FACTSResidueInfo( *this );
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

	inline SelfNeighInfo const &selfneigh( Size const atm ) const { return selfneigh_[atm]; }
	inline Real dG_dCi( Size const atm ) const { return dG_dCi_[atm]; }
	inline Real dSA_dDi( Size const atm ) const { return dSA_dDi_[atm]; }
	inline Real dsolv_dBR( Size const atm ) const { return dsolv_dBR_[atm]; }
	inline Vector elecF2( Size const atm ) const { return elecF2_[atm]; }
	inline Vector solvF2d( Size const atm ) const { return solvF2d_[atm]; }
	inline Vector solvF2BR( Size const atm ) const { return solvF2BR_[atm]; }
	inline Vector sasaF2( Size const atm ) const { return sasaF2_[atm]; }
	inline utility::vector1<Vector> const &xyz() const { return xyz_; }

	inline bool flag_for_calculation( Size const atm ) const { return flag_for_calculation_[atm]; }

public:
	const static Size MAXNEIGH = 300;
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

	utility::vector1<Vector> xyz_; // Store coordinate info
	bool changed_;
	bool enumeration_shell_;

	utility::vector1<SelfNeighInfo> selfneigh_;
	utility::vector1<Real> dG_dCi_;
	utility::vector1<Real> dSA_dDi_;
	utility::vector1<Real> dsolv_dBR_;
	utility::vector1<Vector> elecF2_;
	utility::vector1<Vector> solvF2d_;
	utility::vector1<Vector> solvF2BR_;
	utility::vector1<Vector> sasaF2_;

private:
	FACTSRsdTypeInfoCOP restypeinfo_;
};


/****************************************************************************************************/
/**                                                                                                **/
/**    @brief: The class FACTSPoseInfo                                                             **/
/**                                                                                                **/
/****************************************************************************************************/
class FACTSPoseInfo : public basic::datacache::CacheableData {

public:
	typedef conformation::Residue   Residue;
	typedef conformation::ResidueOP ResidueOP;

public:

	FACTSPoseInfo();
	FACTSPoseInfo( FACTSPoseInfo const & src );

	basic::datacache::CacheableDataOP	clone() const
	{	return new FACTSPoseInfo( *this );}

	Size size() const	{	return residue_info_.size(); }

	FACTSResidueInfo & residue_info( Size const i )
	{ return *residue_info_[i];	}

	FACTSResidueInfo const & residue_info( Size const i ) const
	{	return *residue_info_[i];	}

	bool being_packed( Size const seqpos ) const
	{	return being_packed_[ seqpos ];	}

	void set_placeholder( Size const i, ResidueOP rsd, FACTSResidueInfoOP info );

	FACTSResidueInfo const & placeholder_info( Size const seqpos ) const
	{
		assert( placeholder_info_[ seqpos ] );
		return *placeholder_info_[ seqpos ];
	}

	Residue const & placeholder_residue( Size const seqpos ) const
	{
		assert( placeholder_residue_[ seqpos ] );
		return *placeholder_residue_[ seqpos ];
	}

	void initialize( pose::Pose const & pose, FACTSRsdTypeMap &rsdtypemap );

	void set_repack_list( utility::vector1< bool > const & repacking_residues );

	bool is_changed( pose::Pose const &pose );

	void update_enumeration_shell( pose::Pose const &pose,
																 bool const enumerate_second_shell = false );

	inline bool context_derivative_empty() { return context_derivative_empty_; }

public:
	utility::vector1< FACTSResidueInfoOP > residue_info_; // these are allocated in initialize
	utility::vector1< ResidueOP > placeholder_residue_; // these may be null pointers
	utility::vector1< FACTSResidueInfoOP > placeholder_info_;
	utility::vector1< bool > being_packed_; // stores info from the packertask when setup_for_packing calls set_repack_list
	bool context_derivative_empty_;
};

/****************************************************************************************************/
/**                                                                                                **/
/**    @brief: The  class    FACTSRotamerSetInfo                                                   **/
/**                                                                                                **/
/****************************************************************************************************/
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

	basic::datacache::CacheableDataOP	clone() const
	{ return new FACTSRotamerSetInfo( *this ); }

	Size size() const
	{	return residue_info_.size(); }

	FACTSResidueInfo & residue_info( Size const i )
	{	return *residue_info_[i];	}

	FACTSResidueInfo const & residue_info( Size const i ) const
	{	return *residue_info_[i];	}

	FACTSRotamerSetInfo( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap )
	{	initialize( rotamer_set, rsdtypemap ); }

	/// dont forget to 0 the born_radii
	void initialize( RotamerSet const & rotamer_set, FACTSRsdTypeMap &rsdtypemap );

	private:
		utility::vector1< FACTSResidueInfoOP > residue_info_;
	};

/****************************************************************************************************/
/**                                                                                                **/
/**    @brief: The FACTSPotential class provides all the functions, constants, and parameters      **/
/**            common to all atoms required to calculate the free energy of solvation of a         **/
/** (macro)molecule embedded in a continuum solvent using FACTS method                             **/
/**                                                                                                **/
/****************************************************************************************************/
class FACTSPotential : public utility::pointer::ReferenceCount {

public:
	typedef conformation::Residue Residue;

public:
	FACTSPotential();

	void set_default();

	void setup_for_scoring(pose::Pose & pose, bool const & packing ) const;

	void setup_for_derivatives(pose::Pose & pose) const;

	void setup_for_packing(
												 pose::Pose & pose,
												 utility::vector1< bool > const & repacking_residues ) const;

	void update_residue_for_packing( pose::Pose & pose,
																	 Size const seqpos
																	 ) const;

	void get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const;

	void evaluate_polar_energy( Residue const & rsd1,
															FACTSResidueInfo const & facts1,
															Residue const & rsd2,
															Real & E_elec,
															Real & E_solv
															) const;

	Real evaluate_nonpolar_energy( Residue const & rsd1,
																 FACTSResidueInfo const & facts1,
																 Residue const & rsd2
																 ) const;

	void evaluate_context_change_for_packing(
						 Residue const & rsd1_ref,
						 Residue const & rsd1,
						 FACTSResidueInfo const & facts1,
						 Residue const & rsd2_ref,
						 Residue const & rsd2,
						 FACTSResidueInfo const & facts2,
						 utility::vector1< Real > & dBRi1,
						 utility::vector1< Real > & dBRi2,
						 utility::vector1< Real > & dSAi1,
						 utility::vector1< Real > & dSAi2
						 ) const;

	void evaluate_polar_otf_energy(Residue const & rsd1,
																 FACTSResidueInfo const & facts1,
																 Residue const & rsd2,
																 FACTSResidueInfo const & facts2,
																 utility::vector1< Real > const & dBRi1,
																 utility::vector1< Real > const & dBRi2,
																 Real & E_elec,
																 Real & E_solv,
																 bool do_correction
																 ) const;

	void eval_atom_polar_derivative(
					id::AtomID const & id,
					Real const weight_elec,
					Real const weight_solv,
					pose::Pose const & pose,
					kinematics::DomainMap const & domain_map,
					bool const exclude_DNA_DNA,
					Vector & F1,
					Vector & F2
					) const;

	void eval_atom_nonpolar_derivative(
					id::AtomID const & id,
					Real const weight,
					pose::Pose const & pose,
					kinematics::DomainMap const & domain_map,
					bool const exclude_DNA_DNA,
					Vector & F1,
					Vector & F2
					) const;

	void get_single_rotamer_born_radii(
																	Residue const & rsd1,
																	pose::Pose const & pose,
																	FACTSPoseInfo const & facts_info,
																	FACTSResidueInfo & facts1
																	) const;

	Real polar_energy_pack_corrector(
																	 Residue const & ref_rsd,
																	 Residue const & rsd,
																	 FACTSResidueInfo const & facts_info
																	 ) const;

private:
	void res_res_burial(
											Residue const & rsd1,
											FACTSResidueInfo & facts1,
											Residue const & rsd2,
											FACTSResidueInfo const & facts2
											) const;

	void res_res_burial_for_scoring(
												Residue const & rsd1,
												FACTSResidueInfo & facts1,
												Residue const & rsd2,
												FACTSResidueInfo & facts2
												) const;

	void get_self_terms( FACTSRsdTypeInfoCOP factstype1,
											 FACTSResidueInfo & facts1,
											 bool const packing
											 ) const;

	void calculate_GBpair_apprx(
													Residue const & rsd1,
													Residue const & rsd2,
													FACTSResidueInfo & facts1,
													FACTSResidueInfo & facts2
													) const;

	void calculate_GBpair_exact(
													Residue const & rsd1,
													Residue const & rsd2,
													FACTSResidueInfo & facts1,
													FACTSResidueInfo & facts2
													) const;


	void atom_atom_context_derivative( FACTSResidueInfo & facts1,
																		 FACTSResidueInfo & facts2,
																		 Size const & atm1,
																		 Size const & atm2,
																		 Vector const & dxyz,
																		 Real const & dB_dBdnmtr,
																		 Real const & dB_dBnmtr,
																		 Real const & dBR_dG,
																		 Real const & CutOff_sqr1,
																		 Real const & CutOff_sqr,
																		 bool const full_update
																		 ) const;

	void get_template_born_radii(
														pose::Pose const & pose,
														FACTSPoseInfo & facts_info
														) const;

	// Accessors
	inline Real Tau() const{	return Tau_; }
	inline Real inv_die() const{	return inv_die_; }
	inline Real Kappa() const { return Kappa_; }
	inline Real MultiplicitiveFactor() const { return MultiplicitiveFactor_; };
	inline Real GBPair_cut() const { return GBpair_cut_; };

private: //list of private variables and parameters for the FACTS method common to all atoms

	// Map storing parameters for residue types
	mutable FACTSRsdTypeMap FACTSrsdtypemap_;

	bool options_registered_;

	Real MultiplicitiveFactor_;
	Real inv_die_;
	Real Tau_;
	Real Kappa_;
	Real GBpair_cut_;
	bool do_apprx;
	Real selfenergy_scale_;
	Real intrares_scale_;
	Real min_dis_;

	Real saltbridge_correction_;
	Real dshift2_;

	Real elec_sh_exponent_;

	// Below are not being used
	Real cut_off_born_radius_; //The cut off used for calculating born radius
	Real extra_cut_off_self_;
	Real extra_cut_off_interaction_;
	Real dummy_radius_;
	Real dummy_scale_;
	Real dummy_distance_; // also implicitly defined by the gb placeholder params file

};

} // scoring
} // core

#endif
