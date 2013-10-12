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
// @author: Massih Khorvash (massih.khorvash@gmail.com)
// @author: Hahnbeom Park

#ifndef INCLUDED_core_scoring_facts_FACTSPotential_HH
#define INCLUDED_core_scoring_facts_FACTSPotential_HH

// Unit headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
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

	// Structure for storing Self-Neighbor(not GB-Neighbor) info
	struct SelfNeighInfo {
		Size nneigh;
		utility::vector1<Size> resID;
		utility::vector1<Size> atmID;
	};

	// Forward declaration
	class FACTSRsdTypeInfo;
	typedef utility::pointer::owning_ptr< FACTSRsdTypeInfo > FACTSRsdTypeInfoOP;
	typedef utility::pointer::owning_ptr< FACTSRsdTypeInfo const > FACTSRsdTypeInfoCOP;

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
		Size natoms() const{ return natoms_; }
		Real q( Size const atm ) const{ return q_[ atm ]; }
		Real a0(Size const atm) const{ return a0_[atm]; }
		Real a1(Size const atm) const{ return a1_[atm]; }
		Real a2(Size const atm) const{ return a2_[atm]; }
		Real a3(Size const atm) const{ return a3_[atm]; }
		Real b1(Size const atm) const{ return b1_[atm]; }
		Real b2(Size const atm) const{ return b2_[atm]; }
		Real c0(Size const atm) const{ return c0_[atm]; }
		Real c1(Size const atm) const{ return c1_[atm]; }
		Real c2(Size const atm) const{ return c2_[atm]; }
		Real c3(Size const atm) const{ return c3_[atm]; }
		Real d1(Size const atm) const{ return d1_[atm]; }
		Real d2(Size const atm) const{ return d2_[atm]; }

		Real alpha( Size const atm ) const { return alpha_[atm]; }
		Real COradius2(Size const atm ) const{ return COradius2_[atm]; }
		Real volume(Size const atm ) const{ return volume_[atm]; }
		bool not_using( Size const atm ) const {return not_using_[atm]; }

		bool selfpair( Size const atm1, Size atm2 ) const { return selfpair_[atm1][atm2]; }
		bool is_chargedH( Size const atm1 ) const { return is_chargedH_[atm1]; }

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
			return new FACTSResidueInfo( *this );
		}

		//The constructor with a residue as its input parameter.
		FACTSResidueInfo( Residue const & rsd, bool const is_rotamer = false ){
			initialize( rsd, is_rotamer );
		}

		//This function initializes all the private and public variables of the class FACTSResidueInfo
		void initialize( Residue const & rsd, bool const is_rotamer = false );

		void allocate_selfneigh();

		void store_xyz( Residue const & rsd );

		Size natoms() const { return natoms_; }
		Real Ai( Size const atm ) const { return Ai_[atm]; }
		Real Bi( Size const atm ) const { return Bi_[atm]; }
		Real Ci( Size const atm ) const { return Ci_[atm]; }
		Real Di( Size const atm ) const { return Di_[atm]; }
		Real Ei( Size const atm ) const { return Ei_[atm]; }

		utility::vector1< Real > Ai() const { return Ai_; }
		utility::vector1< Vector > nmtr() const { return nmtr_; }
		utility::vector1< Real > dnmtr() const { return dnmtr_; }

		Real esolvE( Size const atm ) const { return esolvE_[atm]; }
		Vector nmtr( Size const atm ) const { return nmtr_[atm]; }
		Real dnmtr( Size const atm ) const { return dnmtr_[atm]; }
		Real sasa( Size const atm ) const { return sasa_[atm]; }
		Real BR( Size const atm ) const { return BR_[atm]; }
		Real GBpair( Size const res ) const { return GBpair_[res]; }

		SelfNeighInfo selfneigh( Size const atm ) const { return selfneigh_[atm]; }
		Real dG_dCi( Size const atm ) const { return dG_dCi_[atm]; }
		Real dSA_dDi( Size const atm ) const { return dSA_dDi_[atm]; }
		Real dE_dBR( Size const atm ) const { return dE_dBR_[atm]; }
		Vector dE_drij2( Size const atm ) const { return dE_drij2_[atm]; }
		Vector polarF2d( Size const atm ) const { return polarF2d_[atm]; }
		Vector polarF2BR( Size const atm ) const { return polarF2BR_[atm]; }
		Vector nonpolarF2( Size const atm ) const { return nonpolarF2_[atm]; }
		utility::vector1<Vector> xyz() const { return xyz_; }

		bool flag_for_calculation( Size const atm ) const { return flag_for_calculation_[atm]; }

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
		utility::vector1<Real> GBpair_;

		utility::vector1<Vector> xyz_; // Store coordinate info

		utility::vector1<SelfNeighInfo> selfneigh_;
		utility::vector1<Real> dG_dCi_;
		utility::vector1<Real> dSA_dDi_;
		utility::vector1<Real> dE_dBR_;
		utility::vector1<Vector> dE_drij2_;
		utility::vector1<Vector> polarF2d_;
		utility::vector1<Vector> polarF2BR_;
		utility::vector1<Vector> nonpolarF2_;

	};

	////////////////////////////////////////////////////////////////////////////////////////////////////
	class FACTSPoseInfo : public basic::datacache::CacheableData {

	public:
		typedef conformation::Residue   Residue;
		typedef conformation::ResidueOP ResidueOP;

	public:

		FACTSPoseInfo() {};
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

		void initialize( pose::Pose const & pose);

		void set_repack_list( utility::vector1< bool > const & repacking_residues );

		bool is_changed( pose::Pose const &pose );

	public:
		utility::vector1< FACTSResidueInfoOP > residue_info_; // these are allocated in initialize
		utility::vector1< ResidueOP > placeholder_residue_; // these may be null pointers
		utility::vector1< FACTSResidueInfoOP > placeholder_info_;
		utility::vector1< bool > being_packed_; // stores info from the packertask when setup_for_packing calls set_repack_list
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////
	class FACTSRotamerSetInfo : public basic::datacache::CacheableData {

	public:
		typedef conformation::Residue   Residue;
		typedef conformation::ResidueOP ResidueOP;
		typedef conformation::RotamerSetBase RotamerSet;

	public:
		///
		FACTSRotamerSetInfo( FACTSRotamerSetInfo const & src ):
			CacheableData()
		{
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

		FACTSRotamerSetInfo( RotamerSet const & rotamer_set )
		{	initialize( rotamer_set ); }

		/// dont forget to 0 the born_radii
		void initialize( RotamerSet const & rotamer_set );

		private:
			utility::vector1< FACTSResidueInfoOP > residue_info_;
		};

	/**************************************************************************************************/
	/*                                                                                                */
	/*    @brief: The FACTSPotential class provides all the functions, constants, and parameters      */
	/*            common to all atoms required to calculate the free energy of solvation of a         */
	/* (macro)molecule embedded in a continuum solvent using FACTS method               */
	/*                                                                                                */
	/**************************************************************************************************/

	class FACTSPotential {

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

		void setup_for_packing(pose::Pose & pose,	pack::task::PackerTask const & task	) const;

		void update_residue_for_packing( pose::Pose & pose,
																		 Size const seqpos
																		 ) const;

		void create_rsdtypeinfo(pose::Pose & pose) const;

		void get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const;

		Real evaluate_polar_energy( Residue const & rsd1,
																FACTSResidueInfo const & facts1,
																Residue const & rsd2
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

		Real evaluate_polar_otf_energy(Residue const & rsd1,
																		 FACTSResidueInfo const & facts1,
																		 Residue const & rsd2,
																		 FACTSResidueInfo const & facts2,
																		 utility::vector1< Real > const & dBRi1,
																		 utility::vector1< Real > const & dBRi2,
																		 bool do_correction
																		 ) const;

		void eval_atom_polar_derivative(
						id::AtomID const & id,
						Real const weight,
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
																		FACTSPoseInfo const & gb_info,
																		FACTSResidueInfo & gb1
																		) const;

		Real polar_energy_pack_corrector( 
																		 Residue const & ref_rsd,
																		 Residue const & rsd,
																		 FACTSResidueInfo const & facts_info
																		 ) const;

	private:
		void res_res_burial(
												Residue const & rsd1,
												FACTSResidueInfo & gb1,
												Residue const & rsd2,
												FACTSResidueInfo const & gb2
												) const;

		void res_res_burial_for_scoring(
													Residue const & rsd1,
													FACTSResidueInfo & gb1,
													Residue const & rsd2,
													FACTSResidueInfo & gb2
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

		void get_template_born_radii(
															pose::Pose const & pose,
															FACTSPoseInfo & gb_info
															) const;

		void	build_placeholders(
												 pose::Pose const & pose,
												 FACTSPoseInfo & facts_info
												 ) const;

		// Accessors
		Real Tau() const{	return Tau_; }
		Real Kappa() const { return Kappa_; }
		Real MultiplicitiveFactor() const { return MultiplicitiveFactor_; };
		Real GBPair_cut() const { return GBpair_cut_; };

		FACTSRsdTypeInfoCOP get_factstype( core::chemical::ResidueType const &rsdtype ) const;

	private: //list of private variables and parameters for the FACTS method common to all atoms

		// Map storing parameters for residue types
		mutable std::map< chemical::ResidueType const *, FACTSRsdTypeInfoCOP > FACTSrsdtypemap_;
		//mutable std::map< std::string const, FACTSRsdTypeInfoOP > FACTSrsdtypemap_;

		//This variable contains the solvation free energy of the macromolecule
		//which is the sum of FACTSResidueInfo::sasa and FACTSResidueInfo::esolvE_i
		//solvation_free_energy corresponds to deltaG(FACTS) = deltaG(el,FACTS) + Gamma * sigma(Si(FACTS)) i.e. equation 12 on page 706 of FACTS paper

		bool options_registered_;

		Real MultiplicitiveFactor_;
		Real Tau_;
		Real Kappa_;
		Real GBpair_cut_;
		bool do_apprx;
		Real selfenergy_scale_;
		Real intrares_scale_;
		Real min_dis_;

		Real saltbridge_correction_;
		Real dshift2_;

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
