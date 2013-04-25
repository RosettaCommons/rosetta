// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran.hh
/// @brief  Ramachandran potential class delcaration
/// @author Guoli Wang

#ifndef INCLUDED_core_scoring_Ramachandran2B_hh
#define INCLUDED_core_scoring_Ramachandran2B_hh

// Unit Headers
#include <core/scoring/Ramachandran2B.fwd.hh>
//#include <core/scoring/ProteinTorsion.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>

// C++ Headers
#include <map>


namespace core {
namespace scoring {


class Ramachandran2B : public utility::pointer::ReferenceCount
{
public:
	typedef pose::Pose Pose;
	typedef chemical::AA AA;

public:
	Ramachandran2B();
	virtual ~Ramachandran2B() ; // auto-removing definition from header{}

	Real
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi
	) const;

	// Guoli Wang
	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		chemical::AA const left_aa,
		chemical::AA const right_aa,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Lower(
		conformation::Residue const &rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE_Upper(
		conformation::Residue const & rsd,
		chemical::AA const &neighbor
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd,
		Real &drama_dphi,
		Real &drama_dpsi
	) const;

	Real
	RamaE(
		conformation::Residue const & rsd
	) const;
	// finished

	void
	eval_rama_score_residue(
		conformation::Residue const & res,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	eval_rama_score_residue(
		AA const res_aa,
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi
	) const;

	void
	IdealizeRamaEnergy(
		Real const phi,
		Real const psi,
		Real & rama,
		Real & drama_dphi,
		Real & drama_dpsi,
		Real const entropy,
		ObjexxFCL::
		FArray2A< Real > const & rama_for_res
	) const;

	// AS: for neighbor-dependent rama sampling
    void
	random_phipsi_from_rama_left(
                                 AA const left_aa,
                                 AA const pos_aa,
                                 Real & phi,
                                 Real & psi
                                 ) const;
	
    void
	random_phipsi_from_rama_right(
								  AA const pos_aa,
								  AA const right_aa,
								  Real & phi,
								  Real & psi
								  ) const;
	
	/// @brief function for torsion-bin specific but otherwise random phi/psi angles
	/// @author Amelie Stein
    
	void
	random_phipsi_from_rama_by_torsion_bin_left(
												AA const left_aa,
												AA const pos_aa,
												Real & phi,
												Real & psi,
												char const torsion_bin
												) const;
	
	void
	random_phipsi_from_rama_by_torsion_bin_right(
												 AA const pos_aa,
												 AA const right_aa,
												 Real & phi,
												 Real & psi,
												 char const torsion_bin
												 ) const;
	
	
	core::Size get_torsion_bin_index(char torsion_bin) const;
	
	void
	init_rama_sampling_tables_by_torsion_bin(); 
	
	void
	get_entries_per_torsion_bin_left( 
									 AA const left_aa,
									 AA const pos_aa,
									 std::map< char, core::Size > & tb_frequencies ) const;
    
	void
	get_entries_per_torsion_bin_right( 
									  AA const pos_aa,
									  AA const right_aa,
									  std::map< char, core::Size > & tb_frequencies ) const;
    
	

	
	///////////////////////////////
	// unused??
	void
	eval_rama_score_all(
		Pose & pose,
		ScoreFunction const & scorefxn
	) const;

	void
	write_rama_score_all(
		Pose const & pose
	) const;


private:

	void read_rama();
    void init_rama_sampling_table_left( const char torsion_bin ); // AS: for neighbor-dependent rama sampling in KIC/NGK
    void init_rama_sampling_table_right( const char torsion_bin ); 

	//static bool rama_initialized_;
	ObjexxFCL::FArray3D< Real > ram_energ_;
	ObjexxFCL::FArray1D< Real > ram_entropy_;
	ObjexxFCL::FArray4D< Real > ram_energ_left_;
	ObjexxFCL::FArray2D< Real > ram_entropy_left_;
	ObjexxFCL::FArray4D< Real > ram_energ_right_;
	ObjexxFCL::FArray2D< Real > ram_entropy_right_;

	static int const n_phi_ = 36;
	static int const n_psi_ = 36;
	static Real const binw_; // 360 / n_phi_ = 10;
	static int const n_aa_ = 20; // Ramachandran score defined for the cananical AAs only.
	static int const nullaa = 21; // Guoli Wang

	Real const rama_score_limit_;

    static ObjexxFCL::FArray4D< Real > left_ram_probabil_; // AS: probability of each phi/psi combination given the left & actual  AA 
    utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > > left_rama_sampling_table_; // vector of allowed phi/psi pairs for each residue triplet
    utility::vector1< utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > > > left_rama_sampling_table_by_torsion_bin_; 
	
    static ObjexxFCL::FArray4D< Real > right_ram_probabil_; // AS: probability of each phi/psi combination given the actual & right  AA 
    utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > > right_rama_sampling_table_; // vector of allowed phi/psi pairs for each residue triplet
    utility::vector1< utility::vector1< utility::vector1< utility::vector1< utility::vector1< Real > > > > > right_rama_sampling_table_by_torsion_bin_; 
	
	static Real const rama_sampling_thold_;
	static Real const rama_sampling_factor_;
	

};

}
}

#endif
