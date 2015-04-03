// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/chemical_shift/RNA_ChemicalShiftPotential.hh
/// @brief  The real workhorse behind RNA_ChemicalShiftEnergy.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftPotential_HH
#define INCLUDED_core_scoring_rna_chemical_shift_RNA_ChemicalShiftPotential_HH


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>

#include <utility/pointer/ReferenceCount.hh>
///////////////////////////////////
#include <core/chemical/AA.hh>
#include <core/scoring/rna/chemical_shift/RNA_CS_Parameters.hh>

#include <core/scoring/EnergyMap.fwd.hh>

// Utility headers


namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {


/////////////////////////////////////////////////////////////////////////////
class ChemicalShiftData : public utility::pointer::ReferenceCount{

	public:

		ChemicalShiftData( core::Size const in_seq_num, 
										 chemical::AA const in_res_aa,
										 std::string const in_atom_name, 
										 core::Size const in_realatomdata_index,
										 core::Real const in_exp_shift,
										 core::Real const in_ref_shift, 
										 std::string const in_data_line,
										 core::Real const in_accuracy_weight):
			seq_num( in_seq_num ),
			res_aa( in_res_aa ),
			atom_name( in_atom_name ),
			realatomdata_index( in_realatomdata_index ),
			exp_shift( in_exp_shift ),
			ref_shift( in_ref_shift ), 
			//Feb 28, 2012: Warning exp_shift this might be switch due to ambiguity of geminal atoms (H5'/H5'' and etcs) or be a duplicate.
			//Always use get_best_exp_to_calc_chem_shift_mapping() to get best mapping to specific calc_chem_shift to get actual exp_shift!
			data_line( in_data_line ),
			accuracy_weight(  in_accuracy_weight )
		{
		}

		~ChemicalShiftData(){};

	public:

		core::Size seq_num;
		chemical::AA res_aa;
		std::string atom_name;
		core::Size realatomdata_index; //For accessing 	realatomdata in RNA_CS_Parameters	
		core::Real exp_shift;	
		core::Real ref_shift;		
		std::string data_line;
		core::Real accuracy_weight;
};


class RNA_ChemicalShiftPotential : public utility::pointer::ReferenceCount{

	public:

		RNA_ChemicalShiftPotential();

		Size
		get_total_exp_chemical_shift_data_points() const;
		

	private:

		Size
		get_realatomdata_index( std::string const & in_atom_name, chemical::AA const res_aa ) const;

		void
		assert_is_calc_chem_shift_atom( ChemicalShiftData const & CS_data ) const;

		bool
		Is_magnetic_anisotropy_source_atom( core::conformation::Residue const & rsd, Size const atomno ) const;

		bool
		atom_has_exp_chemical_shift_data( core::conformation::Residue const & rsd, Size const atomno ) const;

		utility::vector1 < ChemicalShiftData > const & 
		get_matching_CS_data_entry( Size const seq_num, std::string const in_atom_name ) const;

		utility::vector1< std::string > 
		string_list( std::string const string_one ) const;

		utility::vector1< std::string > 
		string_list( std::string const string_one, std::string const string_two ) const;
		

		void
		import_exp_chemical_shift_data( std::string exp_CS_data_filename, 
																 utility::vector1 < core::Size > include_res_list,
																 utility::vector1 < utility::vector1< std::string > > const & proton_entry_list );

		void
		get_best_exp_to_calc_chem_shift_mapping( utility::vector1 < ChemicalShiftData > const & EXP_chem_shift_data_entry, 
																					utility::vector1 < Real > const & calc_chem_shift_entry, 
														 							utility::vector1 < Real > & actual_exp_chem_shift_entry,
														 							utility::vector1 < bool > & do_include_CS_data ) const; 

		core::Real
		get_calc_chem_shift_value_nuchemics( ChemicalShiftData const & CS_data, pose::Pose const & pose) const;

		core::Real
		get_calc_chem_shift_value_larmord( ChemicalShiftData const & CS_data, pose::Pose const & pose) const;

		core::Real
		get_calc_chem_shift_value( ChemicalShiftData const & CS_data, pose::Pose const & pose) const;

		void
		update_calc_chem_shift_list( pose::Pose const & pose, utility::vector1 < utility::vector1 < Real > > & cal_chem_shift_list ) const;

		core::Real
		get_chemical_shift_energy( utility::vector1 < utility::vector1 < Real > > const & calc_chem_shift_list ) const;


	void
	get_deriv_for_chemical_shift( id::AtomID const & atom_id,
																ChemicalShiftData const & CS_data,
																pose::Pose const & pose,
																Vector & f1,
																Vector & f2,
																utility::vector1 < std::string > atom_names,
																std::string atom_name_in,
																std::string atom_name_whitespace_in,
																bool is_source_atom,
																bool is_neighbor_atom) const;

		void
		get_deriv_for_chemical_shift_data_atom(	pose::Pose const & pose,
																					conformation::Residue const & CS_data_rsd,
																					Size const CS_data_atomno,
																					Vector & f1,
																					Vector & f2 ) const;

		void
		get_ring_current_deriv_for_src_base(	pose::Pose const & pose,
																				conformation::Residue const & rc_source_rsd,
																				Size const chi1_torsion_atomnno,
																				Vector & f1,
																				Vector & f2 ) const;

		void
		get_magnetic_anisotropy_deriv_for_src_base( pose::Pose const & pose,
																							conformation::Residue const & ma_source_rsd,
																							Size const chi1_torsion_atomnno,
																							Vector & f1,
																							Vector & f2 ) const;

		void
		get_magnetic_anisotropy_deriv_for_src_atom( pose::Pose const & pose,
																							conformation::Residue const & ma_source_rsd,
																							Size const ma_source_atomno,
																							Vector & f1,
																							Vector & f2 ) const; 


	public:

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

		void
		finalize_total_energy( pose::Pose const & pose, EnergyMap & totals ) const;
				
		void 
		load_larmord_parameters( std::string  const filename );
		
		void 
		load_larmord_weights ( std::string  const filename );

		void 
		load_larmord_reference_shifts ( std::string  const filename );

		void 
		load_larmord_neighbor_atoms ( std::string  const filename );

    bool 
    get_neighbor_atom( std::string const & key ) const;

    Real 
    get_accuracy_weight( std::string const & key ) const;

    Real 
    get_reference_shift( std::string const & key ) const;
    
    Real 
    get_alpha ( std::string const & key ) const;		
    
		/////////////////////////////////
		void
		eval_atom_derivative(
			id::AtomID const & atom_id,
			pose::Pose const & pose,
			kinematics::DomainMap const & domain_map,
			EnergyMap const & weights,
			Vector & F1,
			Vector & F2 ) const;


		void
		indicate_required_context_graphs(
			utility::vector1< bool > & /*context_graphs_required*/
		) const {}

	private:

		RNA_CS_parameters const rna_cs_params_;
		std::string H5_prime_mode_;
		bool const verbose_;
		bool const include_ring_current_effect_;
		bool const include_magnetic_anisotropy_effect_;
		bool nuchemics_mode_;
		bool cs_verbose_mode_;
		utility::vector1 < utility::vector1 < ChemicalShiftData > > EXP_chem_shift_data_list_;
		core::Size total_exp_chemical_shift_data_points_;
		std::map< std::string, Real > reference_shifts_;
		std::map< std::string, Real > alphas_;	
		std::map< std::string, Real > accuracy_weights_;
		std::map< std::string, bool > neighbor_atoms_;
		std::string path_to_parameter_files_;
		Real larmord_distance_cutoff_;
		Real larmord_beta_;
};


} //chemical_shift
} //rna
} //scoring
} //core

#endif // INCLUDED_core_scoring_methods_RG_Energy_RNA_HH
