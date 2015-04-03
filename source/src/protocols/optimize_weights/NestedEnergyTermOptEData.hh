// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/optimize_weights/NestedEnergyTermOptEData.hh
/// @brief Classes to store by-energy-term rotamer data for optE, when using the unfolded state energy term
/// @author Ron Jacak

#ifndef INCLUDED_protocols_optimize_weights_NestedEnergyTermOptEData_hh
#define INCLUDED_protocols_optimize_weights_NestedEnergyTermOptEData_hh

// Unit headers
#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/NestedEnergyTermOptEData.fwd.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/optimization/types.hh>


#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {

class NestedEnergyTermPNatAAOptEPositionData : public protocols::optimize_weights::PNatAAOptEPositionData {

public:
	//typedef core::chemical::AA AA;

	NestedEnergyTermPNatAAOptEPositionData();
	virtual ~NestedEnergyTermPNatAAOptEPositionData();

	virtual Real get_score(
		Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	virtual void print_score(
		std::ostream & ostr, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	Real process_score(
		std::ostream & ostr, bool print, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	virtual OptEPositionDataType type() const;

	virtual void write_to_file( std::ofstream & outfile ) const ;
	virtual void read_from_file( std::ifstream & infile );
	virtual void write_to_binary_file( std::ofstream & outfile ) const;
	virtual void read_from_binary_file( std::ifstream & infile );

	#ifdef USEMPI
		virtual void send_to_node( int const destination_node, int const tag ) const;
		virtual void receive_from_node( int const source_node, int const tag );
	#endif

	void set_unfolded_energy_emap_vector( utility::vector1< EnergyMap > emap_vector ) { unfolded_energy_emap_vector_ = emap_vector; }
	utility::vector1 < EnergyMap > & unfolded_energy_emap_vector() { return unfolded_energy_emap_vector_; }

	// Use the base class implementation for these methods

	//virtual void range( ScoreTypes const & free_score_list, ScoreTypes const & fixed_score_list, EnergyMap & lower_bound, EnergyMap & upper_bound ) const;
	//virtual Size size() const { return data_.size(); }

	//void set_position( Size pos_in ) { position_ = pos_in; }
	//Size position() const { return position_; }

	//void set_native_aa( AA nat_in ) { native_aa_ = nat_in; }
	//AA native_aa() const { return native_aa_; }

	//void set_neighbor_count( Size nb_in ) { neighbor_count_ = nb_in; }
	//Size neighbor_count() const { return neighbor_count_; }

	//void add_rotamer_line_data( PNatAAOptERotamerDataOP rot_in ) { data_.push_back( rot_in ); }

	//PNatAAOptERotamerDataOPs & data() { return data_; }
	//PNatAAOptERotamerDataOPs const & data() const { return data_; }

	//PNatAAOptERotamerDataOPs::const_iterator rotamer_data_begin() const { return data_.begin(); }
	//PNatAAOptERotamerDataOPs::const_iterator rotamer_data_end() const { return data_.end(); }

	virtual Size memory_use() const;

protected:
	/// @brief used by derived class as well -- finds the energies for the best rotamer for each amino acid
	//void
	//process_rotamers( Multivec const & vars, Size const num_energy_dofs, EnergyMap const & fixed_terms,
	//	ScoreTypes const & score_list, ScoreTypes const & fixed_score_list, Size const aa_range,
	//	utility::vector1< Real > const & dummy_set, utility::vector1< Real > & best_energy_by_aa,
	//	utility::vector1< utility::vector1< Real > > & unweighted_E_dof, Multivec & ref_deriv_weight
	//) const;

private:
	//Size position_;
	//AA native_aa_;
	//Size neighbor_count_;
	//PNatAAOptERotamerDataOPs data_;

	utility::vector1< EnergyMap > unfolded_energy_emap_vector_;

};


class NestedEnergyTermDDGMutationOptEData : public DDGMutationOptEData {

public:
	//typedef core::chemical::AA AA;

public:
	NestedEnergyTermDDGMutationOptEData();
	virtual ~NestedEnergyTermDDGMutationOptEData();

	virtual Real get_score(
		Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	virtual void print_score(
		std::ostream & ostr, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	Real process_score(
		std::ostream & ostr, bool print, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	virtual OptEPositionDataType type() const;

	// Use the base class implementation for these methods. (The read/write methods are all empty anyway.)

	//virtual void range( ScoreTypes const & free_score_list, ScoreTypes const & fixed_score_list, EnergyMap & lower_bound, EnergyMap & upper_bound ) const;
	//virtual Size size() const;

	//virtual void write_to_file( std::ofstream & outfile ) const;
	//virtual void read_from_file( std::ifstream & infile );
	//virtual void write_to_binary_file( std::ofstream & outfile ) const;
	//virtual void read_from_binary_file( std::ifstream & infile );

	virtual Size memory_use() const;

#ifdef USEMPI
	virtual void send_to_node( int const destination_node, int const tag ) const;
	virtual void receive_from_node( int const source_node, int const tag );
#endif

	//void set_wt_aa( AA wt_aa );
	//void set_mut_aa( AA mut_aa );

	//void set_experimental_ddg( Real ddg );
	//void add_wt( SingleStructureDataOP wt );

	//void add_mutant( SingleStructureDataOP mut );

	void set_wt_unfolded_energies_emap( EnergyMap e ) { wt_unfolded_energies_emap_ = e; }
	void set_mut_unfolded_energies_emap( EnergyMap e ) { mut_unfolded_energies_emap_ = e; }

private:
	//Real experimental_ddG_;
	//AA wt_aa_;
	//AA mut_aa_;
	//SingleStructureDataOPs wts_;
	//SingleStructureDataOPs muts_;

	EnergyMap wt_unfolded_energies_emap_;
	EnergyMap mut_unfolded_energies_emap_;

};


} // namespace optimize_weights
} // namespace protocols

#endif
