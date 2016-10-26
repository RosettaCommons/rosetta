// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/optimize_weights/DDGBindOptEData.hh
/// @brief Header file for OptEPositionData class that hold interface ddG information
/// @author Ron Jacak

#ifndef INCLUDED_protocols_optimize_weights_DDGBindOptEData_hh
#define INCLUDED_protocols_optimize_weights_DDGBindOptEData_hh

// Unit headers
#include <protocols/optimize_weights/DDGBindOptEData.fwd.hh>
#include <protocols/optimize_weights/OptEData.hh>

#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {

class DDGBindOptEData : public protocols::optimize_weights::OptEPositionData {

public:

	typedef core::chemical::AA AA;
	enum DDG_Bind_File_Index { WT_COMPLEXES_LIST_FILE = 1, MUT_COMPLEXES_LIST_FILE, WT_UNBOUNDS_LIST_FILE, MUT_UNBOUNDS_LIST_FILE };

	DDGBindOptEData();
	~DDGBindOptEData() override;

	Real get_score(
		Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const override;

	void print_score(
		std::ostream & ostr, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const override;

	Real process_score(
		std::ostream & ostr, bool print, Multivec const & component_weights,
		Multivec const & vars, Multivec & dE_dvars, Size const num_energy_dofs, int const num_ref_dofs, int const num_total_dofs,
		EnergyMap const & fixed_terms, ScoreTypes const & score_list, ScoreTypes const & fixed_score_list ) const;

	OptEPositionDataType type() const override;
	void range( ScoreTypes const & free_score_list, ScoreTypes const & fixed_score_list, EnergyMap & lower_bound, EnergyMap & upper_bound ) const override;
	Size size() const override;
	Size memory_use() const override;

#ifdef USEMPI
	virtual void send_to_node( int const destination_node, int const tag ) const override;
	virtual void receive_from_node( int const source_node, int const tag ) override;
#endif

	void write_to_file( std::ofstream & /* outfile */ ) const override {}
	void read_from_file( std::ifstream & /* infile */ ) override {}
	void write_to_binary_file( std::ofstream & /* outfile */ ) const override {}
	void read_from_binary_file( std::ifstream & /* infile */ ) override {}

	// setters
	void set_experimental_ddg_bind( Real exp_ddg_bind );
	void add_mutation( std::pair< Size, std::pair < AA, AA > > mutation );

	void add_wt_complex( SingleStructureDataOP wt );
	void add_mutant_complex( SingleStructureDataOP mut );
	void add_wt_unbounds( SingleStructureDataOP wt );
	void add_mutant_unbounds( SingleStructureDataOP mut );

private:
	Real experimental_ddG_bind_;
	utility::vector1< std::pair< Size, std::pair < AA, AA > > > mutations_;

	SingleStructureDataOPs wt_complexes_;
	SingleStructureDataOPs mutant_complexes_;
	SingleStructureDataOPs wt_unbounds_;
	SingleStructureDataOPs mutant_unbounds_;

};


} // namespace optimize_weights
} // namespace protocols

#endif // INCLUDED_protocols_optimize_weights_DDGBindOptEData_HH
