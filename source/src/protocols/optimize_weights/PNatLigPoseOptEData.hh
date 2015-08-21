// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/optimize_weights/PNatLigPoseOptEData.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_optimize_weights_PNatLigPoseOptEData_hh
#define INCLUDED_protocols_optimize_weights_PNatLigPoseOptEData_hh

#include <protocols/optimize_weights/PNatLigPoseOptEData.fwd.hh>
#include <protocols/optimize_weights/OptEData.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {


/// @brief
///
/// @details
///
class PNatLigPoseOptEData : public PNatStructureOptEData
{
public:

	PNatLigPoseOptEData();
	virtual ~PNatLigPoseOptEData();

	// my invention to avoid code duplication
	virtual
	Real
	do_score(
		std::ostream & ostr,
		Multivec const & component_weights,
		Multivec const & vars,
		Multivec & dE_dvars,
		/// Basically, turn over all the private data from OptEMultiFunc
		Size const num_energy_dofs,
		int const num_ref_dofs,
		int const num_total_dofs,
		EnergyMap const & fixed_terms,
		ScoreTypes const & score_list,
		ScoreTypes const & fixed_score_list,
		bool const print
	) const;

	// inherited from OptEPositionData
	virtual
	Real
	get_score(
		Multivec const & component_weights,
		Multivec const & vars,
		Multivec & dE_dvars,
		/// Basically, turn over all the private data from OptEMultiFunc
		Size const num_energy_dofs,
		int const num_ref_dofs,
		int const num_total_dofs,
		EnergyMap const & fixed_terms,
		ScoreTypes const & score_list,
		ScoreTypes const & fixed_score_list
	) const
	{ return do_score(std::cout, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list, false /* don't print */); }

	virtual
	void
	print_score(
		std::ostream & ostr,
		Multivec const & component_weights,
		Multivec const & vars,
		Multivec & dE_dvars,
		/// Basically, turn over all the private data from OptEMultiFunc
		Size const num_energy_dofs,
		int const num_ref_dofs,
		int const num_total_dofs,
		EnergyMap const & fixed_terms,
		ScoreTypes const & score_list,
		ScoreTypes const & fixed_score_list
	) const
	{ do_score(ostr, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, score_list, fixed_score_list, true /* do print */); }

	/*
	virtual
	Size
	size() const;
	*/

	virtual
	OptEPositionDataType
	type() const;

	/*
	virtual
	void
	write_to_file( std::ofstream & outfile ) const;

	virtual
	void
	read_from_file( std::ifstream & infile );

	virtual
	void
	write_to_binary_file( std::ofstream & outfile ) const;

	virtual
	void
	read_from_binary_file( std::ifstream & infile );

	virtual
	Size
	memory_use() const;

#ifdef USEMPI
	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );
#endif

	void
	set_total_residue( Size total_residue );

	void
	add_native( SingleStructureDataOP native );

	void
	add_decoy( SingleStructureDataOP decoy );

	private:
	Size total_residue_;
	SingleStructureDataOPs natives_;
	SingleStructureDataOPs decoys_;
	*/

private:
	Real kT_;
	Real multiplier_; //< to make it stand up to the much more numerous sequence recovery entries

}; // PNatLigPoseOptEData


} // namespace optimize_weights
} // namespace protocols

#endif // INCLUDED_protocols_optimize_weights_PNatLigPoseOptEData_HH
