// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core::pack::rotamer_set::OptEData.hh
/// @brief Classes to store by-energy-term rotamer data for optE
/// @author Jim Havranek

#ifndef INCLUDED_protocols_optimize_weights_OptEData_hh
#define INCLUDED_protocols_optimize_weights_OptEData_hh

// Unit headers
#include <protocols/optimize_weights/OptEData.fwd.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/optimization/types.hh>

/// Utility headers
#include <utility/LexicographicalIterator.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

/// C++ headers
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace optimize_weights {


class PNatAAOptERotamerData : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~PNatAAOptERotamerData() override;
	typedef core::Real Real;
	typedef core::chemical::AA AA;
	typedef core::Size Size;

public:
	PNatAAOptERotamerData(
		AA aa_in,
		Size rot_num_in,
		utility::vector1< Real > & data_vec_in,
		utility::vector1< Real > & fixed_data_vec_in )
	:
		this_aa_( aa_in ),
		rot_number_( rot_num_in ),
		data_( data_vec_in ),
		fixed_data_( fixed_data_vec_in )
	{}

	AA this_aa()
	{
		return this_aa_;
	}

	Size rot_number()
	{
		return rot_number_;
	}

	Real
	operator [] ( Size const i ) const
	{
		return data_[ i ];
	}

	utility::vector1< Real > const &
	data()
	{
		return data_;
	}

	utility::vector1< Real > const &
	fixed_data()
	{
		return fixed_data_;
	}

private:
	AA this_aa_;
	Size rot_number_;
	utility::vector1< Real > data_;
	utility::vector1< Real > fixed_data_;
};

class PNatRotOptERotamerData : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~PNatRotOptERotamerData() override;
	typedef core::Real Real;
	typedef core::Size Size;

public:
	PNatRotOptERotamerData(
		utility::vector1< Size > const & rotamer_index,
		utility::vector1< Real > const & chi,
		utility::vector1< Real > const & free_data,
		utility::vector1< Real > const & fixed_data
	) :
		rotamer_index_( rotamer_index ),
		chi_( chi ),
		free_data_( free_data ),
		fixed_data_( fixed_data )
	{}

	utility::vector1< Size > const &
	rotamer_index() const
	{
		return rotamer_index_;
	}

	utility::vector1< Real > const &
	chi() const
	{
		return chi_;
	}

	utility::vector1< Real > const &
	free_data() const
	{
		return free_data_;
	}

	Real
	operator [] ( Size const i ) const
	{
		return free_data_[ i ];
	}

	utility::vector1< Real > const &
	fixed_data() const
	{
		return fixed_data_;
	}

private:
	utility::vector1< Size > rotamer_index_;
	utility::vector1< Real > chi_;
	utility::vector1< Real > free_data_;
	utility::vector1< Real > fixed_data_;

};

class SingleStructureData : public utility::pointer::ReferenceCount {
public:
	typedef core::Real Real;
	typedef core::Size Size;

public:
	SingleStructureData(
		utility::vector1< Real > const & free_data,
		utility::vector1< Real > const & fixed_data )
	:
		rms_(0.0),
		free_data_( free_data ),
		fixed_data_( fixed_data )
	{}

	~SingleStructureData() override = default;

	Real
	operator [] ( Size const i ) const {
		return free_data_[ i ];
	}

	utility::vector1< Real > const &
	free_data() const {
		return free_data_;
	}

	utility::vector1< Real > const &
	fixed_data() const {
		return fixed_data_;
	}

	Real rms() const { return rms_; }
	void rms( Real rms_in ) { rms_ = rms_in; }

	std::string tag() const { return tag_; }
	void tag( std::string const & setting ) { tag_ = setting; }

private:
	Real rms_;
	utility::vector1< Real > free_data_;
	utility::vector1< Real > fixed_data_;
	std::string tag_;
};

std::ostream & operator << ( std::ostream & os, PNatAAOptERotamerDataOP rd );


///////////////////////////////////////////////////////////////////////////////

class OptEPositionData : public utility::pointer::ReferenceCount
{
public:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::optimization::Multivec Multivec;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ScoreTypes ScoreTypes;

public:
	OptEPositionData();

	~OptEPositionData() override;

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
	) const = 0;

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
	) const = 0;

	/// @brief Return the upper and lower bound on the unweighted components at this
	/// position if they are larger (or smaller) than the unweighted values already in
	/// the two input EnergyMaps.
	virtual
	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const = 0;

	virtual
	Size
	size() const = 0;

	virtual
	OptEPositionDataType
	type() const = 0;

	virtual
	void
	write_to_file( std::ofstream & outfile ) const = 0;

	virtual
	void
	read_from_file( std::ifstream & infile ) = 0;

	virtual
	void
	write_to_binary_file( std::ofstream & outfile ) const = 0;

	virtual
	void
	read_from_binary_file( std::ifstream & infile ) = 0;

	virtual
	Size
	memory_use() const = 0;

	void
	tag( std::string const & tag_in );

	std::string const &
	tag() const;

#ifdef USEMPI
	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );
#endif

protected:

	/// @brief Helper function for range(); updates lower/upper_bound as needed
	/// so that score_list scores from structure are included in the range.
	void
	update_range(
		SingleStructureDataCOP structure,
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const;

private:
	std::string tag_;

};

class PNatAAOptEPositionData : public OptEPositionData
{
public:
	typedef core::chemical::AA AA;

public:
	PNatAAOptEPositionData();

	~PNatAAOptEPositionData() override;


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
	) const override;


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
	) const override;


	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const override;


	Size
	size() const override
	{
		return data_.size();
	}


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override ;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;


#ifdef USEMPI

	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );

#endif

	void set_position( Size pos_in ) {
		position_ = pos_in;
	}

	Size position() const
	{
		return position_;
	}

	void set_native_aa( AA nat_in )
	{
		native_aa_ = nat_in;
	}

	AA native_aa() const
	{
		return native_aa_;
	}

	void set_neighbor_count( Size nb_in )
	{
		neighbor_count_ = nb_in;
	}

	Size neighbor_count() const
	{
		return neighbor_count_;
	}

	void add_rotamer_line_data( PNatAAOptERotamerDataOP rot_in )
	{
		data_.push_back( rot_in );
	}

	PNatAAOptERotamerDataOPs &
	data()
	{
		return data_;
	}

	PNatAAOptERotamerDataOPs const &
	data() const
	{
		return data_;
	}

	PNatAAOptERotamerDataOPs::const_iterator rotamer_data_begin() const
	{
		return data_.begin();
	}
	PNatAAOptERotamerDataOPs::const_iterator rotamer_data_end() const
	{
		return data_.end();
	}

protected:
	/// @brief used by derived class as well -- finds the energies for the best rotamer for each amino acid
	void
	process_rotamers(
		Multivec const & vars,
		Size const num_energy_dofs,
		EnergyMap const & fixed_terms,
		ScoreTypes const & score_list,
		ScoreTypes const & fixed_score_list,
		Size const aa_range,
		utility::vector1< Real > const & dummy_set,
		utility::vector1< Real > & best_energy_by_aa,
		utility::vector1< utility::vector1< Real > > & unweighted_E_dof,
		Multivec & ref_deriv_weight
	) const;

private:
	Size position_;
	AA native_aa_;
	Size neighbor_count_;
	PNatAAOptERotamerDataOPs data_;

};

class PSSMOptEPositionData : public PNatAAOptEPositionData
{
public:
	typedef PNatAAOptEPositionData parent;

public:

	PSSMOptEPositionData();

	~PSSMOptEPositionData() override;

	void set_pssm_probabilities( utility::vector1< Real > const & pssm_probs );


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
	) const override;


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
	) const override;


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override ;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;

#ifdef USEMPI
	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );
#endif

private:
	utility::vector1< Real > pssm_probabilities_; // sum to 1

};


class PNatRotOptEPositionData : public OptEPositionData
{
public:
	PNatRotOptEPositionData();

	~PNatRotOptEPositionData() override;


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
	) const override;


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
	) const override;

	Real
	process_score(
		std::ostream & ostr,
		bool print,
		Multivec const & component_weights,
		Multivec const & vars,
		Multivec & dE_dvars,
		/// Basically, turn over all the private data from OptEMultiFunc
		Size const num_energy_dofs,
		int const,
		int const,
		EnergyMap const & fixed_terms,
		ScoreTypes const &,
		ScoreTypes const & fixed_score_list
	) const;


	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const override;


	Size
	size() const override;


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;


#ifdef USEMPI
	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );
#endif

	void
	set_native_rotamer_index( utility::vector1< Size > const & native_rotamer_index );

	void
	set_native_rotamer_chi( utility::vector1< Real > const & native_chi );

	void
	set_native_chi_periodicity( utility::vector1< Real > const & native_chi_periodicity );

	bool
	count_rotamer_as_native( PNatRotOptERotamerDataOP rotamer ) const;

	void
	set_rotamer_well_counts( utility::vector1< Size > const & rotamer_well_counts );

	void add_rotamer_line_data( PNatRotOptERotamerDataOP rot_in );

	PNatRotOptERotamerDataOPs &
	data();

	PNatRotOptERotamerDataOPs const &
	data() const;

	PNatRotOptERotamerDataOPs::const_iterator rotamer_data_begin() const;
	PNatRotOptERotamerDataOPs::const_iterator rotamer_data_end() const;


	core::chemical::AA aa() const;
	core::chemical::AA & aa();

	Real phi() const;
	Real psi() const;
	Real & phi();
	Real & psi();

private:

	Size
	rotamer_index_2_well_id( utility::vector1< Size > const & rotamer_index ) const;

	Size
	rotamer_index_2_well_id( utility::LexicographicalIterator const & lexiter ) const;

	bool
	is_native_rotamer_well(
		utility::vector1< Size > const & rotamer_index
	) const;

	bool
	is_native_rotamer_well(
		utility::LexicographicalIterator const & lexiter
	) const;


private:

	utility::vector1< Size > native_rotamer_index_;
	utility::vector1< Size > rotamer_well_counts_;
	int n_wells_;
	core::chemical::AA aa_;
	Real phi_;
	Real psi_;
	PNatRotOptERotamerDataOPs data_;
	utility::vector1< Real > native_chi_;
	utility::vector1< Real > native_chi_periodicity_;
};

class PNatStructureOptEData : public OptEPositionData
{
public:
	PNatStructureOptEData();
	~PNatStructureOptEData() override;


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
	) const override;


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
	) const override;

	Real
	process_score(
		std::ostream & ostr,
		bool print, // generate output?
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
	) const;


	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const override;


	Size
	size() const override;


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;

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

	void
	n_top_natives_to_score( Size n_top );

	Size
	n_top_natives_to_score() const;

	void
	set_normalize_decoy_stddev( bool setting );

	void
	set_initial_decoy_stddev( Real setting );

	Real
	nativeness( Real rms ) const;

	static
	void
	set_nativeness_low( Real nativeness_rms_low );

	static
	void
	set_nativeness_high( Real nativeness_rms_low );

	static
	Real
	nativeness_low();

	static
	Real
	nativeness_high();


protected:
	Size total_residue_;
	SingleStructureDataOPs natives_;
	SingleStructureDataOPs decoys_;
	Size n_top_natives_to_score_;
	Size n_high_entropy_decoys_;
	bool normalize_decoy_stddev_;
	Real initial_decoy_stddev_;
	static Real nativeness_rms_low_; // Above this rms, nativeness starts to decline
	static Real nativeness_rms_high_; // Above this rms, nativeness is zero
	Real nativeness_sum_;
	static Real const high_entropy_rms_cutoff_;
};


class DDGMutationOptEData : public OptEPositionData
{
public:
	typedef core::chemical::AA AA;

public:
	DDGMutationOptEData();
	~DDGMutationOptEData() override;


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
	) const override;


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
	) const override;

	Real
	process_score(
		std::ostream & ostr,
		bool print, // generate output?
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
	) const;


	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const override;


	Size
	size() const override;


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;

#ifdef USEMPI
	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );
#endif

	void
	set_wt_aa( AA wt_aa );

	void
	set_mut_aa( AA mut_aa );

	void
	set_experimental_ddg( Real ddg );

	void
	add_wt( SingleStructureDataOP wt );

	void
	add_mutant( SingleStructureDataOP mut );

protected:
	Real experimental_ddG_;
	AA wt_aa_;
	AA mut_aa_;
	SingleStructureDataOPs wts_;
	SingleStructureDataOPs muts_;

};

struct WeightRangeConstraint
{
public:

	bool active_;
	core::Real min_weight_;
	core::Real max_weight_;
	core::Real spring_constant_;

public:
	WeightRangeConstraint() :
		active_( true ),
		min_weight_( 0 ),
		max_weight_( 10 ),
		spring_constant_( 1000 )
	{}

};

class ConstraintedOptimizationWeightFunc : public OptEPositionData
{
public:
	ConstraintedOptimizationWeightFunc();

	ConstraintedOptimizationWeightFunc( ScoreTypes const & score_list );

	~ConstraintedOptimizationWeightFunc() override;
	void
	initialize_constraints_from_file( std::ifstream & infile );


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
	) const override;


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
	) const override;


	void
	range(
		ScoreTypes const & free_score_list,
		ScoreTypes const & fixed_score_list,
		EnergyMap & lower_bound,
		EnergyMap & upper_bound
	) const override;


	Size
	size() const override;


	OptEPositionDataType
	type() const override;


	void
	write_to_file( std::ofstream & outfile ) const override ;


	void
	read_from_file( std::ifstream & infile ) override;


	void
	write_to_binary_file( std::ofstream & outfile ) const override;


	void
	read_from_binary_file( std::ifstream & infile ) override;


	Size
	memory_use() const override;


#ifdef USEMPI

	virtual
	void
	send_to_node( int const destination_node, int const tag ) const;

	virtual
	void
	receive_from_node( int const source_node, int const tag );

#endif

private:
	ScoreTypes free_terms_;
	EnergyMap  free_term_map_;
	utility::vector1< WeightRangeConstraint > free_term_constraints_;

};

class OptEPositionDataFactory
{
public:
	static
	OptEPositionDataOP
	create_position_data( OptEPositionDataType const type );

	static
	std::string const &
	optE_type_name( OptEPositionDataType const type );

	static
	bool
	is_optE_type_name( std::string const & name );

	static
	OptEPositionDataType
	optE_type_from_name( std::string const & name );

private:
	static
	void
	initialize_optE_type_name_map();

	static bool optE_type_name_map_initialized_;
	static utility::vector1< std::string > optE_type_2_optE_type_name_;
	static std::map< std::string, OptEPositionDataType > optE_type_name_map_;
};

///////////////////////////////////////////////////////////////////////////////


class OptEData : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	~OptEData() override;
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::optimization::Multivec Multivec;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::ScoreTypes ScoreTypes;

public:
	OptEData() {}

	OptEData(
		ScoreTypes const & fixed_score_list,
		ScoreTypes const & free_score_list )
	:
		fixed_energy_terms_( fixed_score_list ),
		energy_terms_( free_score_list )
	{}

	Size num_positions() const
	{
		return data_.size();
	}

	Size num_rotamers() const;

	void add_position_data( OptEPositionDataOP pos_data_in )
	{
		data_.push_back( pos_data_in );
	}

	OptEPositionDataOPs::const_iterator position_data_begin() const
	{
		return data_.begin();
	}
	OptEPositionDataOPs::const_iterator position_data_end() const
	{
		return data_.end();
	}

	ScoreTypes const & fixed_energy_terms() const
	{
		return fixed_energy_terms_;
	}
	ScoreTypes const & energy_terms() const
	{
		return energy_terms_;
	}

	void write_to_file( std::string filename = "opte.data" ) const;
	void read_from_file( std::string filename );
	void write_to_binary_file( std::string filename = "opte.data" ) const;
	void read_from_binary_file( std::string filename );

private:
	ScoreTypes fixed_energy_terms_;
	ScoreTypes energy_terms_;
	OptEPositionDataOPs data_;
};

} // namespace optimize_weights
} // namespace protocols

#endif // INCLUDE_protocols_optimize_weights_OptEData_HH
