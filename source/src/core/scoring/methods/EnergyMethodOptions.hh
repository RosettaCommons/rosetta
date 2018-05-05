// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley

/// NOTE-- this file includes both string and map, use .fwd.hh if
/// you can!


#ifndef INCLUDED_core_scoring_methods_EnergyMethodOptions_hh
#define INCLUDED_core_scoring_methods_EnergyMethodOptions_hh

// Unit headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>

#include <core/types.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/etable/EtableOptions.fwd.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/FreeDOF_Options.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/chemical/AA.hh>

/// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

// C++ Headers
#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace methods {

/// add more options here
/// NOTE: If you add an option, make sure you also update the constructor,
/// the assignment operator, the == comparison operator, and the show method in the .cc file!
/// right now this class should be pretty light-weight since a copy is held inside ScoreFunctionInfo


class EnergyMethodOptions : public utility::pointer::ReferenceCount {
public:

	/// @brief Default constructor, reads from the global options system
	EnergyMethodOptions();

	/// @brief Initialize an EnergyMethodOptions object from a (possibly local) option collection
	EnergyMethodOptions( utility::options::OptionCollection const & options );

	/// copy constructor
	EnergyMethodOptions( EnergyMethodOptions const & src );

	virtual
	~EnergyMethodOptions();

	/// copy operator
	EnergyMethodOptions &
	operator=( EnergyMethodOptions const & src );

	/// clone
	EnergyMethodOptionsOP
	clone() const{
		return EnergyMethodOptionsOP( new EnergyMethodOptions( *this ) );
	}

	/// @brief Initialize a new EnergyMethodOptions with defaults from the global option collection
	void
	initialize_from_options();

	/// @brief Initialize a new EnergyMethodOptions with defaults from a (possibly local) option collection
	void
	initialize_from_options( utility::options::OptionCollection const & options );

	/// @brief Append the option keys read by the initialize_from_options method to the input option-key list
	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );


	/// @brief Get the nth aa_composition setup file name from the list of setup files.
	///
	inline std::string const & aa_composition_setup_file( core::Size const index ) const {
		runtime_assert_string_msg( index > 0 && index <= aa_composition_setup_files_.size(), "Error in core::scoring::methods::EnergyMethodOptions::aa_composition_setup_file(): The index of the file requested is greater than the number of filenames stored." );
		return aa_composition_setup_files_[index];
	}

	/// @brief Get the number of aa_composition setup files.
	///
	inline core::Size aa_composition_setup_file_count() const {
		return aa_composition_setup_files_.size();
	}

	/// @brief Set the aa_composition setup file names.
	/// @details Overrides existing.
	inline void set_aa_composition_setup_files( utility::vector1 < std::string > const &input_filenames ) {
		aa_composition_setup_files_ = input_filenames;
		return;
	}

	/// @brief Appends additional files to the aa_composition setup file names.
	/// @details Does not override existing.
	inline void append_aa_composition_setup_files( utility::vector1 < std::string > const &input_filenames ) {
		for ( core::Size i=1, imax=input_filenames.size(); i<=imax; ++i ) {
			aa_composition_setup_files_.push_back( input_filenames[i] );
		}
		return;
	}

	/// @brief Get the nth netcharge setup file name from the list of setup files.
	///
	inline std::string const & netcharge_setup_file( core::Size const index ) const {
		runtime_assert_string_msg( index > 0 && index <= netcharge_setup_files_.size(), "Error in core::scoring::methods::EnergyMethodOptions::netcharge_setup_file(): The index of the file requested is greater than the number of filenames stored." );
		return netcharge_setup_files_[index];
	}

	/// @brief Get the number of netcharge setup files.
	///
	inline core::Size netcharge_setup_file_count() const {
		return netcharge_setup_files_.size();
	}

	/// @brief Set the netcharge setup file names.
	/// @details Overrides existing.
	inline void set_netcharge_setup_files( utility::vector1 < std::string > const &input_filenames ) {
		netcharge_setup_files_ = input_filenames;
		return;
	}

	/// @brief Appends additional files to the netcharge setup file names.
	/// @details Does not override existing.
	inline void append_netcharge_setup_files( utility::vector1 < std::string > const &input_filenames ) {
		for ( core::Size i=1, imax=input_filenames.size(); i<=imax; ++i ) {
			netcharge_setup_files_.push_back( input_filenames[i] );
		}
		return;
	}

	/// @brief Get the penalty for each aspartimide-forming two-residue sequence.
	/// @details Used by the aspartimide_penalty score term.
	inline core::Real const &
	aspartimide_penalty_value() const {
		return aspartimide_penalty_value_;
	}

	std::string const &
	etable_type() const;


	void
	etable_type( std::string const & type );

	bool analytic_etable_evaluation() const;
	void analytic_etable_evaluation( bool setting );


	std::string const &
	unfolded_energies_type() const;


	void
	unfolded_energies_type( std::string const & type );

	///
	std::string const &
	split_unfolded_label_type() const;

	///
	void
	split_unfolded_label_type(std::string const & label_type);

	///
	std::string const &
	split_unfolded_value_type() const;

	///
	void
	split_unfolded_value_type(std::string const & value_type);

	///
	bool
	exclude_protein_protein_fa_elec() const;

	void
	exclude_protein_protein_fa_elec( bool const setting );

	bool
	exclude_RNA_RNA_fa_elec() const;

	void
	exclude_RNA_RNA_fa_elec( bool const setting );

	bool
	exclude_RNA_protein_fa_elec() const;

	void
	exclude_RNA_protein_fa_elec( bool const setting );


	bool
	exclude_monomer_fa_elec() const;


	void
	exclude_monomer_fa_elec( bool const setting );

	/// @brief The maximum (all atom) distance at which fa_elec is non-zero
	core::Real
	elec_max_dis() const;

	void
	elec_max_dis( core::Real setting );

	/// @brief The minimium (all atom) distance for which fa_elec changes with distances
	core::Real
	elec_min_dis() const;

	void
	elec_min_dis( core::Real setting );

	/// @brief The dielectric used for the fa_elec term
	core::Real
	elec_die() const;

	void
	elec_die( core::Real setting );

	/// @brief Should fa_elec use a constant (non-distance dependant) dielectric?
	bool
	elec_no_dis_dep_die() const;

	void
	elec_no_dis_dep_die( bool setting );

	/// @brief Should fa_elec/gpelec use a sigmoidal dielectric?
	bool
	elec_sigmoidal_die() const;

	void
	elec_sigmoidal_die( bool const setting );

	void
	elec_sigmoidal_die_params(Real &D, Real &D0, Real &S) const;

	void
	set_elec_sigmoidal_die_params( Real D, Real D0, Real S );

	bool
	smooth_fa_elec() const;

	void
	smooth_fa_elec( bool setting );

	std::string
	grpelec_fade_type() const;

	void
	grpelec_fade_type( std::string setting );

	core::Real
	grpelec_fade_param1() const;

	void
	grpelec_fade_param1( core::Real setting );

	core::Real
	grpelec_fade_param2() const;

	void
	grpelec_fade_param2( core::Real setting );

	bool
	grpelec_fade_hbond() const;

	void
	grpelec_fade_hbond( bool setting );

	bool
	grp_cpfxn() const;

	void
	grp_cpfxn( bool setting );

	std::string
	elec_group_file() const;

	void
	elec_group_file( std::string setting );

	bool
	grpelec_context_dependent() const;

	void
	grpelec_context_dependent( bool setting );

	bool
	use_polarization() const;

	void
	use_polarization( bool setting );

	bool
	use_gen_kirkwood() const;

	void
	use_gen_kirkwood( bool setting );

	Real
	protein_dielectric() const;

	void
	protein_dielectric( Real setting );

	Real
	water_dielectric() const;

	void
	water_dielectric( Real setting );

	bool
	exclude_DNA_DNA() const;


	void
	exclude_DNA_DNA( bool const setting );


	bool
	exclude_intra_res_protein() const;


	void
	exclude_intra_res_protein( bool const setting );

	bool
	count_pair_hybrid() const;

	void
	count_pair_hybrid( bool const setting );


	bool
	put_intra_into_total() const;


	void
	put_intra_into_total( bool const setting );


	core::Size
	geom_sol_interres_path_distance_cutoff() const;


	void
	geom_sol_interres_path_distance_cutoff( core::Size const setting );


	core::Size
	geom_sol_intrares_path_distance_cutoff() const;


	void
	geom_sol_intrares_path_distance_cutoff( core::Size const setting );

	bool
	eval_intrares_elec_ST_only() const;

	void
	eval_intrsres_elec_ST_only( bool setting );

	bool
	envsmooth_zero_negatives() const;

	void
	envsmooth_zero_negatives( bool const setting );

	/// @brief Read access to the hbond options object
	hbonds::HBondOptions const &
	hbond_options() const;

	/// @brief non-const access to the hbond options object
	hbonds::HBondOptions &
	hbond_options();

	/// @brief Set the hbond options object -- makes a deep copy
	void
	hbond_options( hbonds::HBondOptions const & opts );

	/// @brief Read access to the etable options object
	etable::EtableOptions const &
	etable_options() const;

	/// @brief non-const access to the etable options object
	etable::EtableOptions &
	etable_options();

	/// @brief Set the etable options object -- makes a deep copy
	void
	etable_options( etable::EtableOptions const & opts );

	/// @brief Read access to the RNA options object
	rna::RNA_EnergyMethodOptions const &
	rna_options() const;

	/// @brief non-const access to the RNA options object
	rna::RNA_EnergyMethodOptions &
	rna_options();

	/// @brief Set the FreeDOF options object -- makes a deep copy
	void
	rna_options( rna::RNA_EnergyMethodOptions const & opts );

	/// @brief Read access to the FreeDOF options object
	methods::FreeDOF_Options const &
	free_dof_options() const;

	/// @brief non-const access to the FreeDOF options object
	methods::FreeDOF_Options &
	free_dof_options();

	/// @brief Set the FreeDOF options object -- makes a deep copy
	void
	free_dof_options( methods::FreeDOF_Options const & opts );

	std::string const & pb_bound_tag() const;
	std::string & pb_bound_tag();
	void pb_bound_tag( std::string const & tag );

	std::string const & pb_unbound_tag() const;
	std::string & pb_unbound_tag();
	void pb_unbound_tag( std::string const & tag );

	/// @brief Should glyceine's Ramachandran and P_AA_PP tables be symmetrized (e.g. for scoring in a mixed D/L context)?
	/// @details Default false.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool symmetric_gly_tables() const;

	/// @brief Set whether glyceine's Ramachandran and P_AA_PP tables should be symmetrized (e.g. for scoring in a mixed D/L context).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void symmetric_gly_tables( bool const setting );

	/// @brief Get the number of cones in which a voxel must lie in order for that voxel to be considered
	/// to be "buried".
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Size voids_penalty_energy_containing_cones_cutoff() const;

	/// @brief Get the cone distance cutoff for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Real voids_penalty_energy_cone_distance_cutoff() const;

	/// @brief Get the cone dot product cutoff for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Real voids_penalty_energy_cone_dotproduct_cutoff() const;

	/// @brief Get the voxel grid padding for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Real voids_penalty_energy_voxel_grid_padding() const;

	/// @brief Get the voxel size for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Real voids_penalty_energy_voxel_size() const;

	/// @brief Get whether we're prohibiting evaluation of the voids_penalty score term outside
	/// of the context of the packer.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	bool voids_penalty_energy_disabled_except_during_packing() const;

	/// @brief Set the number of cones in which a voxel must lie in order for that voxel to be considered
	/// to be "buried".
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void voids_penalty_energy_containing_cones_cutoff( core::Size const setting );

	/// @brief Set the cone distance cutoff for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void voids_penalty_energy_cone_distance_cutoff( core::Real const &setting );

	/// @brief Set the cone dot product cutoff for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void voids_penalty_energy_cone_dotproduct_cutoff( core::Real const &setting );

	/// @brief Set the voxel grid padding for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void voids_penalty_energy_voxel_grid_padding( core::Real const &setting );

	/// @brief Set the voxel size for the voids penalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void voids_penalty_energy_voxel_size( core::Real const &setting );

	/// @brief Set whether we're prohibiting evaluation of the voids_penalty score term outside
	/// of the context of the packer.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void voids_penalty_energy_disabled_except_during_packing( bool const setting );

	/// @brief Get the bonus function shape for the hbnet energy term.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	std::string const & hbnet_bonus_function_ramping() const;

	/// @brief Set the bonus function shape for the hbnet energy term.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void hbnet_bonus_function_ramping( std::string const & setting );

	/// @brief Get the maximum hydrogen bond network size, beyond which the hbnet score term yields no futher bonus.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	core::Size hbnet_max_network_size() const;

	/// @brief Set the maximum hydrogen bond network size, beyond which the hbnet score term yields no futher bonus.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void hbnet_max_network_size( core::Size const setting );

	bool loop_close_use_6D_potential() const;
	void loop_close_use_6D_potential( bool const setting );

	bool fa_stack_base_all() const;
	void fa_stack_base_all( bool const setting );

	/// @brief Set the angle exponent for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_cone_angle_exponent( core::Real const setting );

	/// @brief Set the angle shift factor for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_cone_angle_shift_factor( core::Real const setting );

	/// @brief Set the distance exponent for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_cone_dist_exponent( core::Real const setting );

	/// @brief Set the distance midpoint for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_cone_dist_midpoint( core::Real const setting );

	/// @brief Set the number of cones in which a point must lie to be considered "buried"
	/// by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_burial_threshold( core::Real const setting );

	/// @brief Set the energy threshold above which a hydrogen bond is not counted.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void buried_unsatisfied_penalty_hbond_energy_threshold( core::Real const setting );

	/// @brief Get the angle exponent for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_cone_angle_exponent() const;

	/// @brief Get the angle shift factor for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_cone_angle_shift_factor() const;

	/// @brief Get the distance exponent for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_cone_dist_exponent() const;

	/// @brief Get the distance midpoint for calculating burial by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_cone_dist_midpoint() const;

	/// @brief Get the number of cones in which a point must lie to be considered "buried"
	/// by the method of sidechain neighbor cones.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_burial_threshold() const;

	/// @brief Get the energy threshold above which a hydrogen bond is not counted.
	/// @details Used by the BuriedUnsatPenalty energy.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	core::Real buried_unsatisfied_penalty_hbond_energy_threshold() const;

	utility::vector1< core::Real > const & get_density_sc_scale_byres() const;
	void set_density_sc_scale_byres(core::Real newscscale);
	void set_density_sc_scale_byres(core::chemical::AA aa, core::Real newscscale);

	/// @brief  This is used in the construction of the VDW_Energy's AtomVDW object
	std::string const &
	atom_vdw_atom_type_set_name() const;


	void
	atom_vdw_atom_type_set_name( std::string const & setting );


	core::Size
	cst_max_seq_sep() const;


	void
	cst_max_seq_sep( Size const setting );

	/// deprecated
	utility::vector1<std::string> const &
	bond_angle_central_atoms_to_score() const;

	/// depricated
	void
	bond_angle_central_atoms_to_score(
		utility::vector1<std::string> const & atom_names);

	scoring::mm::MMBondAngleResidueTypeParamSetOP
	bond_angle_residue_type_param_set();

	scoring::mm::MMBondAngleResidueTypeParamSetCOP
	bond_angle_residue_type_param_set() const;

	void
	bond_angle_residue_type_param_set(
		scoring::mm::MMBondAngleResidueTypeParamSetOP param_set);

	void set_strand_strand_weights(
		int ss_lowstrand,
		int ss_cutoff);


	SecondaryStructureWeights const &
	secondary_structure_weights() const;


	SecondaryStructureWeights &
	secondary_structure_weights();


	bool
	has_method_weights( ScoreType const & type ) const;


	utility::vector1< Real > const &
	method_weights( ScoreType const & type ) const;


	void
	set_method_weights(
		ScoreType const & type,
		utility::vector1< Real > const & wts);

	/// @brief get the harmonic bond angle and bond-length spring constants
	void
	get_cartesian_bonded_parameters( Real &len, Real &ang, Real &tors, Real &proton , Real &imp ) const {
		len=cartbonded_len_;
		ang=cartbonded_ang_;
		tors=cartbonded_tors_;
		proton=cartbonded_proton_;
		imp=cartbonded_improper_;
	}

	/// @brief set the harmonic bond angle and bond-length spring constants
	void
	set_cartesian_bonded_parameters( Real len, Real ang, Real tors, Real proton, Real imp ) {
		cartbonded_len_=len;
		cartbonded_ang_=ang;
		cartbonded_tors_=tors;
		cartbonded_proton_=proton;
		cartbonded_improper_=imp;
	}

	/// @brief get the harmonic bond angle and bond-length spring constants
	bool get_cartesian_bonded_linear() const {
		return cartbonded_linear_;
	}

	/// @brief set the harmonic bond angle and bond-length spring constants
	void set_cartesian_bonded_linear( bool lin_in ) {
		cartbonded_linear_ = lin_in;
	}

	core::Real
	ordered_wat_penalty() const;

	core::Real
	ordered_pt_wat_penalty() const;

	/// used inside ScoreFunctionInfo::operator==
	friend
	bool
	operator==( EnergyMethodOptions const & a, EnergyMethodOptions const & b );

	/// used inside ScoreFunctionInfo::operator==
	friend
	bool
	operator!=( EnergyMethodOptions const & a, EnergyMethodOptions const & b );


	void
	show( std::ostream & out ) const;

	static
	void
	write_score_function_method_options_table_schema(
		utility::sql_database::sessionOP db_session
	);

	void
	insert_score_function_method_options_rows(
		Size batch_id,
		std::string const & score_function_name,
		utility::sql_database::sessionOP db_session) const;

private:
	/// expand this to a class and include ss weights inside
	typedef  std::map< ScoreType, utility::vector1< Real > > MethodWeights;

private:

	/////////////////////////////////////////////////
	// IMPORTANT NOTE!  If you add an option, make sure you also update the constructor,
	// the assignment operator, the == comparison operator, and the show method in the .cc file!
	/////////////////////////////////////////////////
	utility::vector1 < std::string > aa_composition_setup_files_;
	utility::vector1 < std::string > netcharge_setup_files_;
	core::Real aspartimide_penalty_value_;
	std::string atom_vdw_atom_type_set_name_;
	std::string unfolded_energies_type_;
	std::string split_unfolded_label_type_;
	std::string split_unfolded_value_type_;
	MethodWeights method_weights_;
	SecondaryStructureWeights ss_weights_;
	bool exclude_protein_protein_fa_elec_;
	bool exclude_RNA_RNA_fa_elec_;
	bool exclude_RNA_protein_fa_elec_;
	bool exclude_monomer_fa_elec_;
	core::Real elec_max_dis_;
	core::Real elec_min_dis_;
	core::Real elec_die_;
	bool elec_no_dis_dep_die_;
	bool smooth_fa_elec_;
	bool elec_sigmoidal_die_;
	Real elec_sigmoidal_D_, elec_sigmoidal_D0_, elec_sigmoidal_S_;
	std::string grpelec_fade_type_;
	core::Real grpelec_fade_param1_;
	core::Real grpelec_fade_param2_;
	bool grpelec_fade_hbond_;
	bool grp_cpfxn_;
	std::string elec_group_file_;
	bool grpelec_context_dependent_;
	bool use_polarization_;
	bool use_gen_kirkwood_;
	core::Real protein_dielectric_;
	core::Real water_dielectric_;
	bool exclude_DNA_DNA_;
	bool exclude_intra_res_protein_;
	bool count_pair_hybrid_;
	bool put_intra_into_total_;
	core::Size geom_sol_interres_path_distance_cutoff_;
	core::Size geom_sol_intrares_path_distance_cutoff_;
	bool eval_intrares_elec_ST_only_;
	bool envsmooth_zero_negatives_;
	hbonds::HBondOptionsOP hbond_options_;
	core::scoring::etable::EtableOptionsOP etable_options_;
	rna::RNA_EnergyMethodOptionsOP rna_options_;
	methods::FreeDOF_OptionsOP free_dof_options_;
	core::Size cst_max_seq_sep_;
	core::Real cartbonded_len_, cartbonded_ang_, cartbonded_tors_, cartbonded_proton_, cartbonded_improper_;
	bool cartbonded_linear_;
	std::string pb_bound_tag_;
	std::string pb_unbound_tag_;
	utility::vector1< core::Real > fastdens_perres_weights_;
	core::Real ordered_wat_penalty_;    //fpd -> penalty for removing water from bulk (HOH_V->HOH)
	core::Real ordered_pt_wat_penalty_; //fpd -> penalty for removing point water from bulk (PWAT_V->PWAT)
	bool symmetric_gly_tables_;
	bool loop_close_use_6D_potential_;
	bool fa_stack_base_all_;

	//Options for the BuriedUnsatPenalty energy:
	core::Real buried_unsatisfied_penalty_cone_angle_exponent_;
	core::Real buried_unsatisfied_penalty_cone_angle_shift_factor_;
	core::Real buried_unsatisfied_penalty_cone_dist_exponent_;
	core::Real buried_unsatisfied_penalty_cone_dist_midpoint_;
	core::Real buried_unsatisfied_penalty_burial_threshold_;
	core::Real buried_unsatisfied_penalty_hbond_energy_threshold_;

	//Options for the VoidsPenaltyEnergy:
	core::Size voids_penalty_energy_containing_cones_cutoff_;
	core::Real voids_penalty_energy_cone_dotproduct_cutoff_;
	core::Real voids_penalty_energy_cone_distance_cutoff_;
	core::Real voids_penalty_energy_voxel_size_;
	core::Real voids_penalty_energy_voxel_grid_padding_;
	bool voids_penalty_energy_disabled_except_during_packing_;

	//Options for the HBNetEnergy:
	std::string hbnet_bonus_ramping_function_;
	core::Size hbnet_max_network_size_;

	/// deprecated
	utility::vector1<std::string> bond_angle_central_atoms_to_score_;
	core::scoring::mm::MMBondAngleResidueTypeParamSetOP bond_angle_residue_type_param_set_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


std::ostream &
operator<< ( std::ostream & out, EnergyMethodOptions const & options );

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_methods_EnergyMethodOptions )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_ScoreFunction_HH
