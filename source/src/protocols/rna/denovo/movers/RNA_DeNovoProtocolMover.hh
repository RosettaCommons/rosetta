// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Andy Watkins, amw579@stanford.edu


#ifndef INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoProtocolMover_HH
#define INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoProtocolMover_HH

#include <utility/VirtualBase.hh>
#include <protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.fwd.hh>
#include <core/import_pose/RNA_DeNovoParameters.fwd.hh>
#include <core/pose/rna/RNA_SecStruct.fwd.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/types.hh>
#include <core/pose/rna/BasePair.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <tuple>

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

class RNA_DeNovoProtocolMover : public protocols::moves::Mover {

	using Size = core::Size;

public:

	//constructor
	RNA_DeNovoProtocolMover();

	//destructor
	~RNA_DeNovoProtocolMover() override;

public:

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	// void
	// initialize_from_options( utility::options::OptionCollection const & opts );
	// void
	// initialize_from_command_line();

	static
	void
	list_options_read( utility::options::OptionKeyList & opts );

	core::import_pose::options::RNA_DeNovoProtocolOptionsCOP options() const { return  options_; }
	core::import_pose::RNA_DeNovoParametersCOP rna_params() const { return rna_params_; }

	void set_fasta_files( utility::vector1< std::string > const & fasta_files ) { fasta_files_ = fasta_files; }
	void set_input_pdbs( utility::vector1< std::string > const & input_pdbs ); // This must propagate to options
	void set_input_silent_files( utility::vector1< std::string > const & input_silent_files ); // This must propagate to options
	void set_sequence_strings( utility::vector1< std::string > const & sequence_strings ) { sequence_strings_ = sequence_strings; }
	void set_minimize_rna( bool const minimize_rna ) { minimize_rna_ = minimize_rna; minimize_rna_has_been_specified_ = true; }
	void set_helical_substructs( utility::vector1< std::string > const & helical_substructs ) { helical_substructs_ = helical_substructs; }
	void set_dock_chunks( utility::vector1< std::string > const & dock_chunks ) { dock_chunks_ = dock_chunks; }
	void set_input_initialization_pdbs( utility::vector1< std::string > const & input_initialization_pdbs ); // This must propagate to options

	// void set_native_pose( core::pose::PoseOP native_pose ) { native_pose_ = native_pose; }
	void set_align_pdb( std::string const & align_pdb );
	void set_nstruct( Size const nstruct );
	void set_silent_file( std::string const & silent_file );


	void
	de_novo_setup_from_command_line();

	void
	de_novo_setup_from_options( utility::options::OptionCollection const & opts );

	void
	de_novo_setup_from_tag( utility::tag::TagCOP const & tag );

	void
	set_protocol_options( core::import_pose::options::RNA_DeNovoProtocolOptionsOP const & protocol_options ) {
		options_ = protocol_options;
	}

	void
	set_residue_type_set( core::chemical::ResidueTypeSetCOP const & rts ) {
		rsd_set_ = rts;
	}

private:


	void
	initialize_scorefxn( core::pose::Pose & pose );



	/////// Steps ////////
	void
	initialize_sequence_information(
		int const offset,
		bool const edensity_mapfile_provided,
		bool const rna_protein_docking,
		bool const virtual_anchor_provided,
		bool const cutpoint_open_provided,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & possible_cutpoint_open_numbering,
		core::pose::full_model_info::FullModelParametersOP & full_model_parameters,
		utility::vector1< core::Size > & cutpoint_open_in_full_model );

	void
	check_secstructs(
		bool const modeling_with_density,
		bool const rna_protein_docking,
		std::string const & sequence,
		core::pose::rna::RNA_SecStruct const & secstruct,
		core::pose::rna::RNA_SecStruct const & secstruct_general );

	void
	input_numbering_setup(
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_res_resnum_and_chain,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_silent_res_resnum_and_chain,
		std::string const & sequence,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > & input_res,
		utility::vector1< utility::vector1< int > > & resnum_list,
		utility::vector1< core::Size > & input_res_initialize,
		utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
		utility::vector1< core::Size > & dock_chunks_res
	);

	void
	input_pdb_numbering_setup(
		std::string const & sequence,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > & input_res,
		utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
		utility::vector1< core::Size > & dock_chunks_res,
		utility::vector1< utility::vector1< int > > & resnum_list,
		utility::vector1< core::Size > & input_res_user_defined,
		Size & input_res_user_defined_count
	);

	void
	input_initialization_pdb_numbering_setup(
		std::string const & sequence,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > & input_initialization_res
	);
	void
	input_silent_numbering_setup(
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & input_silent_res_resnum_and_chain,
		std::string const & sequence,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > const & input_res_user_defined,
		core::Size const & input_res_user_defined_count,
		utility::vector1< core::Size > & input_res,
		utility::vector1< utility::vector1< core::Size > > & helical_substruct_res,
		utility::vector1< core::Size > & dock_chunks_res,
		utility::vector1< utility::vector1< int > > & resnum_list
	);

	void
	setup_obligate_pair(
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_res,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< utility::vector1< int > > const & resnum_list,
		core::pose::rna::RNA_SecStruct const & secstruct,
		core::pose::rna::RNA_SecStruct const & secstruct_general,
		utility::vector1< core::Size > const & cutpoint_open_in_full_model,
		utility::vector1< core::Size > & obligate_pair,
		utility::vector1< std::string > & obligate_pair_explicit,
		utility::vector1< core::Size > & domain_map );

	void
	initial_pose_setup(
		std::string const & in_path,
		bool const lores_sfxn_provided,
		bool const native_provided,
		bool const edensity_mapfile_provided,
		bool const working_native_provided,
		std::string const & working_native,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > const & working_res,
		core::pose::Pose & full_pose,
		bool & is_rna_and_protein );

	void
	constraint_setup(
		bool const constraint_file_provided,
		std::string const & constraint_file,
		utility::vector1< core::Size > const & working_res,
		core::pose::Pose & full_pose );

	void
	refine_working_obligate_pairs(
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & remove_obligate_pair_rc,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > const & obligate_pair,
		utility::vector1< std::string > & obligate_pair_explicit,
		utility::vector1< core::Size > const & working_res,
		utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs
	);

	void
	setup_working_chain_connections(
		utility::vector1< std::string > const & chain_connections,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters,
		utility::vector1< core::Size > const & working_res,
		utility::vector1< std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > > & working_chain_connections
	);

	void
	setup_pairing_sets(
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & working_stems,
		utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs,
		core::pose::rna::RNA_BasePairList & rna_pairing_list,
		utility::vector1< utility::vector1< core::Size > > & working_obligate_pair_sets,
		utility::vector1< utility::vector1< core::Size > > & working_stem_pairing_sets
	);

	void
	setup_ft_from_silent(
		bool const ft_from_silent_provided,
		std::string const & ft_silent_file,
		bool const tag_provided,
		utility::vector1< std::string > const & tag,
		bool & use_fold_tree_from_silent_file,
		core::kinematics::FoldTree & fold_tree_from_silent_file
	);

	void
	setup_final_res_lists(
		utility::vector1< core::Size > const & working_res,
		utility::vector1< core::Size > const & cutpoint_open_in_full_model,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_res_rc,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & extra_minimize_chi_res_rc,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & output_jump_res_rc,
		bool const minimize_rna_specified,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & syn_chi_rc,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & anti_chi_rc,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & block_stack_above_rc,
		std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > > const & block_stack_below_rc,
		core::pose::full_model_info::FullModelParametersOP const & full_model_parameters // can use a const OP -- not a COP!
	);
	/////// Steps ////////


	utility::vector1< core::Size >
	working_res_map( utility::vector1< core::Size > const & vec,
		utility::vector1< core::Size > const & working_res,
		bool const leave_out_last_working_residue = false ) const;

	std::string
	working_res_map( std::string const & seq,
		utility::vector1< core::Size > const & working_res,
		bool const annotations_in_brackets = true ) const;

	core::pose::rna::RNA_SecStruct
	working_res_map( core::pose::rna::RNA_SecStruct const & seq,
		utility::vector1< core::Size > const & working_res ) const;

	void
	get_seq_and_resnum( std::string const & pdb,
		std::string & seq,
		utility::vector1< int > & resnum,
		utility::vector1< char > & chain,
		utility::vector1< std::string > & segid ) const;

	std::string
	get_silent_seq( std::string const & silent_file ) const;

	std::tuple< utility::vector1< int >, utility::vector1< char >, utility::vector1< std::string > >
	get_silent_resnum( std::string const & silent_file ) const;

	bool
	already_listed_in_obligate_pair( utility::vector1< core::Size > const & new_pair,
		utility::vector1< core::Size > const & obligate_pair ) const;

	bool
	already_listed_in_obligate_pair( utility::vector1< Size > const & new_pair,
		utility::vector1< Size > const & obligate_pair,
		utility::vector1< Size > const & obligate_pair_explicit ) const;

	void
	update_working_obligate_pairs_with_stems(
		utility::vector1< core::pose::rna::BasePair > & working_obligate_pairs,
		utility::vector1< utility::vector1< std::pair< core::Size, core::Size > > > const & working_stems,
		utility::vector1< core::Size > const & working_input_res ) const;

	core::pose::PoseOP
	get_native_pose() { return native_pose_; }

private:

	void
	add_chem_shift_info( core::pose::Pose & not_very_const_pose ) const;


	core::import_pose::options::RNA_DeNovoProtocolOptionsOP options_;
	core::import_pose::RNA_DeNovoParametersOP rna_params_;
	core::chemical::ResidueTypeSetCOP rsd_set_;
	core::pose::PoseOP pose_, native_pose_;

	utility::vector1< std::string > fasta_files_;
	utility::vector1< std::string > input_pdbs_;
	utility::vector1< std::string > input_silent_files_;
	utility::vector1< std::string > sequence_strings_;
	bool minimize_rna_;
	bool minimize_rna_has_been_specified_ = false;
	utility::vector1< std::string > helical_substructs_;
	utility::vector1< std::string > dock_chunks_;
	utility::vector1< std::string > input_initialization_pdbs_;

	protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_grid_;
	// std::string lores_silent_file_;
	RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;

	// std::map< std::string, bool > tag_is_done_;

	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP hires_scorefxn_;

	std::list< core::Real > all_lores_score_final_; // used for filtering.

};

} //movers
} //denovo
} //rna
} //protocols

#endif
