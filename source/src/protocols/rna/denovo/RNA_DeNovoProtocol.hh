// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_DeNovo_Protocol.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH
#define INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH

#include <protocols/moves/Mover.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.fwd.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>

//// C++ headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <core/import_pose/RNA_DeNovoParameters.fwd.hh> // AUTO IWYU For RNA_DeNovoParametersCOP
#include <core/io/silent/SilentStruct.fwd.hh> // AUTO IWYU For SilentStruct
#include <map> // AUTO IWYU For map

namespace protocols {
namespace rna {
namespace denovo {

/// @brief The RNA de novo structure modeling protocol
class RNA_DeNovoProtocol: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_DeNovoProtocol( core::import_pose::options::RNA_DeNovoProtocolOptionsCOP options = nullptr,
		core::import_pose::RNA_DeNovoParametersCOP params = nullptr);

	~RNA_DeNovoProtocol() override;

	/// @brief Clone this object
	protocols::moves::MoverOP clone() const override;
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


	/// @brief Apply the RNA denovo modeling protocol to the input pose
	void apply( core::pose::Pose & pose ) override;

	static void register_options();

	void show(std::ostream & output=std::cout) const override;

	void
	output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag,
		bool const score_only = false ) const;

	void
	align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const;

	void
	set_refine_pose_list( utility::vector1<core::pose::PoseOP> const & setting ) { refine_pose_list_ = setting; };

	core::import_pose::options::RNA_DeNovoProtocolOptionsCOP options() const { return options_; }

	RNA_FragmentMonteCarloCOP rna_fragment_monte_carlo() const { return rna_fragment_monte_carlo_; }

private:

	void
	initialize_lores_silent_file();

	void
	initialize_constraints( core::pose::Pose & pose );

	void
	initialize_scorefxn( core::pose::Pose & pose );

	void
	initialize_tag_is_done();

	void
	output_silent_struct( core::io::silent::SilentStruct & s,
		core::io::silent::SilentFileData & silent_file_data,
		std::string const & silent_file,
		core::pose::Pose & pose,
		std::string const & out_file_tag,
		bool const score_only = false ) const;

	void
	check_for_loop_modeling_case( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map, core::pose::Pose const & pose ) const;

	void
	calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose, std::string const & out_file_tag ) const;

	void
	add_chem_shift_info(core::io::silent::SilentStruct & silent_struct, core::pose::Pose const & const_pose) const;

private:

	core::import_pose::options::RNA_DeNovoProtocolOptionsCOP options_;
	core::import_pose::RNA_DeNovoParametersCOP rna_params_;
	protocols::stepwise::modeler::rna::checker::RNA_VDW_BinCheckerOP vdw_grid_;

	std::string lores_silent_file_;

	RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;

	std::map< std::string, bool > tag_is_done_;

	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP hires_scorefxn_;

	utility::vector1<core::pose::PoseOP> refine_pose_list_;
	core::Size refine_pose_id_ = 1;

}; // class RNA_DeNovoProtocol

std::ostream &operator<< ( std::ostream &os, RNA_DeNovoProtocol const &mover );

// AMW: for some reason, no OP/COP defined/no fwd header...
typedef utility::pointer::shared_ptr< RNA_DeNovoProtocol > RNA_DeNovoProtocolOP;
typedef utility::pointer::shared_ptr< RNA_DeNovoProtocol const > RNA_DeNovoProtocolCOP;

} //denovo
} //rna
} //protocols

#endif
