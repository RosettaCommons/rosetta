// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaInterfacePacker.hh
/// @author ashworth
/// @brief Intended to perform anything one would want to do in a DNA interface that uses a single incarnation of the RotamerSets/InteractionGraph combo.
/// @details For basic packing/designing, the more basic PackRotamersMover can be used instead, provided it receives the appropriate TaskFactory/TaskOperations. This derived class, however, takes advantage of the reusability of packer data to accomplish some higher-level functions, such as rapid estimations of multi-state specificity, reversions, and mutational scanning.

#ifndef INCLUDED_protocols_dna_DnaInterfacePacker_hh
#define INCLUDED_protocols_dna_DnaInterfacePacker_hh

#include <protocols/dna/DnaInterfacePacker.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <protocols/dna/typedefs.hh>
#include <protocols/dna/DnaChains.fwd.hh>
#include <protocols/dna/PDBOutput.fwd.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <map>
#include <string>

#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class ResTypeSequence_lt;

typedef std::map< ResTypeSequence, core::Real, ResTypeSequence_lt > SequenceScores;

class ResTypeSequence_lt {
public:
	bool operator() ( ResTypeSequence const & a, ResTypeSequence const & b ) const;
};

class DnaInterfacePacker : public protocols::simple_moves::PackRotamersMover {

public:
	typedef utility::tag::TagCOP TagCOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	DnaInterfacePacker();

	DnaInterfacePacker(
		ScoreFunctionOP,
		bool minimize = false,
		std::string filename_root = "dnapacker"
	);

	~DnaInterfacePacker() override; // important for properly "releasing" owning pointer data

	moves::MoverOP fresh_instance() const override;
	moves::MoverOP clone() const override;
	// XRW TEMP  std::string get_name() const override;

	void apply( Pose & ) override;

	bool initialized() const;

	/// @brief parse XML (specifically in the context of the parser/scripting scheme)
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;

	core::Real unbound_score( Pose const &, bool output_pdb = false, std::string pdbname = "" );

	std::pair< core::Real, core::Real > measure_specificity( Pose & );
	std::pair< SequenceScores, SequenceScores > measure_bp_specificities( Pose & );
	std::pair< SequenceScores, SequenceScores >
	measure_specificities( Pose &, ResTypeSequences const & );

	void reversion_scan(
		Pose &,
		core::Real bound_score = 0.,
		core::Real binding_score = 0.,
		std::pair< core::Real, core::Real > specificities = std::make_pair(0.,0.) );

	// setters, getters, and book-keeping methods
	void reference_pose( Pose const & );
	PoseCOP reference_pose() const;
	void targeted_dna( DnaDesignDefOPs const & );
	DnaDesignDefOPs const & targeted_dna() const;
	void pdboutput( PDBOutputOP );
	void set_filename_root( std::string const & name ) { filename_root_ = name; }
	std::string dna_seq_tag( Pose const &, ResTypeSequence const & ) const;
	ResTypeSequence get_targeted_sequence( Pose const & ) const;
	ResTypeSequence current_working_sequence( Pose const & ) const;
	std::string current_dna_design_string( Pose const & ) const;
	std::string allowed_types() const;
	void clear_initialization();
	std::string pdbname() { return pdbname_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private: // methods
	void standard_packing( Pose & );
	void post_packing( Pose &, ResTypeSequence const &, Size );
	void protein_scan( Pose & );
	void init_standard( Pose & );
	void make_dna_sequence_combinations( core::chemical::ResidueTypeSet const & );
	void add_complementary_sequence( ResTypeSequence &, core::chemical::ResidueTypeSet const & );
	core::Real calculate_specificity( Pose const &, ResTypeSequence const &, SequenceScores const & );

private: // data
	PoseCOP reference_pose_;
	DnaDesignDefOPs targeted_dna_;
	DnaChainsOP dna_chains_;
	ResTypeSequences dna_sequences_;
	bool minimize_;
	std::string filename_root_;
	bool binding_E_;
	bool probe_specificity_;
	bool reversion_scan_;
	bool protein_scan_;
	std::string allowed_types_;
	bool base_only_;
	bool include_dna_potentials_in_specificity_calculations_;
	core::Size num_repacks_;
	core::Size specificity_repacks_;
	core::optimization::MinimizerOptionsOP minimize_options_;
	core::kinematics::MoveMapOP min_movemap_;
	PDBOutputOP pdboutput_;
	utility::vector1< core::chemical::ResidueTypeCOP > reference_residue_types_;
	bool initialization_state_;
	std::string pdbname_;
};

} // namespace dna
} // namespace protocols

#endif
