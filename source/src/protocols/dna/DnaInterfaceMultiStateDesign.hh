// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DnaInterfaceMultiStateDesign.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaInterfaceMultiStateDesign_hh
#define INCLUDED_protocols_dna_DnaInterfaceMultiStateDesign_hh

#include <protocols/dna/DnaInterfaceMultiStateDesign.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/dna/DnaChains.fwd.hh>
#include <protocols/dna/DnaDesignDef.fwd.hh>

#include <protocols/multistate_design/MultiStatePacker.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
// AUTO-REMOVED #include <protocols/genetic_algorithm/GeneticAlgorithm.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <protocols/genetic_algorithm/GeneticAlgorithm.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace dna {

///@brief wraps DNA-interface specific considerations around the general multistate design / genetic algorithm framework
// could also derive from DnaInterfacePacker(?)
class DnaInterfaceMultiStateDesign : public protocols::simple_moves::PackRotamersMover {
public:
	typedef multistate_design::PosType PosType;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef genetic_algorithm::GeneticAlgorithmOP GeneticAlgorithmOP;

public:
	DnaInterfaceMultiStateDesign();
	virtual ~DnaInterfaceMultiStateDesign();
	virtual void apply( Pose & );
	virtual std::string get_name() const;
	void copy_dna_chains( DnaChainsOP );
	void copy_targeted_dna( DnaDesignDefOPs const & );
	void output_results( Pose & );

	///@brief parse XML (specifically in the context of the parser/scripting scheme)
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & );
	///@brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP fresh_instance() const;
	///@brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP clone() const;

private:
	void initialize( Pose & );

	using protocols::simple_moves::PackRotamersMover::run;

	void run();
	void add_dna_states( Pose const & , PackerTaskCOP );

private:
	GeneticAlgorithmOP gen_alg_;
	// direct use of MultiStatePacker is only for outputting results
	// (GeneticAlgorithm also holds pointer to it and uses it heavily)
	multistate_design::MultiStatePackerOP multistate_packer_;
	DnaChainsOP dna_chains_;
	DnaDesignDefOPs targeted_dna_;
	// option flags/parameters: constructor defaults to command line options
	// parse_my_tag method may change them
	core::Size generations_, pop_size_, num_packs_, pop_from_ss_, numresults_;
	core::Real fraction_by_recombination_, mutate_rate_, boltz_temp_, anchor_offset_;
	// checkpointing options
	std::string checkpoint_prefix_;
	core::Size checkpoint_interval_;
	bool checkpoint_gz_, checkpoint_rename_;

};

} // namespace dna
} // namespace protocols

#endif
