// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/NucleotideMutation.hh
/// @author Christoffer Norn (ch.norn@gmail.com)

#ifndef INCLUDED_protocols_protein_interface_design_movers_NucleotideMutation_hh
#define INCLUDED_protocols_protein_interface_design_movers_NucleotideMutation_hh
#include <protocols/protein_interface_design/movers/NucleotideMutation.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief designs alanine residues in place of the residue identities at the interface. Retains interface glycines and prolines.
class NucleotideMutation : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;
public:
	NucleotideMutation();
	//std::string codon2aa( std::string const codon ) const;
	//std::string aa2randomCodon( char const aa ) const;
	//std::vector<std::string> aa_2_nt( std::string const aa ) const;
	void add_nt_seq_to_pose( core::pose::Pose & pose ); // get the segment names for those segments that are constant in this splice function
	void apply( Pose & pose );
	virtual std::string get_name() const;
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new NucleotideMutation ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
	virtual ~NucleotideMutation();
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP tf );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP sfxn );

	std::string init_sequence() const{ return init_sequence_; }
	void init_sequence( std::string const & s ){ init_sequence_ = s; };

	bool continue_if_silent() const{ return continue_if_silent_; }
	void continue_if_silent( bool const c ){ continue_if_silent_ = c; };

	//core::pose::PoseOP reference_pose() const{ return reference_pose_; }
	//void reference_pose( core::pose::PoseOP const p ){ reference_pose_ = p; };

	bool cache_task() const;
	void cache_task( bool cache );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string init_sequence_;
	bool continue_if_silent_;
	core::pose::PoseOP reference_pose_;
	bool cache_task_;
	core::pack::task::PackerTaskCOP task_;
};


} // movers
} // protein_interface_design
} // protocols


#endif /*INCLUDED_protocols_protein_interface_design_movers_NucleotideMutation_HH*/
