// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/DefineMovableLoops.cc
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_DefineMovableLoops_hh
#define INCLUDED_protocols_seeded_abinitio_DefineMovableLoops_hh

#include <core/pose/Pose.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/string_util.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace seeded_abinitio {

class DefineMovableLoops : public protocols::moves::Mover
{
public:
	DefineMovableLoops();

	~DefineMovableLoops() override;
	void apply( core::pose::Pose & pose ) override;
	// XRW TEMP  std::string get_name() const override;
	protocols::moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new DefineMovableLoops( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new DefineMovableLoops ); }

	void parse_my_tag(  utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private: ///functions

	protocols::loops::LoopsOP find_loops (  core::pose::Pose & pose,
		std::string secstr,
		core::Size offset,
		protocols::loops::Loops seeds );

	bool is_cut( utility::vector1<Size> & cut_points, Size residue);

	bool chainbreakweights();

	void add_chainbreakweights( bool acbw );

	bool use_cutpoints();

	void use_cutpoints( bool uc );


private: /// data

	core::pose::PoseOP template_pdb_;

	std::string secstructure_;

	utility::vector1< core::Size > chains_;

	//gather all loops for specified chains or just the ones that have a cutpoint
	bool use_cutpoints_;

	///residues specifying the seeds
	utility::vector1< std::pair < std::string,std::string > > seed_vector_;

	protocols::loops::LoopsOP loops_;

	///add cutpoint variants for closure
	bool chainbreakweights_;

};
}//end seeded_abinitio
}//end protocols

#endif
