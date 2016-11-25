// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/seeded_abinitio/CloseFold.cc
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_CloseFold_hh
#define INCLUDED_protocols_seeded_abinitio_CloseFold_hh


#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/string_util.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/Frame.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
namespace seeded_abinitio {

class CloseFold : public protocols::moves::Mover
{

public:
	CloseFold();

	~CloseFold() override;
	void apply( core::pose::Pose & pose ) override;
	// XRW TEMP  std::string get_name() const override;
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

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

	core::Size trials();

	void set_trials( core::Size trials_quick_ccd );

	void initialize_fragments();

	//   core::scoring::ScoreFunction cen_scorefxn();
	//   core::scoring::ScoreFunction fa_scorefxn();

private: /// data

	/// fragment set used for ccd
	core::fragment::FragSetOP fragments_;

	///position-indexable list of Frames, populating library
	///void initialize_library();

	///index-based access to the data contained in the FragSet
	boost::unordered_map<core::Size, core::fragment::Frame> library_;
	//FrameMap library_;

	void fast_loopclose( core::pose::Pose &pose, protocols::loops::LoopsOP const loops, bool kic );

	void quick_closure( core::pose::Pose &pose, protocols::loops::LoopsOP const loops );

	core::pose::PoseOP template_pdb_;

	std::string secstructure_;

	utility::vector1< core::Size > chains_;

	//gather all loops for specified chains or just the ones that have a cutpoint
	bool use_cutpoints_;

	///residues specifying the seeds
	utility::vector1< std::pair < std::string,std::string > > seed_vector_;

	//protocols::loops::Loops seeds_;
	//utility::vector1< core::Size > seed_vector_;

	protocols::loops::LoopsOP loops_;

	///add cutpoint variants for closure
	bool chainbreakweights_;

	///option for quick ccd protocol, how many attempts
	core::Size trials_;

	///options for fast_ccd protocol
	bool idealize_;

	///should kinematic loop mover be used after the fast ccd
	bool kic_;
	bool ccd_;
	core::scoring::ScoreFunctionOP cen_scorefxn_;

	core::scoring::ScoreFunctionOP fa_scorefxn_;
};
}//end seeded_abinitio
}//end protocols

#endif
