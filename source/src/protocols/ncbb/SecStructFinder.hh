// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_ncbb_SecStructFinder_hh
#define INCLUDED_protocols_ncbb_SecStructFinder_hh

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/ResidueType.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/ncbb/SecStructFinder.fwd.hh>
#include <protocols/ncbb/SecStructFinderCreator.fwd.hh>

#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <sstream>
#include <cmath>

namespace protocols {
namespace ncbb {

class SecStructFinder : public moves::Mover
{
public:
	SecStructFinder();
	SecStructFinder(
		std::string residue,
		core::Size min_length = 5,
		core::Size max_length = 5,
		core::Real bin_size = 10,
		core::Real dissimilarity = 10,
		core::Real dihedral_min = -180,
		core::Real dihedral_max = 180,
		core::Real dump_threshold = -10000,
		std::string dihedral_pattern = "A",
		std::string alpha_beta_pattern = "A",
		bool min_everything = false,
		bool cart = false,
		bool constrain = false
	);

	~SecStructFinder() = default;
	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "SecStructFinder"; }
	protocols::moves::MoverOP fresh_instance() const override { return SecStructFinderOP( new SecStructFinder ); }
	protocols::moves::MoverOP clone() const override;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	std::string residue_;
	core::Size min_length_;
	core::Size max_length_;
	core::Real bin_size_;
	core::Real dump_threshold_;
	core::Real dihedral_min_;
	core::Real dihedral_max_;
	std::string dihedral_pattern_;
	std::string alpha_beta_pattern_;
	bool min_everything_;
	bool cart_;
	bool constrain_;
	core::Real dissimilarity_;

	std::string alpha_to_beta( std::string alpha );
	std::string expand_pattern_to_fit( std::string pattern, core::Size length );
	bool uniq_refers_to_beta ( char uniq );
	void initialize_rtype_vector( utility::vector1< core::chemical::ResidueType > & restypes );
	std::string make_filename ( core::Size number_dihedrals, utility::vector1< core::Real > dihedrals );
	bool too_similar( core::Size i, core::Size j, utility::vector1< core::Real > dihedrals );
	void show_current_dihedrals( core::Size number_dihedral_sets, utility::vector1< char > uniqs, utility::vector1< core::Real > dihedrals );
};

}
}

#endif
