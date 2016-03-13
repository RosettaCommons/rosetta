// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

using namespace core;
using namespace utility;

namespace protocols {
namespace ncbb {

class SecStructFinder : public moves::Mover
{
public:
	SecStructFinder();
	SecStructFinder(
		std::string residue,
		Size min_length = 5,
		Size max_length = 5,
		Real bin_size = 10,
		Real dissimilarity = 10,
		Real dihedral_min = -180,
		Real dihedral_max = 180,
		Real dump_threshold = -10000,
		std::string dihedral_pattern = "A",
		std::string alpha_beta_pattern = "A",
		bool min_everything = false,
		bool cart = false,
		bool constrain = false
	);

	virtual ~SecStructFinder(){}
	virtual void apply( Pose & pose );
	virtual std::string get_name() const { return "SecStructFinder"; }
	protocols::moves::MoverOP fresh_instance() const { return SecStructFinderOP( new SecStructFinder ); }
	protocols::moves::MoverOP clone() const;
	void parse_my_tag( utility::tag::TagCOP, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );

private:
	core::scoring::ScoreFunctionOP score_fxn_;
	std::string residue_;
	Size min_length_;
	Size max_length_;
	Real bin_size_;
	Real dump_threshold_;
	Real dihedral_min_;
	Real dihedral_max_;
	std::string dihedral_pattern_;
	std::string alpha_beta_pattern_;
	bool min_everything_;
	bool cart_;
	bool constrain_;
	Real dissimilarity_;

	std::string alpha_to_beta( std::string alpha );
	std::string expand_pattern_to_fit( std::string pattern, Size length );
	bool uniq_refers_to_beta ( char uniq );
	void initialize_rtype_vector( utility::vector1< core::chemical::ResidueType > & restypes );
	std::string make_filename ( Size number_dihedrals, utility::vector1< Real > dihedrals );
	bool too_similar( Size i, Size j, utility::vector1< Real > dihedrals );
	void show_current_dihedrals( Size number_dihedral_sets, utility::vector1< char > uniqs, utility::vector1< Real > dihedrals );
};

}
}

#endif
