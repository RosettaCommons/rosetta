// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna/RotamerDNAHBondFilter.hh
/// @brief  when passed to a PackerTask, filters rotamers for contact to DNA in an -ex dependent manner
/// @author ashworth

#ifndef INCLUDED_protocols_dna_RotamerDNAHBondFilter_hh
#define INCLUDED_protocols_dna_RotamerDNAHBondFilter_hh

// Unit Headers
#include <protocols/dna/RotamerDNAHBondFilter.fwd.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class RotamerDNAHBondFilter : public core::pack::rotamer_set::RotamerOperation
{
public:
	RotamerDNAHBondFilter(
		core::Real threshold = -0.5,
		bool base_only = true
	);
	~RotamerDNAHBondFilter() override;


	bool
	operator() (
		core::conformation::ResidueOP rotamer,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & scorefxn,
		core::pack::task::ResidueLevelTask const & rtask,
		utility::graph::GraphCOP packer_neighbor_graph,
		core::pack::dunbrack::ChiSetOP chi_set
	) override;

	void report() const;

private:
	core::Real threshold_;
	bool base_only_;
	core::Size nfiltered_;
	core::Size naccepted_;
	core::scoring::hbonds::HBondDatabaseCOP hb_database_;
	core::scoring::hbonds::HBondOptionsOP hbondoptions_;

};

}
}

#endif
