// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/SameChiBinComboGrouper.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/SameChiBinComboGrouper.hh>

// Package headers
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>
#include <protocols/match/Hit.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

// Utility headers
#include <utility/exit.hh>

#include <utility/OrderedTuple.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

SameChiBinComboGrouper::SameChiBinComboGrouper() : n_geometric_constraints_( 0 ) {}
SameChiBinComboGrouper::SameChiBinComboGrouper( Size ncst ) : n_geometric_constraints_( ncst ) {}

SameChiBinComboGrouper::~SameChiBinComboGrouper() = default;

SameChiBinComboGrouper::Size
SameChiBinComboGrouper::assign_group_for_match(
	match const & m
)
{
	return assign_group_for_match( match_dspos1( m, 1 ) );
}

SameChiBinComboGrouper::Size
SameChiBinComboGrouper::assign_group_for_match(
	match_dspos1 const & m
)
{
	using namespace core::scoring;
	using namespace core::pack::dunbrack;
	using namespace core::pack::rotamers;

	runtime_assert( m.upstream_hits.size() == n_geometric_constraints_ );

	SingleResidueRotamerLibraryFactory const & rotlib( *SingleResidueRotamerLibraryFactory::get_instance() );

	utility::vector1< Size > rot_vector( n_geometric_constraints_ * 3, 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		rot_vector[ ii ] = m.upstream_hits[ ii ].scaffold_build_id();
		core::conformation::ResidueCOP iires = hit_cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] ) );

		/*
		for ( Size jj = 1; jj <= iires->nchi(); ++jj ) {
		core::Vector p1, p2, p3, p4;
		p1 = iires->xyz( iires->chi_atoms( jj )[ 1 ] );
		p2 = iires->xyz( iires->chi_atoms( jj )[ 2 ] );
		p3 = iires->xyz( iires->chi_atoms( jj )[ 3 ] );
		p4 = iires->xyz( iires->chi_atoms( jj )[ 4 ] );

		std::cout << "res chi reported: chi " << jj << " " << iires->chi( jj ) << " real: "
		<< numeric::constants::d::radians_to_degrees * numeric::dihedral_radians(
		p1, p2, p3, p4 ) << std::endl;

		}*/

		rot_vector[ n_geometric_constraints_ + ii ] = iires->aa();

		SingleResidueRotamerLibraryCOP srrl( rotlib.get( iires->type() ) );
		if ( ! srrl ) {
			/// ?!?! What do we without a library?
			rot_vector[ 2*n_geometric_constraints_ + ii ] = 1;
		} else {
			SingleResidueDunbrackLibraryCOP srdl(
				utility::pointer::dynamic_pointer_cast< SingleResidueDunbrackLibrary const > ( srrl ));
			if ( srdl ) {
				RotVector rotvect;
				srdl->get_rotamer_from_chi( iires->chi(), rotvect );
				rot_vector[ 2*n_geometric_constraints_ + ii ] = srdl->rotwell_2_rotno( rotvect );
			} else {
				/// ?!?! What do we do with a non-dunbrack library?
				rot_vector[ 2*n_geometric_constraints_ + ii ] = 1;
			}
		}
	}

	ChiBinComboCountMap::const_iterator iter = chibin_combo_indexer_.find( rot_vector );
	if ( iter == chibin_combo_indexer_.end() ) {
		Size next_index = chibin_combo_indexer_.size() + 1;
		chibin_combo_indexer_[ rot_vector ] = next_index;
		return next_index;
	} else {
		return iter->second;
	}
}

void
SameChiBinComboGrouper::reset()
{
	chibin_combo_indexer_.clear();
}

void
SameChiBinComboGrouper::set_n_geometric_constraints( Size n_csts )
{
	n_geometric_constraints_ = n_csts;
}

void
SameChiBinComboGrouper::set_hit_cacher( UpstreamHitCacherOP cacher )
{
	hit_cacher_ = cacher;
}


}
}
}
