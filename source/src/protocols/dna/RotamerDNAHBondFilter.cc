// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/dna/RotamerDNAHBondFilter.cc
/// @brief  filters rotamers for contact to DNA in an -ex dependent manner
/// @author ashworth

#include <protocols/dna/RotamerDNAHBondFilter.hh>
#include <protocols/dna/util.hh> // close_to_dna

// Package Headers
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#ifdef WIN32
// required for VS2005 build
#include <core/graph/Graph.hh>
#endif

#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.dna.RotamerDNAHBondFilter" );

namespace protocols {
namespace dna {

using namespace core;

RotamerDNAHBondFilter::RotamerDNAHBondFilter(
	core::Real threshold, // = -0.5
	bool base_only // = true
) :
	pack::rotamer_set::RotamerOperation(),
	threshold_( threshold ),
	base_only_( base_only ),
	nfiltered_(0),
	naccepted_(0),
	hb_database_(scoring::hbonds::HBondDatabase::get_database()),
	hbondoptions_(core::scoring::hbonds::HBondOptionsOP( new scoring::hbonds::HBondOptions ) )
{}

RotamerDNAHBondFilter::~RotamerDNAHBondFilter() {}


bool
RotamerDNAHBondFilter::operator() (
	conformation::ResidueOP rotamer,
	pose::Pose const & pose,
	scoring::ScoreFunction const & scorefxn,
	pack::task::ResidueLevelTask const & rtask,
	graph::GraphCOP,
	pack::dunbrack::ChiSetOP chi_set
)
{
	using namespace scoring;
	using namespace hbonds;
	using namespace conformation;

	bool filter(false);

	typedef utility::vector1< Real > Reals;
	Reals const & ex_chi_steps( chi_set->ex_chi_steps );

	Size const nchi( ex_chi_steps.size() );
	if ( ( nchi > 0 && rtask.operate_on_ex1() && ex_chi_steps[1] != 0. ) ||
			( nchi > 1 && rtask.operate_on_ex2() && ex_chi_steps[2] != 0. ) ||
			( nchi > 2 && rtask.operate_on_ex3() && ex_chi_steps[3] != 0. ) ||
			( nchi > 3 && rtask.operate_on_ex4() && ex_chi_steps[4] != 0. ) ) filter = true;

	if ( !filter ) return true;
	++nfiltered_;

	for ( Size pos(1); pos <= pose.total_residue(); ++pos ) {
		Residue const & dnares( pose.residue( pos ) );
		if ( !dnares.is_DNA() ) continue;
		if ( !close_to_dna( *rotamer, dnares, 10*10, base_only_ ) ) continue;

		Real hbE_total(0.0);
		for ( Size ratom_i( rotamer->first_sidechain_atom() ), ratom_end( rotamer->natoms() );
				ratom_i <= ratom_end; ++ratom_i ) {
			Atom const & ratom( rotamer->atom( ratom_i ) );
			Size const datom_start( base_only_ ? dnares.first_sidechain_atom() : 1 );
			for ( Size datom_i( datom_start ), datom_end( dnares.natoms() );
					datom_i <= datom_end; ++datom_i ) {
				Atom const & datom( dnares.atom( datom_i ) );

				// rotamer donor, dna acceptor
				if ( rotamer->atom_type( ratom_i ).is_hydrogen() &&
						rotamer->atom_type( rotamer->atom_base( ratom_i ) ).is_donor() &&
						dnares.atom_type( datom_i ).is_acceptor() ) {
					Real dis2( ratom.xyz().distance_squared( datom.xyz() ) );
					if ( dis2 > MAX_R2 || dis2 < MIN_R2 ) continue;
					HBEvalTuple hbtype(rotamer->atom_base(ratom_i), *rotamer, datom_i, dnares);
					Real hbE;
					hb_energy_deriv_u( *hb_database_, *hbondoptions_, hbtype,
						ratom.xyz(), rotamer->xyz( rotamer->atom_base( ratom_i )), create_don_orientation_vector( *rotamer, ratom_i ),
						datom.xyz(), datom.xyz() /*apl -- acceptor base coordinate goes here, but used only for derivatives */,
						create_acc_orientation_vector( *hbondoptions_, dnares, datom_i ),
						Vector(-1.0, -1.0, -1.0) /* acceptor_base2 atom goes here, this is now wrong*/,
						hbE,
						false /*evaluate_derivative*/, DUMMY_DERIVS );
					hbE_total += hbE * scorefxn.get_weight( hbond_sc );
				} else if (
						// rotamer acceptor, dna donor
						dnares.atom_type( datom_i ).is_hydrogen() &&
						dnares.atom_type( dnares.atom_base( datom_i ) ).is_donor() &&
						rotamer->atom_type( ratom_i ).is_acceptor() ) {
					Real dis2( ratom.xyz().distance_squared( datom.xyz() ) );
					if ( dis2 > MAX_R2 || dis2 < MIN_R2 ) continue;
					HBEvalTuple hbtype( dnares.atom_base(datom_i), dnares, ratom_i, *rotamer );
					Real hbE;
					hb_energy_deriv_u( *hb_database_, *hbondoptions_, hbtype,
						datom.xyz(), dnares.xyz( dnares.atom_base( datom_i ) ), create_don_orientation_vector( dnares, datom_i ),
						ratom.xyz(), ratom.xyz() /* apl -- acceptor-base coordinate goes here, but used only for derivatives */,
						create_acc_orientation_vector( *hbondoptions_, *rotamer, ratom_i ),
						Vector(-1.0, -1.0, -1.0)  /* apl -- acceptor_base2 atom goes here, this is now wrong */,
						hbE,
						false /*evaluate_derivative*/, DUMMY_DERIVS );
					hbE_total += hbE * scorefxn.get_weight( hbond_sc );

				} else continue;

				if ( hbE_total < threshold_ ) {
					//     TR << "rotamer type " << rotamer->name() << " passed hbE threshold" << std::endl;
					++naccepted_;
					return true;
				}
			}
		}
	}
	return false;
}

void
RotamerDNAHBondFilter::report() const
{
	TR << nfiltered_ << " rotamers were filtered, " << naccepted_ << " accepted." << std::endl;
}

}
}

