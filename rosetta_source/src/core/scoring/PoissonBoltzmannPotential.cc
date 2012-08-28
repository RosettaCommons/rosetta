// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/PoissonBoltzmannPotential.cc
/// @brief  PoissonBoltzmann potential class implementation
/// @author Yifan Song (yfsong@uw.edu)

// basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Unit Headers
#include <core/scoring/PoissonBoltzmannPotential.hh>

// Package Headers
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/chemical/ResidueType.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Numeric Headers
// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// utility
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

// option key includes

static basic::Tracer TR("core.scoring.PoissonBoltzmannPotential");

namespace core {
namespace scoring {

// @brief Auto-generated virtual destructor
PoissonBoltzmannPotential::~PoissonBoltzmannPotential() {}

typedef PoissonBoltzmannPotential PB;

PB::PoissonBoltzmannPotential()
{
	potential_is_loaded_ = false;
}

void
dump_pqr(
		 core::pose::Pose const & pose,
		 std::ostream & out,
		 std::string const &,
		 utility::vector1 <Size> const & zero_charge_chains
		 ) {
	Size const nres( pose.total_residue() );

	Size number(0);

	static std::string const chains( " ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890" );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		for ( Size j=1; j<= rsd.natoms(); ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

			//skip outputing virtual atom unless specified.
			//fixed so that the last atom in atom type set can be something other than a virtual atom --steven combs
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::output_virtual ]() &&
				rsd.atom_type(j).is_virtual() ) continue;

			++number;
			runtime_assert( rsd.chain() < chains.size() ); // silly restriction
			char const chain( chains[ rsd.chain() ] );
			if (find(zero_charge_chains.begin(), zero_charge_chains.end(), rsd.chain()) == zero_charge_chains.end()) {
				using namespace ObjexxFCL::fmt;
				out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(8,3,pose.residue_type(i).atomic_charge(j)) <<
				F(8,3,rsd.atom_type(j).lj_radius()) << '\n';
			}
			else {
				using namespace ObjexxFCL::fmt;
				out << "ATOM  " << I(5,number) << ' ' << rsd.atom_name(j) << ' ' <<
				rsd.name3() << ' ' << chain << I(4,rsd.seqpos() ) << "    " <<
				F(8,3,atom.xyz()(1)) <<
				F(8,3,atom.xyz()(2)) <<
				F(8,3,atom.xyz()(3)) <<
				F(8,3,0.0) <<
				F(8,3,rsd.atom_type(j).lj_radius()) << '\n';
			}
		}
	}
}

void
PB::read_APBS_potential(std::string const & apbs_potential_fn)
{
	utility::io::izstream p_stream( apbs_potential_fn );

	std::string line;
	while (p_stream) {
		//	# Comments
		char first_char = p_stream.peek();
		if ( first_char == '#' ) {
			p_stream.getline( line );
			continue;
		}

		//object 1 class gridpositions counts nx ny nz
		std::string buff;
		for (Size i=0;i<5;++i) p_stream >> buff;
		p_stream >> n_grid_[0] >> n_grid_[1] >> n_grid_[2];

		//origin xmin ymin zmin
		p_stream >> buff;
		p_stream >> lower_bound_[0];
		p_stream >> lower_bound_[1];
		p_stream >> lower_bound_[2];

		//delta hx 0.0 0.0
		//delta 0.0 hy 0.0
		//delta 0.0 0.0 hz

		p_stream >> buff >> grid_spacing_[0] >> buff >> buff;
		p_stream >> buff >> buff >> grid_spacing_[1] >> buff;
		p_stream >> buff >> buff >> buff >> grid_spacing_[2];
		i2c_ = numeric::xyzMatrix <core::Real>::rows(
													 grid_spacing_[0],0.,0.,
													 0.,grid_spacing_[1],0.,
													 0.,0.,grid_spacing_[2]);
		c2i_ = numeric::xyzMatrix <core::Real>::rows(
													 1./grid_spacing_[0],0.,0.,
													 0.,1./grid_spacing_[1],0.,
													 0.,0.,1./grid_spacing_[2]);

		//object 2 class gridconnections counts nx ny nz
		//object 3 class array type double rank 0 items n data follows
		p_stream.getline( line );
		p_stream.getline( line );
		p_stream.getline( line );

		//u(0,0,0) u(0,0,1) u(0,0,2)
		potential_.dimension(n_grid_[0],n_grid_[1],n_grid_[2]);

		for (int i=1; i<=potential_.u1(); i++)
			for (int j=1; j<=potential_.u2(); j++)
				for (int k=1; k<=potential_.u3(); k++) {
					p_stream >> buff;
					potential_(i,j,k) = atof(buff.c_str());
					if (potential_(i,j,k) > basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_cap ]) {
						potential_(i,j,k) = basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_cap ];
					}
					if (potential_(i,j,k) < -basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_cap ]) {
						potential_(i,j,k) = -basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_cap ];
					}
				}
		break;
	}


	idx2cart(n_grid_, upper_bound_);

	potential_is_loaded_ = true;
	//attribute "dep" string "positions"
	//object "regular positions regular connections" class field
	//component "positions" value 1
	//component "connections" value 2
	//component "data" value 3
}

core::Real
PB::get_potential(ObjexxFCL::FArray3D< core::Real > const & potential, numeric::xyzVector<core::Real> const & cartX) const {
	if (out_of_bounds(cartX)) return 0.0;

	numeric::xyzVector<core::Real> idxX;
	cart2idx(cartX, idxX);
	return core::scoring::electron_density::interp_linear(potential, idxX);
}


void
PB::eval_PB_energy_residue(
						   core::conformation::Residue const & rsd,
						   Real & PB_energy_residue,
						   Real & PB_energy_backbone,
						   Real & PB_energy_sidechain,
						   Real const & PB_burial_weight
						   ) const {
	PB_energy_residue = 0.0;
	PB_energy_backbone = 0.0;
	PB_energy_sidechain = 0.0;
	if (basic::options::option[ basic::options::OptionKeys::corrections::score::PB_score_residue_range ].user()) {
		if (rsd.seqpos() < (Size) basic::options::option[ basic::options::OptionKeys::corrections::score::PB_score_residue_range ][1]) return;
		if (rsd.seqpos() > (Size) basic::options::option[ basic::options::OptionKeys::corrections::score::PB_score_residue_range ][2])	return;
	}

	for ( Size iatom=1; iatom<=rsd.natoms(); ++iatom ) {
		if (rsd.atomic_charge(iatom) > -1e-6 && rsd.atomic_charge(iatom) < 1e-6) continue;
		core::Real iatom_potential = get_potential(potential_, rsd.xyz(iatom));

		core::Real atom_energy = rsd.atomic_charge(iatom)*iatom_potential * PB_burial_weight;

		if (rsd.atom_is_backbone(iatom)) {
			PB_energy_backbone += atom_energy;
		}
		else {
			PB_energy_sidechain += atom_energy;
		}
		PB_energy_residue += atom_energy;
	}

}

PoissonBoltzmannPotential & get_PB_potential() {
	static PoissonBoltzmannPotential potential;
	if (!potential.isLoaded()) {
		if (!basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_file ].user()) {
			TR.Warning << "[ Warning ] No PB potential map specified." << std::endl;
		} else {
			potential.read_APBS_potential(basic::options::option[ basic::options::OptionKeys::corrections::score::PB_potential_file ]);
		}
	}
	return potential;
}

}
}
