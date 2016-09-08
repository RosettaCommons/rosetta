// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/SuckerEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/methods/SuckerEnergy.hh>
#include <core/scoring/methods/SuckerEnergyCreator.hh>
#include <core/scoring/EnergyMap.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <utility/io/izstream.hh>
#include <sstream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


#define DIST_TH 6.0

namespace core {
namespace scoring {
namespace methods {

SuckerEnergy::~SuckerEnergy() {}


/// @details This must return a fresh instance of the SuckerEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SuckerEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SuckerEnergy );
}

ScoreTypes
SuckerEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( suck );
	return sts;
}


SuckerEnergy::SuckerEnergy() :
	parent( methods::EnergyMethodCreatorOP( new SuckerEnergyCreator ) )
{
	utility::io::izstream izin;
	basic::database::open( izin, basic::options::option[ basic::options::OptionKeys::in::file::sucker_params ]() );

	utility::vector1< utility::vector1< Real > > points;

	char s[999];
	izin.getline(s,999); // read in and ignore header line
	while ( true ) {
		izin.getline(s,999);
		// std::cerr << "string: '" << s << "'" << std::endl;
		std::istringstream iss(s);
		core::Real tmp;
		utility::vector1< Real > point;
		if ( !(iss >> tmp ) ) break; // if no line, done
		point.push_back(tmp);
		iss >> tmp;
		point.push_back(tmp);
		if ( iss >> tmp ) point.push_back(tmp);
		points.push_back(point);
	}

	debug_assert( points.size() >= 2 );

	// std::cerr << "printing points" << std::endl;
	// for( Size i = 1; i <= points.size(); ++i ) {
	//  std::cerr << "point: ";
	//  for( Size j = 1; j <= points[i].size(); ++j ) {
	//   std::cerr << points[i][j] << ",";
	//  }
	//  std::cerr << std::endl;
	// }

	using namespace numeric::interpolation::spline;
	SplineGenerator sg(
		points[      1      ][1], points[      1      ][2], points[      1      ][3],
		points[points.size()][1], points[points.size()][2], points[points.size()][3] );

	for ( Size i = 2; i < points.size(); ++i ) {
		if ( points[i].size() == 2 ) sg.add_known_value( points[i][1], points[i][2] );
		if ( points[i].size() == 3 ) sg.add_known_value( points[i][1], points[i][2], points[i][3] );
	}

	interp_ = sg.get_interpolator();

	// core::Real x,y,dy;
	// for( x = 0; x <= 6.0; x += 0.01 ) {
	//  interp_->interpolate(x,y,dy);
	//  std::cerr << "SUCKERE " << x << " " << y << " " << dy << std::endl;
	// }
}


/// clone
EnergyMethodOP
SuckerEnergy::clone() const
{
	return EnergyMethodOP( new SuckerEnergy() );
}


inline bool count_atom( int const & atype ) {
	// return true;
	return (  (  1 <= atype && atype <=  6 ) // just carbons for now...
		|| ( 18 <= atype && atype <= 19 )
	);
}
/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
SuckerEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace core;
	using namespace conformation;

	if ( rsd1.seqpos() == rsd2.seqpos() ) return;

	std::string name1 = rsd1.name();
	std::string name2 = rsd2.name();

	if ( "SUCK" == name1 && "SUCK" == name2 ) return;
	if ( "SUCK" != name1 && "SUCK" != name2 ) return;

	Residue const & rsd( ( "SUCK" != name1 ) ? rsd1 : rsd2 );
	Residue const & sck( ( "SUCK" == name1 ) ? rsd1 : rsd2 );

	numeric::xyzVector<Real> suck_xyz( sck.xyz(1) );
	for ( int i = 1; i <= (int)rsd.nheavyatoms(); ++i ) {
		conformation::Atom const & atom( rsd.atom(i) );
		if ( !count_atom( atom.type() ) ) continue;
		numeric::xyzVector<Real> atom_xyz( rsd.xyz(i) );
		Real const d = atom_xyz.distance( suck_xyz );
		if ( d > DIST_TH ) continue;

		Real f, df;
		interp_->interpolate( d, f, df );
		// std::cerr << rsd.seqpos() << " " << i << " " << d << " " << f << " " << df << " SUCKER" << std::endl;
		// std::cerr << "SCOREATOMS " << rsd.seqpos() << " " << i << " " << sck.seqpos() << " " << 1 << " " << d << std::endl;
		emap[ suck ] += f;
	}
}

/////////////////////////////////////////////////////////////////////////////
/// probably move this logic to PairEPotential ??
void
SuckerEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using namespace numeric;

	// id is a SUCK atom
	if ( "SUCK" == pose.residue(id.rsd()).name() ) { // loop over atom neighbors and suck
		// conformation::Residue const & sck( pose.residue(id.rsd()) );
		debug_assert( 1 <= id.atomno() && id.atomno() <= 3 );
		if ( id.atomno() > 1 ) return;
		numeric::xyzVector<Real> suck_xyz( pose.xyz(id) );

		for ( int ir = 1; ir <= (int)pose.size(); ++ir ) {
			conformation::Residue const & res( pose.residue(ir) );
			if ( "SUCK" == res.name() ) continue;

			for ( int ia = 1; ia <= (int)res.nheavyatoms(); ++ia ) {
				conformation::Atom const & atom( res.atom( ia ) );
				if ( !count_atom( atom.type() ) ) continue;

				xyzVector<Real> const & atom_xyz( atom.xyz() );

				Real const d = suck_xyz.distance(atom_xyz);
				if ( d > DIST_TH ) continue;

				// std::cerr << "DERIVATOMS " << res.seqpos() << " " << ia << " " << sck.seqpos() << " " << 1 << " " << d << " LPA " << id.rsd() << " " << id.atomno() << std::endl;

				// compute suck F1, F2
				numeric::xyzVector<Real> const f1( suck_xyz.cross( atom_xyz ));
				numeric::xyzVector<Real> const f2( suck_xyz -  atom_xyz );
				Real const dist( f2.length() );
				Real deriv, dummy;
				interp_->interpolate( dist, dummy, deriv );
				// std::cerr << "sck_atm " << deriv << std::endl;
				F1 += ( deriv / dist ) * weights[ suck ] * f1;
				F2 += ( deriv / dist ) * weights[ suck ] * f2;
			}
		}
	} else { // loop over suck residues and suck

		// heavy atoms only
		if ( id.atomno() > pose.residue(id.rsd()).nheavyatoms() ) return;
		conformation::Atom const & atom( pose.residue(id.rsd()).atom(id.atomno()) );
		if ( !count_atom( atom.type() ) ) return;

		numeric::xyzVector<Real> atom_xyz( atom.xyz() );
		// while( "SUCK" == pose.residue(i).name() ) {
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( "SUCK" != pose.residue(i).name() ) continue;

			conformation::Residue const & sck( pose.residue(i) );
			numeric::xyzVector<Real> suck_xyz( sck.xyz(1) );
			Real const d = suck_xyz.distance(atom_xyz);
			if ( d > DIST_TH ) continue;

			// std::cerr << "DERIVATOMS " << id.rsd() << " " << id.atomno() << " " << sck.seqpos() << " " << 1 << " " << d << " LPB " << id.rsd() << " " << id.atomno() << std::endl;

			// compute suck F1, F2
			numeric::xyzVector<Real> const f1( atom_xyz.cross( suck_xyz ));
			numeric::xyzVector<Real> const f2( atom_xyz -  suck_xyz );
			Real const dist( f2.length() );
			Real deriv, dummy;
			interp_->interpolate( dist, dummy, deriv );
			// std::cerr << "reg_atm " << deriv << std::endl;
			F1 += ( deriv / dist ) * weights[ suck ] * f1;
			F2 += ( deriv / dist ) * weights[ suck ] * f2;
		}
	}
}

void
SuckerEnergy::eval_intrares_energy(
	conformation::Residue const & /*rsd*/,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap & /*emap*/
) const
{}


/// @brief SuckerEnergy distance cutoff set to the same cutoff used by EtableEnergy, for now
Distance
SuckerEnergy::atomic_interaction_cutoff() const
{
	return interaction_cutoff();
}

/// @details non-virtual accessor for speed; assumption: SuckerEnergy is not inherrited from.
Distance
SuckerEnergy::interaction_cutoff() const
{
	return DIST_TH;
}

/// @brief SuckerEnergy requires that Energies class maintains a TenANeighborGraph
void
SuckerEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
}
core::Size
SuckerEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
