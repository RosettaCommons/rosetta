// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	core/scoring/methods/GaussianOverlapEnergy.cc
/// @brief  Gaussian Overlap Energy
/// @author Ben Borgo (bborgo@genetics.wustl.edu)



#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/methods/GaussianOverlapEnergy.hh>
#include <core/scoring/methods/GaussianOverlapEnergyCreator.hh>
// AUTO-REMOVED #include <core/scoring/NeighborList.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/NumericTraits.hh>
#include <utility/io/izstream.hh>
#include <sstream>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/AtomType.hh>

//Auto Headers
#include <core/id/AtomID.hh>



#define DIST_TH 6.0
//#define TETHER_WT 10.0

namespace core {
namespace scoring {
namespace methods {

GaussianOverlapEnergy::~GaussianOverlapEnergy() {}


/// @details This must return a fresh instance of the GaussianOverlapEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
GaussianOverlapEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new GaussianOverlapEnergy );
}

ScoreTypes
GaussianOverlapEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( gauss );
	return sts;
}


GaussianOverlapEnergy::GaussianOverlapEnergy() :
	parent( methods::EnergyMethodCreatorOP( new GaussianOverlapEnergyCreator ) )
{

/*

	utility::io::izstream izin;
	basic::database::open( izin, basic::options::option[ basic::options::OptionKeys::in::file::native ]() );

	utility::vector1< utility::vector1< Real > > points;

	char s[999];
	izin.getline(s,999); // read in and ignore header line
	while( true ) {
		izin.getline(s,999);
		// std::cerr << "string: '" << s << "'" << std::endl;
		std::istringstream iss(s);
		core::Real tmp;
		utility::vector1< Real > point;
		if( !(iss >> tmp ) ) break; // if no line, done
		point.push_back(tmp);
		iss >> tmp;
		point.push_back(tmp);
		if( iss >> tmp ) point.push_back(tmp);
		points.push_back(point);
	}

debug_assert( points.size() >= 2 );

	// std::cerr << "printing points" << std::endl;
	// for( Size i = 1; i <= points.size(); ++i ) {
	// 	std::cerr << "point: ";
	// 	for( Size j = 1; j <= points[i].size(); ++j ) {
	// 		std::cerr << points[i][j] << ",";
	// 	}
	// 	std::cerr << std::endl;
	// }

	using namespace numeric::interpolation::spline;
	SplineGenerator sg( points[      1      ][1], points[      1      ][2], points[      1      ][3],
							  points[points.size()][1], points[points.size()][2], points[points.size()][3] );

	for( Size i = 2; i < points.size(); ++i ) {
		if( points[i].size() == 2 ) sg.add_known_value( points[i][1], points[i][2] );
		if( points[i].size() == 3 ) sg.add_known_value( points[i][1], points[i][2], points[i][3] );
	}

	interp_ = sg.get_interpolator();

	// core::Real x,y,dy;
	// for( x = 0; x <= 6.0; x += 0.01 ) {
	// 	interp_->interpolate(x,y,dy);
	// 	std::cerr << "SUCKERE " << x << " " << y << " " << dy << std::endl;
	// }

*/
}


/// clone
EnergyMethodOP
GaussianOverlapEnergy::clone() const
{
	return EnergyMethodOP( new GaussianOverlapEnergy() );
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

///
void
GaussianOverlapEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /*pose*/,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

        using namespace core;
        using namespace conformation;

                Real score(0.0);

		 if ( rsd1.seqpos() == rsd2.seqpos() ) return;
//		 if ( rsd1.is_bonded( rsd2 ) ) return;

	        chemical::AA const aa1( rsd1.aa() );
		chemical::AA const aa2( rsd2.aa() );
           if ( aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
                         rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
                         rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;
	   if ( aa1 == chemical::aa_pro && aa2 == chemical::aa_pro ) return;

                for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {
                        Vector const & i_xyz( rsd1.xyz(i) );
                        //Size const i_type( rsd1.atom_type_index(i) );
                        for ( Size j = 1, j_end = rsd2.nheavyatoms(); j <= j_end; ++j ) {
                                Vector const & j_xyz( rsd2.xyz(j) );
                                //Size const j_type( rsd2.atom_type_index(j) );

                                Real const d2( i_xyz.distance_squared( j_xyz ) );
				
                                if( d2 <= 150.0 ){
                                        core::Real r1 = rsd1.atom_type(i).lj_radius();
                                        core::Real r2 = rsd2.atom_type(j).lj_radius();

                                        r1 = .92*r1; r2 = .92*r2;
					if( d2 <= .5*(r1+r2) ){
						score = 1000;}
					else if( d2 <= .8*(r1+r2) ){
					      score = (r1+r2)/(d2*d2*d2*d2*sqrt(d2));}
					else{ 
                                              score = ((sqrt(numeric::NumericTraits<Real>::pi())*r1*r2)/(sqrt(r1*r1+r2*r2)))*(exp(-d2/(r1*r1+r2*r2)));
                                        emap[ gauss ] += score;}}
				else{
					emap[ gauss ] += 0.0;}
}}
}

void
GaussianOverlapEnergy::eval_atom_derivative(
	id::AtomID const & /*id*/,
	pose::Pose const & /*pose*/,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &,
	EnergyMap const &,
	Vector & /*F1*/,
	Vector & /*F2*/
) const
{
	using namespace numeric;
/*
	// id is a SUCK atom
	if( "SUCK" == pose.residue(id.rsd()).name() ) { // loop over atom neighbors and suck
		// conformation::Residue const & sck( pose.residue(id.rsd()) );
	debug_assert( 1 <= id.atomno() && id.atomno() <= 3 );
		if( id.atomno() > 1 ) return;
		numeric::xyzVector<Real> suck_xyz( pose.xyz(id) );

		for( int ir = 1; ir <= (int)pose.total_residue(); ++ir ) {
			conformation::Residue const & res( pose.residue(ir) );
			if( "SUCK" == res.name() ) continue;

			for( int ia = 1; ia <= (int)res.nheavyatoms(); ++ia ) {
				conformation::Atom const & atom( res.atom( ia ) );
				if( !count_atom( atom.type() ) ) continue;

				xyzVector<Real> const & atom_xyz( atom.xyz() );

				Real const d = suck_xyz.distance(atom_xyz);
				if( d > DIST_TH ) continue;

				// std::cerr << "DERIVATOMS " << res.seqpos() << " " << ia << " " << sck.seqpos() << " " << 1 << " " << d << " LPA " << id.rsd() << " " << id.atomno() << std::endl;

				// compute suck F1, F2
				numeric::xyzVector<Real> const f1( suck_xyz.cross( atom_xyz ));
				numeric::xyzVector<Real> const f2( suck_xyz	-	 atom_xyz );
				Real const dist( f2.length() );
				Real deriv, dummy;
				interp_->interpolate( dist, dummy, deriv );
				// std::cerr << "sck_atm " << deriv << std::endl;
				F1 += ( deriv / dist ) * weights[ suck ] * f1;
				F2 += ( deriv / dist ) * weights[ suck ] * f2;

			}
		}

		// // teathers, orig pos from pose data cache
		// using namespace basic;
		// CacheableData const & cd = pose.data().get(SUCKER_ORIG_POSITIONS);
		// CacheableIntXyzVectorMap const & cxmap = dynamic_cast<CacheableIntXyzVectorMap const &>( cd );
		// numeric::xyzVector<Real> const & orig_xyz = cxmap.map().find( id.rsd() )->second;
		//
		// numeric::xyzVector<Real> const f1( suck_xyz.cross( orig_xyz ));
		// numeric::xyzVector<Real> const f2( suck_xyz	-	 orig_xyz  );
		// std::cerr << "intrares dE " << distance(suck_xyz,orig_xyz) << " " << 2*TETHER_WT*distance(suck_xyz,orig_xyz) << std::endl;
		// std::cerr << "  SUCKER ORIG POS " << id.rsd() << " " << orig_xyz.x() << "," << orig_xyz.y() << "," << orig_xyz.z() << std::endl;
		// std::cerr << "  SUCKER		POS " << id.rsd() << " " << suck_xyz.x() << "," << suck_xyz.y() << "," << suck_xyz.z() << std::endl;
		// Real dist  = distance( suck_xyz, orig_xyz );
		// Real deriv = -dist * 2;
		// F1 += ( deriv / dist ) * TETHER_WT * weights[ suck ] * f1;
		// F2 += ( deriv / dist ) * TETHER_WT * weights[ suck ] * f2;


	}	else { // loop over suck residues and suck

		// heavy atoms only
		if( id.atomno() > pose.residue(id.rsd()).nheavyatoms() ) return;
		conformation::Atom const & atom( pose.residue(id.rsd()).atom(id.atomno()) );
		if( !count_atom( atom.type() ) ) return;

		numeric::xyzVector<Real> atom_xyz( atom.xyz() );
		// while( "SUCK" == pose.residue(i).name() ) {
		for( Size i = 1; i <= pose.total_residue(); ++i ) {
			if( "SUCK" != pose.residue(i).name() ) continue;

			conformation::Residue const & sck( pose.residue(i) );
			numeric::xyzVector<Real> suck_xyz( sck.xyz(1) );
			Real const d = suck_xyz.distance(atom_xyz);
			if( d <= DIST_TH ) {

				// std::cerr << "DERIVATOMS " << id.rsd() << " " << id.atomno() << " " << sck.seqpos() << " " << 1 << " " << d << " LPB " << id.rsd() << " " << id.atomno() << std::endl;

				// compute suck F1, F2
				numeric::xyzVector<Real> const f1( atom_xyz.cross( suck_xyz ));
				numeric::xyzVector<Real> const f2( atom_xyz	-	 suck_xyz );
				Real const dist( f2.length() );
				Real deriv, dummy;
				interp_->interpolate( dist, dummy, deriv );
				// std::cerr << "reg_atm " << deriv << std::endl;
				F1 += ( deriv / dist ) * weights[ suck ] * f1;
				F2 += ( deriv / dist ) * weights[ suck ] * f2;

			}
			// --i;

		}

	}
*/

}

void
GaussianOverlapEnergy::eval_intrares_energy(
	conformation::Residue const & /*rsd*/,
	pose::Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap & /*emap*/
) const {
	// using namespace basic;
	// if( "SUCK" == rsd.name() ) {
	// 	CacheableData const & cd = pose.data().get(SUCKER_ORIG_POSITIONS);
	// 	CacheableIntXyzVectorMap const & cxmap = dynamic_cast<CacheableIntXyzVectorMap const &>( cd );
	// 	numeric::xyzVector<Real> const & orig_xyz = cxmap.map().find( rsd.seqpos() )->second;
	// 	numeric::xyzVector<Real> const & curr_xyz( rsd.xyz(1) );
	// 	Real const d2( distance_squared(orig_xyz,curr_xyz) );
	// 	std::cerr << "intrareas E " << distance(curr_xyz,orig_xyz) << " " << TETHER_WT*d2 << std::endl;
	// 	std::cerr << "  SUCKER ORIG POS " << rsd.seqpos() << " " << orig_xyz.x() << "," << orig_xyz.y() << "," << orig_xyz.z() << std::endl;
	// 	std::cerr << "  SUCKER		POS " << rsd.seqpos() << " " << curr_xyz.x() << "," << curr_xyz.y() << "," << curr_xyz.z() << std::endl;
	// 	emap[ suck ] += TETHER_WT*d2;
	// }
}


/// @brief GaussianOverlapEnergy distance cutoff set to the same cutoff used by EtableEnergy, for now
Distance
GaussianOverlapEnergy::atomic_interaction_cutoff() const
{
	return interaction_cutoff();
}

/// @details non-virtual accessor for speed; assumption: GaussianOverlapEnergy is not inherrited from.
Distance
GaussianOverlapEnergy::interaction_cutoff() const
{
	return DIST_TH;
}

/// @brief GaussianOverlapEnergy requires that Energies class maintains a TenANeighborGraph
void
GaussianOverlapEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required*/ ) const
{
}
core::Size
GaussianOverlapEnergy::version() const
{
	return 1; // Initial versioning
}


}
}
}
