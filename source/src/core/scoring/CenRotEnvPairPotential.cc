// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CenRotEnvPairPotential.cc
/// @brief  CenRot version of cen pair/env
/// @author Yuan Liu


#include <core/scoring/CenRotEnvPairPotential.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/kinematics/Stub.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/prof.hh>
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>


namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////////////////////
CenRotEnvPairPotential::CenRotEnvPairPotential() {
	using namespace numeric::interpolation::spline;
	SIGMOID_SLOPE = 12.0;
	cen_dist_cutoff_12_pad = 13.5*13.5;

	Size const max_aa( 20 ); // just the standard aa for now
	pairsplines_.resize(max_aa, utility::vector1< CubicSpline >(max_aa));
	envsplines_.resize(max_aa, CubicSpline());
	angsplines_.resize(max_aa, utility::vector1< BicubicSpline >(max_aa));
	//dihsplines_.resize(max_aa, utility::vector1< BicubicSpline >(max_aa));

	std::string tag,line;
	chemical::AA aa1, aa2;

	//////////////////////////////
	//read cen_rot_pair_params.txt
	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_rot_pair_params.txt");

	while ( stream.getline( line ) ) {
		std::istringstream l(line);
		l >> tag;

		if ( tag=="PAIR:" ) {
			//read one spline for aa1 and aa2
			l >> aa1 >> aa2;
			CubicSpline &cspline = pairsplines_[std::min(aa1,aa2)][std::max(aa1,aa2)];

			//and train
			Real startp, endp;
			Size nbin;
			l >> startp >> endp >> nbin;

			numeric::MathVector<Real> results(nbin);
			for ( Size i=0; i<nbin; i++ ) l >> results(i);

			std::pair< Real, Real > unused;
			unused = std::make_pair( 0.0, 0.0 );
			cspline.train(
				e_Natural,
				startp, (endp-startp)/(nbin-1),
				results, unused
			);
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format for cen_rot_pair_params.txt");
		}

		if ( l.fail() ) {
			utility_exit_with_message("bad format for cen_rot_pair_params.txt");
		}
	}
	stream.close();

	/////////////////////////////
	//read cen_rot_env_params.txt
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_rot_env_params.txt");
	while ( stream.getline( line ) ) {
		std::istringstream l(line);
		l >> tag;

		if ( tag=="ENV:" ) {
			l >> aa1;
			//std::cout << aa1 << std::endl;
			CubicSpline &cspline = envsplines_[aa1];
			numeric::MathVector<Real> results(40);
			for ( Size i=0; i<40; i++ ) l >> results(i);
			//for (Size i=0; i<40; i++) std::cout << results(i) << std::endl;
			std::pair< Real, Real > unused;
			unused = std::make_pair( 0.0, 0.0 );
			cspline.train(
				e_Natural,
				1.0, 1.0,
				results, unused
			);
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format for cen_rot_env_params.txt");
		}

		if ( l.fail() ) {
			utility_exit_with_message("bad format for cen_rot_env_params.txt");
		}
	}
	stream.close();

	////////////////////////////////
	//read cen_rot_cbeta_params.txt
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_rot_cbeta_params.txt");
	while ( stream.getline( line ) ) {
		std::istringstream l(line);
		l >> tag;

		if ( tag == "CBETA6:" ) {
			// do cbeta stuff
			l >> tag;
			if ( tag == "SHIFT" ) {
				Real shift; l >> shift;
				cbeta6_.shift(shift);
			} else if ( tag == "GAUSSIAN" || tag == "SIGMOID" ) {
				Size ngauss; l >> ngauss;
				numeric::xyzVector<Real> g;
				for ( core::Size i=1; i<=ngauss; ++i ) {
					l >> g[0] >> g[1] >> g[2];
					if ( tag == "GAUSSIAN" ) {
						cbeta6_.add_gaussian(g);
					} else {
						cbeta6_.add_sigmoid(g);
					}
				}
			}
		} else if ( tag == "CBETA12:" ) {
			// do cbeta stuff
			l >> tag;
			if ( tag == "SHIFT" ) {
				Real shift; l >> shift;
				cbeta12_.shift(shift);
			} else if ( tag == "GAUSSIAN" || tag == "SIGMOID" ) {
				Size ngauss; l >> ngauss;
				numeric::xyzVector<Real> g;
				for ( core::Size i=1; i<=ngauss; ++i ) {
					l >> g[0] >> g[1] >> g[2];
					if ( tag == "GAUSSIAN" ) {
						cbeta12_.add_gaussian(g);
					} else {
						cbeta12_.add_sigmoid(g);
					}
				}
			}
		} else if ( tag == "CENPACK:" ) {
			// do cenpack stuff
			l >> tag;
			if ( tag == "SHIFT" ) {
				Real shift; l >> shift;
				cenpack_.shift(shift);
			} else if ( tag == "GAUSSIAN" || tag == "SIGMOID" ) {
				Size ngauss; l >> ngauss;
				numeric::xyzVector<Real> g;
				for ( core::Size i=1; i<=ngauss; ++i ) {
					l >> g[0] >> g[1] >> g[2];
					if ( tag == "GAUSSIAN" ) {
						cenpack_.add_gaussian(g);
					} else {
						cenpack_.add_sigmoid(g);
					}
				}
			}
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format of cen_rot_cbeta_params.txt");
		}

		if ( l.fail() ) {
			utility_exit_with_message("fail to read cen_rot_cbeta_params.txt");
		}
	}
	stream.close();

	//////////////////////////////////
	//read cen_rot_pair_ang_params.txt
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_rot_pair_ang_params.txt");
	while ( stream.getline( line ) ) {
		std::istringstream l(line);
		l >> tag;

		if ( tag=="ANGLE1:" ) {
			//read one spline for aa1 and aa2
			l >> aa1 >> aa2;
			BicubicSpline &cspline = angsplines_[aa1][aa2];

			//and train
			Real start_r, end_r;
			Size nbin_r;
			l >> start_r >> end_r >> nbin_r;
			Real start_a, end_a;
			Size nbin_a;
			l >> start_a >> end_a >> nbin_a;

			numeric::MathMatrix< Real > results( nbin_r, nbin_a );
			for ( Size i=0; i<nbin_r; i++ ) {
				stream.getline( line );
				std::istringstream ll(line);
				for ( Size j=0; j<nbin_a; j++ ) {
					ll >> results(i,j);
				}
			}

			std::pair< Real, Real > unused[2];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			BorderFlag periodic_boundary[2] = { e_Natural, e_Natural };
			//BorderFlag periodic_boundary[2] = { e_FirstDer, e_FirstDer };
			Real start_vals[2] = {start_r, start_a};
			Real deltas[2] = {(end_r-start_r)/(nbin_r-1), (end_a-start_a)/(nbin_a-1)};
			bool lincont[2] = {false,false};

			cspline.train(periodic_boundary, start_vals, deltas, results, lincont, unused);
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format for cen_rot_pair_ang_params.txt");
		}

		if ( l.fail() ) {
			utility_exit_with_message("bad format for cen_rot_pair_ang_params.txt");
		}
	}
	stream.close();
}

void
CenRotEnvPairPotential::evaluate_cen_rot_pair_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real & pair_contribution
) const {
	pair_contribution = 0.0;
	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ < 9 ) return;

	pair_contribution = pairsplines_[std::min(aa1,aa2)][std::max(aa1,aa2)].F(cendist);
}

void
CenRotEnvPairPotential::evaluate_cen_rot_pair_deriv(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real & dpair_dr
) const {
	// only returns d_E/d_cendist
	dpair_dr = 0.0;
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ < 9 ) return;
	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );
	dpair_dr = pairsplines_[std::min(aa1,aa2)][std::max(aa1,aa2)].dF(cendist);
}

void
CenRotEnvPairPotential::evaluate_cen_rot_pair_orientation_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real & ang1_contribution,
	Real & ang2_contribution,
	Real & dih_contribution
) const {
	ang1_contribution = 0.0;
	ang2_contribution = 0.0;
	dih_contribution = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;
	//skip ALA, GLY, PRO
	if ( rsd1.aa()==chemical::aa_gly || rsd1.aa()==chemical::aa_ala || rsd1.aa()==chemical::aa_pro
			|| rsd2.aa()==chemical::aa_gly || rsd2.aa()==chemical::aa_ala || rsd2.aa()==chemical::aa_pro ) return;

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ < 9 ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	//save angle
	core::kinematics::Stub::Vector ra, rb, rc, rd;
	ra = rsd1.atom(rsd1.nbr_atom()).xyz();
	rb = rsd1.atom("CEN").xyz();
	rc = rsd2.atom("CEN").xyz();
	rd = rsd2.atom(rsd2.nbr_atom()).xyz();

	Real ang1 = cos(numeric::angle_radians(ra,rb,rc));
	Real ang2 = cos(numeric::angle_radians(rb,rc,rd));

	ang1_contribution = angsplines_[aa1][aa2].F(cendist, ang1);
	ang2_contribution = angsplines_[aa2][aa1].F(cendist, ang2);
}

void
CenRotEnvPairPotential::evaluate_cen_rot_pair_orientation_deriv(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real const ang,
	Real & dE_dr,
	Real & dE_d_ang,
	Real & dE_d_dih
) const {
	// only retures d_E/d_ang
	dE_d_dih = 0.0;
	dE_d_ang = 0.0;
	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;
	//skip ALA, GLY, PRO
	if ( rsd1.aa()==chemical::aa_gly || rsd1.aa()==chemical::aa_ala || rsd1.aa()==chemical::aa_pro
			|| rsd2.aa()==chemical::aa_gly || rsd2.aa()==chemical::aa_ala || rsd2.aa()==chemical::aa_pro ) return;

	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ < 9 ) return;

	Real cos_ang = cos(ang);
	Real wt_ang = -sin(ang);
	dE_d_ang = angsplines_[rsd1.aa()][rsd2.aa()].dFdy(cendist, cos_ang) * wt_ang;
	dE_dr = angsplines_[rsd1.aa()][rsd2.aa()].dFdx(cendist, cos_ang);
}

void
CenRotEnvPairPotential::evaluate_cen_rot_env_and_cbeta_score(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real & env_contribution,
	Real & cbeta6_contribution,
	Real & cbeta12_contribution
) const {

	env_contribution = 0.0;
	cbeta6_contribution = 0.0;
	cbeta12_contribution = 0.0;

	if ( !rsd.is_protein() ) return;

	SigmoidWeightedCenList<Real> const & cenlist( cenlist_from_pose( pose ));
	int const position ( rsd.seqpos() );
	Real fcen6 ( cenlist.fcen6(position) );
	Real fcen10 ( cenlist.fcen10(position) );
	Real fcen12 ( cenlist.fcen12(position) );

	fcen10 = std::min(40.0, fcen10);
	env_contribution = envsplines_[ rsd.aa() ].F( fcen10 );
	fcen6 = std::min(60.0, fcen6);
	fcen12 = std::min(60.0, fcen12);
	cbeta6_contribution = cbeta6_.func( fcen6 );
	cbeta12_contribution = cbeta12_.func( fcen12 );
}

void
CenRotEnvPairPotential::evaluate_cen_rot_env_and_cbeta_deriv(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	numeric::xyzVector<Real> & f2_cen_env,
	numeric::xyzVector<Real> & f2_cen_cb6,
	numeric::xyzVector<Real> & f2_cen_cb12,
	numeric::xyzVector<Real> & f2_cb_env,
	numeric::xyzVector<Real> & f2_cb_cb6,
	numeric::xyzVector<Real> & f2_cb_cb12
) const {

	f2_cen_env = numeric::xyzVector<Real>(0.0, 0.0, 0.0);
	f2_cen_cb6 = numeric::xyzVector<Real>(0.0, 0.0, 0.0);
	f2_cen_cb12 = numeric::xyzVector<Real>(0.0, 0.0, 0.0);
	f2_cb_env = numeric::xyzVector<Real>(0.0, 0.0, 0.0);
	f2_cb_cb6 = numeric::xyzVector<Real>(0.0, 0.0, 0.0);
	f2_cb_cb12 = numeric::xyzVector<Real>(0.0, 0.0, 0.0);

	if ( !rsd.is_protein() ) return;

	SigmoidWeightedCenList<Real> const & cenlist( cenlist_from_pose( pose ));
	SigmoidWeightedCenList< numeric::xyzVector< Real > > const & dcenlist( dcenlist_from_pose( pose ));

	int const position ( rsd.seqpos() );
	// derivative of centroid count wrt x
	numeric::xyzVector< Real > const dcentroids6_dx  ( dcenlist.fcen6(position) );
	numeric::xyzVector< Real > const dcentroids10_dx ( dcenlist.fcen10(position) );
	numeric::xyzVector< Real > const dcentroids12_dx ( dcenlist.fcen12(position) );

	// derivative of score w.r.t centroid count
	Real fcen6  ( cenlist.fcen6(position) );
	Real fcen10 ( cenlist.fcen10(position) );
	Real fcen12 ( cenlist.fcen12(position) );

	fcen10 = std::min(40.0, fcen10);
	fcen6 = std::min(60.0, fcen6);
	fcen12 = std::min(60.0, fcen12);

	// in the cenrot model, env and cbeta term are define by #CB around a CEN
	// which means for these terms, we need to take care about forces on both CEN and CB
	f2_cen_env = envsplines_[ rsd.aa() ].dF( fcen10 ) * dcentroids10_dx;
	f2_cen_cb6 = cbeta6_.dfunc( fcen6 ) * dcentroids6_dx;
	f2_cen_cb12 = cbeta12_.dfunc( fcen12 ) * dcentroids12_dx;

	// change in others' context:
	Vector atom_x = rsd.atom( rsd.nbr_atom() ).xyz(); //CB for cenrot
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node(position)->const_edge_list_begin(),
			irue = energy_graph.get_node(position)->const_edge_list_end();
			iru != irue; ++iru ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		Size const j( edge->get_other_ind(position) );
		conformation::Residue const & rsd2 ( pose.residue(j) );
		if ( !rsd2.is_protein() ) continue;
		//other CEN affected by this CB
		numeric::xyzVector<Real> cenvec = rsd2.atom("CEN").xyz() - atom_x;
		Real const cendist = cenvec.length_squared();

		if ( cendist > cen_dist_cutoff_12_pad ) continue;

		Real x = sqrt(cendist);
		numeric::xyzVector<Real> gradx = cenvec/x;

		Real d6,d10,d12, z;

		z = SIGMOID_SLOPE*(x-6);
		if ( z>-30 && z<30 ) {
			Real e6 = exp(z);
			d6 = SIGMOID_SLOPE*e6 / ((1+e6)*(1+e6));
		} else {
			d6 = 0;
		}

		z = SIGMOID_SLOPE*(x-10);
		if ( z>-20 && z<20 ) {
			Real e10 = exp(z);
			d10 = SIGMOID_SLOPE*e10 / ((1+e10)*(1+e10));
		} else {
			d10 = 0;
		}

		z = SIGMOID_SLOPE*(x-12);
		if ( z>-20 && z<20 ) {
			Real e12 = exp(z);
			d12 = SIGMOID_SLOPE*e12 / ((1+e12)*(1+e12));
		} else {
			d12 = 0;
		}

		fcen10 = std::min(40.0, cenlist.fcen10(j));
		fcen6 = std::min(60.0, cenlist.fcen6(j));
		fcen12 = std::min(60.0, cenlist.fcen12(j));

		f2_cb_env += envsplines_[ rsd2.aa() ].dF( fcen10 ) * (d10*gradx);
		f2_cb_cb6 += cbeta6_.dfunc( fcen6 ) * (d6*gradx);
		f2_cb_cb12 += cbeta12_.dfunc( fcen12 ) *((d12-d6)*gradx);
	}
}

////////////////////////////////////////////////////////
// copy from smooth env pair potential
////////////////////////////////////////////////////////

void
CenRotEnvPairPotential::fill_smooth_cenlist(
	SigmoidWeightedCenList<Real> & cenlist,
	Size const res1,
	Size const,
	Real const cendist
) const {
	Real interp6, interp10, interp12, z;
	z = SIGMOID_SLOPE*(cendist-6);
	if ( z<-30 ) {
		interp6 = 1;
	} else if ( z>30 ) {
		interp6 = 0;
	} else {
		interp6 = 1 / (1+exp(z));
	}
	cenlist.fcen6(res1) += interp6;

	z = SIGMOID_SLOPE*(cendist-10);
	if ( z<-20 ) {
		interp10 = 1;
	} else if ( z>20 ) {
		interp10 = 0;
	} else {
		interp10 = 1 / (1+exp(z));
	}
	cenlist.fcen10(res1) += interp10;

	z = SIGMOID_SLOPE*(cendist-12);
	if ( z<-20 ) {
		interp12 = 1;
	} else if ( z>20 ) {
		interp12 = 0;
	} else {
		interp12 = 1 / (1+exp(z));
	}
	cenlist.fcen12(res1) += interp12 - interp6;  // from 6->12 A
}

void
CenRotEnvPairPotential::fill_smooth_dcenlist(
	SigmoidWeightedCenList< numeric::xyzVector<Real> > & dcenlist,
	Size const res1,
	Size const,
	numeric::xyzVector<Real> const & cenvec  // vector between centroids
) const {
	Real x = cenvec.length();
	numeric::xyzVector<Real> gradx = cenvec/x;

	Real d6, d12, z;
	z = SIGMOID_SLOPE*(x-6);
	if ( z>-30 && z<30 ) {
		Real e6 = exp(z);
		d6 = SIGMOID_SLOPE*e6 / ((1+e6)*(1+e6));
		dcenlist.fcen6(res1) += d6*gradx;
	} else {
		d6 = 0;
	}

	z = SIGMOID_SLOPE*(x-10);
	if ( z>-20 && z<20 ) {
		Real e10 = exp(z);
		Real d10 = SIGMOID_SLOPE*e10 / ((1+e10)*(1+e10));
		dcenlist.fcen10(res1) += d10*gradx;
	}

	z = SIGMOID_SLOPE*(x-12);
	if ( z>-20 && z<20 ) {
		Real e12 = exp(z);
		d12 = SIGMOID_SLOPE*e12 / ((1+e12)*(1+e12));
	} else {
		d12 = 0;
	}
	dcenlist.fcen12(res1) += (d12-d6)*gradx; // from 6->12 A
}

void
CenRotEnvPairPotential::compute_centroid_environment(
	pose::Pose & pose
) const {
	SigmoidWeightedCenList<Real> & cenlist( nonconst_cenlist_from_pose( pose ));

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Size const nres( energy_graph.num_nodes() );

	/// calculate the cenlist info only if it has not been calculated since the last score evaluation
	if ( cenlist.calculated() ) return;

	cenlist.initialize( pose.total_residue(), 0 );  //different from smoothenvpairpotential, starts from 0

	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;

		//we need to go through all the edges here rather than just upper ones
		//because the edges aren't sym due to the new env def
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_other_ind(i) );
			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;

			//Real const cendist = edge->square_distance();
			numeric::xyzVector<Real> cenvec =
				rsd2.atom( rsd2.nbr_atom() ).xyz() - rsd1.atom("CEN").xyz();
			Real const cendist = cenvec.length_squared();

			if ( cendist <= cen_dist_cutoff_12_pad ) {
				fill_smooth_cenlist( cenlist, i, j, sqrt(cendist) );
			}
		}
	}

	cenlist.calculated() = true;
}

void
CenRotEnvPairPotential::compute_dcentroid_environment(
	pose::Pose & pose
) const {
	nonconst_cenlist_from_pose( pose );
	SigmoidWeightedCenList< numeric::xyzVector< Real > > & dcenlist( nonconst_dcenlist_from_pose( pose ));

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Size const nres( energy_graph.num_nodes() );

	if ( dcenlist.calculated() ) return;

	dcenlist.initialize( pose.total_residue(), numeric::xyzVector< Real >(0,0,0) );

	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_other_ind(i) );
			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;

			numeric::xyzVector<Real> cenvec =
				rsd2.atom( rsd2.nbr_atom() ).xyz() - rsd1.atom("CEN").xyz();
			Real const cendist = cenvec.length_squared();

			if ( cendist <= cen_dist_cutoff_12_pad ) {
				fill_smooth_dcenlist( dcenlist, i, j, cenvec );
			}
		}
	}

	dcenlist.calculated() = true;
}

void
CenRotEnvPairPotential::finalize( pose::Pose & pose ) const {
	SigmoidWeightedCenList<Real> & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;
	SigmoidWeightedCenList<numeric::xyzVector <Real> > & dcenlist( nonconst_dcenlist_from_pose( pose ));
	dcenlist.calculated() = false;
}

/// @details Pose must already contain a cenlist object or this method will fail.
SigmoidWeightedCenList< Real > const &
CenRotEnvPairPotential::cenlist_from_pose( pose::Pose const & pose ) const {
	using namespace core::pose::datacache;
	return *( utility::pointer::static_pointer_cast< SigmoidWeightedCenList<Real> const > ( pose.data().get_const_ptr( CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ) ));

}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenlist object, places it in the pose, and returns
/// a non-const reference to it.
SigmoidWeightedCenList< Real > &
CenRotEnvPairPotential::nonconst_cenlist_from_pose( pose::Pose & pose ) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ) ) {
		return *( utility::pointer::static_pointer_cast< SigmoidWeightedCenList<Real> > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ) ));
	}
	// else
	SigmoidWeightedCenListRealOP cenlist( new SigmoidWeightedCenList< Real > );
	pose.data().set( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST, cenlist );
	return *cenlist;
}

/// @details Pose must already contain a cenlist object or this method will fail.
SigmoidWeightedCenList< numeric::xyzVector< Real > > const &
CenRotEnvPairPotential::dcenlist_from_pose( pose::Pose const & pose ) const {
	using namespace core::pose::datacache;
	return *( utility::pointer::static_pointer_cast< SigmoidWeightedCenList<numeric::xyzVector<Real> > const > ( pose.data().get_const_ptr( CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ) ));

}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
SigmoidWeightedCenList< numeric::xyzVector< Real > > &
CenRotEnvPairPotential::nonconst_dcenlist_from_pose( pose::Pose & pose ) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ) ) {
		return *( utility::pointer::static_pointer_cast< SigmoidWeightedCenList<numeric::xyzVector<Real> > > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ) ));
	}
	// else
	SigmoidWeightedCenListVectorOP cenlist( new SigmoidWeightedCenList< numeric::xyzVector< Real > > );
	pose.data().set( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST, cenlist );
	return *cenlist;
}


}
}
