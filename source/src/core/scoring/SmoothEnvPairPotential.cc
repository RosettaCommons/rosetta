// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SmoothEnvPairPotential.cc
/// @brief  Smooth, differentiable version of centroid env and pair terms
/// @author Frank DiMaio


#include <core/scoring/SmoothEnvPairPotential.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>
#include <utility/io/izstream.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>


#include <utility/vector1.hh>


namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////////////////////

Real SmoothScoreTermCoeffs::func( Real x ) const {
	Real y = shift_;
	for ( core::Size i=1; i<=sigmoid_coeffs_.size(); ++i ) {
		numeric::xyzVector< Real > const &q = sigmoid_coeffs_[i];
		y += q[0] / (1 + exp( -q[2]*x - q[1] ));
	}
	for ( core::Size i=1; i<=gaussian_coeffs_.size(); ++i ) {
		numeric::xyzVector< Real > const &n = gaussian_coeffs_[i];
		y += n[0] * exp( -(x-n[1])*(x-n[1]) / (2*n[2]*n[2]) );
	}
	return y;
}


Real SmoothScoreTermCoeffs::dfunc( Real x ) const {
	Real dy = 0;
	for ( core::Size i=1; i<=sigmoid_coeffs_.size(); ++i ) {
		numeric::xyzVector< Real > const &q = sigmoid_coeffs_[i];
		Real e = exp( -q[2]*x - q[1] );
		Real d = (1 + e);
		dy += q[0]*q[2]*e / (d*d);
	}
	for ( core::Size i=1; i<=gaussian_coeffs_.size(); ++i ) {
		numeric::xyzVector< Real > const &n = gaussian_coeffs_[i];
		dy += n[0] * (n[1]-x) / (n[2]*n[2]) * exp( -(n[1]-x)*(n[1]-x) / (2*n[2]*n[2]) );
	}
	return dy;
}

///////////////////////////////////////////////////////////////////////////////////////////////

SmoothEnvPairPotential::SmoothEnvPairPotential() {
	SIGMOID_SLOPE = 12.0;
	cen_dist_cutoff_12_pad = 13.5*13.5;

	// load the data
	Size const max_aa( 20 ); // just the standard aa's for now
	env_.resize(max_aa);
	pair_.resize(max_aa, utility::vector1< SmoothScoreTermCoeffs >(max_aa) );

	std::string tag,line;
	chemical::AA aa1,aa2;

	utility::io::izstream stream;
	basic::database::open( stream, "scoring/score_functions/centroid_smooth/cen_smooth_params.txt");

	while ( getline( stream, line ) ) {
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
		} else if ( tag == "ENV:" ) {
			// do env stuff
			l >> aa1;
			SmoothScoreTermCoeffs & currEnv = env_[aa1];

			l >> tag;
			if ( tag == "SHIFT" ) {
				Real shift; l >> shift;
				currEnv.shift(shift);
			} else if ( tag == "GAUSSIAN" || tag == "SIGMOID" ) {
				Size ngauss; l >> ngauss;
				numeric::xyzVector<Real> g;
				for ( core::Size i=1; i<=ngauss; ++i ) {
					l >> g[0] >> g[1] >> g[2];
					if ( tag == "GAUSSIAN" ) {
						currEnv.add_gaussian(g);
					} else {
						currEnv.add_sigmoid(g);
					}
				}
			}
		} else if ( tag == "PAIR:" ) {
			// do pair stuff
			// do env stuff
			l >> aa1 >> aa2;
			SmoothScoreTermCoeffs & currPair =
				pair_[std::min(aa1,aa2)][std::max(aa1,aa2)];  // symmetrical; just store one half

			l >> tag;
			if ( tag == "SHIFT" ) {
				Real shift; l >> shift;
				currPair.shift(shift);
			} else if ( tag == "GAUSSIAN" || tag == "SIGMOID" ) {
				Size ngauss; l >> ngauss;
				numeric::xyzVector<Real> g;
				for ( core::Size i=1; i<=ngauss; ++i ) {
					l >> g[0] >> g[1] >> g[2];
					if ( tag == "GAUSSIAN" ) {
						currPair.add_gaussian(g);
					} else {
						currPair.add_sigmoid(g);
					}
				}
			}
		} else if ( tag != "#" ) {
			utility_exit_with_message("bad format for cen_smooth_params.txt");
		}

		if ( l.fail() ) utility_exit_with_message("bad format for cen_smooth_params.txt");
	}
}

void
SmoothEnvPairPotential::fill_smooth_cenlist(
	SigmoidWeightedCenList<Real> & cenlist,
	Size const res1,
	Size const res2,
	Real const cendist
) const {
	Real interp6, interp10, interp12, z;
	z = SIGMOID_SLOPE*(cendist-6);
	if ( z<-20 ) {
		interp6 = 1;
	} else if ( z>20 ) {
		interp6 = 0;
	} else {
		interp6 = 1 / (1+exp(z));
	}
	cenlist.fcen6(res1) += interp6;
	cenlist.fcen6(res2) += interp6;

	z = SIGMOID_SLOPE*(cendist-10);
	if ( z<-20 ) {
		interp10 = 1;
	} else if ( z>20 ) {
		interp10 = 0;
	} else {
		interp10 = 1 / (1+exp(z));
	}
	cenlist.fcen10(res1) += interp10;
	cenlist.fcen10(res2) += interp10;

	z = SIGMOID_SLOPE*(cendist-12);
	if ( z<-20 ) {
		interp12 = 1;
	} else if ( z>20 ) {
		interp12 = 0;
	} else {
		interp12 = 1 / (1+exp(z));
	}
	cenlist.fcen12(res1) += interp12 - interp6;  // from 6->12 A
	cenlist.fcen12(res2) += interp12 - interp6;  // from 6->12 A
}


void
SmoothEnvPairPotential::fill_smooth_dcenlist(
	SigmoidWeightedCenList< numeric::xyzVector<Real> > & dcenlist,
	Size const res1,
	Size const res2,
	numeric::xyzVector<Real> const & cenvec  // vector between centroids
) const {
	Real x = cenvec.length();
	numeric::xyzVector<Real> gradx = cenvec/x;

	Real d6, d12 , z;
	z =  SIGMOID_SLOPE*(x-6);
	if ( z>-20 && z<20 ) {
		Real e6 = exp(z);
		d6 = SIGMOID_SLOPE*e6 / ((1+e6)*(1+e6));
	} else {
		d6 = 0;
	}
	dcenlist.fcen6(res1) += d6*gradx;
	dcenlist.fcen6(res2) -= d6*gradx;

	z = SIGMOID_SLOPE*(x-10);
	if ( z>-20 && z<20 ) {
		Real e10 = exp(z);
		Real d10 = SIGMOID_SLOPE*e10 / ((1+e10)*(1+e10));
		dcenlist.fcen10(res1) += d10*gradx;
		dcenlist.fcen10(res2) -= d10*gradx;
	}

	z = SIGMOID_SLOPE*(x-12);
	if ( z>-20 && z<20 ) {
		Real e12 = exp(z);
		d12 = SIGMOID_SLOPE*e12 / ((1+e12)*(1+e12));
	} else {
		d12 = 0;
	}
	dcenlist.fcen12(res1) += (d12-d6)*gradx;
	dcenlist.fcen12(res2) -= (d12-d6)*gradx;
}

void
SmoothEnvPairPotential::compute_centroid_environment(
	pose::Pose & pose
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
	}

	SigmoidWeightedCenList<Real> & cenlist( nonconst_cenlist_from_pose( pose ));

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Size const nres( energy_graph.num_nodes() );

	/// calculate the cenlist info only if it has not been calculated since the last score evaluation
	if ( cenlist.calculated() ) return;

	cenlist.initialize( pose.total_residue(), 1 );  // every res has 1 neighbor (itself)

	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;

			//Real const cendist = edge->square_distance();
			numeric::xyzVector<Real> cenvec =
				rsd2.atom( rsd2.nbr_atom() ).xyz() - rsd1.atom( rsd1.nbr_atom() ).xyz();
			Real const cendist = cenvec.length_squared();

			if ( cendist <= cen_dist_cutoff_12_pad ) {
				fill_smooth_cenlist( cenlist, i, j, sqrt(cendist) );
			}
		}
	}

	// symetrize cenlist (if necessary)
	if ( symminfo ) {
		for ( Size i = 1; i <= nres; ++i ) {
			conformation::Residue const & rsd1 ( pose.residue(i) );
			if ( !rsd1.is_protein() ) continue;
			if ( !symminfo->bb_is_independent( i ) ) continue; // probably unnecessary...
			utility::vector1< core::Size > bbclones = symminfo->bb_clones( i );
			for ( Size j = 1; j <= bbclones.size(); ++j ) {
				cenlist.fcen6( bbclones[j] ) = cenlist.fcen6(i);
				cenlist.fcen10( bbclones[j] ) = cenlist.fcen10(i);
				cenlist.fcen12( bbclones[j] ) = cenlist.fcen12(i);
			}
		}
	}

	cenlist.calculated() = true;
}

void
SmoothEnvPairPotential::compute_dcentroid_environment(
	pose::Pose & pose
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	core::conformation::symmetry::SymmetryInfoCOP symminfo(0);
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
	}

	nonconst_cenlist_from_pose( pose ); // must call 1st (?)
	SigmoidWeightedCenList< numeric::xyzVector< Real > > & dcenlist( nonconst_dcenlist_from_pose( pose ));

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Size const nres( energy_graph.num_nodes() );

	if ( dcenlist.calculated() ) return;

	dcenlist.initialize( pose.total_residue(), numeric::xyzVector< Real >(0,0,0) );

	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;
		if ( symminfo && !symminfo->bb_is_independent( i ) ) continue;

		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;

			numeric::xyzVector<Real> cenvec =
				rsd2.atom( rsd2.nbr_atom() ).xyz() - rsd1.atom( rsd1.nbr_atom() ).xyz();
			Real const cendist = cenvec.length_squared();

			if ( cendist <= cen_dist_cutoff_12_pad ) {
				fill_smooth_dcenlist( dcenlist, i, j, cenvec );
			}
		}
	}

	// symetrize cenlist (if necessary)
	if ( symminfo ) {
		core::conformation::symmetry::SymmetricConformation const &symmconf =
			dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation());
		for ( Size i = 1; i <= nres; ++i ) {
			conformation::Residue const & rsd1 ( pose.residue(i) );
			if ( !rsd1.is_protein() ) continue;
			if ( !symminfo->bb_is_independent( i ) ) continue; // probably unnecessary...
			utility::vector1< core::Size > bbclones = symminfo->bb_clones( i );
			for ( Size j = 1; j <= bbclones.size(); ++j ) {
				dcenlist.fcen6( bbclones[j] ) = symmconf.apply_transformation_norecompute( dcenlist.fcen6(i), bbclones[j], i, true );
				dcenlist.fcen10( bbclones[j] ) = symmconf.apply_transformation_norecompute( dcenlist.fcen10(i), bbclones[j], i, true );
				dcenlist.fcen12( bbclones[j] ) = symmconf.apply_transformation_norecompute( dcenlist.fcen12(i), bbclones[j], i, true );
			}
		}
	}

	dcenlist.calculated() = true;
}


void
SmoothEnvPairPotential::finalize( pose::Pose & pose ) const {
	SigmoidWeightedCenList<Real> & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;
	SigmoidWeightedCenList<numeric::xyzVector <Real> > & dcenlist( nonconst_dcenlist_from_pose( pose ));
	dcenlist.calculated() = false;
}

void
SmoothEnvPairPotential::evaluate_env_and_cbeta_scores(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real & env_score,
	Real & cb_score6,
	Real & cb_score12
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	env_score = 0.0;
	cb_score6  = 0.0;
	cb_score12 = 0.0;

	if ( !rsd.is_protein() ) return;

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetryInfoCOP symminfo
			= dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
		if ( !symminfo->bb_is_independent( rsd.seqpos() ) ) return;
	}

	SigmoidWeightedCenList<Real> const & cenlist( cenlist_from_pose( pose ));
	int const position ( rsd.seqpos() );

	// derivative of score w.r.t centroid count
	Real const fcen6  ( cenlist.fcen6(position) );
	Real const fcen10 ( cenlist.fcen10(position) );
	Real const fcen12 ( cenlist.fcen12(position) );

	env_score = env_[ rsd.aa() ].func( fcen10 );
	cb_score6  = cbeta6_.func( fcen6 );
	cb_score12 = cbeta12_.func( fcen12 );
}

void
SmoothEnvPairPotential::evaluate_env_and_cbeta_deriv(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	numeric::xyzVector<Real> & d_env_score,
	numeric::xyzVector<Real> & d_cb_score6,
	numeric::xyzVector<Real> & d_cb_score12
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	d_env_score = numeric::xyzVector<Real>(0.0,0.0,0.0);
	d_cb_score6  = numeric::xyzVector<Real>(0.0,0.0,0.0);
	d_cb_score12 = numeric::xyzVector<Real>(0.0,0.0,0.0);

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::conformation::symmetry::SymmetryInfoCOP symminfo
			= dynamic_cast<const core::conformation::symmetry::SymmetricConformation & >( pose.conformation()).Symmetry_Info();
		if ( !symminfo->bb_is_independent( rsd.seqpos() ) ) return;
	}

	if ( !rsd.is_protein() ) return;

	SigmoidWeightedCenList<Real> const & cenlist( cenlist_from_pose( pose ));
	SigmoidWeightedCenList< numeric::xyzVector< Real > > const & dcenlist( dcenlist_from_pose( pose ));

	int const position ( rsd.seqpos() );

	// derivative of centroid count wrt x
	numeric::xyzVector< Real > const dcentroids6_dx  ( dcenlist.fcen6(position) );
	numeric::xyzVector< Real > const dcentroids10_dx ( dcenlist.fcen10(position) );
	numeric::xyzVector< Real > const dcentroids12_dx ( dcenlist.fcen12(position) );

	// derivative of score w.r.t centroid count
	Real const fcen6  ( cenlist.fcen6(position) );
	Real const fcen10 ( cenlist.fcen10(position) );
	Real const fcen12 ( cenlist.fcen12(position) );

	d_env_score = env_[ rsd.aa() ].dfunc( fcen10 ) * dcentroids10_dx;
	d_cb_score6  = cbeta6_.dfunc( fcen6 ) * dcentroids6_dx;
	d_cb_score12 = cbeta12_.dfunc( fcen12 ) * dcentroids12_dx;

	// change in others' context:
	Vector atom_x = rsd.atom( rsd.nbr_atom() ).xyz();
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node(position)->const_edge_list_begin(),
			irue = energy_graph.get_node(position)->const_edge_list_end();
			iru != irue; ++iru ) {
		EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
		Size const j( edge->get_other_ind(position) );
		conformation::Residue const & rsd2 ( pose.residue(j) );
		if ( !rsd2.is_protein() ) continue;
		numeric::xyzVector<Real> cenvec = rsd2.atom( rsd2.nbr_atom() ).xyz() - atom_x;
		Real const cendist = cenvec.length_squared();

		if ( cendist > cen_dist_cutoff_12_pad ) continue;

		Real x = sqrt(cendist);
		numeric::xyzVector<Real> gradx = cenvec/x;

		Real d6,d10,d12, z;
		z =  SIGMOID_SLOPE*(x-6);
		if ( z>-20 && z<20 ) {
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

		d_env_score  += env_[ rsd2.aa() ].dfunc( cenlist.fcen10(j) ) * (d10*gradx);
		d_cb_score6  += cbeta6_.dfunc( cenlist.fcen6(j) ) * (d6*gradx);
		d_cb_score12 += cbeta12_.dfunc( cenlist.fcen12(j) ) * ((d12-d6)*gradx);
	}
}


void
SmoothEnvPairPotential::evaluate_pair_and_cenpack_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist2,
	Real & pair_contribution,
	Real & cenpack_contribution
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	pair_contribution    = 0.0;
	cenpack_contribution = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	//CAR  no pair score if a disulfide
	if ( aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
			rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
			rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ <= 8 ) return;

	Real cendist = sqrt(cendist2);
	SmoothScoreTermCoeffs const & currPair = pair_[std::min(aa1,aa2)][std::max(aa1,aa2)];
	pair_contribution = currPair.func( cendist );

	// smooth pair by an additional sigmoid so pair_score->0 as dist->inf
	Real z = SIGMOID_SLOPE*(cendist-10.5);
	if ( z>20 ) {
		pair_contribution = 0;
	} else if ( z>-20 ) {
		pair_contribution *= 1 / (1+exp(z));
	}
	cenpack_contribution = cenpack_.func( cendist * 10 + 0.5 );
}


void
SmoothEnvPairPotential::evaluate_pair_and_cenpack_deriv(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist2,
	Real & d_pair,
	Real & d_cenpack
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	d_pair = 0.0;
	d_cenpack = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	if ( aa1 == chemical::aa_cys && aa2 == chemical::aa_cys &&
			rsd1.is_bonded( rsd2 ) && rsd1.polymeric_sequence_distance( rsd2 ) > 1 &&
			rsd1.has_variant_type( chemical::DISULFIDE ) && rsd2.has_variant_type( chemical::DISULFIDE ) ) return;
	if ( rsd1.polymeric_sequence_distance( rsd2 ) <= 8 ) return;

	// pair
	Real cendist = sqrt(cendist2);
	SmoothScoreTermCoeffs const & currPair = pair_[std::min(aa1,aa2)][std::max(aa1,aa2)];
	Real pair = currPair.func( cendist );
	Real z = SIGMOID_SLOPE*(cendist-10.5), sigmoid, d_sigmoid;
	if ( z<-20 ) {
		sigmoid=1;
		d_sigmoid=0;
	} else if ( z>20 ) {
		sigmoid=0;
		d_sigmoid =0;
	} else {
		Real e = exp(z);
		Real d = (1 + e);
		sigmoid = 1/d;
		d_sigmoid = -SIGMOID_SLOPE*e / (d*d);
	}
	d_pair = sigmoid * currPair.dfunc( cendist ) + pair * d_sigmoid;

	// cenpack
	d_cenpack = 10 * cenpack_.dfunc( cendist * 10 + 0.5 );
	// 10 x because cenpack_smooth was fit to original function on 0.1A grid
}


/// @details Pose must already contain a cenlist object or this method will fail.
SigmoidWeightedCenList< Real > const &
SmoothEnvPairPotential::cenlist_from_pose( pose::Pose const & pose ) const {
	using namespace core::pose::datacache;
	return static_cast< SigmoidWeightedCenList< Real > const & >
		( pose.data().get( CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
SigmoidWeightedCenList< Real > &
SmoothEnvPairPotential::nonconst_cenlist_from_pose( pose::Pose & pose ) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ) ) {
		return static_cast< SigmoidWeightedCenList< Real > & >
			( pose.data().get( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST ));
	}
	// else
	SigmoidWeightedCenListRealOP cenlist( new SigmoidWeightedCenList< Real > );
	pose.data().set( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_CEN_LIST, cenlist );
	return *cenlist;
}

/// @details Pose must already contain a cenlist object or this method will fail.
SigmoidWeightedCenList< numeric::xyzVector< Real > > const &
SmoothEnvPairPotential::dcenlist_from_pose( pose::Pose const & pose ) const {
	using namespace core::pose::datacache;
	return static_cast< SigmoidWeightedCenList< numeric::xyzVector< Real > > const & >
		( pose.data().get( CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
SigmoidWeightedCenList< numeric::xyzVector< Real > > &
SmoothEnvPairPotential::nonconst_dcenlist_from_pose( pose::Pose & pose ) const {
	if ( pose.data().has( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ) ) {
		return static_cast< SigmoidWeightedCenList< numeric::xyzVector< Real > > & >
			( pose.data().get( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST ));
	}
	// else
	SigmoidWeightedCenListVectorOP cenlist( new SigmoidWeightedCenList< numeric::xyzVector< Real > > );
	pose.data().set( core::pose::datacache::CacheableDataType::SIGMOID_WEIGHTED_D_CEN_LIST, cenlist );
	return *cenlist;
}


}
}
