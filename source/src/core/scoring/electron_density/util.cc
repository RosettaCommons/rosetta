// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/SplineInterp.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/options/option.hh>


#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>  // b factors
#include <core/pose/symmetry/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <iostream>

// option key includes
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/patterson.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace electron_density {

static basic::Tracer TR( "protocols.scoring.electron_density.util" );

/// @brief read density weights from the cmd line into the scorefunction
void add_dens_scores_from_cmdline_to_scorefxn( core::scoring::ScoreFunction &scorefxn ) {
	using namespace basic::options;

	if ( option[ OptionKeys::edensity::fastdens_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_fast,
			option[ OptionKeys::edensity::fastdens_wt ]() );
	}
	if ( option[ OptionKeys::edensity::sliding_window_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_window,
			option[ OptionKeys::edensity::sliding_window_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_ca_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_ca,
			option[ OptionKeys::edensity::whole_structure_ca_wt ]() );
	}
	if ( option[ OptionKeys::edensity::whole_structure_allatom_wt ].user() ) {
		scorefxn.set_weight( core::scoring::elec_dens_whole_structure_allatom,
			option[ OptionKeys::edensity::whole_structure_allatom_wt ]() );
	}
}

bool pose_has_nonzero_Bs( core::pose::Pose const & pose ) {
	if ( !pose.pdb_info() ) return false;

	// to save time check only the first atom of each residue
	for ( uint i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue const & rsd_i ( pose.residue(i) );
		if ( rsd_i.aa() == core::chemical::aa_vrt ) continue;
		Real B = pose.pdb_info()->temperature( i, 1 );
		if ( B > 0 ) {
			return true;
		}
	}
	return false;
}

bool pose_has_nonzero_Bs( poseCoords const & pose ) {
	for ( int i=1 ; i<=(int)pose.size(); ++i ) {
		if ( pose[i].B_ > 0 ) {
			return true;
		}
	}
	return false;
}


/// @brief spline interpolation with periodic boundaries
core::Real interp_spline(
	ObjexxFCL::FArray3D< double > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX,
	bool mirrored
) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp3(&coeffs[0], dims, pt, mirrored);
	return retval;
}

/// @brief spline interpolation with periodic boundaries
numeric::xyzVector<core::Real> interp_dspline(
	ObjexxFCL::FArray3D< double > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX,
	bool mirrored
) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[3] = { 0,0,0 };
	SplineInterp::grad3(&grad[0], &coeffs[0], dims, pt, mirrored);
	return numeric::xyzVector<core::Real>(grad[2],grad[1],grad[0]);
}

/// @brief spline interpolation with periodic boundaries, single precision
core::Real interp_spline(
	ObjexxFCL::FArray3D< float > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX,
	bool mirrored
) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp3< float >(&coeffs[0], dims, pt, mirrored);
	return retval;
}

/// @brief spline interpolation with periodic boundaries, single precision
numeric::xyzVector<core::Real> interp_dspline(
	ObjexxFCL::FArray3D< float > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX,
	bool mirrored
) {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = { idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[3] = { 0,0,0 };
	SplineInterp::grad3< float >(&grad[0], &coeffs[0], dims, pt, mirrored);
	return numeric::xyzVector<core::Real>(grad[2],grad[1],grad[0]);
}

void spline_coeffs(
	ObjexxFCL::FArray3D< double > const & data ,
	ObjexxFCL::FArray3D< double > & coeffs,
	bool mirrored
) {
	int dims[3] = { data.u3(), data.u2(), data.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients3( const_cast<double*>(&coeffs[0]) , dims, mirrored );  // external code wants nonconst even though array is unchanged
}

void spline_coeffs(
	ObjexxFCL::FArray3D< float > const & data ,
	ObjexxFCL::FArray3D< double > & coeffs,
	bool mirrored
) {
	int N = data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray3D< double > data_d(data.u1(),data.u2(),data.u3()) ;
	for ( int i=0; i<N; ++i ) {
		data_d[i] = (double)data[i];
	}
	spline_coeffs( data_d, coeffs, mirrored );
}

void conj_map_times(ObjexxFCL::FArray3D< std::complex<double> > & map_product, ObjexxFCL::FArray3D< std::complex<double> > const & mapA, ObjexxFCL::FArray3D< std::complex<double> > const & mapB) {
	debug_assert(mapA.u1() == mapB.u1());
	debug_assert(mapA.u2() == mapB.u2());
	debug_assert(mapA.u3() == mapB.u3());

	map_product.dimension(mapA.u1(), mapA.u2(), mapA.u3());
	for ( Size i=0; i < mapA.size(); i++ ) {
		map_product[i] = std::conj(mapA[i]) * mapB[i];
	}
}

//∑[A(x) * B(y-x)] over y
ObjexxFCL::FArray3D< double > convolute_maps( ObjexxFCL::FArray3D< double > const & mapA, ObjexxFCL::FArray3D< double > const & mapB) {

	ObjexxFCL::FArray3D< std::complex<double> > FmapA;
	numeric::fourier::fft3(mapA, FmapA);

	ObjexxFCL::FArray3D< std::complex<double> > FmapB;
	numeric::fourier::fft3(mapB, FmapB);

	ObjexxFCL::FArray3D< std::complex<double> > Fconv_map;
	conj_map_times(Fconv_map, FmapB, FmapA );

	ObjexxFCL::FArray3D< double > conv_map;
	numeric::fourier::ifft3(Fconv_map , conv_map);

	return conv_map;
}


/// 4D interpolants

/// @brief spline interpolation with periodic boundaries
core::Real interp_spline(
	ObjexxFCL::FArray4D< double > & coeffs ,
	core::Real slab,
	numeric::xyzVector< core::Real > const & idxX )
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[4] = { slab-1.0, idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real retval = SplineInterp::interp4(&coeffs[0], dims, pt);

	return retval;
}

/// @brief spline interpolation with periodic boundaries
void interp_dspline(
	ObjexxFCL::FArray4D< double > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX ,
	core::Real slab,
	numeric::xyzVector< core::Real > & gradX,
	core::Real & gradSlab )
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[4] = { slab-1.0, idxX[2]-1.0 , idxX[1]-1.0, idxX[0]-1.0 };
	core::Real grad[4] = { 0,0,0,0 };
	SplineInterp::grad4(&grad[0], &coeffs[0], dims, pt);
	gradX = numeric::xyzVector<core::Real>(grad[3],grad[2],grad[1]);
	gradSlab = grad[0];
}

void spline_coeffs(
	ObjexxFCL::FArray4D< double > const & data ,
	ObjexxFCL::FArray4D< double > & coeffs)
{
	int dims[4] = { coeffs.u4(), coeffs.u3(), coeffs.u2(), coeffs.u1() };
	coeffs = data;
	SplineInterp::compute_coefficients4( const_cast<double*>(&coeffs[0]) , dims );
}

void spline_coeffs(
	ObjexxFCL::FArray4D< float > const & data ,
	ObjexxFCL::FArray4D< double > & coeffs)
{
	int N = data.u4()*data.u3()*data.u2()*data.u1();
	ObjexxFCL::FArray4D< double > data_d(data.u1(),data.u2(),data.u3(),data.u4()) ;
	for ( int i=0; i<N; ++i ) {
		data_d[i] = (double)data[i];
	}
	spline_coeffs( data_d, coeffs );
}

void
calculate_density_nbr(
	pose::Pose & pose,
	std::map< Size, Real > & per_rsd_dens,
	std::map< Size, Real > & per_rsd_nbrdens,
	core::conformation::symmetry::SymmetryInfoCOP symminfo,
	bool mixed_sliding_window,
	Size sliding_window_size)
{

	using namespace core::scoring;
	// rescore pose
	core::Size nres_calculate = per_rsd_dens.size();

	core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
	edm.setScoreWindowContext( true );
	edm.setWindow( sliding_window_size );  // smoother to use 3-res window

	// score the pose
	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::elec_dens_window, 1.0 );

	if ( pose.is_fullatom() ) {
		myscore->set_weight( core::scoring::fa_rep, 10e-30 );
	} else {
		myscore->set_weight( core::scoring::vdw, 10e-30 );
	}

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	(*myscore)(pose);

	///////////////////////////////////////////////////////
	// get density correlation score from pose.size(), rather than nres_
	// get zscore for real-space density correlation scores

	core::Real rscc_sum=0, sq_rscc_sum=0;
	for ( Size r=1; r<=pose.size(); ++r ) { // loop over the entire pose

		if ( ! per_rsd_dens.count(r) ) continue;

		if ( mixed_sliding_window ) {
			if ( pose.residue(r).is_protein() ) {
				edm.setWindow( 3 );
			} else {
				edm.setWindow( 1 );
			}
		}

		Real dens_rscc = core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, symminfo , false);
		Size asymm_num_r = r;
		per_rsd_dens[asymm_num_r] = dens_rscc;
		rscc_sum    += per_rsd_dens[asymm_num_r];
		sq_rscc_sum += per_rsd_dens[asymm_num_r]*per_rsd_dens[asymm_num_r];

		TR.Trace << "res: " << asymm_num_r << " symmetric num: " << r << " dens_rscc: " << dens_rscc << std::endl;
	}
	// get mean and stdev for density rscc
	Real per_rsd_dens_mean  = rscc_sum/nres_calculate;
	Real per_rsd_dens_stdev = std::sqrt( sq_rscc_sum/nres_calculate - per_rsd_dens_mean*per_rsd_dens_mean );


	///////////////////////////////////////////////////////
	// for each residue, get neighbors from energy graph,
	// and calculate density-zscore from neighbors
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Size i=1; i <= pose.size(); ++i ) {

		if ( ! per_rsd_dens.count(i) ) continue;

		core::conformation::Residue const &rsd_i( pose.residue(i) );
		//if (rsd_i.name3()=="GLY") continue;

		Size asymm_num_i = i;
		Real i_dens_rscc = per_rsd_dens[asymm_num_i];
		Real sum    = i_dens_rscc;
		Real sq_sum = i_dens_rscc*i_dens_rscc;
		Size n_nbrs = 1;

		TR.Trace << "rsd: " << i << " " << rsd_i.name3() << " i_dens_rscc: " << i_dens_rscc << std::endl;

		// get density score per i
		for ( utility::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(i)->const_edge_list_begin(),
				irue = energy_graph.get_node(i)->const_edge_list_end();
				iru != irue; ++iru ) {

			auto const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_other_ind(i) );
			Size asymm_num_j = j;

			if ( core::pose::symmetry::is_symmetric( pose )  ) {
				if ( ! symminfo->bb_is_independent(j) ) {
					asymm_num_j = symminfo->bb_follows(j);
				}
			}

			conformation::Residue const & rsd_j ( pose.residue(j) );
			//if (rsd_j.name3()=="GLY") continue;

			Real dist = std::pow( edge->square_distance(), 0.5 );
			if ( dist <= 10.0 && per_rsd_dens.count(asymm_num_j) ) {
				Real j_dens_rscc = per_rsd_dens[asymm_num_j];
				sum    += j_dens_rscc;
				sq_sum += j_dens_rscc*j_dens_rscc;
				n_nbrs ++;

				TR.Trace << "computing energy graph: " << j
					//<< " dist: " <<  caled_dist
					<< " "                        << rsd_j.name3()
					<< " dist: "                 << dist
					<< " j_dens_rscc: "           << j_dens_rscc
					<< " sum: "                   << sum
					<< " sq_sum: "               << sq_sum
					<< " n_nbrs: "               << n_nbrs
					//<< " delta: "
					//<< dist-caled_dist
					<< std::endl;
			} else {
				continue;
			}
		} // res j

		TR.Trace << " sum: " << sum
			<< " sq_sum: " << sq_sum
			<< " n_nrbs: " << n_nbrs
			<< std::endl;

		Real nbrdens_mean  = sum/n_nbrs;
		Real nbrdens_stdev = std::sqrt( sq_sum/n_nbrs - nbrdens_mean*nbrdens_mean );

		// z-score for rscc of residue i to rscc of its neighbors
		Real i_nbrdens_zscore = (i_dens_rscc - nbrdens_mean) / nbrdens_stdev;

		TR.Trace << "N Neighbors: " << n_nbrs << std::endl;
		//Protect from NAN, use density zscore if no neighbors to estimate
		if ( n_nbrs <= 3 ) {
			i_nbrdens_zscore = (i_dens_rscc - per_rsd_dens_mean) / per_rsd_dens_stdev;
		}


		per_rsd_nbrdens[asymm_num_i] = i_nbrdens_zscore;

		// z-score for rscc of residue i to rscc of all residues
		Real i_per_rsd_dens_zscore = (i_dens_rscc-per_rsd_dens_mean)/per_rsd_dens_stdev;

		TR.Trace << "rsd: "              << asymm_num_i
			<< " rsn: "             << rsd_i.name3()
			<< " symm_rsd: "        << i
			<< " dens_rscc: "       << i_dens_rscc
			<< " nbrdens_mean: "    << nbrdens_mean
			<< " nbrdens_stdev: "   << nbrdens_stdev
			<< " i_nbr_zscore: "    << i_nbrdens_zscore
			<< " i_per_rsd_zscore: " << i_per_rsd_dens_zscore
			<< " nbrs: "            << n_nbrs
			<< std::endl;

	} // for i in range(pose.size())
}

void
calculate_rama(
	pose::Pose &pose,
	std::map< Size, Real > & rama,
	Size n_symm_subunit,
	Real weight /*=0.2*/
) {


	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( core::scoring::rama_prepro, weight );

	if ( pose.conformation().contains_carbohydrate_residues() ) {
		myscore->set_weight( core::scoring::sugar_bb, weight );
	}
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	(*myscore)(pose);

	for ( Size r=1; r<=pose.size(); ++r ) {
		if ( rama.count(r) == 0 ) continue;

		if ( pose.residue(r).is_protein() ) {
			rama[r] = weight*(pose.energies().residue_total_energies(r)[ core::scoring::rama_prepro ]/n_symm_subunit);
		} else if ( pose.residue(r).is_carbohydrate() ) {
			rama[r] = weight*(pose.energies().residue_total_energies(r)[ core::scoring::sugar_bb ]/n_symm_subunit);
		} else {
			rama[r] = 0.0;
		}
	}
}

void
calculate_geometry(
	pose::Pose & pose,
	std::map< Size, Real > & geometry,
	Size n_symm_subunit,
	Real weight  /*1.0*/)
{

	// clean the container
	calc_per_rsd_score( pose,
		core::scoring::cart_bonded_angle,
		geometry,
		n_symm_subunit,
		weight );

}

void
calc_per_rsd_score(
	pose::Pose &pose,
	scoring::ScoreType const & score_type,
	std::map< Size, Real > & per_rsd_score,
	Size n_symm_subunit,
	Real weight
) {

	runtime_assert( per_rsd_score.size());

	core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
	myscore->set_weight( score_type, weight );

	if ( core::pose::symmetry::is_symmetric(pose) ) {
		myscore = core::scoring::symmetry::symmetrize_scorefunction(*myscore);
	}

	(*myscore)(pose);
	//myscore->show_line( fragbias_tr, pose ); fragbias_tr << std::endl;

	for ( auto score_pair : per_rsd_score ) {

		Size r = score_pair.first;
		per_rsd_score[r] = weight*(pose.energies().residue_total_energies(r)[ score_type ]/n_symm_subunit);
	}
}



} // namespace constraints
} // namespace scoring
} // namespace core
