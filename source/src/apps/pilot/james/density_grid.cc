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
/// @author James Thompson

#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <devel/init.hh>
#include <basic/database/open.hh>
#include <basic/prof.hh>
#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <numeric/kdtree/util.hh>
#include <numeric/kdtree/KDTree.hh>
#include <numeric/kdtree/KDPoint.hh>
#include <numeric/kdtree/KDPointList.hh>
#include <numeric/kdtree/WrappedReal.hh>

#include <apps/pilot/james/james_util.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

// silly typedefs
using namespace numeric::kdtree;
typedef WrappedPrimitive< std::pair< char, char > > SS_Prediction;
typedef utility::pointer::owning_ptr< SS_Prediction > SS_PredictionOP;

static THREAD_LOCAL basic::Tracer tr( "density_grid" );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char* argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );

	using namespace core::id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace ObjexxFCL::format;
	using namespace numeric::kdtree;
	using namespace core::sequence;

	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;
	using utility::file::FileName;

	vector1< string > infiles = option[ in::file::s ]();
	vector1< string >::const_iterator iter, end;

	Size length           (  5   );
	if ( option[ james::debug ]() ) {
		length = 4;
	}
	string const atom_name( "CA" );
	Real const max_dist   ( 10.0 );
	Size const min_seqsep ( 10   );

	vector1< vector1< Real > > query_points;
	vector1< vector1< Real > > db_points;
	vector1< utility::pointer::ReferenceCountOP > errors;

	string closest_line;
	for ( iter = infiles.begin(), end = infiles.end(); iter != end; ++iter ) {
		utility::io::izstream input( *iter );

		string line;
		getline(input,line);
		while ( getline(input,line) ) {
			if ( option[ james::debug ]() ) {
				string q_id, t_id, dummy;
				Size q_resi, t_resi;
				Real e_value, blosum_score, burial_score, dgap_score;
				char ss_pred, ss_nat;

				std::istringstream tokens( line );
				tokens
					>> e_value >> q_resi >> t_resi >> blosum_score >> dgap_score >> burial_score
					>> ss_pred >> ss_nat
					>> q_id >> t_id;
				vector1< Real > vec1;
				vec1.resize( length );
				vec1[1] = std::max( e_value, -400.0 );
				vec1[2] = blosum_score;
				vec1[3] = burial_score;
				vec1[4] = dgap_score;
				db_points.push_back( vec1 );
				//errors.push_back( new WrappedChar( ) );
				errors.push_back( new SS_Prediction( std::make_pair( ss_pred, ss_nat ) ) );
			} else {
				string q_id, t_id;
				Size resi, resj;
				Real e_value, blosum_avg, burial_avg, dgap_avg, query_cover,
					nat_dist, pred_dist;

				std::istringstream tokens( line );
				tokens
					>> q_id >> t_id >> resi >> resj
					>> e_value >> blosum_avg >> burial_avg >> dgap_avg >> query_cover
					>> nat_dist >> pred_dist;

				vector1< Real > vec1;
				vec1.resize( length );
				vec1[1] = std::max( e_value, -400.0 );
				vec1[2] = blosum_avg;
				vec1[3] = burial_avg;
				vec1[4] = dgap_avg;
				vec1[5] = pred_dist;
				db_points.push_back( vec1 );
				errors.push_back( new WrappedReal( std::abs( pred_dist - nat_dist ) ) );
			}
		}
	} // for infiles

	using namespace basic;

	PROF_START( KDTREE_CONSTRUCT );
	HyperRectangleOP bounds( get_percentile_bounds( db_points ) );

	transform_percentile( db_points );
	KDTree tree( db_points, errors );
	PROF_STOP( KDTREE_CONSTRUCT );

	vector1< DistancePrediction > predictions;
	vector1< FileName > aln_files( option[ in::file::alignment ]() );
	for ( vector1< FileName >::const_iterator file_it = aln_files.begin(),
				file_end = aln_files.end(); file_it != file_end; ++file_it
	) {
		vector1< SequenceAlignment > alns(
			read_aln( option[ cm::aln_format ](), *file_it )
		);
		tr.Debug << "read " << alns.size() << " alignments from " << *file_it
			<< std::endl;

		typedef vector1< SequenceAlignment >::iterator align_iter;
		for ( align_iter it = alns.begin(), end = alns.end(); it != end; ++it ) {
			Real ln_e_value = it->score();
			if ( it->identities() == it->length() ) continue;
			tr.Debug << "alignment between " << it->sequence(1)->id() << " and "
				<< it->sequence(2)->id() << std::endl;

			// get the template structure
			core::pose::PoseOP templ_pose = get_pose_by_id( it->sequence(2)->id() );
			SequenceMapping t_aln_to_pdb = map_seq1_seq2(
				it->sequence(2),
				new Sequence(
					templ_pose->sequence(), it->sequence(2)->id() + ".pdb", 1
				)
			);

			vector1< Real > scores ( calc_blosum_scores( *it ) );
			vector1< Real > dgaps( calc_dgaps( *it ) );
			vector1< int > t_burial( calculate_burial( *templ_pose ) );

			if ( option[ james::debug ]() ) {
				vector1< char > templ_ss( get_ss( *templ_pose ) );
				for ( Size ii = 1; ii <= it->length(); ++ii ) {
					if ( it->is_gapped( ii ) ) continue;
					Real const dgap_avg ( dgaps[ ii ] );
					//Size const q_resi_aln( it->sequence(1)->resnum(ii) );
					Size const t_resi_aln( it->sequence(2)->resnum(ii) );
					Size const t_resi_pdb( t_aln_to_pdb[ t_resi_aln ] ); // pdb residue index
					if ( t_resi_pdb == 0 ) continue; // skip residues not in template PDB

					vector1< Real > vec1( length, 0.0 );
					vec1[1] = ln_e_value;
					vec1[2] = scores[ii];
					vec1[3] = t_burial[ii];
					vec1[4] = dgap_avg;
					//char const t_ss( templ_ss[ii] );
				}
			} else {
				DistanceMatrix t_dists ( *templ_pose, atom_name );

				for ( Size ii = 1; ii <= it->length(); ++ii ) {
					if ( it->is_gapped( ii ) ) continue;
					Size const q_resi_aln( it->sequence(1)->resnum(ii) );
					Size const t_resi_aln( it->sequence(2)->resnum(ii) );
					Size const t_resi_pdb( t_aln_to_pdb[ t_resi_aln ] ); // pdb residue index

					for ( Size jj = ii + min_seqsep; jj <= it->length(); ++jj ) {
						if ( it->is_gapped( jj ) ) continue;
						Size const q_resj_aln( it->sequence(1)->resnum(jj) );
						Size const t_resj_aln( it->sequence(2)->resnum(jj) );
						Size const t_resj_pdb( t_aln_to_pdb[ t_resj_aln ] ); // pdb residue index
						Real const pred_dist(
							t_dists.distance( t_resi_pdb, t_resj_pdb )
						);
						Real const dgap_avg ( ( dgaps[ ii ] + dgaps[ jj ] ) / 2 );
						Real const burial_avg( ( t_burial[t_resi_pdb] + t_burial[t_resj_pdb] ) / 2 );

						if ( t_resi_pdb == 0 ) continue; // skip residues not in template PDB
						if ( pred_dist > max_dist ) continue;

						vector1< Real > vec1( length, 0.0 );
						vec1[1] = ln_e_value;
						vec1[2] = scores[ii];
						vec1[3] = burial_avg;
						vec1[4] = dgap_avg;
						vec1[5] = pred_dist;

						transform_percentile_single_pt( vec1, bounds );
						query_points.push_back( vec1 );

						// hack!
						DistancePrediction prediction(
							AtomID( 2, q_resi_aln ),
							AtomID( 2, q_resj_aln ),
							pred_dist,
							vec1
						);
						predictions.push_back( prediction );
					} // resj
				} // for resi
			} // debug
		} // alns
	} // aln_files

	Size const nn( option[ cm::nn ]() );
	PROF_START( KDTREE_SEARCH );
	//std::map< Size, std::map< Size, Real > > sdevs;
	typedef vector1< DistancePrediction >::const_iterator pred_iter;
	for ( pred_iter pred_it = predictions.begin(), pred_end = predictions.end();
				pred_it != pred_end; ++pred_it
	) {
		utility::vector1< core::Real > vars = pred_it->vars();
		KDPointList neighbors = nearest_neighbors( tree, vars, nn );
		vector1< KDPointOP > values = neighbors.sorted_values();

		Real total( 0.0 );
		for ( vector1< KDPointOP >::const_iterator
					it = values.begin(), end = values.end(); it != end; ++it
		) {
			WrappedRealOP err = dynamic_cast< WrappedReal * > ( (*it)->data()() );
			//total += std::abs( err->val() );
			total += err->val();
			//total_sq += err->val() * err->val();
			//std::cout << (*it)->location() << " -> " << err->val() << std::endl;
		}

		//Real const sdev( sqrt( total_sq / values.size() ) );
		Real const sdev( total / values.size() );
		Real pred_dist = pred_it->predicted_distance();
		std::cout << "AtomPair "
			//<< pred_it->atom_i().atomno() << " " << pred_it->atom_i().rsd() << " "
			//<< pred_it->atom_j().atomno() << " " << pred_it->atom_j().rsd() << " "
			<< " CA " << pred_it->atom_i().rsd()
			<< " CA " << pred_it->atom_j().rsd()
			<< " GAUSSIANFUNC " << pred_dist
			<< " " << sdev
			<< std::endl;
	}
	PROF_STOP( KDTREE_SEARCH );

	prof_show();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
