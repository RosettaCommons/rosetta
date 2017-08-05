/*
* MeanField.cc
*
*  Created on: Jan 15, 2013
*      Author: arubenstein
*/
/// @brief   Pilot application source code for comparing specificity profiles using various metrics.
/// @author  arubenstein

// Unit headers

// Project headers
#include <devel/init.hh>
#include <protocols/mean_field/AAMatrix.hh>
#include <protocols/mean_field/jagged_array.hh>
#include <protocols/mean_field/jagged_array.functions.hh>
#include <protocols/mean_field/AAProb.hh>

// util
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh> // file_exists
#include <iostream>

// option key includes

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/mean_field.OptionKeys.gen.hh>
#include <basic/options/option.hh>

static THREAD_LOCAL basic::Tracer TR("apps.pilot.arubenstein.spec_prof_compare") ;

utility::vector1< std::string >
tokenize_line( std::istream & inputstream )
{

	utility::vector1< std::string > tokens;
	std::string input_line;
	std::getline( inputstream, input_line );

	unsigned int llength = input_line.size();
	unsigned int processing = 0;
	unsigned int token_start_char = 0;

	while ( processing < llength ) {
		if ( std::isspace( input_line[ processing ] ) ) { // == ' ') {
			if ( !( std::isspace(input_line[ token_start_char ] ) ) ) { // != ' ' ) {
				std::string token = input_line.substr( token_start_char, processing - token_start_char);
				tokens.push_back( token );
			}
			++processing;
			token_start_char = processing;
		} else {
			++processing;
		}
	}
	if ( processing != token_start_char ) { // last token on the line
		std::string token = input_line.substr( token_start_char, processing - token_start_char + 1 );
		tokens.push_back( token );
	}

	return tokens;
}

int main(int argc, char *argv[])
{


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	devel::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	protocols::mean_field::AAMatrix specprof;
	std::istringstream sp;
	if ( basic::options::option[ basic::options::OptionKeys::mean_field::spec_profile ].active() ) {
		std::string filename = basic::options::option[ basic::options::OptionKeys::mean_field::spec_profile];
		std::string sp_str;
		utility::io::izstream file( filename );
		if ( !file ) {
			TR.Info << "File:" << filename << " not found!\n";
		} else {
			utility::slurp( file, sp_str );
			sp.str( sp_str );
			specprof.build_aa_matrix( sp );
			specprof.show( TR.Info );
		}
	}
	typedef utility::vector1< utility::file::FileName > Filenames;

	utility::vector1 < protocols::mean_field::AAMatrix > specprofiles;
	utility::vector1 < std::string > specprofilenames;
	if ( basic::options::option[ basic::options::OptionKeys::mean_field::bb_list ].active() ) {
		Filenames listnames( basic::options::option[ basic::options::OptionKeys::mean_field::bb_list ]().vector() );
		for ( Filenames::const_iterator filename( listnames.begin() );
				filename != listnames.end(); ++filename ) {
			utility::io::izstream list( (*filename).name().c_str() );
			while ( list ) {
				std::string pdbname;
				list >> pdbname;
				if ( pdbname != "" ) {
					specprofilenames.push_back( pdbname );
					utility::io::izstream file( pdbname );
					if ( !file ) {
						TR.Info << "File:" << pdbname << " not found!\n";
					} else {
						std::istringstream sp;
						std::string sp_str;

						utility::slurp( file, sp_str );
						sp.str( sp_str );
						specprofiles.push_back( protocols::mean_field::AAMatrix( sp ) );

					}
				}
			}
		}

	}


	utility::vector1 < std::string >::const_iterator spn( specprofilenames.begin() );
	for ( utility::vector1 < protocols::mean_field::AAMatrix >::const_iterator sp( specprofiles.begin() );
			sp != specprofiles.end() && spn != specprofilenames.end(); ++sp, ++spn ) {
		utility::vector1 < core::Real > cosine_dist = specprof.cosine_distance( *sp );
		utility::vector1 < core::Real > frob_dist = specprof.frob_distance( *sp );
		utility::vector1 < core::Real > ave_abs_diff = specprof.ave_abs_diff( *sp );

		utility::vector1 < utility::vector1 < core::Real > > vecs;
		vecs.push_back( cosine_dist );
		vecs.push_back( frob_dist );
		vecs.push_back( ave_abs_diff );

		TR.Info << *spn << std::endl;
		for ( utility::vector1 < utility::vector1 < core::Real > >::const_iterator vec = vecs.begin(); vec != vecs.end();
				++vec ) {
			utility::vector1 < core::Real > v = *vec;
			for ( utility::vector1 < core::Real >::const_iterator dist = v.begin(); dist != v.end();
					++dist ) {
				TR.Info << *dist << "\t";
			}
			TR.Info << std::endl;
		}
	}


	if ( basic::options::option[ basic::options::OptionKeys::mean_field::bb_boltz_probs ].active() ) {

		TR.Info << "Just before reading in bb_boltz_probs" << std::endl;

		//input bb_boltz_probs
		protocols::mean_field::jagged_array < core::Real > bb_boltz_probs( specprofiles[1].size(), utility::vector1 < core::Real > ( specprofiles.size(), 0.0 ) );
		utility::vector1 < std::string > filenames;
		TR.Info << "Just before actual reading" << std::endl;

		std::istringstream boltz_p;
		std::string boltz_filename = basic::options::option[ basic::options::OptionKeys::mean_field::bb_boltz_probs ];
		std::string boltz_str;
		utility::io::izstream file( boltz_filename );
		if ( !file ) {
			TR.Info << "File:" << boltz_filename << " not found!\n";
		} else {
			utility::slurp( file, boltz_str );
			boltz_p.str( boltz_str );
			core::Size lineno = 1;
			utility::vector1 < std::string > tokens( tokenize_line ( boltz_p) );
			TR.Info << "Just before while loop" << std::endl;

			while ( ! tokens.empty() ) {

				filenames.push_back( tokens[ 1 ] );
				for ( core::Size i = 2; i <= tokens.size(); ++i ) {
					TR.Info << "in for loop " << i << std::endl;

					std::istringstream token_s( tokens[ i ] );
					core::Real prob;
					token_s >> prob;
					bb_boltz_probs[i-1][lineno] = prob;
				}
				tokens = tokenize_line ( boltz_p) ;
				++lineno;

			}
		}
		TR.Info << "Just before deleting empty pwms" << std::endl;

		//delete empty pwms
		for ( core::Size bb = 1; bb <= specprofiles.size(); ++bb ) {
			if ( specprofiles[bb].empty() ) {
				TR.Info << "deleting empty pwms" << bb << std::endl;
				specprofiles.erase( specprofiles.begin() + bb - 1);
				for ( core::Size pos = 1; pos <= specprofiles[bb].size(); ++pos ) {
					bb_boltz_probs[pos].erase( bb_boltz_probs[pos].begin() + bb - 1 );
				}
				--bb;
			}

		}

		TR.Info << "Just before renormalizing bb boltz probs" << std::endl;

		//renormalize bb_boltz_probs in case any backbones were deleted
		utility::vector1 < core::Real > bb_boltz_probs_totals = bb_boltz_probs.get_totals_columns();
		bb_boltz_probs /= bb_boltz_probs_totals;

		//average pwms of structures
		protocols::mean_field::AAMatrix exp_spec_profile;
		TR.Info << "Just before averagomg pwms" << std::endl;

		for ( core::Size pos = 1; pos <= specprofiles[1].size(); ++pos ) {
			exp_spec_profile.push_back( utility::vector1 < protocols::mean_field::AAProb> ( specprofiles[1][pos].size() ) );
			for ( core::Size bb = 1; bb <= specprofiles.size(); ++bb ) {
				//    if ( specprofilenames[bb] != filenames[bb] )
				//    {
				//     TR.Fatal << "Names do not match" << std::endl;
				//    }
				for ( core::Size aa = 1; aa <= specprofiles[bb][pos].size(); ++aa ) {
					//This RotProb wasn't yet initialized so initialize it to value of current RotProb multiplied by the bb probability
					if ( exp_spec_profile[pos][aa].pos() == 0 ) {
						exp_spec_profile[pos][aa] = protocols::mean_field::AAProb( specprofiles[bb][pos][aa] ) * bb_boltz_probs[pos][bb];
						//exp_spec_profile_[pos][aa] = AAProb( spec_profiles_[bb][pos][aa] ) * bb_boltz_probs_per_bb()[bb];
					} else {
						//Add second RotProb multiplied by its bb probability
						exp_spec_profile[pos][aa] += specprofiles[bb][pos][aa] * bb_boltz_probs[pos][bb];
						//exp_spec_profile_[pos][aa] += spec_profiles_[bb][pos][aa] * bb_boltz_probs_per_bb()[bb];

					}
				} // loops through rot
			} //loops through bb
		} //loops through pos

		TR.Info << "Just before printing distances" << std::endl;

		utility::vector1 < core::Real > cosine_dist = specprof.cosine_distance( exp_spec_profile );
		utility::vector1 < core::Real > frob_dist = specprof.frob_distance( exp_spec_profile );
		utility::vector1 < core::Real > ave_abs_diff = specprof.ave_abs_diff( exp_spec_profile );

		utility::vector1 < utility::vector1 < core::Real > > vecs;
		vecs.push_back( cosine_dist );
		vecs.push_back( frob_dist );
		vecs.push_back( ave_abs_diff );

		TR.Info << "Expected Spec Profile using Following Boltzmann Probs" << std::endl;
		bb_boltz_probs.show( TR.Info );

		for ( utility::vector1 < utility::vector1 < core::Real > >::const_iterator vec = vecs.begin(); vec != vecs.end();
				++vec ) {
			utility::vector1 < core::Real > v = *vec;
			for ( utility::vector1 < core::Real >::const_iterator dist = v.begin(); dist != v.end();
					++dist ) {
				TR.Info << *dist << "\t";
			}
			TR.Info << std::endl;
		}

		if ( basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ].active() ) {
			std::string filename = basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ];
			filename += ".transfac";
			exp_spec_profile.dump_transfac( filename );

		}

	}



}


