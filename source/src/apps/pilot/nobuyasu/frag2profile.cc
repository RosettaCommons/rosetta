// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief
/// @author Nobuyasu Koga

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/scoring/rms_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/excn/Exceptions.hh>

// Utility Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <devel/init.hh>
#include <ObjexxFCL/format.hh>
#include <fstream>

static basic::Tracer TR("frag2profile");

typedef core::Size Size;
typedef std::string String;

using core::fragment::FragmentIO;
using core::fragment::FrameList;
using core::fragment::FragSetOP;
typedef std::map< Size, Size >::iterator M_iter;


class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, f )
OPT_KEY( File, s )
OPT_KEY( File, o )
OPT_KEY( Real, rmsd )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( f, "fragment file", "" );
	NEW_OPT( s, "pdb file name", "" );
	NEW_OPT( o, "output filename", "output" );
	NEW_OPT( rmsd, "output filename", 1.0 );
}


bool compare( M_iter const & a, M_iter const & b )
{
	return ( a->second > b->second );
}


int
main( int argc, char * argv [] )
{
	try{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ThisApplication::register_options();
	devel::init(argc, argv);

	// output file //////////////////////////////////////////////////////////////////////////////////////////////////
	std::ofstream output_;
	std::ostringstream filename;
	if( option[ o ].active() ) {
		filename <<  option[ o ]() << ".profile";
	} else {
		filename <<  "default.profile";
	}
	output_.open( filename.str().c_str() ,std::ios::out );

	std::ofstream resfile_;
	filename.str("");
	if( option[ o ].active() ) {
		filename <<  option[ o ]() << ".resfile";
	} else {
		filename <<  "default.resfile";
	}
	resfile_.open( filename.str().c_str() ,std::ios::out );

	std::ofstream fasta_;
	filename.str("");
	if( option[ o ].active() ) {
		filename <<  option[ o ]() << ".fasta";
	} else {
		filename <<  "defa.fasta";
	}
	fasta_.open( filename.str().c_str() ,std::ios::out );


	// rmsd //////////////////////////////////////////////////////////////////////////////////////////////////
	core::Real rmsd_cutoff_ = option[ rmsd ]();

	// read structure //////////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, option[ s ]() );

	// read fragments //////////////////////////////////////////////////////////////////////////////////////////////////
	using core::fragment::FragSetOP;
	using core::fragment::FragData;
	using core::fragment::FrameList;
	using core::fragment::FragmentIO;
	using core::fragment::FrameIterator;
	using core::scoring::CA_rmsd;

	//
	FragSetOP fragset;
	FrameList frames;
	fragset = FragmentIO().read_data( option[ f ]() );

	//
	core::pose::Pose test_pose( pose );
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	core::util::switch_to_residue_type_set( test_pose, core::chemical::CENTROID );

	// intialize //////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1<Size> > freq( pose.total_residue(), utility::vector1<Size>( 20, 0 ) );
	for( Size i=1; i<=pose.total_residue(); i++ ) {
		for( Size j=1; j<=20; j++ ) {
			freq[ i ][ j ] = 0;
		}
	}

	for ( FrameIterator frame = fragset->begin(); frame != fragset->end(); ++frame ) {

		Size const start ( frame->start() );
		if( ( start + ( frame->length() - 1 ) ) > pose.total_residue() ) continue;
		runtime_assert( start <= pose.total_residue() );
		
		for ( Size i=1; i<=frame->nr_frags(); i++ ) {
			// insert fragment
			
			frame->apply( i, test_pose );
			// calc rmsd
			core::Real rmsd = CA_rmsd( pose, test_pose, start, start + frame->length() - 1 );

			
			if( rmsd <= rmsd_cutoff_ ) {
				FragData fragdat = frame->fragment( i );
				for( Size j=1; j<=fragdat.size(); j++ ) {
					
					Size res = start + j - 1;
					freq[ res ][ core::chemical::aa_from_oneletter_code( fragdat.sequence( j ) ) ] ++;
					
					// std::cout << j << " " << res << " " << fragdat.sequence( j ) << " " << fragdat.sequence() << std::endl;
					
				}
			}
		}
	} // FrameIterator
	
	/// output //////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Size > total;
	total.resize( pose.total_residue() );
	for( Size i=1; i<=pose.total_residue(); i++ ) {
		total[ i ] = 0;
		for( Size j=1; j<=20; j++ ) {
			total[ i ] += freq[ i ][ j ];
		}
	}

	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::RJ;

	output_ << " RES" << " ";
	for( Size i=1; i<=20; i++ ) {
		output_ << RJ( 4, core::chemical::name_from_aa( core::chemical::AA( i ) ) ) << " ";
	}
	output_ << std::endl;

	for( Size i=1; i<=pose.total_residue(); i++ ) {
		output_ << I( 4, i ) << " ";
		for( Size j=1; j<=20; j++ ) {
			output_ << I( 4, freq[ i ][ j ] ) << " ";
		}
		output_ << std::endl;
	}

	output_ << std::endl;
	output_ << " RES" << " ";
	for( Size i=1; i<=20; i++ ) {
		output_ << RJ( 5, core::chemical::name_from_aa( core::chemical::AA( i ) ) ) << " ";
	}
	output_ << std::endl;

	for( Size i=1; i<=pose.total_residue(); i++ ) {
		output_ << I( 4, i ) << " ";
		for( Size j=1; j<=20; j++ ) {
			if( total[ i ] == 0 ) {
				output_ << F( 5, 2, 0.0 ) << " ";
			} else {
				output_ << F( 5, 2, core::Real( freq[ i ][ j ] )/core::Real( total[ i ] ) ) << " ";
			}
		}
		output_ << std::endl;
	}
	output_ << std::endl;

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size l( 0 );
	fasta_ << ">ideal.L" << std::endl;
	resfile_ << " start" << std::endl;
	for( Size i=1; i<=pose.total_residue(); i++ ) {

		std::map< Size, Size > mfreq;
		for( Size j=1; j<=20; j++ ) {
			mfreq.insert( std::map< Size, Size >::value_type( j, freq[ i ][ j ] ) );
		}

		utility::vector1< M_iter > sorted;
		M_iter it = mfreq.begin();
		while( it != mfreq.end() ) {
			sorted.push_back( it );
			++it;
		}

		std::sort( sorted.begin(), sorted.end(), compare );

		Size k( 0 );
		output_ << I( 4, i ) << " " << RJ( 2, pose.sequence().substr( i-1, 1 ) ) << "," ;
		for( utility::vector1< M_iter >::const_iterator it=sorted.begin(); it != sorted.end(); it++ ) {
			M_iter itt( *it );
			output_ << RJ( 2, core::chemical::oneletter_code_from_aa( core::chemical::AA( itt->first ) ) )
							<< I( 4, itt->second ) << ",";
			k++;
			if( k == 1 ) {
				l++;
				resfile_ << I( 7, i ) << " A PIKAA "
								 << core::chemical::oneletter_code_from_aa( core::chemical::AA( itt->first ) ) << std::endl;
				fasta_ << core::chemical::oneletter_code_from_aa( core::chemical::AA( itt->first ) );
				if( l == 60 ) {
					fasta_ << std::endl;
					l = 0;
				}
			}
		}
		output_ << std::endl;
	}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;
}

