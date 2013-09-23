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

// Package headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/BB_Pos.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/fldsgn/NcontactsCalculator.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

// C++ header
#include <fstream>
#include <map>
#include <numeric/xyzVector.hh>
#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

typedef std::string String;


static basic::Tracer TR("pair_distance");

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, output )
OPT_KEY( File, blue )
OPT_KEY( Real, dist )
OPT_KEY( Real, radius )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( output, "output file name for log files", "pair_distance.out" );
	NEW_OPT( blue, "blueprint", "" );
	NEW_OPT( dist, "distance", 10.0 );
	NEW_OPT( radius, "pore radius", 2.0 );
}


/// local mover for testing purposes
class PairDistance : public protocols::moves::Mover {
public: // typedef


	typedef protocols::jd2::parser::BluePrint BluePrint;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;


public: // constructor/deconstructor


	PairDistance()
	{
		// output file
		std::ostringstream filename;
		filename <<  option[ output ]();
		output_.open( filename.str().c_str() ,std::ios::out );

		ssinput_ = false;
		if( option[ blue ].user() ) {
			blueprint_ = new BluePrint( option[ blue ]() );
			ssinfo_ = new SS_Info2( blueprint_->secstruct() );
			ssinput_ = true;
		} else {
			ssinfo_ = new SS_Info2;
		}

	}

	virtual ~PairDistance(){ output_.close(); }

	virtual std::string get_name() const { return ""; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return new PairDistance;
	}

	utility::vector1< Real >
	calc_rsd_sasa( Pose const & pose ) const {

		// define atom_map for main-chain and CB
		core::id::AtomID_Map< bool > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, false );
		for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
			for ( Size j = 1; j<=5; ++j ) {
				core::id::AtomID atom( j, ir );
				atom_map.set( atom, true );
			}
		}

		// calc sasa
		core::id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, option[ radius ](), false, atom_map );

		return rsd_sasa;
	} // calc_residue_sasa



public: // apply


	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using namespace ObjexxFCL::format;
		using namespace protocols::jd2;
		JobOP job_me( JobDistributor::get_instance()->current_job() );
		String me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

		Real cutoff = option[ dist ]();

		if( ssinput_ ) {
			ssinfo_->set_SSorient( pose );
			blueprint_->insert_ss_into_pose( pose );
		} else {
			core::scoring::dssp::Dssp dssp( pose );
			ssinfo_->initialize( pose, dssp.get_dssp_secstruct() );
			dssp.insert_ss_into_pose( pose );
		}

		protocols::fldsgn::topology::BB_Pos bbpos = ssinfo_->bb_pos();
		Size max_ssele = ssinfo_->ss_element_id( pose.total_residue() );
		utility::vector1< utility::vector1< Size > > ncon_sselements( max_ssele, (utility::vector1< Size >(max_ssele, 0)));
		utility::vector1< utility::vector1< bool > > ncon_calc( max_ssele, (utility::vector1< Size >(max_ssele, false)));
		for ( Size i=1 ;i<=max_ssele; ++i ) {
			for ( Size j=1 ;j<=max_ssele; ++j ) {
				ncon_sselements[i][j] = 1;
			}
		}
		utility::vector1< Size > ncon_per_res( pose.total_residue(), 1 );
		utility::vector1< bool > ncon_per_res_calc( pose.total_residue(), false );

		for( Size ii=1; ii<=pose.total_residue(); ii++ ) {

			if( ssinfo_->secstruct( ii ) != 'H' && ssinfo_->secstruct( ii ) != 'E' ) continue;
			Size ii_ssid = ssinfo_->ss_element_id( ii );
			for( Size jj=ii+1; jj<=pose.total_residue(); jj++ ) {

				Size jj_ssid = ssinfo_->ss_element_id( jj );
				if( ssinfo_->secstruct( jj ) != 'H' && ssinfo_->secstruct( jj ) != 'E' ) continue;
				if( ssinfo_->secstruct( ii ) == 'E' && ssinfo_->secstruct( jj ) == 'E' ) continue;
				if( ii_ssid == jj_ssid ) continue;

				Real dist = ( bbpos.CB( ii ) - bbpos.CB( jj ) ).length();
				ncon_calc[ ii_ssid ][ jj_ssid ] = true;

				ncon_per_res_calc[ ii ] = true;
				ncon_per_res_calc[ jj ] = true;

				if( dist <= cutoff ) {
					output_ << me << " " << ii << " " << jj << " " << dist << std::endl;
					ncon_sselements[ ii_ssid ][ jj_ssid ] ++ ;
					ncon_per_res[ ii ] ++;
					ncon_per_res[ jj ] ++;
				}

			}
		}

		//
		Size h( 0 );
		Size total( 0 );
		Size totale( 0 );
		for( Size iaa=1; iaa<=pose.total_residue(); iaa++ ) {
			if( ncon_per_res_calc[ iaa ] ) {
				total ++;
				totale += ncon_per_res[ iaa ];
				// std::cout << iaa << " " << ncon_per_res[ iaa ] << std::endl;
				if( ncon_per_res[ iaa ] > 1 ) {
					h++;
				}
			}
		}

		Real make_sure = 0.0;
		Real entropy( 0.0 );
		Real prob( 0.0 );
		for( Size iaa=1; iaa<=pose.total_residue(); iaa++ ) {
			if( ncon_per_res_calc[ iaa ] ) {
				prob = Real( ncon_per_res[ iaa ] )/Real( totale );
				entropy += -prob*( std::log( prob )/std::log( 2 ) );
				// std::cout << iaa << " " << prob << std::endl;
				make_sure += prob;
			}
		}

		// std::cout << make_sure << std::endl;

		// finalize
		Size tot( 0 );
		for ( Size i=1 ;i<=max_ssele; ++i ) {
			for ( Size j=1 ;j<=max_ssele; ++j ) {
				if( ncon_calc[i][j] ){
					// std::cout << i << " " << j << " " << ncon_sselements[i][j] << std::endl;
					tot += ncon_sselements[i][j] ;
				}
			}
		}


		Real ss_entrpy = 0.0;
		for ( Size i=1 ;i<=max_ssele; ++i ) {
			for ( Size j=1 ;j<=max_ssele; ++j ) {
				if( ncon_calc[i][j] ){
					Real prob = Real(ncon_sselements[i][j])/Real(tot);
					//std::cout << prob << " " << ncon_sselements[i][j] << " " << Real(tot) << std::endl;
					ss_entrpy += -prob*( std::log( prob )/ std::log( 2 ) );
				}
			}
		}

		core::pose::metrics::CalculatorFactory::Instance().remove_calculator( "ncontact" );
		protocols::fldsgn::NcontactsCalculator ncon;
		ncon.contact_distance( cutoff );
		ncon.ignore_loops( true );
		ncon.ignore_same_sselement( true );
		ncon.ignore_same_sheet( true );
		// ncon.use_only_calpha( true );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "ncontact", ncon.clone() );

		basic::MetricValue< Real > score, score1;
		pose.metric( "ncontact", "ss_entropy", score );
		pose.metric( "ncontact", "sidechain_heavy_apolar_atm", score1 );


		Real minv( 35.0 );
		Real maxv( 75.0 );
		Size num( 0 );
		Real total_surface( 0.0 );
		utility::vector1< Real > rsd_sasa( calc_rsd_sasa( pose ) );
		for ( Size iaa=2; iaa<=pose.total_residue()-1; iaa++ ) {
			if( ssinfo_->secstruct( iaa ) == 'E' ) {
				num++;
				Real val( 0.0 );
				if ( rsd_sasa[ iaa ] < minv ) {
					val = minv;
					// val = 0.0;
				} else if ( rsd_sasa[ iaa ] > maxv ) {
					val = maxv;
					// val = 1.0;
				} else {
					val = rsd_sasa[ iaa ];
					// val = ( rsd_sasa[ iaa ] - minv )/( maxv - minv );
				}

				total_surface += val;
				// std::cout << "UNKO" << iaa << " " << rsd_sasa[ iaa ] << std::endl;
			}
		}

		std::cout << "SSE " << me << " " << ss_entrpy << " " << entropy << " "
							<< h << "  " << total << " " << Real( h )/Real( total ) << " "
							<< total_surface/Real( num ) << std::endl;


	}

private: // data

	protocols::fldsgn::topology::SS_Info2_OP ssinfo_;
	protocols::jd2::parser::BluePrintOP blueprint_;
	std::ofstream output_;
	bool ssinput_;

};

typedef utility::pointer::owning_ptr< PairDistance > PairDistanceOP;


int
main( int argc, char * argv [] )
{
	try{
	ThisApplication::register_options();

	// init
  devel::init(argc, argv);

	// mover
	protocols::moves::MoverOP protocol;
	protocol = new PairDistance();

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;

}
