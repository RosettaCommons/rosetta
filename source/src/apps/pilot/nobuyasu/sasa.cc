// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga

// Unit headers

// Project headers
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/dssp/Dssp.hh>

// type headers
#include <core/types.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/string_util.hh>

// C++ header
#include <fstream>
#include <map>
#include <utility/excn/Exceptions.hh>

typedef core::Size Size;
typedef std::string String;

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "sasa" );

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, output )
OPT_KEY( File, rasmol )
OPT_KEY( Real, pore_radius )
OPT_KEY( Real, asa_core )
OPT_KEY( Real, asa_surface )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( output, "output file ", "out.sasa" );
	NEW_OPT( rasmol, "output filename of rasmol script", "out.sasa.ras" );
	NEW_OPT( pore_radius, "pore size", 2.0 );
	NEW_OPT( asa_core, "asa for defining core positions", 20.0 );
	NEW_OPT( asa_surface, "asa for defining core positions", 40.0 );
}


namespace sasa
{
	StringOptionKey output_name( "sasa:output_name" );
}

/// local mover for testing purposes
class Sasa : public protocols::moves::Mover {
public:
	Sasa():
		pore_radius_( option[ pore_radius ]() ),
		core_( option[ asa_core ]() ),
		surface_( option[ asa_surface ]() )
	{
		std::ostringstream filename;
		filename << option[ output ]();
		output_.open( filename.str().c_str() ,std::ios::out );

		filename.str("");
		filename << option[ rasmol ]();
		rasmol_.open( filename.str().c_str() ,std::ios::out );
	}

	virtual ~Sasa(){};

	virtual std::string get_name() const { return "Sasa"; }

	virtual
	void
	apply( core::pose::Pose & pose ){

		using namespace protocols::jd2;
		JobOP job_me( JobDistributor::get_instance()->current_job() );
		String me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

		// define atom_map for main-chain and CB
		core::id::AtomID_Map< bool > atom_map;
		core::pose::initialize_atomid_map( atom_map, pose, false );
		for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
			for ( Size j = 1; j<=5; ++j ) {
				id::AtomID atom( j, ir );
				atom_map.set( atom, true );
			}
		}

		core::scoring::dssp::Dssp dssp( pose );
		String ss = dssp.get_dssp_secstruct();

		// calc sasa
		core::id::AtomID_Map< Real > atom_sasa;
		utility::vector1< Real > rsd_sasa;
		scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius_, false, atom_map );

		String aa( pose.sequence() );

		output_ << me << ' ' << 1 << ' ' << aa.substr( 0, 1 ) << ' ' << ss.substr( 0, 1 ) << ' ' << "999" << std::endl;
		for( Size i=2;i<=pose.total_residue()-1; ++i ){
			output_ << me << ' ' << i << ' ' << aa.substr( i-1, 1 ) << ' ' << ss.substr( i-1, 1 ) << ' ' << rsd_sasa[ i ] << std::endl;
		}
		output_ << me << ' ' << pose.total_residue() << ' ' << aa.substr( pose.total_residue()-1, 1 ) << ' ' << "999" << std::endl;


		for( Size i=1; i<=pose.total_residue(); i++ ) {
			rasmol_ << "select " << i << std::endl;
			if ( core_ > rsd_sasa[ i ] ) {
				rasmol_ << "color blue " << std::endl;
			} else if ( surface_ < rsd_sasa[ i ] ) {
				rasmol_ << "color red " << std::endl;
			} else {
				rasmol_ << "color yellow " << std::endl;
			}
		}

		/*
		job_me->add_string("somestuff");
		job_me->add_string_string_pair("alhpa", "btea");
		job_me->add_string_real_pair("theanswer", 42.0);
		JobDistributor::get_instance()->job_outputter()->other_pose( job_me, pose, "lookit");
		JobDistributor::get_instance()->job_outputter()->file( job_me, "hey here's some junk from job " + me + "\n" );
		*/


		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new Sasa;
	}


private:
	std::ofstream output_;
	std::ofstream rasmol_;
	Real pore_radius_;
	Real core_;
	Real surface_;

};

typedef utility::pointer::owning_ptr< Sasa > SasaOP;


int
main( int argc, char * argv [] )
{
	try{
	//
	// option.add( sasa::output_name, "output name" );
	// init
  //devel::init(argc, argv);

	ThisApplication::register_options();
	devel::init( argc, argv );

	// mover
	protocols::moves::MoverOP protocol;
	protocol = new Sasa();

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
