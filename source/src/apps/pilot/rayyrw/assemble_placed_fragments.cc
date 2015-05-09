#include <devel/init.hh>

#include <core/types.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>

#include <iostream>
#include <string>
#include <map>


#include <basic/Tracer.hh>

#include <apps/pilot/rayyrw/FragMonteCarlo.hh>

//////////////////////////////////////////////////////////////////////////////

OPT_KEY( String, outfile )
OPT_KEY( Integer, nmodels )
OPT_KEY( Integer, runid )
OPT_KEY( Integer, init_type )
OPT_KEY( Boolean, sa_verbose )
// sampling
OPT_KEY( Real, sa_start_temp )
OPT_KEY( Real, sa_end_temp )
OPT_KEY( Integer, sa_nsteps )
OPT_KEY( Integer, mc_nsteps )
// scoring
OPT_KEY( Real, null_frag_score )
OPT_KEY( Real, wt_dens )
OPT_KEY( Real, wt_clash )
OPT_KEY( Real, wt_closab )
OPT_KEY( Real, wt_overlap )
// score files
OPT_KEY( String, fragidx_file)
OPT_KEY( String, densfile)
OPT_KEY( String, overlapfile)
OPT_KEY( String, nonoverlapfile)


void 
register_options(
){
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    NEW_OPT( outfile, "outfile name", "mc_sampling.out" );
    NEW_OPT( nmodels, "nmodels", 1 );
    NEW_OPT( runid, "runid", 0 );
    NEW_OPT( init_type, "init_type: 1=null; 2=random; 3=lowrmsd", 2 );
    NEW_OPT( sa_verbose, "dump frag assignments throughout simulations", false );
    NEW_OPT( sa_start_temp, "sa_start_temp", 1000.0 );
    NEW_OPT( sa_end_temp, "sa_end_temp", 1.0 );
    NEW_OPT( sa_nsteps, "sa_nsteps", 25 );
    NEW_OPT( mc_nsteps, "mc_nsteps", 1500 );
    NEW_OPT( null_frag_score, "null", -150.0 );
    NEW_OPT( wt_dens , "dens", 1.0 );
    NEW_OPT( wt_clash , "clash", 30.0 );
    NEW_OPT( wt_closab , "close", 3.0 );
    NEW_OPT( wt_overlap , "overlap", 3.0 );
    NEW_OPT( fragidx_file, "fragidx_file", "" );
    NEW_OPT( densfile, "densfile", "" );
    NEW_OPT( overlapfile, "overlapfile", "");
    NEW_OPT( nonoverlapfile, "nonoverlapfile", "" );
}


int main( int argc, char* argv[] ) 
{
  try {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    register_options();
    devel::init(argc, argv);

    FragMonteCarlo mc( option[wt_dens], 
                       option[wt_overlap],
                       option[wt_closab], 
                       option[wt_clash],
                       option[null_frag_score] );

    mc.load_scorefiles( option[fragidx_file], 
                        option[densfile],
                        option[overlapfile], 
                        option[nonoverlapfile] );


    utility::io::ozstream outstream;
    std::string outfile_name = option[outfile];

    if( ! utility::file::file_exists( outfile_name ) ){
        outstream.open( outfile_name );
    } else {
        outstream.open_append( outfile_name );
    }

    for ( int i=1; i<=option[nmodels]; ++i ){
        mc.initialize_frag_assignment( option[init_type] );

        mc.run( option[sa_verbose],
                option[sa_start_temp],
                option[sa_end_temp],
                option[sa_nsteps],
                option[mc_nsteps] );

        outstream << mc.report_results( i, option[runid] ) << std::endl;
    }

    outstream.close();

  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}

