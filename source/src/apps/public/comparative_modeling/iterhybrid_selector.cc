// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/mpi_refinement/util.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>

#include <core/types.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/io/ozstream.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <cstdio>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::format;

static basic::Tracer TR( "app.public.iterhybrid_selector" );

void
get_gremlin_d_and_width( std::string const & aa1,
	std::string const & aa2,
	core::Real &d0,
	core::Real &width )
{
	std::vector< std::string > aas_in_order = {"GLY","ALA","SER","VAL","CYS","THR","PRO","ASP","ASN","ILE","LEU","GLU","GLN","MET","HIS","LYS","PHE","TYR","ARG","TRP"};
	std::vector< std::vector< core::Real > > dats( 20 );

	dats[ 0] = {4.467,5.201,5.510,5.671,5.777,5.619,6.140,6.135,6.321,6.413,6.554,7.036,7.297,7.383,7.472,8.216, 7.966, 9.098, 9.166, 8.966};
	dats[ 1] = {0,    5.381,5.829,5.854,6.057,5.982,6.412,6.388,6.766,6.587,6.707,7.124,7.583,7.605,7.591,8.327, 8.162, 9.121, 9.365, 9.252};
	dats[ 2] = {0,    0,    6.190,6.567,6.590,6.450,6.937,6.760,7.081,7.142,7.394,7.483,7.807,8.010,8.051,8.792, 8.694, 9.594, 9.753, 9.770};
	dats[ 3] = {0.017,0,    0,    6.759,6.941,6.791,7.063,6.972,7.219,7.441,7.633,7.404,8.008,8.335,8.179,8.077, 9.057, 9.442, 9.513,10.021};
	dats[ 4] = {0.269,0.262,0,    0,    6.426,6.801,7.157,6.985,7.205,7.476,7.685,7.449,7.962,8.265,8.422,8.494, 9.026, 9.362, 9.460, 9.752};
	dats[ 5] = {0.153,0.291,0.292,0,    0,    6.676,7.062,6.971,7.159,7.442,7.642,7.628,8.055,8.397,8.221,8.715, 9.030, 9.813, 9.764, 9.980};
	dats[ 6] = {0.107,0.312,0.205,0.145,0,    0,    7.288,7.321,7.497,7.554,7.751,7.938,8.308,8.247,8.537,9.198, 8.895, 9.965,10.266, 9.719};
	dats[ 7] = {0.129,0.394,0.240,0.173,0.178,0,    0,    8.001,7.672,7.472,7.696,8.945,8.601,8.401,8.634,9.306, 9.111, 9.979,10.123, 9.867};
	dats[ 8] = {0.120,0.378,0.214,0.138,0.181,0.188,0,    0,    7.682,7.631,7.889,8.485,8.502,8.550,8.672,9.319, 9.168,10.039,10.135, 9.976};
	dats[ 9] = {0.245,0.399,0.321,0.298,0.259,0.320,0.339,0,    0,    8.096,8.342,7.949,8.302,8.874,8.523,8.329, 9.602, 9.719, 9.746,10.470};
	dats[10] = {0.193,0.289,0.323,0.287,0.299,0.307,0.416,0.392,0,    0,    8.522,8.077,8.480,9.122,8.676,8.479, 9.900, 9.889, 9.852,10.707};
	dats[11] = {0.169,0.349,0.305,0.232,0.240,0.262,0.334,0.337,0.249,0,    0,    9.863,9.328,8.870,9.454,9.842, 9.403,10.544,10.713,10.303};
	dats[12] = {0.179,0.214,0.342,0.242,0.295,0.259,0.336,0.341,0.341,0.321,0,    0,    9.074,9.102,9.391,9.667, 9.506,10.534,10.610,10.429};
	dats[13] = {0.125,0.250,0.287,0.179,0.206,0.190,0.317,0.348,0.279,0.261,0.198,0,    0,    9.530,9.396,9.096,10.253,10.400,10.250,11.110};
	dats[14] = {0.249,0.340,0.446,0.510,0.538,0.409,0.475,0.354,0.423,0.453,0.475,0.389,0,    0,   10.607,9.582, 9.602,10.843,10.879,10.661};
	dats[15] = {0.216,0.356,0.408,0.359,0.347,0.378,0.410,0.357,0.373,0.406,0.411,0.450,0.436,0,    0,   10.662, 9.344,10.627,11.322,10.136};
	dats[16] = {0.255,0.394,0.369,0.295,0.439,0.292,0.388,0.361,0.310,0.327,0.318,0.511,0.498,0.457,0,    0,    10.903,10.999,10.577,11.758};
	dats[17] = {0.206,0.380,0.435,0.383,0.203,0.417,0.457,0.325,0.289,0.379,0.401,0.443,0.401,0.342,0.333,0,     0,    11.536,11.615,11.807};
	dats[18] = {0.358,0.550,0.445,0.634,0.521,0.464,0.550,0.343,0.398,0.582,0.591,0.434,0.521,0.611,0.714,0.738, 0,     0,    12.050,11.355};
	dats[19] = {0.219,0.260,0.394,0.246,0.286,0.264,0.425,0.351,0.393,0.347,0.260,0.512,0.451,0.377,0.542,0.441, 0.460, 0,     0,    12.806};
	//dats[20] = {0.267,0.443,0.467,0.535,0.585,0.430,0.506,0.676,0.586,0.589,0.611,0.469,0.547,0.661,0.554,0.704, 0.767, 0.855, 0,     0,   };
	//dats[21] = {0.334,0.485,0.483,0.514,0.491,0.477,0.506,0.327,0.372,0.557,0.578,0.363,0.535,0.641,0.595,0.648, 0.738, 0.822, 0.704, 0,   };
	//dats[22] = {0.239,0.290,0.497,0.271,0.417,0.315,0.462,0.475,0.458,0.397,0.331,0.493,0.490,0.397,0.458,0.470, 0.447, 0.684, 0.889, 0.473};

	core::Size iaa1(99),iaa2(99);
	for ( iaa1 = 0; iaa1 < 20; ++iaa1 ) {
		if ( aa1 == aas_in_order[iaa1] ) break;
	}

	for ( iaa2 = 0; iaa2 < 20; ++iaa2 ) {
		if ( aa2 == aas_in_order[iaa2] ) break;
	}

	if ( iaa1 > 19 || iaa2 > 19 ) return;
	core::Size a = (iaa1 >= iaa2) ? iaa2 : iaa1;
	core::Size b = (iaa1 >= iaa2) ? iaa1 : iaa2;
	d0 = dats[a][b];
	width = dats[b][a]; //1.0/(dats[b][a]+0.100); // original formula
}

bool
check_silent_error( core::io::silent::SilentStructOP ss )
{
	ObjexxFCL::FArray2D< core::Real > xyz = ss->get_CA_xyz();
	bool error( false );

	// just check N-CA/CA-C distance at first 5 residues
	for ( core::Size ires = 1; ires <= 5; ++ires ) {
		core::Real dx = std::abs(xyz(1,ires) - xyz(1,ires+1));
		core::Real dy = std::abs(xyz(2,ires) - xyz(2,ires+1));
		core::Real dz = std::abs(xyz(3,ires) - xyz(3,ires+1));
		core::Real dsum = dx*dx + dy*dy + dz*dz;
		if ( dsum > 100.0 || dsum < 1e-3 || dx > 100 || dy > 100 || dz > 100 ) {
			error = true;
			break;
		}
	}
	return error;
}


core::Size
get_min_index( utility::vector1< core::Real > const &scores,
	utility::vector1< bool > const& redundant )
{
	core::Real scoremin( 1e6);
	core::Size imin( 0 );
	for ( core::Size i = 1; i <= scores.size(); ++i ) {
		if ( redundant[i] ) continue;
		if ( scores[i] < scoremin ) {
			imin = i;
			scoremin = scores[i];
		}
	}

	return imin;
}

core::Real
distance( core::io::silent::SilentStructCOP ss1,
	core::io::silent::SilentStructCOP ss2,
	//core::Size const nres,
	std::string const & mode="Score",
	core::Real const dbase2 = 4.0
)
{
	// assume aligned at the beginning
	// Get Sscore
	core::Real Sscore( 0.0 );
	core::Real rmsd( 0.0 );

	ObjexxFCL::FArray2D< core::Real > const &xyz1 = ss1->get_CA_xyz();
	ObjexxFCL::FArray2D< core::Real > const &xyz2 = ss2->get_CA_xyz();

	core::Size nres( xyz1.size()/3 );

	for ( core::Size ires = 1; ires <= nres; ++ires ) {
		core::Real dx = xyz2(1,ires) - xyz1(1,ires);
		core::Real dy = xyz2(2,ires) - xyz1(2,ires);
		core::Real dz = xyz2(3,ires) - xyz1(3,ires);
		core::Real d2 = dx*dx + dy*dy + dz*dz;
		rmsd += d2;
		Sscore += 1.0/(1.0 + d2/dbase2);
	}

	rmsd /= (core::Real)(nres);
	Sscore /= (core::Real)(nres);
	rmsd = std::sqrt( rmsd );

	core::Real dist( rmsd );
	if ( mode.compare("Sscore") == 0 ) dist = 1.0 - Sscore;

	return dist;
}

core::Size
check_redundant( core::Size iss,
	protocols::wum::SilentStructStore const &library_in,
	utility::vector1< bool > &redundant,
	core::Real const dcut )
{

	core::io::silent::SilentStructCOP ss1 = library_in.get_struct( iss-1 );

	core::Size n( 0 );
	for ( core::Size i = 1; i <= redundant.size(); ++i ) {
		if ( redundant[i] ) {
			n++;
		} else {
			if ( i == iss-1 ) {
				redundant[i] = true;
				n++;
			} else {
				core::io::silent::SilentStructCOP ss2 = library_in.get_struct( i-1 );
				core::Real const dist = distance( ss1, ss2 );
				if ( dist < dcut ) redundant[i] = true;
			}
		}
	}

	return n;
}

// assume already sorted
void
retag( protocols::wum::SilentStructStore &library_inout,
	std::string const & scorename )
{
	library_inout.sort_by( scorename ); // sort again

	// clean columns and keep necessary info only?
	std::string const prefix( option[ out::prefix ]() );
	for ( core::Size iss = 0; iss < library_inout.size(); ++iss ) {
		core::io::silent::SilentStructOP ss = library_inout.get_struct( iss );
		ss->add_string_value( "poolid", std::to_string(iss) );
		ss->set_decoy_tag( prefix+std::to_string(iss)+".pdb" );
		ss->add_energy( "nuse", 0.0 );
		// check if this is call-by-ref
	}
}

protocols::wum::SilentStructStore
read_library_simple( std::string const & silentfile )
{

	protocols::wum::SilentStructStore library;
	core::io::silent::SilentFileOptions sopt;
	core::io::silent::SilentFileData sfd( sopt );
	sfd.read_file( silentfile );

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin();
			iter != sfd.end(); ++iter ) {
		core::io::silent::SilentStructOP ss( *iter );
		library.add( ss );
	}
	return library;
}

protocols::wum::SilentStructStore
read_library_w_simpose( std::string const & silentfile, core::pose::Pose const &pose0 )
{

	protocols::wum::SilentStructStore library;
	core::io::silent::SilentFileOptions sopt;
	core::io::silent::SilentFileData sfd( sopt );
	sfd.read_file( silentfile );

	core::pose::Pose pose;
	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin();
			iter != sfd.end(); ++iter ) {
		core::io::silent::SilentStructOP ss( *iter );
		ss->fill_pose( pose );
		core::scoring::calpha_superimpose_pose( pose, pose0 ); //superimpose to input struct instead...
		ss->fill_struct( pose );
		library.add( ss );
	}
	return library;
}

protocols::wum::SilentStructStore
read_silent_input_as_library( std::string const & silentfile,
	core::pose::Pose const &pose0,
	utility::vector1< core::Real > &scores,
	core::Real const dcut_ref,
	std::string const & scorename,
	core::Size const keep_topn )
{
	protocols::wum::SilentStructStore library_in;

	bool do_score( false );
	if ( option[ score::weights ].user() ) do_score = true;
	utility::vector1< std::string > extags;
	//if( option[ macro::extags ].user() ) extags = option[ macro::extags ]();

	utility::vector1< core::io::silent::SilentStructOP > silents;

	core::scoring::ScoreFunctionOP sfxn
		= core::scoring::ScoreFunctionFactory::create_score_function( option[ score::weights ]() );

	TR << "Read input structure read from " << silentfile << ",";
	do_score ? TR << " rescore silents using input cmd, " : TR << " use scores in silentfile, ";
	TR << " objective function of " << scorename << std::endl;

	core::io::silent::SilentFileOptions sopt;
	core::io::silent::SilentFileData sfd( sopt );
	utility::vector1< bool > redundant;
	sfd.read_file( silentfile );

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(); iter != sfd.end(); ++iter ) {
		core::io::silent::SilentStructOP ss = iter->clone();
		bool error = check_silent_error( ss );

		if ( error ) continue;
		if ( extags.contains( ss->decoy_tag() ) ) continue;
		if ( do_score ) {
			core::pose::Pose pose_tmp;
			ss->fill_pose( pose_tmp );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose_tmp );

			sfxn->score( pose_tmp );
			ss->energies_from_pose( pose_tmp );
		}

		if ( dcut_ref > 0.0 && scorename.compare("penaltysum") == 0 ) {
			protocols::mpi_refinement::add_init_dev_penalty( ss, pose0, "relative",
				dcut_ref, 10.0 );
		}
		silents.push_back( ss );
		scores.push_back( ss->get_energy(scorename) );
		redundant.push_back( false );
	}

	// superimpose into lowest-energy struct and add as SilentStructStore format
	core::pose::Pose pose1, pose2;
	core::Size imin = get_min_index( scores, redundant );
	silents[imin]->fill_pose( pose1 );

	utility::vector1< core::Real > scores_cp( scores );
	std::sort( scores_cp.begin(), scores_cp.end() );

	core::Real const scorecut = scores_cp.size() > keep_topn ? scores_cp[keep_topn] : 1e10;

	TR << "Keep top " << keep_topn << " structures using score cut = " << scorecut
		<< ", superimpose to Emin structure... " << std::endl;

	// reset scores and fill with filtered ones only
	scores.resize( 0 );
	for ( core::Size i = 1; i <= silents.size(); ++i ) {
		core::Real score = silents[i]->get_energy( scorename );
		if ( score > scorecut ) continue;

		silents[i]->fill_pose( pose2 );
		//core::scoring::calpha_superimpose_pose( pose2, pose1 );
		core::scoring::calpha_superimpose_pose( pose2, pose0 ); //superimpose to input struct instead...
		silents[i]->fill_struct( pose2 ); //let's replace structure only to keep the info in silent
		scores.push_back( silents[i]->get_energy( scorename ) );

		library_in.add( silents[i] );
	}

	TR << library_in.size() << " input structures read from " << silentfile << "." << std::endl;

	return library_in;
}

utility::vector1< core::Real >
get_quota_per_silent( core::Size const nsilent,
	core::Size const nout )
{
	bool auto_quota( true );
	utility::vector1< core::Real > nquota;

	if ( option[ cm::quota_per_silent ].user() ) {
		nquota = option[ cm::quota_per_silent ]();
		if ( nquota.size() == nsilent ) auto_quota = false;
	}

	if ( auto_quota ) {
		nquota.resize( 0 );
		core::Size nsum( 0 );
		for ( core::Size isilent = 1; isilent < nsilent; ++isilent ) {
			core::Size n = core::Real(nout/nsilent);
			nquota.push_back( n );
			nsum+=n;
		}
		nquota.push_back( nout -  nsum ); // last silent
	} else {
		nquota = option[ cm::quota_per_silent ]();
	}

	TR << "Assigned quota per silent as: ";
	for ( core::Size isilent = 1; isilent <= nsilent; ++isilent ) TR << " " << nquota[isilent];
	TR << std::endl;

	return nquota;
}

protocols::wum::SilentStructStore
energy_cluster( protocols::wum::SilentStructStore const &library_in,
	core::Size const npick,
	utility::vector1< core::Real > const &scores,
	core::Real const dcut_in,
	std::string const & scorename
)
{

	utility::vector1< core::Size > picked;
	protocols::wum::SilentStructStore library_out;

	if ( library_in.size() < npick ) {
		TR << "Error: Size of input structures is smaller than npick!" << std::endl;
		return library_out;
	}

	// start clustering!
	core::Size const n( library_in.size() );

	TR << "Energy cluster " << n << " -> " << npick << std::endl;

	core::Size it( 0 );
	core::Size nredundant( 1 );
	utility::vector1< bool > redundant( n, false );

	core::Size iss = get_min_index( scores, redundant );
	core::Real minscore = scores[iss];

	core::Real dcut( dcut_in );
	while ( true ) {
		it++;
		iss = get_min_index( scores, redundant );

		if ( dcut > 0.3 && scores[iss]/minscore < 0.7 ) { // increase dcut when energy is too high
			iss = 0;
		} else if ( iss > 0 ) {
			picked.push_back( iss );
			if ( picked.size() >= npick ) break;

			nredundant = check_redundant( iss, library_in, redundant, dcut );

			TR << "picked: " << iss << " (score=" << scores[iss] << ") as " <<  picked.size() << "/" << npick
				<< ", iter=" << it << ", nredun=" << nredundant //- picked.size()
				<< std::endl;
		}

		if ( iss == 0 || nredundant == n ) {
			redundant = utility::vector1< bool >( n, false );
			for ( core::Size i = 1; i <= picked.size(); ++i ) redundant[picked[i]] = true;
			TR << "all redundant, reset and increase dcut from " << dcut << " to " << dcut*0.9 << std::endl;
			dcut *= 0.9;
		}

		// this should not happen, but ensures rare events not iterating forever
		if ( ( dcut < 0.05 ) && (picked.size() < npick ) ) {
			for ( iss = 1; iss <= scores.size(); ++iss ) {
				if ( !picked.has( iss ) ) picked.push_back( iss );
				if ( picked.size() >= npick ) break;
			}
		}
	}

	// convert into SilentStructStore
	for ( core::Size i = 1; i <= picked.size(); ++i ) {
		core::Size iss( picked[i] );
		core::io::silent::SilentStructCOP ss = library_in.get_struct(iss-1);
		library_out.add( ss->clone() ); // let's get copy of it
	}

	library_out.sort_by( scorename );

	return library_out;
}

utility::vector1< std::pair< core::Size, core::Size > >
get_contact_list( core::pose::Pose const &pose,
	utility::vector1< std::pair< core::Real, core::Real > > &d0ws )
{
	// hard coded constant
	core::Real const DTOL( 1.0 );

	utility::vector1< std::pair< core::Size, core::Size > > contacts;

	// explicitly loop over CB
	for ( core::Size res1 = 1; res1 <= pose.total_residue()-9; ++res1 ) {
		core::conformation::Residue const &rsd1 = pose.residue(res1);

		if ( rsd1.is_virtual_residue() ) continue;
		core::Vector crd1 = rsd1.has(" CB ") ? rsd1.xyz(" CB ") : rsd1.xyz(" CA ");

		for ( core::Size res2 = res1+9; res2 <= pose.total_residue(); ++res2 ) {
			core::conformation::Residue const &rsd2 = pose.residue(res2);

			if ( rsd2.is_virtual_residue() ) continue;
			core::Vector crd2 = rsd2.has(" CB ") ? rsd2.xyz(" CB ") : rsd2.xyz(" CA ");

			core::Real d = crd1.distance( crd2 );
			core::Real d0( 8.0 ), width( 0.5 );
			get_gremlin_d_and_width( rsd1.name3(),rsd2.name3(), d0, width );
			if ( d < d0 + DTOL + 2.0*width ) { // some magical formula considering sigmoid func.. perhaps not important
				//std::cout << res1 << " " << res2 << " " << rsd1.name3() << " " << rsd2.name3() << " " << d
				//<< " " << d0 << std::endl;

				contacts.push_back( std::make_pair( res1, res2 ) );
				d0ws.push_back( std::make_pair( d0, width ) );
			}
		}
	}
	return contacts;
}

void
get_contact_devs( core::pose::Pose const &pose,
	utility::vector1< std::pair< core::Size, core::Size > > const &contacts,
	utility::vector1< std::pair< core::Real, core::Real > > const &d0ws,
	utility::vector1< utility::vector1< core::Real > > &contact_devs
)
{
	for ( core::Size ipair = 1; ipair <= contacts.size(); ++ipair ) {
		core::conformation::Residue const &rsd1 = pose.residue( contacts[ipair].first );
		core::conformation::Residue const &rsd2 = pose.residue( contacts[ipair].second );

		core::Vector crd1 = rsd1.has(" CB ") ? rsd1.xyz(" CB ") : rsd1.xyz(" CA ");
		core::Vector crd2 = rsd2.has(" CB ") ? rsd2.xyz(" CB ") : rsd2.xyz(" CA ");
		core::Real d = crd1.distance( crd2 );
		contact_devs[ipair].push_back( d-d0ws[ipair].first );
	}
}


void report_and_dump( protocols::wum::SilentStructStore const & library_inout,
	core::pose::Pose const &pose0,
	bool const dump_cst,
	bool const is_iterative_selection )
{
	core::io::silent::SilentFileOptions sopt;
	core::io::silent::SilentFileData sfd( sopt );
	std::string filename( option[ out::file::silent ]() );
	std::string prefix( "" );
	if ( option[ out::prefix ].user() ) prefix = option[ out::prefix ]();

	utility::vector1< std::pair< core::Size, core::Size > > contacts;
	utility::vector1< std::pair< core::Real, core::Real > > d0ws;
	utility::vector1< utility::vector1< core::Real > > contact_devs;
	if ( dump_cst ) {
		contacts = get_contact_list( pose0, d0ws );
		TR << "contacts: " << contacts.size() << std::endl;
		contact_devs.resize( contacts.size() );
	}

	core::Size i( 0 );
	for ( auto it = library_inout.begin();
			it != library_inout.end(); ++it, ++i ) {
		core::io::silent::SilentStructOP ss( *it );
		if ( prefix != "" ) {
			std::stringstream pdbname;
			pdbname << prefix << "." << i << ".pdb";
			ss->set_decoy_tag( pdbname.str() );
			if ( !is_iterative_selection ) ss->add_energy( "poolid", i ); // add poolid
		}

		if ( dump_cst ) {
			core::pose::Pose pose_tmp; ss->fill_pose( pose_tmp );
			get_contact_devs( pose_tmp, contacts, d0ws, contact_devs );
		}
		sfd.write_silent_struct( *ss, filename, false );
	}

	if ( dump_cst ) {
		TR << "dump cst" << std::endl;
		core::Real const k( 1.0 ), tol( 2.0 ); // hard-coded

		std::string const cstfile1 = option[ constraints::dump_cst_set ]();
		utility::io::ozstream outcst1( cstfile1.c_str() );

		for ( core::Size ipair = 1; ipair <= contacts.size(); ++ipair ) {
			core::Size res1 = contacts[ipair].first, res2 = contacts[ipair].second;
			std::string atm1 = pose0.residue(res1).has("CB") ? "CB":"CA";
			std::string atm2 = pose0.residue(res2).has("CB") ? "CB":"CA";

			core::Vector v1 = pose0.residue(res1).xyz(atm1);
			core::Vector v2 = pose0.residue(res2).xyz(atm2);
			//core::Real d0 = v1.distance(v2);

			core::Real w( 0.0 );
			utility::vector1< core::Real > const &devs = contact_devs[ipair];
			core::Size nstruct( devs.size() );
			for ( core::Size istruct = 1; istruct <= nstruct; ++istruct ) {
				w += 1.0/(1.0 + k*devs[i]*devs[i]);
			}
			w /= (core::Real)(nstruct);

			outcst1 << "AtomPair " << A(3,atm1) << " " << I(4,int(res1)) << " " << A(3,atm2) << " " << I(4,int(res2))
				<< " SCALARWEIGHTEDFUNC " <<  F(8,5,w)
				//<< " BOUNDED 0 " << F(8,3,d0+tol) << " 1.0 0.5" << std::endl;
				<< " BOUNDED 0 " << F(8,3,d0ws[ipair].first+tol) << " 1.0 0.5" << std::endl;
		}
	}
}

int main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		utility::vector1< core::Size > seeds;
		if ( option[ cm::seeds ].user() ) seeds = option[ cm::seeds ]();
		bool const is_iterative_selection( seeds.size() > 0 );

		core::Real dcut( option[ cm::similarity_cut ]() );

		core::Real dcut_ref( 0.0 );
		bool add_penalty_to_ref( false );

		if ( option[ cm::refsimilarity_cut ].user() ) {
			dcut_ref = option[ cm::refsimilarity_cut ]();
			if ( dcut_ref > 0.0 ) add_penalty_to_ref = true;
		}
		std::string scorename = add_penalty_to_ref ? "penaltysum" :  "score";
		std::string prefix = option[ out::prefix ].user()? option[ out::prefix ]() : "";
		bool dump_cst( option[ constraints::dump_cst_set ].user() );

		core::chemical::ResidueTypeSetCOP rsd_set
			= core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		core::pose::Pose pose0;
		core::import_pose::pose_from_file( pose0, *rsd_set,
			option[ in::file::template_pdb ](1), core::import_pose::PDB_file );

		protocols::wum::SilentStructStore library_inout, library_add;

		// read reference structures from in:file:template_silent
		if ( option[ in::file::template_silent ].user() ) {
			//library_inout = read_library_simple( option[ in::file::template_silent ]() );
			library_inout = read_library_w_simpose( option[ in::file::template_silent ](), pose0 );
		}

		core::Size const nout( option[ out::nstruct ]() );
		core::Size const nremain_reset =
			option[ cm::nremain_reset ].user()? option[ cm::nremain_reset ]() : 3; // default

		if ( is_iterative_selection ) {
			// takes only 1 input
			library_add = read_library_simple( option[ in::file::silent ](1) );

			// read template information if available
			protocols::mpi_refinement::MultiObjective fobj;
			fobj.set_init_pose( pose0 );

			// set similarity cut
			if ( pose0.total_residue() > 0 ) fobj.set_iha_cut( dcut_ref );

			//fobj.update_library_NSGAII( library_inout, library_add, nout, false );
			fobj.set_nremain_reset( nremain_reset );
			fobj.update_library_seeds( library_inout, library_add, dcut, seeds, prefix,
				scorename );

		} else {
			// read newly generated structures from in:file:silent
			core::Size nsilent = option[ in::file::silent ]().size();

			utility::vector1< core::Real > nquota = get_quota_per_silent( nsilent, nout );

			core::Size keep_topn = nout*2>100 ? 100 : nout*2; // take twice bigger of noutput or 100 whichever bigger
			for ( core::Size isilent = 1; isilent <= nsilent; ++isilent ) {
				utility::vector1< core::Real > scores;
				library_add = read_silent_input_as_library( option[ in::file::silent ](isilent), pose0, scores,
					dcut_ref, scorename,
					keep_topn ); // keep top100

				protocols::wum::SilentStructStore library_tmp =
					energy_cluster( library_add, core::Size(nquota[isilent]), scores,
					dcut, scorename );

				library_inout.add( library_tmp );
			}
			retag( library_inout, scorename );
		}

		// dump out into silent
		report_and_dump( library_inout, pose0, dump_cst, is_iterative_selection );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
