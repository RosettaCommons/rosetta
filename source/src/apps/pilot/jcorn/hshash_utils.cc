// -*- mode:c++;tab-width:3;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file hshash_utils.cc
/// @brief utilites for working with hotspot hash files
///            Functions:	1) Make a subset of a hash based on a percentage or absolute score cutoff
///							2) Make a list of # of hotspots within 8A of a target residue. Useful for mapping hotspot density onto a target.
///									See RLC's pymol scripts data2bfactor.py and color_b.py
///
/// @author Jacob Corn
/// @created May 2008
/// @usage hshash_utils [-s <existing hash file> ] [-residue <Rosetta Residue name3>] [-scorecut <% cut or absolute score>] [-o <output filename>] [-target <target filename>]
///			-s hotspot.stubs -scorecut 0.01 --> top 1% of each residue type in hotspot.stubs
///			-s hotspot.stubs -target target.pdb --> hotspot density list

#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/conformation/Interface.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/util.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <cmath>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::hotspot_hashing;

basic::Tracer TR( "pilot_apps.jcorn.hshash_utils");

int
main( int argc, char * argv [] )
{
  try {

	OPT(hotspot::target);
	OPT(hotspot::residue);
	OPT(hotspot::hashfile);
	OPT(hotspot::density);
	OPT(hotspot::weighted_density);
	OPT(hotspot::rescore);
	OPT(hotspot::cluster);
	OPT(hotspot::colonyE);
	OPT(hotspot::sc_only);
	OPT(hotspot::rms_target);
	OPT(hotspot::rms_hotspot);
	OPT(hotspot::rms_hotspot_res);
	OPT(hotspot::fxnal_group);
	OPT(score::weights);
	OPT(score::patch);
	OPT(out::scorecut);
	OPT(in::file::s);
	OPT(out::file::o);

	devel::init(argc, argv);
	using namespace core;
	using namespace core::scoring;
	using namespace protocols::hotspot_hashing;

	// Set several variables based on the command line
	std::string tgt="";
	pose::Pose tgt_pose;
	if( option[ hotspot::target ].user() ) {
		tgt= option[ hotspot::target ]();
		core::import_pose::pose_from_pdb( tgt_pose, tgt );
	}

	utility::vector1< std::string > resnames;
	if (option[ hotspot::residue ].user() ) {
		resnames = option[ hotspot::residue ]();
		if (resnames.size() < 1 ) {
			resnames.push_back( "ALL" );
		}
	}
	else {
		resnames.push_back("ALL");
	}

	std::string hashout_fname = "";
	if ( option[ out::file::o ].user() )
	{
		hashout_fname = option[ out::file::o ]();
	}

	protocols::hotspot_hashing::HotspotStubSet stubset; // the main stubset we're going to use throughout
	utility::vector1< std::string > hashfiles( basic::options::start_files() );
	for ( utility::vector1< std::string >::const_iterator it = hashfiles.begin(), end = hashfiles.end();
		it != end; ++it ) {

		stubset.clear();
		TR << "Reading stubset " << *it << "..."<<std::endl;
		stubset.read_data( *it );

		// rescore stubs
		if( option[ hotspot::rescore ].user() ) {
			if( tgt == "" ) utility_exit_with_message( "Must specify a target to rescore against!" );

			ScoreFunctionOP scorefxn;
			if( option[ score::weights ].user() ) scorefxn = get_score_function();
			else {
				ScoreFunctionOP noenvhbond_scorefxn( get_score_function_legacy( "score13" ) );
				methods::EnergyMethodOptions options( noenvhbond_scorefxn->energy_method_options() );
				options.hbond_options().use_hb_env_dep( false );
				noenvhbond_scorefxn->set_energy_method_options( options );
				noenvhbond_scorefxn->set_weight( fa_dun, 0.1 );
				noenvhbond_scorefxn->set_weight( envsmooth, 0 );
				scorefxn = noenvhbond_scorefxn;
			}

			stubset.sc_only( option[ hotspot::sc_only]() );
			stubset = *stubset.rescore( tgt_pose, scorefxn );
		} // if option rescore

		if( option[ hotspot::cluster ].user() ) {
			stubset = *stubset.cluster();
		}

		if( option[ hotspot::colonyE ].user() ) {
			stubset = *stubset.colonyE();
		}

		// method below performs a percentage or absolute score cutoff on a stubset
		if( option[ out::scorecut].user() ) {
			Real score_cutoff = option[ out::scorecut ]();
			TR << "Making a scorecut using " << score_cutoff << "..." << std::endl;
			stubset = *stubset.subset( score_cutoff );
		}

		// calculate rms of each stub to residue in complex
		if( option[ hotspot::rms_target ].user() || option[ hotspot::rms_hotspot ].user() )  {
			using namespace core::conformation;
			if( tgt == "" ) utility_exit_with_message( "Must specify a target complex to calculate rms against!" );

			std::string rms_fname = "";
			if( option[ hotspot::rms_target ].user() ) rms_fname = option[ hotspot::rms_target ]();
			else if( option[ hotspot::rms_hotspot ].user() ) rms_fname = option[ hotspot::rms_hotspot ]();

			utility::io::ozstream ostream;
			if( rms_fname != "" ) ostream.open(rms_fname, std::ios::out);

			//core::Size const jump_num( 1 );
			//Interface iface( jump_num );
			//tgt_pose.update_residue_neighbors();
			//iface.distance( 10 );
			//iface.calculate( tgt_pose );
			bool const fxnal_group( option[hotspot::fxnal_group]() );
			if( option[ hotspot::rms_target ].user() ) { // compare each interface res to stubs of that type
				for( core::Size resnum = 1; resnum <= tgt_pose.total_residue(); ++resnum ) {
					utility::vector1< Real > residue_rms;
					//if( iface.is_interface( resnum ) ) {
						ResidueCOP tgt_res( tgt_pose.residue( resnum ) );
						std::multimap< Real, HotspotStubOP > res_stub_set( stubset.retrieve( tgt_res->name3() ) );

						for (std::multimap<Real, HotspotStubOP >::const_iterator i = res_stub_set.begin(); i != res_stub_set.end(); ++i) {
							ResidueCOP stub = i->second->residue();
							Real const rms = residue_sc_rmsd_no_super( tgt_res, stub, fxnal_group );
							residue_rms.push_back( rms );
						}
					//} // if at interface
					if( residue_rms.size() == 0 ) residue_rms.push_back( 100.0 ); // in case there were no stubs of that type or res was not at the interface
					std::sort( residue_rms.begin(), residue_rms.end() );
					ostream << tgt_pose.pdb_info()->chain( resnum ) << " " << tgt_pose.pdb_info()->number( resnum ) << " " << tgt_pose.residue( resnum ).name3() << " " << residue_rms[1] << std::endl;
				} // if rms_target
			}
			else if ( option[ hotspot::rms_hotspot ].user() ) { // compare each stub to interface residues of that type
				for( HotspotStubSet::const_iterator i = stubset.begin(); i != stubset.end(); ++i ) {
					utility::vector1< Real > stub_rms;
					ResidueCOP stub = i->second.second->residue(); // HSS iterator = pair<string, pair<real, stubOP> >
					for( core::Size resnum = 1; resnum <= tgt_pose.total_residue(); ++resnum ) {
						if( option[ hotspot::rms_hotspot_res ].user() ) {
							if( resnum != (core::Size)option[hotspot::rms_hotspot_res]() ) continue;
						}
						ResidueCOP tgt_res( tgt_pose.residue( resnum ) );
						if( tgt_res->aa() != stub->aa() ) continue; // skip if we have diff resnames
						Real const rms = residue_sc_rmsd_no_super( tgt_res, stub, fxnal_group ); // do rmsd based on functional group only
						stub_rms.push_back( rms );
					}
					std::sort( stub_rms.begin(), stub_rms.end() );
					if( stub_rms.size() > 0 ) ostream << stub->name3() << " " << i->second.second->bonus_value() << " " << stub_rms[1] << std::endl; // residue_rms[1] = best rms
				}
			} // if rms_hotspot
			ostream.close();
		} // if option rms

		// method below outputs a map of hash density relative to the specified target
		if( option[ hotspot::density].user() || option[hotspot::weighted_density].user() ) {
			if( !option[ hotspot::target ].user() ) {
				utility_exit_with_message("You must specify a target for hotspot density calculations!" );
			}
			TR << "Calculating hotspot density..." << std::endl;
			// all aa's except Gly, Cys, Pro
			utility::vector1< std::string > amino_acids;
			amino_acids.push_back( "ALA" );
			amino_acids.push_back( "ARG" );
			amino_acids.push_back( "ASN" );
			amino_acids.push_back( "ASP" );
			amino_acids.push_back( "GLU" );
			amino_acids.push_back( "GLN" );
			amino_acids.push_back( "HIS" );
			amino_acids.push_back( "ILE" );
			amino_acids.push_back( "LEU" );
			amino_acids.push_back( "LYS" );
			amino_acids.push_back( "MET" );
			amino_acids.push_back( "PHE" );
			amino_acids.push_back( "SER" );
			amino_acids.push_back( "THR" );
			amino_acids.push_back( "TRP" );
			amino_acids.push_back( "TYR" );
			amino_acids.push_back( "VAL" );

			utility::vector1< core::Size > neighbors( tgt_pose.total_residue(), 0 );
			utility::vector1< core::Real > weighted_neighbors( tgt_pose.total_residue(), 0 );

			std::string density_fname, weighted_density_fname = "";
			if( option[ hotspot::density].user() ) density_fname = option[ hotspot::density]();
			if( option[ hotspot::weighted_density].user() ) weighted_density_fname = option[ hotspot::weighted_density]();

			utility::io::ozstream ostream, weighted_ostream;
			if( density_fname != "" ) ostream.open(density_fname, std::ios::out);
			if( weighted_density_fname != "" ) weighted_ostream.open(weighted_density_fname, std::ios::out);
			for( core::Size i=1; i <= tgt_pose.total_residue(); ++i ) {
				for( utility::vector1< std::string >::const_iterator aa_it=amino_acids.begin(); aa_it != amino_acids.end(); ++aa_it ) {

					std::multimap< core::Real, protocols::hotspot_hashing::HotspotStubOP > const stubs = stubset.retrieve( *aa_it );
					core::conformation::ResidueCOP res_target( tgt_pose.residue( i ) );

					// for each hotspot of that residue type
					for( std::multimap< core::Real, protocols::hotspot_hashing::HotspotStubOP >::const_iterator stub_it = stubs.begin(); stub_it != stubs.end(); ++stub_it ) {
						core::conformation::ResidueCOP stub( stub_it->second->residue() );
						// get distances from target residue to each hotspot
						core::Real distance( res_target->xyz( res_target->nbr_atom() ).distance( stub->xyz( stub->nbr_atom() )) );
						if( distance <= 8 ) { // take care of unweighted neighbors. They're easy.
							++neighbors[ i ];
						}

						{ // Boltzmann distribution, weighted by LK-like falloff by distance
							core::Real const E = exp( -( stub_it->second->bonus_value() ) );
							core::Real const dist_falloff = 4; // where stubs should stop counting
							core::Real const sigma = 2;
							core::Real const dist_wt = exp( -( pow(distance,2)/(dist_falloff*sigma) ) ); // LK-like falloff e^(-r^2)/8
							weighted_neighbors[i] += dist_wt*E; // keep E positive for B-factor display
						}
					} // for each hotspot of a given residue type
				} // for each residue type
				weighted_neighbors[i] = log(weighted_neighbors[i]) ; // take care of ln part of Bolzmann weighting

				// write to file
				if( density_fname != "" ) ostream << tgt_pose.pdb_info()->chain( i ) << " " << tgt_pose.pdb_info()->number( i ) << " " << tgt_pose.residue( i ).name3() << " " << neighbors[ i ] << std::endl;
				if( weighted_density_fname != "" ) weighted_ostream << tgt_pose.pdb_info()->chain( i ) << " " << tgt_pose.pdb_info()->number( i ) << " " << tgt_pose.residue( i ).name3() << " " << weighted_neighbors[ i ] << std::endl;
			} // for each residue in target
			ostream.close();
			weighted_ostream.close();
		} // if option::density

		// stubset output
		if( hashout_fname != "" ) {
			if ( std::find( resnames.begin(), resnames.end(), "ALL" ) != resnames.end() ) {
				stubset.write_all( hashout_fname );
			}
			else {
				for( utility::vector1<std::string>::const_iterator res_it = resnames.begin(); res_it != resnames.end(); ++res_it ) {
					std::string resname=*res_it;
					HotspotStubSetOP res_stubs ( stubset.subset( resname, 1.0 /*scorecut, and we already cut above*/ ) );
					res_stubs->write_all( hashout_fname );
				}
			}
		} // if hashout_fname
	} // for hashname

  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
}
