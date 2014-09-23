// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/public/scenarios/loops_from_density.cc
/// @brief A simple app that automatically generates a loop file (to be used in the loop
///        modelling protocol) from a pdb file and an electron density map.  It finds regions
///        with the poorest local match to the density, with flags controlling how far
///        loops are allowed to extend into secondary structure elements.
/// @author Frank DiMaio

#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
// AUTO-REMOVED #include <protocols/electron_density/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverStatistics.hh>


static thread_local basic::Tracer TR( "loops_from_density_main" );

namespace loops_from_density { RealOptionKey frac_loop( "loops_from_density:frac_loop" ); }
namespace loops_from_density { RealOptionKey frac_rigid( "loops_from_density:frac_rigid" ); }
namespace loops_from_density { IntegerOptionKey max_helix_melt( "loops_from_density:max_helix_melt" ); }
namespace loops_from_density { IntegerOptionKey max_strand_melt( "loops_from_density:max_strand_melt" ); }

int
main( int argc, char* argv [] )
{
	try {
	option.add( loops_from_density::max_helix_melt, "Max dist to eat into a helix" ).def( -1 );
	option.add( loops_from_density::max_strand_melt, "Max dist to eat into a strand" ).def( -1 );
	option.add( loops_from_density::frac_loop, "Fraction of loop residues" ).def( 0.35 );
	option.add( loops_from_density::frac_rigid, "Fraction of loop residues" ).def( 0.0 );

	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;

	// set parameters
	int max_helix_melt  = option[ loops_from_density::max_helix_melt ]();
	int max_strand_melt = option[ loops_from_density::max_strand_melt ]();
	core::Real frac_loop = option[ loops_from_density::frac_loop ]();
	core::Real frac_rigid = option[ loops_from_density::frac_rigid ]();

	// setup residue types
	core::chemical::ResidueTypeSetCOP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// grab -s and -l
	utility::vector1< std::string > infiles, infilelist;
	if ( option[ in::file::s ].active() )
		infiles = option[ in::file::s ]().vector(); // (-s)
	if ( option[ in::file::l ].active() )
		infilelist = option[ in::file::l ]().vector(); //  (-l)

	for(utility::vector1< std::string >::iterator i = infilelist.begin(), i_end = infilelist.end(); i != i_end; ++i) {
		utility::io::izstream input( i->c_str() );
		if ( !input.good() )
			utility_exit_with_message( "Unable to open file: " + *i + '\n' );
		std::string line;
		while( getline(input, line) )
			infiles.push_back( line );
	}

	for(utility::vector1< std::string >::iterator i = infiles.begin(), i_end = infiles.end(); i != i_end; ++i) {
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *rsd_set, *i );
		int nres = pose.total_residue();

		// stats
		utility::vector1< core::Real > perResCC( nres ), smoothPerResCC( nres );

		// get dssp parse
		core::scoring::dssp::Dssp secstruct( pose );
		ObjexxFCL::FArray1D< char > dssp_pose( nres );
		secstruct.dssp_reduced (dssp_pose);

		// align to map
		protocols::electron_density::SetupForDensityScoringMoverOP dockindens
			( new protocols::electron_density::SetupForDensityScoringMover );
		dockindens->apply( pose );

		// per-res score
		core::scoring::electron_density::getDensityMap().set_nres( nres );

		// score pose to set up graphs
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
		(*scorefxn)(pose);

		//////////////////////
		// now do the matching!
		for (int r=1; r<=nres; ++r) {
			// check for missing density ... how? look for atoms > 40A apart
			bool isMissingDens = false;
			for (int j=1; j<(int)pose.residue(r).natoms(); ++j) {
				for (int i=j+1; j<=(int)pose.residue(r).natoms(); ++j) {
					if ( (pose.residue(r).atom(i).xyz() - pose.residue(r).atom(j).xyz()).length() > 40 ) {
						isMissingDens = true;
						break;
					}
				}
			}

			if (isMissingDens) {
				perResCC[r] = 0.0;
			} else {
				perResCC[r] =
					std::max(
						core::scoring::electron_density::getDensityMap().matchRes( r , pose.residue(r), pose, NULL , false),
						0.001);
			}
		}

		/////////////////////////////////
		/////////////////////////////////
		/// phase 1 LOOPS
		// sort by filtered smoothed CC
		if (frac_loop > 0.0) {
			// outfile
			std::string loopout = i->substr( 0 , i->rfind(".pdb")) + ".loopfile";
			std::ofstream out( loopout.c_str() );

			utility::vector1< bool > loopMarker( nres, false );

			// smooth CCs
			smoothPerResCC[1] = 0.67*perResCC[1] + 0.33*perResCC[2];
			for (int r=2; r<=nres-1; ++r) {
				smoothPerResCC[r] = 0.25*perResCC[r+1] + 0.5*perResCC[r] + 0.25*perResCC[r];
			}
			smoothPerResCC[nres] = 0.67*perResCC[nres] + 0.33*perResCC[nres-1];

			// missing dens
			for (int r=1; r<=nres; ++r) {
				if (perResCC[r] == 0) smoothPerResCC[r] = 0.0;
			}

			// don't eat into sec struct elements too much
			// filter by setting their CCs artificially high
			for (int r=1; r<=nres; ++r) {
				bool in_strand = (max_strand_melt >= 0), in_helix = (max_helix_melt >= 0);
				//bool in_strand = true, in_helix = true;
				for (int i=std::max(1,r-max_strand_melt), i_end=std::min(nres,r+max_strand_melt); i<= i_end; ++i) {
					in_strand &= (dssp_pose(i) == 'E');
				}
				for (int i=std::max(1,r-max_helix_melt), i_end=std::min(nres,r+max_helix_melt); i<= i_end; ++i) {
					in_helix &= (dssp_pose(i) == 'H');
				}
				if ( in_strand || in_helix ) {
					if ( smoothPerResCC[ r ] > 0.0 )  // missing dens
						smoothPerResCC[ r ] = 2;
				}
			}

			utility::vector1< core::Real > sortPerResCC = smoothPerResCC;
			std::sort( sortPerResCC.begin(), sortPerResCC.end() );
			core::Real CCcutoff = sortPerResCC[ (int)(frac_loop*nres+0.5) ];
			for (int r=1; r<=nres; ++r) {
				loopMarker[r] = (smoothPerResCC[r] < CCcutoff);
			}

			// remove singletons
			for (int r=2; r<=nres-1; ++r) {
				if ( loopMarker[r+1] && loopMarker[r-1] ) {
					loopMarker[r] = true;
				}
				if ( !loopMarker[r+1] && !loopMarker[r-1] ) {
					loopMarker[r] = false;
				}
			}

			// fix termini
			if (loopMarker[2]) loopMarker[3] = loopMarker[1] = true;
			if (loopMarker[3]) loopMarker[2] = loopMarker[1] = true;
			if (loopMarker[4]) loopMarker[3] = loopMarker[2] = loopMarker[1] = true;
			if (loopMarker[nres-1]) loopMarker[nres-2] = loopMarker[nres] = true;
			if (loopMarker[nres-2]) loopMarker[nres-1] = loopMarker[nres] = true;
			if (loopMarker[nres-3]) loopMarker[nres-2] = loopMarker[nres-1] = loopMarker[nres] = true;
			if (loopMarker[1] && !loopMarker[2]) loopMarker[1] = false;
			if (loopMarker[nres] && !loopMarker[nres-1]) loopMarker[nres] = false;

			// debugger crap
			for (int r=1; r<=nres; ++r) {
				out << "#  " << r << " " <<  dssp_pose(r) << " " << perResCC[r] << std::endl;
			}

			//////////////////////
			// finally ... write the loopfile
			core::Size loop_start=0, loop_end=0;
			bool inloop=false;
			for (int r=1; r<=nres; ++r) {
				if ( loopMarker[r] && !inloop ) {
					loop_start = r;
					loop_end = r;
					inloop = true;
				} else if ( loopMarker[r] && inloop ) {
					loop_end = r;
				} else if ( !loopMarker[r] && inloop ) {
					out << "LOOP  " << loop_start << " " << loop_end << "  0 0" << std::endl;
					inloop = false;
				}
			}

			if ( inloop ) {
				out << "LOOP  " << loop_start << " " << loop_end << "  0 0" << std::endl;
			}
		}

		/////////////////////////////////
		/////////////////////////////////
		/// phase 2 rigid segments
		if (frac_rigid > 0.0) {
			// outfile
			std::string rigidout = i->substr( 0 , i->rfind(".pdb")) + ".rigid";
			std::ofstream r_out( rigidout.c_str() );

			utility::vector1< bool > rigidMarker( nres, false );

			// smooth CCs
			smoothPerResCC[1] = 0.67*perResCC[1] + 0.33*perResCC[2];
			for (int r=2; r<=nres-1; ++r) {
				smoothPerResCC[r] = 0.25*perResCC[r+1] + 0.5*perResCC[r] + 0.25*perResCC[r];
			}
			smoothPerResCC[nres] = 0.67*perResCC[nres] + 0.33*perResCC[nres-1];

			utility::vector1< core::Real > sortPerResCC = smoothPerResCC;
			std::sort( sortPerResCC.begin(), sortPerResCC.end() );
			core::Real CCcutoff = sortPerResCC[ (int)( (1-frac_rigid) * nres+0.5) ];
			for (int r=1; r<=nres; ++r) {
				rigidMarker[r] = (smoothPerResCC[r] > CCcutoff);
			}

			// remove singletons/doubles/triples
			for (int r=3; r<=nres-2; ++r) {
				if ( rigidMarker[r+2] && rigidMarker[r-2] ) {
					rigidMarker[r-1] = rigidMarker[r] = rigidMarker[r+1] = true;
				}
				if ( !rigidMarker[r+2] && !rigidMarker[r-2] ) {
					rigidMarker[r-1] = rigidMarker[r] = rigidMarker[r+1] = false;
				}
			}

			// fix termini
			if (rigidMarker[2]) rigidMarker[1] = true;
			if (rigidMarker[nres-1]) rigidMarker[nres] = true;

			core::Size r_start=0, r_end=0;
			bool inseg=false;
			for (int r=1; r<=nres; ++r) {
				if ( rigidMarker[r] && !inseg ) {
					r_start = r;
					r_end = r;
					inseg = true;
				} else if ( rigidMarker[r] && inseg ) {
					r_end = r;
				} else if ( !rigidMarker[r] && inseg ) {
					r_out << "RIGID  " << r_start << " " << r_end << std::endl;
					inseg = false;
				}
			}

			if ( inseg ) {
				r_out << "RIGID  " << r_start << " " << r_end << std::endl;
			}
		}
	}
	return 0;
	}
	catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
