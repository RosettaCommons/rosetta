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


// libRosetta headers

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <protocols/loops/ccd_closure.hh>
#include <protocols/relax_protocols.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/id/AtomID.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>
#include <core/types.hh>

#include <core/scoring/sasa.hh>
#include <numeric/xyzVector.hh>

#include <basic/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>

//#include <core/chemical/residue_io.hh>

#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>

#include <core/mm/MMTorsionLibrary.hh>
#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>

#include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>


#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;
using namespace utility;

using utility::vector1;

using core::io::pdb::dump_pdb;


float cost( vector1< vector1< float >  > &quality_vs_ir,
	vector1< core::Size > &model, core::Size slen ){

	float cost = 0;
	//for ( core::Size ir = 1; ir <= slen; ir ++ ) std::cout << model[ir]; std::cout << std::endl;
	for ( core::Size ir = 1; ir <= slen; ir ++ ) {
		cost += quality_vs_ir[model[ir]][ir];
		if ( ir > 1 ) {
			if ( model[ir] != model[ir-1] ) cost+=6.0;
		}
	}
	return cost;
}


struct StructureConfidence{
	float  Confidence;
	bool   AllowLoopEnd;
};

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	try {
		devel::init(argc, argv);

		using namespace protocols::moves;
		using namespace scoring;
		using namespace core::sequence;
		using namespace basic::options;

		vector1< Sequence > seqs = read_fasta_file( option[ OptionKeys::in::file::fasta ]()[1] );
		vector1< core::pose::Pose > structures;

		vector1< Sequence >::const_iterator native = seqs.begin();

		MatrixScoringScheme blosum62;

		blosum62.read_from_file(file::FileName("BLOSUM62",utility::file::PathName("/work/mtyka//bin/")) );

		core::Size slen = native->length();


		for ( vector1< Sequence >::const_iterator it=seqs.begin(), it_end = seqs.end(); it != it_end; it++ ) {
			std::cout << it->id() << std::endl;

			pose::Pose pose;
			core::import_pose::pose_from_file(pose, it->id() + ".pdb.super.pdb" , core::import_pose::PDB_file);

			structures.push_back( pose );
			std::cout << "Adding " <<  it->id() << "  " << structures.size() << " Size" << pose.size() << std::endl;
		}


		//  it->id()
		/*
		vector1< core::Real > scores;
		for ( core::Size ir = 1; ir <= slen; ir ++ ){
		core::Real score = blosum62.score( new Sequence (*native) , new Sequence (*it), ir, ir );
		if( char((*it)[ir]) == '-' ) score = -10;   // gaps are just bad and uninformative
		scores.push_back( score );
		}

		for ( core::Size ir = 1; ir <= slen; ir ++ ){
		core::Real score =
		scores[ std::max((int)1, (int)ir - 05)    ] +
		scores[ std::max((int)1, (int)ir - 04)    ] +
		scores[ std::max((int)1, (int)ir - 03)    ] +
		scores[ std::max((int)1, (int)ir - 02)    ] +
		scores[ std::max((int)1, (int)ir - 01)    ] +
		scores[ ir                                ] +
		scores[ std::min((int)slen, (int)ir + 01) ] +
		scores[ std::min((int)slen, (int)ir + 02) ] +
		scores[ std::min((int)slen, (int)ir + 03) ] +
		scores[ std::min((int)slen, (int)ir + 04) ] +
		scores[ std::min((int)slen, (int)ir + 05) ] +
		+ 0;

		if( score < -10 ) score = -10;
		if( score > +40 ) score =  40;
		score += 50;

		std::cout << "AL "
		<< ir << "  "
		<< (*native)[ir]  << "  "
		<< (*it)[ir] << "  "
		<< score*0.01 << "  "
		<< ((char((*it)[ir]) == '-') ? '0' : '1' )
		<< std::endl;
		*/


		core::pose::Pose newpose = structures[1];

		std::ofstream out("chimera.pdb", std::ios::out | std::ios::binary);
		core::import_pose::FileData fd;
		std::string data;

		vector1< vector1< float >  > quality_vs_ir;

		// now for every sequence position find best PDB
		for ( core::Size pdb = 1; pdb <= structures.size(); pdb ++ ) {
			vector1< float > dists;
			for ( core::Size ir = 1; ir <= slen; ir ++ ) {
				//std::cout << ir << std::endl;
				float dist = numeric::distance( structures[1].xyz(  core::id::AtomID( 3, ir ) ) , structures[pdb].xyz(  core::id::AtomID( 3, ir ) ) );
				if ( dist > 8.0 ) dist = 8.0;
				//std::cout << "Q: " <<  seqs[pdb].id()  << "  " << ir << "  " << dist  << std::endl;
				dists.push_back( dist );
			}
			//smooth it

			std::ofstream confidout( std::string( seqs[pdb].id() + ".confidence").c_str() );

			vector1< float > smoothdists;
			for ( core::Size ir = 1; ir <= slen; ir ++ ) {
				float dist = 0.0;
				float div = 0.0;

				if ( ( ir - 4 ) >  0   ) { dist += 0.1 * dists[  ir - 4  ]; div += 0.1; }
				if ( ( ir - 3 ) >  0   ) { dist += 0.2 * dists[  ir - 3  ]; div += 0.2; }
				if ( ( ir - 2 ) >  0   ) { dist += 0.4 * dists[  ir - 2  ]; div += 0.4; }
				if ( ( ir - 1 ) >  0   ) { dist += 0.8 * dists[  ir - 1  ]; div += 0.8; }
				if ( ( ir - 0 ) >  0   ) { dist += 1.0 * dists[  ir - 0  ]; div += 1.0; }
				if ( ( ir + 1 ) <= slen ) { dist += 0.8 * dists[  ir + 1  ]; div += 0.8; }
				if ( ( ir + 2 ) <= slen ) { dist += 0.4 * dists[  ir + 2  ]; div += 0.4; }
				if ( ( ir + 3 ) <= slen ) { dist += 0.2 * dists[  ir + 3  ]; div += 0.2; }
				if ( ( ir + 4 ) <= slen ) { dist += 0.1 * dists[  ir + 4  ]; div += 0.1; }

				dist/=div;
				if ( dists[  ir - 0  ] > 7.99 ) dist = 10.00;
				confidout << ir << "  " << dist  << std::endl;
				smoothdists.push_back( dist );

			}

			confidout.close();

			quality_vs_ir.push_back(smoothdists);
		}

		vector1< core::Size > model;
		vector1< core::Size > newmodel;

		// now for every sequence position find best PDB
		for ( core::Size ir = 1; ir <= slen; ir ++ ) {
			core::Real bestdist = 1000000.0;
			core::Size bestpdb = 1;
			vector1< float > dists;
			for ( core::Size pdb = 2; pdb <= structures.size(); pdb ++ ) {
				core::Real distance =  (
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 9)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 8)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 7)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 6)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 5)    ]+
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 4)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 3)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 2)    ] +
					quality_vs_ir[pdb][ std::max((int)1, (int)ir - 1)    ] +
					quality_vs_ir[pdb][ ir                                ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 1) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 2) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 3) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 4) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 5) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 6) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 7) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 8) ] +
					quality_vs_ir[pdb][ std::min((int)slen, (int)ir + 9) ]

				);

				if ( distance < bestdist ) {
					bestdist = distance;
					bestpdb = pdb;
				}
			}
			model.push_back( bestpdb );
		}


		float curcost = cost( quality_vs_ir, model, slen );

		for ( int i = 0; i< 10100 ; i++ ) {
			newmodel = model;
			for ( int j = 0; j < 30; j++ ) {
				int g =rand()%(slen-2)+1;
				int a = ((rand()%2)*2) - 1;
				if ( g + a < 1 ) continue;
				if ( (int)(g + a) > (int)(slen) ) continue;
				newmodel[g] = newmodel[g+a];
			}

			float newcost = cost( quality_vs_ir, newmodel, slen );
			if ( newcost < (curcost+10) ) {
				model = newmodel;
				curcost = cost( quality_vs_ir, newmodel, slen );
			}

		}

		// now for every sequence position find best PDB
		for ( core::Size ir = 1; ir <= slen; ir ++ ) {
			std::cout << "P " << ir << "  " << seqs[ model[ir] ].id() << std::endl;
			Size a = 0;
			fd.append_residue( structures[model[ir]].residue(ir), a, structures[model[ir]] );
		}

		data = core::io::pdb::PDB_DReader::createPDBData(fd);
		out.write( data.c_str(), data.size() );

		out.close();

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

