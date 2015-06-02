// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/interface_design/anchored_design/AnchorFinder.cc
/// @brief This code is intended to sift through the entire Protein Data Bank and find PDBs that fit certain criteria; in this case having loops with large numbers of cross-interface neighbors
/// @author Steven Lewis

////////////////////////////////////////////WARNING WARNING WARNING
//This code will not run well without some other changes to Rosetta.  It is intended to run across the entire PDB.
//To use it, you are strongly encouraged to "robustify" Rosetta.  This will cause bad PDBs to be ignored rather than
//causing crashes!  To robustify Rosetta:
// replace all assert statements in the vectorL (vector1) class with runtime_assert statements
// replace all assert statements in the Conformation class with runtime_assert statements
// use the -jd2:delete_old_poses flag to prevent a memory leak in the large -l environment
// use the -in::file::obey_ENDMDL flag to read in only one model from multimodel NMR PDBs
// use the -ignore_unrecognized_res flag to not crash on ligands
// See also the "RobustRosetta" documentation file.
// feel free to contact me for clarification

// Unit Headers


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>

#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/PointGraphData.hh>
#include <core/graph/UpperEdgeGraph.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>

// option key includes
#include <basic/options/keys/pose_metrics.OptionKeys.gen.hh>

// C++ headers
#include <string>

#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "apps.public.interface_design.anchored_design.AnchorFinder" );

basic::options::IntegerOptionKey const window_size("window_size");
basic::options::RealOptionKey const loopness("loopness");
basic::options::RealOptionKey const nbrs_per_residue("nbrs_per_residue");
basic::options::StringOptionKey const bestoutfile("bestoutfile");

/// @brief AnchorFinder mover
class AnchorFinderMover : public protocols::moves::Mover {
public:
	AnchorFinderMover() : out(basic::options::option[ bestoutfile ].value(), std::ios_base::app /*append*/) {
		if ( !out.good() )
			utility_exit_with_message( "Unable to open file: " + basic::options::option[ bestoutfile ].value() + "\n" );
	}

	virtual ~AnchorFinderMover(){
		out.close();
	};

	virtual
	void
	apply( core::pose::Pose & pose ) {

		core::pose::remove_nonprotein_residues( pose );
		core::pose::remove_ligand_canonical_residues( pose );

		core::Size const nres( pose.total_residue() );
		bool quit( false );
		utility::file::FileName const pdbname( protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag() );
		core::Size const num_chains( pose.conformation().num_chains() );

		if ( nres == 0 ) {
			TR << pdbname.base() << " empty, skipping (not a protein?)" << std::endl;
			quit = true;
		} else if ( num_chains < 2 ){
			TR << pdbname.base() << " less than two chains, AnchorFinderMover skipping! " << std::endl;
			quit = true;
		} else if ( nres <= 20 ) {
			TR << pdbname.base() << " too small, skipping" << std::endl;
			quit = true;
		}
		if ( quit ) {
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			return;
		}
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );

		//PointGraph will figure out whose neighbors are whose
		core::Real const distcut( basic::options::option[ basic::options::OptionKeys::pose_metrics::interface_cutoff ].value() );

		core::conformation::PointGraphOP pg( new core::conformation::PointGraph );
		core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );
		core::conformation::find_neighbors< core::conformation::PointGraphVertexData, core::conformation::PointGraphEdgeData >( pg, distcut );

		utility::vector1< utility::vector1< core::Size > > table;
		for( core::Size i( 1 ); i <= nres; ++i) table.push_back( utility::vector1< core::Size >(num_chains, 0) );
		//table is now nres x num_chains filled with zeroes

		//for each residue
		for ( core::Size res1( 1 ); res1 <= nres; ++res1 ) {
			//for each neighbor of that residue
			for ( core::conformation::PointGraph::UpperEdgeListConstIter edge_iter = pg->get_vertex( res1 ).upper_edge_list_begin(),
							edge_end_iter = pg->get_vertex( res1 ).upper_edge_list_end(); edge_iter != edge_end_iter; ++edge_iter ) {

				Size const res2 = edge_iter->upper_vertex(); //a partner
				++table[ res1 ][ pose.chain( res2 ) ];
				++table[ res2 ][ pose.chain( res1 ) ];
			}//for all partners
		}//for all residues

		//print all relevant data to particular file
		std::ostringstream datafile;
		datafile << "Rows are residues, columns are chains, data are neighbors in that chain for each residue" << std::endl;
		datafile << "residue\tchain\tPDBdata\tDSSP";
		for ( core::Size j( 1 ); j <= num_chains; ++j ) datafile << '\t' << j;
		datafile << std::endl;
		for ( core::Size i( 1 ); i <= nres; ++i){
			datafile << i << '\t' << pose.chain(i) << '\t' << pose.pdb_info()->pose2pdb(i)
							 << '\t' << pose.secstruct(i);
			for ( core::Size j( 1 ); j <= num_chains; ++j ) datafile << '\t' << table[i][j];
			datafile << std::endl;
		}
		TR << datafile.str() << std::endl;
		utility::io::ozstream rawdata( pdbname.base() + ".data" );
		rawdata << datafile.str() << std::endl;
		rawdata.close();

		core::Size const windowsize( basic::options::option[ window_size ].value() );
		core::Size const loopness_cutoff( core::Size( basic::options::option[ loopness ].value() ) * windowsize );
		core::Size const nbrs_cutoff( core::Size(basic::options::option[nbrs_per_residue].value() ) * windowsize );

		//now we look for stretches with large numbers of other neighbors AND loop-ish neighbors
		for ( core::Size window( 1 ), windowstart( 1 ), windowend( basic::options::option[ window_size ].value() );
				 windowend <= nres; ++window, ++windowstart, ++windowend ) {

			//skip cross-chain windows
			bool skip( false );
			core::Size const chain( pose.chain( windowstart ) );
			for ( core::Size i( windowstart + 1 ); i <= windowend; ++i ) if ( core::Size( pose.chain( i ) ) != chain ) skip = true;
			if ( skip ) continue;

			core::Size loopness(0);
			utility::vector1< core::Size > interface_nbrs(num_chains, 0);

			for( core::Size i(windowstart); i<=windowend; ++i){
				if( pose.secstruct(i) == 'L' ) ++loopness;
				for( core::Size j(1); j<=num_chains; ++j){
					if( j == chain ) continue; //lots of own neighbors = who cares?
					interface_nbrs[j] += table[i][j];
				}
			}

			core::Size const pdb_windowstart( pose.pdb_info()->number( windowstart ) );
			core::Size const pdb_windowend( pose.pdb_info()->number( windowend ) );

			outcopy << "PDB " << pdbname.base() << " window " << window << " loopness " << loopness << " nbrs";
			for ( core::Size j( 1 ); j <= num_chains; ++j ) outcopy << " " << interface_nbrs[ j ];
			outcopy << " start " << pose.pdb_info()->pose2pdb( windowstart )
							<< " pymol select " << pdbname.base() << " and chain " << pose.pdb_info()->chain( windowstart )
							<< " and resi " << pdb_windowstart << "-" << pdb_windowend;

			bool good_interfacenbrs( false );
			for( core::Size j( 1 ); j <= num_chains; ++j ) if( interface_nbrs[ j ] > nbrs_cutoff ) good_interfacenbrs = true;

			TR << outcopy.str() << std::endl;
			if ( ( loopness >= loopness_cutoff ) && good_interfacenbrs ) {
				out << outcopy.str() << std::endl;
			}

			outcopy.str("");
		}

		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new AnchorFinderMover );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	virtual
	std::string
	get_name() const { return "AnchorFinderMover"; }

private:

	utility::io::ozstream out;

	std::ostringstream outcopy;
};

typedef utility::pointer::shared_ptr< AnchorFinderMover > AnchorFinderMoverOP;

int main( int argc, char* argv[] )
{
	try {
	using basic::options::option;
	option.add( window_size, "window size for loops" ).def( 5 );
	option.add( loopness, "fraction of residues in window that need to be loop" ).def( 0.6 );
	option.add( nbrs_per_residue, "# cross-interface interactions per residue in loop window" ).def( 4.0 );
	option.add( bestoutfile, "file to print best loops in (ALL windows printed to tracer" ).def( "goodfile.out" );
	devel::init( argc, argv );

	protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new AnchorFinderMover ) );

	TR << "************************d**o**n**e**************************************" << std::endl;
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception" << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
