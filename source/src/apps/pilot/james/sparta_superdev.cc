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
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/not_universal_main.hh>
#include <protocols/sparta/Sparta.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/id/AtomID_Map.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <fstream>
#include <algorithm>
#include <numeric>

#include <utility/excn/Exceptions.hh>

class SpartaSuperDeviationMover;
typedef utility::pointer::owning_ptr< SpartaSuperDeviationMover > SpartaSuperDeviationMoverOP;
using core::Size;
using core::Real;
using utility::vector1;
using ObjexxFCL::string_of;

core::Real calc_min(
	utility::vector1< core::Real > const & n
) {
	vector1< Real > numbers( n.begin(), n.end() );
	std::sort( numbers.begin(), numbers.end() );
	Real min = numbers.front();
	return min;
}

core::Real calc_median(
	utility::vector1< core::Real > const & n
) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	Real variance = n.front();
	if ( n.size() != 1 ) {

		vector1< Real > numbers( n.begin(), n.end() );
		std::sort( numbers.begin(), numbers.end() );

		Size const count( numbers.size() );
		if ( count % 2 == 0 ) {
			variance = numbers[ count / 2 ];
		} else {
			variance = numbers[ (count-1) / 2 ]
				+ numbers[ (count+1) / 2 ];
			variance *= 0.5;
		}
	}
	return variance;
}

class SpartaSuperDeviationMover : public protocols::moves::Mover {
public:
	SpartaSuperDeviationMover(
		core::pose::Pose const & native_pose,
		std::string const & cs_file
	) :
		Mover("SpartaSuperDev"),
		dump_pdbs_    ( true ),
		native_pose_  ( native_pose ),
		deviations_   ( native_pose.size(), vector1< Real >() ),
		sparta_scores_( native_pose.size(), vector1< Real >() ),
		sparta_       ( cs_file )
	{}

	virtual std::string get_name() const {
		return "SpartaSuperDev";
	}

	bool dump_pdbs() const {
		return dump_pdbs_;
	}

	void dump_pdbs( bool setting ) {
		dump_pdbs_ = setting;
	}

	void apply( core::pose::Pose & pose ) {
		using core::Real;
		using std::string;
		using core::id::AtomID;

		if ( pose.sequence() != native_pose_.sequence() ) {
			std::cout << "sequences not equal!" << std::endl;
			std::cout << "pose:   " << pose.sequence() << std::endl;
			std::cout << "native: " << native_pose_.sequence() << std::endl;
			std::exit(1);
		} else {
			runtime_assert( pose.sequence() == native_pose_.sequence() );

			// rmsd superposition
			core::scoring::calpha_superimpose_pose( pose, native_pose_ );

			string const tag( core::pose::tag_from_pose(pose) );
			pose_tags_.push_back(tag);
			static string atom_name( "CA" );

			core::id::AtomID_Map< Real > bfactors;
			core::pose::initialize_atomid_map( bfactors, pose, 0.0 );

			for ( Size ii = 1; ii <= pose.size(); ++ii ) {
				Real distance = native_pose_.residue(ii).xyz(atom_name).distance(
					pose.residue(ii).xyz( atom_name )
				);
				deviations_[ii].push_back( distance );
				for ( Size jj = 1, natoms = pose.residue_type(ii).natoms(); jj <= natoms; ++jj ) {
					bfactors[AtomID(jj,ii)] = distance;
				}
			}

			if ( dump_pdbs() ) {
				string output_fn( tag + ".super_dev.pdb" );
				utility::io::ozstream output_stream( output_fn );
				core::io::pdb::dump_bfactor_pdb( pose, bfactors, output_stream );
				output_stream.close();
			}

			// sparta scoring
			vector1< float > scores( sparta_.score_pose_per_residue(pose) );
			core::pose::initialize_atomid_map( bfactors, pose, 0.0 );

			for ( Size ii = 1; ii <= pose.size(); ++ii ) {
				float const score( scores[ii] );
				sparta_scores_[ii].push_back(score);
				for ( Size jj = 1, natoms = pose.residue_type(ii).natoms(); jj <= natoms; ++jj ) {
					bfactors[ AtomID(jj,ii) ] = score;
				}
			}

			if ( dump_pdbs() ) {
				string output_fn( tag + ".sparta_dev.pdb" );
				float const total_sparta_score(
					std::accumulate( scores.begin(), scores.end(), 0.0 ) / 4
				);
				utility::io::ozstream output_stream( output_fn );
				output_stream << "REMARK sparta_score: " << total_sparta_score << std::endl;
				core::io::pdb::dump_bfactor_pdb( pose, bfactors, output_stream );
				output_stream.close();
			}
		}
	} // apply

	void print_stats() {
		using core::Size;
		using core::Real;
		using utility::vector1;
		using core::id::AtomID;

		//Size const n_seen( deviations_.front().size() );
		vector1< Real > variance( deviations_.size(), 0 );

		core::id::AtomID_Map< Real > bfactors;
		core::pose::initialize_atomid_map( bfactors, native_pose_, 0.0 );

		std::string const tag( core::pose::tag_from_pose(native_pose_) );
		utility::io::ozstream data_output( tag + ".stats.txt" );
		//data_output << "resi deviation tag" << std::endl;
		data_output << "tag";
		for ( Size jj = 1; jj <= native_pose_.size(); ++jj ) {
			std::string const col_name( "res_" + string_of(jj) );
			data_output << " " << col_name;
		}
		data_output << std::endl;

		Size const nmodels( deviations_.front().size() );
		Size const nres( deviations_.size() );
		for ( Size jj = 1; jj <= nmodels; ++jj ) {
			data_output << pose_tags_[jj];
			for ( Size ii = 1; ii <= nres; ++ii ) {
				data_output << " " << ObjexxFCL::format::F( 8, 3, deviations_[ii][jj] );
			}
			data_output << std::endl;
		}
		data_output.close();

		for ( Size ii = 1; ii <= nres; ++ii ) {
			Real const median_dev( calc_median( deviations_[ii] ) );
			Real const min_dev( calc_min( deviations_[ii] ) );
			Size const natoms_tot( native_pose_.residue_type(ii).natoms() );

			for ( Size jj = 1, natoms = natoms_tot; jj <= natoms; ++jj ) {
				using namespace basic::options;
				using namespace basic::options::OptionKeys;
				if ( option[ run::debug ]() ) {
					bfactors[AtomID(jj,ii)] = min_dev;
				} else {
					bfactors[AtomID(jj,ii)] = median_dev;
				}
			}
		} // ii

		if ( dump_pdbs() ) {
			std::string output_fn( tag + ".super_dev.pdb" );
			utility::io::ozstream output_stream( output_fn );
			core::io::pdb::dump_bfactor_pdb( native_pose_, bfactors, output_stream );
			output_stream.close();

			vector1< float > scores( sparta_.score_pose_per_residue(native_pose_) );
			core::pose::initialize_atomid_map( bfactors, native_pose_, 0.0 );

			for ( Size ii = 1; ii <= native_pose_.size(); ++ii ) {
				float const score( scores[ii] );
				sparta_scores_[ii].push_back(score);
				for ( Size jj = 1, natoms = native_pose_.residue_type(ii).natoms(); jj <= natoms; ++jj ) {
					bfactors[ AtomID(jj,ii) ] = score;
				}
			}

			output_fn = ( tag + ".sparta_dev.pdb" );
			float const total_sparta_score(
				std::accumulate( scores.begin(), scores.end(), 0.0 ) / 4
			);
			output_stream.open( output_fn );
			output_stream << "REMARK sparta_score: " << total_sparta_score << std::endl;
			core::io::pdb::dump_bfactor_pdb( native_pose_, bfactors, output_stream );
			output_stream.close();
		}

	} // print_stats

private:
	bool dump_pdbs_;
	core::pose::Pose const & native_pose_;
	utility::vector1< std::string > pose_tags_;
	utility::vector1< utility::vector1< core::Real > > deviations_;
	utility::vector1< utility::vector1< core::Real > > sparta_scores_;
	protocols::sparta::Sparta sparta_;
}; // SpartaSuperDeviationMover

int
main( int argc, char * argv [] ) {
	try {

	using namespace protocols::moves;
	using namespace protocols::jobdist;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init( argc, argv );

	// options and usage information
	std::string usage(
		"usage: sparta_superdev "
		"-in:file:native <pdb_file> "
		"-evaluation:chemical_shifts <sparta_file> "
	);
	if ( !option[ in::file::native ].user() ) {
		utility_exit_with_message( usage );
	}
	if ( !option[ evaluation::chemical_shifts ].user() ) {
		utility_exit_with_message( usage );
	}

	// setup residue types
	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	// read in a native pose
	core::pose::Pose native_pose;
	core::import_pose::pose_from_file(
		native_pose, *rsd_set, option[ in::file::native ]()
	);
	core::pose::tag_into_pose( native_pose, option[ in::file::native ]() );
	std::string const cs_file( option[ evaluation::chemical_shifts ]()[1] );

	SpartaSuperDeviationMoverOP mover( new SpartaSuperDeviationMover( native_pose, cs_file ) );
	not_universal_main( *mover );
	mover->print_stats();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
