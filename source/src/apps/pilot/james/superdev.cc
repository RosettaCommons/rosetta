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
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>

#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/id/AtomID_Map.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <fstream>
#include <algorithm>

//Auto Headers
#include <core/pose/util.tmpl.hh>

#include <utility/excn/Exceptions.hh>

class SuperDeviationMover;
typedef utility::pointer::owning_ptr< SuperDeviationMover > SuperDeviationMoverOP;
using core::Size;
using core::Real;
using utility::vector1;
using ObjexxFCL::string_of;

core::Real calc_min(
	utility::vector1< core::Real > const & n
) {
	return *min_element( n.begin(), n.end() );
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

class SuperDeviationMover : public protocols::moves::Mover {
public:
	SuperDeviationMover(
		core::pose::Pose const & native_pose,
		bool superimpose = true
	) :
		Mover("SuperDev"),
		superimpose_(superimpose),
		native_pose_( native_pose ),
		deviations_( native_pose.total_residue(), vector1< Real >() )
	{}


	virtual std::string get_name() const {
		return "SuperDev";
	}

	void apply( core::pose::Pose & pose ) {
		using core::Real;
		using std::string;

		if ( pose.sequence() != native_pose_.sequence() ) {
			std::cout << "sequences not equal!" << std::endl;
			std::cout << "pose:   " << pose.sequence() << std::endl;
			std::cout << "native: " << native_pose_.sequence() << std::endl;
		} else {
			runtime_assert( pose.sequence() == native_pose_.sequence() );

			if ( superimpose_ ) {
				core::scoring::calpha_superimpose_pose( pose, native_pose_ );
			}

			string const tag( core::pose::tag_from_pose(pose) );
			pose_tags_.push_back(tag);
			static string atom_name( "CA" );

			core::id::AtomID_Map< Real > bfactors;
			core::pose::initialize_atomid_map( bfactors, pose, 0.0 );

			for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
				Real distance = native_pose_.residue(ii).xyz(atom_name).distance(
					pose.residue(ii).xyz( atom_name )
				);
				deviations_[ii].push_back( distance );
				for ( Size jj = 1, natoms = pose.residue_type(ii).natoms(); jj <= natoms; ++jj ) {
					core::id::AtomID atom_id(jj,ii);
					bfactors[atom_id] = distance;
				}
			}

			string output_fn( tag + ".superdev.pdb" );
			utility::io::ozstream output_stream( output_fn );
			core::io::pdb::dump_bfactor_pdb( pose, bfactors, output_stream );
			output_stream.close();
		}
	} // apply

	void print_stats() {
		using core::Size;
		using core::Real;
		using utility::vector1;

		//Size const n_seen( deviations_.front().size() );
		vector1< Real > variance( deviations_.size(), 0 );

		core::id::AtomID_Map< Real > bfactors;
		core::pose::initialize_atomid_map( bfactors, native_pose_, 0.0 );

		std::string const tag( core::pose::tag_from_pose(native_pose_) );
		utility::io::ozstream data_output( tag + ".stats.txt" );
		//data_output << "resi deviation tag" << std::endl;
		data_output << "tag";
		for ( Size jj = 1; jj <= native_pose_.total_residue(); ++jj ) {
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
			for ( Size jj = 1, natoms = native_pose_.residue_type(ii).natoms(); jj <= natoms; ++jj ) {
				core::id::AtomID atom_id(jj,ii);
				using namespace basic::options;
				using namespace basic::options::OptionKeys;
				if ( option[ run::debug ]() ) {
					bfactors[atom_id] = min_dev;
				} else {
					bfactors[atom_id] = median_dev;
				}
			}
		} // ii

		std::string output_fn( tag + ".superdev.pdb" );
		utility::io::ozstream output_stream( output_fn );
		core::io::pdb::dump_bfactor_pdb( native_pose_, bfactors, output_stream );
		output_stream.close();
	} // print_stats

private:
	bool superimpose_;
	core::pose::Pose const & native_pose_;
	utility::vector1< std::string > pose_tags_;
	utility::vector1< utility::vector1< core::Real > > deviations_;
};

int
main( int argc, char * argv [] ) {
	try {

	using namespace protocols::moves;
	using namespace protocols::jobdist;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init( argc, argv );

	// setup residue types
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	// read in a native pose
	core::pose::Pose native_pose;
	core::import_pose::pose_from_pdb(
		native_pose, *rsd_set, option[ in::file::native ]()
	);
	core::pose::tag_into_pose( native_pose, option[ in::file::native ]() );

	SuperDeviationMoverOP mover( new SuperDeviationMover( native_pose ) );
	not_universal_main( *mover );
	mover->print_stats();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // int main
