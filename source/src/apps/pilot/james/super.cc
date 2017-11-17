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
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

using core::Size;
using core::Real;
using utility::vector1;

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


class SuperDeviationMover;
typedef utility::pointer::owning_ptr< SuperDeviationMover > SuperDeviationMoverOP;

class SuperDeviationMover : public protocols::moves::Mover {
public:
	SuperDeviationMover(
		core::pose::Pose const & native_pose
	) :
		Mover("Super"),
		native_pose_( native_pose ),
		deviations_( native_pose.size(), vector1< Real >() )
	{}

	virtual std::string get_name() const {
		return "SuperDev";
	}

	void apply( core::pose::Pose & pose ) {
		using core::Real;
		using std::string;
		assert( pose.sequence() == native_pose_.sequence() );

		core::scoring::calpha_superimpose_pose( pose, native_pose_ );

		static string atom_name( "CA" );
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			Real distance = native_pose_.residue(ii).xyz( atom_name ).distance(
				pose.residue(ii).xyz( atom_name )
			);
			deviations_[ii].push_back( distance );
		}
	} // apply

	void print_stats() {
		// for each coordinate
		// measure the median of distances from that center-of-mass

		using core::Size;
		using core::Real;
		using utility::vector1;

		Size const n_seen( deviations_.front().size() );
		vector1< Real > variance( deviations_.size(), 0 );
		for ( Size ii = 1; ii <= deviations_.size(); ++ii ) {
			std::sort( deviations_[ii].begin(), deviations_[ii].end() );

			if ( n_seen % 2 == 0 ) {
				variance[ii] = deviations_[ii][ n_seen / 2 ];
			} else {
				variance[ii]  = deviations_[ii][ (n_seen-1) / 2 ]
					+ deviations_[ii][ (n_seen+1)/2 ];
				variance[ii] *= 0.5;
			}
		} // ii

		std::cout << "resi variance" << std::endl;
		for ( Size ii = 1; ii <= variance.size(); ++ii ) {
			std::cout << ii << ' ' << variance[ii] << std::endl;
		}
	} // print_stats

private:
	core::pose::Pose const & native_pose_;
	utility::vector1< utility::vector1< core::Real > > deviations_;
};

int
main( int argc, char * argv [] ) {
	try {

	using namespace core::chemical;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace protocols::moves;
	using namespace protocols::jobdist;

	devel::init( argc, argv );

	// setup residue types
	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	// read in a native pose
	core::pose::Pose native_pose;
	core::import_pose::pose_from_file(
		native_pose, *rsd_set, option[ in::file::native ]()
	);

	SuperDeviationMoverOP mover ( new SuperDeviationMover( native_pose ) );
	not_universal_main( *mover );
	mover->print_stats();

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
