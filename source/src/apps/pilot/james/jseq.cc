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
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/Energies.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/excn/Exceptions.hh>

using core::Size;
using utility::vector1;

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

class PrintSequenceMover;
typedef utility::pointer::owning_ptr< PrintSequenceMover > PrintSequenceMoverOP;

class PrintSequenceMover : public protocols::moves::Mover {
public:
	PrintSequenceMover() :
		Mover("PrintSequence")
	{}

	void apply( core::pose::Pose & pose ) {
		//if ( !getPoseExtraScore( pose, "score", my_score ) ) my_score = 0;
		core::Real my_score = pose.energies().total_energy();
		seqs_.push_back( std::make_pair( pose.sequence(), my_score ) );
	} // apply

	void print_seqs( std::ostream & out ) {
		using std::string;
		using utility::vector1;
		typedef vector1< std::pair< string, float > >::const_iterator iter;
		for ( iter it = seqs_.begin(), end = seqs_.end(); it != end; ++it ) {
			out << it->first << ' ' << it->second << std::endl;
		}
	} // print_stats

private:
	utility::vector1< std::pair< std::string, float > > seqs_;
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

	option[ out::nooutput ].value( true );

	PrintSequenceMoverOP mover ( new PrintSequenceMover() );
	not_universal_main( *mover );
	mover->print_seqs( std::cout );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
