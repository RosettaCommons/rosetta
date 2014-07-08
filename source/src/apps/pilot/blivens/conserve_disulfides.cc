// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Confirms that disulfides are conserved across the cen<=>fa transition
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created September 2008
/// @details
/// @section cli Command Line
/// @code conserve_disulfides -s input.pdb -o output.pdb -database db @endcode


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>

//Core Chemistry
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/pose/disulfide_util.hh>


//Command line Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Packing
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>

//Scoring
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <apps/pilot/blivens/disulfides.hh>

#include <utility>

using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::chemical;
using namespace core::conformation;

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

basic::Tracer TR( "pilot_apps.blivens.disulfide_handoff" );

int
usage(char* msg)
{
	TR	<< "usage: disulfide_staple -s input.pdb -o output.pdb -database db" << endl
		<< msg << endl;
	exit(1);
}


int main( int argc, char * argv [] )
{
  try {
	//init options system
	option.add_relevant( in::file::s );
	option.add_relevant( out::prefix );

	devel::init(argc, argv);

	if( !option[ in::file::s ].user() ) {
		return usage("No in file given: Use -s to designate a pdb file");
	}
	string infile(basic::options::start_file() );

	if( ! option[ out::prefix ].user() ) {
		return usage("No out file given: Use -prefix to designate an output file");
	}
	string outprefix(option[ out::prefix ]() );
	//done with options

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, infile );

	scoring::ScoreFunctionOP fa_sfxn = scoring::get_score_function_legacy( scoring::PRE_TALARIS_2013_STANDARD_WTS );
	scoring::ScoreFunctionOP cen_sfxn =scoring::ScoreFunctionFactory::create_score_function(scoring::CENTROID_WTS);

	//initialize vectors of all disulf bonds
	vector1< pair<Size,Size> > initial_disulf;
	core::conformation::disulfide_bonds(pose, initial_disulf);
	for(vector1< pair<Size,Size> >::const_iterator ds_it = initial_disulf.begin(); ds_it != initial_disulf.end(); ++ds_it) {
		TR <<"Found Initial disulf at "<<ds_it->first<<" to "<<ds_it->second<<endl;
	}

	//Convert to centroid
	core::util::switch_to_residue_type_set( pose, chemical::CENTROID);

	vector1< pair<Size,Size> > cen_disulf;
	core::conformation::disulfide_bonds(pose, cen_disulf);
	for(vector1< pair<Size,Size> >::const_iterator ds_it = cen_disulf.begin(); ds_it != cen_disulf.end(); ++ds_it) {
		TR <<"Found Centroid disulf at "<<ds_it->first<<" to "<<ds_it->second<<endl;
	}

	string outfile = outprefix;
	outfile.append("centroid.pdb");
	pose.dump_scored_pdb(outfile, *cen_sfxn, "");

	//Convert to fa
	core::util::switch_to_residue_type_set( pose, chemical::FA_STANDARD);

	vector1< pair<Size,Size> > fa_disulf;
	core::conformation::disulfide_bonds(pose, fa_disulf);
	for(vector1< pair<Size,Size> >::const_iterator ds_it = fa_disulf.begin(); ds_it != fa_disulf.end(); ++ds_it) {
		TR <<"Found FA disulf at "<<ds_it->first<<" to "<<ds_it->second<<endl;
	}

	outfile = outprefix;
	outfile.append("fa.pdb");
	pose.dump_scored_pdb(outfile, *fa_sfxn, "");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // end main

