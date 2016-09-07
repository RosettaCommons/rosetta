// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/relax/FastRelax.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>


//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.functions.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("BestBetaBackbones");

// application specific options

inline
bool atoms_within_angle(
	Pose const & pose,
	Size resi,
	Size ai,
	Size aj
) {
	AtomIndices neighbors = pose.residue_type( resi ).bonded_neighbor( ai );
	for ( Size ii = 1; ii <= neighbors.size(); ++ii ) {
		if ( pose.residue_type( resi ).atoms_are_bonded( neighbors[ii], aj ) ) {
			return true;
		}
	}
	return false;
}

bool
bump_check( Pose const & pose, Real const prop_vdw ) {
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size ai = 1; ai <= pose.residue( ii ).natoms(); ++ai ) {

			Real ai_r = pose.residue_type( ii ).atom_type( ai ).lj_radius();

			for ( Size jj = ii; jj <= pose.total_residue(); ++jj ) {
				for ( Size aj = 1; aj <= pose.residue( jj ).natoms(); ++aj ) {
					if ( jj == ii && ai == aj ) continue;
					if ( jj == ii && pose.residue_type( ii ).atoms_are_bonded( ai, aj ) ) continue;
					if ( jj == ii && atoms_within_angle( pose, ii, ai, aj ) ) continue;

					Real aj_r = pose.residue_type( jj ).atom_type( aj ).lj_radius();

					// skip bonded -- explicit for now!
					// Just make it so res 2's C can't clash with residue 3
					if ( ii == 3 && pose.residue( jj ).atom_name( aj ) == " C  " ) continue;
					if ( ii == 1 && pose.residue( jj ).atom_name( aj ) == " N  " ) continue;
					if ( jj == 3 && pose.residue( ii ).atom_name( ai ) == " C  " ) continue;
					if ( jj == 1 && pose.residue( ii ).atom_name( ai ) == " N  " ) continue;

					// Also handle 2 bond inter cases: 3-H 2-C, 3-N 2-O
					if ( ii == 1 && pose.residue_type( ii ).atom_name( ai ) == " C  " && jj == 2 ) {
						AtomIndices neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " N  " ) );
						bool is_within_angle = false;
						for ( Size kk = 1; kk <= neigh.size(); ++kk ) {
							if ( aj == neigh[ kk ] ) is_within_angle = true;
						}
						if ( is_within_angle ) continue;
					}
					if ( ii == 2 && pose.residue_type( ii ).atom_name( ai ) == " N  " && jj == 1 ) {
						AtomIndices neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " C  " ) );
						bool is_within_angle = false;
						for ( Size kk = 1; kk <= neigh.size(); ++kk ) {
							if ( aj == neigh[ kk ] ) is_within_angle = true;
						}
						if ( is_within_angle ) continue;
					}
					if ( ii == 2 && pose.residue_type( ii ).atom_name( ai ) == " C  " && jj == 3 ) {
						AtomIndices neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " N  " ) );
						bool is_within_angle = false;
						for ( Size kk = 1; kk <= neigh.size(); ++kk ) {
							if ( aj == neigh[ kk ] ) is_within_angle = true;
						}
						if ( is_within_angle ) continue;
					}
					if ( ii == 3 && pose.residue_type( ii ).atom_name( ai ) == " N  " && jj == 2 ) {
						AtomIndices neigh = pose.residue_type( jj ).bonded_neighbor( pose.residue_type( jj ).atom_index( " C  " ) );
						bool is_within_angle = false;
						for ( Size kk = 1; kk <= neigh.size(); ++kk ) {
							if ( aj == neigh[ kk ] ) is_within_angle = true;
						}
						if ( is_within_angle ) continue;
					}

					Real const dsq = pose.residue( ii ).xyz( ai ).distance_squared(
						pose.residue( jj ).xyz( aj ) );
					if ( dsq < pow( ai_r + aj_r, 2 ) * prop_vdw * prop_vdw ) {
						//TR << "Bump found: " << ii << " - " << pose.residue_type( ii ).atom_name( ai ) << " and "
						// << jj << " - " << pose.residue_type( jj ).atom_name( aj ) << std::endl;
						//TR << "distance sq " << dsq << " radii " << ai_r << " " << aj_r << std::endl;
						return true; // bump found
					}
				}
			}
		}
	}
	return false;
}

void
report_proportion_accessible(
	std::string const & rtname,
	Real const bump_frac
) {
	TR << "Examining " << rtname << "..." << std::endl;
	scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );

	Pose pose;
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	ResidueTypeCOP ace = ResidueTypeFinder( *restype_set).name3( "ACE" ).get_representative_type();
	ResidueTypeCOP rt = ResidueTypeFinder( *restype_set).name3( rtname ).get_representative_type(); ResidueTypeCOP nme = ResidueTypeFinder( *restype_set).name3( "NME" ).get_representative_type();
	ResidueOP aceres( new Residue( ace, true ) );
	ResidueOP res( new Residue( rt, true ) );
	ResidueOP nmeres( new Residue( nme, true ) );

	pose.append_residue_by_jump( *aceres, true );
	pose.append_residue_by_bond( *res, true );
	pose.append_residue_by_bond( *nmeres, true );

	// There are going to be 72^3 measurements.
	Size const denominator = 72 * 72 * 72;
	Size numerator = 0;

	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_chi( 2, true );
	for ( Size ii = 1; ii <= pose.residue_type( 2 ).nchi(); ++ii ) {
		pose.set_torsion( id::TorsionID( 2, id::CHI, ii ), 180 );
	}
	protocols::simple_moves::MinMoverOP min( new protocols::simple_moves::MinMover( mm, score_fxn, "linmin_iterated", 0.001, true ) );

	for ( Real phi = -175; phi <= 180; phi += 5 ) {
		for ( Real tht = -175; tht <= 180; tht += 5 ) {
			for ( Real psi = -175; psi <= 180; psi += 5 ) {
				pose.set_torsion( id::TorsionID( 2, id::BB, 1 ), phi );
				pose.set_torsion( id::TorsionID( 2, id::BB, 2 ), tht );
				pose.set_torsion( id::TorsionID( 2, id::BB, 3 ), psi );
				min->apply( pose );
				if ( ! bump_check( pose, bump_frac ) ) {
					TR << phi << " " << tht << " " << psi << "\n";
					++numerator;
				}
			}
		}
	}

	TR << "Proportion for " << rtname << " is " << numerator << "/" << denominator << " or " << (Real)(numerator)/(Real)(denominator) << std::endl;
}

void
report_proportion_accessible_alpha(
	std::string const & rtname,
	Real const bump_frac
) {
	TR << "Examining " << rtname << "..." << std::endl;

	scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );

	Pose pose;
	chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	ResidueTypeCOP ace = ResidueTypeFinder( *restype_set).name3( "ACE" ).get_representative_type();
	ResidueTypeCOP rt = ResidueTypeFinder( *restype_set).name3( rtname ).get_representative_type(); ResidueTypeCOP nme = ResidueTypeFinder( *restype_set).name3( "NME" ).get_representative_type();
	ResidueOP aceres( new Residue( ace, true ) );
	ResidueOP res( new Residue( rt, true ) );
	ResidueOP nmeres( new Residue( nme, true ) );

	pose.append_residue_by_jump( *aceres, true );
	pose.append_residue_by_bond( *res, true );
	pose.append_residue_by_bond( *nmeres, true );

	// There are going to be 72^3 measurements.
	Size const denominator = 72 * 72;
	Size numerator = 0;

	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_chi( 2, true );
	protocols::simple_moves::MinMoverOP min( new protocols::simple_moves::MinMover( mm, score_fxn, "linmin_iterated", 0.001, true ) );

	for ( Real phi = -175; phi <= 180; phi += 5 ) {
		for ( Real psi = -175; psi <= 180; psi += 5 ) {
			pose.set_torsion( id::TorsionID( 2, id::BB, 1 ), phi );
			pose.set_torsion( id::TorsionID( 2, id::BB, 2 ), psi );
			min->apply( pose );
			if ( ! bump_check( pose, bump_frac ) ) {
				TR << phi << " " << psi << "\n";
				++numerator;
			}
		}
	}

	TR << "Proportion for " << rtname << " is " << numerator << "/" << denominator << " or " << (Real)(numerator)/(Real)(denominator) << std::endl;
}

void report_proportions_accessible( Real const bump_frac ) {

	utility::vector1< std::string > rts;
	utility::vector1< std::string > alpha_rts;

	rts.push_back( "B3A" );
	rts.push_back( "B3C" );
	rts.push_back( "B3D" );
	rts.push_back( "B3E" );
	rts.push_back( "B3F" );
	rts.push_back( "B3G" );
	rts.push_back( "B3H" );
	rts.push_back( "B3I" );
	rts.push_back( "B3K" );
	rts.push_back( "B3L" );
	rts.push_back( "B3M" );
	rts.push_back( "B3N" );
	rts.push_back( "B3P" );
	rts.push_back( "B3Q" );
	rts.push_back( "B3R" );
	rts.push_back( "B3S" );
	rts.push_back( "B3T" );
	rts.push_back( "B3V" );
	rts.push_back( "B3W" );
	rts.push_back( "B3Y" );

	for ( Size ii = 1; ii <= rts.size(); ++ii ) {
		report_proportion_accessible( rts[ii], bump_frac );
	}

	alpha_rts.push_back( "ALA" );
	alpha_rts.push_back( "CYS" );
	alpha_rts.push_back( "ASP" );
	alpha_rts.push_back( "GLU" );
	alpha_rts.push_back( "PHE" );
	alpha_rts.push_back( "GLY" );
	alpha_rts.push_back( "HIS" );
	alpha_rts.push_back( "ILE" );
	alpha_rts.push_back( "LYS" );
	alpha_rts.push_back( "LEU" );
	alpha_rts.push_back( "MET" );
	alpha_rts.push_back( "ASN" );
	alpha_rts.push_back( "PRO" );
	alpha_rts.push_back( "GLN" );
	alpha_rts.push_back( "ARG" );
	alpha_rts.push_back( "SER" );
	alpha_rts.push_back( "THR" );
	alpha_rts.push_back( "VAL" );
	alpha_rts.push_back( "TRP" );
	alpha_rts.push_back( "TYR" );

	for ( Size ii = 1; ii <= alpha_rts.size(); ++ii ) {
		report_proportion_accessible_alpha( alpha_rts[ii], bump_frac );
	}
}

OPT_KEY( Real, bump_frac )

int
main( int argc, char* argv[] )
{
	try {
		NEW_OPT( bump_frac, "fraction of LJ radius", 0.6 );

		devel::init(argc, argv);
		report_proportions_accessible( option[ bump_frac ].value() );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
