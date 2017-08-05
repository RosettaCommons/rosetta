// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/dump_capped_residue.cc
/// @brief  Useful utility: given an amino/nucleic acid name, dump a polymerically capped version as PDB
/// @author Andy Watkins (amw579@stanford.edu)

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyz.functions.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <string>
#include <sstream>

using namespace core;
using namespace conformation;
using namespace core::chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

static THREAD_LOCAL basic::Tracer TR("DumpCapped");

// application specific options
namespace dumper {
StringOptionKey const residue_name( "dumper::residue_name" );
BooleanOptionKey const nopatch( "dumper::nopatch" );
BooleanOptionKey const fiveprime( "dumper::fiveprime" );
}

int
main( int argc, char* argv[] )
{
	try {
		option.add( dumper::residue_name, "Name of residue to dump." ).def("ALA");
		option.add( dumper::nopatch, "No patches (presumably, part of the name)" ).def(false);
		option.add( dumper::fiveprime, "Add a 7MG with rna_cutpoint_upper and make THIS residue a 5prime capped variant." ).def(false);
		devel::init(argc, argv);

		PoseOP pose( new Pose );

		core::chemical::ResidueTypeSetCOP residue_set_cap = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		std::string name = option[ dumper::residue_name ].value();
		ResidueType const & base_type = residue_set_cap->name_map( name );
		std::string full_name = name;
		if ( option[ dumper::fiveprime ].value() ) {
			full_name += ":5PrimeCap";
			ResidueType const & type = residue_set_cap->name_map( full_name );
			Residue res( type, true );
			pose->conformation().append_residue_by_jump( res, 1 );

			TR << pose->residue(1) << std::endl;
			ResidueType const & type2 = residue_set_cap->name_map( "7MG:rna_cutpoint_upper" );
			Residue res2( type2, true );

			TR << res2 << std::endl;
			pose->append_residue_by_atoms( res2, true, "P", 1, "ZO3'", true );
			pose->conformation().declare_chemical_bond( 1, "ZO3'", 2, "P" );
			core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
			mm->set_chi( true );
			mm->set_bb( true );
			mm->set_jump( true );
			MinMover min;
			min.movemap( mm );
			min.apply( *pose );

		} else {
			if ( !option[ dumper::nopatch ].value() ) {
				if ( base_type.is_protein() ) {
					full_name += ":MethylatedCtermProteinFull:AcetylatedNtermProteinFull";
				} else if ( base_type.is_peptoid() ) {
					full_name += ":AcetylatedNtermDimethylatedCtermPeptoidFull";
				} else if ( base_type.has_property( "RNA" ) ) {
					full_name += ":3prime5prime_methyl_phosphate";
				}
			}
			ResidueType const & type = residue_set_cap->name_map( full_name );
			Residue res( type, true );
			if ( base_type.is_protein() ) {
				for ( Size ii = 1; ii <= res.nchi(); ++ii ) {
					res.set_chi( ii, 180 );
				}
			}
			pose->conformation().append_residue_by_jump( res, 1 );

		}
		pose->dump_pdb( name+".pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
