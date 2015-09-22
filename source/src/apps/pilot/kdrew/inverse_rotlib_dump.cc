// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//awatkins: this is a quick app to dump the backbone-independent rotamers of a single residue to a PDB

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_trials.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/util.hh>
//#include <core/scoring/constraints/Constraint.fwd.hh>
//#include <core/scoring/constraints/CoordinateConstraint.hh>
//#include <core/scoring/constraints/AmbiguousConstraint.hh>
//#include <core/scoring/constraints/BackboneStubConstraint.hh>
//#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/pack/rotamer_set/bb_independent_rotamers.hh>
#include <core/io/pdb/pose_io.hh>

//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
//#include <protocols/simple_moves/oop/OopPatcher.hh>
//#include <protocols/simple_moves/chiral/ChiralMover.hh>
//#include <protocols/hotspot_hashing/HotspotStubSet.hh>
//#include <protocols/hotspot_hashing/HotspotStub.hh>
//#include <protocols/rigid/RigidBodyMover.hh>
//#include <protocols/rigid/RB_geometry.hh>
//#include <protocols/ncbb/oop/OopDockDesignProtocol.hh>

// Filter headers
//#include <basic/MetricValue.hh>
//#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>

//#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
//#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <fstream>

using namespace basic::options;
using namespace protocols;


// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR( "Inverse Rotamer Library Dumper" );

namespace inverse_rotlib_dump {

	StringOptionKey const primary_hs( "inverse_rotlib_dump::primary_hs" );
	FileVectorOptionKey const ancillary_hs("inverse_rotlib_dump::ancillary_hs" );

}

class InverseRotlibDumpMover : public moves::Mover {

	public:

		//default ctor
		InverseRotlibDumpMover(): Mover("InverseRotlibDumpMover"){}
		//default dtor
		virtual ~InverseRotlibDumpMover(){}

		//methods
		virtual void apply( core::pose::Pose& pose );
		virtual std::string get_name() const { return "InverseRotlibDumpMover"; }

};

typedef utility::pointer::owning_ptr< InverseRotlibDumpMover > InverseRotlibDumpMoverOP;
typedef utility::pointer::owning_ptr< InverseRotlibDumpMover const > InverseRotlibDumpMoverCOP;


int
main( int argc, char* argv[] )
{

    try {
		option.add( inverse_rotlib_dump::primary_hs, "The path to primary hotspot stub pdb" ).def( "" );
		option.add( inverse_rotlib_dump::ancillary_hs, "The path to additional hotspot stub pdbs" ).def( "" );

		devel::init(argc, argv);

		//create mover instance
		InverseRotlibDumpMoverOP IRD_mover( new InverseRotlibDumpMover() );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( IRD_mover );

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
				return -1;
    }
    return 0;

}//main


void
InverseRotlibDumpMover::apply(core::pose::Pose& pose)
{

	//core::pose::Pose start_pose = core::pose::Pose(scaffold_pose);

	core::pose::Pose primary_hs_pose;
	core::import_pose::pose_from_pdb(primary_hs_pose, option[ inverse_rotlib_dump::primary_hs].value() );
	core::pose::add_variant_type_to_pose_residue( primary_hs_pose, "SHOVE_BB", primary_hs_pose.total_residue());
	//kdrew: create starting pose by combining the target pose and the scaffold pose

	// the hotspot residue is the final one in the file, i.e. primary_hs_pose.residue(primary_hs_pose.total_residue())

	core::conformation::ResidueCOP primary_hs_residue  =  primary_hs_pose.residue(primary_hs_pose.total_residue()).get_self_ptr();

	// structure sort of cobbled together from Florian's bb independent rotamers and invrottree dumping logic
	utility::vector1< core::conformation::ResidueCOP > primary_rots = core::pack::rotamer_set::bb_independent_rotamers (&primary_hs_residue->type());

	std::string filename("primary.pdb");
	std::vector< core::conformation::ResidueCOP >::const_iterator rot_iterator = primary_rots.begin();
        core::Size atomcounter(0), modelcount(1);
        std::ofstream file_out( filename.c_str() );
	while( rot_iterator != primary_rots.end()){
        	file_out << "MODEL   "+utility::to_string(modelcount++)+"\n";
                core::io::pdb::dump_pdb_residue( **(rot_iterator), atomcounter, file_out );
               	++rot_iterator;
		file_out << "ENDMDL \n";
        } // while( !all_res_iterators_at_end )
	file_out.close();

	//std::vector<protocols::hotspot_hashing::HotspotStubSetOP> ancillary_stubs;
	utility::vector1<std::string> ancillary_locations =  option[ inverse_rotlib_dump::ancillary_hs ]();
	for (core::Size iii = 1; iii <= ancillary_locations.size(); iii++) {

		core::pose::Pose secondary_hs_pose;
		core::import_pose::pose_from_pdb(secondary_hs_pose, ancillary_locations[iii]);
		core::pose::add_variant_type_to_pose_residue( secondary_hs_pose, "SHOVE_BB", secondary_hs_pose.total_residue());

		core::conformation::ResidueCOP secondary_hs_residue  =  &primary_hs_pose.residue(secondary_hs_pose.total_residue());

		// structure sort of cobbled together from Florian's bb independent rotamers and invrottree dumping logic

		utility::vector1< core::conformation::ResidueCOP > secondary_rots = core::pack::rotamer_set::bb_independent_rotamers (&secondary_hs_residue->type());

		std::string filename("secondary_" + utility::to_string(iii) + ".pdb");
		std::vector< core::conformation::ResidueCOP >::const_iterator rot_iterator = secondary_rots.begin();
	        core::Size atomcounter(0), modelcount(1);
	        std::ofstream file_out( filename.c_str() );
		while( rot_iterator != secondary_rots.end()){
        		file_out << "MODEL   "+utility::to_string(modelcount++)+"\n";
	                core::io::pdb::dump_pdb_residue( **(rot_iterator), atomcounter, file_out );
	               	++rot_iterator;
			file_out << "ENDMDL \n";
	        } // while( !all_res_iterators_at_end )
		file_out.close();
	}
}
