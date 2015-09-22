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
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/simple_moves/hbs/HbsRandomSmallMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

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
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>

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
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds hbs patches to the given pdb strucure

// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR( "MikeLinkerMover" );

namespace azide_link_creator {
IntegerOptionKey const lys_one ( "azide_link_creator::lys_one" );
IntegerOptionKey const lys_two ( "azide_link_creator::lys_two" );
StringVectorOptionKey const lys_one_pdb( "azide_link_creator::lys_one_pdb" );
StringVectorOptionKey const lys_two_pdb( "azide_link_creator::lys_two_pdb" );
BooleanOptionKey const min_all( "azide_link_creator::min_all" );
}

class MikeLinkerMover : public Mover {

public:

	//default ctor
	MikeLinkerMover(): Mover("MikeLinkerMover"){}

	//default dtor
	virtual ~MikeLinkerMover(){}

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "MikeLinkerMover"; }

};

typedef utility::pointer::shared_ptr< MikeLinkerMover > MikeLinkerMoverOP;
typedef utility::pointer::shared_ptr< MikeLinkerMover const > MikeLinkerMoverCOP;


class TorsionVectorMover : public Mover {
public:
	//default ctor
	TorsionVectorMover( utility::vector1< utility::vector1< AtomID > > tv ): Mover("TorsionVectorMover"){ tv_ = tv; }

	//default dtor
	virtual ~TorsionVectorMover(){}

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "TorsionVectorMover"; }

	utility::vector1< utility::vector1< AtomID > > tv_;

};

typedef utility::pointer::shared_ptr< TorsionVectorMover > TorsionVectorMoverOP;
typedef utility::pointer::shared_ptr< TorsionVectorMover const > TorsionVectorMoverCOP;


void
TorsionVectorMover::apply( core::pose::Pose & pose ) {

	for ( Size i = 1; i <= 10; ++i ) {
		Real val = numeric::random::rg().uniform();
		Size num_options = tv_.size();
		Size index = Size( val * num_options ) + 1;

		Real start = pose.conformation().torsion_angle( tv_[index][1],
			tv_[index][2],
			tv_[index][3],
			tv_[index][4] );

		if ( numeric::random::rg().uniform() < 0.5 ) {
			start -= numeric::NumericTraits<float>::pi()/3;
		} else {
			start += numeric::NumericTraits<float>::pi()/3;
		}

		pose.conformation().set_torsion_angle( tv_[index][1],
			tv_[index][2],
			tv_[index][3],
			tv_[index][4],
			start );

	}
}


int
main( int argc, char* argv[] )
{
	try {

		utility::vector1< core::Size > empty_vector(0);
		option.add( azide_link_creator::lys_one, "lysine number one" ).def(1);
		option.add( azide_link_creator::lys_two, "lysine number two" ).def(2);
		option.add( azide_link_creator::lys_one_pdb, "lysine number one, pdb notation" ).def("A 1");
		option.add( azide_link_creator::lys_two_pdb, "lysine number two, pdb notation" ).def("A 2");
		option.add( azide_link_creator::min_all, "minimize everything" ).def(false);


		devel::init(argc, argv);

		//create mover instance
		MikeLinkerMoverOP ML_mover( new MikeLinkerMover() );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( ML_mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

void
enforce_atoms_coplanar(
	utility::vector1< id::AtomID > atoms,
	core::pose::Pose & pose
) {

	using namespace core;
	using namespace scoring;
	using namespace constraints;
	using namespace func;

	CircularHarmonicFuncOP dih_func( new CircularHarmonicFunc( 0, 0.02 ) );
	CircularHarmonicFuncOP ang_func( new CircularHarmonicFunc(
		numeric::NumericTraits<float>::pi() * ( 1-float(float(180*(atoms.size()-2))/float(atoms.size()))), 0.02 ) );

	for ( Size i = 1; i <= atoms.size(); ++i ) {
		DihedralConstraintOP d1( new DihedralConstraint( atoms[ ( i - 1 ) % atoms.size() + 1 ], atoms[ ( i ) % atoms.size() + 1 ], atoms[ ( i + 1 ) % atoms.size() + 1 ], atoms[ ( i + 2 ) % atoms.size() + 1 ], dih_func ) );
		TR << "angle constraint on " << atoms[ ( i - 1 ) % atoms.size() + 1 ] << ", " << atoms[ ( i ) % atoms.size() + 1 ] << ", " << atoms[ ( i + 1 ) % atoms.size() + 1 ] << std::endl;
		AngleConstraintOP a1( new AngleConstraint( atoms[ ( i - 1 ) % atoms.size() + 1 ], atoms[ ( i ) % atoms.size() + 1 ], atoms[ ( i + 1 ) % atoms.size() + 1 ], ang_func ) );

		pose.add_constraint( d1 );
		pose.add_constraint( a1 );
	}
}

void
enforce_triazole_distance(
	utility::vector1< id::AtomID > atoms,
	core::pose::Pose & pose
) {

	using namespace core;
	using namespace scoring;
	using namespace constraints;
	using namespace func;

	HarmonicFuncOP harm_func_1_385( new HarmonicFunc( 1.385, .1 ) );
	HarmonicFuncOP harm_func_2_285( new HarmonicFunc( 2.285, .1 ) );

	for ( Size i = 1; i <= atoms.size(); ++i ) {

		AtomPairConstraintOP ap1( new AtomPairConstraint( atoms[ ( i - 1 ) % atoms.size() + 1 ], atoms[ ( i ) % atoms.size() + 1 ], harm_func_1_385 ) );
		AtomPairConstraintOP ap2( new AtomPairConstraint( atoms[ ( i - 1 ) % atoms.size() + 1 ], atoms[ ( i + 1 ) % atoms.size() + 1 ], harm_func_2_285 ) );

		pose.add_constraint( ap1 );
		pose.add_constraint( ap2 );

	}
}

void
add_triazole_constraints(
	core::pose::Pose & pose,
	Size az1,
	Size az2,
	Size linker
) {

	using namespace core;
	using namespace scoring;
	using namespace constraints;
	using namespace func;
	using namespace id;

	HarmonicFuncOP harm_func( new HarmonicFunc( 1.385, .1 ) );
	HarmonicFuncOP harm_func2( new HarmonicFunc( 2.285, .1 ) );

	//HarmonicFuncOP harm_func_0( new HarmonicFunc( 0, std ) );
	//CircularHarmonicFuncOP dih_func_0( new CircularHarmonicFunc( 0, 0.1 ) );
	CircularHarmonicFuncOP dih_func_180( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi(), 0.1 ) );
	CircularHarmonicFuncOP ang_func_120( new CircularHarmonicFunc( numeric::NumericTraits<float>::pi_2_over_3(), 0.1 ) );

	AtomID CE1( pose.residue( az1 ).atom_index("CE"), az1 );
	AtomID NZ1( pose.residue( az1 ).atom_index("NZ"), az1 );
	AtomID NT1( pose.residue( az1 ).atom_index("NT"), az1 );
	AtomID NI1( pose.residue( az1 ).atom_index("NI"), az1 );
	AtomID CK( pose.residue( linker ).atom_index("CK"), linker );
	AtomID CT1( pose.residue( linker ).atom_index("CT1"), linker );
	AtomID CI1( pose.residue( linker ).atom_index("CI1"), linker );


	AtomID CE2( pose.residue( az2 ).atom_index("CE"), az2 );
	AtomID NZ2( pose.residue( az2 ).atom_index("NZ"), az2 );
	AtomID NT2( pose.residue( az2 ).atom_index("NT"), az2 );
	AtomID NI2( pose.residue( az2 ).atom_index("NI"), az2 );
	AtomID CM( pose.residue( linker ).atom_index("CM"), linker );
	AtomID CT2( pose.residue( linker ).atom_index("CT2"), linker );
	AtomID CI2( pose.residue( linker ).atom_index("CI2"), linker );

	utility::vector1< AtomID > r1;
	r1.push_back( NZ1 );
	r1.push_back( NT1 );
	r1.push_back( NI1 );
	r1.push_back( CI1 );
	r1.push_back( CT1 );
	enforce_atoms_coplanar( r1, pose );
	enforce_triazole_distance( r1, pose );

	utility::vector1< AtomID > r2;
	r2.push_back( NZ2 );
	r2.push_back( NT2 );
	r2.push_back( NI2 );
	r2.push_back( CI2 );
	r2.push_back( CT2 );
	enforce_atoms_coplanar( r2, pose );
	enforce_triazole_distance( r2, pose );

	/*
	AtomPairConstraintOP NT1_CT1_atompair( new AtomPairConstraint( NT1, CT1, FuncOP( new TopOutFunc( 10, 2.285, 50 ) ) ) );
	AtomPairConstraintOP NT1_CI1_atompair( new AtomPairConstraint( NT1, CI1, FuncOP( new TopOutFunc( 10, 2.285, 50 ) ) ) );
	AtomPairConstraintOP NT2_CT2_atompair( new AtomPairConstraint( NT2, CT2, FuncOP( new TopOutFunc( 10, 2.285, 50 ) ) ) );
	AtomPairConstraintOP NT2_CI2_atompair( new AtomPairConstraint( NT2, CI2, FuncOP( new TopOutFunc( 10, 2.285, 50 ) ) ) );

	AtomPairConstraintOP NZ1_CT1_atompair( new AtomPairConstraint( NZ1, CT1, FuncOP( new TopOutFunc( 10, 1.385, 50 ) ) ) );
	AtomPairConstraintOP NI1_CI1_atompair( new AtomPairConstraint( NI1, CI1, FuncOP( new TopOutFunc( 10, 1.385, 50 ) ) ) );
	AtomPairConstraintOP NZ2_CT2_atompair( new AtomPairConstraint( NZ2, CT2, FuncOP( new TopOutFunc( 10, 1.385, 50 ) ) ) );
	AtomPairConstraintOP NI2_CI2_atompair( new AtomPairConstraint( NI2, CI2, FuncOP( new TopOutFunc( 10, 1.385, 50 ) ) ) );
	*/

	DihedralConstraintOP exo11( new DihedralConstraint( CE1, NZ1, CT1, CI1, dih_func_180 ) );
	DihedralConstraintOP exo21( new DihedralConstraint( CE2, NZ2, CT2, CI2, dih_func_180 ) );
	DihedralConstraintOP exo12( new DihedralConstraint( CM, CI1, CT1, NZ2, dih_func_180 ) );
	DihedralConstraintOP exo22( new DihedralConstraint( CK, CI2, CT2, NZ1, dih_func_180 ) );

	/*
	pose.add_constraint( NZ1_CT1_atompair );
	pose.add_constraint( NI1_CI1_atompair );
	pose.add_constraint( NZ2_CT2_atompair );
	pose.add_constraint( NI2_CI2_atompair );

	pose.add_constraint( NT1_CT1_atompair );
	pose.add_constraint( NT1_CI1_atompair );
	pose.add_constraint( NT2_CT2_atompair );
	pose.add_constraint( NT2_CI2_atompair );
	*/

	pose.add_constraint( exo11 );
	pose.add_constraint( exo21 );
	pose.add_constraint( exo12 );
	pose.add_constraint( exo22 );
}

utility::vector1< utility::vector1< AtomID > >
get_linker_dihedrals( Pose pose, Size resi ) {
	AtomID CI1( pose.residue( resi ).atom_index("CI1"), resi );
	AtomID CT1( pose.residue( resi ).atom_index("CT1"), resi );
	AtomID CM( pose.residue( resi ).atom_index("CM"), resi );
	AtomID OL( pose.residue( resi ).atom_index("OL"), resi );
	AtomID CK( pose.residue( resi ).atom_index("CK"), resi );
	AtomID CT2( pose.residue( resi ).atom_index("CT2"), resi );
	AtomID CI2( pose.residue( resi ).atom_index("CI2"), resi );

	utility::vector1< AtomID > t1;
	utility::vector1< utility::vector1< AtomID > > ret;

	t1.push_back( CI1 ); t1.push_back( CT1 ); t1.push_back( CM ); t1.push_back( OL );
	ret.push_back( t1 );
	t1.clear();

	t1.push_back( CT1 ); t1.push_back( CM ); t1.push_back( OL ); t1.push_back( CK );
	ret.push_back( t1 );
	t1.clear();

	t1.push_back( CM ); t1.push_back( OL ); t1.push_back( CK ); t1.push_back( CT2 );
	ret.push_back( t1 );
	t1.clear();

	t1.push_back( OL ); t1.push_back( CK ); t1.push_back( CT2 ); t1.push_back( CI2 );
	ret.push_back( t1 );
	t1.clear();

	return ret;
}

core::kinematics::FoldTree double_unlinked_fold_tree(
	core::pose::Pose pose,
	core::Size r1,
	core::Size r2
) {
	core::Size chain_end_1 = pose.conformation().chain_endings()[1];
	core::Size chain_end_2 = pose.total_residue() - 1;
	core::kinematics::FoldTree f = pose.fold_tree();
	f.clear();

	f.add_edge( 1, r1, -1 );
	f.add_edge( r1, chain_end_1, -1 );
	f.add_edge( r1, pose.total_residue(), 1 );
	f.add_edge( pose.total_residue(), r2, 2 );
	f.add_edge( r2, chain_end_1 + 1, -1 );
	f.add_edge( r2, chain_end_2, -1 );

	TR << f;

	return f;
}

core::kinematics::FoldTree double_linked_fold_tree(
	core::pose::Pose pose,
	core::Size r1,
	core::Size r2
) {
	core::Size chain_end_1 = pose.conformation().chain_endings()[1];
	core::Size chain_end_2 = pose.total_residue() - 1;
	core::kinematics::FoldTree f = pose.fold_tree();
	f.clear();

	f.add_edge( 1, r1, -1 );
	f.add_edge( r1, chain_end_1, -1 );
	f.add_edge( r1, pose.total_residue(), "NZ", "CT1" );
	f.add_edge( pose.total_residue(), r2, "CT2", "NZ" );
	f.add_edge( r2, chain_end_1 + 1, -1 );
	f.add_edge( r2, chain_end_2, -1 );

	TR << f;

	return f;
}

void
MikeLinkerMover::apply(
	core::pose::Pose & pose
)
{
	using namespace core;
	using namespace scoring;
	using namespace constraints;
	using namespace func;
	using namespace id;
	using namespace basic;
	using namespace basic::options;

	chemical::ResidueTypeSetCOP rts =
		chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	Size lys1;
	if ( option[ azide_link_creator::lys_one ].user() ) {
		lys1 = option[ azide_link_creator::lys_one ].value();
	} else {
		lys1 = pose.pdb_info()->pdb2pose( *option[ azide_link_creator::lys_one_pdb ].value()[1].c_str(), atoi(option[ azide_link_creator::lys_one_pdb ].value()[2].c_str() ) );
	}
	Size lys2;
	if ( option[ azide_link_creator::lys_two ].user() ) {
		lys2 = option[ azide_link_creator::lys_two ].value();
	} else {
		lys2 = pose.pdb_info()->pdb2pose( *option[ azide_link_creator::lys_two_pdb ].value()[1].c_str(), atoi(option[ azide_link_creator::lys_two_pdb ].value()[2].c_str() ) );
	}


	ResidueOP res1 = core::conformation::ResidueFactory::create_residue( rts->name_map( "ALT" ) );
	ResidueOP res2 = core::conformation::ResidueFactory::create_residue( rts->name_map( "PET" ) );

	pose.replace_residue( lys1, *res1, true );
	pose.replace_residue( lys2, *res1, true );
	//pose.append_residue_by_bond( *res2, false, 3, lys2, 3 );
	pose.append_residue_by_jump( *res2, 1 );


	pose.conformation().declare_chemical_bond( lys1, "NZ", pose.total_residue(), "CT1" );
	pose.conformation().declare_chemical_bond( lys1, "NI", pose.total_residue(), "CI1" );
	pose.conformation().declare_chemical_bond( lys2, "NZ", pose.total_residue(), "CT2" );
	pose.conformation().declare_chemical_bond( lys2, "NI", pose.total_residue(), "CI2" );

	core::Size final = pose.total_residue();
	//protocols::rigid::RigidBodyTransMover trans_mover( pose, 1 );
	//trans_mover.step_size( 1 );
	//trans_mover.apply( pose );

	// We just did some disruptive stuff, so get the foldtree back

	//pose.conformation().detect_bonds();
	//core::import_pose::set_reasonable_fold_tree( pose );

	// create score function
	ScoreFunctionOP score_fxn = get_score_function();//( ScoreFunctionFactory::create_score_function( "talaris2013_cart" ) );

	Real apc_wt = 0.1;

	score_fxn->set_weight( core::scoring::atom_pair_constraint, apc_wt );
	score_fxn->set_weight( core::scoring::angle_constraint, apc_wt );
	score_fxn->set_weight( core::scoring::dihedral_constraint, apc_wt );

	add_triazole_constraints( pose, lys1, lys2, final );

	pose.dump_pdb( "start.pdb" );

	//pose.fold_tree( double_unlinked_fold_tree( pose, lys1, lys2 ) );

	pose.fold_tree( double_linked_fold_tree( pose, lys1, lys2 ) );

	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	if ( option[ azide_link_creator::min_all ].value() ) {
		movemap->set_chi_true_range( 1, final );
		movemap->set_bb_true_range( 1, final );
		//movemap->set_bb( lys2, true );
		//movemap->set_jump( 1, true );
	} else {
		movemap->set_chi( lys1, true );
		movemap->set_chi( lys2, true );
		movemap->set_bb( lys1, true );
		movemap->set_bb( lys2, true );


		movemap->set_chi( final, true );
		movemap->set_bb( final, true );

		movemap->set_jump( 1, true );
		movemap->set_jump( 2, true );
	}
	movemap->set_branches( lys1, true );
	movemap->set_branches( lys2, true );
	movemap->set_branches( pose.total_residue(), true );



	protocols::simple_moves::MinMoverOP minmover( new MinMover( movemap, score_fxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
	minmover->cartesian( true );
	minmover->apply( pose );

	pose.dump_pdb( "min.pdb" );



	/* MC LOOP */
	protocols::jd2::JobOP curr_job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// create a monte carlo object for the full cycle
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *score_fxn, 1 ) );
	moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, .8 ) );
	//rigid::RigidBodyPerturbMoverOP pert_dock_rbpm( new rigid::RigidBodyPerturbMover(1, 1,  .5) );

	simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( movemap, .3, 1 ) );
	pert_pep_small->angle_max( 'H', 2 );
	pert_pep_small->angle_max( 'L', 2 );
	pert_pep_small->angle_max( 'E', 2 );

	simple_moves::ShearMoverOP pert_pep_shear( new simple_moves::ShearMover( movemap, .3 , 1 ) );
	pert_pep_shear->angle_max( 'H', 2 );
	pert_pep_shear->angle_max( 'L', 2 );
	pert_pep_shear->angle_max( 'E', 2 );


	utility::vector1< utility::vector1< AtomID > > tv = get_linker_dihedrals( pose, final );
	TorsionVectorMoverOP pert_linker( new TorsionVectorMover( tv ) );

	moves::RandomMoverOP pert_pep_random( new moves::RandomMover() );
	pert_pep_random->add_mover( pert_pep_small, 1 );
	//pert_pep_random->add_mover( pert_linker, 1 );
	//pert_pep_random->add_mover( pert_pep_shear, 1 );

	moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( pert_pep_random, 100 ) );

	using core::pack::task::operation::TaskOperationCOP;

	// create a task factory and task operations
	TaskFactoryOP pert_tf( new TaskFactory() );
	pert_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

	operation::ReadResfileOP pert_rrop( new operation::ReadResfile() );
	pert_rrop->default_filename();
	pert_tf->push_back( pert_rrop );

	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	// create a rotamer trials mover
	simple_moves::RotamerTrialsMoverOP pert_rt( new simple_moves::EnergyCutRotamerTrialsMover( score_fxn, pert_tf, pert_mc, 0.1 /*energycut*/ ) );

	// create a random mover to hold the docking, and peptide pertubation movers
	moves::RandomMoverOP pert_random( new moves::RandomMover() );
	//pert_random->add_mover( pert_dock_rbpm, 1 );
	pert_random->add_mover( pert_pep_repeat, 0.5 );

	// create a sequence move to hold random and rotamer trials movers
	moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
	pert_sequence->add_mover( pert_random );
	pert_sequence->add_mover( pert_rt );

	// create a TrialMover for the pertubation
	moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

	for ( Size k = 1; k <= 10; ++k ) {
		pert_mc->reset(pose);
		// Perturbation phase - loop
		for ( Size j = 1; j <= 100; ++j ) {
			TR << "PERTURB: " << k << " / "  << j << std::endl;
			pert_trial->apply( pose );
			curr_job->add_string_real_pair( "ENERGY_PERT (pert score)", (*score_fxn)(pose) );
		}
		TR << "Score after perturbation " << ( *score_fxn )( pose ) << std::endl;
		pert_mc->recover_low( pose );
		curr_job->add_string_real_pair( "ENERGY_PERT (pert score) recovered low", (*score_fxn)(pose) );

		minmover->apply( pose );
		TR << "Score after minimization " << ( *score_fxn )( pose ) << std::endl;

		if ( k == 1 ) {
			mc->reset(pose);
			TR << "after mc->reset" << std::endl;
			mc->show_state();
		}

		TR << "pre mc->boltzmann" << std::endl;
		mc->show_state();
		mc->boltzmann( pose );
		TR << "post mc->boltzmann" << std::endl;
		mc->show_state();

		apc_wt += .1;
		score_fxn->set_weight( core::scoring::atom_pair_constraint, apc_wt );
		score_fxn->set_weight( core::scoring::angle_constraint, apc_wt );
		score_fxn->set_weight( core::scoring::dihedral_constraint, apc_wt );
	}

	mc->recover_low( pose );
	curr_job->add_string_real_pair( "ENERGY_FINAL (pert score) ", (*score_fxn)(pose) );
	curr_job->add_string_real_pair( "ENERGY_FINAL (hard score) ", (*score_fxn)(pose) );

	TR << "Ending main loop..." << std::endl;
	TR << "Checking pose energy..." << std::endl;

	// We just did some disruptive stuff, so get the foldtree back

	pose.conformation().detect_bonds();
	//core::import_pose::set_reasonable_fold_tree( pose );
}
