/*
 * CycPepMover.cpp
 *
 *  Created on: Nov 18, 2009
 *      Author: liorz06
 */
#include "CycPepMover.hh"

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.CycPepMover");
using namespace core;
using namespace ObjexxFCL;

namespace protocols {
namespace CycPepMover {
std::string CycPepMover::get_name() const {
	return "CycPepMover";
}
CycPepMover::CycPepMover():Mover() {
	Mover::type( "CycPepMover" );

	//set defaults to scoring function
	_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );

	//set defaults to loop relax mover
	_loop_relax_mover.remodel("perturb_kic");
	_loop_relax_mover.refine("refine_kic");
	_loop_relax_mover.fa_scorefxn(_scorefxn);
	_loop_relax_mover.relax("fullrelax");
	_loop_relax_mover.copy_sidechains(true);

	/*set defaults to task repacker: (commented out, moved to the packRotamers function since segfault is created if object not
	created inside the function (probably due to reference counting issue...)
*/

//	core::pack::task::TaskFactoryOP _packfactory = new pack::task::TaskFactory;
//	_packfactory->push_back( new pack::task::operation::InitializeFromCommandline ); // -ex1,ex2,use_input_sc,etc.
//
//	_packfactory->push_back( new pack::task::operation::IncludeCurrent ); // TODO: since our input is a prepacked structure, I always included its side-chains (these are NOT necessarily the native side-chains). But, maybe this should be left to the user (Barak)
//	_packfactory->push_back( new pack::task::operation::RestrictToRepacking ); // prevents design
}

CycPepMover::~CycPepMover() {

}

void CycPepMover::minimize(core::pose::Pose& workpose){
	core::kinematics::MoveMapOP moveMap = new core::kinematics::MoveMap();
	moveMap->set_bb_true_range(1,workpose.n_residue());

	protocols::moves::MinMover minimizer(moveMap, _scorefxn, "dfpmin_armijo_atol", 0.0001, true /*nb_list*/ );
	minimizer.apply(workpose);
}

void CycPepMover::packRotamers(core::pose::Pose& workpose){
	core::pack::task::TaskFactoryOP _packfactory = new pack::task::TaskFactory;
	_packfactory->push_back( new pack::task::operation::InitializeFromCommandline ); // -ex1,ex2,use_input_sc,etc.

	_packfactory->push_back( new pack::task::operation::IncludeCurrent ); // TODO: since our input is a prepacked structure, I always included its side-chains (these are NOT necessarily the native side-chains). But, maybe this should be left to the user (Barak)
	_packfactory->push_back( new pack::task::operation::RestrictToRepacking ); // prevents design
	_packTask = _packfactory->create_task_and_apply_taskoperations(workpose);
	TR<<"Starting repacking..."<<std::endl;
	_preventer.include_residue(1);
	_preventer.include_residue(workpose.n_residue());
	_preventer.apply(workpose,*_packTask);
	_noRepackDisulf.apply(workpose, *_packTask);
	protocols::moves::PackRotamersMoverOP packer = new protocols::moves::PackRotamersMover(_scorefxn, _packTask);
	packer->apply(workpose);
}

void CycPepMover::rotateUntilCys(core::pose::Pose& workpose, Size untilCys){
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.n_residue());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);

		Size counter = 1;
		pose::Pose tempPose(workpose,untilCys,untilCys);
		core::pose::remove_lower_terminus_type_from_pose_residue(tempPose,tempPose.n_residue());
		core::pose::remove_upper_terminus_type_from_pose_residue(tempPose,1);
	for( Size j=untilCys+1; j%workpose.n_residue()!=untilCys; ++j){

		//handle the case where j % n_residue == 0.. we actually want it to point to the last residue
		Size resi = (j%workpose.n_residue() == 0)?workpose.n_residue():j%workpose.n_residue();
		tempPose.append_polymer_residue_after_seqpos(workpose.residue(resi),counter++,false);
	}
	workpose = tempPose;
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.n_residue());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);
	workpose.conformation().declare_chemical_bond(1,"N",workpose.n_residue(),"C");
}
void CycPepMover::updateSSAtoms(pose::Pose& workpose, Size upNum, utility::vector1_int& vec, Size pepsize) {
	Size toAdd = pepsize+1-upNum;
	utility::vector1<std::pair<Size,Size> > dss;
	for(Size i=1; i<=vec.size(); ++i){
		vec[i] = (vec[i] + toAdd) % pepsize;
	}
	std::sort(vec.begin(),vec.end());
	for (Size i=1; i<=vec.size(); i+=2) {
		std::pair<Size,Size> p(vec[i],vec[i+1]);
		dss.push_back(p);
	}

	workpose.conformation().fix_disulfides(dss);
}

core::scoring::constraints::ConstraintSetOP CycPepMover::createDihedralConstraint(pose::Pose& workpose) {
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(180,3);
	id::AtomID atom1(workpose.residue(workpose.n_residue()).atom_index("CA"),workpose.n_residue());
	id::AtomID atom2(workpose.residue(workpose.n_residue()).atom_index("C"),workpose.n_residue());
	id::AtomID atom3(workpose.residue(1).atom_index("N"),1);
	id::AtomID atom4(workpose.residue(1).atom_index("CA"),1);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
	return cst_set;
}
Real CycPepMover::scoreNoConstraint(pose::Pose& workpose){
	_scorefxn->set_weight(core::scoring::dihedral_constraint,0);
	Real sc =  _scorefxn->score(workpose);
	_scorefxn->set_weight(core::scoring::dihedral_constraint,1);
	return sc;
}
void CycPepMover::modelSSLoop(Size startCys, Size endCys, pose::Pose& workpose){
	utility::vector1< std::pair<Size, Size> > dss;
	std::pair<Size,Size> bond(startCys,endCys);
	dss.push_back(bond);
	workpose.conformation().fix_disulfides(dss);

	protocols::loops::Loops loop;
	TR<<"Starting loop modeling protocol..."<<std::endl;
	loop.add_loop(startCys+1,endCys-1,((startCys+endCys)/2));
	_loop_relax_mover.loops(loop);
	_loop_relax_mover.apply(workpose);
}

void CycPepMover::printEnergies(pose::Pose& workpose) {
	for (Size i=1; i<=workpose.n_residue(); ++i) {
		std::cout<<"=========Residue "<< i<<": "<< workpose.residue(i).name()<<"============"<<std::endl;
		workpose.energies().residue_total_energies(i).print();
	}
}

core::scoring::constraints::ConstraintSetOP CycPepMover::prolineConstraint(pose::Pose& workpose) {

	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(-65.1,4);
	id::AtomID atom1(workpose.residue(6).atom_index("C"),6);
	id::AtomID atom2(workpose.residue(7).atom_index("N"),7);
	id::AtomID atom3(workpose.residue(7).atom_index("CA"),7);
	id::AtomID atom4(workpose.residue(7).atom_index("C"),7);

	core::scoring::constraints::HarmonicFuncOP spring2 = new core::scoring::constraints::HarmonicFunc(160.0,4);
	id::AtomID atom1a(workpose.residue(6).atom_index("N"),6);
	id::AtomID atom2a(workpose.residue(6).atom_index("CA"),6);
	id::AtomID atom3a(workpose.residue(6).atom_index("C"),6);
	id::AtomID atom4a(workpose.residue(7).atom_index("N"),7);
	core::scoring::constraints::HarmonicFuncOP spring3 = new core::scoring::constraints::HarmonicFunc(-81.8,4);
	id::AtomID atom1b(workpose.residue(5).atom_index("C"),5);
	id::AtomID atom2b(workpose.residue(6).atom_index("N"),6);
	id::AtomID atom3b(workpose.residue(6).atom_index("CA"),6);
	id::AtomID atom4b(workpose.residue(6).atom_index("C"),6);

	core::scoring::constraints::HarmonicFuncOP spring4 = new core::scoring::constraints::HarmonicFunc(150.7,4);
	id::AtomID atom1c(workpose.residue(7).atom_index("N"),7);
	id::AtomID atom2c(workpose.residue(7).atom_index("CA"),7);
	id::AtomID atom3c(workpose.residue(7).atom_index("C"),7);
	id::AtomID atom4c(workpose.residue(8).atom_index("N"),8);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1a,atom2a,atom3a,atom4a,spring2));
	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1b,atom2b,atom3b,atom4b,spring3));
	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1c,atom2c,atom3c,atom4c,spring4));

	return cst_set;
}
core::scoring::constraints::ConstraintSetOP CycPepMover::IleConstraint(pose::Pose& workpose) {

	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(-8.95,0.001);
	id::AtomID atom1(workpose.residue(6).atom_index("CA"),7);
	id::AtomID atom2(workpose.residue(7).atom_index("C"),7);
	id::AtomID atom3(workpose.residue(7).atom_index("N"),8);
	id::AtomID atom4(workpose.residue(7).atom_index("CA"),8);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
	return cst_set;
}
void CycPepMover::apply(pose::Pose& workpose) {
	//remove terminal atoms
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,workpose.n_residue());
	utility::vector1 < std::pair< Size, Size > > dss;
	core::conformation::disulfide_bonds(workpose.conformation(),dss);
	//declare a bond between N-C terminals
	workpose.conformation().declare_chemical_bond(1,"N",workpose.n_residue(),"C");
//	minimize(workpose);
//	packRotamers(workpose);

	//vector of Cys atoms
	utility::vector1_int ssAtoms;
	for(Size i=1; i<=dss.size(); ++i){
		ssAtoms.push_back(dss[i].first);
		ssAtoms.push_back(dss[i].second);
	}

	std::sort(ssAtoms.begin(),ssAtoms.end());
	for(Size i=1; i<=ssAtoms.size(); ++i){
		rotateUntilCys(workpose,ssAtoms[i]);
		updateSSAtoms(workpose,ssAtoms[i],ssAtoms,workpose.n_residue());
		core::scoring::constraints::ConstraintSetOP cst_set = createDihedralConstraint(workpose);
		workpose.add_constraints(cst_set->get_all_constraints());
//		core::scoring::constraints::ConstraintSetOP cst_set2;
//		 if (i == 1 ){
//			cst_set2 = IleConstraint(workpose);
//			workpose.add_constraints(cst_set2->get_all_constraints());
//		}
		packRotamers(workpose);
		//if ( i == 1 )
		//	workpose.dump_pdb("after_repack.pdb");
		minimize(workpose);
		//if ( i == 1 )
			//workpose.dump_pdb("after_minimization.pdb");

		Real sc = scoreNoConstraint(workpose);
		TR<<"Score after minimization and repacking: "<<sc<<std::endl;
		dss.clear();
		core::conformation::disulfide_bonds(workpose.conformation(),dss);
		modelSSLoop(ssAtoms[1],ssAtoms[2],workpose);
		if ( i == 1 )
			workpose.dump_pdb("after_modeling1.pdb");
		packRotamers(workpose);

		minimize(workpose);

		sc = scoreNoConstraint(workpose);
		TR<<"Score after loopmodel + minimization + repacking: "<<sc<<std::endl;
		workpose.remove_constraints(cst_set->get_all_constraints());
	}
	_scorefxn->set_weight(core::scoring::dihedral_constraint,0);
	TR<<"Finished with a score of " << _scorefxn->score(workpose)<<std::endl;
	printEnergies(workpose);
}

}
}
