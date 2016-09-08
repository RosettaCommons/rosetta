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

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.CycPepMover" );
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
	_scorefxn = core::scoring::get_score_function();

	//set defaults to loop relax mover
	_loop_relax_mover.remodel("perturb_kic");
	_loop_relax_mover.refine("refine_kic");
	_loop_relax_mover.fa_scorefxn(_scorefxn);
	_loop_relax_mover.relax("fullrelax");
	_loop_relax_mover.copy_sidechains(true);

}

CycPepMover::~CycPepMover() = default;

void CycPepMover::minimize(core::pose::Pose& workpose){
	core::kinematics::MoveMapOP moveMap = new core::kinematics::MoveMap();
	moveMap->set_jump( 1, true );
	protocols::simple_moves::MinMover minimizer(moveMap, _scorefxn, "lbfgs_armijo_atol", 0.0001, true /*nb_list*/ );
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
	_preventer.include_residue(workpose.size());
	_preventer.apply(workpose,*_packTask);
	_noRepackDisulf.apply(workpose, *_packTask);
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover(_scorefxn, _packTask);
	packer->apply(workpose);
}

void CycPepMover::rotateUntilCys(core::pose::Pose& workpose, Size untilCys){
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.size());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);

		Size counter = 1;
		pose::Pose tempPose(workpose,untilCys,untilCys);
		core::pose::remove_lower_terminus_type_from_pose_residue(tempPose,tempPose.size());
		core::pose::remove_upper_terminus_type_from_pose_residue(tempPose,1);
	for( Size j=untilCys+1; j%workpose.size()!=untilCys; ++j){

		//handle the case where j % n_residue == 0.. we actually want it to point to the last residue
		Size resi = (j%workpose.size() == 0)?workpose.size():j%workpose.size();
		tempPose.append_polymer_residue_after_seqpos(workpose.residue(resi),counter++,false);
	}
	workpose = tempPose;
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.size());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);
	workpose.conformation().declare_chemical_bond(1,"N",workpose.size(),"C");
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
    const core::pose::Pose::Residue & lastRes = workpose.residue(workpose.size());
    const core::pose::Pose::Residue & firstRes = workpose.residue(1);

    id::AtomID atom1(lastRes.atom_index("CA"),workpose.size());
	id::AtomID atom2(lastRes.atom_index("C"),workpose.size());
    id::AtomID atom3(firstRes.atom_index("N"),1);
	id::AtomID atom4(firstRes.atom_index("CA"),1);

	//create atompair constraint:
	Real distance = firstRes.xyz("N").distance(lastRes.xyz("C"));
	TR<<"Sequence: " << workpose.sequence()<<std::endl;
	TR<<"Distance is: " <<distance <<std::endl;
	core::scoring::constraints::HarmonicFuncOP springAtomPair = new core::scoring::constraints::HarmonicFunc(distance,0.01);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
//	cst_set->add_constraint(new core::scoring::constraints::AtomPairConstraint(atom3,atom2,springAtomPair));
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
	for (Size i=1; i<=workpose.size(); ++i) {
		std::cout<<"=========Residue "<< i<<": "<< workpose.residue(i).name()<<"============"<<std::endl;
		workpose.energies().residue_total_energies(i).print();
	}
}

void CycPepMover::apply(pose::Pose& workpose) {
	//remove terminal atoms
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,workpose.size());
	utility::vector1 < std::pair< Size, Size > > dss;
	core::conformation::disulfide_bonds(workpose.conformation(),dss);
	//declare a bond between N-C terminals
	workpose.conformation().declare_chemical_bond(1,"N",workpose.size(),"C");
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
		updateSSAtoms(workpose,ssAtoms[i],ssAtoms,workpose.size());
		core::scoring::constraints::ConstraintSetOP cst_set = createDihedralConstraint(workpose);
		workpose.remove_constraints();
		workpose.add_constraints(cst_set->get_all_constraints());
		packRotamers(workpose);
		std::stringstream ss;
		ss << "after_repack" << i <<".pdb";
		workpose.dump_pdb(ss.str());
		Size numCst = workpose.constraint_set()->get_all_constraints().vector().size();
		TR<<"Constraints: " << numCst<< std::endl;
		minimize(workpose);
		ss.str("");
		ss<< "after_minimization" << i << ".pdb";
		workpose.dump_pdb(ss.str());

		Real sc = scoreNoConstraint(workpose);
		TR<<"Score after minimization and repacking: "<<sc<<std::endl;
		dss.clear();
		core::conformation::disulfide_bonds(workpose.conformation(),dss);
		modelSSLoop(ssAtoms[1],ssAtoms[2],workpose);
		if ( i == 1 )
			workpose.dump_pdb("after_modeling1.pdb");

		sc = scoreNoConstraint(workpose);
		TR<<"Score after loopmodel + minimization + repacking: "<<sc<<std::endl;
		workpose.remove_constraints(cst_set->get_all_constraints());
	}
		packRotamers(workpose);

		minimize(workpose);
	
	_scorefxn->set_weight(core::scoring::dihedral_constraint,0);
	TR<<"Finished with a score of " << _scorefxn->score(workpose)<<std::endl;
}

}
}
