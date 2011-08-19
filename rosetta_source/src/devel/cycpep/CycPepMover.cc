/*
 * CycPepMover.cpp
 *
 *  Created on: Nov 18, 2009
 *      Author: liorz06
 */
/*
 * CycDsulf.cc
 *
 *  Created on: Mar 16, 2010
 *      Author: liorz06
 */

#include "CycPepMover.hh"

//Auto Headers
#include <core/pose/util.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.moves.CycPepMover");
using namespace core;
using namespace ObjexxFCL;

namespace protocols {
namespace CycPepMover {
CycPepMover::CycPepMover():Mover() {
	Mover::type( "CycPepMover" );

	//set defaults to scoring function
	_scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );

	//set defaults to loop relax mover
	_loop_relax_mover.remodel("perturb_kic");
	_loop_relax_mover.refine("refine_kic");
	_loop_relax_mover.fa_scorefxn();
	_loop_relax_mover.copy_sidechains(true);

	//set defaults to task repacker:
	task = pack::task::TaskFactory::create_packer_task( workpose );
	task->initialize_from_command_line().restrict_to_repacking();//flags+no-design

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
		std::cout<<resi<<std::endl;
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
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(180,1);
	id::AtomID atom1(workpose.residue(workpose.n_residue()).atom_index("CA"),workpose.n_residue());
	id::AtomID atom2(workpose.residue(workpose.n_residue()).atom_index("C"),workpose.n_residue());
	id::AtomID atom3(workpose.residue(1).atom_index("N"),1);
	id::AtomID atom4(workpose.residue(1).atom_index("CA"),1);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
	return cst_set;
}

void CycPepMover::modelSSLoop(Size startCys, Size endCys, pose::Pose& workpose){
	utility::vector1< std::pair<Size, Size> > dss;
	std::pair<Size,Size> bond(startCys,endCys);
	dss.push_back(bond);
	workpose.conformation().fix_disulfides(dss);

	protocols::loops::Loops loop;
	TR<<"Starting loop modeling protocol...";

	_loop.add_loop(startCys+1,endCys-1,(startCys+endCys)/2);
	_loop_relax_mover.loops(loop);
	_loop_relax_mover.apply(workpose);
}

void CycPepMover::apply(pose::Pose& workpose) {
	//remove terminal atoms
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,1);
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,workpose.n_residue());

	utility::vector1 < std::pair< Size, Size > > dss;
	core::conformation::disulfide_bonds(workpose,dss);
	//declare a bond between N-C terminals
	workpose.conformation().declare_chemical_bond(1,"N",workpose.n_residue(),"C");
	minimize(workpose);
	packRotamers(workpose);
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
		_scorefxn->set_weight(core::scoring::dihedral_constraint,1);
		minimize(workpose);
		dss.clear();
		core::conformation::disulfide_bonds(workpose,dss);

		modelSSLoop(ssAtoms[1],ssAtoms[2],workpose,scorefxn);
		minimize(workpose);
		packRotamers(workpose);
		scorefxn->set_weight(core::scoring::dihedral_constraint,0);
		workpose.remove_constraints(cst_set->get_all_constraints());
	}

}


}
}
