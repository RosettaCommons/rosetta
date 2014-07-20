// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/bder/ZincSiteEvaluator.cc
/// @brief
/// @details
/// @author Bryan Der

#include <devel/init.hh>
#include <devel/metal_interface/MetalSiteResidue.hh>
#include <devel/metal_interface/ZincSiteFinder.hh>
#include <devel/metal_interface/AddZincSiteConstraints.hh>
#include <devel/metal_interface/LigandBurial.hh>
#include <devel/metal_interface/ZincSecondShell.hh>


#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <sstream>


#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh> //degrees-radians
#include <numeric/xyz.functions.hh>




//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.bder.ZincSiteEvaluator");

typedef numeric::xyzVector<core::Real> point;
typedef point axis;
typedef core::pose::Pose Pose;

using namespace core;

namespace local{
basic::options::BooleanOptionKey const evaluate_coordinating_res_types( "evaluate_coordinating_res_types" );
basic::options::BooleanOptionKey const evaluate_coordinating_res_energies( "evaluate_coordinating_res_energies" );
basic::options::BooleanOptionKey const evaluate_zinc_geometry( "evaluate_zinc_geometry" );
basic::options::BooleanOptionKey const evaluate_zinc_neighbors( "evaluate_zinc_neighbors" );
basic::options::BooleanOptionKey const evaluate_ligand_sasa( "evaluate_ligand_sasa" );
basic::options::BooleanOptionKey const evaluate_zinc_second_shell( "evaluate_zinc_second_shell" );
basic::options::BooleanOptionKey const evaluate_mutations_native( "evaluate_mutations_native" );
basic::options::FileOptionKey const evaluate_mutations_native_pdb( "evaluate_mutations_native_pdb" );


basic::options::StringOptionKey const ligand_3_letter_code( "ligand_3_letter_code" );
basic::options::IntegerOptionKey const zinc_resnum( "zinc_resnum" );
basic::options::IntegerOptionKey const n_ligands( "n_ligands" );



}//local


///@brief
class ZincSiteEvaluator : public protocols::moves::Mover {
public:
  ZincSiteEvaluator()
  {
  }
  virtual ~ZincSiteEvaluator(){};


	virtual
	bool
	find_zinc_site( Pose const & pose ) {

		bool found_zinc_site( false );

		if( ! basic::options::option[local::zinc_resnum].user() ) {
			TR << "No zinc_resnum specified by user, finding zinc site" << std::endl;
			devel::metal_interface::ZincSiteFinderOP find_zinc = new devel::metal_interface::ZincSiteFinder();
			find_zinc->set_expecting_n_ligands( basic::options::option[local::n_ligands] );
			utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr( find_zinc->find_zinc_site( pose ) );
			found_zinc_site = !find_zinc->check_for_parse_error();
			msr_ = msr;
		}
		else {
			TR << "User specified zinc_resnum: " << basic::options::option[local::zinc_resnum] << ". Finding zinc site" << std::endl;
			devel::metal_interface::ZincSiteFinderOP find_zinc = new devel::metal_interface::ZincSiteFinder( basic::options::option[local::zinc_resnum] );
			find_zinc->set_expecting_n_ligands( basic::options::option[local::n_ligands] );
			utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr( find_zinc->find_zinc_site( pose ) );
			found_zinc_site = !find_zinc->check_for_parse_error();
			msr_ = msr;
		}

		return found_zinc_site;
	}


	virtual
	void
	score_pose( Pose & pose ) {
		using namespace core::scoring;
		ScoreFunctionOP scorefxn = get_score_function();
		scorefxn_ = scorefxn;
		scorefxn_->score( pose );
		return;
	}



	virtual
	void
	evaluate_zinc_geometry( Pose & pose ) {

		devel::metal_interface::AddZincSiteConstraintsOP zinc_constraints = new devel::metal_interface::AddZincSiteConstraints( msr_ );
		zinc_constraints->add_constraints( pose );
		zinc_constraints->evaluate_constraints( pose );
		zinc_constraints->view_constraints_in_pymol( pose );
		//zinc_constraints->output_constraints_file( pose );

		return;
	}


	virtual
	void
	evaluate_coordinating_res_types( Pose const & pose ) {
		Size num_HIS( 0 );
		Size num_CYS( 0 );
		Size num_ASP( 0 );
		Size num_GLU( 0 );

		for(Size i(2); i <= msr_.size(); ++i) {
			if(msr_[i]->get_resname() == "HIS") { ++num_HIS; }
			if(msr_[i]->get_resname() == "CYS") { ++num_CYS; }
			if(msr_[i]->get_resname() == "ASP") { ++num_ASP; }
			if(msr_[i]->get_resname() == "GLU") { ++num_GLU; }
		}

		TR << pose.pdb_info()->name() << "  zinc coordination number: " << num_HIS + num_CYS + num_ASP + num_GLU << "  zinc coordinating residues: H-C-D-E: " << num_HIS << "-" << num_CYS << "-" << num_ASP << "-" << num_GLU << std::endl;

		return;
	}


	virtual
	void
	evaluate_coordinating_res_energies( Pose const & pose ) {

		//strip out the ZN-containting residue for energy evaluation
		core::pose::Pose protein_only( pose );
		protein_only.conformation().delete_residue_slow( protein_only.total_residue() ); // HIZ is last residue
		scorefxn_->score(protein_only);

		for(Size i(2); i <= msr_.size(); ++i) {

			TR << pose.pdb_info()->name()
			<< "  res1_dun: " << protein_only.energies().residue_total_energies(msr_[i]->get_seqpos())[scoring::fa_dun]
			<< "  res1_rep: " << protein_only.energies().residue_total_energies(msr_[i]->get_seqpos())[scoring::fa_rep]
			<< "  res1_total: " << protein_only.energies().residue_total_energy(msr_[i]->get_seqpos()) << std::endl;
		}


		return;
	}


	virtual
	void
	evaluate_zinc_neighbors( Pose const & pose ) {

		std::string ligand_name( basic::options::option[local::ligand_3_letter_code] );
		devel::metal_interface::LigandBurialOP ligand_burial = new devel::metal_interface::LigandBurial( pose, ligand_name );

		//NeighborsByDistance needs to know the ligand resnum, so find_ligand() comes first
		Size ligand_resnum( ligand_burial->find_ligand() );
		TR << "Ligand resnum: " << ligand_resnum << std::endl;
		ligand_burial->register_calculators();
		ligand_burial->calculate_ligand_neighbors();

		return;
	}

	virtual
	void
	evaluate_ligand_sasa( Pose const & pose ) {

		std::string ligand_name( basic::options::option[local::ligand_3_letter_code] );
		devel::metal_interface::LigandBurialOP ligand_burial = new devel::metal_interface::LigandBurial( pose, ligand_name );

		//NeighborsByDistance needs to know the ligand resnum, so find_ligand() comes first
		Size ligand_resnum( ligand_burial->find_ligand() );
		if(ligand_resnum == 0) {
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return;
		}

		ligand_burial->register_calculators();
		ligand_burial->calculate_ligand_sasa();

		return;
	}


	virtual
	void
	evaluate_zinc_second_shell( Pose const & pose ) {

		devel::metal_interface::ZincSecondShellOP zinc_second_shell = new devel::metal_interface::ZincSecondShell( pose, msr_ );
		zinc_second_shell->register_calculators();
		zinc_second_shell->fill_first_shell_atom_ids();
		zinc_second_shell->fill_second_shell_atom_ids();

		//utility::vector1< id::AtomID > first_shell( zinc_second_shell->get_first_shell_atom_ids() );
		utility::vector1< id::AtomID > second_shell( zinc_second_shell->get_second_shell_atom_ids() );



		for(Size i(1); i <= second_shell.size(); ++i) {
			TR << second_shell[i] << std::endl;
		}

		zinc_second_shell->calculate_hbonds_and_sasa( pose );
		zinc_second_shell->report_buried_unsat( second_shell );


		return;
	}


	virtual
	void
	evaluate_mutations_native( Pose const & pose, Pose const & native_pose ) {

		// //strip out the ZN-containting residue for energy evaluation
		// core::pose::Pose protein_only( pose );
		// protein_only.conformation().delete_residue_slow( protein_only.total_residue() ); // HIZ is last residue

		//could do some asserts to ensure good input

		Size counter(0);
		std::stringstream ss_mutations;
		std::stringstream ss_pymol_mutations;
		utility::file::FileName filename( pose.pdb_info()->name() );
		std::string pdbname_base = filename.base();

		ss_mutations << pose.pdb_info()->name() << "  Mutations: ";
		ss_pymol_mutations << pose.pdb_info()->name() << "  PYMOL_SEL: select " << pdbname_base << "_muts, /" << pdbname_base << "//A/"; //hardcode chain A


		for(Size i(1); i <= native_pose.total_residue(); ++i) {
			char design_res = pose.residue(i).name1();
			char native_res = native_pose.residue(i).name1();

			if( design_res != native_res ) {
				++counter;
				ss_mutations << native_res << i << design_res << " ";
				ss_pymol_mutations << i << "+";
			}
		}

		ss_mutations << std::endl;
		ss_pymol_mutations << std::endl;

		TR << "Number of mutations: " << counter << std::endl;
		TR << ss_mutations.str() << std::endl;
		TR << ss_pymol_mutations.str() << std::endl;


		return;
	}



  virtual
  void
  apply( Pose & pose ){

		TR << std::endl << std::endl;

		bool found_zinc( find_zinc_site( pose ) );


		TR << "FOUND ZINC: true/false: " << found_zinc << std::endl;

		if( !found_zinc ) {
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return;
		}

		score_pose(pose);


		//Cys torsion preference?
		// for(Size i(2); i <= msr_.size(); ++i) {
		// 	std::string atom_name = msr_[i]->get_ligand_atom_name();
		// 	if(atom_name == " SG " || atom_name == "SG") {
		// 		Real dihed = numeric::dihedral_degrees( msr_[1]->get_ligand_atom_xyz(), msr_[i]->get_ligand_atom_xyz(), pose.residue( msr_[i]->get_seqpos() ).atom("CB").xyz(), pose.residue( msr_[i]->get_seqpos() ).atom("CA").xyz() );
		// 		TR << "Dihedral " << pose.residue( msr_[i]->get_seqpos() ).name3() << pose.pdb_info()->number( msr_[i]->get_seqpos() ) << " res " << i << " dihedral=" << dihed << std::endl;
		// 	}
		// }



		//Cross-zinc torsion preference?
		// for(Size i(2); i <= msr_.size() - 1; ++i) {
		// 	for(Size j(i+1); j <= msr_.size(); ++j) {

		// 		core::id::AtomID i_pre_id( msr_[i]->get_pre_ligand_atom_id() );
		// 		core::id::AtomID j_pre_id( msr_[i]->get_pre_ligand_atom_id() );
		// 		point i_pre_xyz = pose.residue( i_pre_id.rsd() ).atom( i_pre_id.atomno() ).xyz();
		// 		point j_pre_xyz = pose.residue( j_pre_id.rsd() ).atom( j_pre_id.atomno() ).xyz();

		// 		point i_xyz( msr_[i]->get_ligand_atom_xyz() );
		// 		point zn_xyz( msr_[1]->get_ligand_atom_xyz() );
		// 		point j_xyz( msr_[j]->get_ligand_atom_xyz() );


		// 		Real dihed1 = numeric::dihedral_degrees( i_xyz, zn_xyz, j_xyz, j_pre_xyz );
		// 		Real dihed2 = numeric::dihedral_degrees( i_pre_xyz, i_xyz, zn_xyz, j_xyz );

		// 		TR << "Cross-zinc Dihedral: " << "i_pre, i, zn, j: " << dihed1 << std::endl;
		// 		TR << "Cross-zinc Dihedral: " << "i, zn, j, j_pre: " << dihed2 << std::endl;
		// 	}
		// }





		if(basic::options::option[local::evaluate_zinc_geometry]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////                 Evaluating Zinc Geometry                 /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_zinc_geometry(pose);
		}

		if(basic::options::option[local::evaluate_coordinating_res_types]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////              Evaluating Coordinating Res Types           /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_coordinating_res_types(pose);
		}

		if(basic::options::option[local::evaluate_coordinating_res_energies]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////              Evaluating Coordinating Res Energies        /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_coordinating_res_energies(pose);
		}

		if(basic::options::option[local::evaluate_zinc_neighbors]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////                 Evaluating Zinc Neighbors                 /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_zinc_neighbors(pose);
		}

		if(basic::options::option[local::evaluate_ligand_sasa]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////                 Evaluating Ligand SASA                   /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_ligand_sasa(pose);
		}

		if(basic::options::option[local::evaluate_zinc_second_shell]) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////                 Evaluating Zinc Second Shell             /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			evaluate_zinc_second_shell(pose);
		}


		if(basic::options::option[local::evaluate_mutations_native] && basic::options::option[local::evaluate_mutations_native_pdb].user() ) {
			TR << std::endl << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			TR << "////////                 Evaluating Mutations from Native         /////////" << std::endl;
			TR << "///////////////////////////////////////////////////////////////////////////" << std::endl;
			//read in native structure
			core::pose::Pose native_pose;
			core::import_pose::pose_from_pdb( native_pose, basic::options::option[local::evaluate_mutations_native_pdb].value() );
			evaluate_mutations_native(pose, native_pose);
			TR << std::endl << std::endl;
		}

    return;
  }


	virtual
	std::string
	get_name() const { return "ZincSiteEvaluator"; }



private:

	utility::vector1< devel::metal_interface::MetalSiteResidueOP > msr_;
	core::scoring::ScoreFunctionOP scorefxn_;

};

typedef utility::pointer::owning_ptr< ZincSiteEvaluator > ZincSiteEvaluatorOP;

int main( int argc, char* argv[] )
{
	try {

	using basic::options::option;
	option.add( local::ligand_3_letter_code, "ligand 3 letter code" ).def(" ZN");
	option.add( local::n_ligands, "number of expected zinc-coordinating residues" ).def(0);
	option.add( local::zinc_resnum, "zinc residue number (rosetta numbering)" ).def(0);

	option.add( local::evaluate_zinc_geometry, "evaluate zinc geometry" ).def(true);
	option.add( local::evaluate_coordinating_res_types, "evaluate coordinating res types" ).def(false);
	option.add( local::evaluate_coordinating_res_energies, "evaluate coordinating res energies" ).def(false);
	option.add( local::evaluate_zinc_neighbors, "evaluate zinc neighbors" ).def(false);
	option.add( local::evaluate_ligand_sasa, "evaluate ligand (HIZ) sasa" ).def(false);
	option.add( local::evaluate_zinc_second_shell, "evaluate zinc second shell" ).def(false);
	option.add( local::evaluate_mutations_native, "evaluate mutations native" ).def(false);
	option.add( local::evaluate_mutations_native_pdb, "evaluate mutations native pdb" ).def("native.pdb");

  devel::init(argc, argv);
  protocols::jd2::JobDistributor::get_instance()->go(new ZincSiteEvaluator);

  TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}

