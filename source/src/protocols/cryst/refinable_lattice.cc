// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/cryst/refinable_lattice_creator.hh>
#include <protocols/cryst/refinable_lattice.hh>

#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/cryst/cryst_rms.hh>
#include <protocols/cryst/spacegroup.hh>
#include <protocols/cryst/util.hh>
#include <protocols/cryst/wallpaper.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/viewer/viewers.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>


using basic::Error;
using basic::Warning;

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

namespace protocols {
namespace cryst {

static basic::Tracer TR("protocols.cryst.refinable_lattice");

using namespace protocols;
using namespace core;
using namespace kinematics;
using namespace scoring;
using namespace scoring::symmetry;
using namespace conformation;
using namespace conformation::symmetry;

core::Real
getMW(core::pose::Pose const & pose) {
	static std::map<std::string,core::Real> MWLOOKUP;
	if ( MWLOOKUP.size() == 0 ) {
		MWLOOKUP["C"] = 12.011;
		MWLOOKUP["F"] = 18.9984032;
		MWLOOKUP["H"] = 1.00794;
		MWLOOKUP["N"] = 14.00674;
		MWLOOKUP["O"] = 15.9994;
		MWLOOKUP["P"] = 30.973762;
		MWLOOKUP["S"] = 32.066;
		MWLOOKUP["F"] = 18.9984;
		MWLOOKUP["CL"] = 35.453;
		MWLOOKUP["BR"] = 79.904;
		MWLOOKUP["I"] = 126.904;
	}


	core::Real mw=0.0;
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		core::conformation::Residue const &rsd = pose.residue(i);
		for ( core::Size j=1; j<=rsd.natoms(); ++j ) {
			core::chemical::AtomTypeSet const &ats = rsd.type().atom_type_set();
			std::string elt = ats[rsd.atom(j).type()].element();
			runtime_assert( MWLOOKUP.find(elt) != MWLOOKUP.end() );
			mw += MWLOOKUP[elt];
		}
	}
	return mw;
}


////////////////////////////////////////////////////

void
UpdateCrystInfo::apply( core::pose::Pose & pose ) {
	auto & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	Size Ajump_=0, Bjump_=0, Cjump_=0, SUBjump_=0;

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	// find lattice jumps
	for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
		if ( jumpname == "A" ) Ajump_=i;
		else if ( jumpname == "B" ) Bjump_=i;
		else if ( jumpname == "C" ) Cjump_=i;
		else if ( jumpname == "SUB" ) SUBjump_=i;
	}
	runtime_assert( Ajump_!=0 &&  Bjump_!=0 ); // && Cjump_ != 0 );

	numeric::xyzVector< core::Real > Axform = pose.jump(Ajump_).get_translation();
	numeric::xyzVector< core::Real > Bxform = pose.jump(Bjump_).get_translation();
	numeric::xyzVector< core::Real > Cxform(-1000,0,0);

	if ( Cjump_ != 0 ) Cxform = pose.jump(Cjump_).get_translation();

	Bxform = numeric::xyzVector< core::Real >( Bxform[1], Bxform[2], Bxform[0] );
	Cxform = numeric::xyzVector< core::Real >( Cxform[2], Cxform[0], Cxform[1] );

	io::CrystInfo ci = pose.pdb_info()->crystinfo();

	bool need_angles=false;
	numeric::xyzVector<Size> grid;
	if ( Cjump_ != 0 ) {
		Spacegroup sg;
		sg.set_spacegroup(ci.spacegroup());
		if ( sg.setting() == TRICLINIC || sg.setting() == MONOCLINIC ) {
			need_angles = true;
		}
		grid = sg.get_nsubdivisions();
	} else {
		WallpaperGroup wg;
		wg.set_wallpaper_group(ci.spacegroup());
		if ( wg.setting() == wgMONOCLINIC ) {
			need_angles = true;
		}
		grid = wg.get_nsubdivisions();
	}

	Real A = Axform.length()*grid[0];
	Real B = Bxform.length()*grid[1];
	Real C = Cxform.length()*grid[2];
	ci.A(A); ci.B(B); ci.C(C);

	if ( need_angles ) {
		Axform = numeric::xyzVector< core::Real >(Axform[0], Axform[1], Axform[2]);
		Bxform = numeric::xyzVector< core::Real >(Bxform[1], Bxform[2], Bxform[0]);
		Cxform = numeric::xyzVector< core::Real >(-Cxform[2], Cxform[0], Cxform[1]);
		Real alpha = RAD2DEG*acos( Bxform.dot(Cxform) / (Bxform.length()*Cxform.length()) );
		Real beta = RAD2DEG*acos( Axform.dot(Cxform) / (Axform.length()*Cxform.length()) );
		Real gamma = RAD2DEG*acos( Axform.dot(Bxform) / (Axform.length()*Bxform.length()) );
		ci.alpha(alpha); ci.beta(beta); ci.gamma(gamma);
	}

	// origin jump root should be at (0,0,0)
	numeric::xyzVector< core::Real > k =  pose.residue( pose.fold_tree().jump_edge( SUBjump_ ).start() ).xyz("ORIG");
	pose_asu.apply_transform_Rx_plus_v( numeric::xyzMatrix<Real>::identity(),-k );

	// copy scores
	if ( pose.data().has( ( core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA ) ) ) {
		pose_asu.data().set(
			core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA,
			basic::datacache::DataCache_CacheableData::DataOP(
			pose.data().get_ptr(core::pose::datacache::CacheableDataType::ARBITRARY_FLOAT_DATA) )
		);
	}

	// now delete the VRT
	pose = pose_asu;
	pose.conformation().delete_residue_slow( pose.size() );
	pose.pdb_info()->set_crystinfo(ci);
}

// parse_my_tag
void
UpdateCrystInfo::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{

}

std::string UpdateCrystInfo::get_name() const {
	return mover_name();
}

std::string UpdateCrystInfo::mover_name() {
	return "UpdateCrystInfo";
}

void UpdateCrystInfo::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Convert back to an asymm pose with a valid CRYST1 line.", attlist );
}

std::string UpdateCrystInfoCreator::keyname() const {
	return UpdateCrystInfo::mover_name();
}

protocols::moves::MoverOP
UpdateCrystInfoCreator::create_mover() const {
	return protocols::moves::MoverOP( new UpdateCrystInfo );
}

void UpdateCrystInfoCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UpdateCrystInfo::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////

DockLatticeMover::DockLatticeMover(core::scoring::ScoreFunctionOP sf_in) {
	sf_ = sf_in->clone();
}

DockLatticeMover::DockLatticeMover() {
	sf_ = core::scoring::get_score_function();
	set_defaults();
}

void
DockLatticeMover::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	SUBjump_ = Ajump_ = Bjump_ = Cjump_ = 0; // will get updated in init(pose);
	rot_mag_ = 5;
	trans_mag_ = 1;
	chi_mag_ = 10;
	lattice_mag_ = 0; // will get updated in init(pose);
	kT_ = 2.0;
	ncycles_ = 1000;
	contact_dist_ = 16.0;
	init_score_cut_ = 0.0;

	randomize_ = false;
	min_ = false;
	min_lattice_ = false;
	final_min_ = false;
	pack_ = false;
	perturb_chi_ = false;
	verbose_ = false;
	recover_low_ = true;

	if ( option[ in::file::native ].user() ) {
		native_ = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_file( *native_, option[ in::file::native ]() , core::import_pose::PDB_file);
		MakeLatticeMover make_lattice;
		make_lattice.set_refinable_lattice( true );
		make_lattice.contact_dist( contact_dist_-1.0 ); //?
	}
}

void
DockLatticeMover::initialize(core::pose::Pose & pose, core::Real scale_contact_dist /*=1.0*/) {
	// resolve spacegroup
	io::CrystInfo ci = pose.pdb_info()->crystinfo();

	// if a spacegroup is specified, override it
	Spacegroup sg;
	std::string sgname;

	if ( spacegroup_.length() > 0 && spacegroup_ != "input" && spacegroup_ != "INPUT" ) {
		if ( spacegroup_ == "random_chiral" ) {
			utility::vector1< std::string > sg_picked( 5 );
			sg_picked[1] = "P212121"; sg_picked[2] = "P1211"; sg_picked[3] = "C121"; sg_picked[4] = "P21212";
			sg_picked[5] = ci.spacegroup();
			Size irand = Size(5.0*numeric::random::uniform())+1;
			sgname = sg_picked[irand];
			sg.set_spacegroup( sgname );
			//core::pose::add_score_line_string( pose, "Spacegrp", sgname );
			TR << "Assigning random_chiral space group of: " << sgname << std::endl;

		} else if ( spacegroup_ == "random_achiral" ) {
			utility::vector1< std::string > sg_picked( 10 );
			sg_picked[1] = "P121/c1"; sg_picked[2] = "P-1"; sg_picked[3] = "C12/c1";
			sg_picked[4] = "Pbca"; sg_picked[5] = "Pna21"; sg_picked[6] = "C1c1";
			sg_picked[7] = "Pbcn"; sg_picked[8] = "Pca21"; sg_picked[8] = "Pccn";
			sg_picked[9] ="231"; sg_picked[10] = ci.spacegroup();

			Size irand = Size(10.0*numeric::random::uniform())+1;
			sgname = sg_picked[irand];

			sg.set_spacegroup( sgname );
			TR << "Assigning random_achiral space group of: " << sg_picked[irand] << std::endl;
			//core::pose::add_score_line_string( pose, "Spacegrp", sgname );

		} else {
			sgname = spacegroup_;
			sg.set_spacegroup(spacegroup_);
		}
		std::string fullSG = sg.pdbname();
		if ( ci.spacegroup() != fullSG ) {
			randomize_ = true;
			ci.spacegroup(fullSG);
		}
	} else {
		sg.set_spacegroup(ci.spacegroup());
		sgname = ci.spacegroup();
	}

	// tgt: 50% solvent

	// make a copy of sf_ to tweak with vdw rep
	core::scoring::ScoreFunctionOP scorevdw (new core::scoring::symmetry::SymmetricScoreFunction());
	scorevdw->set_weight( fa_atr, 1.0 );
	scorevdw->set_weight( fa_rep, 0.2 );

	// now make lattice
	MakeLatticeMover make_lattice;
	make_lattice.set_refinable_lattice( true );
	make_lattice.contact_dist( scale_contact_dist*contact_dist_ );


	bool good = false;
	core::Size iter = 0;
	core::Size scorecut_increment_step = 10;
	core::Real scorecut_increment_size = 100.0;
	core::Real init_score_cut_loc( init_score_cut_ );

	while ( !good ) {
		good = true;
		if ( randomize_ ) {
			// randomize volume
			sg.set_parameters(1.0,1.0,1.0,90.0,90.0,90.0);
			core::Real mw = getMW( pose );
			core::Real tgtVol = 1.23/(0.65) * mw * sg.nsymmops();
			core::Real volscale = sg.volume();
			core::Real A = 0.2+0.8*numeric::random::rg().uniform();
			core::Real B = 0.2+0.8*numeric::random::rg().uniform();
			core::Real C = 0.2+0.8*numeric::random::rg().uniform();
			core::Real alpha = 60+60*numeric::random::rg().uniform();
			core::Real beta  = 60+60*numeric::random::rg().uniform();
			core::Real gamma = 60+60*numeric::random::rg().uniform();
			core::Real scale = std::pow( volscale*tgtVol/(A*B*C), 1.0/3.0);
			ci.A(scale*A);
			ci.B(scale*B);
			ci.C(scale*C);
			ci.alpha(alpha);
			ci.beta(beta);
			ci.gamma(gamma);
			pose.pdb_info()->set_crystinfo(ci);

			// randomize position
			numeric::xyzVector< core::Real > x( numeric::random::random_normal() );
			numeric::xyzVector< core::Real > ytmp;
			do {
				ytmp = numeric::random::random_normal();
			} while ( x.cross( ytmp ).length_squared() < 0.1 );
			numeric::xyzVector< core::Real > z( x.cross( ytmp ).normalized() ), y( z.cross(x) );
			numeric::xyzMatrix< core::Real > const R( numeric::xyzMatrix< core::Real >::cols( x, y, z ) );
			numeric::xyzVector< core::Real > const v(
				scale*A*numeric::random::uniform(), scale*B*numeric::random::uniform(), scale*C*numeric::random::uniform());
			pose.apply_transform_Rx_plus_v( R,v);

			// chis?
			if ( perturb_chi_ ) randomize_chis(pose);
		}


		make_lattice.apply( pose );

		if ( randomize_ ) {
			slide_lattice( pose );
			core::Real vdwscore = 999;
			if ( regenerate_lattice( pose ) ) {
				//vdwscore = (*sf_)(pose);
				vdwscore = (*scorevdw)(pose);
				good = (vdwscore < init_score_cut_loc); //????
				if ( !good ) {
					protocols::cryst::UpdateCrystInfo update_cryst1;
					update_cryst1.apply( pose );
				}
			} else {
				good = false; //????
			}

			TR << "Attempt " << iter++ << " score = " << vdwscore << std::endl;
		}

		if ( iter%scorecut_increment_step == 0 ) {
			init_score_cut_loc += scorecut_increment_size;
			TR << "increase score cut, now " << init_score_cut_loc << std::endl;
		}
	}

	core::Real scale = std::pow( make_lattice.symmops().size() , 1.0/3.0 );
	numeric::xyzVector<core::Size> nsub = make_lattice.sg().get_nsubdivisions();
	lattice_mag_ = scale * numeric::xyzVector<core::Real>(  1.0/nsub[0],1.0/nsub[1],1.0/nsub[2] );

	// finally, some stuff we will need later
	auto & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	symdofs_ = SymmConf.Symmetry_Info()->get_dofs();

	for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
		if ( jumpname == "SUB" ) SUBjump_=i;
		if ( jumpname == "A" ) Ajump_=i;
		if ( jumpname == "B" ) Bjump_=i;
		if ( jumpname == "C" ) Cjump_=i;
	}
	runtime_assert( SUBjump_ != 0 );

	// make sure put info of randomly selected SG
	core::pose::add_score_line_string( pose, "Spacegrp", sgname );
}

void
DockLatticeMover::slide_lattice( core::pose::Pose & pose ) {
	core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
	mm->set_jump(true); mm->set_chi(false); mm->set_bb(false);
	core::pose::symmetry::make_symmetric_movemap( pose, *mm );
	protocols::minimization_packing::MinMoverOP min_rb( new protocols::minimization_packing::MinMover(mm, sf_, "lbfgs_armijo", 0.01, true) );
	min_rb->apply( pose );
}

void
DockLatticeMover::perturb_lattice( core::pose::Pose & pose ) {
	std::map< Size, SymDof >::iterator it, it_begin = symdofs_.begin(), it_end = symdofs_.end();

	for ( it = it_begin; it != it_end; ++it ) {
		if ( it->first == SUBjump_ ) continue;

		core::Real mag = 0;
		if ( it->first == Ajump_ ) mag = lattice_mag_[0];
		if ( it->first == Bjump_ ) mag = lattice_mag_[1];
		if ( it->first == Cjump_ ) mag = lattice_mag_[2];
		runtime_assert( mag!=0 );

		core::kinematics::Jump jump_i = pose.jump( it->first );
		numeric::xyzVector< core::Real > transNorm = jump_i.get_translation().normalized();
		numeric::xyzVector< core::Real > transdel (
			mag * numeric::random::rg().gaussian() * transNorm[0],
			mag * numeric::random::rg().gaussian() * transNorm[1],
			mag * numeric::random::rg().gaussian() * transNorm[2]
		);

		jump_i.set_translation( jump_i.get_translation() - transdel );  // neg. direction expands lattice

		pose.set_jump( it->first, jump_i );
	}
}

void
DockLatticeMover::randomize_chis( core::pose::Pose & pose ) {
	core::Size nres = pose.total_residue();
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		SymmetricConformation & symm_conf( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		nres = symm_conf.Symmetry_Info()->num_independent_residues();
	}

	for ( core::Size i=1; i<=nres; ++i ) {
		for ( core::Size j=1; j<=pose.residue(i).nchi(); ++j ) {
			pose.set_chi( j, i, 360*numeric::random::rg().uniform() );
		}
	}
}

void
DockLatticeMover::perturb_chis( core::pose::Pose & pose ) {
	core::Size nres = pose.total_residue();
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		SymmetricConformation & symm_conf( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		nres = symm_conf.Symmetry_Info()->num_independent_residues();
	}

	for ( core::Size i=1; i<=nres; ++i ) {
		for ( core::Size j=1; j<=pose.residue(i).nchi(); ++j ) {
			core::Real mag = chi_mag_*numeric::random::rg().gaussian();
			pose.set_chi( j, i, pose.chi( j, i) + mag );
		}
	}
}

void
DockLatticeMover::perturb_rb( core::pose::Pose & pose ) {
	std::map< Size, SymDof >::iterator it, it_begin = symdofs_.begin(), it_end = symdofs_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		core::kinematics::Jump flexible_jump = pose.jump( it->first );
		SymDof const &dof = it->second;

		for ( Size i = 1; i<= 3; ++i ) {
			if ( dof.allow_dof(i) ) {
				flexible_jump.gaussian_move_single_rb( 1, trans_mag_, i );
			}
		}

		for ( Size i = 4; i<= 6; ++i ) {
			if ( dof.allow_dof(i) ) {
				flexible_jump.gaussian_move_single_rb( 1, rot_mag_, i );
			}
		}

		pose.set_jump( it->first, flexible_jump );
	}
}

bool
DockLatticeMover::regenerate_lattice( core::pose::Pose & pose ) {
	protocols::cryst::UpdateCrystInfo update_cryst1;
	protocols::cryst::MakeLatticeMover setup;
	setup.contact_dist( contact_dist_ );
	setup.set_refinable_lattice( true );

	update_cryst1.apply( pose );
	io::CrystInfo ci = pose.pdb_info()->crystinfo();

	// check if regeneration will mess things up
	//if ( ci.A()<0.5 || ci.B()<0.5 || ci.C()<0.5 ) return false;
	//if ( ci.alpha()<10 || ci.beta()<10 || ci.gamma()<10 ) return false;
	//if ( ci.alpha()>170 || ci.beta()>170 || ci.gamma()>170 ) return false;
	//if ( ci.A()>50.0 || ci.B()>50.0 || ci.C()>50.0 ) return false;
	//if ( ci.A()*ci.B()*ci.C() > 20000.0 ) return false;

	// fd unfortunately above checks are not general
	// instead make sure that at least 10% of volume is occupied and no more than 200%
	core::Real const deg2rad( 57.29577951308232 );
	core::Real ca = cos(deg2rad*ci.alpha()), cb = cos(deg2rad*ci.beta()), cg = cos(deg2rad*ci.gamma());
	core::Real volume = ci.A()*ci.B()*ci.C() * std::sqrt( 1-ca*ca-cb*cb-cg*cg+2*ca*cb*cg );

	// monomeric structure at this point
	core::Real MW = getMW( pose );
	core::Real occ = 1.23*MW / volume;
	if ( occ < 0.1 || occ > 2.0 ) {
		TR << "Regenerated lattice fails occupancy check [occ=" << occ << "]: ";
		TR << ci.A()<<","<<ci.B()<<","<<ci.C() <<" , "<< ci.alpha()<<","<<ci.beta()<<","<<ci.gamma()<< std::endl;
		return false;
	}

	// check valid
	core::pose::Pose pose_tmp( pose );
	//if ( !setup.check_valid( pose_tmp ) ) {
	// return false;
	//}
	setup.apply( pose );


	SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	symdofs_ = SymmConf.Symmetry_Info()->get_dofs();

	for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
		std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
		if ( jumpname == "SUB" ) SUBjump_=i;
		if ( jumpname == "A" ) Ajump_=i;
		if ( jumpname == "B" ) Bjump_=i;
		if ( jumpname == "C" ) Cjump_=i;
	}

	return true;
}

void
DockLatticeMover::repack( core::pose::Pose & pose, bool rottrials ) {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::simple_task_operations;

	TaskFactoryOP tf (new TaskFactory);
	tf->push_back( TaskOperationCOP(new InitializeFromCommandline) );
	tf->push_back( TaskOperationCOP(new IncludeCurrent) );
	tf->push_back( TaskOperationCOP(new RestrictToRepacking) );
	tf->push_back( TaskOperationCOP(new NoRepackDisulfides) );
	tf->push_back( TaskOperationCOP(new RestrictToInterface( SUBjump_ ) ) );

	if ( rottrials ) {
		protocols::minimization_packing::RotamerTrialsMoverOP pack_interface_rtrials( new protocols::minimization_packing::RotamerTrialsMover( sf_, tf ) );
		pack_interface_rtrials->apply( pose );
	} else {
		protocols::minimization_packing::PackRotamersMoverOP pack_interface_repack( new protocols::minimization_packing::PackRotamersMover( sf_ ) );
		pack_interface_repack->task_factory(tf);
		pack_interface_repack->apply( pose );
	}
}


void
DockLatticeMover::apply( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");

	// main loop!
	std::string movetag;

	// MC mover
	bool good = false;
	core::Size nfail( 0 );
	std::string sgname("");

	core::Real contact_dist_scale = 1.0;

	while ( !good ) {
		// randomize and symmetrize pose
		initialize(pose, contact_dist_scale);

		core::pose::get_score_line_string(pose, "Spacegrp", sgname );
		//TR << "Spacegrp value? " << sgname << std::endl;

		// set up minimizers
		core::kinematics::MoveMapOP mmrb(new core::kinematics::MoveMap);
		mmrb->set_jump(true); mmrb->set_chi(false); mmrb->set_bb(false);
		core::pose::symmetry::make_symmetric_movemap( pose, *mmrb );
		protocols::minimization_packing::MinMoverOP min_rb
			= protocols::minimization_packing::MinMoverOP ( new protocols::minimization_packing::MinMover(mmrb, sf_, "lbfgs_armijo", 0.01, true) );

		core::kinematics::MoveMapOP mmall(new core::kinematics::MoveMap);
		mmall->set_jump(true); mmall->set_chi(true); mmall->set_bb(true);
		core::pose::symmetry::make_symmetric_movemap( pose, *mmall );
		protocols::minimization_packing::MinMoverOP min_all
			= protocols::minimization_packing::MinMoverOP ( new protocols::minimization_packing::MinMover(mmall, sf_, "lbfgs_armijo", 0.01, true) );
		min_all->max_iter(20);

		core::Real temp = kT_ ;
		protocols::cryst::UpdateCrystInfo update_cryst1;

		(*sf_)(pose);
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(pose, *sf_, temp ) );

		for ( int i=1; i<=(int)ncycles_; ++i ) {
			if ( i % 10 == 0 ) {
				perturb_lattice(pose);
				movetag="lattice";
			} else {
				if ( perturb_chi_ && (i % 2)==0 ) {
					perturb_chis(pose);
					perturb_rb(pose);
					movetag="chi+rb";
				} else {
					perturb_rb(pose);
					movetag="rb";
				}
			}

			if ( pack_ ) {
				repack( pose, i % 10 == 0 );
			}

			if ( min_ ) {
				min_all->apply(pose);
			} else if ( min_lattice_ ) {
				min_rb->apply(pose);
			}

			(*sf_)(pose);
			mc->boltzmann( pose, movetag );

			if ( verbose_ && i%10==0 ) {
				mc->show_scores();
				mc->show_counters();
				// update_cryst1.report( pose );
			}
		}

		if ( recover_low_ ) mc->recover_low( pose );

		core::Real score_pre = (*sf_)(pose);

		if ( !regenerate_lattice(pose) ) {
			contact_dist_scale += 1.1;
			TR << "Error regenerating lattice.  Increase contact dist by factor of " << contact_dist_scale << std::endl;
			nfail++;
			if ( nfail > 10 ) { // don't try more than 10 times
				TR << "Terminate due to repeating failures in regenerate_lattices" << std::endl;
				core::pose::setPoseExtraScore( pose, "crystScore", 0.0);
				core::pose::setPoseExtraScore( pose, "crystRMS", 10.0);
				return;
			}

		} else {
			core::Real score_post = (*sf_)(pose);

			good = (std::fabs( score_post-score_pre ) < 1);

			TR << "After refinement score = " << score_pre << std::endl;
			TR << "     after regen score = " << score_post << std::endl;

			if ( !good ) {
				TR << "Rerunning!" << std::endl;
				update_cryst1.apply( pose );
			} else if ( verbose_ ) {
				core::Real RMS=-1.0;
				if ( native_ ) {
					io::CrystInfo cin = native_->pdb_info()->crystinfo();
					TR << "native Xtal: " << cin.A() << " " << cin.B() << " " << cin.C() << std::endl;
					RMS = crystRMS( pose, *native_, true);
				}
				TR << "           rms = " << RMS << std::endl;
			}
		}
	}

	if ( min_ || final_min_ ) {
		core::kinematics::MoveMapOP mmall(new core::kinematics::MoveMap);
		mmall->set_jump(true); mmall->set_chi(true); mmall->set_bb(true);
		core::pose::symmetry::make_symmetric_movemap( pose, *mmall );
		protocols::minimization_packing::MinMoverOP min_all = protocols::minimization_packing::MinMoverOP ( new
			protocols::minimization_packing::MinMover(mmall, sf_, "lbfgs_armijo", 0.01, true) );
		min_all->max_iter(20);

		min_all->apply(pose);
	}

	core::Real RMS=-1.0;
	if ( native_ ) {
		RMS = crystRMS( pose, *native_, true);
		core::Real score = (*sf_)(pose);
		core::pose::setPoseExtraScore( pose, "crystScore", score);
		core::pose::setPoseExtraScore( pose, "crystRMS", RMS);
	}

	TR << "adding score line Spacegrp: " << sgname << std::endl;
	core::pose::add_score_line_string( pose, "Spacegrp", sgname );
}


// parse_my_tag
void
DockLatticeMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	if ( tag->hasOption( "scorefxn" ) ) {
		sf_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	if ( tag->hasOption( "ncycles" ) ) {
		ncycles_ = tag->getOption<core::Size>( "ncycles" );
	}
	if ( tag->hasOption( "trans_step" ) ) {
		trans_mag_ = tag->getOption<core::Real>( "trans_step" );
	}
	if ( tag->hasOption( "kT" ) ) {
		kT_ = tag->getOption<core::Real>( "kT" );
	}
	if ( tag->hasOption( "rot_step" ) ) {
		rot_mag_ = tag->getOption<core::Real>( "rot_step" );
	}
	if ( tag->hasOption( "chi_step" ) ) {
		chi_mag_ = tag->getOption<core::Real>( "chi_step" );
	}
	if ( tag->hasOption( "init_score_cut" ) ) {
		init_score_cut_ = tag->getOption<core::Real>( "init_score_cut" );
	}

	if ( tag->hasOption( "spacegroup" ) ) {
		spacegroup_ = tag->getOption<std::string>( "spacegroup" );
	}
	if ( tag->hasOption( "randomize" ) ) {
		randomize_ = tag->getOption<bool>( "randomize" );
	}
	if ( tag->hasOption( "recover_low" ) ) {
		recover_low_ = tag->getOption<bool>( "recover_low" );
	}
	if ( tag->hasOption( "perturb_chi" ) ) {
		perturb_chi_ = tag->getOption<bool>( "perturb_chi" );
	}

	if ( tag->hasOption( "verbose" ) ) {
		verbose_ = tag->getOption<bool>( "verbose" );
	}

	if ( tag->hasOption( "min" ) ) {
		min_ = tag->getOption<bool>( "min" );
	}

	if ( tag->hasOption( "lattice_min" ) ) {
		min_lattice_ = tag->getOption<bool>( "lattice_min" );
	}

	if ( tag->hasOption( "final_min" ) ) {
		final_min_ = tag->getOption<bool>( "final_min" );
	}

}

std::string DockLatticeMover::get_name() const {
	return mover_name();
}

std::string DockLatticeMover::mover_name() {
	return "DockLatticeMover";
}

void DockLatticeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "fullatom" , xsct_rosetta_bool, "If fullatom=1 it does fullatom docking with repacks and min steps; if fullatom=0 you need to provide a centroid energy function (e.g. score4_smooth) and it will also add lattice slide moves.", "false" )
		+ XMLSchemaAttribute( "ncycles" , xsct_non_negative_integer , "# of steps to run." )
		+ XMLSchemaAttribute( "trans_step" , xsct_real, "Magnitude of translational perturbations." )
		+ XMLSchemaAttribute( "rot_step" , xsct_real, "Magnitude of rotational perturbations." )
		+ XMLSchemaAttribute( "chi_step", xsct_real, "Magnitude of torsion angle perturbations." )
		+ XMLSchemaAttribute( "kT" , xsct_real, "Simulation temperature" )
		+ XMLSchemaAttribute( "randomize", xsct_rosetta_bool,"Use a random starting point?" )
		+ XMLSchemaAttribute( "perturb_chi", xsct_rosetta_bool,"Do we also randomize the torsions?" )
		+ XMLSchemaAttribute( "spacegroup", xs_string, "What spacegroup to run? Either a spacegroup name (no spaces or one of: input, random_chiral, or random_achiral)" )
		+ XMLSchemaAttribute( "verbose", xsct_rosetta_bool,"Be verbose in output?" )
		+ XMLSchemaAttribute( "init_score_cut", xsct_real,"the initial score cut to use for lattice generation" )
		+ XMLSchemaAttribute( "recover_low", xsct_rosetta_bool,"Recover the low energy sampled conformation?" )
		+ XMLSchemaAttribute( "min", xsct_rosetta_bool,"Minimize the lattice after each perturbation?" )
		+ XMLSchemaAttribute( "final_min", xsct_rosetta_bool,"Minimize the lattice at the end of each simulation?" );

	protocols::rosetta_scripts::attributes_for_parse_score_function ( attlist ) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Docking within a crystal lattice.", attlist );
}

std::string DockLatticeMoverCreator::keyname() const {
	return DockLatticeMover::mover_name();
}

protocols::moves::MoverOP
DockLatticeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DockLatticeMover );
}

void DockLatticeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DockLatticeMover::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////

// bool
// MakeLatticeMover::check_valid( core::pose::Pose & pose ) {
//  // initialize sg_ from pose CRYST1 line
//  //bool valid = true;
//  io::CrystInfo ci = pose.pdb_info()->crystinfo();
//  runtime_assert(ci.A()*ci.B()*ci.C() != 0);  // TODO: allow these to be randomized
//
//  sg_.set_spacegroup(ci.spacegroup());
//  sg_.set_parameters(ci.A(),ci.B(),ci.C(), ci.alpha(), ci.beta(), ci.gamma());
//
//  // get principal axes
//  numeric::xyzVector< core::Real > Ax, Bx, Cx;
//  Ax = sg_.f2c()*numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
//  Bx = sg_.f2c()*numeric::xyzVector< core::Real >(0,1,0); Bx.normalize();
//  Cx = sg_.f2c()*numeric::xyzVector< core::Real >(0,0,1); Cx.normalize();
//
//  //Size rootres = place_near_origin( pose );
//  place_near_origin( pose );
//
//  core::pose::Pose posebase;
//  utility::vector1<Size> Ajumps, Bjumps, Cjumps, monomer_jumps, monomer_anchors;
//  //Size base_monomer;
//
//  numeric::xyzVector< core::Real > max_extent(0,0,0);
//  for ( Size i=1; i<= pose.size(); ++i ) {
//   if ( !pose.residue(i).is_protein() ) continue;
//   for ( Size j=1; j<=pose.residue(i).natoms(); ++j ) {
//    numeric::xyzVector< core::Real > f_ij = sg_.c2f()*pose.residue(i).xyz(j);
//    for ( int k=0; k<3; ++k ) max_extent[k] = std::max( std::fabs(f_ij[k]) , max_extent[k] );
//   }
//  }
//  max_extent += sg_.c2f()*numeric::xyzVector< core::Real >(6,6,6);
//
//  if ( max_extent[0] > 1 || max_extent[1] > 1 || max_extent[2] > 1 ||
//    max_extent[0] < 0 || max_extent[1] < 0 || max_extent[2] < 0 ) return false;
//  return true;
// }

void
MakeLatticeMover::apply( core::pose::Pose & pose ) {
	// initialize sg_ from pose CRYST1 line
	io::CrystInfo ci = pose.pdb_info()->crystinfo();
	runtime_assert(ci.A()*ci.B()*ci.C() != 0);  // TODO: allow these to be randomized

	sg_.set_spacegroup(ci.spacegroup());
	sg_.set_parameters(ci.A(),ci.B(),ci.C(), ci.alpha(), ci.beta(), ci.gamma());

	// get principal axes
	numeric::xyzVector< core::Real > Ax, Bx, Cx;
	Ax = sg_.f2c()*numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
	Bx = sg_.f2c()*numeric::xyzVector< core::Real >(0,1,0); Bx.normalize();
	Cx = sg_.f2c()*numeric::xyzVector< core::Real >(0,0,1); Cx.normalize();

	Size rootres = place_near_origin( pose );

	core::pose::Pose posebase;
	utility::vector1<Size> Ajumps, Bjumps, Cjumps, monomer_jumps, monomer_anchors;
	Size base_monomer;

	numeric::xyzVector< core::Real > max_extent(0,0,0);
	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( !pose.residue(i).is_protein() ) continue;
		for ( Size j=1; j<=pose.residue(i).natoms(); ++j ) {
			numeric::xyzVector< core::Real > f_ij = sg_.c2f()*pose.residue(i).xyz(j);
			for ( int k=0; k<3; ++k ) max_extent[k] = std::max( std::fabs(f_ij[k]) , max_extent[k] );
		}
	}
	max_extent += sg_.c2f()*numeric::xyzVector< core::Real >(6,6,6);

	// in each direction, find closest symmcenter we are not considering
	numeric::xyzVector< Size > grid = sg_.get_nsubdivisions();
	numeric::xyzVector<core::Real> closest_nongen_symmcenter(1.0/grid[0], 1.0/grid[1], 1.0/grid[2]);

	numeric::xyzVector<int> EXTEND(
		std::max( 1, (int)std::ceil( 2*max_extent[0]-closest_nongen_symmcenter[0] )),
		std::max( 1, (int)std::ceil( 2*max_extent[1]-closest_nongen_symmcenter[1] )),
		std::max( 1, (int)std::ceil( 2*max_extent[2]-closest_nongen_symmcenter[2] )) );

	TR.Debug << "using extent " << EXTEND[0] << "," << EXTEND[1] << "," << EXTEND[2] << std::endl;
	TR.Debug << "with max_extent = " << max_extent[0] << "," << max_extent[1] << "," << max_extent[2] << std::endl;
	TR.Debug << "with closest_nongen_symmcenter = " << closest_nongen_symmcenter[0] << "," << closest_nongen_symmcenter[1] << "," << closest_nongen_symmcenter[2] << std::endl;

	build_lattice_of_virtuals( posebase, EXTEND, Ajumps, Bjumps, Cjumps, monomer_anchors, base_monomer);

	Size nvrt = posebase.size();
	Size nres_monomer = pose.size();

	// only connecting ones!
	// updates monomer_anchors
	detect_connecting_subunits( pose, posebase, monomer_anchors, base_monomer );
	add_monomers_to_lattice( pose, posebase, monomer_anchors, monomer_jumps, rootres );

	Size nsubunits = monomer_anchors.size();

	core::pose::PDBInfoOP pdbinfo_old=pose.pdb_info(), pdbinfo_new;
	pose = posebase;

	conformation::symmetry::SymmetryInfo syminfo;
	setup_xtal_symminfo( pose, nsubunits, nvrt, base_monomer, nres_monomer, Ajumps, Bjumps, Cjumps, monomer_jumps, syminfo );

	bool symmdetectdisulf = basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ]();
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(false);
	pose::symmetry::make_symmetric_pose( pose, syminfo );
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(symmdetectdisulf);

	pdbinfo_new = pose::PDBInfoOP( new pose::PDBInfo( pose, true ) );
	core::pose::symmetry::make_symmetric_pdb_info( pose, pdbinfo_old, pdbinfo_new );
	pdbinfo_new->set_crystinfo(ci);
	pose.pdb_info( pdbinfo_new );

	//fpd foldtree above assumes no monomer jumps
	//fpd here, we reinitialize the foldtree to let the SymmetryInfo machinery take care of
	//    monomer jumps
	//fpd this also has the advantage of renumbering jumps in a consistent way so that
	//    all symmetric foldtree manipulations work well with symm poses made from this function
	core::kinematics::FoldTree ft = core::conformation::symmetry::get_asymm_unit_fold_tree(pose.conformation());
	core::conformation::symmetry::symmetrize_fold_tree(pose.conformation(), ft);
	pose.fold_tree( ft );

	// force symmetrization
	core::conformation::symmetry::SymmetryInfoCOP symminfo_new =
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & >( pose.conformation()).Symmetry_Info();
	for ( core::Size i=pose.fold_tree().num_jump(); i>=1; --i ) {
		if ( symminfo_new->jump_is_independent(i) && symminfo_new->jump_clones(i).size() > 0 ) {
			core::kinematics::Jump j_i = pose.jump( i );
			pose.set_jump( i, j_i );
		}
	}

	// if we are mirror symmetric, update restypes
	if ( core::conformation::symmetry::is_mirror_symmetric( pose.conformation() ) ) {
		auto & mirror_conf(
			dynamic_cast< core::conformation::symmetry::MirrorSymmetricConformation& >( pose.conformation() ) );
		mirror_conf.update_residue_identities();
	}

	// update disulf info
	pose.conformation().detect_disulfides();

}

Size
MakeLatticeMover::place_near_origin (
	Pose & pose
) {
	Size rootpos=0;
	Size nres = pose.size();

	numeric::xyzVector< core::Real > com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) {
		core::conformation::Residue const &rsd = pose.residue(i);
		com += rsd.atom(rsd.nbr_atom()).xyz();
		if ( pose.residue(i).is_upper_terminus() ) break;
	}
	com /= nres;

	Real mindis2(1e6);
	for ( Size i=1; i<= nres; ++i ) {
		core::conformation::Residue const &rsd = pose.residue(i);
		Real const dis2( com.distance_squared(  rsd.atom(rsd.nbr_atom()).xyz() ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
			rootpos = i;
		}
		if ( pose.residue(i).is_upper_terminus() ) break;
	}

	Size nsymm = sg_.nsymmops();
	Size bestxform=0;
	numeric::xyzVector< core::Real > bestoffset(0,0,0);
	mindis2=1e6;
	com = sg_.c2f()*com;
	for ( Size i=1; i<=nsymm; ++i ) {
		numeric::xyzVector< core::Real > foffset = sg_.symmop(i).get_rotation()*com + sg_.symmop(i).get_translation(), rfoffset;
		rfoffset[0] = min_mod( foffset[0], 1.0 );
		rfoffset[1] = min_mod( foffset[1], 1.0 );
		rfoffset[2] = min_mod( foffset[2], 1.0 );
		Real dist = (sg_.f2c()*rfoffset).length_squared();
		if ( dist<mindis2 ) {
			mindis2=dist;
			bestxform=i;
			bestoffset = foffset - rfoffset;
		}
	}

	numeric::xyzMatrix<Real> R = sg_.f2c()*sg_.symmop(bestxform).get_rotation()*sg_.c2f();
	numeric::xyzVector<Real> T = sg_.f2c()*(sg_.symmop(bestxform).get_translation() - bestoffset);
	pose.apply_transform_Rx_plus_v( R,T );

	return rootpos;
}

void
MakeLatticeMover::detect_connecting_subunits(
	Pose const & monomer_pose,
	Pose const & pose,
	utility::vector1<Size> & monomer_anchors,
	Size &basesubunit
) {
	utility::vector1<Size> new_monomer_anchors;
	Size new_basesubunit=0;

	utility::vector1< numeric::xyzMatrix<core::Real> > new_allRs;
	utility::vector1< numeric::xyzVector<core::Real> > new_allTs;

	// get pose radius
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.size () );

	// as interaction centers, use:
	// a) all protein CB's
	// b) all ligand heavyatoms
	core::Size nIntCtrs=0;
	for ( Size i=1; i<= nres_monomer; ++i ) {
		if ( monomer_pose.residue(i).is_protein() ) {
			nIntCtrs++;
		} else {
			nIntCtrs += monomer_pose.residue(i).nheavyatoms();
		}
	}

	utility::vector1<numeric::xyzVector< core::Real >> monomer_cas(nIntCtrs);

	numeric::xyzVector< core::Real > T0 = pose.residue(monomer_anchors[basesubunit]).xyz("ORIG");
	runtime_assert( T0.length() < 1e-6);

	core::Size ctr=1;
	for ( Size i=1; i<= nres_monomer; ++i ) {
		core::conformation::Residue const &rsd = monomer_pose.residue(i);
		if ( rsd.is_protein() ) {
			if ( rsd.aa() == core::chemical::aa_gly ) {
				monomer_cas[ctr++] = rsd.xyz(2);  //CA
			} else {
				monomer_cas[ctr++] = rsd.xyz(5);  //CB
			}
		} else {
			for ( Size j=1; j<= rsd.nheavyatoms(); ++j ) {
				monomer_cas[ctr++] = rsd.xyz(j);
			}
		}
	}

	Real radius = 0;
	for ( Size i=1; i<=nIntCtrs; ++i ) {
		radius = std::max( monomer_cas[i].length_squared() , radius );
	}
	radius = sqrt(radius);

	// make master first
	new_monomer_anchors.push_back(monomer_anchors[basesubunit]);
	new_allRs.push_back( allRs_[basesubunit] );
	new_allTs.push_back( allTs_[basesubunit] );
	new_basesubunit = 1;

	for ( Size i=1; i<=num_monomers; ++i ) {
		if ( i==basesubunit ) continue;

		// pass 1 check vrt-vrt dist to throw out very distant things
		numeric::xyzVector< core::Real > T = pose.residue(monomer_anchors[i]).xyz("ORIG");
		Real disVRT = T.length();
		if ( disVRT>contact_dist_+2*radius ) continue;

		// pass 2 check ca-ca dists
		numeric::xyzVector< core::Real > X = pose.residue(monomer_anchors[i]).xyz("X") - T;
		numeric::xyzVector< core::Real > Y = pose.residue(monomer_anchors[i]).xyz("Y") - T;
		numeric::xyzVector< core::Real > Z = X.cross( Y );

		numeric::xyzMatrix<core::Real> R = numeric::xyzMatrix<core::Real>::cols(X,Y,Z);

		bool contact=false;
		for ( Size j=1; j<=nIntCtrs && !contact; ++j ) {
			numeric::xyzVector< core::Real > Ri = R*monomer_cas[j] + T;
			for ( Size k=1; k<= nres_monomer && !contact; ++k ) {
				contact = ((Ri-monomer_cas[k]).length_squared() < contact_dist_*contact_dist_);
			}
		}
		if ( contact ) {
			new_monomer_anchors.push_back(monomer_anchors[i]);
			new_allRs.push_back( allRs_[i] );
			new_allTs.push_back( allTs_[i] );
		}
	}
	basesubunit = new_basesubunit;
	monomer_anchors = new_monomer_anchors;
	allRs_ = new_allRs;
	allTs_ = new_allTs;
}


void
MakeLatticeMover::add_monomers_to_lattice(
	Pose const & monomer_pose,
	Pose & pose,
	utility::vector1<Size> const & monomer_anchors,
	utility::vector1<Size> & monomer_jumps,
	Size rootpos
) {
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.size () );

	monomer_jumps.clear();
	Size n_framework_jumps = pose.fold_tree().num_jump();

	Size nres_protein(0);
	for ( Size i=1; i<= num_monomers; ++i ) {
		//std::cerr << "inserting " << i << " of " << num_monomers << std::endl;
		Size const anchor( monomer_anchors[i] + nres_protein ); // since we've already done some insertions
		Size const old_nres_protein( nres_protein );
		pose.insert_residue_by_jump( monomer_pose.residue(rootpos), old_nres_protein+1, anchor ); ++nres_protein;
		for ( Size j=rootpos-1; j>=1; --j ) {
			pose.prepend_polymer_residue_before_seqpos( monomer_pose.residue(j), old_nres_protein+1, false ); ++nres_protein;
		}
		for ( Size j=rootpos+1; j<= monomer_pose.size(); ++j ) {
			if ( monomer_pose.residue(j).is_lower_terminus() || !monomer_pose.residue(j).is_polymer() ) {
				pose.insert_residue_by_jump( monomer_pose.residue(j), nres_protein+1, nres_protein ); ++nres_protein;
			} else {
				pose.append_polymer_residue_after_seqpos( monomer_pose.residue(j), nres_protein, false ); ++nres_protein;
			}
		}
		monomer_jumps.push_back( n_framework_jumps+i );
	}
	for ( Size i=1; i<= num_monomers; ++i ) {
		pose.conformation().insert_chain_ending( nres_monomer*i );
	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	f.reorder( nres_protein+1 );
	pose.fold_tree(f);
}


void
MakeLatticeMover::build_lattice_of_virtuals(
	core::pose::Pose & posebase,
	numeric::xyzVector<int> EXTEND,
	utility::vector1<Size> &Ajumps,
	utility::vector1<Size> &Bjumps,
	utility::vector1<Size> &Cjumps,
	utility::vector1<Size> &subunit_anchors,
	Size &basesubunit
) {
	numeric::xyzVector<int> nvrts_dim(2*EXTEND[0]+1, 2*EXTEND[1]+1, 2*EXTEND[2]+1);

	posebase.clear();
	allRs_.clear();
	allTs_.clear();

	// get gridspace in each dimension
	numeric::xyzVector< Size > grid = sg_.get_nsubdivisions();
	numeric::xyzVector< Size > trans_dofs = sg_.get_trans_dofs();

	numeric::xyzVector< core::Real > O, Ax, Bx, Cx, Ay, By, Cy;
	if ( sg_.setting() == HEXAGONAL ) {
		Ax = sg_.f2c()*numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
		Bx = sg_.f2c()*numeric::xyzVector< core::Real >(0,1,0); Bx.normalize();
		Cx = sg_.f2c()*numeric::xyzVector< core::Real >(0,0,1); Cx.normalize();
	} else {
		Ax = numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
		Bx = numeric::xyzVector< core::Real >(0,1,0); Bx.normalize();
		Cx = numeric::xyzVector< core::Real >(0,0,1); Cx.normalize();
	}

	// we don't care where y is as long as it is perpendicular
	Ay = Bx - Ax.dot(Bx)*Ax; Ay.normalize();
	By = Cx - Bx.dot(Cx)*Bx; By.normalize();
	Cy = Ax - Cx.dot(Ax)*Cx; Cy.normalize();

	ObjexxFCL::FArray3D<int> vrtX, vrtY, vrtZ;
	vrtX.dimension(nvrts_dim[0]*grid[0]+1,1,1); vrtX=0;
	vrtY.dimension(nvrts_dim[0]*grid[0]+1,nvrts_dim[1]*grid[1]+1,1); vrtY=0;
	vrtZ.dimension(nvrts_dim[0]*grid[0]+1,nvrts_dim[1]*grid[1]+1,nvrts_dim[2]*grid[2]+1); vrtZ=0;

	// 1 expand A (j==k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		numeric::xyzVector< core::Real > fX( (Real)(i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], -EXTEND[1], -EXTEND[2] );
		O = sg_.f2c()*fX;

		// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
		ResidueOP vrt_x = make_vrt(O,Ax,Ay);
		ResidueOP vrt_y = make_vrt(O,Bx,By);
		ResidueOP vrt_z = make_vrt(O,Cx,Cy);

		if ( i==1 ) {
			posebase.append_residue_by_bond( *vrt_x );    vrtX(1,1,1) = 1;
			posebase.append_residue_by_jump( *vrt_y, 1);  vrtY(1,1,1) = 2;
			posebase.append_residue_by_jump( *vrt_z, 2);  vrtZ(1,1,1) = 3;
		} else {
			posebase.append_residue_by_jump( *vrt_x, vrtX(i-1,1,1));  vrtX(i,1,1) = posebase.size();
			Ajumps.push_back(posebase.fold_tree().num_jump());
			posebase.append_residue_by_jump( *vrt_y, vrtX(i  ,1,1));  vrtY(i,1,1) = posebase.size();
			posebase.append_residue_by_jump( *vrt_z, vrtY(i  ,1,1));  vrtZ(i,1,1) = posebase.size();
		}
	}

	// expand B (k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		for ( int j=2; j<=(int)(nvrts_dim[1]*grid[1]+1); ++j ) {
			numeric::xyzVector< core::Real > fX( (i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], (j-1-(Real)(EXTEND[1]*grid[1]))/(Real)grid[1], -EXTEND[2] );
			O = sg_.f2c()*fX;

			// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
			ResidueOP vrt_y = make_vrt(O,Bx,By);
			ResidueOP vrt_z = make_vrt(O,Cx,Cy);

			posebase.append_residue_by_jump( *vrt_y, vrtY(i,j-1,1)); vrtY(i,j,1) = posebase.size();
			Bjumps.push_back(posebase.fold_tree().num_jump());
			posebase.append_residue_by_jump( *vrt_z, vrtY(i,j,1));   vrtZ(i,j,1) = posebase.size();
		}
	}

	// expand C
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		for ( int j=1; j<=(int)(nvrts_dim[1]*grid[1]+1); ++j ) {
			for ( int k=2; k<=(int)(nvrts_dim[2]*grid[2]+1); ++k ) {
				numeric::xyzVector< core::Real > fX( (i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], (j-1-(Real)(EXTEND[1]*grid[1]))/(Real)grid[1], (k-1-(Real)(EXTEND[2]*grid[2]))/(Real)grid[2] );
				O = sg_.f2c()*fX;

				// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
				ResidueOP vrt_z = make_vrt(O,Cx,Cy);

				posebase.append_residue_by_jump( *vrt_z, vrtZ(i,j,k-1));  vrtZ(i,j,k) = posebase.size();
				Cjumps.push_back(posebase.fold_tree().num_jump());
			}
		}
	}

	// add "hanging" virtuals
	for ( int s=1; s<=(int)sg_.nsymmops(); ++s ) {
		numeric::xyzMatrix<Real> R_i = sg_.symmop(s).get_rotation(), R_i_cart;
		numeric::xyzVector<Real> T_i = sg_.symmop(s).get_translation();

		// T_i -> indices
		for ( int i=-(int)EXTEND[0]; i<=(int)EXTEND[0]; ++i ) {
			for ( int j=-(int)EXTEND[1]; j<=(int)EXTEND[1]; ++j ) {
				for ( int k=-(int)EXTEND[2]; k<=(int)EXTEND[2]; ++k ) {
					// find lattice anchor
					auto x_i = (int)std::floor( (i+EXTEND[0]+T_i[0])*grid[0] + 1.5 );
					auto y_i = (int)std::floor( (j+EXTEND[1]+T_i[1])*grid[1] + 1.5 );
					auto z_i = (int)std::floor( (k+EXTEND[2]+T_i[2])*grid[2] + 1.5 );

					O = posebase.residue(vrtZ(x_i,y_i,z_i)).xyz("ORIG");

					if ( sg_.setting() == HEXAGONAL ) {
						R_i_cart = sg_.f2c()*R_i*sg_.c2f();
						Ax = R_i_cart*numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
						Ay = R_i_cart*numeric::xyzVector< core::Real >(0,1,0); Ay.normalize();
					} else {
						R_i_cart = R_i;
						Ax = R_i*numeric::xyzVector< core::Real >(1,0,0); Ax.normalize();
						Ay = R_i*numeric::xyzVector< core::Real >(0,1,0); Ay.normalize();
					}
					ResidueOP vrt_z = make_vrt(O,Ax,Ay, (R_i_cart.det()<0));

					posebase.append_residue_by_jump( *vrt_z, vrtZ(x_i,y_i,z_i));

					allRs_.push_back(R_i);
					allTs_.push_back(numeric::xyzVector<Real>(i+T_i[0],j+T_i[1],k+T_i[2]));

					subunit_anchors.push_back(posebase.size());
					if ( s==1 && i==0 && j==0 && k==0 ) {
						basesubunit = subunit_anchors.size();
					}
				}
			}
		}
	}
}


void
MakeLatticeMover::setup_xtal_symminfo(
	Pose & pose,
	Size const num_monomers,
	Size const num_virtuals,
	Size const base_monomer,
	Size const nres_monomer,
	utility::vector1<Size> const &Ajumps,
	utility::vector1<Size> const &Bjumps,
	utility::vector1<Size> const &Cjumps,
	utility::vector1<Size> const &monomer_jumps,
	conformation::symmetry::SymmetryInfo & symminfo
) {

	// bb clones
	for ( Size i=1; i<= num_monomers; ++i ) {
		if ( i != base_monomer ) {
			Size const offset( (i-1)*nres_monomer ), base_offset( (base_monomer-1)*nres_monomer);
			for ( Size j=1; j<= nres_monomer; ++j ) {
				symminfo.add_bb_clone ( base_offset+j, offset+j );
				symminfo.add_chi_clone( base_offset+j, offset+j );
			}
		}
	}

	// subunit base jump clones
	Size const base_monomer_jump( monomer_jumps[ base_monomer ] );
	for ( Size i=1; i<=monomer_jumps.size(); ++i ) {
		if ( monomer_jumps[i]!= base_monomer_jump ) {
			symminfo.add_jump_clone( base_monomer_jump, monomer_jumps[i], 0.0 );
		}
	}

	// unit cell clones
	using core::conformation::symmetry::SymDof;
	std::map< Size, SymDof > symdofs;

	Size Amaster=Ajumps[1], Bmaster=Bjumps[1], Cmaster=Cjumps[1];
	numeric::xyzVector<core::Size> linked_dofs=sg_.get_trans_dofs();
	if ( linked_dofs[1]==1 ) { Bmaster=Ajumps[1]; /*std::cerr << "clone B->A" << std::endl;*/ }
	if ( linked_dofs[2]==1 ) { Cmaster=Ajumps[1]; /*std::cerr << "clone C->A" << std::endl;*/ }

	for ( Size i=2; i<=Ajumps.size(); ++i ) {
		symminfo.add_jump_clone( Amaster, Ajumps[i], 0.0 );
	}
	for ( Size i=1; i<=Bjumps.size(); ++i ) {
		if ( Bmaster!=Bjumps[i] ) {
			symminfo.add_jump_clone( Bmaster, Bjumps[i], 0.0 );
		}
	}
	for ( Size i=1; i<=Cjumps.size(); ++i ) {
		if ( Cmaster!=Cjumps[i] ) {
			symminfo.add_jump_clone( Cmaster, Cjumps[i], 0.0 );
		}
	}

	SymDof symdof_a;
	SymDof symdof_b;
	SymDof symdof_c;

	if ( refinable_lattice_ ) {
		core::Size nrot_dofs=sg_.get_nrot_dofs();
		if ( nrot_dofs == 3 ) {
			symdof_c.read( "x y z" );
			symdof_b.read( "x z" );
		} else if ( nrot_dofs == 1 ) {
			symdof_c.read( "x y" );
			symdof_b.read( "x" );
		} else {
			symdof_c.read( "x" );
			symdof_b.read( "x" );
		}
		symdof_a.read( "x" );

		symdofs[ Ajumps[1] ] = symdof_a;
		if ( Bmaster==Bjumps[1] ) {
			symdofs[ Bjumps[1] ] = symdof_b;
		}
		if ( Cmaster==Cjumps[1] ) {
			symdofs[ Cjumps[1] ] = symdof_c;
		}
	}

	// jump names
	TR << "Initializing " << pose.num_jump() << " jumps." << std::endl;
	for ( Size v=1; v<=pose.num_jump(); ++v ) symminfo.set_jump_name(v, "v_"+utility::to_string(v));

	symminfo.set_jump_name(Ajumps[1], "A");
	symminfo.set_jump_name(Bjumps[1], "B");
	symminfo.set_jump_name(Cjumps[1], "C");

	SymDof symdof_m;
	symdof_m.read( sg_.get_moveable_dofs()+" angle_x angle_y angle_z");
	symdofs[ base_monomer_jump ] = symdof_m;
	symminfo.set_jump_name(base_monomer_jump, "SUB");

	symminfo.set_dofs( symdofs );

	symminfo.num_virtuals( num_virtuals );
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( pose.size(), 0 );
	symminfo.set_nres_subunit( nres_monomer );

	Size const nres_protein( num_monomers * nres_monomer );
	for ( Size i=1; i<= nres_protein; ++i ) {
		if ( symminfo.bb_is_independent( i ) ) symminfo.set_score_multiply( i, 2 );
		else symminfo.set_score_multiply( i, 1 );
	}
	symminfo.update_score_multiply_factor();
}

// parse_my_tag
void
MakeLatticeMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	if ( tag->hasOption( "contact_dist" ) ) {
		contact_dist_ = tag->getOption<core::Real>( "contact_dist" );
	}
}

std::string MakeLatticeMover::get_name() const {
	return mover_name();
}

std::string MakeLatticeMover::mover_name() {
	return "MakeLatticeMover";
}

void MakeLatticeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "contact_dist", xsct_real, "Command line '-interaction_shell ##' OR XML 'contact_dist=##' specifies the distance (in Angstrom) away from the input structure to generate symmetric partners." ) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "If you wish to model a structure in its crystal lattice, a symmetry definition file is not needed. Rather, one can use the flag -symmetry_definition CRYST1.", attlist );
}

std::string MakeLatticeMoverCreator::keyname() const {
	return MakeLatticeMover::mover_name();
}

protocols::moves::MoverOP
MakeLatticeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeLatticeMover );
}

void MakeLatticeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeLatticeMover::provide_xml_schema( xsd );
}


}
}
