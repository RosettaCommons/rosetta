// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/znhash/ZnHash.cc
/// @brief  Implementation of zinc-match hash for use in optimizing zinc coordination
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

// Unit headers
#include <devel/znhash/SymmZnMoversAndTaskOps.hh>
#include <devel/znhash/SymmZnMoversAndTaskOpsCreators.hh>

// Package headers
#include <devel/znhash/ZnHash.hh>

// Protocols headers
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

// Core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <numeric/polynomial.hh>
#include <core/scoring/hbonds/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>

// Protocols headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesMovers.hh>

// Utility headers
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>

// DUMP INCLUDES and clean up later
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/graph/Graph.hh>

#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/pack/make_symmetric_task.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>

#include <protocols/symmetric_docking/SymDockProtocol.hh>
#include <protocols/symmetric_docking/SymDockingHiRes.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RampingMover.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesMovers.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>

#include <utility/string_util.hh>

// option keys
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <sstream>


namespace devel {
namespace znhash {

/// @details Auto-generated virtual destructor
FindZnCoordinatingResidues::~FindZnCoordinatingResidues() {}

static thread_local basic::Tracer TR( "devel.znhash.SymmZnMoversAndTaskOps" );

InitializeZNCoordinationConstraintMover::InitializeZNCoordinationConstraintMover() :
	parent( "InitializeZNCoordinationConstraintMover" ),
	znreach_( 3.0 ),
	orbital_dist_( 1.0 ),
	orbital_reach_( 1.0 ),
	znwelldepth_( 3.0 ),
	clash_weight_( 0.0 ),
	require_3H_( false ),
	idealize_input_virtual_atoms_( true )
{}

void InitializeZNCoordinationConstraintMover::set_zn_reach( core::Real reach ) { znreach_ = reach; }
void InitializeZNCoordinationConstraintMover::set_orbital_dist( core::Real dist ) { orbital_dist_ = dist; }
void InitializeZNCoordinationConstraintMover::set_orbital_reach( core::Real reach ) { orbital_reach_ = reach; }
void InitializeZNCoordinationConstraintMover::set_zn_well_depth( core::Real depth ) { znwelldepth_ = depth; }
void InitializeZNCoordinationConstraintMover::set_clash_weight( core::Real weight ) { clash_weight_ = weight; }
void InitializeZNCoordinationConstraintMover::require_3H( bool setting ) { require_3H_ = setting; }
void InitializeZNCoordinationConstraintMover::set_idealize_input_virtual_atoms( bool setting ) { idealize_input_virtual_atoms_ = setting; }
void InitializeZNCoordinationConstraintMover::set_matcher_constraint_file_name( std::string const & fname ) { matcher_constraint_file_name_ = fname; }
void InitializeZNCoordinationConstraintMover::set_reference_pdb( std::string const & fname ) { reference_pdb_ = fname; }
void InitializeZNCoordinationConstraintMover::set_match_pdb_listfilename( std::string const & fname ) { match_pdb_listfilename_ = fname; }

InitializeZNCoordinationConstraintMover::MoverOP
InitializeZNCoordinationConstraintMover::clone() const {
	return new InitializeZNCoordinationConstraintMover( *this );
}

std::string
InitializeZNCoordinationConstraintMover::get_name() const
{
	return "InitializeZNCoordinationConstraintMover";
}

void InitializeZNCoordinationConstraintMover::apply( core::pose::Pose & p )
{
	recover_sidechains_ = new protocols::simple_moves::ReturnSidechainMover( p );
	if ( zn_score_ == 0 ) {
		zn_score_ = new devel::znhash::ZnCoordinationScorer;

		zn_score_->set_zn_reach( znreach_ );
		zn_score_->set_orbital_dist( orbital_dist_ );
		zn_score_->set_orbital_reach( orbital_reach_ );
		zn_score_->set_zn_well_depth( znwelldepth_ );
		zn_score_->set_clash_weight( clash_weight_ );
		zn_score_->require_3H( require_3H_ );
		zn_score_->set_idealize_input_virtual_atoms( idealize_input_virtual_atoms_ );
		zn_score_->set_matcher_constraint_file_name( matcher_constraint_file_name_ );
		zn_score_->set_reference_pdb( reference_pdb_ );

		// initialize the assym_resid and the third_resid
		if (p.total_residue() < 3 ) {
			utility_exit_with_message("Cannot use ZnCoordinationConstraint on a pose with fewer than three residues" );
		}
		zn_score_->set_asymm_resid(1);
		if ( p.residue(1).chain() == p.residue(2).chain() ) {
			zn_score_->set_third_resid(2);
		} else {
			utility_exit_with_message("Expected chain 1 to have at least 2 residues");
		}

		// look for the first residue on chain b -- it should not be either residue 1 or residue 2
		for ( core::Size ii = 3; ii <= p.total_residue(); ++ii ) {
			if ( p.residue(ii).chain() == 2 ) {
				zn_score_->set_symm_resid( ii );
				break;
			}
			if ( ii == p.total_residue()) {
				utility_exit_with_message("Did not find chain B in input structure" );
			}
		}

		std::ifstream listfile( match_pdb_listfilename_.c_str() );
		core::Size linenum(0), count_matches_added(0);
		while ( listfile.good() ) {
			std::string line;
			getline( listfile, line );
			++linenum;
			if ( line.size() == 0 ) { continue; }

			std::istringstream linestream( line );
			if ( ! linestream.good() ) {
				utility_exit_with_message("!linestream.good()\nWhile reading line " + utility::to_string(linenum) + " of " + match_pdb_listfilename_ + " could not read matchfilename\n" + line );
			}
			std::string matchfilename;
			if ( linestream.peek() == '#' ) continue;
			if ( linestream.peek() == '\n') continue;
			linestream >> matchfilename;
			if ( matchfilename == "" ) {
				utility_exit_with_message("While reading line " + utility::to_string(linenum) + " of " + match_pdb_listfilename_ + " could not read matchfilename\n" + line );
			}
			zn_score_->add_match_from_file( matchfilename );
			++count_matches_added;
		}
		zn_score_->finalize_after_all_matches_added();
		//TR << "Added " << count_matches_added << " Zn matches to the znscore" << std::endl;
	}

	ZnCoordinationConstraintOP zncst = new ZnCoordinationConstraint( zn_score_ );
	p.add_constraint( zncst );

}

ZnCoordinationScorerOP
InitializeZNCoordinationConstraintMover::zn_score() const { return zn_score_; }

protocols::simple_moves::ReturnSidechainMoverOP InitializeZNCoordinationConstraintMover::recover_sidechains() const { return recover_sidechains_; }

/////////////////////////////////////////////////////////////////////

ZNCoordinationConstraintReporterMover::ZNCoordinationConstraintReporterMover(
	InitializeZNCoordinationConstraintMoverOP init_zn
) : parent( "ZNCoordinationConstraintReporterMover" ), init_zn_( init_zn ) {}

ZNCoordinationConstraintReporterMover::MoverOP
ZNCoordinationConstraintReporterMover::clone() const {
	return new ZNCoordinationConstraintReporterMover( init_zn_ );
}

std::string
ZNCoordinationConstraintReporterMover::get_name() const
{
	return "ZNCoordinationConstraintReporterMover";
}

void ZNCoordinationConstraintReporterMover::apply( core::pose::Pose & p )
{
	// send to the log file info on which residues have been identified by the zncoordcst as best for coordinating zinc
	devel::znhash::ZnCoordinationScorerOP zn_score = init_zn_->zn_score();
	std::pair< core::Size, core::Size > result = zn_score->best_match( p );
	if ( result.first == 0 || result.second == 0 ) {
		//TR << "No zn constraint satisfied" << std::endl;
		return;
	}

	//devel::znhash::ZnMatchData const & m1 = zn_score->zn_matches()[ result.first ];
	//devel::znhash::ZnMatchData const & m2 = zn_score->zn_matches()[ result.second ];

	//core::Size symm_res_resid( zn_score->focused_clone_atids()[1].rsd() );

	//using namespace core::conformation::symmetry;

	//SymmetricConformation const & SymmConf (
	//	dynamic_cast<SymmetricConformation const &> ( p.conformation()) );
	//SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	//Size nsubunits = symm_info->subunits();
	//Size nres_asu = symm_info->num_independent_residues();
	//Size nres_monomer = symm_info->num_independent_residues() - 1;


	TR << "Best zn coordination: score = " << zn_score->score( p ) << std::endl;
	//TR << "Coordinated by residues " << m1.res1() << " and " << m1.res2() << " on chain A" << std::endl;
	//TR << "Coordinated by residues " << m2.res1() << " and " << m2.res2() << " on chain B" << std::endl;
	//TR << "Pymol: hide lines, elem H; show cartoon; show sticks, res " << m1.res1() << "+" << m1.res2()
	//	<< " and chain A; show sticks, res " << nres_asu + m2.res1() << "+" << nres_asu + m2.res2() << " and chain B; show sticks, resn ZNX" << std::endl;
}


/////////////////////////////////////////////////////////////////////

matchfile_header
read_match_header_line_from_pdb(
	std::string const & fname,
	core::Size linenum,
	std::istringstream & matchline_stream
)
{
	matchfile_header header;

	matchline_stream >> header.remark_string_;
	if ( header.remark_string_ != "REMARK" || ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read REMARK in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.remark_number_;
	if ( header.remark_number_ != "666"  || ! matchline_stream.good()) {
		utility_exit_with_message( "Expected to read 666 in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.match_string_;
	if ( header.match_string_ != "MATCH" || ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read MATCH in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.template_string_;
	if ( header.template_string_ != "TEMPLATE" || ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read TEMPLATE in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.lig_chain_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a chain for the ligand residue in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.lig_resname_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read ligand resname in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.lig_resnum_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a residue number for the ligand residue in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.matchstr2_;
	if ( header.matchstr2_ != "MATCH" || ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read MATCH in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.motifstr_;
	if ( header.motifstr_ != "MOTIF" || ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read MOTIF in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.prot_chain_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a chain for the protein residue in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.aa_str_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read an amino acid 3-letter code in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.aares_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a residue index for the protein residue in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.geocst_indstr_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a geometric constraint index (1) in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	matchline_stream >> header.geocst_indstr2_;
	if ( ! matchline_stream.good() ) {
		utility_exit_with_message( "Expected to read a geometric constraint index (2) in file " + fname + " on line " + utility::to_string( linenum ) );
	}
	return header;
}

/////////////////////////////////////////////////////////////////////
void
add_znx_coordination_remark_lines_to_pose(
	core::pose::Pose & p,
	symdes_znx_coordination_data const & coordination_data
)
{
	// REMARK 666 MATCH TEMPLATE X ZNX    0 MATCH MOTIF A HIS  140  1  1

	core::pose::Remarks match_remarks;
	for ( core::Size ii = 1; ii <= 6; ++ii ) {
		core::Size znind = ( ii + 1 ) / 2;

		core::pose::RemarkInfo rinfo;
		rinfo.num = 666;
		std::ostringstream ss;
		ss << "MATCH TEMPLATE "
			<< coordination_data.zinc_data_[ znind ].chain_
			<< " ZNX "
			<< coordination_data.zinc_data_[ znind ].resind_
			<< " MATCH MOTIF "
			<< coordination_data.protein_data_[ ii ].chain_
			<< " "
			<< coordination_data.protein_data_[ ii ].name3_
			<< " "
			<< coordination_data.protein_data_[ ii ].resindex_
			<< " " << ii << " "
			<< coordination_data.protein_data_[ ii ].exgeom_index_ << "\n";
		rinfo.value = ss.str();
		TR << "remark " << ii << ": " << rinfo.value;
		match_remarks.push_back( rinfo );
	}
	TR << std::endl;
	core::pose::PDBInfoOP info = p.pdb_info();
	if ( ! info ) { info = new core::pose::PDBInfo; }
	info->remarks() = match_remarks;
	p.pdb_info( info );

}

/////////////////////////////////////////////////////////////////////
ZNCoordinationConstraintPlacerMover::ZNCoordinationConstraintPlacerMover(
	InitializeZNCoordinationConstraintMoverOP init_zn
) :
	parent( "ZNCoordinationConstraintPlacerMover" ),
	init_zn_( init_zn ),
	constraint_energy_cutoff_( 25.0 )
{
}

void ZNCoordinationConstraintPlacerMover::set_constraint_energy_cutoff(
	core::Real setting
)
{
	constraint_energy_cutoff_ = setting;
}

void ZNCoordinationConstraintPlacerMover::set_four_residue_cst_fname(
	std::string const & fname
)
{
	cst_fname_ = fname;
}

ZNCoordinationConstraintPlacerMover::MoverOP
ZNCoordinationConstraintPlacerMover::clone() const
{
	return new ZNCoordinationConstraintPlacerMover( init_zn_ );
}

std::string
ZNCoordinationConstraintPlacerMover::get_name() const {
	return "ZNCoordinationConstraintPlacerMover";
}

void ZNCoordinationConstraintPlacerMover::apply( core::pose::Pose & p )
{

	// remove the existing ZnHash constraint from the Pose.
	p.constraint_set( core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ) );

	// Record which residues have been identified by the zncoordcst as best for coordinating zinc
	// devel::znhash::ZnCoordinationScorerOP zn_score = init_zn_->zn_score();
	devel::znhash::ZnCoordinationScorerOP zn_score = init_zn_->zn_score();
	std::pair< core::Size, core::Size > result = zn_score->best_match( p );
	if ( result.first == 0 || result.second == 0 ) {
		//TR << "No zn constraint satisfied" << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );

		return;
	}
	m1_ = zn_score->zn_matches()[ result.first ];
	m2_ = zn_score->zn_matches()[ result.second ];

	/// restore the native sidechains and now
	core::pose::Pose recovered_sidechain_pose;
	init_zn_->recover_sidechains()->apply( p );
	recovered_sidechain_pose = p;

	mutate_the_interface_to_alanine( p );

	//TR << "BEFORE ZINC RESIDUE ADDITION" << std::endl;
	//TR << "fold tree: " << std::endl;
	//TR << p.fold_tree() << std::endl;

	//core::conformation::symmetry::SymmetricConformation const & symm_conf (
	//	dynamic_cast< core::conformation::symmetry::SymmetricConformation const & > ( p.conformation()) );
	//core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	//TR << "SymmInfo" << std::endl;
	//TR << *symm_info << std::endl;


	zn_score->insert_match_onto_pose( p, result.first, 1 );
	zn_score->insert_match_onto_pose( p, result.second, 2 );

	insert_zn_residues_into_pose( p );
	add_matcher_remark_lines_for_zn_coordination( p );
	minimize_zinc_coordination( p );
	quick_hires_symdock( p );
	restore_alanine_interface_residues_to_wtconf( recovered_sidechain_pose, p );

	filter_by_constraint_score( p );
}

void ZNCoordinationConstraintPlacerMover::mutate_the_interface_to_alanine( core::pose::Pose & p )
{
	// mutate the interface to alanine.
	core::scoring::ScoreFunctionOP fa_sfxn = core::scoring::get_score_function();
	(*fa_sfxn)( p );

	protocols::simple_moves::symmetry::SymPackRotamersMover sympack;
	mutate_interface_residues_to_alanine_task_ = core::pack::task::TaskFactory::create_packer_task( p );
	core::pack::make_symmetric_PackerTask_by_truncation( p, mutate_interface_residues_to_alanine_task_ );
	protocols::toolbox::task_operations::RestrictToInterface rti_taskop;
	for ( Size ii = 1; ii <= p.fold_tree().num_jump(); ++ii ) { rti_taskop.add_jump( ii ); }
	rti_taskop.distance( 12 );
	rti_taskop.apply( p, *mutate_interface_residues_to_alanine_task_ );
	core::pack::task::operation::RestrictAbsentCanonicalAAS make_ala;
	make_ala.keep_aas( "A" ); // keep only alanine.
	make_ala.include_residue(0); // 0 means apply this to all residues
	make_ala.apply( p, *mutate_interface_residues_to_alanine_task_ );
	core::pack::task::operation::NoRepackDisulfides nodisulf;
	nodisulf.apply( p, *mutate_interface_residues_to_alanine_task_ );

	//std::cout << "packer task: " << *task << std::endl;
	sympack.task( mutate_interface_residues_to_alanine_task_ );
	sympack.score_function( fa_sfxn );
	sympack.apply( p );

}

void ZNCoordinationConstraintPlacerMover::insert_zn_residues_into_pose( core::pose::Pose & p ) {
	// Now, let's create a zinc residue and attach it by a jump to m1.res1 on all four chains?
	devel::znhash::ZnCoordinationScorerOP zn_score = init_zn_->zn_score();

	assert( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->has_name( "ZNX" ));

	core::chemical::ResidueType const & znx_restype =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->
		name_map( "ZNX" );

	HTReal match_frame(  m1_.res1conf().xyz(1), m1_.res1conf().xyz(2), m1_.res1conf().xyz(3) );

	//utility::vector1< Size > new_zn_resids(4, 0);
	// after Frank's fix, should only need one append_residue_by_jump call.
	for ( core::Size ii = 1; ii <= 1; ++ii ) {

		core::Size m1resid = (ii-1) * ( zn_score->r2() - zn_score->r1() ) + m1_.res1();
		core::conformation::ResidueOP znres = core::conformation::ResidueFactory::create_residue( znx_restype );
		core::conformation::Residue const & res1chA = p.residue( m1resid );
		HTReal dock_frame( res1chA.xyz(1), res1chA.xyz(2), res1chA.xyz(3) );
		HTReal transform_match_coords_to_docked_coords = dock_frame * match_frame.inverse();

		for ( core::Size ii = 1; ii <= 5; ++ii ) {
			znres->set_xyz( ii, transform_match_coords_to_docked_coords * m1_.znconf().xyz( ii ) );
		}
		//new_zn_resids.push_back( p.total_residue() + 1 );

		p.append_residue_by_jump( *znres, m1resid, "CA", "ZN", true );

		//core::conformation::Conformation c;
		//c.append_residue_by_jump( *znres, 1, "", "ZN", true );
		//p.conformation().insert_conformation_by_jump( c, 1, 1, m1resid, 0, "CA", "ZN" );

	}

	//core::pose::symmetry::make_symmetric_pdb_info()
	core::pose::PDBInfoOP newinfo = new core::pose::PDBInfo( p ); // fake new info
	p.pdb_info( newinfo );

	//for ( Size ii = 1; ii <= p.total_residue(); ++ii ) {
	//	std::cout << "Residue " << ii << " chain: " << p.residue(ii).chain() << " ";
	//	std::cout << p.residue(ii).name() << std::endl;
	//}
	//core::kinematics::FoldTree ft( p.fold_tree() );
	//p.fold_tree( ft ); // see if this works?!
}

void
ZNCoordinationConstraintPlacerMover::add_matcher_remark_lines_for_zn_coordination(
	core::pose::Pose &  p
)
{
	// Let's try and add constraints using the existing enzdes constraint addition code
	// we'll have to create remarks that match the match file and that are tailored for this pose.
	// Here's an example:
	// REMARK 666 MATCH TEMPLATE X ZNX    0 MATCH MOTIF A HIS  140  1  1
	// REMARK 666 MATCH TEMPLATE X ZNX    0 MATCH MOTIF A ASP   93  2  2
	// We'll need to add two sets of constraints: one between the asymmetric unit
	// and the first zinc residue; the second between the assymetric unit and the Nth
	// zinc residue -- the symmetric clone of the first zinc residue that has wrapped
	// around back to

	// Symmetry info
	using namespace core::conformation::symmetry;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( p.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	Size nsubunits = symm_info->subunits();
	Size nres_asu = symm_info->num_independent_residues();
	// Size nres_monomer = symm_info->num_independent_residues() - 1; // Unused variable causes warning.

	symdes_znx_coordination_data coordination_data;
	// Read in the original remark lines from the m1 match constraint file and translate them
	// into new
	std::ifstream m1matchfile(m1_.match_pdb_file().c_str());
	core::pose::Remarks match_remarks;

	coordination_data.zinc_data_[ 1 ].chain_ = 'B';
	coordination_data.zinc_data_[ 1 ].resind_ = nres_asu;
	for ( Size ii = 1; ii <= 2; ++ii ) {
		std::string match_line;
		std::getline( m1matchfile, match_line );
		std::cout << "original line: " << match_line << std::endl;
		std::istringstream matchline_stream( match_line );

		matchfile_header header = read_match_header_line_from_pdb( m1_.match_pdb_file(), ii, matchline_stream );

		coordination_data.protein_data_[ ii ].chain_ = 'A';
		coordination_data.protein_data_[ ii ].resindex_ = ( ii == 1 ? m1_.res1() : m1_.res2() );
		coordination_data.protein_data_[ ii ].name3_ = ( ii == 1 ?
			p.residue( m1_.res1()).name3() :
			p.residue( m1_.res2()).name3() );
		coordination_data.protein_data_[ ii ].exgeom_index_ = header.geocst_indstr2_;

		// OK: now what?

		//core::pose::RemarkInfo rinfo;
		//rinfo.num = 666;
		//std::ostringstream ss;
		//ss << "MATCH TEMPLATE B ZNX " << nres_asu << " MATCH MOTIF A "
		//	<< ( ii == 1 ?
		//		p.residue( m1_.res1()).name3() :
		//		p.residue( m1_.res2()).name3() )
		//	<< " "
		//	<< ( ii == 1 ? m1_.res1() : m1_.res2() )
		//	<< " " << ii << " " << header.geocst_indstr2_ << "\n";  // I don't think a \n is required.
		//rinfo.value = ss.str();
		//std::cout << "remark " << ii << ": " << rinfo.value << std::endl;
		//match_remarks.push_back( rinfo );
	}


	// OK, now prepare the match remark lines for the constraints to the other zinc that's chemically bound
	// to chain A.

	std::ifstream m2matchfile(m2_.match_pdb_file().c_str());

	/// v3v4: Does residue 1 coordinate v3 and residue 2 coordinate v4?
	/// if not, then residue 1 coordinates v4 and residue 2 coordinates v3.
	bool const v3v4 = ! init_zn_->zn_score()->optimal_coordination_is_reversed( p );

	coordination_data.zinc_data_[ 2 ].chain_ = 'B';
	coordination_data.zinc_data_[ 2 ].resind_ = nres_asu;
	coordination_data.zinc_data_[ 3 ].chain_ = char( int('A') + 2*nsubunits - 1 );
	coordination_data.zinc_data_[ 3 ].resind_ = nres_asu*nsubunits;

	for ( Size ii = 1; ii <= 2; ++ii ) {
		std::string match_line;
		std::getline( m2matchfile, match_line );
		std::cout << "original line: " << match_line << std::endl;
		std::istringstream matchline_stream( match_line );

		Size const which_cst = ( v3v4 == (ii==1) ) ? 3 : 4;

		matchfile_header header = read_match_header_line_from_pdb( m2_.match_pdb_file(), ii, matchline_stream );

		//core::pose::RemarkInfo rinfo;
		//rinfo.num = 666;
		//std::ostringstream ss;
		//ss << "MATCH TEMPLATE " << char( int('A') + 2*nsubunits - 1 )
		//	<< " ZNX " << nsubunits * nres_asu << " MATCH MOTIF A "
		//	<< ( ii == 1 ?
		//		p.residue( m2_.res1()).name3() :
		//		p.residue( m2_.res2()).name3() )
		//	<< " "
		//	<< ( ii == 1 ? m2_.res1() : m2_.res2() )
		//	<< " " << which_cst << " " << header.geocst_indstr2_<< "\n";  // I don't think a \n is required.
		//rinfo.value = ss.str();
		//std::cout << "remark " << which_cst << ": " << rinfo.value << std::endl;
		//match_remarks.push_back( rinfo );

		coordination_data.protein_data_[ which_cst ].chain_ = 'C';
		coordination_data.protein_data_[ which_cst ].resindex_ = nres_asu + ( ii == 1 ? m2_.res1() : m2_.res2() );
		coordination_data.protein_data_[ which_cst ].name3_ = ( ii == 1 ?
			p.residue( m2_.res1()).name3() :
			p.residue( m2_.res2()).name3() );
		coordination_data.protein_data_[ which_cst ].exgeom_index_ = header.geocst_indstr2_;

		coordination_data.protein_data_[ which_cst + 2].chain_ = 'A';
		coordination_data.protein_data_[ which_cst + 2].resindex_ = ( ii == 1 ? m2_.res1() : m2_.res2() );
		coordination_data.protein_data_[ which_cst + 2 ].name3_ = ( ii == 1 ?
			p.residue( m2_.res1()).name3() :
			p.residue( m2_.res2()).name3() );
		coordination_data.protein_data_[ which_cst + 2 ].exgeom_index_ = header.geocst_indstr2_;

	}

	//core::pose::PDBInfoOP info = p.pdb_info();
	//if ( ! info ) { info = new core::pose::PDBInfo; }
	//info->remarks() = match_remarks;
	//p.pdb_info( info );
	add_znx_coordination_remark_lines_to_pose( p, coordination_data );

	protocols::enzdes::AddOrRemoveMatchCsts cst_adder;
	cst_adder.cstfile( cst_fname_ );
	cst_adder.set_cst_action( protocols::enzdes::ADD_NEW );
	cst_adder.apply( p );

	core::scoring::ScoreFunctionOP fa_sfxn = fa_sfxn_w_cstterms();

	//TR << "Score after adding first set of constraints" << std::endl;
	//fa_sfxn->show( TR, p );
	//TR << std::endl;

	//std::string tag = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
	//std::string nstruct_tag = utility::to_string( protocols::jd2::JobDistributor::get_instance()->current_job()->nstruct_index() );
	//p.dump_pdb( tag + "_" + nstruct_tag + "_before_min.pdb" );

}


core::scoring::ScoreFunctionOP
ZNCoordinationConstraintPlacerMover::fa_sfxn_w_cstterms() const {
	core::scoring::ScoreFunctionOP fa_sfxn = core::scoring::get_score_function();
	fa_sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	fa_sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
	fa_sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );
	return fa_sfxn;
}

void
ZNCoordinationConstraintPlacerMover::minimize_zinc_coordination( core::pose::Pose &  p )
{
	core::scoring::ScoreFunctionOP fa_sfxn = fa_sfxn_w_cstterms();
	fa_sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	fa_sfxn->set_weight( core::scoring::angle_constraint, 1.0 );
	fa_sfxn->set_weight( core::scoring::dihedral_constraint, 1.0 );

	//TR << "fold tree: " << std::endl;
	//TR << p.fold_tree() << std::endl;

	core::conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast< core::conformation::symmetry::SymmetricConformation const & > ( p.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	//TR << "SymmInfo" << std::endl;
	//TR << *symm_info << std::endl;


	//TR << "Before minimization score: ";
	//fa_sfxn->show( TR, p  );
	//TR << std::endl;

	core::kinematics::MoveMapOP mm = new core::kinematics::MoveMap;
	mm->set_jump( true );
	mm->set_chi( m1_.res1(), true );
	mm->set_chi( m1_.res2(), true );
	mm->set_chi( m2_.res1(), true );
	mm->set_chi( m2_.res2(), true );
	protocols::simple_moves::symmetry::SymMinMoverOP symminmover =
		new protocols::simple_moves::symmetry::SymMinMover( mm, fa_sfxn, "dfpmin_armijo_nonmonotone", 1e-4, true );

	core::scoring::EnergyMap start_weights = fa_sfxn->weights();
	core::scoring::EnergyMap end_weights = start_weights;

	start_weights[ core::scoring::atom_pair_constraint ] = 200.0;
	start_weights[ core::scoring::angle_constraint     ] = 200.0;
	start_weights[ core::scoring::dihedral_constraint  ] = 200.0;
	start_weights[ core::scoring::fa_rep               ] = 0.04;

	end_weights[ core::scoring::atom_pair_constraint ] = 50.0;
	end_weights[ core::scoring::angle_constraint     ] = 50.0;
	end_weights[ core::scoring::dihedral_constraint  ] = 50.0;

	protocols::moves::RampingMoverOP ramp =
		new protocols::moves::RampingMover(symminmover, fa_sfxn, start_weights, end_weights, 10, 1, 0 );
	ramp->apply( p );

	//TR << "After minimization score: ";
	//fa_sfxn->show( TR, p  );
	//TR << std::endl;

	//std::string tag = protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag();
	//std::string nstruct_tag = utility::to_string( protocols::jd2::JobDistributor::get_instance()->current_job()->nstruct_index() );
	//p.dump_pdb( tag + "_" + nstruct_tag + "_after_min.pdb" );


}

void
ZNCoordinationConstraintPlacerMover::quick_hires_symdock( core::pose::Pose & p ) {
	core::scoring::ScoreFunctionOP fa_sfxn = fa_sfxn_w_cstterms();
	fa_sfxn->set_weight( core::scoring::atom_pair_constraint, 50.0 );
	fa_sfxn->set_weight( core::scoring::angle_constraint, 50.0 );
	fa_sfxn->set_weight( core::scoring::dihedral_constraint, 50.0 );

	protocols::symmetric_docking::SymDockingHiRes docking_high( fa_sfxn, fa_sfxn );
	docking_high.apply( p );
	//TR << "Score after hi-res dicking" << std::endl;
	//fa_sfxn->show( TR, p );
}

void
ZNCoordinationConstraintPlacerMover::restore_alanine_interface_residues_to_wtconf(
	core::pose::Pose const & pose_w_scs,
	core::pose::Pose & p
)
{
		// Symmetry info
	using namespace core::conformation::symmetry;

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( p.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	Size nres_asu = symm_info->num_independent_residues();

	for ( Size ii = 1; ii <= nres_asu - 1; ++ii ) {
		if ( ii == m1_.res1() || ii == m1_.res2() || ii == m2_.res1() || ii == m2_.res2() ) continue;
		if ( ! mutate_interface_residues_to_alanine_task_->being_packed( ii ) ) continue;
		p.replace_residue( ii, pose_w_scs.residue(ii), true );
	}
}


void
ZNCoordinationConstraintPlacerMover::filter_by_constraint_score( core::pose::Pose const & p ) {
	using namespace core::scoring;
	EnergyMap tots = p.energies().total_energies();
	tots *= p.energies().weights();
	core::Real sum = tots[ atom_pair_constraint ] + tots[ dihedral_constraint ] + tots[ angle_constraint ];
	TR << "Final score: " << sum << " from dist: " << tots[ atom_pair_constraint ] << " ang: " << tots[ angle_constraint ] << " dihe: " << tots[ dihedral_constraint ];
	if ( sum < constraint_energy_cutoff_ ) {
		TR << " PASSED FILTER" << std::endl;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << " FAILED FILTER" << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}
}

/////////////////////////////////////////////////////////////////////

FindZnCoordinatingResidues::FindZnCoordinatingResidues() : fail_on_absent_coordinators_( true ) {}

void FindZnCoordinatingResidues::fail_on_absent_coordinators( bool setting )
{
	fail_on_absent_coordinators_ = setting;
}

// @brief Try to find, based solely on the structure of the input pose,
// which residues were meant to be those residues coordinating the two
// zinc residues coordinated by the asymmetric unit.
// Necessary only because the silent file doesn't hold remarks, it seems.
void FindZnCoordinatingResidues::find_coordinating_residues(
	core::pose::Pose const & p
)
{
	using core::chemical::aa_his;
	using core::chemical::aa_asp;
	using core::chemical::aa_glu;

	resinds_.clear(); resinds_.reserve( 4 );
	atomids_.clear(); atomids_.reserve( 4 );

	core::pose::Pose copy_pose( p );
	core::scoring::symmetry::SymmetricScoreFunction sfxn;
	sfxn.set_weight( core::scoring::fa_atr, 1.0 );
	sfxn( copy_pose );

	core::conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast< core::conformation::symmetry::SymmetricConformation const & > ( p.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	core::Size nres_asu = symm_info->num_independent_residues();
	core::Size nsubunits = symm_info->subunits();

	utility::vector1< std::string > his_coordinating_atoms; his_coordinating_atoms.push_back("ND1"); his_coordinating_atoms.push_back("NE2");
	utility::vector1< std::string > asp_coordinating_atoms; asp_coordinating_atoms.push_back("OD1"); asp_coordinating_atoms.push_back("OD2");
	utility::vector1< std::string > glu_coordinating_atoms; glu_coordinating_atoms.push_back("OE1"); glu_coordinating_atoms.push_back("OE2");

	// ok -- look at residues in the vacinity of the 1st znx residue
	for ( Size ii = 1; ii <= 4; ++ii ) {
		Size znxind = ii < 3 ? nres_asu : nsubunits*nres_asu;
		assert( copy_pose.residue(znxind).name() == "ZNX" );


		std::string vname( std::string("V") + std::string(1, char( ii + '1' - 1 )));
		core::conformation::Residue const & znx = copy_pose.residue( znxind );
		core::Vector zncoord = znx.xyz(znx.atom_index("ZN"));
		core::Vector vcoord = znx.xyz(znx.atom_index(vname));
		core::Real best_dist = 0;
		core::Size ind_w_best_dist( 0 ), atomind_w_best_dist( 0 );
		for ( core::graph::Graph::EdgeListConstIter
				eiter    = copy_pose.energies().energy_graph().get_node(znxind)->const_edge_list_begin(),
				eiterend = copy_pose.energies().energy_graph().get_node(znxind)->const_edge_list_end();
				eiter != eiterend; ++eiter ) {
			core::conformation::Residue const & nbr( copy_pose.residue( (*eiter)->get_other_ind( znxind ) ));

			std::pair< core::Real, core::Size > closest_dist_info;
			switch ( nbr.aa() ) {
				case aa_his :
					closest_dist_info = closest_distance_to_desired_vrt( zncoord, vcoord, nbr, his_coordinating_atoms );
				break;
				case aa_asp :
					closest_dist_info = closest_distance_to_desired_vrt( zncoord, vcoord, nbr, asp_coordinating_atoms );
				break;
				case aa_glu :
					closest_dist_info = closest_distance_to_desired_vrt( zncoord, vcoord, nbr, glu_coordinating_atoms );
				break;
			default:
				// Only looking for EDH coordination of the zinc -- this might need to change in the future
				continue;
			}
			if ( closest_dist_info.first < 0 ) continue;
			if ( ind_w_best_dist == 0 ||  closest_dist_info.first < best_dist ) {
				best_dist = closest_dist_info.first;
				ind_w_best_dist = nbr.seqpos();
				atomind_w_best_dist = closest_dist_info.second;
			}
		}
		if ( ind_w_best_dist == 0 && fail_on_absent_coordinators_ ) {
			utility_exit_with_message( "Unable to find a protein residue coordinating ZNX residue, virt atom #" + utility::to_string( ii ) );
		}
		if ( ind_w_best_dist != 0 ) {
			resinds_.push_back( ind_w_best_dist );
			atomids_.push_back( core::id::AtomID( atomind_w_best_dist, ind_w_best_dist ));
			TR << "Determined residue " << ind_w_best_dist << " coordinates Zinc at virtual atom " << ii << std::endl;
		}
	}
}

std::pair< core::Real, core::Size >
FindZnCoordinatingResidues::closest_distance_to_desired_vrt(
	core::Vector const & zncoord,
	core::Vector const & vcoord,
	core::conformation::Residue const & nbr,
	utility::vector1< std::string > const & coord_atnames
) const
{
	core::Real closest_d2(-1); core::Size ind = 0;
	for ( core::Size ii = 1; ii <= coord_atnames.size(); ++ii ) {
		assert( nbr.has( coord_atnames[ii] ));
		core::Size iiind = nbr.atom_index( coord_atnames[ii] );
		core::Vector iipos = nbr.xyz(iiind);
		core::Real d2zn = iipos.distance_squared( zncoord );
		if ( d2zn > 9 ) continue; // cutoff of 3 angstroms between coordinating atom and ZN
		core::Real d2v = iipos.distance_squared( vcoord );
		if ( closest_d2 < 0 || d2v < closest_d2 ) {
			closest_d2 = d2v;
			ind = iiind;
		}
	}
	if ( closest_d2 > 0 ) {
		return std::make_pair( std::sqrt( closest_d2 ), ind );
	} else {
		return std::make_pair( core::Real(-1.0), core::Size(0) );
	}
}
/////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
InsertZincCoordinationRemarkLinesCreator::create_mover() const
{
	return new InsertZincCoordinationRemarkLines;
}

std::string InsertZincCoordinationRemarkLinesCreator::keyname() const
{
	return "InsertZincCoordinationRemarkLines";
}


/////////////////////////////////////////////////////////////////////

InsertZincCoordinationRemarkLines::InsertZincCoordinationRemarkLines() {}
InsertZincCoordinationRemarkLines::~InsertZincCoordinationRemarkLines() {}

protocols::moves::MoverOP
InsertZincCoordinationRemarkLines::clone() const { return new InsertZincCoordinationRemarkLines; }

std::string
InsertZincCoordinationRemarkLines::get_name() const { return "InsertZincCoordinationRemarkLines"; }

void InsertZincCoordinationRemarkLines::apply( core::pose::Pose & p )
{
	using namespace core::conformation::symmetry;

	if ( ! dynamic_cast<SymmetricConformation const * > ( & p.conformation()) ) {
		utility_exit_with_message( "InsertZincCoordinationRemarkLines requires a symmetric pose" );
	}

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( p.conformation()) );

	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	Size nsubunits = symm_info->subunits();
	Size nres_asu = symm_info->num_independent_residues();
	// Size nres_monomer = symm_info->num_independent_residues() - 1; // Unused variable causes warning.

	FindZnCoordinatingResidues finder;
	finder.find_coordinating_residues( p );

	// Now create the needed coordination data
	symdes_znx_coordination_data coordination_data;
	coordination_data.zinc_data_[ 1 ].chain_ = 'B';
	coordination_data.zinc_data_[ 1 ].resind_ = nres_asu;
	coordination_data.zinc_data_[ 2 ].chain_ = 'B';
	coordination_data.zinc_data_[ 2 ].resind_ = nres_asu;
	coordination_data.zinc_data_[ 3 ].chain_ = char( int('A') + 2*nsubunits - 1 );
	coordination_data.zinc_data_[ 3 ].resind_ = nres_asu*nsubunits;
	for ( Size ii = 1; ii <= 6; ++ii ) {

		if ( ii <= 2 ) {
			coordination_data.protein_data_[ ii ].chain_ = 'A';
			coordination_data.protein_data_[ ii ].resindex_ = finder.resinds()[ ii ];
			coordination_data.protein_data_[ ii ].name3_ =	p.residue( finder.resinds()[ ii ] ).name3();
		} else if ( ii <= 4 ) {
			coordination_data.protein_data_[ ii ].chain_ = 'C';
			coordination_data.protein_data_[ ii ].resindex_ = finder.resinds()[ ii ] + nres_asu;
			coordination_data.protein_data_[ ii ].name3_ =	p.residue( finder.resinds()[ ii ] ).name3();
		} else {
			coordination_data.protein_data_[ ii ].chain_ = 'A';
			coordination_data.protein_data_[ ii ].resindex_ = finder.resinds()[ ii - 2 ]; // ii - 2 == 3 or 4
			coordination_data.protein_data_[ ii ].name3_ =	p.residue( finder.resinds()[ ii - 2 ] ).name3();
		}
		coordination_data.protein_data_[ ii ].exgeom_index_ = utility::to_string( coordination_data.protein_data_[ ii ].name3_ == "HIS" ? 1 : 2); // UGLY HACK!  Assumes that HIS is described in the first block of a variable constraint and that ASP/GLU are described in the second block. This is true for the one particular input file I'm working with, but may not be true generally.

	}

	add_znx_coordination_remark_lines_to_pose( p, coordination_data );

}

void InsertZincCoordinationRemarkLines::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	//noop
}


/////////////////////////////////////////////////////////////////////

core::pack::task::operation::TaskOperationOP
DisableZnCoordinationResiduesTaskOpCreator::create_task_operation() const
{ return new DisableZnCoordinationResiduesTaskOp;}

std::string DisableZnCoordinationResiduesTaskOpCreator::keyname() const
{
	return "DisableZnCoordinationResiduesTaskOp";
}

////////////////////////////////////////////////////////////////////

DisableZnCoordinationResiduesTaskOp::DisableZnCoordinationResiduesTaskOp() {}
DisableZnCoordinationResiduesTaskOp::~DisableZnCoordinationResiduesTaskOp() {}

DisableZnCoordinationResiduesTaskOp::TaskOperationOP
DisableZnCoordinationResiduesTaskOp::clone() const
{
	return new DisableZnCoordinationResiduesTaskOp;
}

void DisableZnCoordinationResiduesTaskOp::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	//FindZnCoordinatingResidues finder;
	//finder.find_coordinating_residues( pose );
	//for ( Size ii = 1; ii <= 4; ++ii ) {
	//	task.nonconst_residue_task( finder.resinds()[ ii ] ).prevent_repacking();
	//}

	using namespace core::conformation::symmetry;

	if ( ! dynamic_cast<SymmetricConformation const * > ( & pose.conformation()) ) {
		utility_exit_with_message( "InsertZincCoordinationRemarkLines requires a symmetric pose" );
	}

	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );

	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	Size nsubunits = symm_info->subunits();
	Size nres_asu = symm_info->num_independent_residues();
	// Size nres_monomer = symm_info->num_independent_residues() - 1; // Unused varaible causes warning.

	task.nonconst_residue_task( nres_asu ).prevent_repacking();

	protocols::enzdes::EnzdesConstraintReporter reporter1;
	reporter1.ligand_resno( nres_asu );
	reporter1.find_constraints_to_ligand( pose );

	for ( core::Size ii = 1; ii <= reporter1.constrained_nonligand_atoms().size(); ++ii ) {
		core::id::AtomID iiatid = reporter1.constrained_nonligand_atoms()[ ii ];
		TR << "Disabling residue " << iiatid.rsd() << std::endl;
		task.nonconst_residue_task( iiatid.rsd() ).prevent_repacking();
	}

	protocols::enzdes::EnzdesConstraintReporter reporter2;
	reporter2.ligand_resno( nres_asu*nsubunits );
	reporter2.find_constraints_to_ligand( pose );

	task.nonconst_residue_task( nres_asu*nsubunits ).prevent_repacking();

	for ( core::Size ii = 1; ii <= reporter2.constrained_nonligand_atoms().size(); ++ii ) {
		core::id::AtomID iiatid = reporter2.constrained_nonligand_atoms()[ ii ];
		TR << "Disabling residue " << iiatid.rsd() << std::endl;
		task.nonconst_residue_task( iiatid.rsd() ).prevent_repacking();
	}

}

void DisableZnCoordinationResiduesTaskOp::parse_tag( TagCOP, DataMap & )
{
}

/////////////////////////////////////////////////////////////////////////////////

ZnCoordNumHbondCalculator::ZnCoordNumHbondCalculator()
{
	finder_.fail_on_absent_coordinators( false );
}

core::pose::metrics::PoseMetricCalculatorOP
ZnCoordNumHbondCalculator::clone() const {
	return new ZnCoordNumHbondCalculator;
}

void ZnCoordNumHbondCalculator::notify_structure_change() {
	core::pose::metrics::EnergyDependentCalculator::notify_structure_change();
	nhbcalc_.notify_structure_change();
}

void ZnCoordNumHbondCalculator::notify_energy_change() {
	core::pose::metrics::EnergyDependentCalculator::notify_energy_change();
	nhbcalc_.notify_energy_change();
}



void ZnCoordNumHbondCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const
{
   if ( key == "all_Hbonds" ) {
     basic::check_cast( valptr, &all_Hbonds_, "all_Hbonds expects to return a Size" );
     (static_cast< basic::MetricValue< core::Size > * >(valptr))->set( all_Hbonds_ );

	 } else if ( key == "atom_Hbonds" ) {
     basic::check_cast( valptr, &atom_Hbonds_, "atom_Hbonds expects to return a id::AtomID_Map< Size >" );
     (static_cast< basic::MetricValue< core::id::AtomID_Map< Size > > * >(valptr))->set( atom_Hbonds_ );

   } else if ( key == "residue_Hbonds" ) {
     basic::check_cast( valptr, &residue_Hbonds_, "residue_Hbonds expects to return a utility::vector1< Size >" );
     (static_cast<basic::MetricValue<utility::vector1< Size > > * >(valptr))->set( residue_Hbonds_ );

   } else {
     basic::Error() << "NumberHbondsCalculator cannot compute the requested metric " << key << std::endl;
     utility_exit();
   }

}

std::string
ZnCoordNumHbondCalculator::print( std::string const & key ) const
{
  if ( key == "all_Hbonds" ) {
    return utility::to_string( all_Hbonds_ );
  } else if ( key == "atom_Hbonds" ) {
    basic::Error() << "id::AtomID_Map< Size > has no output operator, for metric " << key << std::endl;
    utility_exit();
  } else if ( key == "residue_Hbonds" ) {
    return utility::to_string( residue_Hbonds_ );
  }

  basic::Error() << "ZnCoordNumHbondCalculator cannot compute metric " << key << std::endl;
  utility_exit();
  return "";
}

void
ZnCoordNumHbondCalculator::recompute( core::pose::Pose const & this_pose )
{
	TR << "ZnCoordNumHbondCalculator::recompute" << std::endl;
	//nhbcalc_.notify_structure_change();
	basic::MetricValue< core::Size > mv_all_Hbonds;
	nhbcalc_.get( "all_Hbonds", mv_all_Hbonds, this_pose );
	all_Hbonds_ = mv_all_Hbonds.value();

	basic::MetricValue< core::id::AtomID_Map< core::Size > > mv_atom_Hbonds;
	nhbcalc_.get( "atom_Hbonds", mv_atom_Hbonds, this_pose );
	atom_Hbonds_ = mv_atom_Hbonds.value();

	basic::MetricValue< utility::vector1< core::Size > > mv_residue_Hbonds;
	nhbcalc_.get( "residue_Hbonds", mv_residue_Hbonds, this_pose );
	residue_Hbonds_ = mv_residue_Hbonds.value();

	finder_.find_coordinating_residues( this_pose );
	for ( Size ii = 1; ii <= finder_.atomids().size(); ++ii ) {
		TR << "Incrementing hbonds for " << finder_.atomids()[ ii ].rsd() << " atom " << this_pose.residue(finder_.atomids()[ ii ].rsd() ).atom_name( finder_.atomids()[ ii ].atomno() ) << std::endl;
		atom_Hbonds_[ finder_.atomids()[ ii ] ] += 1; // add one hydrogen bond to every atom coordinating zinc
		residue_Hbonds_[ finder_.resinds()[ ii ] ] += 1; // add one hydrogen bond to every residue coordinating zinc
	}
	all_Hbonds_ += 4;

}

/////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
LoadZnCoordNumHbondCalculatorMoverCreator::create_mover() const
{
	return new LoadZnCoordNumHbondCalculatorMover;
}
std::string LoadZnCoordNumHbondCalculatorMoverCreator::keyname() const
{
	return "LoadZnCoordNumHbondCalculatorMover";
}

/////////////////////////////////////////////////////////////////////////////////

LoadZnCoordNumHbondCalculatorMover::LoadZnCoordNumHbondCalculatorMover() :
	protocols::moves::Mover( "LoadZnCoordNumHbondCalculatorMover" )
{}

LoadZnCoordNumHbondCalculatorMover::~LoadZnCoordNumHbondCalculatorMover()
{}

protocols::moves::MoverOP
LoadZnCoordNumHbondCalculatorMover::clone() const {
	return new LoadZnCoordNumHbondCalculatorMover;
}

std::string
LoadZnCoordNumHbondCalculatorMover::get_name() const { return "LoadZnCoordNumHbondCalculatorMover"; }

void
LoadZnCoordNumHbondCalculatorMover::apply( core::pose::Pose & )
{
	TR << "LoadZnCoordNumHbondCalculatorMover::apply" << std::endl;
	using core::pose::metrics::CalculatorFactory;
	if ( CalculatorFactory::Instance().check_calculator_exists( "bur_unsat_calc_default_hbond_calc" ) ) {
		CalculatorFactory::Instance().remove_calculator( "bur_unsat_calc_default_hbond_calc" );
	}
	CalculatorFactory::Instance().register_calculator( "bur_unsat_calc_default_hbond_calc", core::pose::metrics::PoseMetricCalculatorOP( new ZnCoordNumHbondCalculator ) );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void LoadZnCoordNumHbondCalculatorMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{}




}
}


