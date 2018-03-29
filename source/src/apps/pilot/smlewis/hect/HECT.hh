// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/smlewis/HECT/HECT.hh
/// @brief This mover is a base mover for my pair of HECT apps.  It contains shared methods, mostly for setting up the odd chemical linkage.  See the other HECT files.
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/task_operations/RestrictByCalculatorsOperation.hh>
#include <protocols/task_operations/RestrictToNeighborhoodOperation.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/pose_metric_calculators/InterGroupNeighborsCalculator.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh> //BoundFunc
#include <core/scoring/constraints/AmbiguousConstraint.hh>

#include <core/scoring/rms_util.hh>

//movers
#include <protocols/simple_moves/BackboneMover.hh> //SmallMover
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/simple_moves/TorsionDOFMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/moves/MonteCarlo.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <numeric/xyz.functions.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <set>
#include <utility>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.pilot.smlewis.HECT" );

namespace basic { namespace options { namespace OptionKeys {
basic::options::IntegerOptionKey const e3_hinge_start_resnum("e3_hinge_start_resnum");
basic::options::IntegerOptionKey const e3_hinge_stop_resnum("e3_hinge_stop_resnum");
basic::options::StringOptionKey const e3_hinge_chain("e3_hinge_chain");
basic::options::IntegerOptionKey const e3_catalytic_resnum("e3_catalytic_resnum");

basic::options::IntegerOptionKey const UBQ_tail_start_resnum("UBQ_tail_start_resnum");
basic::options::StringOptionKey const UBQ_tail_chain("UBQ_tail_chain");

basic::options::StringOptionKey const e2_chain("e2_chain");
basic::options::IntegerOptionKey const e2_catalytic_resnum("e2_catalytic_resnum");

basic::options::BooleanOptionKey const debug_skip_fragment_generation("debug_skip_fragment_generation");
basic::options::BooleanOptionKey const debug("debug");

basic::options::IntegerOptionKey const cycles("cycles");
basic::options::IntegerOptionKey const repack_cycles("cycles");

basic::options::RealOptionKey const catalytic_cst_sd("catalytic_cst_sd");
basic::options::RealOptionKey const ubq_cst_sd("ubq_cst_sd");

}}}//basic::options::OptionKeys

/// @brief utility function to register options
void register_options(){
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	option.add( e3_hinge_start_resnum, "" );
	option.add( e3_hinge_stop_resnum, "" );
	option.add( e3_hinge_chain, "" );
	option.add( e3_catalytic_resnum, "" );

	option.add( UBQ_tail_start_resnum, "" );
	option.add( UBQ_tail_chain, "" );

	option.add( e2_chain, "" );
	option.add( e2_catalytic_resnum, "" );

	option.add( debug_skip_fragment_generation, "" ).def(false);
	option.add( debug, "" ).def(false);

	option.add( cycles, "run refinement for this many cycles" ).def(50);
	option.add( repack_cycles, "repack every this many cycles" ).def(50);

	option.add( catalytic_cst_sd, "sd (weight) of catalytic constraint").def(1.0);
	option.add( ubq_cst_sd, "sd (weight) of UBQ/E3 constraints").def(3.0);

	return;
}

/// @brief HECT mover
class HECTMover : public protocols::moves::Mover {

protected: //enum for atomID vector
	enum atomID {
		CYS_C  = 1,
		CYS_CA = 2,
		CYS_CB = 3,
		CYS_SG = 4,
		GLY_C  = 5,
		GLY_CA = 6,
		GLY_N  = 7,
		GLY2_C = 8,
		E3_cat_SG = 9,
		atomID_tot = E3_cat_SG
	};

	//putting the data members first - for ease of organizing ctor
protected: //protected does not require writing getters/setters for all these members for this use-once code

	/// @brief flag for when init function has run.
	bool init_for_input_yet_;

	/// @brief variables for points of interest in the protein (brackets for movers) - this may be overkill
	char e3chain_, e2chain_, ubqchain_;
	core::Size e3chain_num_, e3chain_start_, e3chain_end_, e3_hinge_start_, e3_hinge_stop_, e3_catalytic_res_,
		ubqchain_num_, ubqchain_start_, ubqchain_end_, utail_start_, utail_stop_,
		e2chain_num_, e2chain_start_, e2chain_end_, e2_catalytic_res_;

	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP movemap_;

	/// @brief fragment sets for hinge region and tail region
	core::fragment::FragSetCOP fragset3mer_e3_hinge_, fragset3mer_ubq_tail_;

	/// @brief this mover ignores the apply-passed PDB, instead returning to the starting pose set up by the init function.  Under the assumption that the starting pose is not changing this lets up skip running the setup code every time.
	core::pose::PoseCOP fixed_starting_pose_, xtal_pose_;

	/// @brief vector contains atomIDs for thioester bond and atoms before/after bond to determine various torsions - indexed by enum
	utility::vector1< core::id::AtomID > atomIDs;

public:
	HECTMover() :  init_for_input_yet_(false),
		e3chain_(0), e2chain_(0), ubqchain_(0),
		e3chain_num_(0), e3chain_start_(0), e3chain_end_(0), e3_hinge_start_(0), e3_hinge_stop_(0), e3_catalytic_res_(0),
		ubqchain_num_(0), ubqchain_start_(0), ubqchain_end_(0), utail_start_(0), utail_stop_(0),
		e2chain_num_(0), e2chain_start_(0), e2chain_end_(0), e2_catalytic_res_(0),
		fullatom_scorefunction_(NULL), task_factory_(NULL), movemap_(NULL), fragset3mer_e3_hinge_(NULL), fragset3mer_ubq_tail_(NULL),
		fixed_starting_pose_(NULL), xtal_pose_(NULL),
		atomIDs(atomID_tot, core::id::AtomID::BOGUS_ATOM_ID() )

	{
		//set up fullatom scorefunction
		using namespace core::scoring;
		fullatom_scorefunction_ = get_score_function();
		TR << "Using fullatom scorefunction (TALARIS_2013), \n" << *fullatom_scorefunction_ << std::flush;

		//read native structure for CA RMSD
		core::pose::Pose xtal_pose;
		core::import_pose::pose_from_file(xtal_pose, basic::options::option[basic::options::OptionKeys::in::file::native].value(), core::import_pose::PDB_file);
		xtal_pose_ = new core::pose::Pose(xtal_pose);
	}

	/// @brief helper code for fragments generation
	core::fragment::FragSetCOP make_frags(core::Size const start, core::Size const stop, std::string const & seq){

		core::Size const frags_length(3); //magic number: 3mer fragments!!
		core::fragment::FragSetOP fragset(new core::fragment::ConstantLengthFragSet( frags_length ));
		std::string ss_string(frags_length, 'L');
		core::fragment::FragDataOPs list;

		if ( !basic::options::option[basic::options::OptionKeys::debug_skip_fragment_generation].value() ) {
			for ( core::Size j = start; j <= stop-frags_length+1; ++j ) {
				std::string const seqsubstr(seq, j-1, frags_length); //j-1 accounts for string [] from 0
				TR << "adding frame, start at " << j << " go for " << frags_length << " to " << j+frags_length << " seq " << seqsubstr << std::endl;
				list =  core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa( ss_string, seqsubstr, 200, false ); //magic number: 200 fragments per position (not duplicated - this will be like robetta server fragments)
				core::fragment::FrameOP frame;
				frame = new core::fragment::Frame( j );
				frame->add_fragment( list );
				fragset->add( frame );
			}//hinge residues
		}//debugging skip
		return fragset;
	}

	/// @brief parse_options will grab things from the option system and store them in local data - the stored data are used in the later setup functions
	void
	parse_options(core::pose::Pose const & pose) {

		//determine where the flexible regions are
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//there should only be 3 chains in the input pose
		runtime_assert(pose.conformation().num_chains() == 3);

		//e3 points of interest
		e3chain_ = (option[e3_hinge_chain].value()[0]); //just take the first one
		e3_hinge_start_ = pose.pdb_info()->pdb2pose().find(e3chain_, option[e3_hinge_start_resnum].value());
		e3_hinge_stop_  = pose.pdb_info()->pdb2pose().find(e3chain_, option[e3_hinge_stop_resnum ].value());
		e3chain_num_ = (pose.chain(e3_hinge_start_));
		e3chain_start_ = (pose.conformation().chain_begin(e3chain_num_));
		e3chain_end_ = (pose.conformation().chain_end(e3chain_num_));
		e3_catalytic_res_ = (pose.pdb_info()->pdb2pose().find(e3chain_, option[e3_catalytic_resnum ].value()));

		//ubq points of interest
		ubqchain_ = (option[UBQ_tail_chain].value()[0]); //just take the first one
		utail_start_ = pose.pdb_info()->pdb2pose().find(ubqchain_, option[UBQ_tail_start_resnum].value());
		ubqchain_num_ = (pose.chain(utail_start_));
		utail_stop_ = pose.conformation().chain_end(ubqchain_num_);
		ubqchain_start_ = (pose.conformation().chain_begin(ubqchain_num_));
		ubqchain_end_ = utail_stop_;

		//E2 start/end/points of interest
		e2chain_ = (option[e2_chain].value()[0]); //just take the first one
		e2_catalytic_res_ = (pose.pdb_info()->pdb2pose().find(e2chain_, option[e2_catalytic_resnum].value()));
		e2chain_num_ = (pose.chain(e2_catalytic_res_));
		e2chain_start_ = (pose.conformation().chain_begin(e2chain_num_));
		e2chain_end_ = (pose.conformation().chain_end(e2chain_num_));

		TR << "resnums:" << std::endl
			<< "e3 chain start catalytic hingestart hingeend stop " << e3chain_num_ << " " << e3chain_start_ << " "
			<< e3_catalytic_res_ << " " << e3_hinge_start_ << " " << e3_hinge_stop_ << " " << e3chain_end_ << std::endl
			<< "e2 chain start catalytic stop " << e2chain_num_ << " " << e2chain_start_ << " "
			<< e2_catalytic_res_ << " " << e2chain_end_ << std::endl
			<< "ubq chain start flexstart stop " << " " << ubqchain_num_ << " " << ubqchain_start_ << " "
			<< utail_start_ << " " << ubqchain_end_ << std::endl;
		TR << "num chains " << pose.conformation().num_chains() << std::endl;
		return;
	}

	/// @brief set_up_foldtree
	void
	set_up_foldtree(core::pose::Pose & pose) {

		core::Size const E2_E3_JUMP(1);

		using namespace core::kinematics;
		FoldTree foldtree(pose.size());
		foldtree.clear();
		foldtree.add_edge( Edge(e3chain_start_, e3chain_end_, Edge::PEPTIDE));
		foldtree.add_edge( Edge(e2chain_start_, e2_catalytic_res_, Edge::PEPTIDE));
		foldtree.add_edge( Edge(e2_catalytic_res_, e2chain_end_, Edge::PEPTIDE));
		foldtree.add_edge( Edge(e2chain_start_, e3chain_start_, E2_E3_JUMP));
		foldtree.add_edge( Edge(ubqchain_end_, ubqchain_start_, Edge::PEPTIDE)); //entirety of ubiquitin
		foldtree.add_edge( Edge(e2_catalytic_res_, ubqchain_end_, Edge::CHEMICAL, "SG", "C", true ));

		foldtree.reorder(e3chain_start_);
		pose.fold_tree(foldtree);

		TR << "new fold tree " << pose.fold_tree() << std::endl;
		return;
	}

	/// @brief some shared setup of task factory
	void
	set_up_taskfactory(){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		task_factory_ = new TaskFactory;
		task_factory_->push_back( new InitializeFromCommandline );
		if ( option[ packing::resfile ].user() ) {
			task_factory_->push_back( new ReadResfile );
		}
		task_factory_->push_back( new IncludeCurrent );
		task_factory_->push_back( new RestrictToRepacking );

		//prevent repacking at chemical linkage!
		operation::PreventRepackingOP PROP(new operation::PreventRepacking);
		PROP->include_residue(e2_catalytic_res_);
		task_factory_->push_back( PROP );


		//Sharing all this code for both HECT movers is slightly inefficient - you could remove some group pairs for the ubq only case - but I don't expect it to make a big difference

		//first we build residue sets for the different domains and regions of interest
		std::set< core::Size > E3C, E3N, E3hinge, E2, UBQ, UBQtail;
		for ( core::Size i(e3chain_start_);   i<=e3_hinge_start_; ++i ) E3N.insert(i);
		for ( core::Size i(e3_hinge_start_); i<=e3_hinge_stop_;  ++i ) E3hinge.insert(i);
		for ( core::Size i(e3_hinge_stop_);  i<=e3chain_end_;     ++i ) E3C.insert(i);
		for ( core::Size i(e2chain_start_);   i<=e2chain_end_;     ++i ) E2.insert(i);
		for ( core::Size i(ubqchain_start_);  i<=utail_start_;    ++i ) UBQ.insert(i);
		for ( core::Size i(utail_start_);    i<=ubqchain_end_;    ++i ) UBQtail.insert(i);


		//borrowing IGNC's typedefs to set up the group pairs
		using namespace protocols::pose_metric_calculators;
		InterGroupNeighborsCalculator::group_set groups;
		//push back pairs of sets for the regions we need:
		// E3hinge_self
		// E3N_E3hinge
		// E3hinge_E3C
		// E3N_E3C
		// E3C_UBQ
		// UBQ_UBQtail
		// UBQtail_UBQtail
		// UBQtail_E2
		// UBQ_E2
		// E3C_E2
		groups.push_back(std::make_pair(E3hinge, E3hinge));
		groups.push_back(std::make_pair(E3N, E3hinge));
		groups.push_back(std::make_pair(E3hinge, E3C));
		groups.push_back(std::make_pair(E3N, E3C));
		groups.push_back(std::make_pair(E3C, UBQ));
		groups.push_back(std::make_pair(UBQ, UBQtail));
		groups.push_back(std::make_pair(UBQtail, UBQtail));
		groups.push_back(std::make_pair(UBQtail, E2));
		groups.push_back(std::make_pair(UBQ, E2));
		groups.push_back(std::make_pair(E3C, E2));

		//make the calculator
		std::string const calc_g("IGNC_g");
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc_g, new protocols::pose_metric_calculators::InterGroupNeighborsCalculator(groups) );

		//this is the constructor parameter for the TaskOperation - pairs of calculators and calculations to perform
		utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
		calcs_and_calcns.push_back(std::make_pair(calc_g, "neighbors"));

		using protocols::task_operations::RestrictByCalculatorsOperation;
		task_factory_->push_back( new RestrictByCalculatorsOperation( calcs_and_calcns ) );

		return;
	}

	/// @brief ubq_constraints will set up some ubiquitin to E3 constraints based on the input crystal structures.  To refresh your memory, this code runs on the PDBs 3JW0 and 3JVZ, from reference Kamadurai HB, Souphron J, Scott DC, Duda DM, Miller DJ, Stringer D, Piper RC, Schulman BA.  Insights into ubiquitin transfer cascades from a structure of a UbcH5B approximately ubiquitin-HECT(NEDD4L) complex.  Mol Cell. 2009 Dec 25;36(6):1095-102.  Anyway, the constraints are drawn from some mutational experiments in the paper, which showed that these interactions matter for functionality.
	void
	ubq_constraints(core::pose::Pose & pose) {
		using namespace core::scoring::constraints;
		using core::id::AtomID;

		//build the BoundFuncs we are going to use
		using namespace basic::options;
		BoundFuncOP chargepair( new BoundFunc( 0, 3, option[OptionKeys::ubq_cst_sd].value(), "chargepair") );
		BoundFuncOP hbond( new BoundFunc( 0, 3, option[OptionKeys::ubq_cst_sd].value(), "hbond") );

		//here we are hardcoding some constraints between the UBQ and E3C regions - this is based on contacts observed in the publication
		//E3C     UBQ
		//K915 NZ D39 OD2 (but switchable to OD1?)
		//K940 NZ G35 O
		//L916 O  Q40 NE2
		//L916 N  Q40 OE1

		core::chemical::ResidueTypeSetCAP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		core::chemical::ResidueType const &
			lys_rt( fa_standard->name_map("LYS") ),
			leu_rt( fa_standard->name_map("LEU") ),
			asp_rt( fa_standard->name_map("ASP") ),
			gln_rt( fa_standard->name_map("GLN") ),
			gly_rt( fa_standard->name_map("GLY") );

		//E3C     UBQ                //this is an ion pair
		//K915 NZ D39 OD2 (but switchable to OD1?)
		core::Size const NZ(lys_rt.atom_index("NZ"));
		core::Size const asp_OD2(asp_rt.atom_index("OD2"));
		core::Size const asp_OD1(asp_rt.atom_index("OD1"));

		core::id::AtomID const k915_NZ(NZ, pose.pdb_info()->pdb2pose(e3chain_, 915));
		core::Size D39(pose.pdb_info()->pdb2pose(ubqchain_, 39));
		AmbiguousConstraintOP ambig_ionpair( new AmbiguousConstraint() );
		ambig_ionpair->add_individual_constraint(new AtomPairConstraint( AtomID(asp_OD2, D39), k915_NZ, chargepair ) );
		ambig_ionpair->add_individual_constraint(new AtomPairConstraint( AtomID(asp_OD1, D39), k915_NZ, chargepair ) );
		pose.add_constraint( ambig_ionpair );

		//E3C     UBQ                //this is an hbond
		//K940 NZ G35 O
		core::id::AtomID const k940_NZ(NZ, pose.pdb_info()->pdb2pose(e3chain_, 940));
		core::Size const gly_O(gly_rt.atom_index("O"));
		core::Size const G35(pose.pdb_info()->pdb2pose(ubqchain_, 35));
		AtomPairConstraintOP hbond1( new AtomPairConstraint( AtomID(gly_O, G35), k940_NZ, hbond));
		pose.add_constraint( hbond1 );

		//E3C     UBQ                //this is a pair of hbonds
		//L916 O  Q40 NE2
		//L916 N  Q40 OE1
		core::Size const Q40(pose.pdb_info()->pdb2pose(ubqchain_, 40));
		core::Size const L916(pose.pdb_info()->pdb2pose(e3chain_, 916));
		AtomPairConstraintOP hbond2a( new AtomPairConstraint( AtomID(gln_rt.atom_index("NE2"), Q40), AtomID(leu_rt.atom_index("O"), L916), hbond));
		AtomPairConstraintOP hbond2b( new AtomPairConstraint( AtomID(gln_rt.atom_index("OE1"), Q40), AtomID(leu_rt.atom_index("N"), L916), hbond));
		pose.add_constraint( hbond2a );
		pose.add_constraint( hbond2b );

		fullatom_scorefunction_->set_weight( core::scoring::atom_pair_constraint, 1.0 ); //option weight control?

		return;
	}

	/// @brief build_AtomID_vec builds the atomID vector used for setting up the TorsionDOFMovers and used to calculate atom distances for later scoring
	void
	build_AtomID_vec(){
		core::chemical::ResidueTypeSetCAP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		core::chemical::ResidueType const &
			cyx_rt( fa_standard->name_map("CYX") ),
			cyz_rt( fa_standard->name_map("CYZ") ),
			gly_rt( fa_standard->name_map("GLY") );

		using core::id::AtomID;
		atomIDs[ CYS_C     /*1*/] = AtomID(cyx_rt.atom_index("C"), e2_catalytic_res_);
		atomIDs[ CYS_CA    /*2*/] = AtomID(cyx_rt.atom_index("CA"), e2_catalytic_res_);
		atomIDs[ CYS_CB    /*3*/] = AtomID(cyx_rt.atom_index("CB"), e2_catalytic_res_);
		atomIDs[ CYS_SG    /*4*/] = AtomID(cyx_rt.atom_index("SG"), e2_catalytic_res_);
		atomIDs[ GLY_C     /*5*/] = AtomID(gly_rt.atom_index("C"), ubqchain_end_);
		atomIDs[ GLY_CA    /*6*/] = AtomID(gly_rt.atom_index("CA"), ubqchain_end_);
		atomIDs[ GLY_N     /*7*/] = AtomID(gly_rt.atom_index("N"), ubqchain_end_);
		atomIDs[ GLY2_C    /*8*/] = AtomID(gly_rt.atom_index("C"), ubqchain_end_-1);
		atomIDs[ E3_cat_SG /*9*/] = AtomID(cyz_rt.atom_index("SG"), e3_catalytic_res_);
		return;
	}

	virtual ~HECTMover(){};

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return true; }

	/// @details This function is specific to the original system for which this code was written - if you are not trying to duplicate the initial results you should remove it!
	void create_extra_output( core::pose::Pose const & pose ){

		//find Job
		using protocols::jd2::JobDistributor;
		protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );

		//print catalytic distance
		core::Real const cat_dist(distance( pose.xyz(atomIDs[GLY_C]), pose.xyz(atomIDs[E3_cat_SG])));
		job_me->add_string_real_pair("catalytic_distance", cat_dist);

		//print RMSD to start pose
		job_me->add_string_real_pair("RMSD_start", core::scoring::CA_rmsd(pose, *fixed_starting_pose_ ));

		//print RMSD to native pose
		job_me->add_string_real_pair("RMSD_native", core::scoring::CA_rmsd(pose, *xtal_pose_ ));

		return;
	}

	/// @brief this apply function holds the main apply for BOTH HECT_UBQ and HECT_ALL.  It has a function call to a function which will handle the different aspect for those two peices of code; the remainder of the differences are encoded in the setup function defined by each subclass.
	virtual
	void
	apply( core::pose::Pose & pose ){

		if ( !init_for_input_yet_ ) init_on_new_input(pose);
		pose = *fixed_starting_pose_; //copy operation

		/*
		movers we'll need:
		bb cycles:
		small/shear/fragment movers for hinge region
		small/shear/fragment movers for UBQ tail
		TorsionDOF movers for chemical linkage
		RotamerTrialsMover
		MinMover
		packing cycles:
		PackRotamersMover
		TaskAwareMinMover
		*/

		using namespace protocols::moves;

		protocols::moves::RandomMoverOP random_mover( new protocols::moves::RandomMover() );

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *fullatom_scorefunction_, 0.8 ) );

		/////////////////////////////////rotamer trials mover///////////////////////////////////////////
		using protocols::minimization_packing::RotamerTrialsMoverOP;
		using protocols::minimization_packing::EnergyCutRotamerTrialsMover;
		protocols::minimization_packing::RotamerTrialsMoverOP rt_mover(new protocols::minimization_packing::EnergyCutRotamerTrialsMover(
			fullatom_scorefunction_,
			task_factory_,
			mc,
			0.01 /*energycut*/ ) );

		//////////////////////////TorsionDOFMovers for the chemical linkage////////////////////////////////
		//utility::vector1< protocols::simple_moves::TorsionDOFMoverOP > DOF_movers;
		//we are going to iterate over the ordered quadruplets describing the torsional degrees of freedom involved in the thioester bond, like so: 1,2,3,4; 2,3,4,5; 3,4,5,6
		//This iterates from the cysteine backbone C to the second glycine backbone C.  See the enum names for which atoms are which
		for ( core::Size i(1); i <= GLY2_C-3; ++i ) {
			protocols::simple_moves::TorsionDOFMoverOP DOF_mover(new protocols::simple_moves::TorsionDOFMover);
			DOF_mover->set_DOF(atomIDs[i], atomIDs[i+1], atomIDs[i+2], atomIDs[i+3]);
			DOF_mover->check_mmt(true);
			DOF_mover->temp(0.4);
			DOF_mover->set_angle_range(-180, 180);
			DOF_mover->tries(1000);
			//DOF_movers.push_back(DOF_mover); //we are really keeping this handle only as a backup / for debugging
			random_mover->add_mover(DOF_mover, 1.0);
		}

		//  //DOF movers test
		//   core::pose::Pose copy(pose);
		//   std::ostringstream outputfilename;
		//   for(core::Size i(1); i<=DOF_movers.size(); ++i){
		//    DOF_movers[i]->apply(copy);
		//    outputfilename.str(""); //clears stream
		//    outputfilename << "torsionDOF_" << i << ".pdb";
		//    copy.dump_pdb(outputfilename.str());
		//    copy = pose;
		//   }

		///////////////////////////////////fragments////////////////////////////////////////////////
		using protocols::simple_moves::ClassicFragmentMover;
		using protocols::simple_moves::ClassicFragmentMoverOP;
		add_frag_mover(random_mover); //adds hinge fragment mover, if appropriate

		ClassicFragmentMoverOP frag_mover_tail = new ClassicFragmentMover(fragset3mer_ubq_tail_, movemap_);
		frag_mover_tail->enable_end_bias_check(false);
		random_mover->add_mover(frag_mover_tail, 1.0);

		/////////////////////////////////////SmallMover, protocols::simple_moves::ShearMover//////////////////////////////////
		protocols::simple_moves::BackboneMoverOP small_mover = new protocols::simple_moves::SmallMover(movemap_, 0.8, 0);
		small_mover->angle_max( 'H', 4.0 );
		small_mover->angle_max( 'E', 4.0 );
		small_mover->angle_max( 'L', 4.0 );
		random_mover->add_mover(small_mover, 1.0);

		protocols::simple_moves::BackboneMoverOP shear_mover = new protocols::simple_moves::ShearMover(movemap_, 0.8, 0);
		shear_mover->angle_max( 'H', 4.0 );
		shear_mover->angle_max( 'E', 4.0 );
		shear_mover->angle_max( 'L', 4.0 );
		random_mover->add_mover(shear_mover, 1.0);

		/////////////////////////////generate full repack mover////////////////////////////////////////
		protocols::minimization_packing::PackRotamersMoverOP pack_mover = new protocols::minimization_packing::PackRotamersMover;
		pack_mover->task_factory( task_factory_ );
		pack_mover->score_function( fullatom_scorefunction_ );

		///////////////////////////////////////////Minimizer mover and TAMinmover///////////////////////////
		protocols::minimization_packing::MinMoverOP min_mover = new protocols::minimization_packing::MinMover(
			movemap_,
			fullatom_scorefunction_,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true /*use_nblist*/ );

		using protocols::minimization_packing::TaskAwareMinMoverOP;
		using protocols::minimization_packing::TaskAwareMinMover;
		protocols::minimization_packing::TaskAwareMinMoverOP TAmin_mover = new protocols::minimization_packing::TaskAwareMinMover(min_mover, task_factory_);

		///////////////////////package RT/min for JumpOutMover////////////////////////////////////////
		protocols::moves::SequenceMoverOP RT_min_seq( new protocols::moves::SequenceMover );
		RT_min_seq->add_mover(rt_mover);
		RT_min_seq->add_mover(min_mover);

		protocols::moves::JumpOutMoverOP bb_if_RT_min( new protocols::moves::JumpOutMover(
			random_mover,
			RT_min_seq,
			fullatom_scorefunction_,
			20.0)); //20 score units

		//run loop
		using basic::options::option;
		using namespace basic::options::OptionKeys;
		core::Size const applies = option[ cycles ].value(); //default 5
		core::Size const repackcycles = option[ repack_cycles ].value();
		TR << "   Current     Low    total cycles =" << applies << std::endl;
		for ( core::Size i(1); i <= applies; ++i ) {
			if ( (i % repackcycles == 0) || (i == applies) ) { //full repack
				pack_mover->apply(pose);
				TAmin_mover->apply(pose);
			} else {
				bb_if_RT_min->apply(pose);
			}

			mc->boltzmann(pose);
			if ( mc->mc_accepted() ) TR << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;

		}//end the exciting for loop
		mc->recover_low( pose );

		//let's store some energies/etc of interest
		if ( true ) create_extra_output(pose);

		(*fullatom_scorefunction_)(pose);
		set_last_move_status(protocols::moves::MS_SUCCESS); //this call is unnecessary but let's be safe
		return;
	}

	/// @brief default version does nothing
	virtual
	void
	add_frag_mover(protocols::moves::RandomMoverOP /*random_mover*/) {return;}

	/// @brief pure virtual version, this is where most differences between ubq and all version lie
	virtual
	void
	init_on_new_input(core::pose::Pose & pose) = 0;

};
