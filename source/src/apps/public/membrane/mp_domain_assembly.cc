// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_domain_assembly.cc
/// @brief   Build a full-length membrane protein model from a bunch of poses and a fasta file
/// @author  JKLeman (julia.koehler.leman@gmail.com)

// App headers
#include <devel/init.hh>

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/AddMembraneMover.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/pose/annotated_sequence.hh>

#include <protocols/relax/membrane/MPRangeRelaxMover.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/conformation/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <protocols/jd2/util.hh>

#include <utility/string_util.hh>
#include <utility/minmax.hh>
#include <numeric/random/random.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::pose;
using namespace core::conformation;
using namespace core::conformation::membrane;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::membrane;
using namespace protocols::simple_moves;
using namespace core::chemical;
using namespace protocols::moves;

static basic::Tracer TR( "apps.public.mp_domain_assembly" );

class MPDomainAssembly : public Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Mover
	MPDomainAssembly();

	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MPDomainAssembly( MPDomainAssembly const & src );

	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MPDomainAssembly &
	operator=( MPDomainAssembly const & src );

	/// @brief Destructor
	~MPDomainAssembly();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	virtual std::string get_name() const;

	/// @brief Get Membrane protein interface statistics
	virtual void apply( Pose & pose );

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	// /// @brief Create a Clone of this mover
	// virtual protocols::moves::MoverOP clone() const;
	//
	// /// @brief Create a Fresh Instance of this Mover
	// virtual protocols::moves::MoverOP fresh_instance() const;
	//
	// /// @brief Pase Rosetta Scripts Options for this Mover
	// void parse_my_tag(
	//       utility::tag::TagCOP tag,
	//       basic::datacache::DataMap &,
	//       protocols::filters::Filters_map const &,
	//       protocols::moves::Movers_map const &,
	//       core::pose::Pose const &
	//       );

private: // methods

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();

	// create residue
	ResidueOP create_residue_from_olc( char olc );

	// create residue
	Residue create_residue_from_resn( Pose & pose, core::Size resnumber );

	// set random torsion
	void set_random_torsion( Pose & pose, Size resn );

	/// @brief Print score to cout
	void print_score( Pose & pose, core::scoring::ScoreFunctionOP sfxn );


private: // data

	/// @brief Fragment filename
	std::string frag9_fn_;

	/// @brief Fragment filename
	std::string frag3_fn_;

	/// @brief Sequence of the full-length protein: given as fasta
	std::string full_seq_;

	/// @brief Number of the pose in the -in:file:l list that is the TM domain
	core::Size tmpdb_;

	/// @brief Vector of pose input files
	utility::vector1< std::string > infiles_;

	/// @brief Vector of poses_ that are read in from PDB files
	utility::vector1< core::pose::Pose > poses_;

	/// @brief Vector of sequences_ in poses_
	utility::vector1< std::string > sequences_;

	/// @brief Vector of offsets_, i.e. residue numbers along the full-length sequence
	///   for which the PDB poses_ start
	utility::vector1< core::Size > offsets_;

	/// @brief Vector of length (nres of full-length pose), describing whether the
	///   residue exists in a pose (0 = no, 1 = yes, >1 = overlap)
	utility::vector1< core::Size > in_pdb_;

	/// @brief Vectors of linker sequences_ to model
	utility::vector1< std::string > linkers_;

	/// @brief Movemap of which residues can move (i.e. are in linkers_ +/-1 residue)
	core::kinematics::MoveMap movemap_;

	/// @brief scorefunction (i.e. membrane sfxn_)
	core::scoring::ScoreFunctionOP sfxn_;

};

//////////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MPDomainAssembly::MPDomainAssembly() :
	Mover(),
	frag9_fn_(),
	frag3_fn_(),
	full_seq_(),
	tmpdb_(),
	infiles_(),
	poses_(),
	sequences_(),
	offsets_(),
	in_pdb_(),
	linkers_(),
	movemap_(),
	sfxn_()
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MPDomainAssembly::MPDomainAssembly( MPDomainAssembly const & src ) :
	Mover( src ),
	frag9_fn_( src.frag9_fn_ ),
	frag3_fn_( src.frag9_fn_ ),
	full_seq_( src.full_seq_ ),
	tmpdb_( src.tmpdb_ ),
	infiles_( src.infiles_ ),
	poses_( src.poses_ ),
	sequences_( src.sequences_ ),
	offsets_( src.offsets_ ),
	in_pdb_( src.in_pdb_ ),
	linkers_( src.linkers_ ),
	movemap_( src.movemap_ ),
	sfxn_( src.sfxn_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MPDomainAssembly &
MPDomainAssembly::operator=( MPDomainAssembly const & src )
{
	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new MPDomainAssembly( *this ) );
}

/// @brief Destructor
MPDomainAssembly::~MPDomainAssembly() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MPDomainAssembly::get_name() const {
	return "MPDomainAssembly";
}

/// @brief run domain assembly protocol for membrane proteins
void MPDomainAssembly::apply( Pose & pose ) {

	using namespace protocols::membrane;

	TR << "Calling MPDomainAssembly" << std::endl;

	register_options();
	init_from_cmd();

	// clear vector contents, needed for running things in parallel (not sure why)
	pose.clear();
	poses_.clear();
	sequences_.clear();
	offsets_.clear();

	////////////////////////////////////////////////////////////////////////////////////////
	/// SETUP

	// read in poses_
	utility::vector1< core::pose::Pose > poses_ = core::import_pose::poses_from_files( infiles_ );
	Pose tmdomain = poses_[ tmpdb_ ];

	// setup poses_
	TR << "switching residue typeset and removing cutpoint variants" << std::endl;
	for ( core::Size i = 1; i <= poses_.size(); ++i ) {

		// switch poses_ into full-atom
		protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap( core::chemical::FA_STANDARD );
		typeset_swap.apply( poses_[i] );

		// remove cutpoint variants for all subposes_ except for N- and C-terminus
		if ( i != 1 ) {
			remove_lower_terminus_type_from_pose_residue( poses_[i], 1 );
		}
		if ( i != poses_.size() ) {
			remove_upper_terminus_type_from_pose_residue( poses_[i], poses_[i].total_residue() );
		}
	}

	// get sequences_ from poses_
	for ( core::Size i = 1; i <= poses_.size(); ++i ) {
		sequences_.push_back( poses_[i].sequence() );
	}

	TR << "POSES " << poses_.size() << std::endl;
	TR << "SEQUENCES " << sequences_.size() << std::endl;


	////////////////////////////////////////////////////////////////////////////////////////
	/// OFFSETS

	// check for gaps in poses_: error for mutations or missing density
	// go through all pose sequences_ (ATOM lines)
	for ( core::Size i = 1; i <= sequences_.size(); ++i ) {

		std::string this_seq = sequences_[i];
		core::Size max_matches = 0;
		core::Size offset = 0;
		TR << "sequence: " << sequences_[i] << std::endl;

		// go through the sequence itself
		for ( core::Size j = 0; j <= this_seq.size(); ++j ) {

			// go through the full-length sequence
			for ( core::Size k = 0; k <= full_seq_.size(); ++k ) {

				core::Size matches = 0;

				// if characters match, go through the rest of the sequence using an offset
				if ( full_seq_[k] == this_seq[j] ) {

					for ( core::Size m = 0; m <= this_seq.size() - j; ++m ) {
						if ( k+m <= full_seq_.size() && full_seq_[ k+m ] == this_seq[ j+m ] ) {
							matches++;
						}
					}

					// replace the maximum number of matches
					if ( matches > max_matches ) {
						max_matches = matches;
						offset = k;
					}
				}
			}

			// if the number of matches for each sequence is shorter than the length of the sequence itself,
			// then there must either be a mutation or missing density
			if ( max_matches < this_seq.length() ) {
				TR << "max_matches: " << max_matches << ", this_seq.size(): " << this_seq.size() << std::endl;
				utility_exit_with_message( "The sequence in one of the PDBs doesn't match the sequence in the fasta file! There could either be a mutation or missing density. Please close loops within PDB files before running this protocol (i.e. connecting the pieces with linkers)! The sequence in question is " + sequences_[i] + "\n" );
			}
		}

		offsets_.push_back( offset );
	}

	// make sure that the offsets_ are in ascending order, i.e. the provided poses_ are
	// given in order from the N-terminus to the C-terminus
	for ( core::Size i = 2; i <= offsets_.size(); ++i ) {
		if ( offsets_[i] < offsets_[ i-1 ] ) {

			for ( core::Size k = 1; k <= offsets_.size(); ++k ) {
				TR << "offsets " << k << " " << offsets_[k] << std::endl;
			}
			utility_exit_with_message( "The provided PDBs are not in order: they should be in order from the N-terminus to the C-terminus. Quitting!" );
		}
	}

	TR << "offsets_:" << std::endl;
	for ( core::Size i = 1; i <= offsets_.size(); ++i ) {
		TR << offsets_[i] << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	/// LINKERS

	// figure out which linkers_ need to be modeled
	// create a vector of booleans: 0 means 'residue not in any PDBs', 1
	// means 'needs to be modeled', and >1 means overlap between different poses_ and error message
	utility::vector1< core::Size > in_pdb_( full_seq_.size(), 0 );

	// go through domains and set boolean to true, if residue exists in any PDB
	for ( core::Size i = 1; i <= sequences_.size(); ++i ) {

		// get offset for that sequence
		core::Size this_offset = offsets_[i];

		for ( core::Size j = 1 + this_offset; j < 1 + this_offset + sequences_[i].size(); ++j ) {
			in_pdb_[ j ] += 1;
		}
	}

	// if in_pdb_ has a value >1, then there is overlap between poses_, crash
	if ( in_pdb_[ static_cast< core::Size > ( utility::argmax( in_pdb_ ) ) ] > 1 ) {
		utility_exit_with_message( "There is overlap between the sequences_ to model. Can't determine which coordinates to use. Please make sure there is always at least a 1-residue linker in between PDBs. Quitting!" );
	}

	// create poses_ for linkers_, use pose_from_sequence
	// vectors full_seq_vec and in_pdb_ correspond to each other
	utility::vector1< std::string > linkers_;
	std::string linker = "";
	for ( core::Size i = 1; i <= in_pdb_.size(); ++i ) {

		if ( in_pdb_[i] == 0 ) {

			// start of new linker, pushback to vector of linkers_
			if ( i > 1 && in_pdb_[ i-1 ] == 1 && linker.size() > 0 ) {
				linkers_.push_back( linker );
				linker = "";
			}

			// string is zero-indexed!
			linker += full_seq_[ i-1 ];
		}
	}

	// pushback last linker
	linkers_.push_back( linker );

	// do we have an N-terminal linker?
	std::string Nterm;
	if ( in_pdb_[1] == 0 ) {
		Nterm = linkers_[1];
		linkers_.pop( Nterm );
	}

	TR << "Nterminal sequence? " << Nterm << std::endl;
	TR << "linker vector: " << std::endl;
	for ( core::Size i = 1; i <= linkers_.size(); ++i ) {
		TR << i << " " << linkers_[i] << std::endl;
	}

	////////////////////////////////////////////////////////////////////////////////////////
	/// RUN DOMAIN ASSEMBLY STARTING WITH TM DOMAIN

	// create a pose from the TM domain
	Pose full_pose( poses_[ tmpdb_ ] );

	TR << "walking C-terminally..." << std::endl;

	// start with TM domain, then add residues C-terminally:
	// linker, then next pose,linker,pose,linker etc
	for ( core::Size i = tmpdb_; i <= poses_.size(); ++i ) {

		// adding pose
		if ( i != tmpdb_ ) {

			TR << "adding pose " << i << " " << poses_[i].sequence() << std::endl;

			// add two residues at the beginning
			Residue rsd( create_residue_from_resn( poses_[i], 1 ) );
			full_pose.conformation().safely_append_polymer_residue_after_seqpos( rsd, full_pose.total_residue(), true );

			Residue rsd1( create_residue_from_resn( poses_[i], 2 ) );
			full_pose.conformation().safely_append_polymer_residue_after_seqpos( rsd1, full_pose.total_residue(), true );

			// superimpose first two residues in ref pose with first two in full pose
			// superimpose requires two residues!!!
			SuperimposeMoverOP align( new SuperimposeMover( full_pose, full_pose.total_residue()-1, full_pose.total_residue(), 1, 2, false ) );
			align->apply( poses_[i] );

			// remove the added residues
			full_pose.conformation().delete_polymer_residue( full_pose.total_residue() );
			full_pose.conformation().delete_polymer_residue( full_pose.total_residue() );

			// append pose to pose
			append_pose_to_pose( full_pose, poses_[i], false );

		}

		// only continue if we have a linker C-terminally of the TM domain
		if ( i > linkers_.size() ) {
			continue;
		}

		// go through linker
		for ( core::Size j = 1; j <= linkers_[i].size(); ++j ) {

			TR << "adding linker " << i << " " << linkers_[i] << std::endl;

			ResidueOP rsd( create_residue_from_olc( static_cast< char >( linkers_[i][j] ) ) );
			full_pose.conformation().safely_append_polymer_residue_after_seqpos( *rsd, full_pose.total_residue(), true );
			set_random_torsion( full_pose, full_pose.total_residue() );

		}
	}

	TR << "walking N-terminally..." << std::endl;

	// now add residues N-terminally:
	// linker, then previous pose,linker,pose,linker etc
	for ( core::Size i = tmpdb_-1; i >= 1; --i ) {

		TR << "adding linker " << i << " " << linkers_[i] << std::endl;

		// go through linker
		for ( core::Size j = linkers_[i].size(); j >= 1; --j ) {

			ResidueOP rsd( create_residue_from_olc( static_cast< char >( linkers_[i][j-1] ) ) );
			full_pose.conformation().safely_prepend_polymer_residue_before_seqpos( *rsd, 1, true );
			set_random_torsion( full_pose, 2 );

		}

		// adding pose
		TR << "adding pose " << i << " " << poses_[i].sequence() << std::endl;

		// add two residues at end
		Residue rsd( create_residue_from_resn( poses_[i], poses_[i].total_residue() ) );
		full_pose.conformation().safely_prepend_polymer_residue_before_seqpos( rsd, 1, true );

		Residue rsd1( create_residue_from_resn( poses_[i], poses_[i].total_residue()-1 ) );
		full_pose.conformation().safely_prepend_polymer_residue_before_seqpos( rsd1, 1, true );

		// superimpose last two residues in ref pose with first two in full pose
		// superimpose requires two residues!!!
		SuperimposeMoverOP align( new SuperimposeMover( full_pose, 1, 2, poses_[i].total_residue()-1, poses_[i].total_residue(), false ) );
		align->apply( poses_[i] );

		// remove the added residues
		full_pose.conformation().delete_polymer_residue( 1 );
		full_pose.conformation().delete_polymer_residue( 1 );

		// append pose to pose
		Pose dummy( poses_[i] );
		append_pose_to_pose( dummy, full_pose, false );
		full_pose = *dummy.clone();

	}

	// superimpose full pose to TM domain again
	SuperimposeMoverOP align( new SuperimposeMover( tmdomain, 1, tmdomain.total_residue(), offsets_[tmpdb_]+1, offsets_[tmpdb_]+tmdomain.total_residue(), false ) );
	align->apply( full_pose );

	// full_pose.dump_pdb( "full_pose1.pdb" );

	////////////////////////////////////////////////////////////////////////////////////////
	/// CREATE MOVEMAP

	// create a MoveMap and set the backbone moveable to the linkers_ plus 1 flanking residue on either side
	// the +1 is for the membrane residue
	TR << "finding movable residues" << std::endl;
	utility::vector1< bool > move_bb( in_pdb_.size()+1, false );
	for ( core::Size i = 1; i <= in_pdb_.size(); ++i ) {

		// only flexible linkers_ will be movable
		if ( in_pdb_[ i ] == 0 ) {
			move_bb[ i ] = true;

			// and the residues before and after
			if ( i > 0 ) {
				move_bb[ i-1 ] = true;
			}
			if ( i < in_pdb_.size() ) {
				move_bb[ i+1 ] = true;
			}
		}
	}

	// finally create the MoveMap
	TR << "and creating a MoveMap" << std::endl;
	core::kinematics::MoveMapOP movemap_( new core::kinematics::MoveMap() );
	movemap_->set_bb( move_bb );
	//  movemap_->show();

	////////////////////////////////////////////////////////////////////////////////////////
	/// SPANNING TOPOLOGY

	// add the membrane to the pose but create SpanningTopology from TM pose first
	TR << "getting chain and z" << std::endl;
	std::pair< utility::vector1< core::Real >, utility::vector1< core::Real > > chain_and_z( get_chain_and_z( tmdomain ) );
	utility::vector1< char > secstruct( get_secstruct( tmdomain ) );
	core::Real thickness( 15 );

	// create SpanningTopology
	TR << "topology:" << std::endl;
	SpanningTopologyOP topo( new SpanningTopology( chain_and_z.first, chain_and_z.second, secstruct, thickness ) );
	TR << "topology before shifting: " << std::endl;
	topo->show( TR );
	SpanningTopologyOP topo_tm( new SpanningTopology( chain_and_z.first, chain_and_z.second, secstruct, thickness ) );

	// shift the spans according to their position in the full-length pose
	for ( core::Size i = 1; i <= topo->nspans(); ++i ) {

		topo->span( i )->shift( offsets_[ tmpdb_ ] );

	}
	TR << "topology after shifting: " << std::endl;
	topo->show( TR );

	// write out spanfile
	TR << "writing spanfile" << std::endl;
	topo->write_spanfile( "domain_assembly.span" );

	////////////////////////////////////////////////////////////////////////////////////////
	/// MEMBRANE, FOLDTREE, SCOREFUNCTION

	// find anchor point from TM domain, then shift in full pose
	AddMembraneMoverOP addmem( new AddMembraneMover( topo_tm ) );
	addmem->apply( tmdomain );
	TR << "nres: " << tmdomain.total_residue() << " " << nres_protein( tmdomain ) << std::endl;
	TR << "finding anchor point from TM domain" << std::endl;
	core::Size anchor_point = create_membrane_foldtree_anchor_pose_tmcom( tmdomain );

	// add the membrane
	core::Size anchor_full = anchor_point + offsets_[ tmpdb_ ];
	TR << "adding membrane at residue " << anchor_full << std::endl;
	AddMembraneMoverOP addmem1( new AddMembraneMover( topo, anchor_full, 0 ) );
	addmem1->apply( full_pose );

	// find anchor residue: this should be TM COM but shifted according to offset of the TM domain
	TR << "======= setting foldtree" << std::endl;
	core::Size anchor_point1 = create_membrane_foldtree_anchor_pose_tmcom( full_pose );
	TR << "anchor point for TM pose " << anchor_point1 << std::endl;
	full_pose.fold_tree().show( TR );

	// create scorefxn
	TR << "sfxn_" << std::endl;
	using namespace core::scoring;
	ScoreFunctionOP sfxn_ = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

	// set chains in PDBInfo
	PDBInfoOP new_info( new PDBInfo( full_pose ) );
	new_info->set_chains( 'A' );
	full_pose.pdb_info( new_info );

	// full_pose.dump_pdb( "full_pose2.pdb" );

	// create MC object
	using namespace protocols::moves;
	core::Real kT = 2.0;
	MonteCarloOP mc( new MonteCarlo( full_pose, *sfxn_, kT ) );

	// getting fragment set
	using namespace core::fragment;
	FragSetOP frags9 = FragmentIO().read_data( frag9_fn_ );
	FragSetOP frags3 = FragmentIO().read_data( frag3_fn_ );

	// setup fragment mover
	ClassicFragmentMoverOP frag9_mover( new ClassicFragmentMover( frags9, movemap_ ) );
	ClassicFragmentMoverOP frag3_mover( new ClassicFragmentMover( frags3, movemap_ ) );

	// figure out number of fragment insertions we want to make
	// depends on the number of residues in linkers
	core::Size n_linker_res( 0 );
	for ( core::Size i = 1; i <= in_pdb_.size(); ++i ) {
		if ( in_pdb_[i] == 0 ) {
			n_linker_res++;
		}
	}

	// initial fragment insertion cycle
	core::Size iterations( 1 );
	core::Size stage1_cycles( n_linker_res * 8 );

	// Small and ShearMover for linkers_
	core::Size nmoves( full_pose.total_residue() );
	SmallMoverOP small( new SmallMover( movemap_, kT, nmoves ) );
	ShearMoverOP shear( new ShearMover( movemap_, kT, nmoves ) );
	small->angle_max( 180.0 );
	shear->angle_max( 180.0 );

	// initialize AtomTreeMinimizer
	// tried CartesianMinimizer, but it made breaks into the pose
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	min_opts.max_iter( 2000 );
	core::optimization::AtomTreeMinimizer atm;

	// create packer task - will be re-used
	using namespace core::pack::task;
	PackerTaskOP repack = TaskFactory::create_packer_task( full_pose );
	repack->restrict_to_residues( move_bb );
	repack->restrict_to_repacking();

	for ( core::Size j = 1; j <= iterations; ++j ) {

		TR << "Fragment insertion for " << n_linker_res << " linker residues in " << stage1_cycles << " cycles" << std::endl;
		for ( Size i = 1; i <= stage1_cycles; ++i ) {

			// 9-residue fragments
			frag9_mover->apply( full_pose );
			( *sfxn_ )( full_pose );
			mc->boltzmann( full_pose , "frag9" );

			// 3-residue fragments
			frag3_mover->apply( full_pose );
			( *sfxn_ )( full_pose );
			mc->boltzmann( full_pose , "frag3" );

			mc->recover_low( full_pose );

			// Small and ShearMover for linkers_
			TR << "small and shearmover" << std::endl;
			small->apply( full_pose );
			core::pack::pack_rotamers( full_pose, *sfxn_, repack );
			( *sfxn_ )( full_pose );
			mc->boltzmann( full_pose , "small_pack" );
			mc->recover_low( full_pose );
			print_score( full_pose, sfxn_ );

			shear->apply( full_pose );
			core::pack::pack_rotamers( full_pose, *sfxn_, repack );
			( *sfxn_ )( full_pose );
			mc->boltzmann( full_pose , "shear_pack" );
			mc->recover_low( full_pose );
			print_score( full_pose, sfxn_ );

			if ( i % 100 == 0 ) {
				mc->show_scores();
				mc->show_counters();
			}
		}

		// Small and ShearMover for linkers_
		TR << "FINAL small and shearmover" << std::endl;
		small->apply( full_pose );
		print_score( full_pose, sfxn_ );
		mc->boltzmann( full_pose , "small" );
		mc->recover_low( full_pose );

		shear->apply( full_pose );
		print_score( full_pose, sfxn_ );
		mc->boltzmann( full_pose , "shear" );
		mc->recover_low( full_pose );

		// pack rotamers
		core::pack::pack_rotamers( full_pose, *sfxn_, repack );
		( *sfxn_ )( full_pose );
		mc->boltzmann( full_pose , "pack" );
		mc->recover_low( full_pose );
		//  full_pose.dump_pdb( "full_pose3.pdb" );

		//  full_pose.dump_pdb( "full_pose4.pdb" );

		// run AtomTreeMinimizer
		TR << "AtomTree Minimization" << std::endl;
		atm.run( full_pose, *movemap_, *sfxn_, min_opts );

		print_score( full_pose, sfxn_ );
		mc->boltzmann( full_pose , "AT_min" );
		mc->recover_low( full_pose );

	}

	// another MPRangeRelax step if angle max flag is given
	if ( option[relax::range::angle_max].user() ) {

		using namespace protocols::relax::membrane;
		MPRangeRelaxMoverOP mprr( new MPRangeRelaxMover() );
		mprr->optimize_membrane( false );
		mprr->apply( full_pose );

	}

	// do a final alignment onto the TM domain
	align->apply( full_pose );

	pose = full_pose;

} // apply

/// @brief Register Options with JD2
void MPDomainAssembly::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::fasta );
	option.add_relevant( OptionKeys::mp::assembly::poses );
	option.add_relevant( OptionKeys::in::file::frag3 );
	option.add_relevant( OptionKeys::in::file::frag9 );
	option.add_relevant( OptionKeys::mp::assembly::TM_pose_number );
	option.add_relevant( OptionKeys::relax::range::angle_max );

} // register options

/// @brief Initialize Mover options from the comandline
void MPDomainAssembly::init_from_cmd() {

	using namespace basic::options;

	// cry if PDB list not given
	if ( option[OptionKeys::in::file::fasta].user() ) {
		full_seq_ = core::sequence::read_fasta_file( option[ OptionKeys::in::file::fasta ]()[1] )[1]->sequence();
		TR << "Read in fasta file " << option[OptionKeys::in::file::fasta]()[1] << std::endl;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide fasta file with -in:file:fasta!");
	}

	// read in PDB list
	if ( option[mp::assembly::poses].user() ) {
		infiles_ = basic::options::option[mp::assembly::poses]();
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide a list of PDB files with -in:file:l!");
	}

	// figure out which of the poses_ is the TM span
	if ( option[OptionKeys::mp::assembly::TM_pose_number].user() ) {
		tmpdb_ = option[OptionKeys::mp::assembly::TM_pose_number]();
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Which number of the PDBs in the file list is located in the membrane?");
	}

	// read in fragments
	if ( option[OptionKeys::in::file::frag3].user() ) {
		frag3_fn_ = option[OptionKeys::in::file::frag3]();
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide fragments with -in:file:frag3!");
	}
	if ( option[OptionKeys::in::file::frag9].user() ) {
		frag9_fn_ = option[OptionKeys::in::file::frag9]();
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please provide fragments with -in:file:frag9!");
	}

} // init from commandline

////////////////////////////////////////////////////////////////////////////////
// create residue
Residue MPDomainAssembly::create_residue_from_resn( Pose & pose, core::Size resnumber ) {

	using namespace core::conformation;

	std::string name3 = pose.residue( resnumber ).name3();

	ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	ResidueType const & rsd_type( residue_set->name_map( name3 ) );

	// true: preserve Cbeta
	ResidueOP rsd( ResidueFactory::create_residue( rsd_type, pose.residue( resnumber ), pose.conformation(), true ) );

	Residue rsd1( *rsd );

	return rsd1;

} // create residue from tlc

////////////////////////////////////////////////////////////////////////////////
// create residue
ResidueOP MPDomainAssembly::create_residue_from_olc( char olc ) {

	using namespace core::chemical;
	ResidueTypeSetCOP const & residue_set( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	ResidueTypeCOP rsd_type( residue_set->get_representative_type_name1( olc ) );
	ResidueType const & res( *rsd_type );
	ResidueOP rsd( ResidueFactory::create_residue( res ) );

	return rsd;
} // create residue from OLC

////////////////////////////////////////////////////////////////////////////////
// set random torsion
void MPDomainAssembly::set_random_torsion( Pose & pose, Size resn ) {

	// get random phi and psi
	core::Real random_phi( numeric::random::uniform() * 360 - 180 );
	core::Real random_psi( numeric::random::uniform() * 360 - 180 );

	// set phi and psi to random and omega to 180
	pose.set_phi( resn, random_phi );
	pose.set_phi( resn, random_psi );
	pose.set_omega( resn, 180 );

} // set random torsion

////////////////////////////////////////////////////////////////////////////////
/// @brief Print score to cout
void MPDomainAssembly::print_score( Pose & pose, core::scoring::ScoreFunctionOP sfxn ) {

	// print energies and iteration
	Real tot_score = ( *sfxn )( pose );
	Real fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
	Real fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];

	pose.energies().show_total_headers( TR );
	TR << std::endl;
	pose.energies().show_totals( TR );
	TR << std::endl;

	TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

} // print score

/////////////////////////////////
///// Rosetta Scripts Methods ///
/////////////////////////////////
//
///// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//MPDomainAssembly::clone() const {
// return ( protocols::moves::MoverOP( new MPDomainAssembly( *this ) ) );
//}
//
///// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//MPDomainAssembly::fresh_instance() const {
// return protocols::moves::MoverOP( new MPDomainAssembly() );
//}
//
///// @brief Pase Rosetta Scripts Options for this Mover
//void
//MPDomainAssembly::parse_my_tag(
//            utility::tag::TagCOP tag,
//            basic::datacache::DataMap &,
//            protocols::filters::Filters_map const &,
//            protocols::moves::Movers_map const &,
//            core::pose::Pose const &
//            ) {
//  // TODO
//
//}
//
///// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//MPDomainAssemblyCreator::create_mover() const {
// return protocols::moves::MoverOP( new MPDomainAssembly );
//}
//
///// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//MPDomainAssemblyCreator::keyname() const {
// return MPDomainAssemblyCreator::mover_name();
//}
//
///// @brief Mover name for Rosetta Scripts
//std::string
//MPDomainAssemblyCreator::mover_name() {
// return "MPDomainAssembly";
//}

//////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< MPDomainAssembly > MPDomainAssemblyOP;

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MAIN ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{
	try {

		// initialize option system, RNG, and all factory-registrators
		devel::init(argc, argv);

		//  protocols::jd2::register_options();

		// Create and kick off a new load membrane mover
		core::pose::Pose pose;
		MPDomainAssemblyOP mp_domain_assembly( new MPDomainAssembly() );
		protocols::jd2::JobDistributor::get_instance()->go( mp_domain_assembly );

		return 0;

	}
catch (utility::excn::Exception const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
