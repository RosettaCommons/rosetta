// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/NucleotideMutation.cc
/// @brief Substitutes amino acids in pose based on mutations in nucleotide space
/// @details Translate amino acid sequence to nucleotide sequence, introduce
/// random mutation, and affectuates the amino acid substitution in pose. Stop
/// codons are avoided and with default setting nucleotide mutations will continue
/// untill a non-silent mutation occurs.
/// @author Christoffer Norn (ch.norn@gmail.com)

// Unit headers
#include <protocols/protein_interface_design/movers/NucleotideMutation.hh>
#include <protocols/protein_interface_design/movers/NucleotideMutationCreator.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/NucleotideTools.hh>

#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#include <sstream>
#include <core/chemical/AA.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>

#include <utility/vector0.hh>
#include <map>
#include <math.h>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.NucleotideMutation" );

std::string
NucleotideMutationCreator::keyname() const
{
	return NucleotideMutationCreator::mover_name();
}

protocols::moves::MoverOP
NucleotideMutationCreator::create_mover() const {
	return protocols::moves::MoverOP( new NucleotideMutation );
}

std::string
NucleotideMutationCreator::mover_name()
{
	return "NucleotideMutation";
}

NucleotideMutation::NucleotideMutation() :
	Mover( NucleotideMutationCreator::mover_name() ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	init_sequence_(""),
	continue_if_silent_( true ),
	reference_pose_( /* NULL */ )
{
}


NucleotideMutation::~NucleotideMutation() {}

void
NucleotideMutation::add_nt_seq_to_pose( core::pose::Pose & pose ){
	using namespace core::chemical;

	string aa_sequence = pose.sequence();
	string nt_sequence;
	if ( init_sequence() != "" ) {
		TR << "Initializing nucleotide sequence with sequence " << std::endl;
		TR << init_sequence() << std::endl;
		runtime_assert_string_msg( aa_sequence.length()*3 == init_sequence().length(), "Your nucleotide sequence does not have the same length as your protein sequence * 3" );
		core::pose::add_comment(pose, "nt_seq", init_sequence());
	} else {
		TR << "Initializing nucleotide sequence with random codons " << std::endl;
		for ( size_t i = 0; i < aa_sequence.length(); i++ ) {
			string nt = NucleotideTools::aa2randomCodon( aa_sequence[i] );
			nt_sequence.append( nt );
		}
		core::pose::add_comment(pose, "nt_seq", nt_sequence);
	}
}

void
NucleotideMutation::apply( core::pose::Pose & pose )
{

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;


	/////////////////////////////////////////////////////////////////
	// translate aa to nt if the protein doesn't already have a nt //
	// sequence in comments                                        //
	/////////////////////////////////////////////////////////////////
	std::map<std::string, std::string> comments = core::pose::get_all_comments( pose );
	bool has_nt_sequence = ( comments.find( "nt_seq" ) != comments.end() );
	if ( !has_nt_sequence ) {
		add_nt_seq_to_pose( pose );
	} else {
		TR << "Nucleotide sequence initialized from comments " << std::endl;
	}

	/////////////////////////////////////////////////////////////////
	// find out which residues are designable and pick one randomly//
	/////////////////////////////////////////////////////////////////
	PackerTaskCOP task;
	if ( cache_task_ && task_ ) {
		if ( pose.total_residue() == task_->total_residue() ) {
			task = task_;
		} else {
			task_.reset(); // Invalidate cached task.
		}
	}
	if ( ! task ) {
		task = task_factory()->create_task_and_apply_taskoperations( pose );
	}
	if ( cache_task_ && !task_ ) {
		task_ = task;
	}

	utility::vector1< core::Size > being_designed;
	being_designed.clear();
	for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
		if ( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() ) {
			being_designed.push_back( resi );
		}
	}
	if ( being_designed.empty() ) {
		TR.Warning << "WARNING: No residues are listed as designable." << std::endl;
		return;
	}


	//////////////////////////////////////////////////////////////////////////////////////
	// Make mutation in nucleotide space, ignore stop codons, continue if silent mutation
	//////////////////////////////////////////////////////////////////////////////////////
	bool cont = true;
	while ( cont ) {
		core::Size const target_aa_no = being_designed[ (core::Size) floor( numeric::random::rg().uniform() * being_designed.size() ) + 1 ];

		string nt_sequence = core::pose::get_all_comments(pose)[ "nt_seq" ];
		string target_codon = nt_sequence.substr( ( target_aa_no - 1) * 3, 3); // substrating by 1 here as aa seq starts from 1.

		TR << "nt seq before mut " << nt_sequence << std::endl;
		TR << "aa seq before mut " << pose.sequence() << std::endl;
		TR << "mutation will happen at amino acid " << pose.sequence()[ target_aa_no - 1 ] << target_aa_no << " with codon " << target_codon << std::endl;

		core::Size nt_no_in_codon = numeric::random::random_range(0,2);
		//string target_nt(1, target_codon[nt_no_in_codon]);
		char target_nt = target_codon[nt_no_in_codon];

		// Make mutation according to the K80 model with kappa=2, distance=1
		core::Real kappa = 2.0;
		core::Real d = 1.0;
		// Transition probability matrix for K80 with, d=1 and kappa=2
		//   A  T  C  G
		// A p0 p1 p2 p2
		// T p1 p0 p2 p2
		// C p2 p2 p0 p1
		// G p2 p2 p1 p0
		//
		// where
		core::Real p0 = 0.25 + 0.25 * exp(-4.0 * d/(kappa+2.0)) + 0.5 * exp(-2.0*(kappa+1.0)/(kappa+2.0));
		core::Real p1 = 0.25 + 0.25 * exp(-4.0 * d/(kappa+2.0)) - 0.5 * exp(-2.0*(kappa+1.0)/(kappa+2.0));
		core::Real p2 = 0.25 - 0.25 * exp(-4.0 * d/(kappa+2.0));

		// if kappa=2 and d=1 this works out to be
		//         A              T              C              G
		// A 0.453534940367 0.230404780219 0.158030139707 0.158030139707
		// T 0.230404780219 0.453534940367 0.158030139707 0.158030139707
		// C 0.158030139707 0.158030139707 0.453534940367 0.230404780219
		// G 0.158030139707 0.158030139707 0.230404780219 0.453534940367
		//
		// The non-changing nucleotide mutations will only add a constant
		// scaling to the branch time, equal for all nts. We shouldn't
		// care about this. To normalize we do:

		core::Real norm = 1/(1-p0);

		// and multiply this on p1 and p2. For kappa=2 and d=1 the new table is
		//
		//          A              T              C              G
		// A        0       0.421627652   0.289186174   0.289186174
		// T  0.421627652          0      0.289186174   0.289186174
		// C  0.289186174   0.289186174          0      0.421627652
		// G  0.289186174   0.289186174   0.421627652        0
		//
		// the cummulative probabilities are then
		//          A
		// A        0
		// T  0.421627652
		// C  0.710813826       ect ...
		// G        1

		core::Real cum_p_transition = p1*norm;
		core::Real cum_p_transversion_type_1 = p1*norm+p2*norm;

		// I'll use these probabilities instead of the transition
		// probabilities, so it is K80, but speeded up.

		core::Real const random_num = numeric::random::rg().uniform();

		std::string mut_nt;
		//TR << "mutating nt position " << nt_no_in_codon << " in codon "<< std::endl;
		if ( random_num < cum_p_transition /* 0.421627652 for d=1 k=2 */ ) {
			//TR << "making transition mutation " << std::endl;
			if      ( target_nt == 'T' ) { mut_nt = 'A'; }
			else if ( target_nt == 'A' ) { mut_nt = 'T'; }
			else if ( target_nt == 'G' ) { mut_nt = 'C'; }
			else if ( target_nt == 'C' ) { mut_nt = 'G'; }
		} else if ( random_num < cum_p_transversion_type_1 /* 0.710813826 for d=1 k=2 */ ) {
			//TR << "making transversion type 1 mutation " << std::endl;
			if      ( target_nt == 'T' ) { mut_nt = 'C'; }
			else if ( target_nt == 'A' ) { mut_nt = 'C'; }
			else if ( target_nt == 'G' ) { mut_nt = 'A'; }
			else if ( target_nt == 'C' ) { mut_nt = 'A'; }
		} else {
			//TR << "making transversion type 2 mutation " << std::endl;
			if      ( target_nt == 'T' ) { mut_nt = 'G'; }
			else if ( target_nt == 'A' ) { mut_nt = 'G'; }
			else if ( target_nt == 'G' ) { mut_nt = 'T'; }
			else if ( target_nt == 'C' ) { mut_nt = 'T'; }
		}

		string old_aa = NucleotideTools::codon2aa(target_codon);
		TR << "old codon was " << target_codon << " encoding " << old_aa << std::endl;

		string new_codon = target_codon.replace(nt_no_in_codon, 1, mut_nt);
		string new_aa = NucleotideTools::codon2aa(new_codon);

		TR << "new codon is " << new_codon << " encoding " << new_aa << std::endl;

		ostringstream mutation;
		mutation << old_aa << target_aa_no << new_aa;

		string new_aa_seq = pose.sequence().replace( target_aa_no - 1, 1, new_aa);
		core::pose::add_comment(pose, "aa_seq", new_aa_seq);

		if ( (new_aa != "Stop") && ( new_aa == old_aa ) ) {
			nt_sequence.replace( ( target_aa_no - 1 ) * 3, 3, new_codon);
			core::pose::add_comment(pose, "nt_seq", nt_sequence);
			core::pose::add_comment(pose, "last_nt_mut_was_non_silent", "false");
			core::pose::add_comment(pose, "mutation", mutation.str());
			if ( continue_if_silent() == true ) {
				TR << "Made a silent mutation and will attempt another mutation." << std::endl;
				cont = true;
			} else {
				TR << "Made a silent mutation" << std::endl;
				cont = false;
			}
		} else if ( new_aa != "Stop" ) { // then make the mutation
			nt_sequence.replace( ( target_aa_no - 1 ) * 3, 3, new_codon);
			core::pose::add_comment(pose, "nt_seq", nt_sequence);
			core::pose::add_comment(pose, "last_nt_mut_was_non_silent", "true");
			core::pose::add_comment(pose, "mutation", mutation.str());
			TR << "Made a mutation resulting in an amino acid substitution." << std::endl;
			cont = false;

			std::map<std::string, std::string> final_comments = core::pose::get_all_comments( pose );

			// It is important to reset the pose to the reference pose here,
			// instead of carrying on the minimized pose... otherwise the pose
			// get more and more minimized throughout the trajectory. That
			// would cause a bunch of wierd artifacts including locked
			// trajectories, where mutations can't be accepted for long streches...
			// Also, in evolution, the new sequence does not remember
			// the conformation of ancestral structure (it folds anew).
			pose = *reference_pose_;

			// First add back the comments
			for ( std::map<string, string>::const_iterator i = final_comments.begin(); i != final_comments.end(); ++i ) {
				std::string const key(i->first);
				std::string const val(i->second);
				core::pose::add_comment(pose, key, val);
			}

			// Now thread on the sequence we have in comments.
			using namespace protocols::toolbox::task_operations;
			using namespace core::pack::task;
			using namespace core::pack::task::operation;

			ThreadSequenceOperationOP tso( new ThreadSequenceOperation );
			tso->target_sequence( final_comments["aa_seq"] );
			tso->allow_design_around(false);

			TaskFactoryOP tf;
			tf = TaskFactoryOP( new TaskFactory );
			tf->push_back(tso);
			PackerTaskOP ptask = tf->create_task_and_apply_taskoperations(pose);

			TR<<"Effectuating mutation of residue " << pose.residue( target_aa_no ).name3() << " " << target_aa_no <<" to ";
			protocols::simple_moves::PackRotamersMoverOP pack;
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				pack = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::symmetry::SymPackRotamersMover( scorefxn(), ptask ) );
			} else {
				pack = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( scorefxn(), ptask ) );
			}
			pack->apply( pose );
			TR << pose.residue( target_aa_no ).name3() << " in pose. " << std::endl;
			(*scorefxn())(pose);
		} else {
			TR << "Mutation resulted in stop codon. Retrying. " << std::endl;
		}
	}

	TR << "nt seq after mut " << core::pose::get_all_comments(pose)[ "nt_seq" ] << std::endl;
	TR << "aa seq after mut " << core::pose::get_all_comments(pose)[ "aa_seq" ] << std::endl;

	if ( pose.sequence() != core::pose::get_all_comments(pose)[ "aa_seq" ] ) {
		TR << "The pose sequence is different from the target sequence. This shouldn't happen." << std::endl;
		TR << "Pose seq " << pose.sequence() << std::endl;
		TR << "target seq" << core::pose::get_all_comments(pose)[ "aa_seq" ] << std::endl;
		utility_exit_with_message("This should not be happening!");
	}
}

std::string
NucleotideMutation::get_name() const {
	return NucleotideMutationCreator::mover_name();
}

void
NucleotideMutation::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	init_sequence( tag->getOption< std::string >( "init_sequence", "" ) );
	continue_if_silent( tag->getOption< bool >( "continue_if_silent", true ) );

	// read the reference pose only if it doesn't exist already
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose( tag, data );
	} else throw utility::excn::EXCN_RosettaScriptsOption("Need to specify name under which to save pose.");

	if ( ( reference_pose_ != NULL /* if already read, don't reread */) && ( tag->hasOption( "reference_pdb_file" )) ) {
		std::string const template_pdb_fname( tag->getOption< std::string >( "reference_pdb_file" ));
		core::import_pose::pose_from_file( *reference_pose_, template_pdb_fname , core::import_pose::PDB_file);
		TR <<"reading in " << template_pdb_fname << " pdb with " << reference_pose_->total_residue() <<" residues"<<std::endl;
	}

	cache_task_ = tag->getOption< bool >( "cache_task", false );
}

protocols::moves::MoverOP
NucleotideMutation::clone() const {
	return( protocols::moves::MoverOP( new NucleotideMutation( *this ) ));
}

core::scoring::ScoreFunctionOP
NucleotideMutation::scorefxn() const{
	return scorefxn_;
}

void
NucleotideMutation::scorefxn( core::scoring::ScoreFunctionOP scorefxn )
{
	scorefxn_ = scorefxn;
}

core::pack::task::TaskFactoryOP
NucleotideMutation::task_factory() const{
	return( task_factory_ );
}

void
NucleotideMutation::task_factory( core::pack::task::TaskFactoryOP task_factory){
	task_factory_ = task_factory;
}

bool NucleotideMutation::cache_task() const {
	return( cache_task_ );
}

void NucleotideMutation::cache_task( bool cache ) {
	cache_task_ = cache;
}

} //movers
} //protein_interface_design
} //protocols
