// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/evolution/NucleotideMutation.cc
/// @brief Substitutes amino acids in pose based on mutations in nucleotide space
/// @details Translate amino acid sequence to nucleotide sequence, introduce
/// random mutation, and affectuates the amino acid substitution in pose. Stop
/// codons are avoided and with default setting nucleotide mutations will continue
/// untill a non-silent mutation occurs.
/// @author Christoffer Norn (ch.norn@gmail.com)

// Unit headers
#include <protocols/evolution/NucleotideMutation.hh>
#include <protocols/evolution/NucleotideMutationCreator.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/kinematics/MoveMap.hh>


// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
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
#include <cmath>

#include <ctime>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace evolution {

using namespace std;
using namespace core::scoring;
using namespace core::kinematics;

using bools = utility::vector1<bool>;


static basic::Tracer TR( "protocols.evolution.NucleotideMutation" );

// XRW TEMP std::string
// XRW TEMP NucleotideMutationCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return NucleotideMutation::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP NucleotideMutationCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new NucleotideMutation );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP NucleotideMutation::mover_name()
// XRW TEMP {
// XRW TEMP  return "NucleotideMutation";
// XRW TEMP }

NucleotideMutation::NucleotideMutation() :
	Mover( NucleotideMutation::mover_name() ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	init_sequence_(""),
	continue_if_silent_( false ),
	flexbb_( true ),
	bbnbrs_( 0 ),
	reference_pose_( /* NULL */ )
{
}


NucleotideMutation::~NucleotideMutation() = default;

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
		for ( char i : aa_sequence ) {
			string nt = NucleotideTools::aa2randomCodon( i );
			nt_sequence.append( nt );
		}
		core::pose::add_comment(pose, "nt_seq", nt_sequence);
	}
}

void
NucleotideMutation::find_neighbors(
	bools const & is_mutated,
	core::pose::Pose const & pose,
	bools & is_flexible,
	core::Real const heavyatom_distance_threshold = 6.0 )
{
	Size const nres( pose.size() );
	is_flexible = is_mutated;

	for ( Size i=1; i<= nres; ++i ) {
		core::conformation::Residue const & rsd1( pose.residue(i) );
		if ( rsd1.is_virtual_residue() ) continue;
		for ( Size j=1; j<= nres; ++j ) {
			if ( !is_mutated[j] ) continue;
			if ( is_flexible[i] ) break;
			core::conformation::Residue const & rsd2( pose.residue(j) );
			if ( rsd2.is_virtual_residue() ) continue;
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) {
				is_flexible[i] = true;
			}
		}
	}
}

void
NucleotideMutation::compute_folding_energies(
	ScoreFunctionOP fa_scorefxn,
	Pose & pose,
	bools const & is_flexible,
	bools const & is_mutpos,
	Size bbnbrs=1,
	bool flexbb=true)
{
	protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
	fastrelax.cartesian( true );
	MoveMapOP movemap(new MoveMap);
	movemap->set_bb( false );
	movemap->set_chi( false );
	movemap->set_jump( false );

	for ( Size i=1; i<= pose.size(); ++i ) {
		if ( is_flexible[i] ) {
			movemap->set_chi( i, true );
			if ( flexbb > 0 ) {
				for ( Size j=0; j<=bbnbrs; j++ ) {
					if ( is_mutpos[i] || ( i+j <= pose.size() && is_mutpos[i+j] ) || ( j < i && is_mutpos[i-j] ) ) {
						movemap->set_bb ( i, true );
					}
				}
			}
		}
	}

	fastrelax.set_movemap( movemap );
	fastrelax.apply(pose);
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
		TR << "Nucleotide sequence will be initialized from comments " << std::endl;
	}

	/////////////////////////////////////////////////////////////////
	// find out which residues are designable and pick one randomly//
	/////////////////////////////////////////////////////////////////
	PackerTaskCOP task;
	if ( cache_task_ && task_ ) {
		if ( pose.size() == task_->total_residue() ) {
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
	for ( core::Size resi = 1; resi <= pose.size(); ++resi ) {
		if ( task->residue_task( resi ).being_designed() && pose.residue(resi).is_protein() ) {
			being_designed.push_back( resi );
		}
	}
	if ( being_designed.empty() ) {
		TR.Warning << "No residues are listed as designable." << std::endl;
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

		// Write the energy of the current pose
		std::ostringstream e_curr;
		e_curr.precision(16);
		core::Real energy_current = (*scorefxn())(pose);
		e_curr << energy_current;
		core::pose::add_comment(pose, "energy_current", e_curr.str()); //TR << "score current is " << score_current << std::endl;

		nt_sequence.replace( ( target_aa_no - 1 ) * 3, 3, new_codon);
		core::pose::add_comment(pose, "nt_seq", nt_sequence);
		core::pose::add_comment(pose, "mutation", mutation.str());


		if ( new_aa == "*" ) {
			// EvolutionaryDynamicsMover looks for stop_codon and sets fitness to 0 if found.
			core::pose::add_comment(pose, "stop_codon", "1");
			core::pose::add_comment(pose, "silent_mutation", "0");

			// Also, for analysis, it will probably be useful to set the energy of stops to something high.
			std::ostringstream e_proposed;
			e_proposed.precision(16);
			core::Real energy_proposed = 10000.0; // energy for chain break
			e_proposed << energy_proposed;
			core::pose::add_comment(pose, "energy_proposed", e_proposed.str());

			cont = false;
		} else if ( new_aa == old_aa ) { // silent mutation
			core::pose::add_comment(pose, "stop_codon", "0");
			core::pose::add_comment(pose, "silent_mutation", "1");

			core::pose::add_comment(pose, "energy_proposed", e_curr.str()); // if silent the energy is just the same as the current.

			if ( continue_if_silent() == false ) {
				TR << "Made a silent mutation and will attempt another mutation." << std::endl;
				cont = true;
			} else {
				TR << "Made a silent mutation" << std::endl;
				cont = false;
			}
		} else { // non-silent, non-stop codon, then make the mutation
			core::pose::add_comment(pose, "stop_codon", "0");
			core::pose::add_comment(pose, "silent_mutation", "0");
			TR << "Made a mutation resulting in an amino acid substitution." << std::endl;
			cont = false;

			// std::map<std::string, std::string> final_comments = core::pose::get_all_comments( pose );


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Let's mutate the residue in question
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			AA const target_aa( core::chemical::aa_from_oneletter_code( new_aa[0] ) );
			utility::vector1< bool > allowed_aas;
			allowed_aas.clear();
			allowed_aas.assign( num_canonical_aas, false );
			allowed_aas[ target_aa ] = true;

			PackerTaskOP mutate_residue( task->clone() );
			mutate_residue->initialize_from_command_line().or_include_current( true );
			for ( core::Size resi = 1; resi <= pose.total_residue(); ++resi ) {
				if ( resi != target_aa_no ) {
					mutate_residue->nonconst_residue_task( resi ).prevent_repacking(); // We could also do restrict_to_repacking instead of prevent_repacking
				} else {
					mutate_residue->nonconst_residue_task( resi ).restrict_absent_canonical_aas( allowed_aas );
				}
			}
			TR<<"Effectuating mutation of residue " << pose.residue( target_aa_no ).name3() << " " << target_aa_no <<" to ";
			protocols::minimization_packing::PackRotamersMoverOP pack;
			if ( core::pose::symmetry::is_symmetric( pose ) ) {
				pack = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::symmetry::SymPackRotamersMover( scorefxn(), mutate_residue ) );
			} else {
				pack = protocols::minimization_packing::PackRotamersMoverOP( new protocols::minimization_packing::PackRotamersMover( scorefxn(), mutate_residue ) );
			}
			pack->apply( pose );
			TR << pose.residue( target_aa_no ).name3() << " in pose. " << std::endl;



			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Now lets relax the pose
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// suboptions:
			// It can be full bb flexibility then bbnbrs just needs to be > pose length. Takes 80 seconds per trial (87 residues, cutoff=8).
			// It can be with localized bb flex (bbnbrs = 1). Takes 20 seconds per trial (87 residues, cutoff=8).
			// It can be without bb flexi then flexbb=false, and no definition for bbnbrs. Takes 1.6 second per trial (87 residues, cutoff=8).


			bools is_mutated( pose.size(), false );   //mutated position
			bools is_flexible( pose.size(), false );  //repackable
			is_mutated[ target_aa_no ] = true;
			core::Real cutoff = 8.0;
			find_neighbors( is_mutated, pose, is_flexible, cutoff );

			NucleotideMutation::compute_folding_energies( scorefxn(), pose, is_flexible, is_mutated, bbnbrs(), flexbb() ); // this takes 33 seconds for 1PTF (87 residues), with bbnbr=1

			std::ostringstream e_prop;
			e_prop.precision(16);
			core::Real score_proposed = (*scorefxn())(pose);
			e_prop << score_proposed;
			core::pose::add_comment(pose, "energy_proposed", e_prop.str()); //TR << "score_proposed is " << score_proposed << std::endl;

		}
	}

	TR << "nt seq after mut " << core::pose::get_all_comments(pose)[ "nt_seq" ] << std::endl;
	TR << "aa seq after mut " << core::pose::get_all_comments(pose)[ "aa_seq" ] << std::endl;
}

// XRW TEMP std::string
// XRW TEMP NucleotideMutation::get_name() const {
// XRW TEMP  return NucleotideMutation::mover_name();
// XRW TEMP }

void
NucleotideMutation::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	init_sequence( tag->getOption< std::string >( "init_sequence", "" ) );
	continue_if_silent( tag->getOption< bool >( "continue_if_silent", false ) );
	flexbb( tag->getOption< bool >("flexbb", true ) );
	bbnbrs( tag->getOption< Size >("bbnbrs", 0 ) );
	runtime_assert_string_msg( !(bbnbrs()>0 && flexbb()==false), "You have bbnbrs>0, but flexbb=false?! The bbnbrs sets which backbone neighbor residues that should be flexible" );

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

std::string NucleotideMutation::get_name() const {
	return mover_name();
}

std::string NucleotideMutation::mover_name() {
	return "NucleotideMutation";
}

void NucleotideMutation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute( "init_sequence", xs_string, "Initial sequence" )
		+ XMLSchemaAttribute::attribute_w_default( "continue_if_silent", xsct_rosetta_bool, "Make another mutation if the first mutation is silent", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "flexbb", xsct_rosetta_bool, "Let backbone atoms be flexible during relax", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "bbnbrs", xsct_non_negative_integer, "Let the +/- n (dflt=1) residues from the mutated residue be bb flexible", "1" )
		// AMW XRW TODO: do we have a helper for reference pose parsing?
		//+ XMLSchemaAttribute( "reference_name", xs_string, "Saved reference pose" )
		//+ XMLSchemaAttribute( "reference_pdb_file", xs_string, "Saved reference pdb" )
		+ XMLSchemaAttribute::attribute_w_default( "cache_task", xsct_rosetta_bool, "Cache the initially calculated packer task", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string NucleotideMutationCreator::keyname() const {
	return NucleotideMutation::mover_name();
}

protocols::moves::MoverOP
NucleotideMutationCreator::create_mover() const {
	return protocols::moves::MoverOP( new NucleotideMutation );
}

void NucleotideMutationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NucleotideMutation::provide_xml_schema( xsd );
}


} //evolution
} //protocols
