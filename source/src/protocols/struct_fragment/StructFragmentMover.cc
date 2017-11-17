// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   StructFragmentMover.cc
/// @brief  Creating fragments from supplied structures.
/// @author andreas scheck (graziano.bud@gmail.com), Correia's LPDI/EPFL


// Unit Headers
#include <protocols/struct_fragment/StructFragmentMover.hh>
#include <protocols/struct_fragment/StructFragmentMoverCreator.hh>

// Package Headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/BoundedCollector.hh>

// Project Headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/moves/DsspMover.hh>

// Utility Headers
#include <basic/database/open.hh>
#include <basic/datacache/ConstDataMap.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/Tag.hh>

// XSD XRW Includes
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace struct_fragment {

namespace cf = core::fragment;
namespace cp = core::pose;
namespace pfp = protocols::frag_picker;

static basic::Tracer TR("protocols.struct_fragment.StructFragmentMover");

StructFragmentMover::StructFragmentMover():
	structPicker_( new pfp::FragmentPicker() ),
	loop_angle_conf_( default_value_for_loop_angle_conf() ),  // default confidence for phi/psi angles instead of 1 (which would be correct because we know them from the acutal structure) to allow for more varied search results
	small_frag_size_( default_value_for_small_frag_size() ),
	large_frag_size_( default_value_for_large_frag_size() ),
	small_frag_file_( default_value_for_empty_file() ),
	large_frag_file_( default_value_for_empty_file() ),
	output_frag_files_( default_value_for_output_frag_files() ),
	steal_small_frags_( default_value_for_steal_small_frags() ),
	steal_large_frags_( default_value_for_steal_large_frags() ),
	changed_frags_( false ),
	frag_weight_file_( default_value_for_empty_file() ),
	sequence_profile_( default_value_for_empty_file() ),
	vall_file_( default_value_for_empty_file() )
{
	structPicker_->n_candidates_ = default_value_for_n_candidates();
	structPicker_->n_frags_      = default_value_for_n_frags();
	structPicker_->prefix_       = default_value_for_prefix();
}

StructFragmentMover::~StructFragmentMover()
{}

core::Size
StructFragmentMover::evaluate_job() {
	if (  not ( small_frag_file_.empty() || large_frag_file_.empty() ) ) { // We have defined fragment files (both!)
		if ( utility::file::file_exists( small_frag_file_ ) && utility::file::file_exists( large_frag_file_ ) ) { // And they exist!
			return 0; // STATUS 0: Load Fragments From File.
		} else if ( !vall_file_.empty() && utility::file::file_exists( vall_file_ ) ) { // If they don't exist we are making for them!
			return 1; // STATUS 1: Create Fragments.
		} else { // We don't have it and cannot create it
			if ( steal_small_frags_ && steal_large_frags_ ) { // We can steal if we want to steal all.
				return 2; // STATUS 2: Only Steal Fragments.
			} else {
				throw( CREATE_EXCEPTION(utility::excn::Exception, "With the provided data Fragments cannot be created, not readed." ) );
			}
		}
	}
	return 2;
}

void StructFragmentMover::apply( cp::Pose & pose ) {
	if ( not pose.empty() ) {

		TR.Trace << TR.Green << "small fragments pointer address " << smallF_ << " size: " << smallF_->size() << TR.Reset << std::endl;
		TR.Trace << TR.Green << "large fragments pointer address " << largeF_ << " size: " << largeF_->size() << TR.Reset << std::endl;

		if ( smallF_->size() == 0 || largeF_->size() == 0 ) { //  With this we make sure that this is only runned once in the RosettaScripts
			utility::vector1< cf::FragSetOP > frag_result;
			core::Size status = evaluate_job();

			if (  status == 0 ) {
				TR << "Small fragment file " << small_frag_file_ << " and large fragment file " << large_frag_file_ << " are provided." << std::endl;
				TR << "Reading fragments from provided fragment files." << std::endl;
				frag_result = read_frag_files( small_frag_file_, large_frag_file_ );
			} else if ( status == 1 ) {
				TR << "[Calculating Fragments]" << std::endl;
				frag_result = get_fragments( pose );
			}

			frag_result = steal_fragments( pose, frag_result);

			for ( core::Size i = 1; i <= frag_result.size(); ++i ) {
				if ( frag_result[i]->min_pos() != 1 ) {
					frag_result[i]->shift_by( 1 - frag_result[i]->min_pos() );
				}
			}

			if ( output_frag_files_ && changed_frags_ ) {
				TR.Debug << "Printing fragment files." << std::endl;
				cf::FragmentIO().write_data( structPicker_->prefix_ + "." + utility::to_string(structPicker_->n_frags_) + "." + utility::to_string(small_frag_size_) + "mers", *frag_result[1] );
				cf::FragmentIO().write_data( structPicker_->prefix_ + "." + utility::to_string(structPicker_->n_frags_) + "." + utility::to_string(large_frag_size_) + "mers", *frag_result[2] );
			}

			smallF_->add( *frag_result[1] );
			largeF_->add( *frag_result[2] );
		}

		TR.Trace << TR.Green << "small fragments pointer address " << smallF_ << " size: " << smallF_->size() << TR.Reset << std::endl;
		TR.Trace << TR.Green << "large fragments pointer address " << largeF_ << " size: " << largeF_->size() << TR.Reset << std::endl;
	} else {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Pose was empty.");
	}
} // end apply method

utility::vector1< cf::FragSetOP > StructFragmentMover::get_fragments( cp::Pose const & pose ) {
	TR.Debug << TR.Bold << "[MANAGING TEMPLATE FRAGMENTS]" << TR.Reset << std::endl;

	// Pose Length
	core::Size pose_length = pose.total_residue();
	TR.Trace << "Pose length is " << pose_length << std::endl;

	// First residue
	core::Size first_res = pose.pdb_info()->number(1);


	// Compute new frag files
	// Declare empty fragment sets for small and large fragments
	cf::FragSetOP mfrag = cf::FragSetOP( new cf::ConstantLengthFragSet( small_frag_size_ ) );
	cf::FragSetOP Mfrag = cf::FragSetOP( new cf::ConstantLengthFragSet( large_frag_size_ ) );

	TR.Debug << "Calculating new fragments." << std::endl;
	if ( not ( steal_small_frags_ && steal_large_frags_ ) && vall_file_.empty() ) {
		throw ( CREATE_EXCEPTION(utility::excn::BadInput,  "with no frag_file and no steal, a vall file must be provided." ) );
	}
	if ( not vall_file_.empty() && utility::file::file_exists( vall_file_ ) ) {
		pfp::FragmentPickerOP pickIt = make_fragment_picker( pose, vall_file_ );
		TR.Trace << "FragmentPicker ready!" << std::endl;
		collector_to_picker( pickIt, small_frag_size_, pose_length );
		collector_to_picker( pickIt, large_frag_size_, pose_length );

		pickIt->pick_candidates();

		mfrag = get_fragset( pickIt, first_res, small_frag_size_ );
		Mfrag = get_fragset( pickIt, first_res, large_frag_size_ );
	} else {
		throw ( CREATE_EXCEPTION(utility::excn::BadInput,  "The provided vall file could not be found!" ) );
	}

	// if ( steal_small_frags_ || steal_large_frags_ ) {
	//     cp::Pose ideal_pose;
	//     ideal_pose = cp::Pose( pose );
	//     protocols::idealize::IdealizeMover idealizer;
	//     idealizer.fast( false );
	//     idealizer.apply( ideal_pose );
	//     if ( steal_small_frags_ ) mfrag = steal_fragments_by_size( ideal_pose, mfrag, small_frag_size_);
	//     if ( steal_large_frags_ ) Mfrag = steal_fragments_by_size( ideal_pose, Mfrag, large_frag_size_);
	// }

	TR.Debug << TR.Bold << "[            END            ]" << TR.Reset << std::endl;
	utility::vector1< cf::FragSetOP > result;
	result.push_back( mfrag );
	result.push_back( Mfrag );

	changed_frags_ = true;
	return result;
} // end get_fragments

utility::vector1< cf::FragSetOP > StructFragmentMover::steal_fragments( cp::Pose const & pose, utility::vector1< cf::FragSetOP > fsets ) {
	if ( steal_small_frags_ || steal_large_frags_ ) {
		cp::Pose ideal_pose;
		ideal_pose = cp::Pose( pose );
		protocols::idealize::IdealizeMover idealizer;
		idealizer.fast( false );
		idealizer.apply( ideal_pose );
		// cp:PDBInfo info = cp::PDBInfo( pose.pdb_info() );
		ideal_pose.pdb_info()->name( pose.pdb_info()->name() );
		if ( steal_small_frags_ ) fsets[1] = steal_fragments_by_size( ideal_pose, fsets[1], small_frag_size_);
		if ( steal_large_frags_ ) fsets[2] = steal_fragments_by_size( ideal_pose, fsets[2], large_frag_size_);
	}
	changed_frags_ = true;
	return fsets;
}

cf::FragSetOP StructFragmentMover::steal_fragments_by_size( cp::Pose const & pose, cf::FragSetOP fset, core::Size size ) {
	TR.Debug << TR.Underline << "Stealing fragments of size " << size << TR.Reset << std::endl;
	cf::SingleResidueFragDataOP rfdata = cf::SingleResidueFragDataOP( new cf::BBTorsionSRFD );
	cf::FragDataCOP fdata = cf::FragDataCOP( cf::FragDataOP( new cf::FragData( rfdata, size ) ) );
	cf::FrameOP frame;
	for ( core::Size pos = 1; pos <= pose.total_residue() - size + 1; ++pos ) {
		frame = cf::FrameOP( new cf::Frame( pos, fdata ) );
		frame->steal( pose );
		fset->add( frame );
	}
	return fset;
} // end steal_fragments

pfp::FragmentPickerOP StructFragmentMover::make_fragment_picker( cp::Pose pose, std::string vall_file_name ) {
	core::Size num_res = pose.total_residue();

	// get secondary structure
	protocols::moves::DsspMover dssp;
	dssp.apply( pose );

	// get sequence and structure; DSSP has to be executed before the SS can be set
	std::string sequence = pose.sequence();
	std::string structure = pose.secstruct();

	// get accessible surface area
	core::scoring::sasa::SasaCalcOP sasa_calc( new core::scoring::sasa::SasaCalc() );
	sasa_calc->calculate( pose );   // calculate sasa for whole pose first before extracting the per residue sasa
	utility::vector1< core::Real > res_sasa = sasa_calc->get_residue_sasa();

	// angles
	utility::vector1< core::Real > psi_angles;
	utility::vector1< core::Real > phi_angles;
	// get angles for all residues
	for ( core::Size pos = 1; pos <= num_res; ++pos ) {
		psi_angles.push_back( pose.psi( pos ) );
		phi_angles.push_back( pose.phi( pos ) );
	}

	// set default confidence for phi and psi angles
	utility::vector1<core::Real> angle_conf( num_res, loop_angle_conf_ );

	// set all values needed for fragment picking
	core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );
	if ( !sequence_profile_.empty() ) {
		if ( utility::file::file_exists( sequence_profile_ ) ) {
			q_prof->read_from_file( sequence_profile_ );
		} else {
			throw( CREATE_EXCEPTION(utility::excn::Exception, "The provided sequence profile file could not be found!" ) );
		}
	} else {
		q_prof->generate_from_sequence( core::sequence::Sequence( sequence, "sequence" ) );
	}
	//  q_prof->convert_profile_to_probs(1.0);
	structPicker_->set_query_profile(q_prof);

	structPicker_->add_query_ss( structure, "native" );
	structPicker_->set_query_sa( res_sasa );
	structPicker_->set_query_psi( psi_angles );
	structPicker_->set_query_phi( phi_angles );
	structPicker_->set_query_psi_conf( angle_conf );
	structPicker_->set_query_phi_conf( angle_conf );

	structPicker_->read_vall( vall_file_name );
	TR << "Loading vall file " << vall_file_name << std::endl;

	// Setup Scores
	pfp::scores::FragmentScoreManagerOP fscore = structPicker_->get_score_manager();
	if ( frag_weight_file_.empty() ) {
		fscore->create_scoring_method( "SecondaryIdentity", 350, 1.0, 0.0, false, structPicker_, "native" );
		TR << "SecondaryIdentity score added" << std::endl;
		std::string blosum_file = basic::database::full_name( "sequence/substitution_matrix/BLOSUM62" );
		fscore->create_scoring_method( "ProfileScoreSubMatrix", 200, 1.0, 0.0, false, structPicker_, blosum_file );
		// fscore->create_scoring_method( "ProfileScoreL1", 200, 1.0, 0.0, false, structPicker_, "sequence" );
		// TR << "ProfileScoreSubMatrix score added" << std::endl;
	} else {
		if ( utility::file::file_exists( frag_weight_file_ ) ) {
			fscore->create_scores( frag_weight_file_, structPicker_ );
		} else {
			throw ( CREATE_EXCEPTION(utility::excn::BadInput,  "The provided weight file could not be found!" ) );
		}
	}

	pfp::FragmentSelectingRuleOP selector( new pfp::BestTotalScoreSelector( structPicker_->n_frags_, fscore ) );
	structPicker_->selector_ = selector;
	return structPicker_;
} // end make_fragment_picker

void StructFragmentMover::collector_to_picker( pfp::FragmentPickerOP pickIt, core::Size size, core::Size seqlength ) {
	pfp::CompareTotalScore comparator( pickIt->get_score_manager() );
	pickIt->frag_sizes_.push_back( size );
	pfp::CandidatesCollectorOP collector( new pfp::BoundedCollector<pfp::CompareTotalScore> (
		seqlength,             // collector must know the size of query
		pickIt->n_candidates_, // how many candidates to collect
		comparator,            // yes, here the comparator comes to sort fragments within the collector
		pickIt->get_score_manager()->count_components()
		) );
	pickIt->set_candidates_collector( size, collector );
} //end collector_to_picker

cf::FragSetOP StructFragmentMover::get_fragset ( pfp::FragmentPickerOP pickIt, core::Size position, core::Size size ) {
	core::Size residueInPose = position;
	cf::FragSetOP myFragSet  = cf::FragSetOP( new cf::ConstantLengthFragSet( size ) );
	pfp::CandidatesCollectorOP storage = pickIt->get_candidates_collector( size );
	for ( core::Size qPos = 1; qPos <= pickIt->size_of_query(); ++qPos ) {
		if ( storage->get_candidates( qPos ).size() == 0 ) continue;
		pfp::Candidates out;
		pickIt->selector_->select_fragments( storage->get_candidates( qPos ), out );
		cf::FrameOP frame( new cf::Frame(residueInPose++) );
		for ( core::Size fi = 1; fi <= out.size(); ++fi ) {
			cf::FragDataOP current_fragment( NULL );

			for ( core::Size i = 1; i <= out[1].first->get_length(); ++i ) {
				pfp::VallResidueOP r = out[fi].first->get_residue(i);
				std::string pdbid    = out[fi].first->get_pdb_id();
				core::Size index     = r->resi();
				char aa              = toupper(r->aa());
				char ss              = r->ss();
				core::Real phi       = r->phi();
				core::Real psi       = r->psi();
				core::Real omega     = r->omega();

				if ( i == 1 ) {
					current_fragment = cf::FragDataOP( new cf::AnnotatedFragData( pdbid, index ) );
				}

				cf::BBTorsionSRFDOP res_torsions( new cf::BBTorsionSRFD(3,ss,aa) );
				res_torsions->set_torsion   ( 1, phi   );
				res_torsions->set_torsion   ( 2, psi   );
				res_torsions->set_torsion   ( 3, omega );
				res_torsions->set_secstruct ( ss );

				// Add residue to fragment
				current_fragment->add_residue( res_torsions );

			} // End VallResidue loop

			if ( current_fragment ) { // != NULL) {
				current_fragment->set_valid(); //it actually containts data

				// Add fragment to frame
				if ( !frame->add_fragment(current_fragment) ) {
					current_fragment->show( TR.Fatal );
					exit( 1111 );
				}
			}
		} // End FragmentCandidate loop
		// Add frame to myFragSet.
		myFragSet->add(frame);
	}
	return myFragSet;
} // end get_fragset



cf::FragSetOP StructFragmentMover::read_frag_file( std::string frag_file) {
	if ( utility::file::file_exists( frag_file ) ) {
		cf::FragSetOP frags = cf::FragmentIO().read_data( frag_file );
		core::Size frag_length = frags->max_frag_length();
		TR.Debug << "Fragments (" << frag_length << ") read from " << frag_file <<std::endl;
		return frags;
	} else {
		throw ( CREATE_EXCEPTION(utility::excn::BadInput,  "The provided fragment file could not be found!" ) );
	}
}


utility::vector1< cf::FragSetOP > StructFragmentMover::read_frag_files( std::string small_frag_file, std::string large_frag_file ) {
	cf::FragSetOP mfrags = read_frag_file( small_frag_file );
	cf::FragSetOP Mfrags = read_frag_file( large_frag_file );

	utility::vector1< cf::FragSetOP > result;
	result.push_back( mfrags) ;
	result.push_back( Mfrags );
	return result;
}

void StructFragmentMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	StructFragmentMover::provide_xml_schema( xsd );
}

std::string StructFragmentMoverCreator::keyname() const {
	return StructFragmentMover::mover_name();
}

protocols::moves::MoverOP
StructFragmentMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StructFragmentMover );
}

void
StructFragmentMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// -- Attributes with default values -- //
	small_frag_size_             = tag->getOption< core::Size >( "small_frag_size", default_value_for_small_frag_size() );
	large_frag_size_             = tag->getOption< core::Size >( "large_frag_size", default_value_for_large_frag_size() );
	output_frag_files_           = tag->getOption< bool >( "output_frag_files", default_value_for_output_frag_files() );
	steal_small_frags_           = tag->getOption< bool >( "steal_small_frags", default_value_for_steal_small_frags() );
	steal_large_frags_           = tag->getOption< bool >( "steal_large_frags", default_value_for_steal_large_frags() );
	structPicker_->n_candidates_ = tag->getOption< core::Size >( "n_candidates", default_value_for_n_candidates() );
	structPicker_->n_frags_      = tag->getOption< core::Size >( "n_frags", default_value_for_n_frags() );
	structPicker_->prefix_       = tag->getOption< std::string >( "prefix", default_value_for_prefix() );

	// -- Linked Attributes (depend of one another) -- //
	small_frag_file_ = tag->getOption< std::string >( "small_frag_file", default_value_for_empty_file() );
	large_frag_file_ = tag->getOption< std::string >( "large_frag_file", default_value_for_empty_file() );
	vall_file_       = tag->getOption< std::string >( "vall_file", default_value_for_empty_file() );

	// Link1: small_frag_file & large_frag_file need to be either both specified or none.
	if ( ( small_frag_file_.empty() && not large_frag_file_.empty() ) || ( not small_frag_file_.empty() && large_frag_file_.empty() ) ) {
		throw ( CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Please make sure to either provide the small and large fragment files or none of them!" ) );
	}
	// Link2: vall_file is only optional if small_frag_file & large_frag_file are specified, is mandatory otherwise
	if ( small_frag_file_.empty() && vall_file_.empty() ) {
		throw ( CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "To create fragments a vall file needs to be provided: vall_file option" ) );
	}

	// -- Optional Attributes -- //
	frag_weight_file_ = tag->getOption< std::string >( "frag_weight_file", default_value_for_empty_file() );
	sequence_profile_ = tag->getOption< std::string >( "sequence_profile", default_value_for_empty_file() );

	smallF_ = core::fragment::FragSetOP(new core::fragment::ConstantLengthFragSet );
	largeF_ = core::fragment::FragSetOP(new core::fragment::ConstantLengthFragSet );
	if ( !data.has( "FragSet", "small" ) && !data.has( "FragSet", "large" ) ) {
		data.add( "FragSet", "small", smallF_ );
		data.add( "FragSet", "large", largeF_ );
		TR.Trace << TR.Green << "small fragments pointer address " << smallF_ << TR.Reset << std::endl;
		TR.Trace << TR.Green << "large fragments pointer address " << largeF_ << TR.Reset << std::endl;
	}

}

void StructFragmentMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "small_frag_size", xsct_non_negative_integer, "Set the size of the small fragments", utility::to_string( default_value_for_small_frag_size() ) )
		+ XMLSchemaAttribute::attribute_w_default( "large_frag_size", xsct_non_negative_integer, "Set the size of the large fragments", utility::to_string( default_value_for_large_frag_size() ) )
		+ XMLSchemaAttribute( "small_frag_file", xs_string, "Path to input small fragment file (requires large_frag_file too)" )
		+ XMLSchemaAttribute( "large_frag_file", xs_string, "Path to input large fragment file (requires small_frag_file too)" )
		+ XMLSchemaAttribute::attribute_w_default( "output_frag_files", xsct_rosetta_bool, "Write fragment files", utility::to_string( default_value_for_output_frag_files() ) )
		+ XMLSchemaAttribute::attribute_w_default( "steal_small_frags", xsct_rosetta_bool, "Steal small fragments", utility::to_string( default_value_for_steal_small_frags() ) )
		+ XMLSchemaAttribute::attribute_w_default( "steal_large_frags", xsct_rosetta_bool, "Steal large fragments", utility::to_string( default_value_for_steal_large_frags() ) )
		+ XMLSchemaAttribute( "frag_weight_file", xs_string, "Path to input weight file" )
		+ XMLSchemaAttribute( "sequence_profile", xs_string, "Path to input sequence profile" )
		+ XMLSchemaAttribute( "vall_file", xs_string, "Path to vall file (required when no fragment files are provided)" )
		+ XMLSchemaAttribute::attribute_w_default( "n_candidates", xsct_non_negative_integer, "Set the number of candidates", utility::to_string( default_value_for_n_candidates() ) )
		+ XMLSchemaAttribute::attribute_w_default( "n_frags", xsct_non_negative_integer, "Set the number of fragments", utility::to_string( default_value_for_n_frags() ) )
		+ XMLSchemaAttribute::attribute_w_default( "prefix", xs_string, "Prefix of output files", default_value_for_prefix() );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Create fragments from a supplied pose.", attlist );
}

} //end namespace struct_fragment
} // end namespace protocols
