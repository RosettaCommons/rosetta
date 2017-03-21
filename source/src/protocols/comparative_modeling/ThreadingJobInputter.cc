// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/comparative_modeling/ThreadingJobInputter.cc

///Unit headers
#include <protocols/comparative_modeling/ThreadingJobInputter.hh>
#include <protocols/comparative_modeling/ThreadingJobInputterCreator.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

///C++ headers
#include <string>

///Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

///Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>

#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/comparative_modeling/ExtraThreadingMover.hh>

namespace protocols {
namespace comparative_modeling {

static THREAD_LOCAL basic::Tracer tr( "protocols.comparative_modeling.ThreadingJobInputter" );

std::map< std::string, utility::vector1< core::Size > > read_extra_residues(
	utility::vector1< std::string > const & fns
) {
	std::map< std::string, utility::vector1< core::Size > > extra_residues;
	using core::Size;
	for ( Size ii = 1; ii <= fns.size(); ++ii ) {
		utility::io::izstream input( fns[ii] );
		if ( !input.good() ) {
			std::string const & msg( "Error reading file " + (std::string)fns[ii] );
			utility_exit_with_message(msg);
		}
		std::string line;
		while ( getline(input,line) ) {
			if ( line.substr(0,1) == "#" ) continue;
			std::istringstream ss(line);
			std::string aln_id;
			ss >> aln_id;
			aln_id = ObjexxFCL::uppercased(aln_id);
			extra_residues[aln_id] = utility::vector1< core::Size >();

			core::Size res(0);
			ss >> res;
			if ( res != 0 ) {
				extra_residues[aln_id].push_back(res);
			}
			while ( ss.good() ) {
				ss >> res;
				if ( res != 0 ) {
					extra_residues[aln_id].push_back(res);
				}
				res = 0;
			}
		}
		input.close();
	}

	return extra_residues;
}

ThreadingJobInputter::ThreadingJobInputter() :
	input_source_( protocols::jd2::JobInputterInputSource::NONE )
{
	using namespace core;
	using namespace core::pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	tr.Debug << "Instantiate ThreadingJobInputter" << std::endl;

	/// read alignments from command-line
	utility::vector1< std::string > const & aln_fns( option[ in::file::alignment ]() );
	for ( Size ii = 1; ii <= aln_fns.size(); ++ii ) {
		utility::vector1< core::sequence::SequenceAlignment > alns(
			sequence::read_aln( option[ cm::aln_format ](), aln_fns[ii] )
		);
		for ( Size jj = 1; jj <= alns.size(); ++jj ) {
			alignments_.push_back( alns[jj] );
		}
	}

	/// get template-pdbs from files
	if ( option[ in::file::template_pdb ].user() ) {
		FileList template_pdb_filenames = option[ in::file::template_pdb ]();
		typedef utility::vector1< pose::PoseOP > PoseOPvec;
		PoseOPvec poses = core::import_pose::poseOPs_from_files( template_pdb_filenames , core::import_pose::PDB_file);

		/// put template-pdbs into map --- use filename as key --- this is used to match pdb and alignment
		for ( PoseOPvec::const_iterator it = poses.begin(); it != poses.end(); ++it ) {
			utility::file::FileName fn( (*it)->pdb_info()->name() );
			std::string const base_fn( static_cast< std::string > (fn.base()) );
			std::string const match( ObjexxFCL::uppercased( base_fn.substr(0,5) ) );
			tr.Trace << "add template " << match << std::endl;
			template_poses_[ match ].push_back( *it );
		}
		input_source_ = protocols::jd2::JobInputterInputSource::PDB_FILE;
	} else if ( option[ in::file::template_silent ].user() ) { // get template-pdbs from silent-file
		core::io::silent::SilentFileOptions opts;
		io::silent::SilentFileData sfd( opts );
		sfd.read_file( option[ in::file::template_silent ]() );
		for ( io::silent::SilentFileData::iterator it = sfd.begin(); it != sfd.end(); ++it ) {
			PoseOP pose( new Pose );
			it->fill_pose( *pose );
			std::string const match( ObjexxFCL::uppercased( it->decoy_tag().substr(2,5) ) );
			tr.Trace << "add template " << match << std::endl;
			template_poses_[ match ].push_back( pose );
		}
		input_source_ = protocols::jd2::JobInputterInputSource::SILENT_FILE;
	} else {
		//no -in:file:template_xxx option
		utility_exit_with_message("ThreadingJobInputter needs parent-pdbs either as in:file:template_pdb or as in:file:template_silent");
	}

	// read in extra residues to steal
	if ( option[ cm::steal_extra_residues ].user() ) {
		utility::vector1< utility::file::FileName > const & fns( option[ cm::steal_extra_residues ]() );
		extra_residues_ = read_extra_residues(fns);
	} // steal_extra_residues
} // ThreadingJobInputter()

/// @details This function will first see if the pose already exists in the Job.
/// If not, it will read it into the pose reference, and hand a COP cloned from
/// that pose to the Job. If the pose pre-exists it just copies the COP's pose
/// into it.
void ThreadingJobInputter::pose_from_job(
	core::pose::Pose & pose,
	protocols::jd2::JobOP job
) {
	tr.Debug << "ThreadingJobInputter::pose_from_job" << std::endl;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using core::Real;
	using core::Size;

	///cast to ThreadingJob ... to access alignment and template pdb
	ThreadingJobCOP tjob = utility::pointer::dynamic_pointer_cast< protocols::comparative_modeling::ThreadingJob const > ( job->inner_job() );

	pose = core::pose::Pose();  //fpd  symmetry-safe
	std::string sequence;
	if ( option[ in::file::fasta ].user() ) {
		utility::vector1< core::sequence::SequenceOP > input_fasta
			= core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] );
		if ( input_fasta.size() == 0 ) {
			utility_exit_with_message(
				"ERROR: Fasta file specified by -in::file::fasta contains no valid sequence"
			);
		}

		if ( input_fasta.size() > 1 ) {
			utility_exit_with_message(
				"ERROR: Fasta file specified by -in::file::fasta should contain a single sequence"
			);
		}

		sequence = input_fasta[1]->sequence();
	} else {
		sequence = tjob->alignment().sequence(1)->ungapped_sequence();
	}

	core::pose::make_pose_from_sequence(
		pose,
		sequence,
		*( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ))
	);

	Real const alignment_coverage( tjob->alignment().length() - tjob->alignment().gapped_positions() );
	Real const alignment_identities( (core::Real)tjob->alignment().identities() / pose.size() );
	Real const alignment_perc ( alignment_coverage / pose.size() );

	// Add the alignment length, perc coverage and total length
	job->add_string_real_pair( "aln_len", alignment_coverage  );
	job->add_string_real_pair( "aln_perc", alignment_perc );
	job->add_string_real_pair( "aln_ident", alignment_identities );
	job->add_string_real_pair( "nres", pose.size() );

	// thread to get starting model
	comparative_modeling::ThreadingMover mover( tjob->alignment(), *(tjob->get_pose()) );
	mover.build_loops(false);

	// broken-chain folding from a sequence alignment and template pdb requires
	// that missing loop coordinates be randomized to avoid clashes. the option's
	// default value is false to match existing behavior.
	mover.randomize_loop_coords( option[OptionKeys::nonlocal::randomize_missing]() );
	mover.repack_query(false);
	mover.apply( pose );
	core::sequence::alignment_into_pose( tjob->alignment(), pose );

	// add extra residues from template here
	utility::vector1< core::Size > extra_res( tjob->extra_residues_to_steal() );
	if ( extra_res.size() > 0 ) {
		protocols::comparative_modeling::ExtraThreadingMover thief(
			tjob->alignment(), *tjob->get_pose(), extra_res);
		thief.apply(pose);
	}
} // pose_from_job

/// @details this function determines what jobs
void ThreadingJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs ) {
	tr.Debug << "ThreadingJobInputter::fill_jobs" << std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::Size;
	using core::Real;

	jobs.clear(); //should already be empty anyway

	//read command line
	Size const nstruct( get_nstruct() );
	Real filter_threshold = -1;

	if ( option[ cm::aln_length_filter_quantile ].user() ) {
		Real quantile = option[ cm::aln_length_filter_quantile ]();

		// make list of lengths
		//create jobs for each alignment
		std::vector < int > length_list;
		for ( Alignments::const_iterator align_it = alignments_.begin(); align_it != alignments_.end(); ++align_it ) {
			length_list.push_back(  align_it->length() - align_it->gapped_positions() );
			tr << "Len " <<  align_it->length() - align_it->gapped_positions() << std::endl;
		}

		auto  i = length_list.begin();
		std::vector< int >::size_type m = (size_t)( length_list.size() * quantile );

		std::nth_element(i, i + m, length_list.end());

		filter_threshold = length_list.at(m);

		tr << "Quantile filter threshold = " << filter_threshold << std::endl;
	}

	// create jobs for each alignment
	for ( Alignments::const_iterator align_it = alignments_.begin(); align_it != alignments_.end(); ++align_it ) {
		// alignment id
		std::string const alignment_id( align_it->alignment_id() );
		std::string const template_id( alignment_id.substr(0,5));
		tr.Debug << "creating job for alignment " << alignment_id << " on template " << template_id << std::endl;

		Real const alignment_coverage( align_it->length() - align_it->gapped_positions() );

		if ( option[ cm::aln_length_filter ].user() ) {
			filter_threshold = option[ cm::aln_length_filter ]();
		}

		if ( ( filter_threshold > 0 ) && ( alignment_coverage < filter_threshold ) ) {
			tr << "Skipping alignment " << alignment_id << ": length = "
				<< int( alignment_coverage )
				<< "  threshold = " << int( filter_threshold ) << std::endl;
			continue;
		}

		// find matching template pdb
		PoseMap::const_iterator iter = template_poses_.find( template_id );

		if ( iter != template_poses_.end() ) {
			// found template
			PoseOPs template_poses( iter->second );
			for ( PoseOPs::const_iterator it = template_poses.begin(); it != template_poses.end(); ++it ) {
				// create inner job
				ThreadingJobOP ijob( new ThreadingJob(
					*it, align_it->clone(), "S_" + alignment_id, nstruct
					) );
				// find extra residues
				ExtraResidues::const_iterator extra_res(extra_residues_.find(alignment_id));
				if ( extra_res != extra_residues_.end() ) {
					ijob->extra_residues_to_steal( extra_res->second );
				}

				// make nstruct outer jobs
				for ( Size index = 1; index <= nstruct; ++index ) {
					jobs.push_back( protocols::jd2::JobOP( new protocols::jd2::Job( ijob, index ) ) );
					jobs.back()->add_string_string_pair( "aln_id", alignment_id );
				} // loop over nstruct
			}
		} else { // report error if template pdb not found
			//utility_exit_with_message( "ERROR: no template_pdb provided for alignment " +  alignment_id );
			tr.Error << "no template pdb provided for alignment " << alignment_id << std::endl;
		}
	} // for alignments
} // fill_jobs

/// @brief Return the type of input source that the ThreadingJobInputter is
/// currently using for template structures.
protocols::jd2::JobInputterInputSource::Enum ThreadingJobInputter::input_source() const {
	return input_source_;
}

size_t ThreadingJobInputter::num_templates() const {
	size_t num_templates = 0;
	for ( auto i = template_poses_.begin();
			i != template_poses_.end(); ++i ) {
		++num_templates;
	}
	return num_templates;
}

//CREATOR SECTION
std::string
ThreadingJobInputterCreator::keyname() const
{
	return "ThreadingJobInputter";
}

protocols::jd2::JobInputterOP
ThreadingJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new ThreadingJobInputter );
}

} // comparative_modeling
} // protocols
