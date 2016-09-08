// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Interface ddG bind of mutation prediction protocol
/// @file   src/apps/pilot/andrew/interface_ddg_bind.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// core headers
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/ResidueIndexDescription.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/import_pose/import_pose.hh>

// Protocols headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/rigid/RigidBodyMover.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/Tracer.hh>

// Devel headers
#include <devel/init.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

// C++ headers
#include <map>
#include <string>


///// forward declarations //////
class CloseContactWithResidue;
class InterfaceDDGMutationTask;
class InterfaceDDGBindInnerJob;
class InterfaceDDGBindJobInputter;
class InterfaceDDGMover;

typedef utility::pointer::shared_ptr< CloseContactWithResidue > CloseContactWithResidueOP;
typedef utility::pointer::shared_ptr< CloseContactWithResidue const > CloseContactWithResidueCOP;
typedef utility::pointer::shared_ptr< InterfaceDDGMutationTask > InterfaceDDGMutationTaskOP;
typedef utility::pointer::shared_ptr< InterfaceDDGMutationTask const > InterfaceDDGMutationTaskCOP;
typedef utility::pointer::shared_ptr< InterfaceDDGBindInnerJob > InterfaceDDGBindInnerJobOP;
typedef utility::pointer::shared_ptr< InterfaceDDGBindInnerJob const > InterfaceDDGBindInnerJobCOP;
typedef utility::pointer::shared_ptr< InterfaceDDGBindJobInputter > InterfaceDDGBindJobInputterOP;
typedef utility::pointer::shared_ptr< InterfaceDDGBindJobInputter const > InterfaceDDGBindJobInputterCOP;
typedef utility::pointer::shared_ptr< InterfaceDDGMover > InterfaceDDGMoverOP;
typedef utility::pointer::shared_ptr< InterfaceDDGMover const > InterfaceDDGMoverCOP;

// Macro headers
OPT_1GRP_KEY( File, interface_ddg, jobs )
OPT_1GRP_KEY( Boolean, interface_ddg, ddg_sep_from_complex )

using namespace core;

class CloseContactWithResidue : public core::select::residue_selector::ResidueSelector
{
public:
	CloseContactWithResidue( core::Size seqpos ) : seqpos_( seqpos ) {}

	core::select::residue_selector::ResidueSelectorOP
	clone() const {
		return core::select::residue_selector::ResidueSelectorOP( new CloseContactWithResidue( seqpos_ ));
	}

	std::string get_name() const { return "CloseContactWithResidue"; }

	virtual
	core::select::residue_selector::ResidueSubset
	apply(
		core::pose::Pose const & pose
	) const
	{
		core::select::residue_selector::ResidueSubset subset( pose.size(), false );
		subset[ seqpos_ ] = true;

		core::conformation::Residue const & focused_rsd( pose.residue( seqpos_ ) );
		core::Vector focused_nbr_atom = focused_rsd.xyz( focused_rsd.nbr_atom() );

		core::Real proximity_cutoff = 4.5;
		core::Real cut2 =  proximity_cutoff * proximity_cutoff;
		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( ii == seqpos_ ) continue;

			core::conformation::Residue const & iirsd( pose.residue(ii) );
			Real cutoff = focused_rsd.nbr_radius() + proximity_cutoff + iirsd.nbr_radius();
			if ( iirsd.xyz( iirsd.nbr_atom() ).distance_squared( focused_nbr_atom ) > cutoff * cutoff ) continue; // too far

			for ( core::Size jj = 1; jj <= iirsd.natoms(); ++jj ) {
				for ( core::Size kk = 1; kk <= focused_rsd.natoms(); ++kk ) {
					if ( iirsd.xyz( jj ).distance_squared( focused_rsd.xyz( kk ) ) < cut2 ) {
						subset[ ii ] = true;
						break;
					}
				}
				if ( subset[ ii ] ) break;
			}
		}
		return subset;
	}



private:
	core::Size seqpos_;
};

struct input_mutation {
	core::pose::ResidueIndexDescriptionFromFileOP resind;
	core::chemical::AA mutaa;
};

class InterfaceDDGMutationTask : public basic::datacache::CacheableData {
public:
	InterfaceDDGMutationTask();
	InterfaceDDGMutationTask( core::Size seqpos, bool complex );
	InterfaceDDGMutationTask( core::Size seqpos, bool complex, core::chemical::AA new_aa );
	virtual ~InterfaceDDGMutationTask();

	virtual basic::datacache::CacheableDataOP clone() const;

	/// @brief should the complex be modeled or should it be separated?
	bool complex() const;

	bool unininitialized() const;

	/// @brief Does this mutation perscribe to maintain the wild type amino acid?
	bool wt() const;

	/// @brief Return the pose-indexed sequence position that should be mutated
	core::Size seqpos() const;

	core::chemical::AA new_aa() const;

private:
	core::Size seqpos_;
	bool wt_;
	bool complex_;
	core::chemical::AA new_aa_;
};


class InterfaceDDGBindInnerJob : public protocols::jd2::InnerJob
{
public:
	InterfaceDDGBindInnerJob(
		std::string const & input_tag,
		core::Size nstruct
	);

	virtual ~InterfaceDDGBindInnerJob();

	//InterfaceDDGMutationTask const & mutation() const;
	//void mutation( InterfaceDDGMutationTask const & mut );

	bool complex() const;

	/// @brief Does this mutation perscribe to maintain the wild type amino acid?
	bool wt() const;

	utility::file::FileName const & input_pdb() const;
	void input_pdb( utility::file::FileName const & pdb_name );

	void complex( bool setting );
	void residue_index_description( core::pose::ResidueIndexDescriptionOP setting );
	void new_aa( core::chemical::AA setting );

	//core::pose::ResidueIndexDescription & residue_index_description();
	core::pose::ResidueIndexDescription const & residue_index_description() const;
	core::chemical::AA new_aa() const;

	virtual bool same( protocols::jd2::InnerJob const & other ) const;
	InterfaceDDGBindInnerJobOP reference_job() const;
	void reference_job( InterfaceDDGBindInnerJobOP setting );

private:
	utility::file::FileName input_pdb_;
	InterfaceDDGBindInnerJobOP reference_job_;
	// InterfaceDDGMutationTask mutation_;

	core::pose::ResidueIndexDescriptionOP res_ind_desc_;
	bool wt_;
	bool complex_;
	core::chemical::AA new_aa_;

};

class InterfaceDDGBindJobInputter : public protocols::jd2::JobInputter
{
public:
	InterfaceDDGBindJobInputter();
	virtual ~InterfaceDDGBindJobInputter();

	virtual void pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job );
	virtual void fill_jobs( protocols::jd2::JobsContainer & jobs );
	virtual protocols::jd2::JobInputterInputSource::Enum input_source() const;
};

// class InterfaceDDGBindJobInputterCreator : public protocols::jd2::JobInputterCreator {
// public:
//  virtual JobInputterOP create_JobInputter() const;
//  virtual std::string keyname() const;
// };


class InterfaceDDGMover : public protocols::moves::Mover
{
public:

	virtual void apply( core::pose::Pose & pose );
	virtual protocols::moves::MoverOP clone() const;
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	void sfxn( core::scoring::ScoreFunctionOP sfxn );

private:

	core::Size determine_central_seqpos( core::pose::Pose const & pose );
	bool determine_if_wt( core::pose::Pose const & pose );
	core::chemical::AA determine_mutation( core::pose::Pose const & pose );
	bool determine_if_complex( core::pose::Pose const & pose );

	void separate_complex( core::pose::Pose & pose );
	void mutate_and_relax(
		core::pose::Pose & pose,
		core::Size seqpos,
		bool wt,
		core::chemical::AA mutaa,
		core::select::residue_selector::ResidueSubset const & subset
	);

private:
	core::scoring::ScoreFunctionCOP sfxn_;

};

////////////////// InterfaceDDGMutationTask /////////////////////////


InterfaceDDGMutationTask::InterfaceDDGMutationTask() :
	seqpos_( 0 ),
	wt_( false ),
	complex_( false ),
	new_aa_( core::chemical::aa_unk )
{}

InterfaceDDGMutationTask::InterfaceDDGMutationTask( core::Size seqpos, bool complex ) :
	seqpos_( seqpos ),
	wt_( true ),
	complex_( complex ),
	new_aa_( core::chemical::aa_unk )
{}

InterfaceDDGMutationTask::InterfaceDDGMutationTask( core::Size seqpos, bool complex, core::chemical::AA new_aa ) :
	seqpos_( seqpos ),
	wt_( false ),
	complex_( complex ),
	new_aa_( new_aa )
{}

InterfaceDDGMutationTask::~InterfaceDDGMutationTask() {}

basic::datacache::CacheableDataOP
InterfaceDDGMutationTask::clone() const {
	return basic::datacache::CacheableDataOP( new InterfaceDDGMutationTask( *this ));
}

/// @brief should the complex be modeled or should it be separated?
bool InterfaceDDGMutationTask::complex() const { return complex_; }

bool InterfaceDDGMutationTask::unininitialized() const { return seqpos_ == 0; }

/// @brief Does this mutation perscribe to maintain the wild type amino acid?
bool InterfaceDDGMutationTask::wt() const { return wt_; }

/// @brief Return the pose-indexed sequence position that should be mutated
core::Size
InterfaceDDGMutationTask::seqpos() const {
	return seqpos_;
}

core::chemical::AA
InterfaceDDGMutationTask::new_aa() const {
	return new_aa_;
}

////////////////// InterfaceDDGBindInnerJob /////////////////////////

InterfaceDDGBindInnerJob::InterfaceDDGBindInnerJob(
	std::string const & input_tag,
	core::Size nstruct
) :
	protocols::jd2::InnerJob( input_tag, nstruct ),
	wt_( true ),
	complex_( true ),
	new_aa_( core::chemical::aa_unk )
{}

InterfaceDDGBindInnerJob::~InterfaceDDGBindInnerJob() {}

// InterfaceDDGMutationTask const &
// InterfaceDDGBindInnerJob::mutation() const {
//  return mutation_;
// }
//
// void InterfaceDDGBindInnerJob::mutation( InterfaceDDGMutationTask const & mut ) {
//  mutation_ = mut;
// }

bool InterfaceDDGBindInnerJob::complex() const
{
	return complex_;
}

/// @brief Does this mutation perscribe to maintain the wild type amino acid?
bool InterfaceDDGBindInnerJob::wt() const
{
	return wt_;
}

utility::file::FileName const & InterfaceDDGBindInnerJob::input_pdb() const
{
	return input_pdb_;
}

void InterfaceDDGBindInnerJob::input_pdb( utility::file::FileName const & pdb_name )
{
	input_pdb_ = pdb_name;
}


void InterfaceDDGBindInnerJob::complex( bool setting ) {
	complex_ = setting;
}

void InterfaceDDGBindInnerJob::residue_index_description( core::pose::ResidueIndexDescriptionOP setting )
{
	res_ind_desc_ = setting;
}

void InterfaceDDGBindInnerJob::new_aa( core::chemical::AA setting )
{
	wt_ = false;
	new_aa_ = setting;
}

//core::pose::ResidueIndexDescription &
//InterfaceDDGBindInnerJob::residue_index_description()
//{
// return *res_ind_desc_;
//}

core::pose::ResidueIndexDescription const &
InterfaceDDGBindInnerJob::residue_index_description() const
{
	return *res_ind_desc_;
}

core::chemical::AA
InterfaceDDGBindInnerJob::new_aa() const
{
	return new_aa_;
}

bool
InterfaceDDGBindInnerJob::same( protocols::jd2::InnerJob const & other ) const
{
	InterfaceDDGBindInnerJob const * other_ptr = dynamic_cast< InterfaceDDGBindInnerJob const * > ( &other );
	return other_ptr && InnerJob::same( other );;
}

InterfaceDDGBindInnerJobOP InterfaceDDGBindInnerJob::reference_job() const
{
	return reference_job_;
}

void InterfaceDDGBindInnerJob::reference_job( InterfaceDDGBindInnerJobOP setting )
{
	reference_job_ = setting;
}


////////////////// InterfaceDDGBindJobInputter /////////////////////////

InterfaceDDGBindJobInputter::InterfaceDDGBindJobInputter() {}

InterfaceDDGBindJobInputter::~InterfaceDDGBindJobInputter() {}

void
InterfaceDDGBindJobInputter::pose_from_job( core::pose::Pose & pose, protocols::jd2::JobOP job )
{
	InterfaceDDGBindInnerJobCOP inner_job(
		utility::pointer::dynamic_pointer_cast< InterfaceDDGBindInnerJob const > ( job->inner_job() ) );

	if ( !job->inner_job()->get_pose() ) {
		InterfaceDDGBindInnerJobOP reference_job( inner_job->reference_job() );
		if ( ! reference_job || ! reference_job->get_pose() ) {
			utility::file::FileName const & fname( reference_job ? reference_job->input_pdb() : inner_job->input_pdb() );
			core::import_pose::pose_from_file( pose, fname , core::import_pose::PDB_file);
			core::pose::PoseCOP pose_cop( core::pose::PoseOP( new core::pose::Pose( pose ) ));
			if ( reference_job ) {
				protocols::jd2::JobOP dummy_job( new protocols::jd2::Job( reference_job, 1 ));
				load_pose_into_job( pose_cop, dummy_job ); // point both inner_job and refernce_job at the same Pose
			}
			load_pose_into_job( pose_cop, job );       // to slighty reduce the number of Poses that are held in memory
		} else {
			load_pose_into_job( reference_job->get_pose(), job );
			pose = *reference_job->get_pose();
		}
	} else {
		pose = *(job->inner_job()->get_pose());
	}

	// duplicating vikram's code.
	// if ( option[in::auto_setup_metals].user()
	//   && (option[in::metals_distance_constraint_multiplier]() > 1.0e-10 ||
	//   option[in::metals_angle_constraint_multiplier]() > 1.0e-10) ) {
	//
	//  //If the user has specified the auto_setup_metals option, we need to add back the metal constraints.
	//  TR << "setting up metal constraints for saved pose copy " << job->input_tag() << std::endl;
	//  core::util::auto_setup_all_metal_constraints( pose,
	//   option[in::metals_distance_constraint_multiplier](),
	//   option[in::metals_angle_constraint_multiplier]() );
	// }

	// resolve the index of the residue that's changing, and store that data in the inner job
	core::Size resid = inner_job->residue_index_description().resolve_index( pose );
	InterfaceDDGMutationTask mut;
	if ( ! inner_job->wt() ) {
		mut = InterfaceDDGMutationTask( resid, inner_job->complex(), inner_job->new_aa() );
	} else {
		mut = InterfaceDDGMutationTask( resid, inner_job->complex() );
	}

	pose.data().set( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION, mut.clone() );
}

void
InterfaceDDGBindJobInputter::fill_jobs( protocols::jd2::JobsContainer & jobs )
{
	basic::Tracer TR( "apps.pilot.andrew.InterfaceDDGBindJobInputter" );

	jobs.clear();

	std::list< InterfaceDDGBindInnerJobOP > inner_jobs;
	core::Size const nstruct( get_nstruct() );

	// open -interface_ddg::jobs file
	if ( basic::options::option[ basic::options::OptionKeys::interface_ddg::jobs ].user() ) {

		std::string jobfname( basic::options::option[ basic::options::OptionKeys::interface_ddg::jobs ]() );
		utility::io::izstream ddg_jobs( jobfname );
		if ( ! ddg_jobs.good() ) {
			// could not open jobs file
			throw utility::excn::EXCN_Msg_Exception( "Could not open interface ddg jobs file: '" + jobfname );
		}

		typedef std::pair< std::string, std::string > PDBAndSeqposPair;
		utility::vector1< PDBAndSeqposPair > pdbs_in_order;
		std::map< PDBAndSeqposPair, utility::vector1< input_mutation > > mutations_for_pdb;
		core::Size count_line( 0 );
		while ( ddg_jobs.good() ) {
			std::string line, pdb, chain, resstring;
			char newaa;
			ddg_jobs.getline( line );
			++count_line;
			if ( ddg_jobs.good() ) {
				if ( line == "" ) continue;
				std::stringstream sstream( line );
				if ( sstream.peek() == '#' ) continue;
				sstream >> pdb;
				if ( ! sstream.good() ) {
					// throw an exception
					throw utility::excn::EXCN_Msg_Exception( "Problem reading '" + jobfname + "' on line " +
						utility::to_string( count_line ) + ", which reads:\n" + line );
				}
				sstream >> chain;
				if ( ! sstream.good() ) {
					// throw an exception
					throw utility::excn::EXCN_Msg_Exception( "Problem reading '" + jobfname + "' on line " +
						utility::to_string( count_line ) + ", which reads:\n" + line );
				}
				if ( chain.size() != 1 ) {
					throw utility::excn::EXCN_Msg_Exception( "Expected to read a chain identifier (or an underscore) "
						"after reading the pdb name in file '" + jobfname + "' on line " +
						utility::to_string( count_line ) + ", which reads:\n" + line );
				}
				sstream >> resstring;
				if ( ! sstream.good() ) {
					// throw an exception
					throw utility::excn::EXCN_Msg_Exception( "Problem reading '" + jobfname + "' on line " +
						utility::to_string( count_line ) + ", which reads:\n" + line );
				}
				sstream >> newaa;
				if ( ! sstream.good() ) {
					// throw an exception
					throw utility::excn::EXCN_Msg_Exception( "Problem reading '" + jobfname + "' on line " +
						utility::to_string( count_line ) + ", which reads:\n" + line );
				}
				if ( ! core::chemical::oneletter_code_specifies_aa( newaa ) ) {
					throw utility::excn::EXCN_Msg_Exception( "Could not interpret the input character '" +
						utility::to_string( newaa ) + "' as an amino acid 1-letter code.\n" + "Problem encountered on line " +
						utility::to_string( count_line ) + " of file " + jobfname + " which reads " + line );
				}

				std::string pos_string = utility::to_string( chain ) + utility::to_string( resstring );
				PDBAndSeqposPair pdb_pos_string = std::make_pair( pdb, pos_string );

				if ( mutations_for_pdb.find( pdb_pos_string ) == mutations_for_pdb.end() ) {
					// first encounter with this pdb
					pdbs_in_order.push_back( pdb_pos_string );
				}
				int resid; char insertion_code;
				try {
					core::pose::parse_PDBnum_icode( resstring, jobfname, count_line, resid, insertion_code );
				} catch ( utility::excn::EXCN_Msg_Exception e ) {
					throw utility::excn::EXCN_Msg_Exception( "Could not interpret the string '" + resstring + "' as a residue identifier\n" + e.msg() );
				}
				input_mutation mutation;
				mutation.resind = core::pose::ResidueIndexDescriptionFromFileOP(
					new core::pose::ResidueIndexDescriptionFromFile( jobfname, count_line, chain[0], resid, insertion_code ));
				mutation.mutaa  = core::chemical::aa_from_oneletter_code( newaa );

				mutations_for_pdb[ pdb_pos_string ].push_back( mutation );
			}
		}

		std::map< std::string, std::string > base_names;
		for ( core::Size ii = 1; ii <= pdbs_in_order.size(); ++ii ) {
			std::string const & iipdb = pdbs_in_order[ii].first;
			utility::file::FileName iiname( iipdb );
			std::string iibase = iiname.base();
			if ( base_names.find( iibase ) == base_names.end() ) {
				base_names[ iibase ] = pdbs_in_order[ii].first;
			} else {
				std::string existing_pdb = base_names[ iibase ];
				if ( existing_pdb != iipdb ) {
					throw utility::excn::EXCN_Msg_Exception( "Two different PDBs, '" + existing_pdb + "' and '" +
						iipdb + " have the same base name, '" + iibase + "' and so the output structure names will collide." +
						"\nPlease rename one of them." );
				}
			}
		}

		// ok, now construct the inner jobs
		for ( core::Size ii = 1; ii <= pdbs_in_order.size(); ++ii ) {
			std::string const & iipdb = pdbs_in_order[ ii ].first;
			utility::file::FileName iipdb_fname( iipdb );
			std::string iipdb_base = iipdb_fname.base();
			std::string const & iipos = pdbs_in_order[ ii ].second;
			TR << "Creating mutations for " << iipdb << " at position " << iipos << std::endl;
			utility::vector1< input_mutation > const & ii_mutations = mutations_for_pdb[ pdbs_in_order[ ii ] ];
			InterfaceDDGBindInnerJobOP wt_job_complex( new InterfaceDDGBindInnerJob( iipdb_base + "_" + iipos + "_wt_com", nstruct ));
			wt_job_complex->complex( true );
			wt_job_complex->input_pdb( iipdb_fname );
			wt_job_complex->residue_index_description( ii_mutations[ 1 ].resind );
			InterfaceDDGBindInnerJobOP wt_job_separated( new InterfaceDDGBindInnerJob(  iipdb_base + "_" + iipos + "_wt_sep", nstruct ));
			wt_job_separated->complex( false );
			wt_job_separated->reference_job( wt_job_complex );
			wt_job_separated->residue_index_description( ii_mutations[ 1 ].resind );
			inner_jobs.push_back( wt_job_complex );
			inner_jobs.push_back( wt_job_separated );
			for ( core::Size jj = 1; jj <= ii_mutations.size(); ++jj ) {
				TR << "  Mutation " << iipdb << " pos " << iipos << " " << ii_mutations[jj].mutaa << std::endl;
				for ( int kk = 0; kk <= 1; ++kk ) {
					std::string complex = kk ? "_com" : "_sep";
					InterfaceDDGBindInnerJobOP jj_job( new InterfaceDDGBindInnerJob(
						iipdb_base + "_" + iipos + "_mut" +
						utility::to_string( core::chemical::oneletter_code_from_aa( ii_mutations[jj].mutaa ))
						+ complex, nstruct ) );
					jj_job->complex( kk );
					jj_job->residue_index_description( ii_mutations[ jj ].resind );
					jj_job->reference_job( wt_job_complex );
					jj_job->new_aa( ii_mutations[ jj ].mutaa );
					inner_jobs.push_back( jj_job );
				}
			}
		}
	} else {
		utility::excn::EXCN_Msg_Exception( "InterfaceDDGJobInputter requires the interface_ddg::jobs flag" );
	}

	for ( std::list< InterfaceDDGBindInnerJobOP >::const_iterator
			iter = inner_jobs.begin(); iter != inner_jobs.end(); ++iter ) {
		for ( core::Size ii = 1; ii <= nstruct; ++ii ) {
			jobs.push_back( protocols::jd2::JobOP( new protocols::jd2::Job( *iter, ii )));
		}
	}

	// else if ( basic::options::option[ interface_ddg::pdb_file ].user() &&
	//   basic::options::option[ interface_ddg::mutation ].user() ) {
	//  // a single job
	//
	// }
}

protocols::jd2::JobInputterInputSource::Enum
InterfaceDDGBindJobInputter::input_source() const { return protocols::jd2::JobInputterInputSource::PDB_FILE; }

////////////////// InterfaceDDGMover /////////////////////////

void InterfaceDDGMover::apply( core::pose::Pose & pose ) {

	core::Size seqpos        = determine_central_seqpos( pose );
	bool wt                  = determine_if_wt( pose );
	core::chemical::AA mutaa = determine_mutation( pose );
	bool complex             = determine_if_complex( pose );

	CloseContactWithResidue selector( seqpos );
	core::select::residue_selector::ResidueSubset subset = selector.apply( pose );

	if ( ! complex ) separate_complex( pose );
	mutate_and_relax( pose, seqpos, wt, mutaa, subset );

}

protocols::moves::MoverOP InterfaceDDGMover::clone() const
{
	return protocols::moves::MoverOP( new InterfaceDDGMover( *this ));
}

std::string InterfaceDDGMover::get_name() const
{
	return "InterfaceDDGMover";
}

protocols::moves::MoverOP InterfaceDDGMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new InterfaceDDGMover );
}

void InterfaceDDGMover::sfxn( core::scoring::ScoreFunctionOP setting )
{
	sfxn_ = setting;
}


core::Size InterfaceDDGMover::determine_central_seqpos( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) );
	InterfaceDDGMutationTaskCOP ddg_task( utility::pointer::dynamic_pointer_cast< InterfaceDDGMutationTask const > (
		pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) ));
	return ddg_task->seqpos();
}

bool
InterfaceDDGMover::determine_if_wt( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) );
	InterfaceDDGMutationTaskCOP ddg_task( utility::pointer::dynamic_pointer_cast< InterfaceDDGMutationTask const > (
		pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) ));
	return ddg_task->wt();
}


core::chemical::AA
InterfaceDDGMover::determine_mutation( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) );
	InterfaceDDGMutationTaskCOP ddg_task( utility::pointer::dynamic_pointer_cast< InterfaceDDGMutationTask const > (
		pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) ));
	return ddg_task->new_aa();
}

bool InterfaceDDGMover::determine_if_complex( core::pose::Pose const & pose )
{
	assert( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) );
	InterfaceDDGMutationTaskCOP ddg_task( utility::pointer::dynamic_pointer_cast< InterfaceDDGMutationTask const > (
		pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_DDG_MUTATION ) ));
	return ddg_task->complex();
}

void InterfaceDDGMover::separate_complex( core::pose::Pose & pose )
{
	// uhhh... let's just separate jump 1 for now.
	protocols::rigid::RigidBodyTransMover trans( pose, 1 );
	trans.step_size( 3000 );
	trans.apply( pose );
}

void InterfaceDDGMover::mutate_and_relax(
	core::pose::Pose & pose,
	core::Size seqpos,
	bool wt,
	core::chemical::AA mutaa,
	core::select::residue_selector::ResidueSubset const & subset
)
{
	assert( wt || mutaa <= core::chemical::num_canonical_aas );
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( subset[ ii ] || ii == seqpos ) {
			if ( ii == seqpos && ! wt ) {
				utility::vector1< bool > aa_vect( 20, false );
				aa_vect[ mutaa ] = true;
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( aa_vect );
			} else {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
			}
		} else {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}

	core::pack::pack_rotamers( pose, *sfxn_, task );

}



int main( int argc, char * argv [] )
{
	try {
		NEW_OPT( interface_ddg::jobs, "The list of point mutations to run through the interface_ddg_bind protocol", "" );

		devel::init( argc, argv );
		InterfaceDDGMoverOP iddgm( new InterfaceDDGMover );
		iddgm->sfxn( core::scoring::get_score_function() );

		protocols::jd2::JobDistributor::get_instance()->go( iddgm, protocols::jd2::JobInputterOP( new InterfaceDDGBindJobInputter ));

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
