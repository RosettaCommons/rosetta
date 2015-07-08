// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  A filter that, for each dimer in a pose, outputs the peptide which contributes most to the interface.

/// @author Nir London
/// @author Yuval Sedan (yuval.sedan@mail.huji.ac.il)
/// @author Orly Marcu (orly.marcu@mail.huji.ac.il)
/// @date   Jul. 30, 2014

// TODO : consider whether JD2 has infrastructure to replace PeptideDeriverJobOutputter

// PeptideDeriverFilter


// Unit headers
#include <protocols/analysis/PeptideDeriverFilter.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/util/disulfide_util.hh>
#include <numeric/xyzVector.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/DisulfideInsertionMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

//Utility headers
#include <basic/MetricValue.hh>
#include <basic/options/keys/peptide_deriver.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/ocstream.hh>
#include <utility/tag/Tag.hh>

// RosettaScripts headers
#include <protocols/rosetta_scripts/util.hh>

// External headers
#include <boost/format.hpp>
#include <boost/foreach.hpp>

// C++ headers
#include <cassert>

namespace protocols {
namespace analysis {

static thread_local basic::Tracer tracer( "protocols.analysis.PeptideDeriverFilter" );

// BEGIN PeptideDeriverOutputterContainer implementation

void PeptideDeriverOutputterContainer::push_back(PeptideDeriverOutputterOP const item) {
	list_.push_back(item);
}

void PeptideDeriverOutputterContainer::clear() {
	list_.clear();
}

void PeptideDeriverOutputterContainer::begin_structure(core::pose::Pose const & pose, std::string const &name) {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->begin_structure(pose, name);
	}
}

void PeptideDeriverOutputterContainer::chain_pair_pose_prepared(core::pose::Pose const & pose) {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->chain_pair_pose_prepared(pose);
	}
}

void PeptideDeriverOutputterContainer::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & options_string) {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->begin_receptor_partner_pair(receptor_chain_letter, partner_chain_letter, total_isc, options_string);
	}
}

void PeptideDeriverOutputterContainer::peptide_length(core::Size const pep_length) {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->peptide_length(pep_length);
	}
}

void PeptideDeriverOutputterContainer::peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const pep_start,
	core::Real const linear_isc, std::string const & disulfide_info, bool const was_cyclic_pep_modeled, core::pose::Pose const & cyclic_pose, core::Real const cyclic_isc) {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->peptide_entry(pose, entry_type, pep_start, linear_isc, disulfide_info, was_cyclic_pep_modeled, cyclic_pose, cyclic_isc);
	}
}

void PeptideDeriverOutputterContainer::end_receptor_partner_pair() {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->end_receptor_partner_pair();
	}
}

void PeptideDeriverOutputterContainer::end_structure() {
	BOOST_FOREACH( PeptideDeriverOutputterOP const outputter, list_ ) {
		outputter->end_structure();
	}
}

// END PeptideDeriverOutputterContainer implementation


// BEGIN PeptideDeriverBasicStreamOutputter implementation

PeptideDeriverBasicStreamOutputter::PeptideDeriverBasicStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;
}

void PeptideDeriverBasicStreamOutputter::begin_structure(core::pose::Pose const &, std::string const &name) {
	(*out_p_) << prefix_ << "# structure: " << name << std::endl;
}

void PeptideDeriverBasicStreamOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & options_string) {
	(*out_p_) << prefix_ << "# options: " << options_string << std::endl;
	(*out_p_) << prefix_ << "> chain_pair: receptor= " << receptor_chain_letter << " partner= " << partner_chain_letter << " total_isc= " << total_isc << std::endl;
}

void PeptideDeriverBasicStreamOutputter::peptide_length(core::Size const pep_length) {
	(*out_p_) << prefix_ << ">> peptide_length: " << pep_length << std::endl;
}

void PeptideDeriverBasicStreamOutputter::peptide_entry(core::pose::Pose const & , PeptideDeriverEntryType const entry_type, core::Size const pep_start,
	core::Real const linear_isc, std::string const & disulfide_info, bool const was_cyclic_pep_modeled, core::pose::Pose const &, core::Real const cyclic_isc) {

	(*out_p_) << prefix_
		<< "| "
		<< entry_type << " "
		<< pep_start << " "
		<< linear_isc
		<< (disulfide_info.length()>0? " " : "") << disulfide_info
		<< (was_cyclic_pep_modeled? (boost::format(" %1$.3f") % cyclic_isc).str() : "")
		<< std::endl;
}
void PeptideDeriverBasicStreamOutputter::end_receptor_partner_pair() {
	(*out_p_) << prefix_ << "# end chain pair" << std::endl;
}
void PeptideDeriverBasicStreamOutputter::end_structure() {
	(*out_p_) << prefix_ << "# end structure" << std::endl;
}

// END PeptideDeriverBasicStreamOutputter implementation

// BEGIN PeptideDeriverMarkdownStreamOutputter implementation

void PeptideDeriverMarkdownStreamOutputter::clear_buffers() {
	header_.str("");
	best_linear_peptides_.str("");
	best_cyclic_peptides_.str("");
	all_peptides_.str("");
	footer_.str("");
}

PeptideDeriverMarkdownStreamOutputter::PeptideDeriverMarkdownStreamOutputter(utility::io::orstream & out, std::string prefix) {
	out_p_ = &out;
	prefix_ = prefix;
}

void PeptideDeriverMarkdownStreamOutputter::begin_structure(core::pose::Pose const & pose, std::string const &name) {
	clear_buffers();

	header_ << prefix_ << "# Peptiderive report file " << name << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "- Chains: " << pose.conformation().num_chains() << std::endl
		<< prefix_ << "- Fold tree: " << pose.fold_tree() << std::endl
		<< prefix_ << std::endl;

	best_cyclic_peptides_ << prefix_ << "## Best cyclizable peptides for all chain pairs" << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "| Receptor | Partner | Peptide length | Position | Interface score | Disulfide info               | Cyclized interface score | Sequence             |" << std::endl
		<< prefix_ << "|----------|---------|----------------|----------|-----------------|------------------------------|--------------------------|----------------------|" << std::endl;

	best_linear_peptides_ << prefix_ << "## Best linear peptides for all chain pairs" << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "| Receptor | Partner | Peptide length | Position | Interface score | Disulfide info               | Cyclized interface score | Sequence             |" << std::endl
		<< prefix_ << "|----------|---------|----------------|----------|-----------------|------------------------------|--------------------------|----------------------|" << std::endl;

	all_peptides_ << prefix_ << "## All peptides" << std::endl
		<< prefix_ << std::endl;

}

void PeptideDeriverMarkdownStreamOutputter::end_structure() {
	footer_ << prefix_ << "*end of report*";

	( * out_p_ ) << header_.str() << std::endl << std::endl
		<< best_cyclic_peptides_.str() << std::endl << std::endl
		<< best_linear_peptides_.str() << std::endl << std::endl
		<< all_peptides_.str() << std::endl << std::endl
		<< footer_.str() << std::endl;

	clear_buffers();
}

void PeptideDeriverMarkdownStreamOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const total_isc,
	std::string const & options_string) {

	current_receptor_chain_letter_ = receptor_chain_letter;
	current_partner_chain_letter_ = partner_chain_letter;
	current_total_isc_ = total_isc;

	all_peptides_ << "- Options: " << options_string << std::endl;
}

void PeptideDeriverMarkdownStreamOutputter::peptide_length(core::Size const pep_length) {
	current_pep_length_ = pep_length;

	all_peptides_ << ( boost::format("### Receptor= %1% Partner= %2% Peptide_length= %3%")
			% current_receptor_chain_letter_
			% current_partner_chain_letter_
			% current_pep_length_ ) << std::endl
		<< prefix_ << "- Total interface score: " << current_total_isc_ << std::endl
		<< prefix_ << std::endl
		<< prefix_ << "| Position | Interface score | Disulfide info               | Cyclized interface score |" << std::endl
		<< prefix_ << "|----------|-----------------|------------------------------|--------------------------|" << std::endl;
}

void PeptideDeriverMarkdownStreamOutputter::peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const pep_start,
	core::Real const linear_isc, std::string const & disulfide_info, bool const was_cyclic_pep_modeled, core::pose::Pose const &, core::Real const cyclic_isc ) {

	std::string const cyclic_isc_string = (was_cyclic_pep_modeled? (boost::format("%1$.3f") % cyclic_isc).str() : "N/A");

	core::Size const PEPTIDE_CHAIN = 2;

	switch (entry_type) {
		case ET_BEST_LINEAR:
			best_linear_peptides_ << prefix_ << ( boost::format( "| %1$-8c | %2$-7c | %3$-14d | %4$-8d | %5$-15.3f | %6$-28s | %7$-24s | %8$-20s |" )
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_pep_length_
					% pep_start
					% linear_isc
					% disulfide_info
					% cyclic_isc_string
					% pose.chain_sequence(PEPTIDE_CHAIN) ) << std::endl;

			// the best linear peptide signals the end of a run for a certain chain pair in certain roles
			// add a separator at the end of the table for the 'general' peptides
			all_peptides_ << prefix_ << "---" << std::endl
				<< prefix_ << std::endl;
			break;

		case ET_BEST_CYCLIC:
			best_cyclic_peptides_ << prefix_ << ( boost::format( "| %1$-8c | %2$-7c | %3$-14d | %4$-8d | %5$-15.3f | %6$-28s | %7$-24s | %8$-20s |" )
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_pep_length_
					% pep_start
					% linear_isc
					% disulfide_info
					% cyclic_isc_string
					% pose.chain_sequence(PEPTIDE_CHAIN) ) << std::endl;

			break;

		case ET_GENERAL:
			all_peptides_ << prefix_ << ( boost::format( "| %1$-8d | %2$-15.3f | %3$-28s | %4$-24s |" )
					% pep_start
					% linear_isc
					% disulfide_info
					% cyclic_isc_string ) << std::endl;
			break;

		default: // should never happen; means code wasn't written well
			utility_exit_with_message("Invalid enum value provided for entry_type");
	}
}


void PeptideDeriverMarkdownStreamOutputter::end_receptor_partner_pair() {
	(*out_p_).flush();
}

// END PeptideDeriverMarkdownStreamOutputter implementation


// BEGIN PeptideDeriverPoseOutputter implementation

void PeptideDeriverPoseOutputter::output_pose( core::pose::Pose & pose, std::string const & pose_name ) {
	(*scorefxn_)(pose);
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobOP current_job( protocols::jd2::JobDistributor::get_instance()->current_job() );
		protocols::jd2::JobDistributor::get_instance()->job_outputter()->other_pose( current_job, pose, pose_name );
	}
	else {
		std::string const file_name(pose_name + ".pdb");
	    utility::io::ozstream out_stream( file_name );
		core::io::pdb::FileData::dump_pdb( pose, out_stream );
	}
}

PeptideDeriverPoseOutputter::PeptideDeriverPoseOutputter( bool const is_dump_best_peptide_pose,
		bool const is_dump_prepared_pose, bool const is_dump_cyclic_poses,
		core::scoring::ScoreFunctionCOP scorefxn) :
		is_dump_best_peptide_pose_(is_dump_best_peptide_pose),
		is_dump_prepared_pose_(is_dump_prepared_pose),
		is_dump_cyclic_poses_(is_dump_cyclic_poses),
		scorefxn_(scorefxn) { }

void PeptideDeriverPoseOutputter::chain_pair_pose_prepared(core::pose::Pose const & pose) {
	// note: we use the term 'prepared' and 'minimized' interchangeably
	//       anticipating that this step mgith involve other things
	//       in the meantime, it makes sense to talk about the complex just before
	//       peptiderive goes over it

	// save a copy of the prepared pose aside, but we output only on begin_receptor_partner_pair()
	// since we want to include chain letter in output, and we don't get them from here
	current_chain_pair_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	is_chain_pair_new_ = true;
}

void PeptideDeriverPoseOutputter::begin_receptor_partner_pair(char const receptor_chain_letter,
	char const partner_chain_letter, core::Real const,
	std::string const &) {
	current_receptor_chain_letter_ = receptor_chain_letter;
	current_partner_chain_letter_ = partner_chain_letter;

	// NOTE : this assumes the PeptideDeriver uses the same prepared
	//        structure when switching roles (receptor and partner)
	//        between partners in a chain pair, in case both roles are
	//        indeed evaluated.
	if ( is_dump_prepared_pose_ && is_chain_pair_new_ ) {
		std::string pose_name = ( boost::format("%1%%2%")
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_receptor_chain_letter_ : current_partner_chain_letter_  )
			% ( (current_receptor_chain_letter_ < current_partner_chain_letter_ )? current_partner_chain_letter_  : current_receptor_chain_letter_ ) ).str();

		output_pose( *current_chain_pair_pose_, pose_name );
		is_chain_pair_new_ = false;
	}
}

void PeptideDeriverPoseOutputter::peptide_length(core::Size const pep_length) {
	current_peptide_length_ = pep_length;
}

void PeptideDeriverPoseOutputter::peptide_entry(core::pose::Pose const & pose, PeptideDeriverEntryType const entry_type, core::Size const,
	core::Real const, std::string const & disulfide_info, bool const was_cyclic_pep_modeled, core::pose::Pose const & cyclic_pose, core::Real const) {

	if (is_dump_best_peptide_pose_ && (entry_type == ET_BEST_LINEAR || (entry_type == ET_BEST_CYCLIC && was_cyclic_pep_modeled))) {
		std::string const lin_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_%4%_linear_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_
					% (entry_type == ET_BEST_LINEAR? "linear" : "cyclic") ).str();
		core::pose::Pose temp_pose(pose);
		output_pose(temp_pose, lin_pose_name);

		if ( was_cyclic_pep_modeled ) {
			std::string const cyc_pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_best_%4%_cyclic_peptide")
						% current_receptor_chain_letter_
						% current_partner_chain_letter_
						% current_peptide_length_
						% (entry_type == ET_BEST_LINEAR? "linear" : "cyclic") ).str();
			core::pose::Pose temp_cyclic_pose(cyclic_pose);
			output_pose(temp_cyclic_pose, cyc_pose_name);
		}

	}
	else if (is_dump_cyclic_poses_ && entry_type == ET_GENERAL && was_cyclic_pep_modeled) {
		std::string const pose_name = ( boost::format("receptor%1%_partner%2%_%3%aa_%4%_cyclic_peptide")
					% current_receptor_chain_letter_
					% current_partner_chain_letter_
					% current_peptide_length_
					% disulfide_info ).str();
		core::pose::Pose temp_pose(cyclic_pose);
		output_pose(temp_pose, pose_name);
	}

}

// END PeptideDeriverPoseOutputter implementation


// -------------------------------------------------------------------------------- //
// ---                                                                          --- //
// ---                 PeptideDeriverFilter implementation                      --- //
// ---                                                                          --- //
// -------------------------------------------------------------------------------- //

// NOTE : based on LoopMover version
PeptideDeriverFilter & PeptideDeriverFilter::operator=( PeptideDeriverFilter const & rhs ) {
	if (this == &rhs) return *this;

	// NOTE : if Filter will have an operator=, call it:
	// Filter::operator=(rhs);
	PeptideDeriverFilter::assign(*this, rhs);
	return *this;
}

void PeptideDeriverFilter::parse_options() {
	// read user-defined options

	// the default value for pep_lengths is in options_rosetta.py
	// NOTE : the leading :: is to disambiguate peptide_deriver from the protocols:: namespace with the top-level one
	set_pep_lengths(basic::options::option[basic::options::OptionKeys::peptide_deriver::pep_lengths]());

	set_is_skip_zero_isc(basic::options::option[basic::options::OptionKeys::peptide_deriver::skip_zero_isc]());
	set_is_dump_peptide_pose(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_peptide_pose]());
	set_is_dump_report_file(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_report_file]());
	set_is_report_gzip(basic::options::option[basic::options::OptionKeys::peptide_deriver::report_gzip]());
	set_is_dump_prepared_pose(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_prepared_pose]());
	set_is_dump_cyclic_poses(basic::options::option[basic::options::OptionKeys::peptide_deriver::dump_cyclic_poses]());
	set_is_do_minimize(basic::options::option[basic::options::OptionKeys::peptide_deriver::do_minimize]());
	set_optimize_cyclic_threshold(basic::options::option[basic::options::OptionKeys::peptide_deriver::optimize_cyclic_threshold]());
	set_report_format( PeptideDeriverFilter::parse_report_format_string(basic::options::option[basic::options::OptionKeys::peptide_deriver::report_format]()) );

	utility::vector1<char> restrict_receptors_to_chains;

	BOOST_FOREACH(std::string const chain_string, basic::options::option[basic::options::OptionKeys::peptide_deriver::restrict_receptors_to_chains]()){
		assert(chain_string.size() == 1);
		restrict_receptors_to_chains.push_back(chain_string[0]);
	}

	utility::vector1<char> restrict_partners_to_chains;
	BOOST_FOREACH(std::string const chain_string, basic::options::option[basic::options::OptionKeys::peptide_deriver::restrict_partners_to_chains]()){
		assert(chain_string.size() == 1);
		restrict_partners_to_chains.push_back(chain_string[0]);
	}

	set_restrict_receptors_to_chains(restrict_receptors_to_chains);
	set_restrict_partners_to_chains(restrict_partners_to_chains);
}

PeptideDeriverFilter::PeptideDeriverFilter() :
		Filter("PeptideDeriverFilter") {
	// this is currently redundant, since we always read the options (parse_options())
	// and this is the default specified in options_rosetta.py, but I left it here
	// so that if options are moved out of the ctor, it would be clear that it's probably
	// a good idea to do this!
	pep_lengths_.push_back(10);

	parse_options();

	// the deriver uses the default score function
	scorefxn_deriver_ = core::scoring::get_score_function();

	// the minimizer uses the default score function, but with a full weight
	// on coordinate constraints (see PeptideDeriverFilter::minimize())
	scorefxn_minimizer_ = core::scoring::get_score_function();
	scorefxn_minimizer_->set_weight(core::scoring::coordinate_constraint, 1.0);
}

void
PeptideDeriverFilter::assign(PeptideDeriverFilter & lhs, PeptideDeriverFilter const & rhs) {
	lhs.scorefxn_deriver_          = rhs.scorefxn_deriver_         ;
	lhs.scorefxn_minimizer_        = rhs.scorefxn_minimizer_       ;
	lhs.pep_lengths_               = rhs.pep_lengths_              ;
	lhs.is_skip_zero_isc_          = rhs.is_skip_zero_isc_         ;
	lhs.is_dump_peptide_pose_      = rhs.is_dump_peptide_pose_     ;
	lhs.is_dump_report_file_       = rhs.is_dump_report_file_      ;
	lhs.is_report_gzip_            = rhs.is_report_gzip_           ;
	lhs.is_dump_prepared_pose_     = rhs.is_dump_prepared_pose_    ;
	lhs.is_dump_cyclic_poses_      = rhs.is_dump_cyclic_poses_     ;
	lhs.is_do_minimize_            = rhs.is_do_minimize_           ;
	lhs.optimize_cyclic_threshold_ = rhs.optimize_cyclic_threshold_;
	lhs.report_format_             = rhs.report_format_            ;
}

PeptideDeriverFilter::PeptideDeriverFilter(PeptideDeriverFilter const & rhs) :
		Filter(rhs) {
	PeptideDeriverFilter::assign(*this, rhs);
}

/// Minimize the pose, while using native constraints.
void
PeptideDeriverFilter::minimize(core::pose::Pose & pose) const {

	// set native constraints prior to minimization

	core::scoring::constraints::ConstraintSetOP cst_set(
			new core::scoring::constraints::ConstraintSet());
	core::scoring::func::HarmonicFuncOP spring(
			new core::scoring::func::HarmonicFunc(
				0, // mean
				1  // std. dev
			));
	core::conformation::Conformation const & conformation(pose.conformation());
	for (core::Size i = 1; i <= pose.total_residue(); ++i) {
		if ( ! pose.residue(i).is_protein() ) {
			continue;
		}
		core::id::AtomID CAi( pose.residue(i).atom_index("CA"), i );
		cst_set->add_constraint( core::scoring::constraints::ConstraintCOP(
				new core::scoring::constraints::CoordinateConstraint( CAi, CAi,
					conformation.xyz( CAi ), spring ) ) );
	}
	pose.constraint_set( cst_set );

	// minimize bb and sc
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_bb( true );
	movemap->set_chi( true );

	// TODO : consider using whatever is mentioned in -run:min_type and -run:min_tolerance
	//        and maybe change the defaults to this for Peptiderive app
	protocols::simple_moves::MinMover minimizer( movemap, scorefxn_minimizer_,
			"dfpmin_armijo_atol", 0.01, true /*nb_list*/ );

	// minimize pose
	minimizer.apply( pose );

}

// this filter doesn't do any filtering, only reporting
bool
PeptideDeriverFilter::apply(
		core::pose::Pose const & // pose - unused
	) const {

	// we do nothing here
	// the interesting stuff happens in PeptideDeriverFilter::report()

	return true;
}

PeptideDeriverReportFormat
PeptideDeriverFilter::parse_report_format_string(
		std::string const value) {
	if ( value == "markdown" ) {
		return PRF_MARKDOWN;
	}
	if ( value == "basic" ) {
		return PRF_BASIC;
	}
	throw utility::excn::EXCN_KeyError( "PeptideDeriverFilter::parse_report_format_string(std::string const value): value must be either 'markdown' or 'basic'." );
}

/// Verifies the pose has at least two chains. Minimizes the pose and detects
/// disulfides.
///
/// @param pose the pose to prepare.
void
PeptideDeriverFilter::prepare_pose(
		PeptideDeriverOutputter & ,
		core::pose::Pose & pose) const {
	// check number of chains
	if (pose.conformation().num_chains() < 2) {
		throw utility::excn::EXCN_Msg_Exception(
				"Input PDB contains less than 2 chains.");
	}

	if ( is_do_minimize_ ) {
		// minimize pose
		tracer << "Start: minimization" << std::endl;
		minimize(pose);
		tracer << "End: minimization" << std::endl;
	}

	// detect SS
	tracer << "Start: detect disulfides" << std::endl;
	pose.conformation().detect_disulfides();
	tracer << "End: detect disulfides" << std::endl;
}

// the report() method is where PeptideDeriver analyzes the pose
// and spits out information about derived peptides
void
PeptideDeriverFilter::report(std::ostream & out, core::pose::Pose const & orig_pose) const {

	tracer << "Start PeptideDeriver report" << std::endl;

	// determine what to output, where and how
	utility::io::ozstream file_out;
	utility::io::ocstream log_out( out );

	std::string const pose_name = ( protocols::jd2::jd2_used()? protocols::jd2::JobDistributor::get_instance()->current_output_name() : "" );

	//std::string prefix;
	if ( is_dump_report_file_ ) {
		std::string path = ( basic::options::option[ basic::options::OptionKeys::out::path::all ].user()? basic::options::option[ basic::options::OptionKeys::out::path::all ]().path() : "" );
		std::string file_name = path + ( pose_name.length()? pose_name + "." : "" ) + "peptiderive.txt";
		if ( is_report_gzip_ ) {
			file_name += ".gz";
		}
		tracer << "Setting up report file output to " << file_name << std::endl;
		file_out.open( file_name );
	}

	PeptideDeriverOutputterContainer output;

	utility::io::orstream & output_stream = ( is_dump_report_file_?
                (utility::io::orstream &)file_out :
                (utility::io::orstream &)log_out );

 	// this is simpler than creating a full-fledged factory. Sorry, design pattern people. --yuvals
	switch ( report_format_ ) {
		case PRF_MARKDOWN:
			output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverMarkdownStreamOutputter( output_stream , /*prefix=*/"" ) ) );
			break;
		case PRF_BASIC:
			output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverBasicStreamOutputter( output_stream , /*prefix=*/"" ) ) );
			break;
	}


	if ( is_dump_peptide_pose_ || is_dump_prepared_pose_ || is_dump_cyclic_poses_ ) {
		tracer << "Setting up pose output ("
			<< "peptides=" << is_dump_peptide_pose_ << ","
			<< "prepared=" << is_dump_prepared_pose_ << ","
			<< "cyclic=" << is_dump_cyclic_poses_ << ")"
			<< std::endl;
		output.push_back( PeptideDeriverOutputterOP( new PeptideDeriverPoseOutputter( is_dump_peptide_pose_,
				is_dump_prepared_pose_, is_dump_cyclic_poses_, scorefxn_deriver_ ) ) );
	}

	// prepare the PDB

	tracer << "Begin structure " << pose_name << std::endl;
	output.begin_structure(orig_pose, pose_name);

	// derive

	// if user didn't specify receptor/partner chains, we consider all chains as possible
	// receptor/partner chains (respectively).
	utility::vector1<core::Size> receptor_chain_indices = get_chain_indices(orig_pose, restrict_receptors_to_chains_);
	utility::vector1<core::Size> partner_chain_indices = get_chain_indices(orig_pose, restrict_partners_to_chains_);

	// since we have the 'both_ways' argument in PeptideDeriverFilter::derive_peptide()
	// we do some effort here to call the derivation function as little as possible.
	// as there is some set up done while preparing the two-monomer pose
	// so this figures out which poses can be run with both_ways = true, and which can't
	utility::vector1<core::Size> already_done_both_ways;

	// go through all requested chains
	BOOST_FOREACH( core::Size const receptor_chain, receptor_chain_indices ){

		bool is_receptor_also_partner = partner_chain_indices.contains(receptor_chain);

		BOOST_FOREACH( core::Size const partner_chain, partner_chain_indices ){
			if ( receptor_chain == partner_chain ) {
				continue;
			}

			// when both chains should be both receptors and partners, run only once with both_ways = true
			bool both_ways = (is_receptor_also_partner && receptor_chain_indices.contains(partner_chain));

			if ( both_ways ) {
				// check we've not been here already with the roles switched (see below, where already_done_both_ways is populated)
				if ( !already_done_both_ways.contains(partner_chain) ) {
					derive_peptide( output, orig_pose, receptor_chain, partner_chain, /* both_ways = */ true);
				}
			}
			else {
				derive_peptide( output, orig_pose, receptor_chain, partner_chain, /* both_ways = */ false );
			}
		}

		// avoid deriving from this chain when the time comes and other chains ask for it to be the partner
		if ( is_receptor_also_partner ) {
			already_done_both_ways.push_back(receptor_chain);
		}
	}


	tracer << "End structure " << pose_name << std::endl;
	output.end_structure();


}

utility::vector1<core::Size>
PeptideDeriverFilter::get_chain_indices(core::pose::Pose const & pose, utility::vector1<char> const & restrict_to_chains) {
	utility::vector1<core::Size> chain_indices;

	if ( restrict_to_chains.empty() ) {
		// no restriction mean return all chains
		core::Size chain_count = pose.conformation().num_chains();
		for ( core::Size i = 1; i <= chain_count; ++i ) {
			chain_indices.push_back(i);
		}
	}
	else {
		BOOST_FOREACH( char const chain_char, restrict_to_chains ){
			chain_indices.push_back( core::pose::get_chain_id_from_chain( chain_char, pose ) );
		}
	}
	return chain_indices;
}


/// Pushes the residue indices belonging to the specified chain index into a given vector.
void
PeptideDeriverFilter::push_chain_residue_indices(core::pose::Pose const & many_chain_pose,
		core::Size const chain_index, utility::vector1<core::Size> & residue_indices ) {

	core::Size chain_begin = many_chain_pose.conformation().chain_begin( chain_index );
	core::Size chain_end = many_chain_pose.conformation().chain_end( chain_index );

	for (core::Size i = chain_begin; i <= chain_end; ++i) {
		residue_indices.push_back( i );
	}

}

/// cuts the specified chains from the original pose using core::io::pdb::pose_from_pose()
/// and calculates the complex (receptor and partner)
///
/// @param pose a many-chain pose
/// @param first_chain_index chain index of the receptor (and also partner, if both_ways is true)
/// @param second_chain_index chain index of the partner (and also receptor, if both_ways is true)
/// @param both_ways whether the two specified chains may play the role of both the partner and
///        the receptor. When false, the first chain is considered to be the receptor, and the
///        second one is considered to be the partner

void
PeptideDeriverFilter::derive_peptide(
		PeptideDeriverOutputter & output,
		core::pose::Pose const & pose,
		core::Size const first_chain_index,
		core::Size const second_chain_index,
		bool const both_ways ) const {

	// create a subpose with only the two chains
	utility::vector1<core::Size> residue_indices;

	PeptideDeriverFilter::push_chain_residue_indices( pose, first_chain_index, residue_indices );
	PeptideDeriverFilter::push_chain_residue_indices( pose, second_chain_index, residue_indices );

	core::pose::Pose chain_pair_pose;

	core::io::pdb::pose_from_pose( chain_pair_pose, pose, residue_indices );

	// prepare the chain pair pose
	prepare_pose(output, chain_pair_pose);

	// calculate the total energy of the chain pair pose
	core::Size const jump_id = chain_pair_pose.fold_tree().get_jump_that_builds_residue( chain_pair_pose.conformation().chain_begin(2) );
	core::Real const total_isc( calculate_interface_score(chain_pair_pose, jump_id) );

	tracer << "Chain pair prepared " << chain_pair_pose.pdb_info()->chain(chain_pair_pose.conformation().chain_begin(1)) << chain_pair_pose.pdb_info()->chain(chain_pair_pose.conformation().chain_begin(2)) << std::endl;
	output.chain_pair_pose_prepared(chain_pair_pose);

	// derive the peptide
	utility::vector1<core::pose::PoseOP> chains = chain_pair_pose.split_by_chain();
	core::pose::PoseOP first_chain_pose = chains[1];
	core::pose::PoseOP second_chain_pose = chains[2];

	derive_peptide(output, *first_chain_pose, *second_chain_pose, total_isc);

	if ( both_ways ) {
		derive_peptide(output, *second_chain_pose, *first_chain_pose, total_isc);
	}
}

/// Goes through the residues in partner_pose with a sliding window (the size
/// of which is determined by PeptideDeriverFilter::pep_lengths_), and
/// calculates the energetic contribution of each such linear segment. This is
/// output for each position in the partner_pose, together with additional
/// data, such as the fraction of the energetic contribution to the total
/// interface score (@see PeptideDeriverFilter::calculate_interface_score()),
/// and putative SS bridge sites).
///
/// Finally, the best position is output again, with a special mark ('best')
/// to indicate that it is the best position.
///
/// Note: for lack of a better term, we call the chain from which the
///       peptide will be derived, the partner.
void
PeptideDeriverFilter::derive_peptide(
		PeptideDeriverOutputter & output,
		core::pose::Pose const & receptor_pose,
		core::pose::Pose const & partner_pose,
		core::Real const total_isc ) const {

	core::Size receptor_start = receptor_pose.conformation().chain_begin(1);
	//core::Size receptor_end = receptor_pose.conformation().chain_end(1);

	core::Size partner_start = partner_pose.conformation().chain_begin(1);
	core::Size partner_end = partner_pose.conformation().chain_end(1);

	char receptor_chain_letter = receptor_pose.pdb_info()->chain(receptor_start);
	char partner_chain_letter = partner_pose.pdb_info()->chain(partner_start);

	std::ostringstream options_string;

	// show options which have an effect on the way we're running
	options_string << "pep_lengths= "        << utility::join(pep_lengths_, ",") << " "
				<< "is_skip_zero_isc= "      << is_skip_zero_isc_                << " "
				<< "is_do_minimize= "        << is_do_minimize_                  ;

	tracer << "Begin receptor-partner pair " << receptor_chain_letter << partner_chain_letter << std::endl;
	output.begin_receptor_partner_pair(receptor_chain_letter, partner_chain_letter, total_isc, options_string.str());

	// this is used to skip discontinuities in the chain (and not derive discontinuous peptides)
	utility::vector1<core::Size> partner_cutpoints = partner_pose.fold_tree().cutpoints();
	if ( ! partner_cutpoints.empty() ) {
		tracer << "Partner cutpoints: " << utility::join(partner_cutpoints, ",") << std::endl;
	}

	// will be applied for each peptide with significant contribution to binding
	// where residues immediately before and after can be mutated to cysteins and form a disulfide bridge
	const core::Size peptide_chain = 2;

	// TODO : perhaps we want the jump in the movemap that DisulfideInsertionMover uses to be user-defined (command-line option), to prevent the peptide from escaping the binding pocket
	protocols::simple_moves::DisulfideInsertionMoverOP disulfide_inserter( new protocols::simple_moves::DisulfideInsertionMover(peptide_chain) );

	BOOST_FOREACH( core::Size const pep_length, pep_lengths_ ){

		if (pep_length > (partner_end - partner_start + 1)) {
			tracer << "Skip peptide length " << pep_length << " as it is larger than the chain" << std::endl;
			continue;
		}
		tracer << "Peptide length " << pep_length << std::endl;
		output.peptide_length(pep_length);

		const core::Real UNLIKELY_ISC_VALUE = 2e6;
		const core::Real CHECK_CYCLIZABLE_THRESHOLD_ISC = -0.01;
		const core::Real PRACTICALLY_ZERO_ISC = 1e-7;

		// the following list of variables, which hold information about the best peptides (in terms of
		// cyclic isc or linear isc), is a bit verbose, and could perhaps be shortened by creating a
		// PeptideDeriverEntry class which holds all the information.
		core::Real lin_isc_of_best_lin = UNLIKELY_ISC_VALUE;
		core::Real cyc_isc_of_best_lin = UNLIKELY_ISC_VALUE;
		core::Size pep_start_of_best_lin = 0;
		std::string disulfide_info_of_best_lin;
		core::pose::PoseOP lin_pose_of_best_lin;
		core::pose::PoseOP cyc_pose_of_best_lin;
		bool was_best_lin_cyclic_model_created = false;

		core::Real lin_isc_of_best_cyc = UNLIKELY_ISC_VALUE;
		core::Real cyc_isc_of_best_cyc = UNLIKELY_ISC_VALUE;
		core::Size pep_start_of_best_cyc = 0;
		std::string disulfide_info_of_best_cyc;
		core::pose::PoseOP lin_pose_of_best_cyc;
		core::pose::PoseOP cyc_pose_of_best_cyc;
		bool was_best_cyc_cyclic_model_created = false;

		bool any_peptide_cyclizable = false;
		bool any_peptide_cyclic_model_created = false;
		bool current_peptide_cyclizable = false;
		bool current_peptide_cyclic_model_created = false;

		// keep track of the next chain break
		utility::vector1<core::Size>::const_iterator next_cutpoint_it = partner_cutpoints.begin();
		/*
			pseudo-code for better handling of non-protein residues :
			bool[partner_length] valid_start_points;
			foreach ( residue in partner_pose ) {
				if ( (! residue.is_protein()) || residue.is_cutpoint() ) {
					valid_start_points.set_range(residue.index-pep_length, pep_length, false); // make sure indexes work out
				}
			}
		*/

		for (core::Size pep_start = partner_start; pep_start <= partner_end - pep_length + 1; ++pep_start) {
			core::Size pep_end = pep_start + pep_length - 1;

			// TODO : skip windows with non-protein residues
			/*
				pseudo-code for better handling of non-protein residues :
				if ( ! valid_start_point[pep_start] ) continue;
			*/

			// skip sequences which have cutpoints
			// NOTE : for this to work, Rosetta must generate cutpoints where there is a chain discontinuity
			//        This is enforced by using the -in:missing_density_to_jump command line option, which is
			//        turned on by default in the Peptiderive app
			if ( (next_cutpoint_it != partner_cutpoints.end()) && (*next_cutpoint_it >= pep_start) && (*next_cutpoint_it < pep_end) ) {
				tracer << "Skip " << pep_start << "-" << *next_cutpoint_it << ", contains cutpoint" << std::endl;
				pep_start = *next_cutpoint_it; // it will be incremented further after 'continue'
				++next_cutpoint_it;
				continue;
			}

			// assemble a pose with the receptor and a cut-out peptide, starting at pep_start and ending at pep_end
			// here, the complex is the receptor and the derived peptide
			// NOTE : calculate_interface_score() needs to know the jump number, and for that we need to know the residue of the jump
			//        options were to use receptor_peptide_pose->fold_tree().get_residue_edge(receptor_pose.total_residue() + 1), but we would have to
			//        know whether to look at the Edge's start() or stop() for the jump, and the depends on how build_receptor_peptide_pose() builds
			//        the pose. So it seemed preferable to just output it from build_receptor_peptide_pose().
			core::Size linear_jump_id;
			core::pose::PoseOP receptor_peptide_pose = build_receptor_peptide_pose(receptor_pose, partner_pose, pep_start, pep_end, linear_jump_id);

			// evaluate the newly created pose
			const core::Real linear_isc = calculate_interface_score(*receptor_peptide_pose, linear_jump_id);

			// skip zero linear_isc if requested
			if ((is_skip_zero_isc_) && (fabs(linear_isc) < PRACTICALLY_ZERO_ISC)) {
				tracer << "Skip " << pep_start << "-" << pep_end << ", zero (linear) isc" << std::endl;
				// we skipped a couple of residues, so make sure that the next_cutpoint_it is updated
				while ( (next_cutpoint_it != partner_cutpoints.end()) && (*next_cutpoint_it < pep_end) ) {
					++next_cutpoint_it;
				}
				pep_start = pep_end;
				continue;
			}

			// This section checks and prints out closability of the peptide by a proposed putative wrapper disulphide bridge
			core::Size n_putative_cyd = pep_start - 1;
			core::Size c_putative_cyd = pep_end + 1;
			std::stringstream disulfide_info("");

			current_peptide_cyclizable = false;
			current_peptide_cyclic_model_created = false;

			core::Real cyclic_isc = UNLIKELY_ISC_VALUE;
			core::pose::PoseOP receptor_cyclic_peptide_pose = receptor_peptide_pose; // NOTE : this is a dummy value so that we don't dereference NULL; see below
			core::pose::PoseOP receptor_pre_cyclization_peptide_pose = receptor_peptide_pose; // NOTE : this is a dummy value so that we don't dereference NULL; see below
			if ((linear_isc<CHECK_CYCLIZABLE_THRESHOLD_ISC) && (n_putative_cyd>=partner_start) && (c_putative_cyd<=partner_end)) {
				protocols::simple_moves::DisulfideCyclizationViability cyclization_viable(disulfide_inserter->determine_cyclization_viability(partner_pose, n_putative_cyd, c_putative_cyd));
				current_peptide_cyclizable = ( (cyclization_viable == protocols::simple_moves::DCV_CYCLIZABLE) || (cyclization_viable == protocols::simple_moves::DCV_ALREADY_CYCLIZED) );

				if (current_peptide_cyclizable) {
					any_peptide_cyclizable = true;

					disulfide_info << partner_chain_letter << "_" << n_putative_cyd << "-" << c_putative_cyd;

					// For peptides that can be closed and contribute more then third of the binding energy, mutate to cysteins and re-evaluate energy
					if ((linear_isc / total_isc) >= optimize_cyclic_threshold_) {
						// NOTE : see note on linear_jump_id
						core::Size cyclic_jump_id;
						receptor_pre_cyclization_peptide_pose = build_receptor_peptide_pose(receptor_pose, partner_pose, n_putative_cyd, c_putative_cyd, cyclic_jump_id);
						// create a copy of the linear pose which we are going to mutate and model
						// that is, for peptides for which we model the cyclic peptide, we also want to keep the expanded (+-1 residue) linear version
						receptor_cyclic_peptide_pose = core::pose::PoseOP( new core::pose::Pose( *receptor_pre_cyclization_peptide_pose ) );
						//define and create disulfide bond between the edges of the derived peptide using the DisulfideInsertionMover
						disulfide_inserter->apply(*receptor_cyclic_peptide_pose);
						// any_peptide_cyclic_model_created checks if any peptide in the pose has been cyclized
						// current_peptide_cyclic_model_created checks if the current peptide has been cyclized
						any_peptide_cyclic_model_created = true;
						current_peptide_cyclic_model_created = true;
						cyclic_isc = calculate_interface_score(*receptor_cyclic_peptide_pose, cyclic_jump_id);
				   }
				}
			}

			// output data for each cut-out peptide
			output.peptide_entry(*receptor_peptide_pose, ET_GENERAL, pep_start, normalize_isc(linear_isc), disulfide_info.str(), /*was_cyclic_pep_modeled=*/cyclic_isc<0, *receptor_cyclic_peptide_pose, normalize_isc(cyclic_isc));

			// save details for the linear peptide that scored best
			if (linear_isc < lin_isc_of_best_lin) {
				pep_start_of_best_lin = pep_start;
				lin_isc_of_best_lin = linear_isc;
				cyc_isc_of_best_lin = (current_peptide_cyclic_model_created? cyclic_isc : UNLIKELY_ISC_VALUE);
				disulfide_info_of_best_lin = disulfide_info.str();
				lin_pose_of_best_lin = (current_peptide_cyclic_model_created? receptor_pre_cyclization_peptide_pose : receptor_peptide_pose);
				// TODO : solve the need for a dummy value for receptor_cyclic_peptide_pose
				// NOTE : receptor_peptide_pose is used here as dummy, since we can't pass NULL as reference
				cyc_pose_of_best_lin = (current_peptide_cyclic_model_created? receptor_cyclic_peptide_pose : receptor_peptide_pose);
				was_best_lin_cyclic_model_created = current_peptide_cyclic_model_created;
			}

			// save details for the cyclizable peptide that scored best
			//
			// We only model cyclic peptides out of peptides that contribute enough to the interactions energy
			// (more than the fraction specified by optimize_cyclic_threshold_).
			//
			// There are three ways to define the scoring or the set from which the 'best' peptide may be taken:
			// 1. best scoring linear peptide pose
			// 2. best scoring linear peptide pose among the cyclization candidates (if any candidates exist)
			// 3. best scoring cyclic pose (if such were modeled among the candidates)
			// Right now we merge and 2nd and 3rd into
			// 2/3. the best cyclizable pose, which is the best scoring cyclic pose, if one is avaiable, otherwise, the best linear pose among cyclization candidates
			//
			// Now, consider the following scenario:
			// Total interface score: -50
			// Position    Linear peptide interface score     Cyclic peptide interface score
			// 10          -10                                -4
			// 20          -20                                -3
			// 30          -30                                -3.5
			//
			// Let's assume all listed positions are cyclization candidates, and that they are the only cyclization candidates for this chain pair.
			// We might get different results here depending on the value of optimize_cyclic_threshold
			//
			// optimize_cyclic_threshold    Positions that produce cyclic peptides       Best cyclizable peptide       Based on
			// 0  -0.2                      10, 20, 30                                   10                            Cyclic peptide interface score
			// 0.2-0.4                      20, 30                                       30                            Cyclic peptide interface score
			// 0.4-0.6                      30                                           30                            Cyclic peptide interface score
			// 0.6-1                        N/A                                          30                            Linear peptide interface score
			//
			// Note that position 20 will never be considered best cyclizable, since if no cyclic model was produced, position 30 has a better linear peptide interface score.
			// If position 20 is compared using the cyclic peptide interface score, this means a cyclic model was produced for it, which means it would be produced for position 30 as well (has a better linear peptide interface score), and position 30 has a better cyclic interface score.
			// So, if a peptide B has both a better linear peptide interface score and a better cyclic peptide interface score than peptide A, peptide A will never be the best cyclizable peptide.

			if (current_peptide_cyclizable &&
				(
					(current_peptide_cyclic_model_created && (cyclic_isc < cyc_isc_of_best_cyc)) ||
					((!any_peptide_cyclic_model_created) && (linear_isc < lin_isc_of_best_cyc))
				) ) {
				pep_start_of_best_cyc = pep_start;
				lin_isc_of_best_cyc = linear_isc;
				cyc_isc_of_best_cyc = (current_peptide_cyclic_model_created? cyclic_isc : UNLIKELY_ISC_VALUE);
				disulfide_info_of_best_cyc = disulfide_info.str();
				lin_pose_of_best_cyc = (current_peptide_cyclic_model_created? receptor_pre_cyclization_peptide_pose : receptor_peptide_pose);
				// TODO : solve the need for a dummy value for receptor_cyclic_peptide_pose
				// NOTE : receptor_peptide_pose is used here as dummy, since we can't pass NULL as reference
				cyc_pose_of_best_cyc = (current_peptide_cyclic_model_created? receptor_cyclic_peptide_pose : receptor_peptide_pose);
				was_best_cyc_cyclic_model_created = current_peptide_cyclic_model_created;
			}
		} // for pep_start

		//since the best disulfide position is not always the position of the best linear position and to make things less confusing
		//best cyclic information will be printed in a line of its own
		tracer << "Outputting best peptides" << std::endl;
		output.peptide_entry(*lin_pose_of_best_lin, ET_BEST_LINEAR, pep_start_of_best_lin, normalize_isc(lin_isc_of_best_lin), disulfide_info_of_best_lin, was_best_lin_cyclic_model_created, *cyc_pose_of_best_lin, normalize_isc(cyc_isc_of_best_lin));
		if ( any_peptide_cyclizable ) {
			// otherwise, cyc_pose_of_best_cyc is unset
			output.peptide_entry(*lin_pose_of_best_cyc, ET_BEST_CYCLIC, pep_start_of_best_cyc, normalize_isc(lin_isc_of_best_cyc), disulfide_info_of_best_cyc, was_best_cyc_cyclic_model_created, *cyc_pose_of_best_cyc, normalize_isc(cyc_isc_of_best_cyc));
		}

	} // for each peptide_length

	tracer << "End receptor-partner pair" << std::endl;
	output.end_receptor_partner_pair();
}

/// Given a two-monomer pose and a score function - the interface score is the
/// score difference between the given complex and the unbound monomers,
/// 'unbound' meaning that they are spatially seperated.
///
/// @param pose A two-monomer pose.
/// TODO : this might have code in common with DdgScan. Consider using that.
core::Real
PeptideDeriverFilter::calculate_interface_score(core::pose::Pose const & pose, core::Size const jump_id) const {
	//make a new pose in which the monomers are far apart from one another
	core::pose::Pose unbound_pose(pose);

	core::Real const VERY_LARGE_TRANSLATION = 2e6;

	//create a Rigid body mover - the mover works on a rigid body jump
	protocols::rigid::RigidBodyTransMoverOP translate_away(
			new protocols::rigid::RigidBodyTransMover(unbound_pose, jump_id));

	//set the size of the move
	translate_away->step_size(VERY_LARGE_TRANSLATION);

	//calculate scores
	core::Real bound_energy = (*scorefxn_deriver_)(unbound_pose);  // before the move - as a complex.
	translate_away->apply(unbound_pose);
	core::Real unbound_energy = (*scorefxn_deriver_)(unbound_pose);  // after the move - unbound monomers.

	//answer
	return (bound_energy - unbound_energy);
}


/// Given a two-monomer pose and a score function - the interface score (delta between energy of complex and of the seperated monomers) is calculated for each residue.
/// This signifies the energetic contribution of each residue to the complex binding.
/// @param pose a two-monomer pose.
void
PeptideDeriverFilter::calculate_per_residue_interface_score(
		core::pose::Pose & chain_pair_pose,
		std::set<core::Size> interface_residues,
		std::ostream & report_out,
		core::Size const jump_id) const {

	//seperate the given structure in space using a translate_away mover and score each of the states
	core::Real const VERY_LARGE_TRANSLATION = 2e6;
	core::pose::Pose unbound_pair(chain_pair_pose);
	protocols::rigid::RigidBodyTransMoverOP translate_away(new protocols::rigid::RigidBodyTransMover(unbound_pair, jump_id));
	translate_away->step_size(VERY_LARGE_TRANSLATION);
	translate_away->apply(unbound_pair);
	// make sure hbonds are included in the scoring function's per residue energies
	core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn_deriver_->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
	scorefxn_deriver_->set_energy_method_options( *emopts );
	(*scorefxn_deriver_)(chain_pair_pose);
	(*scorefxn_deriver_)(unbound_pair);

	//go over each interface residue and calculate its energetical contribution to binding
	for(std::set< core::Size >::const_iterator it(interface_residues.begin()), end(interface_residues.end()); it != end; ++it) {
		 // delta between total energy in the bound and unbound states is in fact the interface score
		core::Real res_isc = chain_pair_pose.energies().residue_total_energy(*it) - unbound_pair.energies().residue_total_energy(*it);
		report_out << "residue " << *it << " in chain " << chain_pair_pose.residue(*it).chain() << " has interface score " << res_isc << std::endl;
	}

}

/// for a two-monomer pose, returns a set of interface residues
/// @param pose a two-monomer pose.
std::set<core::Size>
PeptideDeriverFilter::find_interface_residues(core::pose::Pose const & chain_pair_pose) {
	//find the interface residues using an InterfaceNeighborDefinition calculator and store them in the interface_residues set
	core::Size chain1_idx =1, chain2_idx =2;
	basic::MetricValue< std::set< core::Size > > mv_interface_set;
	core::pose::metrics::PoseMetricCalculatorOP if_residues_calculator(new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator(chain1_idx, chain2_idx));
	std::string InterfaceNeighborDefinition = "InterfaceNeighborDefinition_";
	core::pose::metrics::CalculatorFactory::Instance().register_calculator(InterfaceNeighborDefinition, if_residues_calculator);
	chain_pair_pose.metric( InterfaceNeighborDefinition, "interface_residues", mv_interface_set);
	std::set<core::Size> interface_residues = mv_interface_set.value();
	return interface_residues;
}

core::Real
PeptideDeriverFilter::normalize_isc(core::Real const isc) {
	const core::Real NEAR_ZERO_ISC = 2e-2;
	return (fabs(isc) < NEAR_ZERO_ISC? 0 : isc);
}

/// For two given full proteins poses, builds a new pose with the full receptor pose connected to the peptide found between start and end positions in the partner pose
/// @param receptor_pose the receptor chain to be maintained
/// @param partner_pose  the chain from which a peptide should be cut and connected to the receptor in a new pose
/// @param peptide_start starting position in the partner pose to cut from
/// @param peptide_end   last position in the partner pose to be cut out
core::pose::PoseOP PeptideDeriverFilter::build_receptor_peptide_pose(core::pose::Pose const & receptor_pose,
		core::pose::Pose const & partner_pose,
		core::Size peptide_start,
		core::Size peptide_end,
		core::Size & jump_id) const {
	// once created these will be the positions of the peptide start and end in the receptor_peptide_pose
	core::Size two_chain_peptide_start = receptor_pose.total_residue() + 1;
	core::Size two_chain_peptide_end = two_chain_peptide_start + peptide_end - peptide_start;

	//copy receptor and start appending the peptide from its center of mass backward and forward
	core::pose::PoseOP receptor_peptide_pose( new core::pose::Pose(receptor_pose) );
	core::Size pep_cen_mass = core::pose::residue_center_of_mass(partner_pose, peptide_start, peptide_end);

	// tracer messages might be useful for debugging, but I'd rather comment than tracer.Debug
	// from Flexpepdock's flags file; choose the residue in the receptor that is closest to the peptide center-of-mass
	// TODO : handle a case where the closest residue is not a protein -- we needn't ask for "CA"
	core::Size receptor_anchor_pos = core::pose::return_nearest_residue(receptor_pose, 1, receptor_pose.total_residue(), partner_pose.residue(pep_cen_mass).atom("CA").xyz());
	// tracer << "Appending jump from " << receptor_peptide_pose->residue(receptor_anchor_pos) << " to " << partner_pose.residue(pep_cen_mass) << std::endl;

	receptor_peptide_pose->append_residue_by_jump(partner_pose.residue(pep_cen_mass), receptor_anchor_pos, "", "", true);

	for (core::Size i = pep_cen_mass ; i > peptide_start; --i) {
		// tracer << "Prepending at " << two_chain_peptide_start << " " << partner_pose.residue(i-1) << std::endl;
		receptor_peptide_pose->conformation().safely_prepend_polymer_residue_before_seqpos(partner_pose.residue(i-1), two_chain_peptide_start, false);
	}

	for (core::Size i = pep_cen_mass ; i < peptide_end; ++i) {
		// tracer << "Appending at " << receptor_peptide_pose->total_residue() << " " << partner_pose.residue(i+1) << std::endl;
		receptor_peptide_pose->conformation().safely_append_polymer_residue_after_seqpos(partner_pose.residue(i+1), receptor_peptide_pose->total_residue(), false);
	}

	// see note for linear_jump_id
	jump_id = receptor_peptide_pose->fold_tree().get_jump_that_builds_residue(receptor_pose.total_residue() + pep_cen_mass - peptide_start + 1);

	// NOTE : when enabling -in:missing_density_to_jump (which is enabled by default for the
	//        PeptideDeriver app), chain termini which are missing hydrogens will receive a
	//        UPPER/LOWERTERM_TRUNC_VARIANT. Without removing those before trying to add the
	//        UPPER/LOWER_TERMINUS_VARIANT the following error occurs:
	//        ERROR: unable to find desired variant residue: MET:NtermTruncation MET LOWER_TERMINUS_VARIANT
	if ( receptor_peptide_pose->residue(two_chain_peptide_start).has_variant_type(core::chemical::LOWERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( *receptor_peptide_pose, core::chemical::LOWERTERM_TRUNC_VARIANT, two_chain_peptide_start );
	}
	if ( receptor_peptide_pose->residue(two_chain_peptide_end).has_variant_type(core::chemical::UPPERTERM_TRUNC_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue( *receptor_peptide_pose, core::chemical::UPPERTERM_TRUNC_VARIANT, two_chain_peptide_end );
	}

	core::pose::add_variant_type_to_pose_residue(*receptor_peptide_pose, core::chemical::LOWER_TERMINUS_VARIANT, two_chain_peptide_start);
	core::pose::add_variant_type_to_pose_residue(*receptor_peptide_pose, core::chemical::UPPER_TERMINUS_VARIANT, two_chain_peptide_end);
	// if a cysteine in the peptide was disulfide bonded to a residue that no longer exists in the protein-peptide pose, it is converted back to a regular cysteine
	for (core::Size i = two_chain_peptide_start; i <= two_chain_peptide_end ; ++i) {
		if (receptor_peptide_pose->conformation().residue(i).has_variant_type( core::chemical::DISULFIDE )) {
			core::conformation::change_cys_state(i, "CYS", receptor_peptide_pose->conformation());
		}
	}
	receptor_peptide_pose->conformation().detect_disulfides();
	return receptor_peptide_pose;
}


void
PeptideDeriverFilter::parse_my_tag( utility::tag::TagCOP tag,
			basic::datacache::DataMap & data,
			protocols::filters::Filters_map const & , // filters - unused
			protocols::moves::Movers_map const & , // movers - unused
			core::pose::Pose const & // pose - unused
		) {

	// NOTE : we rely on the fact that all the options already have values
	//        from when parse_options() was called. This way command-line
	//        defaults are also effectively RosettaScripts defaults.

	if ( tag->hasOption("pep_lengths") ) {
		set_pep_lengths(utility::string_split( tag->getOption<std::string>("pep_lengths"), ',', core::Size() ));
	}

	if ( tag->hasOption("skip_zero_isc") ) {
		set_is_skip_zero_isc(tag->getOption<bool>("skip_zero_isc"));
	}

	if ( tag->hasOption("dump_peptide_pose") ) {
		set_is_dump_peptide_pose(tag->getOption<bool>("dump_peptide_pose"));
	}

	if ( tag->hasOption("dump_cyclic_poses") ) {
		set_is_dump_cyclic_poses(tag->getOption<bool>("dump_cyclic_poses"));
	}

	if ( tag->hasOption("dump_report_file") ) {
		set_is_dump_report_file(tag->getOption<bool>("dump_report_file"));
	}

	if ( tag->hasOption("report_gzip") ) {
		set_is_report_gzip(tag->getOption<bool>("report_gzip"));
	}

	if ( tag->hasOption("dump_prepared_pose") ) {
		set_is_dump_prepared_pose(tag->getOption<bool>("dump_prepared_pose"));
	}

	if ( tag->hasOption("do_minimize") ) {
		set_is_do_minimize(tag->getOption<bool>("do_minimize"));
	}

	if ( tag->hasOption("scorefxn_deriver") ) {
		set_scorefxn_deriver( protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_deriver", data, "" ) );
	}

	if ( tag->hasOption("optimize_cyclic_threshold") ) {
		set_optimize_cyclic_threshold( tag->getOption<core::Real>( "optimize_cyclic_threshold" ) );
	}

	if ( tag->hasOption("report_format") ) {
		set_report_format( PeptideDeriverFilter::parse_report_format_string(tag->getOption<std::string>( "report_format" )) );
	}

	utility::vector1<char> restrict_receptors_to_chains(
		utility::string_split( tag->getOption<std::string>("restrict_receptors_to_chains", ""), ',', char() ));

	utility::vector1<char> restrict_partners_to_chains(
		utility::string_split( tag->getOption<std::string>("restrict_partners_to_chains", ""), ',', char() ));

	set_restrict_receptors_to_chains(restrict_receptors_to_chains);
	set_restrict_partners_to_chains(restrict_partners_to_chains);
}

} // namespace analysis
} // namespace protocols
