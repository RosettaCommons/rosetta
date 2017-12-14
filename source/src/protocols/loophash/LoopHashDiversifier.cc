// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/loophash/movers/LoopHashDiversifier.cc
/// @brief Simple mover that uses loophash to replace randomly chosen fragments in a given pose.
/// Heavy inspiration taken by LoophashMoverWrapper.
/// @author Tim Jacobs

// Unit headers
#include <protocols/loophash/LoopHashDiversifier.hh>
#include <protocols/loophash/LoopHashDiversifierCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/string_util.hh>


#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/sort_predicates.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/AtomID.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace loophash {

using core::pose::Pose;
using namespace utility;
using namespace protocols::moves;
using core::Real;
using core::Size;
using std::string;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;

static basic::Tracer TR( "protocols.loophash.LoopHashDiversifier" );

///****Creator Methods****///
// XRW TEMP std::string
// XRW TEMP LoopHashDiversifierCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return LoopHashDiversifier::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LoopHashDiversifierCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LoopHashDiversifier );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LoopHashDiversifier::mover_name()
// XRW TEMP {
// XRW TEMP  return "LoopHashDiversifier";
// XRW TEMP }
///****End Creator Methods****///

LoopHashDiversifier::~LoopHashDiversifier() = default;

LoopHashDiversifier::LoopHashDiversifier() :
	protocols::moves::Mover( LoopHashDiversifier::mover_name() ),
	library_( /* NULL */ ),
	min_inter_ss_bbrms_( 0.0 ),
	max_inter_ss_bbrms_( 100000.0 ),
	min_intra_ss_bbrms_( 0.0 ),
	max_intra_ss_bbrms_( 100000.0 ),
	min_rms_( 0.0 ),
	max_rms_( 100.0 ),
	start_res_( "2" ),
	stop_res_( "0" ),
	window_size_(4),
	max_radius_(4),
	max_struct_(10),
	num_iterations_(100),
	num_try_div_(100),
	diversify_loop_only_( true ),
	ideal_( false ),
	filter_by_phipsi_( false ),
	cenfilter_( /* NULL */ ),
	ranking_cenfilter_( /* NULL */ ),
	scorefxn_cen_cst_(/* NULL */),
	scorefxn_rama_cst_(/* NULL */)
{
	loop_sizes_.clear();
	loop_sizes_.push_back(window_size_);

	library_ = LoopHashLibraryOP( new LoopHashLibrary( loop_sizes() , 1 , 0 ) );
	library_->load_mergeddb();
	library_->mem_foot_print();
}

LoopHashDiversifier::LoopHashDiversifier(
	LoopHashLibraryOP library,
	core::Real min_inter_ss_bbrms,
	core::Real max_inter_ss_bbrms,
	core::Real min_intra_ss_bbrms,
	core::Real max_intra_ss_bbrms,
	core::Real min_rms,
	core::Real max_rms,
	std::string const & start_res,
	std::string const & stop_res,
	core::Size window_size,
	core::Size max_radius,
	core::Size max_struct,
	core::Size num_iterations,
	core::Size num_try_div,
	bool diversify_loop_only,
	bool ideal,
	bool filter_by_phipsi,
	protocols::filters::FilterOP cenfilter,
	protocols::filters::FilterOP ranking_cenfilter,
	core::scoring::ScoreFunctionOP scorefxn_cen_cst,
	core::scoring::ScoreFunctionOP scorefxn_rama_cst
) :
	protocols::moves::Mover( LoopHashDiversifier::mover_name() ),
	library_(std::move( library )),
	min_inter_ss_bbrms_( min_inter_ss_bbrms ),
	max_inter_ss_bbrms_( max_inter_ss_bbrms ),
	min_intra_ss_bbrms_( min_intra_ss_bbrms ),
	max_intra_ss_bbrms_( max_intra_ss_bbrms ),
	min_rms_( min_rms ),
	max_rms_( max_rms ),
	start_res_( start_res ),
	stop_res_( stop_res ),
	window_size_(window_size),
	max_radius_(max_radius),
	max_struct_(max_struct),
	num_iterations_(num_iterations),
	num_try_div_(num_try_div),
	diversify_loop_only_(diversify_loop_only),
	ideal_( ideal ),
	filter_by_phipsi_( filter_by_phipsi ),
	cenfilter_(std::move( cenfilter )),
	ranking_cenfilter_(std::move( ranking_cenfilter )),
	scorefxn_cen_cst_(std::move(scorefxn_cen_cst)),
	scorefxn_rama_cst_(std::move(scorefxn_rama_cst))
{
	loop_sizes_.clear();
	loop_sizes_.push_back(window_size_);

	library_ = LoopHashLibraryOP( new LoopHashLibrary( loop_sizes() , 1 , 0 ) );
	library_->load_mergeddb();
	library_->mem_foot_print();
}


void
LoopHashDiversifier::apply( Pose & pose )
{
	using namespace core::io::silent;
	runtime_assert( library_ != nullptr );
	Pose const saved_pose( pose );

	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );

	core::pose::set_ss_from_phipsi( pose );

	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	if ( scorefxn_cen_cst_ ) { simple_inserter->scorefxn_cen_cst(*scorefxn_cen_cst_);}
	if ( scorefxn_rama_cst_ ) { simple_inserter->scorefxn_rama_cst(*scorefxn_rama_cst_);}

	LoopHashSampler lsampler( library_, simple_inserter );

	//Configure iteration-indepenent loophash options

	//Max rms of decoys
	lsampler.set_min_rms( min_rms() );
	lsampler.set_max_rms( max_rms() );

	//Run Loophash
	lsampler.set_nonideal( !ideal_ );
	lsampler.set_max_radius( max_radius_ );
	lsampler.set_max_struct( max_struct_ );
	lsampler.set_filter_by_phipsi( filter_by_phipsi_ );

	core::Size start_res = core::pose::parse_resnum( start_res_, pose );
	core::Size stop_res = pose.size();
	if ( !stop_res_.empty() ) {
		stop_res = core::pose::parse_resnum( stop_res_, pose );
	}

	for ( core::Size cur_iter=1; cur_iter<=num_iterations_; ++cur_iter ) {
		core::Size cur_num_try_div = 1 ;
		bool div_success = false;
		while ( !div_success && cur_num_try_div <= num_try_div_ ) // it tries upto num_try_div_ times to diversify
				{
			//Choose a random window-size of residues to run loophash on
			core::Size lh_start = numeric::random::random_range(start_res, stop_res-window_size_+1);
			core::Size lh_stop = lh_start+window_size_-1;
			TR << "lh_start: " << lh_start << ", lh_stop: " << lh_stop << std::endl;

			bool ok_to_diversify = true;
			if ( diversify_loop_only_ ) {
				for ( core::Size res=lh_start; res<=lh_stop; ++res ) {
					if ( pose.conformation().secstruct(res) != 'L' ) {
						ok_to_diversify=false;
					}
				}
			}
			if ( !ok_to_diversify ) {
				TR<<"Randomly selected residues have non-loop, so select other residues!"<<std::endl;
				cur_num_try_div++;
				TR << "cur_num_try_div: " << cur_num_try_div << std::endl;
				continue;
			}


			lsampler.set_start_res( lh_start );
			lsampler.set_stop_res ( lh_start );


			//Determine min and max torsion RMSD based on secondary structure
			char sec_struct = pose.conformation().secstruct(lh_start);
			bool inter_ss=false;
			for ( core::Size res=lh_start+1; res<=lh_stop; ++res ) {
				if ( pose.conformation().secstruct(res) != sec_struct ) {
					inter_ss=true;
				}
			}

			if ( inter_ss ) {
				lsampler.set_min_bbrms( min_inter_ss_bbrms() );
				lsampler.set_max_bbrms( max_inter_ss_bbrms() );
			} else {
				lsampler.set_min_bbrms( min_intra_ss_bbrms() );
				lsampler.set_max_bbrms( max_intra_ss_bbrms() );
			}


			std::vector< SilentStructOP > lib_structs;
			Size starttime = time( nullptr );
			lsampler.build_structures( pose, lib_structs );
			Size endtime = time( nullptr );
			Size nstructs = lib_structs.size();
			TR << "Found " << nstructs << " alternative states in time: " << endtime - starttime << std::endl;

			//Shuffle the loophash structures
			//  numeric::random::random_permutation(lib_structs.begin(), lib_structs.end(), numeric::random::rg());

			std::vector< std::pair< Real, SilentStructOP > > cen_scored_structs;
			for ( std::vector< SilentStructOP >::const_iterator struct_it = lib_structs.begin();
					struct_it != lib_structs.end(); ++struct_it ) {
				Pose rpose;
				(*struct_it)->fill_pose( rpose );

				// apply selection criteria
				bool passed_i = cenfilter_->apply( rpose );
				if ( passed_i ) {
					core::Real score_i = ranking_cenfilter()->report_sm( rpose );
					core::io::silent::SilentFileOptions opts;
					core::io::silent::SilentStructOP new_struct = ideal_?
						core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( opts ) :
						core::io::silent::SilentStructFactory::get_instance()->get_silent_struct("binary", opts );
					new_struct->fill_struct( rpose );
					cen_scored_structs.emplace_back(-score_i,new_struct );
				}
			}

			// sort by centroid criteria
			std::sort( cen_scored_structs.begin(), cen_scored_structs.end(), utility::SortFirst<Real, SilentStructOP>() );
			all_structs_ = cen_scored_structs;
			TR << "After centroid filter: " << all_structs_.size() << " of " << cen_scored_structs.size() << " structures" << std::endl;

			if ( all_structs_.size() == 0 ) {
				TR<<"No structures survived centroid filter. Consider relaxing filters"<<std::endl;
				cur_num_try_div++;
				TR << "cur_num_try_div: " << cur_num_try_div << std::endl;
				//set_last_move_status( protocols::moves::FAIL_RETRY );
				//return;
			} else {
				div_success = true;
			}
		} //while (!div_success && cur_num_try_div <= num_try_div_)

		if ( !div_success ) {
			TR<<"diversification failed after " << num_try_div_ << " trials"<<std::endl;
			set_last_move_status( protocols::moves::FAIL_RETRY );
			return;
		}

		//Success!
		set_last_move_status(protocols::moves::MS_SUCCESS);

		// make best from list the next starting structure
		std::pair< Real, SilentStructOP > currbest = all_structs_.back();
		all_structs_.pop_back();
		TR << "Best score after round " << cur_iter << ": " << -currbest.first << std::endl;

		currbest.second->fill_pose( pose );

	}

	//Change to FA
	core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
}

//core::pose::PoseOP
//LoopHashDiversifier::get_additional_output() {
// if ( all_structs_.size() == 0)
//  return NULL;
//
// // pop best
// core::pose::PoseOP pose( new core::pose::Pose() );
// std::pair< core::Real, core::io::silent::SilentStructOP > currbest = all_structs_.back();
// all_structs_.pop_back();
// TR << "Returning score = " << -currbest.first << std::endl;
// currbest.second->fill_pose( *pose );
//
// return pose;
//}

// XRW TEMP std::string
// XRW TEMP LoopHashDiversifier::get_name() const {
// XRW TEMP  return LoopHashDiversifier::mover_name();
// XRW TEMP }

void
LoopHashDiversifier::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &filters,
	Movers_map const & /*movers*/,
	Pose const &
){
	num_iterations_ = tag->getOption< Size >( "num_iterations", 100 );
	num_try_div_ = tag->getOption< Size >( "num_try_div", 100 );

	diversify_loop_only_ = tag->getOption< bool >( "diversify_loop_only", false ); // specify in xml 'diversify_loop_only=1/0'

	min_inter_ss_bbrms_ = tag->getOption< Real >( "min_inter_ss_bbrms", 0 );
	max_inter_ss_bbrms_ = tag->getOption< Real >( "max_inter_ss_bbrms", 100000 );

	min_intra_ss_bbrms_ = tag->getOption< Real >( "min_intra_ss_bbrms", 0 );
	max_intra_ss_bbrms_ = tag->getOption< Real >( "max_intra_ss_bbrms", 100000 );

	min_rms_ = tag->getOption< Real >( "min_rms",   0.0 );
	max_rms_ = tag->getOption< Real >( "max_rms",   100.0 );

	max_radius_ = tag->getOption< Size >( "max_radius", 4 );

	max_struct_ = tag->getOption< Size >( "max_struct", 10 );

	ideal_ = tag->getOption< bool >( "ideal",  false );  // by default, assume structure is nonideal
	filter_by_phipsi_ = tag->getOption< bool >( "filter_by_phipsi", false );

	start_res_ = core::pose::get_resnum_string( tag, "start_", start_res_ );
	stop_res_ = core::pose::get_resnum_string( tag, "stop_", "" );

	window_size_ = tag->getOption< Size >( "window_size", 4 );

	if ( tag->hasOption("scorefxn_cen_cst") ) {
		std::string scorefxn_name = tag->getOption<string>( "scorefxn_cen_cst" );
		scorefxn_cen_cst_ = data.get_ptr<ScoreFunction>( "scorefxns", scorefxn_name );
	}
	if ( tag->hasOption("scorefxn_rama_cst") ) {
		std::string scorefxn_name = tag->getOption<string>( "scorefxn_rama_cst" );
		scorefxn_rama_cst_ = data.get_ptr<ScoreFunction>( "scorefxns", scorefxn_name );
	}

	//Currently use only window_size fragment sizes
	add_loop_size( window_size_ ) ;

	// path to DB -- if not specified then command-line flag is used
	library_ = LoopHashLibraryOP( new LoopHashLibrary( loop_sizes() , 1 , 0 ) );
	if ( tag->hasOption( "db_path" ) ) {
		std::string db_path = tag->getOption< string >( "db_path" );
		library_->set_db_path( db_path );
	}
	library_->load_mergeddb();
	library_->mem_foot_print();

	// centroid filter
	string const centroid_filter_name( tag->getOption< string >( "centroid_filter", "true_filter" ) );
	auto find_cenfilter( filters.find( centroid_filter_name ) );
	if ( find_cenfilter == filters.end() ) {
		utility_exit_with_message( "Filter " + centroid_filter_name + " not found in LoopHashDiversifier" );
	}
	cenfilter( find_cenfilter->second );
	ranking_cenfilter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "ranking_cenfilter", centroid_filter_name ), filters ) );
}


//Min RMS
Real
LoopHashDiversifier::min_rms() const{ return min_rms_; }

void
LoopHashDiversifier::min_rms( Real const min_rms ){
	min_rms_ = min_rms;
}

//Max RMS
Real
LoopHashDiversifier::max_rms() const{ return max_rms_; }

void
LoopHashDiversifier::max_rms( Real const max_rms ){
	max_rms_ = max_rms;
}

//Min inter-ss bbrms
Real
LoopHashDiversifier::min_inter_ss_bbrms() const{ return min_inter_ss_bbrms_; }

void
LoopHashDiversifier::min_inter_ss_bbrms( Real const min_inter_ss_bbrms ){
	min_inter_ss_bbrms_ = min_inter_ss_bbrms;
}

//Max inter-ss bbrms
Real
LoopHashDiversifier::max_inter_ss_bbrms() const{ return max_inter_ss_bbrms_; }

void
LoopHashDiversifier::max_inter_ss_bbrms( Real const max_inter_ss_bbrms ){
	max_inter_ss_bbrms_ = max_inter_ss_bbrms;
}


//Min intra-ss bbrms
Real
LoopHashDiversifier::min_intra_ss_bbrms() const{ return min_intra_ss_bbrms_; }

void
LoopHashDiversifier::min_intra_ss_bbrms( Real const min_intra_ss_bbrms ){
	min_intra_ss_bbrms_ = min_intra_ss_bbrms;
}

//Max intra-ss bbrms
Real
LoopHashDiversifier::max_intra_ss_bbrms() const{ return max_intra_ss_bbrms_; }

void
LoopHashDiversifier::max_intra_ss_bbrms( Real const max_intra_ss_bbrms ){
	max_intra_ss_bbrms_ = max_intra_ss_bbrms;
}


//Number of iterations
Size LoopHashDiversifier::num_iterations() const { return num_iterations_; }
void LoopHashDiversifier::num_iterations( Size const num_iterations ){
	num_iterations_ = num_iterations;
}

//Number of trys in each iteration
Size LoopHashDiversifier::num_try_div() const { return num_try_div_; }
void LoopHashDiversifier::num_try_div( Size const num_try_div ){
	num_try_div_ = num_try_div;
}

void
LoopHashDiversifier::cenfilter( protocols::filters::FilterOP cenfilter ) {
	cenfilter_ = cenfilter;
}


utility::vector1< Size >
LoopHashDiversifier::loop_sizes() const{
	return loop_sizes_;
}

void
LoopHashDiversifier::add_loop_size( Size const loop_size ){
	loop_sizes_.push_back( loop_size );
}

std::string LoopHashDiversifier::get_name() const {
	return mover_name();
}

std::string LoopHashDiversifier::mover_name() {
	return "LoopHashDiversifier";
}

void LoopHashDiversifier::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("num_iterations", xsct_non_negative_integer, "Number of loophash runs to execute", "100")
		+ XMLSchemaAttribute::attribute_w_default("num_try_div", xsct_non_negative_integer, "Number of loophash runs to try in each execution", "100")
		+ XMLSchemaAttribute::attribute_w_default("diversify_loop_only", xsct_rosetta_bool, "XRW TO DO", "false")
		+ XMLSchemaAttribute::attribute_w_default("min_inter_ss_bbrms", xsct_real,
		"Min torsion-angle rmsd for residue windows that span two pieces of secondary structure", "0")
		+ XMLSchemaAttribute::attribute_w_default("max_inter_ss_bbrms", xsct_real,
		"Max torsion-angle rmsd for residue windows that span two pieces of secondary structure", "100000")
		+ XMLSchemaAttribute::attribute_w_default("min_intra_ss_bbrms", xsct_real,
		"Min torsion-angle rmsd for residue windows that are a single pieces of secondary structure", "0")
		+ XMLSchemaAttribute::attribute_w_default("max_intra_ss_bbrms", xsct_real,
		"Max torsion-angle rmsd for residue windows that are a single pieces of secondary structure", "100000")
		+ XMLSchemaAttribute::attribute_w_default("min_rms", xsct_real, "Min Anstrom Rmsd cuttofs for loophash generated structures", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("max_rms", xsct_real, "Max Anstrom Rmsd cuttofs for loophash generated structures", "100.0")
		+ XMLSchemaAttribute::attribute_w_default("max_radius", xsct_non_negative_integer, "max loophash radius", "4")
		+ XMLSchemaAttribute::attribute_w_default("max_struct", xsct_non_negative_integer,
		"Max models to create in loophash (number of fragment tried is 200x this)", "10")
		+ XMLSchemaAttribute::attribute_w_default("ideal", xsct_rosetta_bool, "Save space and assume structure is ideal?", "false")
		+ XMLSchemaAttribute::attribute_w_default("filter_by_phipsi", xsct_rosetta_bool, "Filter-out non-ideal phipsi.", "false");

	core::pose::attributes_for_get_resnum_string( attlist, "start_" );
	core::pose::attributes_for_get_resnum_string( attlist, "stop_" );

	attlist
		+ XMLSchemaAttribute::attribute_w_default("window_size", xsct_non_negative_integer, "loophash window size and loophash fragment size", "4")
		+ XMLSchemaAttribute("scorefxn_cen_cst", xs_string, "Score function for constraint centroid mode sampling")
		+ XMLSchemaAttribute("scorefxn_rama_cst", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute("db_path", xs_string, "path to DB -- if not specified then command-line flag is used" )
		+ XMLSchemaAttribute::attribute_w_default("centroid_filter", xs_string, "Filter after centroid stage using this filter.", "true_filter" )
		+ XMLSchemaAttribute("ranking_cenfilter", xs_string, "Prune decoys based on the filter's apply function." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Simple mover that uses loophash to replace randomly chosen fragments in a given pose.", attlist );
}

std::string LoopHashDiversifierCreator::keyname() const {
	return LoopHashDiversifier::mover_name();
}

protocols::moves::MoverOP
LoopHashDiversifierCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopHashDiversifier );
}

void LoopHashDiversifierCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopHashDiversifier::provide_xml_schema( xsd );
}


} //loophash
} //protocols

