// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AlignEndsMover.cc
/// @brief

// Unit headers
#include <protocols/splice/AlignEndsMover.hh>
#include <protocols/splice/AlignEndsMoverCreator.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.splice.AlignEndsMover" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <string>
#include <core/import_pose/import_pose.hh>
#include <protocols/splice/FindEndpointsOperation.hh>
#include <protocols/toolbox/superimpose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace splice {

using namespace::protocols;

// XRW TEMP std::string
// XRW TEMP AlignEndsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return AlignEndsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AlignEndsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AlignEndsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AlignEndsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "AlignEnds";
// XRW TEMP }

AlignEndsMover::AlignEndsMover(): moves::Mover("AlignEnds"),
	distance_threshold_( 18.0 ),
	strand_length_( 3 ),
	neighbors_( 6 ),
	N_terminal_count_( 3 ),
	max_strands_( 10 ),
	odd_( true ),
	even_( true ),
	template_pose_( /* NULL */ ),
	stagger_( 0 ),
	parallel_( true ),
	chain_( 1 ),
	sequence_separation_( 15 )
{
	residues_to_align_on_pose_.clear();
	residues_to_align_on_template_.clear();
}

AlignEndsMover::~AlignEndsMover() = default;

void
AlignEndsMover::template_pose( core::pose::PoseOP p ){ template_pose_ = p; }

core::pose::PoseOP
AlignEndsMover::template_pose() const{ return template_pose_; }

utility::vector1< core::Size >
AlignEndsMover::reference_positions( core::pose::Pose const & pose ) const{
	using namespace core::scoring::dssp;

	Dssp dssp( pose );
	dssp.dssp_reduced();

	utility::vector1< core::Size > strand_positions;
	strand_positions.clear();
	bool odd_strand( true ); // if antiparallel we need to intermittently connect odd and even strands...
	core::Size const resi_end( chain() == 0 ? pose.size() : pose.conformation().chain_end( chain() ) );
	core::Size const resi_start( chain() == 0 ? 1 : pose.conformation().chain_begin( chain() ) );
	TR<<"Chain start: "<<resi_start<<" resi_end: "<<resi_end<<std::endl;
	// note the fast-forwarding of resi in the inner-loop below! should be fine, but if there are bugs, this is a good place to dig
	for ( core::Size resi = resi_start; resi <= resi_end - strand_length() + 1; ++resi ) {  // at least n-strand positions in a row
		if ( dssp.get_dssp_secstruct( resi ) == 'E' ) {
			core::Size strand_count = 1;
			//find end of beta strand
			for ( ; strand_count < strand_length() && dssp.get_dssp_secstruct( resi + strand_count ) == 'E'; ++strand_count ) {}
			if ( strand_count == strand_length() ) { // we have n-length beta strand
				if ( !parallel() && !odd_strand ) {
					// move to the end of the strand (even strands will be represented by the C-termini)
					for ( ; resi <= resi_end && dssp.get_dssp_secstruct( resi + strand_length() ) == 'E'; ++resi ) {}
				}  // fi !parallel() && !odd()
				core::Size position_count = 0; // how many strand positions were added to the strand_positions array
				for ( ; resi <= resi_end && dssp.get_dssp_secstruct( resi ) == 'E'; ++resi ) { // this causes 'fast-forwarding' of resi, but it's okay, since we don't want to double count the strands
					// we'll add only strand_length() positions to the array, but will forward the resi index to the end of the strand.
					if ( position_count < strand_length() ) {
						strand_positions.push_back( resi );
						++position_count;
					}
				}  // for resi inner loop
				odd_strand = !odd_strand;  // flip odd->even->odd...
			}  // fi strand_count
		}  // fi dssp.get_dssp_secstruct(resi)
	}  // for resi outer loop
	TR<<"DEBUG: strand positions: ";
	for ( core::Size const sp : strand_positions ) {
		TR<<sp<<'+';
	}
	TR<<std::endl;
	utility::vector1< core::Size > strand_ntermini;
	strand_ntermini.clear();
	odd_strand = true; // antiparallel strands have neighboring c/n terminal positions whereas parallel strands have neighboring n terminal positions.
	for ( utility::vector1< core::Size >::const_iterator resi = strand_positions.begin(); resi != strand_positions.end(); resi += strand_length() ) { /// find strand ntermini
		if ( !odd_strand && !parallel() ) {
			strand_ntermini.push_back( *(resi + strand_length() - 1 ) );
		} else {
			strand_ntermini.push_back( *resi );
		}
		odd_strand = !odd_strand;/// flip odd->even->odd...
	}
	TR<<"DEBUG: nterminal positions: ";
	for ( core::Size const resi : strand_ntermini ) {
		TR<<resi<<'+';
	}
	TR<<std::endl;
	utility::vector1< core::Size > ntermini_w_neighbors;
	ntermini_w_neighbors.clear();
	/// In the following we ensure that each of the ntermini has neighbors within the ntermini array. We're looking for tightly packed barrel structures.
	for ( core::Size const res : strand_ntermini ) { /// strand_ntermini with enough neighbors
		core::Size const resi_neighbors( neighbors_in_vector( pose, res, strand_ntermini, distance_threshold(), dssp, sequence_separation() ));
		if ( resi_neighbors >= neighbors() ) {
			ntermini_w_neighbors.push_back( res );
			TR<<res<<" has "<<resi_neighbors<<" neighbors"<<std::endl;
		}//fi resi_neighbors
	}//foreach res

	utility::vector1< core::Size > ntermini_w_close_neighbors; // ntermini with 2 close neighbors
	ntermini_w_close_neighbors.clear();
	if ( parallel() ) {
		for ( core::Size const res : ntermini_w_neighbors ) {
			core::Size const close_neighbors( neighbors_in_vector( pose, res, strand_ntermini, 10.0, dssp, sequence_separation() ));
			if ( close_neighbors >= 2 ) {
				ntermini_w_close_neighbors.push_back( res );
			} else {
				core::Size const close_neighbors_one_down( neighbors_in_vector( pose, res+1, strand_ntermini, 10.0, dssp, sequence_separation() ));
				if ( close_neighbors_one_down >= 2 ) {
					ntermini_w_close_neighbors.push_back( res + 1);
				}
			}
			TR<<res<<" has "<<close_neighbors<<" close neighbors"<<std::endl;
		}
	} else { //fi parallel /// antiparallel beta barrels are less well packed; ignoring close neighbors condition
		ntermini_w_close_neighbors = ntermini_w_neighbors;
	}
	TR<<"Found "<<ntermini_w_close_neighbors.size()<<" strands that fit all of the criteria"<<std::endl;
	if ( ntermini_w_close_neighbors.size() > max_strands() ) {
		TR<<"max_strands = "<<max_strands()<<std::endl;
		ntermini_w_close_neighbors.resize( max_strands() );
		TR<<"Only considering the first (n-terminal) "<<ntermini_w_close_neighbors.size()<<" strands. The c-terminal strands are ignored"<<std::endl;
	}

	utility::vector1< core::Size > positions_for_alignment;
	positions_for_alignment.clear();
	odd_strand = true;
	for ( core::Size const res : ntermini_w_close_neighbors ) {
		if ( ( odd_strand && odd() ) || ( !odd_strand && even() ) ) {
			if ( !odd_strand && !parallel() ) {
				for ( core::Size count = 0; count < N_terminal_count(); ++count ) {
					positions_for_alignment.push_back( res - count ); // count residues backwards if antiparallel and even
				}
			} else {
				for ( core::Size count = 0; count < N_terminal_count(); ++count ) {
					positions_for_alignment.push_back( res + count );
				}
			}
		}//fi
		odd_strand = !odd_strand;
	}//foreach
	return positions_for_alignment;
}

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	for ( core::Size const pos : positions ) {
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

void
AlignEndsMover::apply( Pose & pose ){
	using namespace protocols::toolbox;
	using namespace core::scoring::dssp;

	runtime_assert( residues_to_align_on_pose_.size() == residues_to_align_on_template_.size() );
	utility::vector1< core::Size > template_positions, pose_positions;
	template_positions.clear(); pose_positions.clear();
	if ( residues_to_align_on_pose_.size() == 0 ) {
		template_positions = reference_positions( *template_pose() );
		pose_positions = reference_positions( pose );
	} else {
		template_positions = residues_to_align_on_template_;
		pose_positions = residues_to_align_on_pose_;
	}

	TR<<"Aligning positions on template: ";
	for ( core::Size const p : template_positions ) {
		TR<<p<<'+';
	}
	TR<<std::endl;
	TR<<"To pose positions: ";
	for ( core::Size const p : pose_positions ) {
		TR<<p<<'+';
	}
	TR<<std::endl;


	utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coords( pose, pose_positions ) ), ref_coords( Ca_coords( *template_pose(), template_positions ) );

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	runtime_assert( init_coords.size() >= N_terminal_count() * stagger() );
	std::rotate( init_coords.begin(), init_coords.begin() + ( N_terminal_count() * stagger()), init_coords.end());

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( pose, rotation, to_init_center, to_fit_center );
}

// XRW TEMP std::string
// XRW TEMP AlignEndsMover::get_name() const {
// XRW TEMP  return AlignEndsMover::mover_name();
// XRW TEMP }

moves::MoverOP
AlignEndsMover::clone() const
{
	return moves::MoverOP( new AlignEndsMover( *this ) );
}

moves::MoverOP
AlignEndsMover::fresh_instance() const
{
	return moves::MoverOP( new AlignEndsMover );
}

void
AlignEndsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	parallel( tag->getOption< bool >( "parallel", true ) );
	distance_threshold( tag->getOption< core::Real >( "distance_threshold", 18.0 ));
	neighbors( tag->getOption< core::Size >( "neighbors", 6 ) );
	N_terminal_count( tag->getOption< core::Size >( "N_terminal_count", 3 ) );
	odd( tag->getOption< bool >( "odd", true ) );
	even( tag->getOption< bool >( "even", true ) );
	std::string const template_fname( tag->getOption< std::string >( "template_pose" ) );
	template_pose( core::import_pose::pose_from_file( template_fname ) );
	stagger( tag->getOption< core::Size >( "stagger", 0 ) );
	strand_length( tag->getOption< core::Size >( "strand_length", 3 ) );
	max_strands( tag->getOption< core::Size >( "max_strands", 10 ) );
	chain( tag->getOption< core::Size >( "chain", 1 ) );
	sequence_separation( tag->getOption< core::Size >( "sequence_separation", 15 ) );
	if ( tag->hasOption( "residues_to_align_on_pose" ) && tag->hasOption( "residues_to_align_on_template" ) ) {
		TR<<"Will respect user defined residues for alignment, ignoring beta-strand calculations"<<std::endl;
		residues_to_align_on_pose_ = utility::string_split(tag->getOption<std::string>("residues_to_align_on_pose"),'+',core::Size());
		residues_to_align_on_template_ = utility::string_split(tag->getOption<std::string>("residues_to_align_on_template"),'+',core::Size());
		runtime_assert( residues_to_align_on_pose_.size() == residues_to_align_on_template_.size() );
		TR<<"Aligning pose residues: ";
		for ( core::Size const s : residues_to_align_on_pose_ ) {
			TR<<s<<',';
		}
		TR<<"\nwith template residues";
		for ( core::Size const s : residues_to_align_on_template_ ) {
			TR<<s<<',';
		}
		TR<<std::endl;
	}//fi tag->hasOption

	TR<<"paralell: "<<parallel()<<" distance_threshold: "<<distance_threshold()<<" max_strands: "<<max_strands()<<" strand_length: "<<strand_length()<<" neighbors: "<<neighbors()<<" N_terminal_count: "<<N_terminal_count()<<" odd: "<<odd()<<" even: "<<even()<<" template_pose: "<<template_fname<<" stagger: "<<stagger()<<" chain: "<<chain()<<" sequence_separation: "<<sequence_separation()<<std::endl;
}

std::string AlignEndsMover::get_name() const {
	return mover_name();
}

std::string AlignEndsMover::mover_name() {
	return "AlignEnds";
}

void AlignEndsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// AMW TODO
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "parallel", xsct_rosetta_bool, "Align ends in a parallel orientation", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_threshold", xsct_real, "XRW TO DO", "18.0" )
		+ XMLSchemaAttribute::attribute_w_default( "neighbors", xsct_non_negative_integer, "XRW TO DO", "6" )
		+ XMLSchemaAttribute::attribute_w_default( "N_terminal_count", xsct_non_negative_integer, "XRW TO DO", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "odd", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "even", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "template_pose", xs_string, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "stagger", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "strand_length", xsct_non_negative_integer, "XRW TO DO", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "max_strands", xsct_non_negative_integer, "XRW TO DO", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "chains", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "sequence_separation", xsct_non_negative_integer, "XRW TO DO", "15" )
		+ XMLSchemaAttribute( "residues_to_align_on_pose", xs_string, "+ separated list of non-negative integers" )
		+ XMLSchemaAttribute( "residues_to_align_on_template", xs_string, "+ separated list of non-negative integers" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string AlignEndsMoverCreator::keyname() const {
	return AlignEndsMover::mover_name();
}

protocols::moves::MoverOP
AlignEndsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AlignEndsMover );
}

void AlignEndsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AlignEndsMover::provide_xml_schema( xsd );
}

} // simple_moves
} // protocols

