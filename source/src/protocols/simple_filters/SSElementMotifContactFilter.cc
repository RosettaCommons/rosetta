// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/SSElementMotifContactFilter
/// @brief filter structures by IntraRepeatContacts
/// @details
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/SSElementMotifContactFilter.hh>
#include <protocols/simple_filters/SSElementMotifContactFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/motif/reference_frames.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/scoring/motif/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoringManager.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
// Parser headers
#include <protocols/filters/Filter.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <set>
#include <map>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer tr("protocols.filters.SSElementMotifContactFilter");

namespace protocols {
namespace simple_filters {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;
using namespace core::scoring::motif;
// @brief default constructor
SSElementMotifContactFilter::SSElementMotifContactFilter():
	Filter( "SSDegree" ),
	filtered_value_( 99 )
{}

// @brief copy constructor
SSElementMotifContactFilter::SSElementMotifContactFilter( SSElementMotifContactFilter const & rval ):
	Super( rval ),
	filtered_value_( rval.filtered_value_ )
{}

// @brief destructor
SSElementMotifContactFilter::~SSElementMotifContactFilter() = default;

// @brief set filtered value
void SSElementMotifContactFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
SSElementMotifContactFilter::Real
SSElementMotifContactFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
}

/// @brief
void
SSElementMotifContactFilter::report( std::ostream & out, Pose const & pose ) const
{
	if ( report_avg_ ) {
		out << "Avg ss contact" <<  compute( pose ) << std::endl;
	} else {
		out << "Worst ss contact" <<  compute( pose ) << std::endl;
	}
}

/// @brief get ss_elements
protocols::loops::Loops SSElementMotifContactFilter::get_ss_elements(const Pose & pose) const{
	using protocols::loops::Loop;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char lastSecStruct = dssp.get_dssp_secstruct( 1 );
	Size startSS = 0;
	Size endSS = 0;
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
		if ( option[OptionKeys::score::motif_ignore_symmmetry]() ) {
			nres1= core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
	}
	protocols::loops::Loops ss_elements;
	if ( dssp.get_dssp_secstruct(1)=='H' || dssp.get_dssp_secstruct(1)=='E' ) {
		startSS=  1;
	}
	for ( core::Size ii = 2; ii <= nres1; ++ii ) {
		if ( (dssp.get_dssp_secstruct(ii) == 'H' || dssp.get_dssp_secstruct(ii)) && lastSecStruct == 'L' ) {
			startSS = ii;
		}
		if ( dssp.get_dssp_secstruct(ii)!=lastSecStruct && (lastSecStruct ==  'H' || lastSecStruct ==  'E') ) {
			endSS = ii-1;
			if ( endSS-startSS >= 2 ) {
				ss_elements.add_loop(Loop(startSS,endSS));
			}
		}
		lastSecStruct = dssp.get_dssp_secstruct(ii);
	}
	return(ss_elements);
}

Size SSElementMotifContactFilter::which_ssElement(Size res,protocols::loops::Loops ssElements) const{
	for ( Size ii=1; ii<=ssElements.size(); ++ii ) {
		if ( res>= ssElements[ii].start() && res<=ssElements[ii].stop() ) {
			return(ii);
		}
	}
	return(0);
}

Size SSElementMotifContactFilter::get_ssElements_in_contact_w_threshold(std::multiset<Size> ssElements_in_contact) const {
	std::map<Size,Size> element_ct;
	std::multiset<Size>::iterator itr;
	std::map<Size,Size>::iterator map_itr;
	for ( itr = ssElements_in_contact.begin(); itr!=ssElements_in_contact.end(); ++itr ) {
		Size counts = ssElements_in_contact.count(*itr);
		element_ct.insert(std::pair<Size,Size>(*itr,counts));
	}
	Size above_threshold_ct =0;
	for ( map_itr=element_ct.begin(); map_itr!=element_ct.end(); ++map_itr ) {
		if ( map_itr->second>=contacts_between_ssElement_threshold_ ) {
			above_threshold_ct++;
		}
	}
	return(above_threshold_ct);
}


Size SSElementMotifContactFilter::get_SSelements_in_contact(Size element,protocols::loops::Loops ssElements, const Pose & pose) const{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::kinematics::FoldTree;
	std::multiset<Size> ssElements_in_contact;
	core::scoring::dssp::Dssp dssp( pose );
	const FoldTree& tree = pose.fold_tree();
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
		if ( option[OptionKeys::score::motif_ignore_symmmetry]() ) {
			nres1= core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
	}
	for ( size_t ir = ssElements[element].start(); ir <= ssElements[element].stop(); ++ir ) {
		Xform const ibb_stub = core::pose::motif::get_backbone_reference_frame(pose,ir);
		char ss1 = dssp.get_dssp_secstruct( ir );
		char aa1 = pose.residue(ir).name1();
		for ( size_t jr = 1; jr <= nres1; ++jr ) {
			if ( !tree.is_jump_point(ir) && !tree.is_jump_point(jr) ) {
				Real dist = pose.residue(ir).xyz("CA").distance(pose.residue(jr).xyz("CA"));
				if ( dist < 12 ) {
					char ss2 = dssp.get_dssp_secstruct( jr );
					char aa2 = pose.residue(jr).name1();
					Xform const jbb_stub = core::pose::motif::get_backbone_reference_frame(pose,jr);
					Xform const Xbb = ibb_stub.inverse() * jbb_stub;
					core::scoring::motif::XformScoreCOP xs_bb_fxn1(mman_->get_xform_score_BB_BB(ss1,ss2,aa1,aa2));
					core::scoring::motif::XformScoreCOP xs_bb_fxn2(mman_->get_xform_score_BB_BB(ss2,ss1,aa2,aa1));
					Real tmpScore = 0;
					if ( xs_bb_fxn1 != nullptr ) {
						tmpScore += xs_bb_fxn1->score_of_bin(Xbb);
						tmpScore += xs_bb_fxn2->score_of_bin(Xbb.inverse());
					}
					if ( tmpScore>0 ) {
						Size resHit = which_ssElement(jr,ssElements);
						if ( resHit != 0 && resHit !=element ) { //must be in SS element and not the current element
							ssElements_in_contact.insert(resHit);
						}
					}
				}
			}
		}
	}
	return(get_ssElements_in_contact_w_threshold(ssElements_in_contact));
}


/// @brief
Real
SSElementMotifContactFilter::compute( const Pose & pose ) const
{
	using protocols::loops::Loop;
	protocols::loops::Loops ss_elements = get_ss_elements(pose);
	Size startElement = 1;
	Size endElement = ss_elements.size();
	if ( ignore_terminal_SS_>0 ) {
		startElement+=ignore_terminal_SS_;
		endElement-=ignore_terminal_SS_;
	}
	Real score = 0;
	if ( only_n_term_ ) {
		return(get_SSelements_in_contact(1,ss_elements,pose));
	}
	if ( only_c_term_ ) {
		return(get_SSelements_in_contact(ss_elements.size(),ss_elements,pose));
	}
	if ( report_avg_ ) {
		Real tmpScore = 0;
		Size ct = 0;
		for ( Size ii=startElement; ii<=endElement; ++ii ) {
			tmpScore+=get_SSelements_in_contact(ii,ss_elements,pose);
			ct+=1;
		}
		if ( ct!=0 ) {
			score = tmpScore/Real(ct);
		}
	} else {
		Real tmpScore = 9999;
		for ( Size ii=startElement; ii<=endElement; ++ii ) {
			Real tmp=get_SSelements_in_contact(ii,ss_elements,pose);
			if ( tmp<tmpScore ) {
				tmpScore=tmp;
			}
		}
		score= tmpScore;
	}
	return(score);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SSElementMotifContactFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	tr << "value" << value << "filtered_value_" << threshold_ << std::endl;
	if ( value >= threshold_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << threshold_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
SSElementMotifContactFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set threshold
	threshold_ = tag->getOption<Real>("threshold",0);
	contacts_between_ssElement_threshold_ = tag->getOption<Real>("contacts_between_ssElement_threshold",3);
	report_avg_ = tag->getOption<bool>( "report_avg", true );
	ignore_terminal_SS_ = tag->getOption<Size>("ignore_terminal_ss",0);
	only_n_term_ = tag->getOption<bool>("only_n_term",false);
	only_c_term_ = tag->getOption<bool>("only_c_term",false);
	mman_ = core::scoring::motif::MotifHashManager::get_instance();

}

// XRW TEMP filters::FilterOP
// XRW TEMP SSElementMotifContactFilterCreator::create_filter() const { return protocols::filters::FilterOP(new SSElementMotifContactFilter); }

// XRW TEMP std::string
// XRW TEMP SSElementMotifContactFilterCreator::keyname() const { return "SSDegree"; }

std::string SSElementMotifContactFilter::name() const {
	return class_name();
}

std::string SSElementMotifContactFilter::class_name() {
	return "SSDegree";
}

void SSElementMotifContactFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("contacts_between_ssElement_threshold", xsct_real, "XRW TO DO", "3")
		+ XMLSchemaAttribute::attribute_w_default("report_avg", xsct_rosetta_bool, "XRW TO DO", "true")
		+ XMLSchemaAttribute::attribute_w_default("ignore_terminal_ss", xsct_non_negative_integer, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("only_n_term", xsct_rosetta_bool, "XRW TO DO", "false")
		+ XMLSchemaAttribute::attribute_w_default("only_c_term", xsct_rosetta_bool, "XRW TO DO", "false");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SSElementMotifContactFilterCreator::keyname() const {
	return SSElementMotifContactFilter::class_name();
}

protocols::filters::FilterOP
SSElementMotifContactFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SSElementMotifContactFilter );
}

void SSElementMotifContactFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SSElementMotifContactFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
