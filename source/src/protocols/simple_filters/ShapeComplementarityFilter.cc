// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/ShapeComplementarityFilter.cc
/// @brief  Filter structures by shape complementarity and/or interface area
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

// Unit Headers
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>
#include <protocols/simple_filters/ShapeComplementarityFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


//// C++ headers
static thread_local basic::Tracer tr( "protocols.simple_filters.ShapeComplementarityFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
ShapeComplementarityFilter::ShapeComplementarityFilter():
	Filter( "ShapeComplementarity" ),
	filtered_sc_( 0.50 ),
	filtered_area_( 250 ),
	jump_id_( 1 ),
	quick_( false ),
	verbose_( false ),
	residues1_( ),
	residues2_( ),
	write_int_area_(false),
	multicomp_(false),
	sym_dof_name_("")
{}


// @brief constructor with arguments
ShapeComplementarityFilter::ShapeComplementarityFilter( Real const & filtered_sc, Real const & filtered_area,
	Size const & jump_id, Size const & quick, Size const & verbose):
	Filter( "ShapeComplementarity" ),
	filtered_sc_( filtered_sc ),
	filtered_area_( filtered_area ),
	jump_id_( jump_id ),
	quick_( quick ),
	verbose_( verbose ),
	residues1_( ),
	residues2_( ),
	sym_dof_name_("")
{}

// @brief copy constructor
ShapeComplementarityFilter::ShapeComplementarityFilter( ShapeComplementarityFilter const & rval ):
	Super( rval ),
	filtered_sc_( rval.filtered_sc_ ),
	filtered_area_( rval.filtered_area_ ),
	jump_id_( rval.jump_id_ ),
	quick_( rval.quick_ ),
	verbose_( rval.verbose_ ),
	residues1_( rval.residues1_ ),
	residues2_( rval.residues2_ ),
	write_int_area_( rval.write_int_area_ ),
	multicomp_( rval.multicomp_ ),
	sym_dof_name_( rval.sym_dof_name_ )
{}

void ShapeComplementarityFilter::filtered_sc( Real const & filtered_sc ) { filtered_sc_ = filtered_sc; }
void ShapeComplementarityFilter::filtered_area( Real const & filtered_area ) { filtered_area_ = filtered_area; }
void ShapeComplementarityFilter::jump_id( Size const & jump_id ) { jump_id_ = jump_id; }
void ShapeComplementarityFilter::quick( Size const & quick ) { quick_ = quick; }
void ShapeComplementarityFilter::verbose( Size const & verbose ) { verbose_ = verbose; }
void ShapeComplementarityFilter::residues1( utility::vector1< core::Size > const & residues ) { residues1_ = residues; }
void ShapeComplementarityFilter::residues2( utility::vector1< core::Size > const & residues ) { residues2_ = residues; }
void ShapeComplementarityFilter::sym_dof_name( std::string const & sym_dof_name ) { sym_dof_name_ = sym_dof_name; }
std::string ShapeComplementarityFilter::sym_dof_name() const { return sym_dof_name_; }
void ShapeComplementarityFilter::write_int_area( bool write_int_area ) { write_int_area_ = write_int_area; }
bool ShapeComplementarityFilter::write_int_area() const { return write_int_area_; }
void ShapeComplementarityFilter::multicomp( bool multicomp ) { multicomp_ = multicomp; }
bool ShapeComplementarityFilter::multicomp() const { return multicomp_; }

/// @brief
core::Size ShapeComplementarityFilter::compute( Pose const & pose ) const
{
	if(scc_.GetResults().valid)
		return 1;

	if(!scc_.Init())
		return 0;
	if(quick_)
		scc_.settings.density = 5.0;
	scc_.Reset();

	bool symm = core::pose::symmetry::is_symmetric( pose );
	Real nsubs_scalefactor = 1.0;

	if (!residues1_.empty() && !residues2_.empty()) {
		for(utility::vector1<Size>::const_iterator r = residues1_.begin();
			r != residues1_.end(); ++r)
				scc_.AddResidue(0, pose.residue(*r));

		for(utility::vector1<Size>::const_iterator r = residues2_.begin();
			r != residues2_.end(); ++r)
				scc_.AddResidue(1, pose.residue(*r));

		if(!scc_.Calc())
			return 0;

	} else if (!symm) {

		if(!scc_.Calc( pose, jump_id_ ))
			return 0;

	} else {
		if ( multicomp_ ) {
			// MULTI COMPONENT SYMM
			runtime_assert( sym_dof_name() != "" );
			utility::vector1<std::string> sym_dof_name_list = utility::string_split( sym_dof_name() , ',' );

			Size sym_dof_index = 1;

			if( sym_dof_name_list.size() > 1) {
				Size intracontact_count = 0;
 				for (core::Size i=1; i<=sym_dof_name_list.size(); i++) {
					if( core::pose::symmetry::intracomponent_contact(pose, sym_dof_name_list[i], 12.0)) {
						intracontact_count++;
						sym_dof_index = i;
					}
				}
				if( intracontact_count > 1) {
					tr.Warning << "Intracontacts detected between multiple components.  Calculating sc based off of the last component with intracontacts.  Separate sc calculations are recommended for each component with intracontacts." << std::endl;
				}
			}

			utility::vector1<Size> full_intracomp_resis = core::pose::symmetry::get_full_intracomponent_resis(pose, sym_dof_name_list[sym_dof_index]);
			utility::vector1<Size> full_intracomp_neighbor_resis = core::pose::symmetry::get_full_intracomponent_neighbor_resis(pose, sym_dof_name_list[sym_dof_index], 12.0 );

			//core::pose::Pose full_intracomp_subpose = core::pose::symmetry::get_full_intracomponent_subpose(pose, sym_dof_name_list[sym_dof_index]);
			//core::pose::Pose full_intracomp_neighbor_subpose = core::pose::symmetry::get_full_intracomponent_neighbor_subpose(pose, sym_dof_name_list[sym_dof_index], 12.0);
			//full_intracomp_subpose.dump_pdb("full_intracomp_subpose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");
			//full_intracomp_neighbor_subpose.dump_pdb("full_intracomp_neighbor_subpose_" + protocols::jd2::JobDistributor::get_instance()->current_output_name() + ".pdb");

 			for (core::Size i=1; i<=full_intracomp_resis.size(); i++) {
				scc_.AddResidue(0,pose.residue(full_intracomp_resis[i]));
			}
 			for (core::Size i=1; i<=full_intracomp_neighbor_resis.size(); i++) {
				scc_.AddResidue(1,pose.residue(full_intracomp_neighbor_resis[i]));
			}

			//tr << "Using jump_id " << sym_aware_jump_id << " to partition pose" << std::endl;
			//if(!scc_.Calc( pose, sym_aware_jump_id ))
			if(!scc_.Calc())
				return 0;

			utility::vector1<Size> subs = core::pose::symmetry::get_jump_name_to_subunits( pose, sym_dof_name_list[sym_dof_index] );
			nsubs_scalefactor = (Real) subs.size() ;
		} else {
			// SINGLE COMPONENT SYMM
			ObjexxFCL::FArray1D_bool is_upstream ( pose.total_residue(), false );
			utility::vector1<Size> sym_aware_jump_ids;

			if ( sym_dof_name() != "" ) {
				sym_aware_jump_ids.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() ) );
			} else {
				// all slidable jumps
				Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
				for (Size j = 1; j <= nslidedofs; j++)
					sym_aware_jump_ids.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
			}

			// partition & fill residueX_ vectors
			core::pose::symmetry::partition_by_symm_jumps( sym_aware_jump_ids, pose.fold_tree(), core::pose::symmetry::symmetry_info(pose), is_upstream );
			Size ndownstream=0;
			for (core::Size i=1; i<=pose.total_residue(); ++i) {
				if (pose.residue(i).aa() == core::chemical::aa_vrt) continue;
				scc_.AddResidue(is_upstream(i)?1:0, pose.residue(i));
				if (!is_upstream(i)) ndownstream++;
			}
			// scalefactor
			nsubs_scalefactor = (Real)( ndownstream / core::pose::symmetry::symmetry_info(pose)->get_nres_subunit() );

			if(!scc_.Calc())
				return 0;
		}
	}

	core::scoring::sc::RESULTS const &r = scc_.GetResults();
	if(verbose_) {

		// Verbose view
		tr << "==================================================" << std::endl;
		tr << std::endl;
		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << ":" << std::endl;
			else
				tr << "Total/Average for both molecules:" << std::endl;

			tr << "          Total Atoms: " << r.surface[i].nAtoms << std::endl;
			tr << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << std::endl;
			tr << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << std::endl;
			tr << "           Total Dots: " << r.surface[i].nAllDots << std::endl;
			tr << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << std::endl;
			tr << "         Trimmed Area: " << r.surface[i].trimmedArea << " (avg) " << std::endl;
			tr << std::endl;
    }
		tr << std::endl;

		for(int i = 0; i <= 2; i++) {
			if(i < 2)
				tr << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << std::endl;
			else
				tr << "Average for both molecules:" << std::endl;
			tr << "      Mean Separation: " << r.surface[i].d_mean << std::endl;
			tr << "    Median Separation: " << r.surface[i].d_median << std::endl;
			tr << "    Mean Shape Compl.: " << r.surface[i].s_mean << std::endl;
			tr << "  Median Shape Compl.: " << r.surface[i].s_median << std::endl;
			tr << std::endl;
		}

	}

	tr << "Shape complementarity: " << r.sc << std::endl;
	tr << "Interface area: " << r.area << std::endl;
	if ( nsubs_scalefactor != 1) {
		tr << "Area per monomer: " << ( (core::Real) r.area / nsubs_scalefactor ) << std::endl ;
	}
	tr << "Interface seperation: " << r.distance << std::endl;

	return 1;
}

/// @brief
core::Real ShapeComplementarityFilter::report_sm( Pose const & pose ) const
{
	scc_.Reset(); // Unfortunately, this line had to be added. While reducing
								// efficiency in normal use cases by forcing recalculation of
								// presumably the same value, it is necessary for greedy
								// optimization using the GreedyOptMutationMover, which calls
								// the report_sm() function of filters directly. -Neil King
	if(compute( pose )) {
		if ( write_int_area_ ) {
			protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
			std::string column_header = this->get_user_defined_name() + "_int_area";
			core::Real int_area = scc_.GetResults().area ;

			// symmetric scalefactor
			if (core::pose::symmetry::is_symmetric( pose )) {
				if ( multicomp_ ) {
					utility::vector1<Size> subs = core::pose::symmetry::get_jump_name_to_subunits( pose, sym_dof_name() );
					int_area /= (Real) subs.size() ;
				} else {
					ObjexxFCL::FArray1D_bool is_upstream ( pose.total_residue(), false );
					utility::vector1<Size> sym_aware_jump_ids;
					if ( sym_dof_name() != "" ) {
						sym_aware_jump_ids.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name() ) );
					} else {
						Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
						for (Size j = 1; j <= nslidedofs; j++) sym_aware_jump_ids.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
					}
					core::pose::symmetry::partition_by_symm_jumps( sym_aware_jump_ids, pose.fold_tree(), core::pose::symmetry::symmetry_info(pose), is_upstream );
					Size ndownstream=0;
					for (Size i=1; i<=pose.total_residue(); ++i) {
									if (pose.residue(i).aa() == core::chemical::aa_vrt) continue;
									if (!is_upstream(i)) ndownstream++;
					}
					int_area /= (Real)( ndownstream / core::pose::symmetry::symmetry_info(pose)->get_nres_subunit() );
 				}
			}
			job->add_string_real_pair(column_header, int_area );
		}
		return scc_.GetResults().sc;
	}
	tr.Error << "Issue computing shape complementarity value - returning -1 instead." << std::endl;
	if( write_int_area_ ) {
		// Need to add placeholder so that all structures in a run have the same number of scorefile headers
		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		std::string column_header = this->get_user_defined_name() + "_int_area";
		job->add_string_real_pair(column_header, -1 );
	}
	return -1;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose has high enough shape
// complementarity.
bool ShapeComplementarityFilter::apply( Pose const & pose ) const
{
	scc_.Reset();

	if(!compute( pose )) {
		tr.Error << "Issue computing shape complementarity value - failing filter." << std::endl;
		return false;
	}

	Real sc = scc_.GetResults().sc;
	Real area = scc_.GetResults().area;

	if( sc < filtered_sc_ ) {
		tr << "Filter failed current < threshold sc: " << sc << " < " << filtered_sc_ << std::endl;
		return false;
	}

	if( area < filtered_area_ ) {
		tr << "Filter failed current < threshold interface area: " << area << " < " << filtered_area_ << std::endl;
		return false;
	}

	tr << "Successfully filtered: " << sc << std::endl;
	return true;
} // apply_filter

/// @brief parse xml
void
ShapeComplementarityFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap &,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & pose )
{
	filtered_sc_ = tag->getOption<Real>( "min_sc", 0.50 );
	filtered_area_ = tag->getOption<Real>( "min_interface", 0 );
	verbose_ = tag->getOption<Size>( "verbose", false );
	quick_ = tag->getOption<Size>( "quick", false );
	jump_id_ = tag->getOption<Size>( "jump", 1 );
	write_int_area_ = tag->getOption<bool>( "write_int_area", false );
	sym_dof_name(tag->getOption<std::string>( "sym_dof_name", "" ));
	multicomp( tag->getOption< bool >("multicomp", false) );

	if(tag->hasOption("residues1")) {
		residues1_ = core::pose::get_resnum_list(tag, "residues1", pose);
		if(residues1_.empty())
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues1") << ". Using default." << std::endl;
	}
	if(tag->hasOption("residues2")) {
		residues2_ = core::pose::get_resnum_list(tag, "residues2", pose);
		if(residues2_.empty())
			tr.Warning << "Failed to parse residue range: " << tag->getOption<std::string> ("residues2") << ". Using default." << std::endl;
	}

	tr.Info << "Structures with shape complementarity < " << filtered_sc_ << ", interface area < " <<
		filtered_area_ << " A^2 will be filtered." << std::endl;

	if(quick_)
		tr.Info << "Calculating shape complementarity in quick mode with less accuracy." << std::endl;
	if(!residues1_.empty() && !residues2_.empty()) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		for(utility::vector1<Size>::const_iterator r = residues1_.begin(); r != residues1_.end(); ++r)
	                tr.Info << (r == residues1_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		for(utility::vector1<Size>::const_iterator r = residues2_.begin(); r != residues2_.end(); ++r)
	                tr.Info << (r == residues2_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
	} else {
		if(!residues1_.empty() || !residues2_.empty())
			tr.Warning << "Ignoring residue range selection since residues" << (residues1_.empty() ? 1 : 2) << " is empty." << std::endl;
		if(jump_id_ != 1)
			tr.Info << "Using Jump ID " << jump_id_ << " to define surfaces." << std::endl;
	}
}

void ShapeComplementarityFilter::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & /*score_fxns*/,
		utility::lua::LuaObject const & /*tasks*/ ) {
	filtered_sc_ = def["min_sc"] ? def["min_sc"].to<Real>() : 0.50;
	filtered_area_ = def["min_interface"] ? def["min_interface"].to<Real>() : 0;
	verbose_ = def["verbose"] ? def["verbose"].to<Size>() : 0;
	quick_ = def["quick"] ? def["quick"].to<Size>() : 0;
	jump_id_ = def["jump"] ? def["jump"].to<Size>() : 1;
	write_int_area_ = def["write_int_area"] ? def["write_int_area"].to<bool>() : false;

	if(def["residues1"] ) {
		for (utility::lua::LuaIterator i=def["residues1"].begin(), end; i != end; ++i) {
			residues1_.push_back( (*i).to<core::Size>() );
		}
		if(residues1_.empty())
			tr.Warning << "Failed to parse residue1 range, using default." << std::endl;
	}
	if(def["residues2"] ) {
		for (utility::lua::LuaIterator i=def["residues2"].begin(), end; i != end; ++i) {
			residues2_.push_back( (*i).to<core::Size>() );
		}
		if(residues2_.empty())
			tr.Warning << "Failed to parse residue2 range, using default." << std::endl;
	}

	tr.Info << "Structures with shape complementarity < " << filtered_sc_ << ", interface area < " <<
		filtered_area_ << " A^2 will be filtered." << std::endl;

	if(quick_)
		tr.Info << "Calculating shape complementarity in quick mode with less accuracy." << std::endl;
	if(!residues1_.empty() && !residues2_.empty()) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		for(utility::vector1<Size>::const_iterator r = residues1_.begin(); r != residues1_.end(); ++r)
	                tr.Info << (r == residues1_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		for(utility::vector1<Size>::const_iterator r = residues2_.begin(); r != residues2_.end(); ++r)
	                tr.Info << (r == residues2_.begin() ? "" : ", ") << *r;
		tr.Info << std::endl;
	} else {
		if(!residues1_.empty() || !residues2_.empty())
			tr.Warning << "Ignoring residue range selection since residues" << (residues1_.empty() ? 1 : 2) << " is empty." << std::endl;
		if(jump_id_ != 1)
			tr.Info << "Using Jump ID " << jump_id_ << " to define surfaces." << std::endl;
	}
}

filters::FilterOP
ShapeComplementarityFilterCreator::create_filter() const { return filters::FilterOP( new ShapeComplementarityFilter ); }

std::string
ShapeComplementarityFilterCreator::keyname() const { return "ShapeComplementarity"; }


} // filters
} // protocols
