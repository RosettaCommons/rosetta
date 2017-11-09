// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/RmsdFilter.hh>
#include <protocols/protein_interface_design/filters/RmsdFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueSpanSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <algorithm>
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <protocols/toolbox/superimpose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace protein_interface_design {
namespace filters {


using namespace::protocols::toolbox;
using namespace ObjexxFCL;

RmsdFilter::RmsdFilter() :
	protocols::filters::Filter( "Rmsd" ),
	superimpose_( true ),
	symmetry_( false ),
	threshold_( 5.0 ),
	reference_pose_( /* NULL */ ),
	selection_from_segment_cache_(false),
	superimpose_on_all_( false ),
	specify_both_spans_( false ),
	CA_only_ ( true ),
	begin_native_ (0),
	end_native_ (0),
	begin_pose_ (0),
	end_pose_ ( 0 ),
	by_aln_ ( false )
{
	aln_files_.clear();
	template_names_.clear();
	query_names_.clear();
}

RmsdFilter::RmsdFilter(
	utility::vector1<core::Size> const selection,
	bool const superimpose,
	core::Real const threshold,
	core::pose::PoseOP reference_pose
) : protocols::filters::Filter( "Rmsd" ),
	selection_( core::select::residue_selector::ResidueSelectorOP( new core::select::residue_selector::ResidueIndexSelector(selection) ) ),
	superimpose_(superimpose),
	threshold_(threshold),
	reference_pose_(reference_pose),
	selection_from_segment_cache_(false),
	specify_both_spans_( false ),
	CA_only_ ( true ),
	begin_native_ (0),
	end_native_ (0),
	begin_pose_ (0),
	end_pose_ ( 0 ),
	by_aln_ (false )
{
	aln_files_.clear();
	template_names_.clear();
	query_names_.clear();
}


RmsdFilter::~RmsdFilter() {}

protocols::filters::FilterOP
RmsdFilter::clone() const {
	return protocols::filters::FilterOP( new RmsdFilter( *this ) );
}

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.RmsdFilter" );
core::Real
RmsdFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core;
	using namespace core::scoring;
	core::pose::Pose copy_pose = pose;
	core::pose::Pose native = *reference_pose_;
	core::Real rmsd( 0.0 );
	if ( by_aln_ ) {
		utility::vector1<std::string> template_seq;
		utility::vector1<std::string> query_seq;
		// read sequences from alignment files
		for ( core::Size iter = 1; iter <= aln_files_.size(); ++iter ) {
			utility::io::izstream data( aln_files_[iter] );
			query_seq.push_back(""); template_seq.push_back("");
			runtime_assert( data );
			std::string line;
			getline(data, line);
			while ( data ) {
				if ( line.length() == 0 ) {
					continue;
				}
				if ( line.substr(1, query_names_[iter].length()) == query_names_[iter] ) {
					while ( data ) {
						getline( data, line );
						if ( line[0] == '>' ) {
							break;
						}
						query_seq[iter] += line;
					}
				}
				if ( line.length() == 0 ) {
					continue;
				}
				if ( line.substr(1, template_names_[iter].length()) == template_names_[iter] ) {
					while ( data ) {
						getline( data, line );
						if ( line[0] == '>' ) {
							break;
						}
						template_seq[iter] += line;
					}
				}
				if ( line.substr(1, query_names_[iter].length()) == query_names_[iter] ) {
					while ( data ) {
						getline( data, line );
						if ( line[0] == '>' ) {
							break;
						}
						query_seq[iter] += line;
					}
				}
				while ( data ) {
					getline( data, line);
					if ( line[0] == '>' ) {
						break;
					}
				}
			}
		}
		utility::vector1<core::Size> template_segm;
		utility::vector1<core::Size> query_segm;
		core::Size template_length = 0;
		core::Size query_length = 0;
		core::Size template_prev_chain;
		core::Size query_prev_chain;
		// read through aligned sequences, count every "aligned" residue into segments
		for ( core::Size i=1; i <= aln_files_.size(); ++i ) {
			template_prev_chain = template_length;
			query_prev_chain = query_length;
			template_length = 0;
			query_length = 0;
			core::Size template_gaps = 0;
			core::Size query_gaps = 0;
			for ( core::Size pos = 0; pos < template_seq[i].size(); ++pos ) {
				if ( template_seq[i][pos] == '-' && query_seq[i][pos] == '-' ) {
					template_gaps ++;
					query_gaps ++;
				} else if ( template_seq[i][pos] != '-' && query_seq[i][pos] == '-' ) {
					query_gaps ++;
					template_length ++;
				} else if ( template_seq[i][pos] == '-' && query_seq[i][pos] != '-' ) {
					template_gaps ++;
					query_length ++;
				} else if ( template_seq[i][pos] != '-' && query_seq[i][pos] != '-' ) {
					template_segm.add_back(template_prev_chain + pos - template_gaps + 1);
					query_segm.add_back(query_prev_chain + pos - query_gaps + 1);
					template_length ++;
					query_length ++;
				}
			}
		}
		TR << "template segments to align and calculate RMSD on:\n" << template_segm << "\nquery segments to align and calculate RMSD on:\n" << query_segm << std::endl;
		TR << "Found " << template_segm.size() << " residues to calculate RMSD over" << std::endl;
		runtime_assert_msg( template_segm.size() > 10 , "there must be more that 10 residues to calculate RMSD over" );
		// align and calculate RMSD over the above calculated segments
		utility::vector1< numeric::xyzVector< core::Real > >  query_coords( Ca_coords( pose, template_segm ) );
		utility::vector1< numeric::xyzVector< core::Real > >  template_coords( Ca_coords( *reference_pose_, query_segm ) );
		numeric::xyzMatrix< core::Real > rotation;
		numeric::xyzVector< core::Real > to_init_center, to_fit_center;
		using namespace protocols::toolbox;
		superposition_transform( query_coords, template_coords, rotation, to_init_center, to_fit_center );
		apply_superposition_transform( copy_pose, rotation, to_init_center, to_fit_center );
		rmsd = res_rmsd( template_segm, query_segm, copy_pose, *reference_pose_ );
		return rmsd;
	}

	if ( specify_both_spans_ ) {

		if ( CA_only_ ) {
			if ( superimpose_ ) {
				FArray2D_double pose_coordinates_;
				FArray2D_double ref_coordinates_;
				pose_coordinates_.redimension(3, 1*(end_native_-begin_native_), 0.0);
				ref_coordinates_.redimension(3, 1*(end_native_-begin_native_), 0.0);

				core::Size n_at = 1;
				for ( core::Size i = 0; i < end_native_-begin_native_; ++i ) {
					for ( core::Size d = 1; d <= 3; ++d ) {
						pose_coordinates_(d, n_at) = copy_pose.xyz(core::id::NamedAtomID("CA", begin_pose_+i))[d - 1];
						ref_coordinates_(d, n_at) = native.xyz(core::id::NamedAtomID("CA", begin_native_+i))[d - 1];
					}
					n_at++;
				}

				rmsd = numeric::model_quality::rms_wrapper((end_native_-begin_native_), pose_coordinates_, ref_coordinates_);
			} else {
				core::Real sum2( 0.0 );
				core::Size natoms( 0 );
				for ( core::Size i = 0; i < end_native_-begin_native_; ++i ) {
					core::Vector diff = copy_pose.xyz(core::id::NamedAtomID("CA", begin_pose_+i))-native.xyz(core::id::NamedAtomID("CA", begin_native_+i));
					sum2 += diff.length_squared();
					natoms += 1;
				}
				rmsd = std::sqrt(sum2 / natoms);
			}

		} else {

			if ( superimpose_ ) {
				FArray2D_double pose_coordinates_;
				FArray2D_double ref_coordinates_;
				pose_coordinates_.redimension(3, 4*(end_native_-begin_native_), 0.0);
				ref_coordinates_.redimension(3, 4*(end_native_-begin_native_), 0.0);

				core::Size n_at = 1;
				for ( core::Size i = 0; i < end_native_-begin_native_; ++i ) {

					for ( core::Size d = 1; d <= 3; ++d ) {
						pose_coordinates_(d, n_at) = copy_pose.xyz(core::id::NamedAtomID("N", begin_pose_+i))[d - 1];
						ref_coordinates_(d, n_at) = native.xyz(core::id::NamedAtomID("N", begin_native_+i))[d - 1];
					}
					n_at++;

					for ( core::Size d = 1; d <= 3; ++d ) {
						pose_coordinates_(d, n_at) = copy_pose.xyz(core::id::NamedAtomID("CA", begin_pose_+i))[d - 1];
						ref_coordinates_(d, n_at) = native.xyz(core::id::NamedAtomID("CA", begin_native_+i))[d - 1];
					}
					n_at++;

					for ( core::Size d = 1; d <= 3; ++d ) {
						pose_coordinates_(d, n_at) = copy_pose.xyz(core::id::NamedAtomID("C", begin_pose_+i))[d - 1];
						ref_coordinates_(d, n_at) = native.xyz(core::id::NamedAtomID("C", begin_native_+i))[d - 1];
					}
					n_at++;

					for ( core::Size d = 1; d <= 3; ++d ) {
						pose_coordinates_(d, n_at) = copy_pose.xyz(core::id::NamedAtomID("O", begin_pose_+i))[d - 1];
						ref_coordinates_(d, n_at) = native.xyz(core::id::NamedAtomID("O", begin_native_+i))[d - 1];
					}
					n_at++;
				}

				rmsd = numeric::model_quality::rms_wrapper(4*(end_native_-begin_native_), pose_coordinates_, ref_coordinates_);
			} else {
				core::Real sum2( 0.0 );
				core::Size natoms( 0 );
				for ( core::Size i = 0; i < end_native_-begin_native_; ++i ) {
					core::Vector diff1 = copy_pose.xyz(core::id::NamedAtomID("N", begin_pose_+i))-native.xyz(core::id::NamedAtomID("N", begin_native_+i));
					sum2 += diff1.length_squared();
					natoms += 1;
					core::Vector diff2 = copy_pose.xyz(core::id::NamedAtomID("CA", begin_pose_+i))-native.xyz(core::id::NamedAtomID("CA", begin_native_+i));
					sum2 += diff2.length_squared();
					natoms += 1;
					core::Vector diff3 = copy_pose.xyz(core::id::NamedAtomID("C", begin_pose_+i))-native.xyz(core::id::NamedAtomID("C", begin_native_+i));
					sum2 += diff3.length_squared();
					natoms += 1;
					core::Vector diff4 = copy_pose.xyz(core::id::NamedAtomID("O", begin_pose_+i))-native.xyz(core::id::NamedAtomID("O", begin_native_+i));
					sum2 += diff4.length_squared();
					natoms += 1;
				}
				rmsd = std::sqrt(sum2 / natoms);
			}
		}


		return rmsd;
	}

	if ( !symmetry_ ) {
		runtime_assert_msg( copy_pose.size() == native.size(), "the reference pose must be the same size as the working pose" );
	}

	// generate temporary FArray
	FArray1D_bool selection_array( pose.size(), false ); // on which residues to check rmsd
	FArray1D_bool superimpose_array( pose.size(), false ); // which residues to superimpose

	if ( selection_from_segment_cache_ ) {
		core::pose::datacache::SpecialSegmentsObserver::set_farray_from_sso( superimpose_array, pose, true );
	} else {
		core::select::residue_selector::ResidueSubset selection(pose.size(), false);
		if ( selection_ ) {
			selection = selection_->apply( *reference_pose_ ); // In parse_my_tag this was the pose reference -- Is there a way to apply the selection independently to each pose?
		}
		if ( selection_ && superimpose_on_all() ) {
			core::pose::datacache::SpecialSegmentsObserver::set_farray_from_sso( superimpose_array, pose, true );
		}
		for ( core::Size ii(1); ii <= selection.size(); ++ii ) {
			if ( selection[ ii ] ) {
				selection_array[ii-1] = true; // FArray1D is 0 indexed
				superimpose_array[ii-1] = true;
			}
		}
	}

	if ( superimpose_ && !superimpose_on_all() ) {
		if ( symmetry_ ) { // SJF I haven't changed symmetry selection_array b/c I don't know which tests to use
			rmsd = sym_rmsd_with_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
		} else {
			rmsd = rmsd_with_super_subset( copy_pose, native, superimpose_array, core::scoring::is_protein_CA );
		}
	} else if ( superimpose_ && superimpose_on_all() ) {
		calpha_superimpose_pose( copy_pose, native );
		rmsd = core::scoring::rmsd_no_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
	} else {
		rmsd = core::scoring::rmsd_no_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
	}
	return rmsd;
}

bool
RmsdFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const rmsd( compute( pose ));
	TR << "RMSD over selected residues: " << rmsd ;
	if ( rmsd <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else {
		TR<<" failing." << std::endl;
		return( false );
	}
}

void
RmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd<<'\n';
}

core::Real
RmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
RmsdFilter::set_selection( core::select::residue_selector::ResidueSelectorOP const & sele ) {
	selection_ = sele;
}

void
RmsdFilter::add_selector( core::select::residue_selector::ResidueSelectorOP const & sele  ) {
	selection_ = core::select::residue_selector::OR_combine( selection_, sele );
}

void
RmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
{
	using namespace core::select::residue_selector;

	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	///against the ref pose
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map );
		TR<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<reference_pose_->size()<<" residues"<<std::endl;
	} else {
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( reference_pose ) );
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( *reference_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}
	}

	symmetry_ = tag->getOption<bool>( "symmetry", 0 );
	superimpose_on_all( tag->getOption< bool >( "superimpose_on_all", false ) );

	std::string chains = tag->getOption<std::string>( "chains", "" );
	if ( chains != "" ) {
		for ( char const & chain: chains ) {
			ResidueSelectorOP chain_select( new ChainSelector( chain ) );
			add_selector( chain_select );
		}
	}

	utility::vector0< utility::tag::TagCOP > const rmsd_tags( tag->getTags() );
	for ( utility::vector0< utility::tag::TagCOP >::const_iterator it=rmsd_tags.begin(); it!=rmsd_tags.end(); ++it ) {
		utility::tag::TagCOP const rmsd_tag = *it;
		if ( rmsd_tag->getName() == "residue" ) {
			ResidueSelectorOP res_select( new ResidueIndexSelector( core::pose::get_resnum_string( rmsd_tag ) ) );
			add_selector( res_select );
		}

		if ( rmsd_tag->getName() == "span" ) {
			std::string const & begin( core::pose::get_resnum_string( rmsd_tag, "begin_" ) );
			std::string const & end( core::pose::get_resnum_string( rmsd_tag, "end_" ) );
			ResidueSelectorOP span_select( new ResidueSpanSelector( begin, end ) );
			add_selector( span_select );
		}

		if ( rmsd_tag->getName() == "span_two" ) {
			specify_both_spans_=true;
			begin_native_ = rmsd_tag->getOption<core::Size>( "begin_native", 0 );
			end_native_ = rmsd_tag->getOption<core::Size>( "end_native", 0 );
			begin_pose_ = rmsd_tag->getOption<core::Size>( "begin_pose", 0 );
			end_pose_ = rmsd_tag->getOption<core::Size>( "end_pose", 0 );
			CA_only_ = rmsd_tag->getOption<bool>( "CA_only", true );
			runtime_assert( end_pose_ > begin_pose_ );
			runtime_assert( begin_pose_>=1);
			runtime_assert( end_native_ > begin_native_ );
			runtime_assert( begin_native_>=1);
			runtime_assert( end_native_ - begin_native_ == end_pose_ - begin_pose_ );
		}
	}

	superimpose_ = tag->getOption<bool>( "superimpose", 1 );
	threshold_ = tag->getOption<core::Real>( "threshold", 5 );
	///
	if ( tag->hasOption("by_aln") ) {
		by_aln_ = (tag->getOption<bool>("by_aln"));
		if ( by_aln_ ) {
			aln_files_ = utility::string_split(tag->getOption<std::string>("aln_files"), ',', std::string());
			template_names_ = utility::string_split(tag->getOption<std::string>("template_names"), ',', std::string());
			query_names_ = utility::string_split(tag->getOption<std::string>("query_names"), ',', std::string());
		}
	}

	if ( tag->hasOption("rms_residues_from_pose_cache") ) {
		selection_from_segment_cache_ = tag->getOption<bool>( "rms_residues_from_pose_cache", 1 );
		if ( selection_ != nullptr ) std::cerr << "Warning: in rmsd filter tag, both a span selection and the instruction to set the residues from the pose cache is given. Incompatible, defined span will be ignored." << std::endl;
	}

	TR<<"RMSD filter with superimpose=" << superimpose_ << " and threshold="<< threshold_ << std::endl;
	if ( superimpose_on_all() ) {
		TR<<" superimpose_on_all set to true. Any spans defined will be used only to measure RMSd but the pose will be supreimposed on the reference pose through all residues ";
	}
	if ( selection_from_segment_cache_ ) TR << " that are in pose segment observer cache at apply time." << std::endl;
	else if ( selection_ == nullptr ) {
		TR << "ALL" << std::endl;
		selection_ = ResidueSelectorOP( new TrueResidueSelector );
	} else {
		TR << "selected residues" << std::endl;
	}
}

core::Real atom_distance( core::conformation::Residue const & r1, std::string const & a1, core::conformation::Residue const & r2, std::string const & a2 ){
	return( r1.xyz( a1 ).distance( r2.xyz( a2 ) ) );
}

core::Real res_rmsd( utility::vector1< core::Size > const pose_res_for_rmsd, utility::vector1< core::Size > const template_res_for_rmsd, core::pose::Pose const & copy_pose, core::pose::Pose const & template_pose ){
	//two_res_rmsd( core::conformation::Residue const & r1, core::conformation::Residue const & r2 ){
	core::Real rmsd( 0.0 );
	core::Size count = 0.0;
	runtime_assert( pose_res_for_rmsd.size() == template_res_for_rmsd.size() );

	for ( core::Size i = 1; i <= pose_res_for_rmsd.size(); i++ ) {
		core::conformation::Residue r1 = copy_pose.conformation().residue( pose_res_for_rmsd[i] );
		core::conformation::Residue r2 = template_pose.conformation().residue( template_res_for_rmsd[i] );
		rmsd += pow( atom_distance( r1, "CA", r2, "CA" ), 2.0 );
		count = count + 1;
		//if ( (r1.name3()!="GLY") && (r2.name3()!="GLY") ){
		//rmsd += pow( atom_distance( r1, "CB", r2, "CB" ), 2.0 );
		//count++;
		//}
	}
	rmsd /= (core::Real) count;
	rmsd = pow( rmsd, 0.5 );
	return rmsd;
}

utility::vector1< numeric::xyzVector< core::Real > >Ca_coords( core::pose::Pose const & pose, utility::vector1< core::Size > const & positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;
	for ( core::Size const pos : positions ) {
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP RmsdFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RmsdFilter ); }

// XRW TEMP std::string
// XRW TEMP RmsdFilterCreator::keyname() const { return "Rmsd"; }

std::string RmsdFilter::name() const {
	return class_name();
}

std::string RmsdFilter::class_name() {
	return "Rmsd";
}

void RmsdFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_saved_reference_pose( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "symmetry", xsct_rosetta_bool, "Should we be symmetry-aware?", "0" )
		+ XMLSchemaAttribute( "chains", xs_string, "Chain characters without separation" )
		+ XMLSchemaAttribute::attribute_w_default( "superimpose_on_all", xsct_rosetta_bool, "Do we want to superimpose on all residues?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "superimpose", xsct_rosetta_bool, "Do we want to superimpose?", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold in RMSD above which the filter fails", "5" )
		+ XMLSchemaAttribute( "by_aln", xsct_rosetta_bool, "Perform RMSD calculation as guided by a sequence alignment file" )
		+ XMLSchemaAttribute( "aln_files", xs_string, "Comma-separated list of alignment files" )
		+ XMLSchemaAttribute( "template_names", xs_string, "Comma-separated list of template names" )
		+ XMLSchemaAttribute( "query_names", xs_string, "Comma-separated list of query names" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_residues_from_pose_cache", xs_string, "Are the RMS residues stored in the Pose's cache?", "1" );

	AttributeList rmsd_subtag_attributes;
	rmsd_subtag_attributes + XMLSchemaAttribute( "resnum", xsct_refpose_enabled_residue_number, "Residue number for the beginning of the alignment region" )
		+ XMLSchemaAttribute( "pdb_num", xsct_refpose_enabled_residue_number, "Residue number for the beginning of the alignment region" );

	AttributeList span_subtag_attributes;
	span_subtag_attributes + XMLSchemaAttribute( "begin_res_num", xsct_refpose_enabled_residue_number, "Residue number for the beginning of the alignment region" )
		+ XMLSchemaAttribute( "begin_pdb_num", xsct_refpose_enabled_residue_number, "Residue number for the beginning of the alignment region" )
		+ XMLSchemaAttribute( "end_res_num", xsct_refpose_enabled_residue_number, "Residue number for the end of the alignment region" )
		+ XMLSchemaAttribute( "end_pdb_num", xsct_refpose_enabled_residue_number, "Residue number for the end of the alignment region" );

	AttributeList spn2_subtag_attributes;
	spn2_subtag_attributes + XMLSchemaAttribute( "begin_native", xsct_non_negative_integer, "Residue number for the beginning of the alignment region, on the native" )
		+ XMLSchemaAttribute( "end_native", xsct_non_negative_integer, "Residue number for the end of the alignment region, on the native" )
		+ XMLSchemaAttribute( "begin_pose", xsct_non_negative_integer, "Residue number for the beginning of the alignment region, on the pose" )
		+ XMLSchemaAttribute( "end_pose", xsct_non_negative_integer, "Residue number for the end of the alignment region, on the pose" );

	XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "rmsd", rmsd_subtag_attributes, "XRW TODO"/*, 0 minoccurs*/ );
	ssl.add_simple_subelement( "span", span_subtag_attributes, "XRW TODO"/*, 0 minoccurs*/ );
	ssl.add_simple_subelement( "span_two", spn2_subtag_attributes, "XRW TODO"/*, 0 minoccurs*/ );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "Filter based on the C-alpha RMSD to a reference structure.", attlist, ssl );
}

std::string RmsdFilterCreator::keyname() const {
	return RmsdFilter::class_name();
}

protocols::filters::FilterOP
RmsdFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new RmsdFilter );
}

void RmsdFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RmsdFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // devel


