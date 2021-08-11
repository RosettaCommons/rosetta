// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.cc
/// @brief a vector of Poses sorted by the secondary structure of the termini
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/data_storage/TerminalDSSPSortedPoseVector.hh>
#include <protocols/pose_sewing/util.hh>

#include <core/pose/subpose_manipulation_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/pose/variant_util.hh>

#include <core/conformation/Conformation.hh>
#include <core/chemical/VariantType.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueEnergyMetric.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
//#include <core/chemical/ResidueProperties.hh>
#include <core/simple_metrics/SimpleMetric.hh>

#include <protocols/pose_sewing/data_storage/SegmentEnvelope.hh>
#include <protocols/pose_sewing/data_storage/PoseSegment.hh>
#include <protocols/pose_sewing/data_storage/PoseWithTerminalSegmentsOfKnownDSSP.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/io/GeneralFileManager.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <time.h>
#include <sstream>

static basic::Tracer TR( "protocols.pose_sewing.data_storage.TerminalDSSPSortedPoseVector" );


namespace protocols {
namespace pose_sewing {
namespace data_storage {

using core::Size;

TerminalDSSPSortedPoseVector::TerminalDSSPSortedPoseVector():
	utility::VirtualBase()
{
	utility::vector1<SegmentEnvelopeOP> seg_1;
	segment_envelopes_.push_back(seg_1);
}

TerminalDSSPSortedPoseVector::~TerminalDSSPSortedPoseVector(){}

TerminalDSSPSortedPoseVector::TerminalDSSPSortedPoseVector( TerminalDSSPSortedPoseVector const & ) = default;



TerminalDSSPSortedPoseVectorOP
TerminalDSSPSortedPoseVector::clone() const {
	return TerminalDSSPSortedPoseVectorOP( new TerminalDSSPSortedPoseVector( *this ) );
}
void
TerminalDSSPSortedPoseVector::parse_motif_file(std::string motif_filename){
	clear_segment_envelopes();
	utility::io::izstream motif_file(motif_filename);
	std::string line;
	std::string motifs;
	while ( getline( motif_file, line) ) {
		add_segment_envelopes_from_string(line);
	}

}
void
TerminalDSSPSortedPoseVector::add_segment_envelopes_from_string(std::string motifs){
	utility::vector1<std::string> tokens;
	core::Size current_token = 1;
	utility::vector1<SegmentEnvelopeOP> current_segment_envelope_vector;
	tokens = utility::split(motifs);
	if ( tokens.size()%3!=0 ) {
		utility_exit_with_message("Malformed motif file!");
	}
	while ( current_token < tokens.size() ) {
		SegmentEnvelopeOP new_envelope = SegmentEnvelopeOP(new SegmentEnvelope());
		new_envelope->set_permissible_secondary_structures(tokens[current_token]);
		new_envelope->set_minimum_length(utility::string2Size(tokens[current_token+1]));
		new_envelope->set_maximum_length(utility::string2Size(tokens[current_token+2]));
		current_segment_envelope_vector.push_back(new_envelope);
		current_token = current_token+3;
	}
	this->add_to_segment_envelopes(current_segment_envelope_vector);
	if ( segment_envelopes_[1].size() == 0 ) {
		utility_exit_with_message("No motif file created");
	}
}

void
TerminalDSSPSortedPoseVector::store_pose(PoseSegmentOP n_terminal, PoseSegmentOP c_terminal){
	if ( n_terminal->get_source_pose() != c_terminal->get_source_pose() ) {
		utility_exit_with_message("terminal segments come from different poses!");
	}
	PoseWithTerminalSegmentsOfKnownDSSPOP pose_to_store = PoseWithTerminalSegmentsOfKnownDSSPOP(new PoseWithTerminalSegmentsOfKnownDSSP());
	pose_to_store->set_N_term_DSSP(n_terminal->get_dssp_code());
	pose_to_store->set_C_term_DSSP(c_terminal->get_dssp_code());
	pose_to_store->set_N_term_length(n_terminal->get_length());
	pose_to_store->set_C_term_length(c_terminal->get_length());

	core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose());

	utility::vector1< core::Size > residues;
	for ( core::Size i = n_terminal->get_starting_residue(); i <= c_terminal->get_ending_residue(); ++i ) {
		residues.push_back(i);
	}

	core::pose::pdbslice(*pose, *n_terminal->get_source_pose(), residues);

	core::pose::set_reasonable_fold_tree(*pose);

	core::pose::PoseOP dssp_pose = core::pose::PoseOP(new core::pose::Pose(*(n_terminal->get_source_pose())));
	core::scoring::dssp::Dssp dssp(*dssp_pose);
	dssp.insert_ss_into_pose(*dssp_pose);

	std::string secstruct;
	for ( core::Size current_residue = 1; current_residue <= pose->size(); ++current_residue ) {
		pose->set_secstruct(current_residue, dssp_pose->secstruct(current_residue-1+n_terminal->get_starting_residue()));
		secstruct+=dssp_pose->secstruct(current_residue-1+n_terminal->get_starting_residue());
	}
	pose_to_store->set_secstruct(secstruct);
	for ( core::Size conf_resnum=1; conf_resnum<pose->size(); ++conf_resnum ) {
		pose->conformation_ptr()->update_polymeric_connection(conf_resnum);
	}


	core::pose::remove_upper_terminus_type_from_pose_residue(*pose,1);
	core::pose::remove_lower_terminus_type_from_pose_residue(*pose,pose->size());

	pose_to_store->set_source_pose(pose);
	this->store_pose(pose_to_store);
}
void
TerminalDSSPSortedPoseVector::store_pose(PoseWithTerminalSegmentsOfKnownDSSPOP pose_to_store){

	if ( nterm_secstructs_.find(pose_to_store->get_N_term_DSSP())==std::string::npos ) {
		nterm_secstructs_ = nterm_secstructs_+pose_to_store->get_N_term_DSSP();
	}
	if ( cterm_secstructs_.find(pose_to_store->get_C_term_DSSP())==std::string::npos ) {
		cterm_secstructs_ = cterm_secstructs_+pose_to_store->get_C_term_DSSP();
	}

	poseOP_map_[pose_to_store->get_N_term_DSSP()][pose_to_store->get_C_term_DSSP()].push_back(pose_to_store);
}

PoseWithTerminalSegmentsOfKnownDSSPOP
TerminalDSSPSortedPoseVector::get_random_pose(std::string nterm_secstructs,std::string cterm_secstructs){
	core::Size sum_of_vector_sizes = 0;
	for ( char n_secstruct:nterm_secstructs ) {
		for ( char c_secstruct:cterm_secstructs ) {
			sum_of_vector_sizes+=poseOP_map_[n_secstruct][c_secstruct].size();
		}
	}
	if ( sum_of_vector_sizes == 0 ) {
		TR << " N secstructs are: " << nterm_secstructs << std::endl;
		TR << " C secstructs are: " << cterm_secstructs << std::endl;
		utility_exit_with_message("No permissible segments found. Check permissible segment ends.");
	}
	core::Size selected_entry = numeric::random::random_range(1,sum_of_vector_sizes);

	for ( char n_secstruct:nterm_secstructs ) {
		for ( char c_secstruct:cterm_secstructs ) {
			if ( selected_entry > poseOP_map_[n_secstruct][c_secstruct].size() ) {
				selected_entry -= poseOP_map_[n_secstruct][c_secstruct].size();
			} else {
				return poseOP_map_[n_secstruct][c_secstruct].at(selected_entry);
			}
		}
	}
	utility_exit_with_message("Miscalculated pose index!");
	return nullptr;
}

utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
TerminalDSSPSortedPoseVector::get_vector(std::string nterm_secstructs,std::string cterm_secstructs){
	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP> out_vector;

	for ( char n_secstruct:nterm_secstructs ) {
		for ( char c_secstruct:cterm_secstructs ) {
			out_vector=out_vector.append(poseOP_map_[n_secstruct][c_secstruct]);
		}
	}

	return out_vector;
}
void
TerminalDSSPSortedPoseVector::populate_from_pdb_list(std::string pdb_list){
	if ( segment_envelopes_[1].size() == 0 ) {
		utility_exit_with_message("No motif file loaded");
	}
	utility::io::izstream pdbs(pdb_list);
	std::string line;
	core::pose::PoseOP current_pose;
	while ( getline( pdbs, line) ) {
		utility::vector1<PoseSegmentOP> pose_segments;
		current_pose = core::import_pose::pose_from_file(line);
		core::scoring::dssp::Dssp dssp(*current_pose);
		dssp.insert_ss_into_pose(*current_pose);
		//split dssp string into blocks
		core::Size last_residue = 1;
		char last_dssp = current_pose->secstruct(last_residue);
		for ( core::Size current_residue=2; current_residue<current_pose->size(); ++current_residue ) { // not a foreach because we need to manipulate current_residue numerically
			if ( current_pose->secstruct(current_residue)!=last_dssp ) {
				pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_residue-1,last_dssp,current_pose)));
				last_residue = current_residue;
				last_dssp = current_pose->secstruct(current_residue);
			}

		}
		pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_pose->size(),last_dssp,current_pose)));
		//now compare blocks to segment envelopes
		for ( core::Size current_motif_number=1; current_motif_number <= segment_envelopes_.size(); ++current_motif_number ) {
			utility::vector1<SegmentEnvelopeOP> current_motif = segment_envelopes_[current_motif_number];
			for ( core::Size offset = 0; offset <= pose_segments.size() - current_motif.size(); ++offset ) {
				core::Size segment_to_compare=1;
				bool good_match = true;
				while ( good_match ) {
					if ( current_motif[segment_to_compare]->is_valid(pose_segments[segment_to_compare+offset]) ) {
						if ( segment_to_compare == current_motif.size() ) {
							//match is good!
							this->store_pose(pose_segments[offset+1],pose_segments[offset+segment_to_compare]);
							good_match = false;
						} else { //keep going
							++segment_to_compare;
						}
					} else { // advance forward
						good_match = false;
					}
				}
			}
		}

	}

}
void
TerminalDSSPSortedPoseVector::populate_from_pdb_list(std::string motif_filename, std::string pdb_list){

	this->parse_motif_file(motif_filename);
	this->populate_from_pdb_list(pdb_list);
}
void
TerminalDSSPSortedPoseVector::set_segment_envelopes(utility::vector1<utility::vector1<SegmentEnvelopeOP>> segment_envelopes){
	segment_envelopes_ = segment_envelopes;
}
void
TerminalDSSPSortedPoseVector::add_to_segment_envelopes(utility::vector1<SegmentEnvelopeOP> segment_envelope){
	segment_envelopes_.push_back(segment_envelope);
}

void
TerminalDSSPSortedPoseVector::clear_segment_envelopes(){
	segment_envelopes_.clear();
}

void
TerminalDSSPSortedPoseVector::populate_from_segment_file(std::string file_path, bool clear /*false*/) {
	if ( clear ) {
		poseOP_map_.clear();
	}

	std::string segment_file_path = file_path + "/segment.txt";
	if ( ! utility::file::file_exists(segment_file_path) ) {
		utility_exit_with_message("Segment file "+segment_file_path+" does not exist!");
	}

	utility::io::izstream segment_file( segment_file_path );
	std::string line;
	utility::vector1<std::string> tokens;
	while ( getline( segment_file, line) ) {
		tokens = utility::split(line);
		PoseWithTerminalSegmentsOfKnownDSSPOP current_pose_segment = PoseWithTerminalSegmentsOfKnownDSSPOP( new PoseWithTerminalSegmentsOfKnownDSSP);
		current_pose_segment->set_N_term_DSSP(tokens[1].at(0));
		current_pose_segment->set_N_term_length(utility::string2Size(tokens[2]));
		current_pose_segment->set_C_term_DSSP(tokens[3].at(0));
		current_pose_segment->set_C_term_length(utility::string2Size(tokens[4]));
		current_pose_segment->set_secstruct(tokens[5]);
		current_pose_segment->set_filename(file_path + "/" + tokens[6]);
		current_pose_segment->set_segfile_path(file_path);
		this->store_pose(current_pose_segment);
	}
	segment_file.close();
}

utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
TerminalDSSPSortedPoseVector::populate_from_segment_file_and_get_random_vector_set(
	std::string const & path,
	std::string const & allowed_n_ter_ss,
	std::string const & allowed_c_ter_ss,
	core::Size max_vector_size /* 5000 */
) const {
	utility::vector1< std::string > paths;
	paths.push_back(path);
	return populate_from_segment_file_and_get_random_vector_set( paths, allowed_n_ter_ss, allowed_c_ter_ss, max_vector_size);
}

utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
TerminalDSSPSortedPoseVector::populate_from_segment_file_and_get_random_vector_set(
	utility::vector1< std::string > const & file_paths,
	std::string const & allowed_n_ter_ss,
	std::string const & allowed_c_ter_ss,
	core::Size max_vector_size /* 5000 */
) const {

	std::string line;
	utility::vector1<std::string> tokens;
	utility::vector1< std::string > passed_lines;
	//TR << "Allowed nter_ss " << allowed_n_ter_ss << std::endl;
	//TR << "Allowed cter_ss " << allowed_c_ter_ss << std::endl;
	//No way to randomize this, so need to do this twice.
	for ( std::string const & file_path : file_paths ) {
		std::string const segment_file_path = file_path + "/segment.txt";
		//TR << "Reading " <<segment_file_path << std::endl;
		if ( ! utility::file::file_exists(segment_file_path) ) {
			utility_exit_with_message("Segment file "+segment_file_path+" does not exist!");
		}


		//Only read the file once per thread.
		utility::vector1< std::string > const & lines = utility::io::GeneralFileManagerVector::get_instance()->get_file_contents(segment_file_path);

		//utility::io::izstream segment_file( segment_file_path );
		for ( auto const & line : lines ) {
			tokens = utility::split(line);
			char const N_term_DSSP = tokens[1].at(0);
			char const C_term_DSSP = tokens[3].at(0);
			if ( (allowed_n_ter_ss.find(N_term_DSSP) != std::string::npos) && allowed_c_ter_ss.find(C_term_DSSP) != std::string::npos ) {
				passed_lines.push_back( line+" "+file_path);
			}
		}
		//segment_file.close();
	}

	numeric::random::random_permutation(passed_lines);
	TR << "Total Segments read from seg file: " << passed_lines.size() << std::endl;
	utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP> populated_segments;

	for ( std::string const & current_line : passed_lines ) {
		if ( populated_segments.size() >= max_vector_size ) break;

		tokens = utility::split(current_line);
		if ( tokens.size() < 6 ) continue;
		PoseWithTerminalSegmentsOfKnownDSSPOP current_pose_segment = PoseWithTerminalSegmentsOfKnownDSSPOP( new PoseWithTerminalSegmentsOfKnownDSSP);
		current_pose_segment->set_N_term_DSSP(tokens[1].at(0));
		current_pose_segment->set_N_term_length(utility::string2Size(tokens[2]));
		current_pose_segment->set_C_term_DSSP(tokens[3].at(0));
		current_pose_segment->set_C_term_length(utility::string2Size(tokens[4]));
		current_pose_segment->set_secstruct(tokens[5]);
		current_pose_segment->set_filename(tokens[tokens.size()] + "/" + tokens[6]);
		current_pose_segment->set_segfile_path(tokens[tokens.size()]);
		populated_segments.push_back(current_pose_segment);
	}

	TR << "Total Segments returned: " << populated_segments.size() << std::endl;
	return populated_segments;
}

void
TerminalDSSPSortedPoseVector::populate_all_OP_vector(){

	all_poseOPs_.clear();

	for ( char n_secstruct:nterm_secstructs_ ) {
		for ( char c_secstruct:cterm_secstructs_ ) {
			for ( PoseWithTerminalSegmentsOfKnownDSSPOP current_poseOP : poseOP_map_[n_secstruct][c_secstruct] ) {
				all_poseOPs_.push_back(current_poseOP);
			}
		}
	}

}

utility::vector1<PoseWithTerminalSegmentsOfKnownDSSPOP>
TerminalDSSPSortedPoseVector::get_all_poseOPs(){
	this->populate_all_OP_vector();
	return all_poseOPs_;
}
void
TerminalDSSPSortedPoseVector::clear_all_vectors(){
	all_poseOPs_.clear();
	this->clear_terminal_vectors();
}
void
TerminalDSSPSortedPoseVector::clear_terminal_vectors(){
	for ( char n_secstruct:nterm_secstructs_ ) {
		for ( char c_secstruct:cterm_secstructs_ ) {
			poseOP_map_[n_secstruct][c_secstruct].clear();
		}
	}
}


void
TerminalDSSPSortedPoseVector::create_segment_file(std::string file_path) {
	this->populate_all_OP_vector();
	utility::io::ozstream out_segment_stream;
	out_segment_stream.open(file_path+"/segment.txt", std::ios::out);
	core::Size current_OP_out = 1;
	std::string pdb_filename;
	for ( PoseWithTerminalSegmentsOfKnownDSSPOP current_poseOP : all_poseOPs_ ) {
		pdb_filename = std::to_string(current_OP_out) + ".pdb";
		out_segment_stream << current_poseOP->get_N_term_DSSP() << " " << current_poseOP->get_N_term_length() << " " << current_poseOP->get_C_term_DSSP() << " " << current_poseOP->get_C_term_length() << " " << current_poseOP->get_secstruct() << " " << pdb_filename << std::endl;
		current_poseOP->store_source_pose_for_segment(file_path+"/"+pdb_filename);
		++current_OP_out;
	}
	out_segment_stream.close();
}

void
TerminalDSSPSortedPoseVector::simultaneously_populate_and_write_segment_file_pdb (
	core::pose::Pose const & pose,
	std::string const & motif_filename,
	std::string const & file_path,
	bool store_reference_pdbs,
	bool use_pdbs_only,
	utility::vector1< protocols::filters::FilterOP > const & filters,
	utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
	bool output_elements
) {
	this->parse_motif_file(motif_filename);
	this->simultaneously_populate_and_write_segment_file_pdb_using_stored_motifs(pose, file_path, store_reference_pdbs, use_pdbs_only, filters, metrics, output_elements);
}

void
TerminalDSSPSortedPoseVector::simultaneously_populate_and_write_segment_file_pdb_using_stored_motifs (
	core::pose::Pose const & pose,
	std::string const & file_path,
	bool store_reference_pdbs,
	bool use_pdbs_only,
	utility::vector1< protocols::filters::FilterOP > const & filters,
	utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
	bool output_elements
) {

	core::pose::PoseOP current_pose = utility::pointer::make_shared< core::pose::Pose >(pose); //We require OPs here and we will be copying here or every place we need an op.

	utility::vector1<PoseSegmentOP> pose_segments;
	utility::vector1<PoseSegmentOP> segments_to_store;
	PoseSegmentOP segment_to_consider;

	//Wipe PDBInfo, insert DSSP, fix DSSP, etc.
	current_pose->pdb_info() = nullptr;
	core::scoring::dssp::Dssp dssp(*current_pose);
	dssp.insert_ss_into_pose(*current_pose);


	//Prepare outputs
	std::string pose_fname = pose.pdb_info()->name();
	utility::file::FileName fn(pose_fname);

	std::string nstruct_tag = "_0001";
	std::string pose_basename = utility::remove_from_string(fn.base(), nstruct_tag);

	//Here we open the file and append it if needed.  This allows us to run in MPI mode.
	utility::io::ozstream out_segment_stream;
	std::stringstream header; //Currently empty string, but in future, we want columns here as column and ignoring of comments.



	core::Size current_OP_out = 1;
	std::string pdb_filename;
	std::string mmtf_filename;
	if ( segment_envelopes_[1].size() == 0 ) {
		utility_exit_with_message("No motif file loaded");
	}

	//split dssp string into blocks
	core::Size last_residue = 1;
	char last_dssp = current_pose->secstruct(last_residue);
	for ( core::Size current_residue=2; current_residue <= current_pose->size(); ++current_residue ) { // not a foreach because we need to manipulate current_residue numerically
		if ( current_pose->secstruct(current_residue)!=last_dssp ) {
			pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_residue-1,last_dssp,current_pose)));
			last_residue = current_residue;
			last_dssp = current_pose->secstruct(current_residue);
		}

	}
	pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_pose->size(),last_dssp,current_pose)));
	//now compare blocks to segment envelopes
	for ( core::Size current_motif_number=1; current_motif_number <= segment_envelopes_.size(); ++current_motif_number ) {
		utility::vector1<SegmentEnvelopeOP> current_motif = segment_envelopes_[current_motif_number];
		if ( pose_segments.size() >= current_motif.size() ) {
			for ( core::Size offset = 0; offset < (pose_segments.size() - current_motif.size()+1); ++offset ) {
				segments_to_store.clear();
				core::Size segment_to_compare=1;
				bool good_match = true;
				while ( good_match ) {
					segment_to_consider = pose_segments[segment_to_compare+offset];
					if ( segment_to_consider->get_length() > current_motif[segment_to_compare]->get_maximum_length() ) {
						if ( segment_to_compare == 1 ) {
							segment_to_consider = PoseSegmentOP(new PoseSegment(*segment_to_consider,segment_to_consider->get_ending_residue()-(current_motif[segment_to_compare]->get_maximum_length()-1),segment_to_consider->get_ending_residue()));
						} else if ( segment_to_compare == current_motif.size() ) {
							segment_to_consider = PoseSegmentOP(new PoseSegment(*segment_to_consider,segment_to_consider->get_starting_residue(), segment_to_consider->get_starting_residue()+(current_motif[segment_to_compare]->get_maximum_length()-1)));
						}
					}


					if ( current_motif[segment_to_compare]->is_valid(segment_to_consider) ) {
						segments_to_store.push_back(segment_to_consider);
						if ( segment_to_compare == current_motif.size() ) {
							TR << "FOUND GOOD MATCH" << std::endl;
							//match is good!
							this->store_pose(segments_to_store[1],segments_to_store[segments_to_store.size()]);
							good_match = false;
							segments_to_store.clear();
						} else { //keep going
							++segment_to_compare;
						}
					} else { // advance forward
						TR << "FOUND BAD MATCH" << std::endl;
						good_match = false;
						segments_to_store.clear();
					}
				}
			}
		}
	}
	this->purge_missing_density();
	this->populate_all_OP_vector();
	utility::vector1< std::string > lines;

	core::Size failed_filter_counts = 0;
	core::Size total_segments = 0;
	for ( PoseWithTerminalSegmentsOfKnownDSSPOP current_poseOP : all_poseOPs_ ) {

		core::pose::PoseOP out_pose = current_poseOP->create_source_pose_for_segment(output_elements);
		if ( filters.size() > 0 ) {
			bool filter_failed = false;
			for ( const std::shared_ptr<protocols::filters::Filter> & filter: filters ) {
				if ( ! filter->apply(*out_pose) ) {
					failed_filter_counts+=1;
					filter_failed = true;
					TR << "Filter failed!" << std::endl;
					break; //Break out of filters
				}
			}
			if ( filter_failed ) {
				TR << "Filter failed, onto next segment." << std::endl;
				continue; //Next segment
			}
		}
		if ( metrics.size() > 0 ) {
			for ( const std::shared_ptr<core::simple_metrics::SimpleMetric> & metric: metrics ) {
				metric->apply(*out_pose);
			}
		}
		total_segments+=1;

		std::string segment_pdb_path;
		std::string filename_string;

		mmtf_filename = pose_basename + "_" + std::to_string(current_OP_out) + ".mmtf";
		pdb_filename = pose_basename + "_" + std::to_string(current_OP_out) + ".pdb";

		if ( use_pdbs_only ) {
			segment_pdb_path = pdb_filename;
			filename_string = pdb_filename;
		} else {
			segment_pdb_path = mmtf_filename;
			filename_string = mmtf_filename;
		}

		TR << "Storing pose "<< pdb_filename << std::endl;
		std::stringstream line;
		line << current_poseOP->get_N_term_DSSP() << " " << current_poseOP->get_N_term_length() << " " << current_poseOP->get_C_term_DSSP() << " " << current_poseOP->get_C_term_length() << " " << current_poseOP->get_secstruct() << " " << segment_pdb_path << " " << fn.path() << std::endl;

		out_pose->dump_pdb(file_path+"/"+filename_string );


		if ( store_reference_pdbs && ! use_pdbs_only ) {
			current_poseOP->store_reference_pdb(file_path+"/"+pdb_filename);
		}

		lines.push_back(line.str());
		++current_OP_out;
	}
	out_segment_stream.open_append_if_existed(file_path+"/segment.txt", header);
	for ( std::string const & line : lines ) {
		out_segment_stream << line;
	}
	out_segment_stream.close();
	this->clear_all_vectors();
	if ( filters.size() ) {
		TR << "A total of " << failed_filter_counts << " segments failed at least one filter." << std::endl;
	}
	TR << " A total of " << total_segments << " segments were output from this pdb." <<std::endl;
}

void
TerminalDSSPSortedPoseVector::simultaneously_populate_and_write_segment_file (
	std::string const & motif_filename,
	std::string const & pdb_list,
	std::string const & file_path,
	bool store_reference_pdbs,
	bool use_pdbs_only,
	utility::vector1< protocols::filters::FilterOP > const & filters,
	utility::vector1< core::simple_metrics::SimpleMetricOP > const & metrics,
	bool output_elements) {

	utility::io::izstream pdbs(pdb_list);
	std::string line;
	core::pose::PoseOP current_pose;
	while ( getline( pdbs, line) ) {
		try {
			current_pose = core::import_pose::pose_from_file(line);
		}
catch ( utility::excn::Exception& excn ) {
	std::cerr << "Exception : " << std::endl;
	excn.show( std::cerr );
	excn.show( TR );
	TR << "caught exception.  skipping graft. printing structures." << std::endl;
}

		simultaneously_populate_and_write_segment_file_pdb(*current_pose, motif_filename, file_path, store_reference_pdbs, use_pdbs_only, filters, metrics, output_elements);
	}

}

core::Size
TerminalDSSPSortedPoseVector::get_total_size () {
	core::Size total_size = 0;
	for ( char n_secstruct:nterm_secstructs_ ) {
		for ( char c_secstruct:cterm_secstructs_ ) {
			total_size += poseOP_map_[n_secstruct][c_secstruct].size();
		}
	}
	return total_size;
}

void
TerminalDSSPSortedPoseVector::purge_missing_density() {

	this->populate_all_OP_vector();
	this->clear_terminal_vectors();
	for ( PoseWithTerminalSegmentsOfKnownDSSPOP current_posesegOP : all_poseOPs_ ) {
		bool store_pose = true;
		core::pose::PoseCOP current_pose = current_posesegOP->get_source_pose();
		for ( core::Size current_resnum = 1; current_resnum < current_pose->size(); ++current_resnum ) {
			if ( !store_pose ) {
				break;
			}
			numeric::xyzVector<core::Real> this_CA = current_pose->residue(current_resnum).xyz(2);
			numeric::xyzVector<core::Real> next_CA = current_pose->residue(current_resnum+1).xyz(2);
			if ( this_CA.distance(next_CA) > 5.0 ) {
				TR << "break found, purging!" << std::endl;
				store_pose = false;
				break;
			}
			if ( current_pose->residue(current_resnum).name3() == "CYD" ) {
				TR << "disulfide found, purging!" << std::endl;
				store_pose = false;
				break;
			}
			if ( !(current_pose->residue(current_resnum).is_bonded(current_pose->residue(current_resnum+1))) ) {
				TR << "broken bond found, purging!" << std::endl;
				store_pose = false;
				break;
			} else if ( !(current_pose->residue(current_resnum+1).is_bonded(current_pose->residue(current_resnum))) ) {
				TR << "unsymmetrical bond found, purging!" << std::endl;
				store_pose = false;
				break;
			}
		}
		for ( core::Size current_resnum = 1; current_resnum <= current_pose->size(); ++current_resnum ) {
			if ( !store_pose ) {
				break;
			}
			if ( !(current_pose->residue(current_resnum).type().is_canonical_aa()) ) {
				TR << "noncanonical found! PURGING" << std::endl;
				store_pose = false;
				break;
			}
			if ( current_pose->residue(current_resnum).annotated_name(true)!=current_pose->residue(current_resnum).annotated_name(false) ) {
				TR << "VARIANT FOUND PURGING" << std::endl;
				store_pose = false;
				break;
			}
			//hbond
			TR <<"Running hbond-dependant purge." << std::endl;
			core::scoring::ScoreFunctionOP scorefxn_;
			scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("beta_nov16");
			core::pose::PoseOP copy_pose = current_pose->clone();
			core::scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn_->energy_method_options() );
			energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
			scorefxn_->set_energy_method_options( energymethodoptions );
			scorefxn_->score( *copy_pose );
			core::scoring::hbonds::HBondSet hbond_set;
			hbond_set.setup_for_residue_pair_energies( *copy_pose, false, false );

			std::set< Size > donor_set;
			std::set< Size > acceptor_set;

			TR << "Ready for scoring" << std::endl;
			for ( core::Size i = 1; i <= hbond_set.nhbonds(); ++i ) {
				core::scoring::hbonds::HBond const & hbond = hbond_set.hbond(i);
				if ( hbond.don_hatm_is_backbone() && hbond.acc_atm_is_backbone() ) {
					acceptor_set.insert(hbond.acc_res());
					donor_set.insert(hbond.don_res());
				}
			}

			if ( current_pose->size() > 4 ) {
				for ( core::Size current_resnum = 5; current_resnum <= current_pose->size()-4; ++current_resnum ) {
					if ( current_pose->secstruct(current_resnum) == 'H' ) {
						if ( current_pose->secstruct(current_resnum-4)=='H' && donor_set.find(current_resnum) == donor_set.end() ) {
							TR << "helix break found, purging!" << std::endl;
							store_pose = false;
							break;
						}
						if ( current_pose->secstruct(current_resnum+4)=='H' && acceptor_set.find(current_resnum) == acceptor_set.end() ) {
							TR << "helix break found, purging!" << std::endl;
							store_pose = false;
							break;
						}
					}
				}
			}
			TR <<"Sidechain stuff" << std::endl;
			if ( current_pose->residue(current_resnum).has_variant_type(core::chemical::SIDECHAIN_CONJUGATION) ) {
				TR << "SIDECHAIN CONJUGATION FOUND PURGING" << std::endl;
				store_pose = false;
				break;
			}
			if ( current_pose->residue(current_resnum).has_variant_type(core::chemical::DIMETHYLATION) ) {
				TR << "DIMETHYLATION FOUND PURGING" << std::endl;
				store_pose = false;
				break;
			}
			if ( current_pose->residue(current_resnum).has_variant_type(core::chemical::ACETYLATION) ) {
				TR << "DIMETHYLATION FOUND PURGING" << std::endl;
				store_pose = false;
				break;
			}
			if ( current_pose->residue(current_resnum).name3()=="UDP" ) {
				TR << "UDP FOUND PURGING" << std::endl;
				store_pose = false;
				break;
			}

			//-hbond
			// core::simple_metrics::per_residue_metrics::PerResidueEnergyMetricOP rama_trip = core::simple_metrics::per_residue_metrics::PerResidueEnergyMetricOP(new core::simple_metrics::per_residue_metrics::PerResidueEnergyMetric());
			// rama_trip->set_scorefunction(scorefxn_);
			// rama_trip->set_scoretype(core::scoring::score_type_from_name("rama_prepro"));
			// std::map< core::Size, core::Real > ramas = rama_trip->calculate(*copy_pose);

			// for(auto current_pair : ramas){
			//  if(current_pair.second > 1) {
			//   TR << "bad rama " << current_pair.second << " found at residue " << current_pair.first << ", purging" << std::endl;
			//   store_pose = false;
			//   break;
			//  }
			// }
			// rama_trip->set_scoretype(core::scoring::score_type_from_name("omega"));
			// std::map< core::Size, core::Real > omegas = rama_trip->calculate(*copy_pose);

			// for(auto current_pair : omegas){
			//  if(current_pair.second > 1 && current_pair.first != 1) {
			//   TR << "bad omega " << current_pair.second << " found at residue " << current_pair.first << ", purging" << std::endl;
			//   store_pose = false;
			//   break;
			//  }
			// }

		}
		if ( store_pose ) {
			this->store_pose(current_posesegOP);
		}
	}
	all_poseOPs_.clear();
}
utility::vector1<core::pose::PoseCOP>
TerminalDSSPSortedPoseVector::convert_pose_into_segment_vector(core::pose::PoseOP current_pose){
	utility::vector1<core::pose::PoseCOP> return_vector;
	this->clear_all_vectors();
	utility::vector1<PoseSegmentOP> pose_segments;
	core::scoring::dssp::Dssp dssp(*current_pose);
	dssp.insert_ss_into_pose(*current_pose);
	//split dssp string into blocks
	core::Size last_residue = 1;
	char last_dssp = current_pose->secstruct(last_residue);
	for ( core::Size current_residue=2; current_residue <= current_pose->size(); ++current_residue ) { // not a foreach because we need to manipulate current_residue numerically
		if ( current_pose->secstruct(current_residue)!=last_dssp ) {
			pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_residue-1,last_dssp,current_pose)));
			last_residue = current_residue;
			last_dssp = current_pose->secstruct(current_residue);
		}

	}
	pose_segments.push_back(PoseSegmentOP(new PoseSegment(last_residue,current_pose->size(),last_dssp,current_pose)));
	//now compare blocks to segment envelopes
	for ( core::Size current_motif_number=1; current_motif_number<=segment_envelopes_.size(); ++current_motif_number ) {
		utility::vector1<SegmentEnvelopeOP> current_motif = segment_envelopes_[current_motif_number];
		for ( core::Size offset = 0; offset < pose_segments.size()-current_motif.size()+1; ++offset ) {
			core::Size segment_to_compare=1;
			bool good_match = true;
			while ( good_match ) {
				if ( current_motif[segment_to_compare]->is_valid(pose_segments[segment_to_compare+offset]) ) {
					if ( segment_to_compare == current_motif.size() ) {
						//match is good!
						this->store_pose(pose_segments[offset+1],pose_segments[offset+segment_to_compare]);
						good_match = false;
					} else { //keep going
						++segment_to_compare;
					}
				} else { // advance forward
					good_match = false;
				}
			}
		}
	}
	this->purge_missing_density();
	this->populate_all_OP_vector();

	for ( PoseWithTerminalSegmentsOfKnownDSSPOP current_posesegOP : all_poseOPs_ ) {
		return_vector.push_back(current_posesegOP->get_source_pose());
	}
	return return_vector;
}

} //protocols
} //pose_sewing
} //data_storage






