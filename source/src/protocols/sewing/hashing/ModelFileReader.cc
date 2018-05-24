// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/ModelFileReader.cc
/// @brief a reader of SEWING model files
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <core/conformation/Atom.hh>
static basic::Tracer TR( "protocols.sewing.hashing.ModelFileReader" );


namespace protocols {
namespace sewing {
namespace hashing {

ModelFileReader::ModelFileReader():
	utility::pointer::ReferenceCount()
{

}

ModelFileReader::~ModelFileReader(){}

ModelFileReader::ModelFileReader( ModelFileReader const & ) {

}



ModelFileReaderOP
ModelFileReader::clone() const {
	return ModelFileReaderOP( new ModelFileReader( *this ) );
}

SegmentVectorOP
ModelFileReader::read_model_file(std::string filename){

	core::Size starttime = time(NULL);
	utility::io::izstream file(filename);
	if ( !file.good() ) {
		utility_exit_with_message("Could not find segment file with name: " + filename);
	}
	TR << "Reading segments from file: " << filename << std::endl;
	//Starting to declare the vars we need
	std::string line;
	core::Size version = 0;
	SegmentVectorOP segment_list = SegmentVectorOP(new SegmentVector());

	data_storage::SmartSegmentOP current_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment);
	data_storage::SmartSegmentOP last_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment);
	data_storage::SmartSewingResidueOP current_residue = data_storage::SmartSewingResidueOP(new data_storage::SmartSewingResidue);
	utility::vector1<data_storage::SmartSewingResidueOP> residues_vector;
	bool c_fixed = false;
	utility::vector1<core::conformation::Atom> atoms_vector;
	utility::vector1<core::Real> chi_angles_;
	//core::Size current_segmentID = 0;
	//So the first line will be the version--unless there are possibly comments above it?
	// if( getline( file, line ) ) {
	//  utility::trim( line );
	while ( getline( file, line ) ) {

		//skip comment lines and empty lines
		utility::trim(line);
		if ( line.length() == 0 || line[0] == '#' ) {
			continue;
		}
		utility::vector1<std::string> tokens = utility::string_split(line);
		runtime_assert(tokens.size() > 0);
		if ( tokens[1]=="VERSION" ) {
			version += utility::string2int(tokens[2]);
			//  }
			//  else if ( tokens[1]=="MODEL" ) {
			////   TR << "Found model" << std::endl;
			//   //our current segment has no c-terminal neighbor
			//   current_segment -> set_c_terminal_neighbor(nullptr);
			//   if(residues_vector.size() > 0){
			//    current_segment->set_residue_vector(residues_vector);
			//    current_segment->set_length(residues_vector.size());
			//    residues_vector.clear();
			//    segment_list->push_back(current_segment);//add it to the list
			//   }
			//   current_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment);
			//   current_segment->set_segment_id(segment_list->size()+1);// do we really need this anymore?
			//
		} else if ( tokens[1]=="SEGMENT" ) {
			// if (atoms_vector.size() > 0){
			//  residues_vector.push_back(current_residue);
			//  atoms_vector.clear();
			// }
			//   TR << "Found segment" << std::endl;
			if ( residues_vector.size() > 0 ) {
				current_segment->set_residue_vector(residues_vector);
				current_segment->set_length(residues_vector.size());
				residues_vector.clear();
				segment_list->push_back(current_segment);//add it to the list
				last_segment = current_segment; //hand over the new segment
				current_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment); //make a new segment
				if ( tokens[3] == "N_FIXD" ) {
					if ( c_fixed ) {
						data_storage::SmartSegment::link_to( *(segment_list->rbegin()),current_segment); // set up the dual links to it
					} else {
						utility_exit_with_message( "ERROR: Bad connectivity in input segment file!" );
					}
				}
			}
			if ( tokens[4] == "C_FIXD" ) {
				c_fixed = true;
			} else {
				c_fixed = false;
			}
			current_segment->set_segment_id(segment_list->size()+1);// do we really need this anymore?
			current_segment->set_dssp_code(tokens[2][0]);
		} else if ( tokens[1]=="RESIDUE" ) {
			//   TR << "Found residue" << std::endl;
			current_residue->set_atom_vector(atoms_vector);
			atoms_vector.clear();
			current_residue = data_storage::SmartSewingResidueOP(new data_storage::SmartSewingResidue);
			residues_vector.push_back(current_residue);
			current_residue->set_amino_acid_type(tokens[2]);// we need a dictionary here
			current_residue->set_full_type_name(tokens[3]);
			for ( core::Size this_chi = 4; this_chi <= tokens.size(); this_chi++ ) {
				chi_angles_.push_back(utility::string2float(tokens[this_chi]));
			}
			current_residue->set_chis(chi_angles_);
			chi_angles_.clear();
		} else if ( tokens[1]=="ATOM" ) {
			//   TR << "Found atom" << std::endl;
			core::Real x = utility::string2float(tokens[2]);
			core::Real y = utility::string2float(tokens[3]);
			core::Real z = utility::string2float(tokens[4]);
			core::conformation::Atom current_atom;
			current_atom.xyz(numeric::xyzVector<core::Real>(x,y,z));
			atoms_vector.push_back(current_atom);
		} else {
			TR << line << std::endl;
		}
	}
	//Done looping
	//Now add back in the last segment
	if ( atoms_vector.size() > 0 ) {
		current_residue->set_atom_vector( atoms_vector );
	}
	if ( residues_vector.size() > 0 ) {
		current_segment->set_residue_vector(residues_vector);
		current_segment->set_length(residues_vector.size());
		//residues_vector.clear();
		//last_segment = current_segment; //hand over the new segment
		//current_segment = data_storage::SmartSegmentOP(new data_storage::SmartSegment); //make a new segment
		//if (tokens[3] == "N_FIXD"){
		//if( c_fixed ){ //The last segment's C-terminus was fixed
		//This should already be connected
		//data_storage::SmartSegment::link_to( *(segment_list->rbegin()),current_segment); // set up the dual links to it
		//}
		segment_list->push_back(current_segment);//add it to the list
	}
	core::Size endtime = time(NULL);
	TR << "Read " << segment_list->size() << " segments in " << endtime - starttime << " seconds" << std::endl;
	//for(core::Size segcounter = 1; segcounter <= segment_list->size(); segcounter++){
	// current_segment = segment_list->at(segcounter);
	// TR << " Current segment is " << segcounter << std::endl;
	// for(core::Size rescounter = 1; rescounter <= current_segment->get_residue_vector().size(); rescounter++){
	//  current_residue = current_segment->get_residue_vector().at(rescounter);
	//  TR << " Current residue is " << rescounter << std::endl;
	//  for(core::Size atomcounter = 1; atomcounter <= current_residue->get_atom_vector().size(); atomcounter++){
	//   TR << "Atom " << atomcounter << " of residue " << rescounter << " of segment " << segcounter << " has coordinates " << current_residue->get_atom(atomcounter).xyz().at(0) << " " << current_residue->get_atom(atomcounter).xyz().at(1) << " " << current_residue->get_atom(atomcounter).xyz().at(2) << std::endl;
	//  }
	// }
	//}

	segment_list->set_version( version );
	return segment_list;
}

} //protocols
} //sewing
} //hashing






