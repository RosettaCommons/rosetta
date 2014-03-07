// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)
///
// Unit Headers
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <core/types.hh>

// Boost Headers
#include <boost/foreach.hpp>



namespace protocols {
namespace ligand_docking {
namespace ligand_options {

Interface::Interface(
		core::Size num,
		InterfaceInfo info
)	: utility::vector1<InterfaceInfo>(num, info)
{}

void Interface::enforce_minimum_length(
		core::Size const chain_begin,
		core::Size const chain_end,
		core::Size const window
){
	for(core::Size i=chain_begin; i<=chain_end; ++i){
		if( (*this)[i].type == InterfaceInfo::is_interface){
			for(core::Size j = 1; j <= window; ++j){
				if(
						i>j	// otherwise the next line fails since we use unsigned ints
						&& i-j >= chain_begin // bounds check
						&& (*this)[i-j].type == InterfaceInfo::non_interface
				){
					(*this)[i-j].type = InterfaceInfo::near_interface;
					(*this)[i-j].chain=(*this)[i].chain;
					(*this)[i-j].chain_id=(*this)[i].chain_id;
				}
				if(
						i+j <= chain_end //bounds check
						&& (*this)[i+j].type == InterfaceInfo::non_interface
				){
					(*this)[i+j].type = InterfaceInfo::near_interface;
					(*this)[i+j].chain=(*this)[i].chain;
					(*this)[i+j].chain_id=(*this)[i].chain_id;
				}
			}
		}
	}
}

core::Size Interface::find_first_interface_residue(core::Size chain_begin, core::Size const chain_end) const{
	return (*this)[chain_begin].type == InterfaceInfo::non_interface?
			find_start_of_next_interface_region(chain_begin, chain_end): chain_begin;
}

core::Size Interface::find_start_of_next_interface_region(
		core::Size start_from,
		core::Size const chain_end
)const{
	assert( (*this)[start_from].type == InterfaceInfo::non_interface);

	for(++start_from; start_from <= chain_end; ++start_from){
		if((*this)[start_from].type != InterfaceInfo::non_interface) {
			return start_from;
		}
	}
	return 0;
}

core::Size Interface::find_stop_of_this_interface_region(
		core::Size start_from,
		core::Size const chain_end
)const{
	assert( (*this)[start_from].type != InterfaceInfo::non_interface );
	for(++start_from ;start_from<= chain_end; ++start_from){
		if((*this)[start_from].type == InterfaceInfo::non_interface){
			return start_from-1;
		}
	}
	return chain_end;

}

utility::vector1<core::Size>
Interface::get_interface_residues() const{
	utility::vector1<core::Size> interface_residues;

	for(core::Size position=1; position<= this->size(); ++position){
		if(this->at(position).type == InterfaceInfo::is_interface){
			interface_residues.push_back(position);
		}
	}
	return interface_residues;
}

utility::vector1<core::Size>
Interface::get_near_interface_residues() const{
	utility::vector1<core::Size> interface_residues;

	for(core::Size position=1; position<= this->size(); ++position){
		if(this->at(position).type == InterfaceInfo::near_interface){
			interface_residues.push_back(position);
		}
	}
	return interface_residues;
}

std::string Interface::get_python_string() const{
	std::stringstream python_stream;

	python_stream<<"interface residues: ";
	BOOST_FOREACH(core::Size res_id, get_interface_residues()){
		python_stream<< res_id << '+';
	}
	python_stream<< std::endl;
	python_stream<<"near interface residues: ";

	utility::vector1<core::Size> near_residues= get_near_interface_residues();
	utility::vector1<core::Size>::const_iterator j= near_residues.begin();
	for(; j != near_residues.end(); ++j){
		python_stream<< *j << '+';
	}
	python_stream<< std::endl;



	return python_stream.str();
}

std::ostream & operator<<(std::ostream& output, Interface const & interface){
	output<<interface.get_python_string();
	return output;
}

} //namespace ligand_options
} //namespace ligand_docking
} //namespace protocols
