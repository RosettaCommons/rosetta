/*
 * FastRelaxCreator.hh
 *
 *  Created on: Jun 21, 2010
 *      Author: zellnehj
 */


#ifndef INCLUDED_protocols_relax_FastRelaxCreator_hh
#define INCLUDED_protocols_relax_FastRelaxCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace relax {

class FastRelaxCreator : public moves::MoverCreator
{
public:
        virtual moves::MoverOP create_mover() const;
        virtual std::string keyname() const;
        static  std::string mover_name();
};

}
}


#endif
