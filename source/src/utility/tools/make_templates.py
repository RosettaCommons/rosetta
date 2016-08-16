#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#
## @file   utility/tools/make_templates.py
## @brief  Program that generate make_vector.hh and make_map.hh files
## @author Sergey Lyskov

NFunctions = 25  # number of functional argument that will be genrated

Header = """// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file   utility/tools/make_vector.hh
/// @brief  Common function to build vector, vector0, vector1, map.
/// @author Sergey Lyskov

#ifndef INCLUDED_utility_tools_%(iguard)s
#define INCLUDED_utility_tools_%(iguard)s

%(include)s

namespace utility {
namespace tools {

"""

Footer = """
} // namespace utility
} // namespace tools

#endif // INCLUDED_utility_tools_%(iguard)s

"""

def getVectorFunction(n, className, functionName):
    r = "template<typename T>\n%s<T> %s(const T & i0" % (className, functionName)
    for i in range(1, n): r += ', const T & i%d' % i
    r += ')\n{\n  %s<T> v;\n ' % className
    for i in range(n): r += ' v.push_back(i%d);' % i
    r += '\n  return(v);\n}\n\n'
    return r

def getMapFunction(n, className, functionName):
    r = "template<typename T1, typename T2>\n%s<T1, T2> make_map(const T1 &f0, const T2 & s0" % className
    for i in range(1, n): r += ', const T1 & f%d, const T2 & s%d' % (i,i)
    r += ')\n{\n  %s<T1, T2> m;\n '% className
    for i in range(n): r += ' m[f%d]=s%d;' % (i, i)
    r += '\n  return(m);\n}\n\n'
    return r


def makeFile(fileName, include, function, className, cppName):
    iguard = fileName.replace('.', '_')
    s = Header % dict(include=include, iguard=iguard)
    for i in range(1, NFunctions):
        s += function(i, className, cppName)
    s += Footer % dict(iguard=iguard)
    f = file(fileName, 'wb');  f.write(s);  f.close()


tasks = [
    dict(fileName='make_vector.hh', include='#include <vector>\n',
         function=getVectorFunction, className='std::vector', cppName='make_vector'),

    dict(fileName='make_vector0.hh', include='#include <utility/vector0.fwd.hh>\n',
         function=getVectorFunction, className='utility::vector0', cppName='make_vector0'),

    dict(fileName='make_vector1.hh', include='#include <utility/vector1.fwd.hh>\n',
         function=getVectorFunction, className='utility::vector1', cppName='make_vector1'),

    dict(fileName='make_map.hh', include='#include <map>\n',
         function=getMapFunction, className='std::map', cppName='make_map'),
]


def main():
    map(lambda x: makeFile(**x), tasks)

    """s = Header
    for i in range(1, NFunctions):
        s += getVectorFunction(i, 'std::vector',      'make_vector')
        s += getVectorFunction(i, 'utility::vector0', 'make_vector0')
        s += getVectorFunction(i, 'utility::vector1', 'make_vector1')
        s += getMapFunction(i, 'std::map')
    s += Footer

    f = file('make_vector.hh', 'wb');  f.write(s);  f.close()
    """

if __name__ == "__main__": main()
