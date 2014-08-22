# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   basic/options/options_class.py
## @brief  Program options generation script classes
## @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)


import sys
import time

KnownTypes=['Boolean', 'Integer', 'Real', 'String', 'File', 'Path', 'BooleanVector', 'IntegerVector', 'ResidueChainVector', 'RealVector', 'StringVector', 'FileVector', 'PathVector']

# argument should be string or list<string>
#  - wrap strings with "".
def wrapCStrig(s):
    if s is not None:
        if type(s) == list:
            for i, e in enumerate(s): s[i] = wrapCStrig(e)
            return s

        if len(s) < 2 :
            return '"' + s + '"'
        else:
            if s[0] == '"' and s[-1] == '"': return s
            else: return '"' + s + '"'
    return s


class Option:
    def __init__(self, name=None, ctype=None, group=None, desc="No description", short="",
                 oldName="-",
                 lower=None, upper=None, default=None, legal=None, n=None, n_lower=None, n_upper=None, restrict_access=False):

        if ctype not in KnownTypes:
            print 'Unknown type:%s!!!' % ctype
            sys.exit()
        if default == 'none' or default == 'None':
            print '*** Option %(name)s will default to the *string* "%(default)s"' % vars()
            print '    If you want no default, write default=None (no quotes)'

        self.ctype = ctype;  self.name = name;    self.group = group
        self.desc = desc;    self.short = short;  self.oldName = oldName
        self.lower = lower;  self.upper = upper;  self.default = default;  self.legal=legal
        self.restrict_access = restrict_access;
        self.n = n;  self.n_lower = n_lower;  self.n_upper = n_upper
        #Is this the synonymous member of an option group (e.g. -in:file:file)
        self.is_option_group = False
        # Wraping c-strings in "
        if ctype == 'String' or ctype == 'Path' or ctype == 'File' or \
           ctype == 'StringVector' or ctype == 'PathVector' or ctype == 'FileVector':
            if ctype == 'StringVector':
              self.default = self.default
            else:
              self.default = wrapCStrig( self.default )
            self.lower = wrapCStrig( self.lower )
            self.upper = wrapCStrig( self.upper )
            if type( self.legal ) == type([]):
                for i in range(0, len(self.legal) ):
                    self.legal[i] = wrapCStrig( self.legal[i] )
            else: self.legal = wrapCStrig( self.legal )


    def get_namespace(self,level):
            if self.group:
                namespaces = self.group.split(':')
                #print namespaces
                if len(namespaces) > level:
                    return namespaces[level]

    # return C++ name of option object
    def getCName(self):
        if str( self.name[:1] ).isdigit() : return 'n' + self.name
        else: return self.name


    def getOptionCC(self):
        s = ''
        s += 'option.add( basic::options::OptionKeys::'
        if self.group:
            s += self.group.replace(':','::')+'::'
        s += self.getCName()+', "'+self.desc+'" )'
        if self.short:  s+= '.shortd( "' + self.short + '" )'
        if self.lower : s+= '.lower(' + self.lower + ')'
        if self.upper : s+= '.upper(' + self.upper + ')'
        if self.n: s+= '.n(' + self.n + ')'
        if self.n_lower: s+= '.n_lower(' + self.n_lower + ')'
        if self.n_upper: s+= '.n_upper(' + self.n_upper + ')'
        if self.legal :
            if type(self.legal) == type(''):
                s+= '.legal(' + self.legal + ')'
            else:
                for l in self.legal:
                    s+= '.legal(' + l + ')'
        if self.default is not None :
            if type(self.default) == type(''): s+= '.def(' + self.default + ')'
            else:
                if len( self.default ) > 0:
                    for d in self.default:
                        s+= '.def(' + d + ')'
                else:
                    s+='.def()'
        if self.restrict_access: s+= '.restrict_access(true)'
        #These have to be last, as they're in the base class
        if self.is_option_group: s+= '.is_group(true)'
        return s + ';\n'


    def getOptionKeysHH(self):
        s = '';  se = ''
        if self.group:
            for ns in self.group.split(':'):
                s += 'namespace ' + ns + ' { '
                se += ' }'
        s += 'extern '+ self.ctype+'OptionKey const '+self.getCName()+';'
        return s + se + '\n'


    def getOptionKeysCC(self):
        s = '';  se = ''
        if self.group:
            for ns in self.group.split(':'):
                s += 'namespace ' + ns + ' { '
                se += ' }'

        s += self.ctype+'OptionKey const '+self.getCName()+'( "'

        if self.group:
            s += self.group
            if self.name != self.group.split(':')[-1] : s += ':' + self.name
        else: s += self.name
        s+= '" ); ' + se + '\n'
        return s


    def getWikiTableRow(self):
        def smStr(s): return s or ''

        s =  ' |-\n'
        s += ' | -%(name)s <%(ctype)s>\n' % {'name':self.name, 'ctype':self.ctype}
        s += ' | ' + self.desc + '\n'
        s += ' | ' + smStr(self.lower) + '-' + smStr(self.upper) + '\n'
        if self.legal=='true' and self.default=='true': s += ' |\n'
        else:
            if type(self.default) == type( [] ):
                s += ' | ' + str( self.default ) + '\n'
            else: s += ' | ' + smStr(self.default) + '\n'
        s += ' | -' + self.oldName + '\n'
        s += ' |-\n'
        return s


    def getDoxygenRow(self):
        def smStr(s): return s or ''

        s = "<dl>\n"
        s += '<dt><b>-%(name)s</b> \\<%(ctype)s\\><dt>\n' % {'name':self.name, 'ctype':self.ctype}
        s += '<dd>' + self.desc + '</dd><br>\n'
        if self.lower or self.upper:
            s += '<dd>Range: ' + smStr(self.lower) + '-' + smStr(self.upper) + '</dd><br>\n'
        if self.legal=='true' and self.default=='true': pass #s += ' |\n'
        else:
            if type(self.default) == type( [] ): df = str( self.default )
            else: df = smStr(self.default)

            if df: s += '<dd>Default: ' + df + '</dd><br>\n'

        #s += ' | -' + self.oldName + '\n'
        s += '</dl>\n'
        return s


    def getMarkdownRow(self):
        #Unfortunately, pure markdown doesn't have a defintion list element - steal HTML
        def smStr(s): return s or ''

        s =  '<dt><b>-%(name)s</b> \\<%(ctype)s\\></dt>\n' % {'name':self.name, 'ctype':self.ctype}
        s += '<dd>' + self.desc + '<br/>'
        if self.lower or self.upper:
            s += 'Range: ' + smStr(self.lower) + '-' + smStr(self.upper) + '<br/>'
        if self.legal=='true' and self.default=='true': pass #s += ' |\n'
        else:
            if type(self.default) == type( [] ): df = str( self.default )
            else: df = smStr(self.default)

            if df: s += 'Default: ' + df + '<br/>'

        s += '</dd>\n'
        return s



def Option_Group(group, *args):
    res = []

    for o in args:  # first concat all lists
        if type(o) == type([]): res += o
        else: res.append( o )

    if group:
        found_option_for_group = False
        for o in res:
            if o.group: o.group = group + ':' + o.group
            else: o.group = group
            if o.name == group:
                found_option_for_group = True
                o.is_option_group = True
        # In order to nest option specifiers in a flags file (passed with @flags.txt),
        # there must be a boolean option with the same name as each option group.
        if not found_option_for_group:
            o = Option( group, 'Boolean', group=group, desc=group+" option group", legal='true', default='true' )
            o.is_option_group = True
            res.insert( 0, o )
    return res


def writeToFile(opt_list, fileName, mapFunction):
    l = map(mapFunction, opt_list)
    f = file(fileName, 'wb');  f.write( "".join(l) );  f.close()

def printWikiTable(opt_list):
    s = ""
    prevGroup = None
    for o in opt_list:
        if prevGroup != o.group:  # Generating new table
            if prevGroup : s += ' |}\n' # Closing previos table if any
            s += """{| border="1" cellpadding="10" width="100%"\n |+ '''""" + (o.group or '')
            s += " Option Group'''\n"
            s += ' ! Option name\n'                    # width="10%" |
            s += ' ! Description\n'
            s += ' ! Range\n'
            s += ' ! Default\n'
            s += ' ! Old name\n'
            s += ' |-\n'
        s += o.getWikiTableRow()
        prevGroup = o.group
    return s + ' |}\n'


def getDoxygenPage(opt_list):
    s = "/*!\n@page full_options_list Rosetta command line option descriptions.\n"
    s += "<i>(This is an automatically generated file, do not edit!)</i> Generated: "+time.strftime("%Y-%m-%d")+"\n"
    s += "<ul>\n"
    prevGroup = None
    for o in opt_list:
        if prevGroup != o.group:  # Generating new table
            if prevGroup : s += ' </li>\n' # Closing previos table if any
            s += "<li><h2>" + (o.group or '') + "</h2>\n"
        s += o.getDoxygenRow()
        prevGroup = o.group
    return s + "</ul>\n */\n"


def getMarkdownPage(opt_list):
    s =  "# List of Rosetta command line options.\n\n"
    s += "_(This is an automatically generated file, do not edit!)_\n"
    s += "Generated: " + time.strftime("%Y-%m-%d") + "\n\n"
    s += "_Note that some application specific options may not be present in this list._\n\n"
    s += "[[_TOC_]]\n"
    in_dl = False
    prevGroup = []
    for o in opt_list:
        if prevGroup != o.group:  # Generating new group
            hlevel=min(6, len(o.group.split(":"))+1)
            if in_dl:
                s += "</dl>\n"
            s += "+ <h"+str(hlevel)+">-" + (o.group or '') + "</h"+str(hlevel)+">\n"
            s += "<dl>\n"
            in_dl = True
        s += o.getMarkdownRow()
        prevGroup = o.group
    if in_dl:
        s += "</dl>"
    return s + "\n"
