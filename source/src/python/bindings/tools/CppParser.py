#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
# This code made available under dual license: RosettaCommons license and GPLv3

## @file   CppParser.py
## @brief  Parse C++ code and wrap it in Python
## @author Sergey Lyskov


import os, sys, gc, glob

import xml.dom.minidom

import doxygen

class CD:
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = '|'
        for i in dir(self):
            if not i.startswith('__'): r += '%s --> %s, ' % (i, getattr(self, i))
        return r[:-2]+'|'



class ReferenceSection:
    def __init__(self, ReferenceSection):
        # do we even need this??? self.ReferenceSection = ReferenceSection
        self.TypeNodes,    self.Types    = {}, {}
        self.ContextNodes, self.Contexts = {}, {}
        self.Files = {}

        for i in ReferenceSection.childNodes:
            if i.nodeName == 'Types':
                for j in i.childNodes:
                    if j.nodeName != '#text':  self.TypeNodes[j.getAttribute('id')] = j

            if i.nodeName == 'Contexts':
                for j in i.childNodes:
                    if j.nodeName != '#text':  self.ContextNodes[j.getAttribute('id')] = j

            if i.nodeName == 'Files':
                for j in i.childNodes:
                    if j.nodeName != '#text':  self.Files[j.getAttribute('id')] = j.getAttribute('name')


    def getContext(self, context_id):
        if self.Contexts.get(context_id) is None:
            node = self.ContextNodes[context_id]
            if node.nodeName == 'TranslationUnit': self.Contexts[context_id] = '::'
            else: self.Contexts[context_id] = self.getContext(node.getAttribute('context') ) + self.ContextNodes[context_id].getAttribute('name') +'::'

        #print 'getContext', context_id, '-->', self.Contexts[context_id]
        return self.Contexts[context_id]


    def getType(self, type_id):
        if self.Types.get(type_id) is None:
            node = self.TypeNodes[type_id]

            if node.nodeName == 'FundamentalType':  self.Types[type_id] = node.getAttribute('kind')
            else: self.Types[type_id] = 'Unknow type: ' + node.nodeName


        return self.Types[type_id]


    def getFile(self, file_id):
        pass



def normalize_context(context):  # making sure that context have form ::namespace(s)::
    if not context.startswith('::'): context = '::'+ context
    if not context.endswith('::'): context = context+'::'
    return context

#
# T Y P E S
#
class CppType(object):
    def getPureTypeName(self): return self.T()  # return name of a type without const or *, & modifiers...
    def getFile(self, simplify=False): return self.file_


#    def __repr__(self):  return  self.type_

class CppType_Fundamental(CppType):  # 'FundamentalType'
    def __init__(self, type_, kind): self.type_ = type_;  self.kind = kind;  self.file_ = None
    def getContext(self): return ''
    def getKind(self): return 'Fundamental'
    def T(self, simplify=False): return self.type_

class CppType_Simple(CppType):  # 'PointerType', 'ReferenceType', 'CvQualifiedType'
    def __init__(self, type_, postfix, kind): self.type_ = type_;  self.postfix = postfix;  self.kind = kind
    def getContext(self): return self.type_.getContext()
    def getKind(self): return self.type_.getKind()
    def getPureTypeName(self): return self.type_.getPureTypeName()
    def T(self, simplify=False): return self.type_.T(simplify=simplify) + ' ' + self.postfix
    def getFile(self): return self.type_.getFile()
    #def __repr__(self):  return  'CppType_Simple_'+self.postfix

class CppType_Composite(CppType):  #, 'FunctionType', 'ArrayType' - have some additionla information and non trivial composition
    def __init__(self, type_, kind): self.type_ = type_;  self.kind = kind
    def getContext(self): return ''
    def getKind(self): return 'Composite'
    def T(self, simplify=False): return '___XQWERTY___' + self.type_




def resolve_context_and_name(context, name):  # see if we can map self.context and self.name to more general type
    d = { ('::std::', '_Ios_Openmode') : ('::std::ios_base::', 'openmode'),
          ('::std::', '_Rb_tree_const_iterator<unsigned int>')     : ('::std::', 'set<unsigned int>::const_iterator'),
          ('::std::', '_Rb_tree_const_iterator<core::id::DOF_ID>') : ('::std::', 'set<core::id::DOF_ID>::const_iterator'),
          ('::std::', '_Rb_tree_const_iterator<std::pair<const double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > > >') : ('::std::', 'multimap< double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > >::const_iterator '),
          ('::utility::', 'vector1<std::_Rb_tree_const_iterator<protocols::match::downstream_hit>,std::allocator<std::_Rb_tree_const_iterator<protocols::match::downstream_hit> > >') : ('::utility::', 'vector1< std::set< protocols::match::downstream_hit >::const_iterator >'),
          ('::core::fragment::picking_old::concepts::', 'Extent<__gnu_cxx::__normal_iterator<const core::fragment::picking_old::vall::VallResidue*, std::vector<core::fragment::picking_old::vall::VallResidue, std::allocator<core::fragment::picking_old::vall::VallResidue> > > >') : ('::core::fragment::picking_old::concepts::', 'Extent< core::fragment::picking_old::vall::VallSection::PageConstIterator >'),

          ('::std::', '_Rb_tree_iterator<std::pair<const double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > > >') : ('::std::', 'multimap< double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > >::iterator '),
    }

    #while context.startswith('::'): context = context[2:]
    #while context.endswith('::'): context = context[:-2]
    #if not context.startswith('::'): context = '::'+ context
    #if not context.endswith('::'): context = context+'::'

    # TODO: Port some procedural mapping from T() below
    context_, name_ = d.get( (context, name), (context, name) )

    #if name.startswith('Extent<__gnu_cxx::'): print '~~ %s ~~ %s ~~' % (context, name), '--> ', (context_, name_)

    #return '::'+context_+'::', name_
    return context_, name_


class CppType_Complex(CppType):  # 'Class', 'Struct', 'Union', 'Typedef', 'Enumeration' - these one have context and definition location
    #def __init__(self, type_, context): self.type_ = type_ ;  self.context = context
    def __init__(self, name, context, kind):
        context = normalize_context(context)

        assert context.startswith('::'), 'CppType_Complex type must have context that starts with "::"!, got: {0}'.format(context)
        assert context.endswith('::'), 'CppType_Complex type must have context that ends with "::"!, got: {0}'.format(context)

        assert not context.startswith('::::'), 'CppType_Complex type must have context that starts with "::"!, got: {0}'.format(context)
        assert not context.endswith('::::'), 'CppType_Complex type must have context that ends with "::"!, got: {0}'.format(context)

        self.name = name;  self.context = context;  self.kind = kind;

    def getContext(self): return self.context

    def getKind(self): return self.kind
    def T(self, simplify=False):
        if self.context == '::' : return '::' + self.name


        # WARNING: Code below is mostly deprecated, for adjusting name maps see resolve_context_and_name above

        # work around for some strange types generated by GCCXML
        #if self.context == '__gnu_cxx':
        # d = { '::std::_Ios_Openmode' : '::std::ios_base::openmode',

        #       #'__normal_iterator<const core::chemical::Adduct*,std::vector<core::chemical::Adduct, std::allocator<core::chemical::Adduct> > >' : 'utility::vector1< core::chemical::Adduct >::const_iterator',
        #       #'_List_const_iterator<core::chemical::AA>' : 'std::list< core::chemical::AA >::const_iterator',
        #       #'_Rb_tree_const_iterator<std::pair<const std::basic_string<char, std::char_traits<char>, std::allocator<char> >, utility::pointer::access_ptr<const core::chemical::ResidueType> > >' : 'core::chemical::ResidueTypeSet::const_residue_iterator',

        #       #'__normal_iterator<const utility::pointer::owning_ptr<numeric::kdtree::KDPoint>*,std::vector<utility::pointer::owning_ptr<numeric::kdtree::KDPoint>, std::allocator<utility::pointer::owning_ptr<numeric::kdtree::KDPoint> > > >' : 'numeric::kdtree::KDPointList::const_iterator',
        #       #'__normal_iterator<utility::pointer::owning_ptr<numeric::kdtree::KDPoint>*,std::vector<utility::pointer::owning_ptr<numeric::kdtree::KDPoint>, std::allocator<utility::pointer::owning_ptr<numeric::kdtree::KDPoint> > > >' : 'numeric::kdtree::KDPointList::iterator',

        #       #'__normal_iterator<const core::conformation::PseudoBond*,std::vector<core::conformation::PseudoBond, std::allocator<core::conformation::PseudoBond> > >' : 'core::conformation::PseudoBondCollection::PBIter',

        #       '::core::fragment::picking_old::concepts::Extent<__gnu_cxx::__normal_iterator<const core::fragment::picking_old::vall::VallResidue*, std::vector<core::fragment::picking_old::vall::VallResidue, std::allocator<core::fragment::picking_old::vall::VallResidue> > > >' : 'core::fragment::picking_old::concepts::Extent< core::fragment::picking_old::vall::VallSection::PageConstIterator >',

        #       '::std::_Rb_tree_const_iterator<unsigned int>' : '::std::set<unsigned int>::const_iterator',
        #       '::utility::vector1<std::_Rb_tree_const_iterator<protocols::match::downstream_hit>,std::allocator<std::_Rb_tree_const_iterator<protocols::match::downstream_hit> > >' : 'utility::vector1< std::set< protocols::match::downstream_hit >::const_iterator >',
        #       #'::utility::vector1<std::_Rb_tree_const_iterator<protocols::match::upstream_hit>' : 'utility::vector1< std::set< protocols::match::downstream_hit >::const_iterator >',

        #       '::std::_Rb_tree_const_iterator<std::pair<const double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > > >' : '::std::multimap< double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > >::const_iterator ',
        #       '::std::_Rb_tree_iterator<std::pair<const double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > > >' : '::std::multimap< double, std::pair<unsigned int, std::pair<utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStubSet>, utility::pointer::owning_ptr<protocols::hotspot_hashing::HotspotStub> > > >::iterator ',

        # }
        # arg = '::'+self.context + '::' + self.name
        # r = d.get(arg, arg)
        context, name = resolve_context_and_name(self.context, self.name)
        r = context + name

        if   r == '::std::_Bit_reference': r = '::std::vector<bool>::reference'
        elif r.startswith('::std::_List_const_iterator<'): r = '::std::list< ' + r.partition('<')[2][:-1] + '> ::const_iterator'
        elif r.startswith('::std::_List_iterator<'): r = '::std::list< ' + r.partition('<')[2][:-1] + '> ::iterator'

        elif r.startswith('::std::_Rb_tree_const_iterator<std::pair<const '): r = '::std::map< ' + r.partition('<const ')[2][:-1] + '::const_iterator'
        elif r.startswith('::std::_Rb_tree_iterator<std::pair<const '): r = '::std::map< ' + r.partition('<const ')[2][:-1] + '::iterator'

        #elif r.startswith('::std::_Rb_tree_const_iterator<'): r = '::std::set< ' + r.partition('<')[2][:-1] + '::const_iterator'

        elif r.startswith('::__gnu_cxx::__normal_iterator<const '): r = r.partition('*,')[2][:-1] + ':: const_iterator'

        elif r.startswith('::__gnu_cxx::__normal_iterator<'):
            #print '_____________________', r
            #part_c =  r.partition(',')
            part_s =  r.partition('*,')
            if part_s[0].endswith('const'): r = part_s[2][:-1] + ':: const_iterator'
            else: r = part_s[2][:-1] + ':: iterator'

        #if r != arg: print '_____________ %s --> %s' % (arg, r)
        #if not r.startswith('::'): print'__________________', r, self.context, self.name; a = 1/0

        # normalizing type so it guarantee to have form ::context:: name
        #if not context.startswith('::'): context = '::'+ context
        #if not context.endswith('::'): context = context+'::'

        return r
        #return '::'+self.context + '::' + self.name


class CppType_Typedef(CppType_Complex):
    def __init__(self, name, context, type_):
        #self.name = name;  self.context = context;  self.type_ = type_; self.kind = 'Typedef'
        #CppType_Complex.__init__(self, name, context, kind='Typedef')  #
        super(CppType_Typedef, self).__init__(name, context, kind='Typedef')
        self.type_ = type_

    def T(self, simplify=False):
        if simplify: return self.type_.T(simplify=simplify)
        else: return CppType_Complex.T(self)



class CppNamespace:
    ''' data holder for C++ namespace. name represent namespace name and context is full context in which it recides
    '''
    def __init__(self, name, context):
        self.name, self.context = name, context  #normalize_context(context)

    def getChildrenContext(self):
        if self.context:  return self.context + self.name + '::'
        else:  return self.name

    def wrap(self): pass


class CppEnum:
    def __init__(self, node, context, refSection, name=None):
        self.name = name or node.getAttribute('name')
        self.context = context
        self.public  = (not node.hasAttribute('access')) or  (node.getAttribute('access') == 'public')
        self.values = []

        self.file_ =  refSection.Files[node.getAttribute('file')]
        self.line  =  node.getAttribute('line')

        for c in node.childNodes:
            if c.nodeName == 'EnumValue':  # GCC XML
                self.values.append( CD(name=c.getAttribute('name'), init=c.getAttribute('init') ) )


    def wrap_prefix_code(self, indent=''): return ''


    def wrap(self, indent='', wrappingScope='boost::python::'):
        print 'WrappingEnum %s...' % self.name

        if not self.public: return ''
        r = '  %senum_< %s%s >("%s", "%s")\n' % (wrappingScope, self.context, self.name, self.name, doxygen.getDoxygenComment(self.file_, self.line))
        for v in self.values:
            r += '    .value("%s", %s%s)\n' % (v.name, self.context, v.name)
        r += '    .export_values();\n'
        return  ('\n' + indent).join( r.split('\n') ) + '\n'


    def getFileList(self):  # return list of files where types was created (use it for includes!)
        return [self.file_]



class CppVariable:
    def __init__(self, node, context, refSection):
        self.name = node.getAttribute('name')
        self.context = context
        self.public  = (not node.hasAttribute('access')) or  (node.getAttribute('access') == 'public')
        self.type_ = refSection.getType( node.getAttribute('type') )

        self.file_ =  refSection.Files[node.getAttribute('file')]
        self.line  =  node.getAttribute('line')


    def wrap(self, indent='', wrappingScope='boost::python::'):
        if not self.public: return ''
        #print '________________ self.type_:', self.type_
        #if isinstance(self.type_, CppType_Fundamental): readwrite = 'readwrite'
        #else: readwrite = 'readonly'
        readwrite = 'readwrite'

        r = '  %sdef_%s("%s", &%s%s);\n' % (wrappingScope, readwrite, self.name, self.context, self.name)
        return  ('\n' + indent).join( r.split('\n') ) + '\n'



class CppFunction:
    OperatorToCppNameMap = {'=':'assign', '()':'__call__', '[]':'__getitem__',
        '*':'__mul__', '+':'__add__', '-':'__sub__', '/':'__div__', '%':'__mod__', '^':'__pow__',}
    def __init__(self, node, context, returnType, refSection, memberFunction=False, const=False, constructor=False, const_overloaded=False):
        self.name = node.getAttribute('name')
        self.context = context
        self.demangled = node.getAttribute('demangled')
        self.returnType = returnType
        self.memberFunction = memberFunction
        self.const = const
        self.const_overloaded = const_overloaded
        self.constructor = constructor
        self.reference = refSection

        self.operator      =  node.nodeName == 'OperatorMethod'  or  node.nodeName == 'OperatorFunction'
        self.static        =  node.hasAttribute('static')        # GCC XML
        self.virtual       =  node.hasAttribute('virtual')       # GCC XML
        self.pure_virtual  =  node.hasAttribute('pure_virtual')  # GCC XML
        self.public        =  (not node.hasAttribute('access')) or  (node.getAttribute('access') == 'public')    # to do - create 'access' datamember and
        self.protected     =  (not node.hasAttribute('access')) or  (node.getAttribute('access') == 'protected') #    turn this in to a propertys
        self.private       =  (not node.hasAttribute('access')) or  (node.getAttribute('access') == 'private')
        self.artificial    =  node.hasAttribute('artificial')
        self.file_         =  refSection.Files[node.getAttribute('file')]
        self.line          =  node.getAttribute('line')

        self.argsTypes = []
        for c in node.childNodes:
            if c.nodeName == 'ParmVar':  # Clang XML
                self.argsTypes.append( CD(name=c.getAttribute('name'), type_=refSection.getType( c.getAttribute('type') ) ) )

            if c.nodeName == 'Argument':  # GCC XML
                self.argsTypes.append( CD(name=c.getAttribute('name'), type_=refSection.getType( c.getAttribute('type') ) , default=c.getAttribute('default') ) )


    def __str__(self):
        ''' Generate string representing function name and types (but no context information for function itself), like this: int foo(char a, bool b)
        '''
        return ( '%s %s(%s)' % (self.returnType.T(simplify=True), self.name , self.getSimpleArgsType(simplify=True) ) )
        '''    .replace('::platform::Size', '::core::Size')
                        .replace('::core::scoring::rna::Vector', '::numeric::xyzVector<double>') \
                        .replace('::core::scoring::methods::EnergyMethodOP', '::utility::pointer::owning_ptr<core::scoring::methods::EnergyMethod>') \
                        .replace('::core::Vector', '::numeric::xyzVector<double>').replace('::core::kinematics::DomainMap', '::ObjexxFCL::FArray1D<int>') \
                        .replace('::std::basic_string<char,std::char_traits<char>,std::allocator<char> >', '::std::string') \
                        .replace('double', '::core::Real').replace('long unsigned int', '::core::Size')
                        '''

    def isDefalutArgWrappable(self, x):
        ''' check if default arg can be wrapped and return string for it, return epmty string castable to False otherwise
        '''
        _r = ''
        if x.default:
            # GCC gives non-working default arg here, so what we do is just bind without default... - probably could be fixed later
            # in future we can just locate enums in all parsed code but for now lets not guess...
            known_bad_defaults = ['(((utility::tag::Tag*)operator new(', 'utility::vector0<int, std::allocator<int>', 'MATCH_YR', 'std::cout', 'typename ', 'std::make_pair [with ',
                'core::chemical::ChemicalManager::residue_type_set(const std::string&)(((const std::string&)(& core::chemical::FA_STANDARD',
                'core::scoring::hbonds::DUMMY_DERIVS', 'core::fragment::BBTorsionSRFD', 'std::ios_base::in', 'protocols::forge::build::SegmentInsertConnectionScheme::RANDOM_SIDE',
              'protocols::stepwise::sampling::rna::local_count_data', 'numeric::xyzVector<T>',
            ]
            for bd in known_bad_defaults:
                if x.default.startswith(bd): return _r

            if len( x.default ) > 64: return _r  # Golden rule: if default arg description more then 64 - gcc got it wrong!

            #print '______', x.default, x.type_.getKind()

            if x.type_.getKind() == 'Enumeration':
                if x.default.startswith(x.type_.getContext())  or x.default.startswith(x.type_.getContext()[2:]) :  # starts with '::namespace' or just 'namespace'
                    _r += '=%s' % x.default
                else:
                    _r += '=%s%s' % (x.type_.getContext(), x.default)  # x.type_.getContext(),
            else:
                _r += '=(%s)' % (x.default)
        # Now, we need to find if all args to the right is also wrappable...
        for a in self.argsTypes[ self.argsTypes.index(x)+1 : ]:
            if not self.isDefalutArgWrappable(a): return ''
        return _r

    def writeFunctionDeclarationTypes(self, withArgs=True, withDefaultArgs=False):
        #def write_def(x):
        #    if self.isDefalutArgWrappable(x): return '=%s' % self.isDefalutArgWrappable(x)
        #    #if x.default: return '=%s' % x.default
        #    else: return ''
        return ', '.join( map(lambda (i, x): x.type_.T() + ' __a%s%s' % (i, self.isDefalutArgWrappable(x)), enumerate(self.argsTypes)) )


    def getSimpleArgsType(self, constructor=False, withArgs=True, simplify=False):
        if not constructor: return ', '.join( map(lambda (i, x): x.type_.T(simplify=simplify) + ' __a%s' % i, enumerate(self.argsTypes)) )
        else:
            def foo(x):
                if x.default and self.isDefalutArgWrappable(x):
                    return x.type_.T()
                else: return x.type_.T()
            r, default_state = '', False
            for x in self.argsTypes:
                if (not default_state) and  self.isDefalutArgWrappable(x):
                    default_state = True
                    r += 'boost::python::optional< %s ' % x.type_.T()
                else: r += x.type_.T()
                r += ', '

            if default_state: r = r[:-2] + ' >, '
            return r[:-2]


    def getArgTypes(self):
        if self.argsTypes:
            def wrap_argument(x):
                _r = 'boost::python::arg("%s")' % x.name
                _r += self.isDefalutArgWrappable(x)

                return _r

            #return '    , ( %s )\n' % ', '.join( map(wrap_argument, self.argsTypes) )
            return ', '.join( map(wrap_argument, self.argsTypes) )
        else: return ''


    def isWrapable(self):
        ''' Check if we can really wrap this object now, not theoretical but practically...
        '''
        #print '_____ %s number of args: %s' % (self.name, len(self.argsTypes) )
        if len(self.argsTypes) > 70:
            print '\033[31m\033[1m%s\033[0m' % (  'Too many arguments for function:%s... Skipping...' % self.demangled )
            return False

        tp = [self.returnType] + map(lambda x: x.type_, self.argsTypes)


        for a in tp:  # check if function contain types that we don't know how to deal with yet...
            #print self.name, a.T()
            if a.T().find('___XQWERTY___') >= 0: return False
            if a.T().find('::std::pair<boost::unordered_detail::hash_iterator_equivalent_keys<std::allocator<std::pair<') >= 0: return False

        #if self.pure_virtual or (not self.public): return False
        if (not self.public): return False

        #  check if return result is sane. Some types (like int *) mean that we actully have an iterator functions and this one should be taked care at class abstraction level
        for t in ['void *', 'void const *', 'char *', 'int *', 'int const *', 'double *', 'double const *', 'Real *', 'Real const *', 'Size *',
                  'kiss_fft_cpx *', 'kiss_fft_cfg', 'kiss_fftnd_cfg', 'kiss_fftr_cfg', ]:
            if (self.returnType.T() or '').endswith(t): return False

        if self.returnType.T().startswith('::boost::unordered_detail::hash_const_iterator_unique_keys<std::allocator') : return False

        #if self.operator  and  self.name in ['<', '>', '+', '-', '->', '*', '&', '++', '--', '+=', '-=', '*=', '==', '!=']: return False  # will deal with this later...

        if self.operator  and  self.name=='='  and   self.artificial:  return False
            #print self.context, self.name, self.returnType.T(), self.argsTypes[0].type_.T()
            #if self.returnType.T() != self.context[:-2] + ' &'  or self.argsTypes[0].type_.T() != self.context[:-2] + ' const &': return False
            #print self.name, '= C=%s, N=%s, T=%s' % (self.context, self.name, self.argsTypes[0].type_.T())

        if self.operator:
            if self.name not in self.OperatorToCppNameMap: return False
            if not self.memberFunction: return False

        return True


    def getReturnValuePolicyString(self):
        primitive_types = ['bool', 'char', 'int', 'float', 'double', 'size_t', 'string', 'Real', 'Size', 'Length', 'core::PackerEnergy', 'core::Energy']

        def endsWith(str_, suffix, list_=primitive_types):  # return true if str_ end with any element on list_
            for x in [ e+suffix for e in list_]:
                if str_.endswith(x+' &') or str_.endswith(x+' *'): return True
            return False

        if self.returnType.getKind() == 'Enumeration':
            if endsWith(self.returnType.T(), ' const', ['']):  # enum_type const &, enum_type const *,
                return '    , boost::python::return_value_policy< boost::python::copy_const_reference >()\n'

            elif self.returnType.T().endswith(' *') or self.returnType.T().endswith(' &'): # enum_type &, enum_type *,
                return '    , boost::python::return_value_policy< boost::python::copy_non_const_reference>()\n'

            return ''

        elif endsWith(self.returnType.T(), ' const'):  # int const &, int const *
            return '    , boost::python::return_value_policy< boost::python::copy_const_reference >()\n'

        elif endsWith(self.returnType.T(), ''):  # int &, int *
            return '    , boost::python::return_value_policy< boost::python::copy_non_const_reference>()\n'

        #elif endsWith(self.returnType.T(), ' const', ['']):  # const &, const *, ... &, ... *
        #    return '    , boost::python::return_value_policy< boost::python::reference_existing_object >()\n'

	# The boost::graph vertex descriptors and edge descriptors are implemented with void pointers
	# which gives boost::python a headache. Tell boost::python to treat them as opaque objects.
        elif    self.returnType.T().endswith('core::chemical::VD') \
             or self.returnType.T().endswith('core::chemical::VD const &') \
             or self.returnType.T().endswith('core::chemical::ED'):
            return '    , boost::python::return_value_policy< boost::python::return_by_value >()\n'


        elif self.returnType.T().endswith(' *') or self.returnType.T().endswith(' &'):  # ... &, ... *
            return '    , boost::python::return_value_policy< boost::python::reference_existing_object >()\n'



        return ''


    def isCallbackStructNeeded(self):
        #if self.name == '=': return False  # virtual assign operator? No way we can wrap it...
        if self.virtual and self.public: return True  # and (not f.const)
        return False

    def writeCallbackCode(self, class_, callback, indent=''):
        if not self.isWrapable(): return ''
        r = ''
        #if (self.virtual and self.public and not self.const_overloaded) or  (self.constructor and self.public):  # for now we allow to overload only non-const functions in Python...
        if self.isCallbackStructNeeded() or  (self.constructor and self.public):  # for now we allow to overload only non-const functions in Python...
            def foo( (i, a) ):
                #if (not self.constructor) and a.type_.T()[-1] == '&'  and False: return 'utility::PyAP( __a%s )' % i
                #if (not self.constructor) and a.type_.T()[-1] == '&': return 'utility::py::to_py_object( __a%s )' % i
                #if (not self.constructor) and a.type_.T()[-1] == '&' and False: return 'boost::python::ptr( & __a%s )' % i
                # and  not a.type_.T().endswith(' const &') \
                if (not self.constructor) and a.type_.T()[-1] == '&' \
                   and a.type_.T() not in 'bool &  ::platform::Size &  ::platform::Size const &  ::numeric::Real &  ::core::Size &  ::core::Size const &  ::core::Real &  ::core::chemical::AA const &  ::numeric::Size &':
                    return 'boost::python::ptr( & __a%s )' % i
                else: return '__a%s' % i

            args = ', '.join( map(lambda (i, x): '__a%s' % i, enumerate(self.argsTypes)) )
            args_AP = ', '.join( map(foo, enumerate(self.argsTypes)) )
            if self.constructor: r += '\n%s(%s) : %s(%s) {' % (callback, self.writeFunctionDeclarationTypes(withDefaultArgs=True), class_, args)  # Writing constructors...  # and not self.artificial
            else:
                if self.operator:
                    py_name =  self.OperatorToCppNameMap [self.name]
                    cpp_name = 'operator' + self.name
                else:
                    cpp_name = self.name
                    py_name = self.name

                r += '\n%s %s(%s) %s{\n' % (self.returnType.T(), cpp_name, self.writeFunctionDeclarationTypes(withDefaultArgs=True), ['', ' const '][self.const])
                return_ = ['return', ''][self.returnType.T() == 'void']
                if self.pure_virtual:
                    r += '  #ifndef _MSC_VER\n'
                    r += '    %s this->get_override("%s")(%s);\n' % (return_, py_name, args_AP)
                    r += '  #else\n'
                    if args_AP:
                        r += '      return boost::python::call< %s >(this->get_override("%s").ptr(), %s );\n' % (self.returnType.T(), py_name, args_AP)
                    else:
                        r += '      return boost::python::call< %s >(this->get_override("%s").ptr() );\n' % (self.returnType.T(), py_name)
                    r += '  #endif\n'
                else:
                    r += '  if( boost::python::override f = this->get_override("%s") ) {\n' % py_name
                    r += '    #ifndef _MSC_VER\n'
                    r += '      %s f(%s);\n' % (return_, args_AP)
                    r += '    #else\n'
                    if args_AP:
                        r += '      return boost::python::call< %s >(f.ptr(), %s );\n' % (self.returnType.T(), args_AP)  # Should we double check for void return type here???
                    else:
                        r += '      return boost::python::call< %s >(f.ptr());\n' % (self.returnType.T(), )  # Should we double check for void return type here???
                    r += '    #endif\n'
                    r += '  } else %s %s(%s);\n' % (return_, self.context+cpp_name, args)
                    r += '}\n'
                    r += '\n%s __default_%s(%s) %s{\n  %s this->%s(%s);\n' % (self.returnType.T(), py_name, self.getSimpleArgsType(), ['', ' const '][self.const], return_, self.context+cpp_name, args)

            r += '}\n'
            r = ('\n' + indent).join( r.split('\n') ) + '\n'
        return r


    def wrap_prefix_code(self, indent=''): return ''

    def wrap(self, indent='', wrappingScope='boost::python::', useCallbackStruct=''):#overloaded=False):
        if not self.isWrapable(): return ''
        args = self.getSimpleArgsType()

        if self.operator:
            py_name = self.OperatorToCppNameMap [self.name]
            cpp_name = 'operator' + self.name
        else:
            cpp_name = self.name
            py_name = self.name

        if useCallbackStruct:  useCallbackStruct += '::'

        D = dict(py_name=py_name, cpp_name=cpp_name, context=self.context, returnType=self.returnType.T(), args=args, CallbackStruct=useCallbackStruct,
                 wrappingScope=wrappingScope, const=' const' if self.const else '', doc=doxygen.getDoxygenComment(self.file_, self.line), file_=self.file_, line=self.line)
        r = '\n{ // %(context)s%(cpp_name)s%(const)s file:%(file_)s line:%(line)s\n' % D

        if self.memberFunction and not self.static:
            r += '  typedef %(returnType)s ( %(context)s * %(py_name)s_function_type)(%(args)s) %(const)s;\n\n' % D
            if useCallbackStruct: r += '  typedef %(returnType)s ( %(CallbackStruct)s * %(py_name)s_callback_function_type)(%(args)s) %(const)s;\n\n' % D
        else:
            r += '  typedef %(returnType)s ( * %(py_name)s_function_type)(%(args)s);\n\n' % D

        r += '  %(wrappingScope)sdef("%(py_name)s"\n' % D

        if self.pure_virtual:  r += '    , ::boost::python::pure_virtual( %(py_name)s_function_type( &%(context)s%(cpp_name)s ) )\n' % D
        else:                  r += '    , %(py_name)s_function_type( &%(context)s%(cpp_name)s )\n' % D  # r += '    , %s_function_type( &%s%s )\n' % (self.name, self.context, self.name )

        #if useCallbackStruct  and  self.virtual  and  not self.pure_virtual  and not self.const_overloaded:
        if useCallbackStruct and self.isCallbackStructNeeded() and not self.pure_virtual  and not self.const_overloaded:
            r += '    , %(py_name)s_callback_function_type( &%(CallbackStruct)s__default_%(py_name)s )\n' % D
            #r += '    , &%(CallbackStruct)s::__default_%(py_name)s\n' % D

        if self.argsTypes:  r+= '    , ( %s )\n' % self.getArgTypes()

        #print self.returnType.T(), self.name, '-->', self.getReturnValuePolicyString()
        r += self.getReturnValuePolicyString()

        # Bug in boost: it can't hanlde more then 5 Template args here, so we omit doc string if this is the case...
        if (useCallbackStruct  and  self.virtual  and  not self.pure_virtual  and not self.const_overloaded) and \
           self.argsTypes and self.getReturnValuePolicyString(): r += '    );\n' % D
        else: r += '    , "%(doc)s" );\n' % D

        #if self.static: r += '  %(wrappingScope)sstaticmethod("%(py_name)s");\n' % D
        r += '}\n'

        return  ('\n' + indent).join( r.split('\n') ) + '\n'

    def getFileList(self): # , exclude=[] # return list of files where types was created (use it for includes!)
        '''if self in exclude:
            return []

        exclude.append(self)
        '''
        r = [self.file_]
        for a in [self.returnType] + map(lambda x: x.type_, self.argsTypes):
            r.append(a.getFile())
            '''
            base_type = a.T().split(' ')[0]
            if base_type in self.reference.Objects:
                #print '~~~~~ ', base_type
                self.reference.Objects[ base_type ].getFileList(exclude)
                '''

        return r

    def __repr__(self):
        return 'CppFunction: %s Context=%s isWrapable=%s, %s %s' % (self.name, self.context, self.isWrapable(), self.returnType, self.argsTypes)



class CppClass:
    def __init__(self, node, context, refSection):
        self.reference = refSection

        #print 'CppClass________:', context, node.getAttribute('name')

        #self.name = node.getAttribute('name')
        #self.context = context
        self.context, self.name = resolve_context_and_name(context, node.getAttribute('name'))

        self.bases = []
        self.constructors = []
        self.destructor = None
        self.functions = []
        self.enums = []
        self.dataMembers = []
        self.requiredTypes = []  # types that need to be declared for class to be parsable (ie includes needed)

        self.abstract   = node.hasAttribute('abstract')
        self.incomplete = node.hasAttribute('incomplete')
        self.file_      = refSection.Files[node.getAttribute('file')]
        self.line       = node.getAttribute('line')

        # work-around GCC bug when templates put in to wrong files...
        #if self.name == 'UpperEdgeGraph<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>':
        #    print '~~~~~~~~~~ ', self.file_
        #    self.file_ = 'core/conformation/PointGraph.hh'

        for ch in node.childNodes:
            #if ch.nodeName == '#text': continue  # why on earth this node type is everywhere???
            if ch.nodeName == 'Base':
                self.bases.append( CD( type_=refSection.getType(ch.getAttribute('type')), access=ch.getAttribute('access'), virtual=bool(int(ch.getAttribute('virtual')))  ) )

    def getVirtualBases(self):  # return list of virtual bases classes
        r = []
        for b in self.bases:
            if b.type_.T() in self.reference.Objects:
                o = self.reference.Objects[ b.type_.T() ]
                if b.virtual  and  o.context != '::std::': r.append( o.context + o.name )
                r += o.getVirtualBases()
        #print 'VirtualBases:', r
        return r


    def isBase(self, base):
        for b in self.bases:
            if b.type_.T() == base.context+base.name: return True
            if b.type_.T() in self.reference.Objects and  (self.reference.Objects[ b.type_.T() ].isBase(base)): return True
        return False


    #def _isDepend(self, other): # is SELF depend on OTHER ← too confusing, renaming...
    def isNeeded(self, other): # return True if other needed to construct self, and False otherwise
        if self.context+self.name == other.context+other.name: return True

        # check if other is part of our name as template specilisation, for example: Xforms... [ ::utility::vector1<numeric::xyzTransform<double>,std::allocator<numeric::xyzTransform<double> > > ]
        left, _, right = (self.context+self.name).partition( (other.context+other.name)[2:] )
        if left and right:
            if left[-1] in ' ,<'  and  right[0] in ' ,>': return True

        for b in self.bases:
            if b.type_.T() == (other.context+other.name)[2:] : return True
            #if b.type_.T().find( (other.context+other.name)[2:] ) >= 0 :
            #    left, _, right = b.type_.T().partition( (other.context+other.name)[2:] )
            #    if left and right: return True
            if b.type_.T() in self.reference.Objects:
                o = self.reference.Objects[ b.type_.T() ]
                if o.context+o.name == other.context+other.name: return True
            if b.type_.T() in self.reference.Objects and  (self.reference.Objects[ b.type_.T() ].isNeeded(other)): return True
        return False

    # def isNeeded(self, other):
    #     res = self._isNeeded(other)
    #     check_listA, check_listB = ['RotamerSets', 'FixbbRotamerSets', 'Xforms', 'xyzTransform<double>'], ['FixbbRotamerSets', 'RotamerSets']
    #     if self.name in check_listA  and  other.name in check_listB:
    #         print '____', self.name, other.name, res, other._isNeeded(self)
    #         for b in self.bases:
    #             # if b.type_.T().find( (other.context+other.name)[2:] ) >= 0 :
    #             #     left, _, right = b.type_.T().partition( (other.context+other.name)[2:] )
    #             #     if left and right: print 'b.type_.T():', b.type_.T(), '   other.context+other.name:', other.context+other.name
    #             #     #print 'left:{0}, _:{1}, right:{2}'.format(left, _, right)
    #             if b.type_.T() in self.reference.Objects and  (self.reference.Objects[ b.type_.T() ].isNeeded(other)):
    #                 print '__REF: {0}-{1}'.format(self.reference.Objects[ b.type_.T() ].name, other.name)
    #     return res


    def getChildrenContext(self):  return self.context + self.name + '::'

    def getAllVirtualFunctions(self, D, checkAll=True):
        ''' fill dict with all virtual functions that need to be implemented before class could be created. {full cpp name : CppFunction}
        - return non-empty dict if all function that needed (pure_virtual) are wrapable, and false ({}) otherwise
        '''
        for f in self.functions:
            if f.public and f.virtual  and (checkAll or f.pure_virtual):
                name =  str(f)

                if name not in D  and f.virtual:
                    D[name] = f
                    if f.pure_virtual and (not f.isWrapable()): return {}


            elif f.pure_virtual: return {}

        for b in self.bases:
            if b.type_.T() in self.reference.Objects:
                o = self.reference.Objects[ b.type_.T() ]
                #if o.isWrapable():
                if b.access == 'public' and  (not o.getAllVirtualFunctions(D)): return {}
                elif not o.getAllVirtualFunctions(D, checkAll=False): return {}

        return D



    def isCreatable(self, asBase=False, with_callback_struct=False):
        if with_callback_struct:
            return self.getAllVirtualFunctions({})
            '''
            for f in self.functions:
                if f.pure_virtual  and  (not f.isWrapable()): return False

            for b in self.bases:
                if b.type_.T() in self.reference.Objects:
                    o = self.reference.Objects[ b.type_.T() ]
                    if o.isWrapable() and not o.isCreatable():
                        #print o, ' --> False'
                        return False

            return True
            '''
        else:

            if self.abstract: return False
            for f in self.functions:
                if f.pure_virtual: return False

            if self.constructors:
                for c in self.constructors:
                    if len(c.argsTypes) == 1  and  c.argsTypes[0].type_.T().find(self.context+self.name+' ') >= 0: continue  # skip copy constructor
                    if c.public: return True
                    if asBase and c.protected: return True
                #print self.name, 'Is NOT creatable!'
                return False

            return True


    def isCopyable(self, asBase=False):
        ''' check if there is private constructor and assign operator
        '''
        if self.abstract: return False
        if (self.destructor and (not self.destructor.public)) or (not self.isCreatable(asBase)): return False
        for c in self.constructors:  # check if copy construtor is present and private...
            if (c.private)  and  len( c.argsTypes ) == 1  and  c.argsTypes[0].type_.T() == self.getChildrenContext()[:-2] + ' const &' :  return False

        for b in self.bases:  # lests check if all bases are copyable too...
            if b.type_.T() in self.reference.Objects and  (not self.reference.Objects[ b.type_.T() ].isCopyable(asBase=True)): return False

        for b in self.bases:
           if b.type_.T() == '::boost::noncopyable_::noncopyable': return False

        #for f in self.functions:  # check if assign opperator is private...
        #    if (not f.public)  and  f.name == '=': return False  # not exactly exact mucth, but good enought for now
        return True




    def isCallbackStructNeeded(self):
        ''' Return True if class have public virtual functions, so it make sense to allow callback from Python...
        '''
        #print 'isCallbackStructNeeded[0]', self.context+self.name
        if self.getVirtualBases(): return False  # we cant auto create virtual base classes

        #print 'isCallbackStructNeeded[1]', self.context+self.name
        for f in self.functions:
            if f.pure_virtual and not f.public: return False

        #print 'isCallbackStructNeeded[2]', self.context+self.name
        if self.constructors: # check if there is a way to construct subclass... (ie constructors that public or no constructors)
            for c in self.constructors:
                if len(c.argsTypes) == 1  and  c.argsTypes[0].type_.T().find(self.context+self.name+' ') >= 0: continue  # skip copy constructor
                if c.public: break
            else: return False

        #print 'isCallbackStructNeeded[3]', self.context+self.name
        if not self.isCreatable(with_callback_struct=True): return False  # there is no point to wrap it if we can't create it

        #print 'isCallbackStructNeeded[4]', self.context+self.name
        for f in self.functions:
            #if f.virtual and f.public: return True  # and (not f.const)
            if f.isCallbackStructNeeded() and f.isWrapable(): return True  # and (not f.const)

        return False




    def isWrapable(self):
        ''' Check if we can really wrap this object now, not theoretical but practically...
        '''
        if self.context+self.name == '::utility::pointer::ReferenceCount': return False
        if self.incomplete: return False
        #for b in self.bases:
        #    if b.type_.T() == '::boost::noncopyable_::noncopyable': return False

        # Some exceptions...
        if self.context+self.name == '::protocols::neighbor::Neighborhood': return False  # temporary disabling wrapping of protocols::neighbor::Neighborhood because we could not compile it on Linux GCC-4.1 due to function pointer in constructor


        # some classes to avoid wrapping (we will do it by hand in utility, numeric)
        skip_list = [
            '::utility::vector1<int,std::allocator<int> >',
            '::utility::vector1<long unsigned int,std::allocator<long unsigned int> >',
            '::utility::vector1<double,std::allocator<double> >',
            '::utility::vector1<std::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::allocator<std::basic_string<char, std::char_traits<char>, std::allocator<char> > > >',
            '::numeric::xyzMatrix<double>',
            '::boost::noncopyable_::noncopyable'
                     ]

        for a in skip_list:
            if (self.context+self.name) == a: return False

        return True


    def isPublicMembersWrappable(self):
        ''' Check if public data members are wrappable or not and return list of wrappable one...
        '''
        if self.context+self.name == '::core::scoring::etable::AtomPairEnergy': return self.dataMembers
        else: return []


    def isHeldTypeOP(self):
        for b in self.bases:
            if b.type_.T() == '::utility::pointer::ReferenceCount': return True
            if b.type_.T() in self.reference.Objects and  (self.reference.Objects[ b.type_.T() ].isHeldTypeOP()): return True
        return False


    def isSelfInherit(self):
        for b in self.bases:
            #if b.type_.T() == '::core::scoring::etable::BaseEtableEnergy<core::scoring::etable::EtableEnergy>': continue
            _case1 = '::core::scoring::etable::BaseEtableEnergy<%s>' % (self.context+self.name)[2:]
            # ^^ special case, template use it self for initiate, we can't list as base class in Python because we will never be able to create it
            if b.type_.T() == _case1 : return _case1;
        return ''


    def getHeldType(self, heldTypeBase=None):
        if not heldTypeBase: heldTypeBase = self.context+self.name

        r =  heldTypeBase  #self.context+self.name

        bases = ''
        for b in self.bases:
            #print b.type_.T()
            # Special cases...
            #if b.type_.T() == '::core::scoring::etable::BaseEtableEnergy<core::scoring::etable::EtableEnergy>': continue
            if b.type_.T() == self.isSelfInherit(): continue

            if b.type_.T() == '::utility::pointer::ReferenceCount': continue
            if b.type_.T() == '::boost::noncopyable_::noncopyable': continue
            if b.type_.T().startswith('::utility::vector'): continue
            if b.type_.T().startswith('::std::iterator'): continue
            if b.type_.T().startswith('::std::'): continue
            if b.type_.T().startswith('::boost::enable_shared_from_this<'): continue
            #print b.type_.T()

            bases += b.type_.T() + ', '

        if bases: r += ', boost::python::bases< %s >' % bases[:-2]

        #if self.isHeldTypeOP(): r+=  ', ::utility::pointer::shared_ptr< %s >' % heldTypeBase #(self.context, self.name)
        #r+=  ', ::utility::pointer::shared_ptr< %s >' % heldTypeBase  # now we want all objects to held in SP
        r+=  ', ::boost::shared_ptr< %s >' % heldTypeBase  # now we want all objects to held in SP
        return r


    def markConstOverloaded(self):
        ''' function.const_overloaded → True/False for all const member functions if the same non-const function present
        '''
        for f in self.functions:
            if f.const:
                for h in self.functions:
                    if f.demangled[:-6] == h.demangled:  # we don't need to compare return types, just type signature with removed 'const'
                        print 'MakringConstOverloaded %s' % f.demangled, '→ True'
                        f.const_overloaded = True; break
                else: f.const_overloaded = False

    def write_implicitly_convertible_code(self, use_callback_struct, D):
        r = ''
        if self.isHeldTypeOP():
            convert_to = ['%(context)s%(name)s const' % D]
            for b in self.bases:
                if b.type_.T() != '::utility::pointer::ReferenceCount' and \
                   b.type_.T() in  self.reference.Objects  and  self.reference.Objects[ b.type_.T() ].isHeldTypeOP():
                     convert_to.append( b.type_.T() )
            r += '\n'

            # Adding COP → OP conversion
            # r += '  boost::python::to_python_converter< utility::pointer::owning_ptr< %(context)s%(name)s const >, utility::COP_to_Python_converter< %(context)s%(name)s >, false >();\n\n' % D

            # for i in convert_to:
            #     r += '  boost::python::implicitly_convertible< utility::pointer::owning_ptr< %(context)s%(name)s >\n' % D
            #     r += '                                       , utility::pointer::owning_ptr< %s > >();\n\n' % i

            # if use_callback_struct:
            #     r += '  boost::python::implicitly_convertible< utility::pointer::owning_ptr< %(callback)s >\n' % D
            #     r += '                                       , utility::pointer::owning_ptr< %(context)s%(name)s > >();\n\n' % D


            # Adding COP → OP conversion
            r += '  boost::python::to_python_converter< ::boost::shared_ptr< %(context)s%(name)s const >, utility::py::COP_to_Python_converter< %(context)s%(name)s >, false >();\n\n' % D

            for i in convert_to:
                r += '  boost::python::implicitly_convertible< ::boost::shared_ptr< %(context)s%(name)s >\n' % D
                r += '                                       , ::boost::shared_ptr< %s > >();\n\n' % i

            if use_callback_struct:
                r += '  boost::python::implicitly_convertible< ::boost::shared_ptr< %(callback)s >\n' % D
                r += '                                       , ::boost::shared_ptr< %(context)s%(name)s > >();\n\n' % D

        return r


    #def getCallbackClassName(self): return 'PY_'+ self.name
    def getCallbackClassName(self): return self.name

    def getMangledName(self, name): return name.replace(' ', '_').replace('<', '_T_').replace('>', '_T').replace('::', '_').replace(',', '_').replace('*', '_st_')

    def getExposerName(self, name): return self.getMangledName(name) + '_exposer'

    def wrap_prefix_code(self, indent=''):
        r = ''
        # Callback Struct →
        if self.isCallbackStructNeeded():
            callback= self.getExposerName(self.getCallbackClassName()) + '_callback'

            D = dict(name=self.name, context=self.context, callback=callback)

            r += 'struct %(callback)s : public %(context)s%(name)s, boost::python::wrapper< %(context)s%(name)s >\n{\n' % D

            #for f in self.constructors+self.functions:
            #    r += f.writeCallbackCode(class_=self.context+self.name, callback=callback, indent='  ')

            for f in self.constructors + self.getAllVirtualFunctions({}).values():
                r += f.writeCallbackCode(class_=self.context+self.name, callback=callback, indent='  ')

            '''
            # now lets iterate thrugh all public bases and see if there is pure virtual functions that we need to implement...
            for b in self.bases:
                if b.type_.T() in self.reference.Objects:
                    o = self.reference.Objects[ b.type_.T() ]
                    for f in o.functions:
                        if f.pure_virtual and  ('%s(%s)' % (f.name , f.getSimpleArgsType()) ) not in written:
                            print written
                            r += f.writeCallbackCode(class_=self.context+self.name, callback=callback, indent='  ')
                            written.append( '%s(%s)' % (f.name , f.getSimpleArgsType()))
            '''

            r += '};\n'
        # ← Callback Struct

        return ('\n' + indent).join( r.split('\n') ) + '\n'

    def wrap(self, indent=''):
        #r = self.wrap_(indent, use_callback_struct=False)
        #if self.isCallbackStructNeeded(): r += self.wrap_(indent, use_callback_struct=True)
        r = self.wrap_(indent, use_callback_struct=self.isCallbackStructNeeded())
        return r

    def wrap_(self, indent='', use_callback_struct=False):
        if not self.isWrapable(): return ''

        print 'WrappingClass %s %s... [' % (self.context, self.name),
        for b in self.bases: print b.type_.T(),
        print ']', '[I]' if use_callback_struct else ''

        self.markConstOverloaded()
        if use_callback_struct: class_name = self.getCallbackClassName()
        else: class_name = self.name

        default_constructors = filter(lambda x: x.default, self.constructors)
        default_constructor = default_constructors and default_constructors[0]  # must be always... at least one... - not really...
        exposer = self.getExposerName(class_name) #self.name.replace(' ', '_').replace('<', '_T_').replace('>', '_T').replace('::', '_').replace(',', '_') + '_exposer'
        callback=exposer+'_callback'

        if use_callback_struct: heldTypeBase = callback
        else: heldTypeBase = self.context+self.name

        D = dict(name=self.name, class_name=class_name, context=self.context, mangled_class_name=self.getMangledName(class_name), exposer=exposer,
                 callback=callback, heldType=self.getHeldType(heldTypeBase), doc=doxygen.getDoxygenComment(self.file_, self.line), file_=self.file_, line=self.line)

        r = '\n{ // %(context)s%(class_name)s file:%(file_)s line:%(line)s\n' % D  # making comments emacs friendly

        r += self.write_implicitly_convertible_code(use_callback_struct, D)

        r += '  utility::py::wrap_access_pointer< %(context)s%(name)s >("%(mangled_class_name)s");\n' % D

        if use_callback_struct: r += '  boost::python::class_< %s, boost::noncopyable >("__CPP_%s__", "", boost::python::no_init);\n\n' % (self.getHeldType(), self.name)

        if self.isCopyable():  # or self.isHeldTypeOP():
            r += '  typedef boost::python::class_< %(heldType)s > %(exposer)s_type;\n' % D
        else:
            r += '  typedef boost::python::class_< %(heldType)s, boost::noncopyable > %(exposer)s_type;\n' % D

        if default_constructor:
            if default_constructor.argsTypes:
                default_constructor_args = '< %s >( (%s) , "%s" )' % (default_constructor.getSimpleArgsType(constructor=True), default_constructor.getArgTypes(), doxygen.getDoxygenComment(default_constructor.file_, default_constructor.line))
            else:
                default_constructor_args = '< %s >()' % default_constructor.getSimpleArgsType(constructor=True)

        #print ' default_constructor and self.isCreatable()...........', default_constructor, self.isCreatable()
        #r += '// mark 1 \n'
        if (default_constructor and self.isCreatable())  or (default_constructor and self.isCreatable(with_callback_struct=use_callback_struct)) :
            r += '  %(exposer)s_type %(exposer)s("%(mangled_class_name)s", "%(doc)s", boost::python::init %(default_constructor_args)s );\n' %  dict(D.items(), default_constructor_args=default_constructor_args )
        elif self.isCreatable() or self.isCreatable(with_callback_struct=use_callback_struct):
            r += '  %(exposer)s_type %(exposer)s("%(mangled_class_name)s", "%(doc)s" );\n' %  D
        else:
            r += '  %(exposer)s_type %(exposer)s("%(mangled_class_name)s", "%(doc)s", boost::python::no_init );\n' %  D
        #r += '// mark 2 \n'

        if default_constructor and self.isCreatable():
            for c in self.constructors:
                if (c is not default_constructor) and  c.public:  #  and (not c.artificial):
                    if c.argsTypes:
                        r += '  %(exposer)s.def( boost::python::init< %(args)s > ( (%(e_args)s) , "%(constr_doc)s" ) );\n' % dict(D.items(), args=c.getSimpleArgsType(constructor=True), e_args=c.getArgTypes(), constr_doc=doxygen.getDoxygenComment(c.file_, c.line))
                    else:
                        r += '  %(exposer)s.def( boost::python::init< %(args)s > () );\n' % dict(D.items(), args=c.getSimpleArgsType() )



        if self.enums:
            r+= '    boost::python::scope %(name)s_scope( %(exposer)s );\n' % D
            for e in self.enums:
                r+= e.wrap(indent='  ')

        wrapped = []
        statics = {}

        all_functions = self.functions
        if self.isSelfInherit():
            print '\033[35m\033[1mSelf inheritance detected! [%s]\033[0m' % self.isSelfInherit()
            # if b.type_.T() in self.reference.Objects:  - skipping test for presense for now - objetc _must_ be there...
            o = self.reference.Objects[ self.isSelfInherit() ]
            all_functions = self.functions + o.functions

        for f in all_functions:  #self.functions:
            r+= f.wrap(indent='  ', wrappingScope=exposer + '.', useCallbackStruct=['', callback][use_callback_struct])  #, overloaded=self.isOverloaded(f) )
            if f.static and f.isWrapable():  statics[f.name] = '  %(wrappingScope)sstaticmethod("%(py_name)s");\n' % dict(wrappingScope=exposer + '.', py_name=f.name)
            wrapped.append( str(f) )

        if use_callback_struct:
            F = self.getAllVirtualFunctions({})
            for f in F:
                if f not in wrapped:
                    r+= F[f].wrap(indent='  ', wrappingScope=exposer + '.', useCallbackStruct=['', callback][use_callback_struct])  #, overloaded=self.isOverloaded(f) )
                    wrapped.append( f )


        for v in self.isPublicMembersWrappable():
            r += v.wrap(indent='  ', wrappingScope=exposer + '.')
            #for v in self.dataMembers:  r += v.wrap(indent='  ', wrappingScope=exposer + '.')
            # commenting out for now because there is no reliable way check if object is wrappable or not

        for i in statics.values():
            r += '\n' + i + '\n'

        #r += self.wrapSpecialMethods()
        # str, looking for 'ostream& operator<<(ostream &, T)
        str_op_name = self.context[2:] + 'operator<<(std::ostream&, %s const&)' % (self.context[2:]+self.name)  # core::conformation::operator<<(std::ostream&, core::conformation::Residue const&)
        #print str_op_name
        if str_op_name in self.reference.Objects:
            if self.reference.Objects[str_op_name].returnType.T()== '::std::ostream &':
                #_o = self.reference.Objects[str_op_name]
                #print '.........', _o.artificial, _o.demangled, _o.file_, _o.line
                r += '  %(exposer)s.def( boost::python::self_ns::str( boost::python::self ) );\n' % D

        #for o in self.reference.Objects.values():
        #    if o.name.startswith('<<'): print '~~~~ ', o

        r += '}\n'
        return ('\n' + indent).join( r.split('\n') ) + '\n'

    def getFileList(self):  # return list of files where types was created (use it for includes!)
        r = [self.file_]
        for f in self.constructors+self.functions: r += f.getFileList()
        for t in self.requiredTypes: r.append( t.getFile() )

        for b in self.bases:
            if b.type_.T() in self.reference.Objects: r += self.reference.Objects[ b.type_.T() ].getFileList()

        return r

    def __repr__(self):
        return 'CppClass: %s %s\n' % (self.name, self.context) + \
               '  constructors: %s\n' % self.constructors + \
               '  function: %s\n' % self.functions


    '''  To do:
             .staticmethod( "get_instance" );              core/pack/task/operation/_ResFilterFactory.cc
             vars: core/id/_TorsionID.cc
        .def( bp::self != bp::self )
        .def( bp::self < bp::self )
        .def( bp::self_ns::str( bp::self ) )
        .def( bp::self == bp::self )
        .def( bp::self_ns::str( bp::self ) );

    bp::scope().attr("BOGUS_TORSION_ID") = core::id::BOGUS_TORSION_ID;


'''

#
#
#
class GccXML:
    def __init__(self, dom):
        self.dom = dom  # Do we even need this???

        self.Objects = {}     # storage of all objects indexed by C++ names
        self.Contexts = {}    # storage of all objects indexed by C++ namespaces
        self.Namespaces = {}  # storage of CppNamespace's by XML id's
        self.Files = {}       # storage of file names (strings) indexed by XML id's
        self.TypeNodes, self.Types = {}, {}  # storage of types nodes and types (strings) indexed by XML id's
        self.Nodes = {}       # storage for everything else...

    def getType(self, type_id):
        if type_id not in self.Types:
            if type_id not in self.TypeNodes: self.Types[type_id] = CppType_Fundamental( 'UnknowType_%s' % type_id, 'unknow')
            else:
                node = self.TypeNodes[type_id]
                if   node.nodeName == 'FundamentalType':   self.Types[type_id] = CppType_Fundamental( node.getAttribute('name'), node.nodeName )

                elif node.nodeName == 'PointerType':       self.Types[type_id] = CppType_Simple( self.getType( node.getAttribute('type')), '*', node.nodeName )
                elif node.nodeName == 'ReferenceType':     self.Types[type_id] = CppType_Simple( self.getType( node.getAttribute('type')), '&', node.nodeName )
                elif node.nodeName == 'CvQualifiedType':   self.Types[type_id] = CppType_Simple( self.getType( node.getAttribute('type')), 'const', node.nodeName )

                elif node.nodeName in ['ArrayType', 'FunctionType']:
                    #print 'Creating composite', type_id, node.nodeName
                    self.Types[type_id] = CppType_Composite( node.nodeName, node.nodeName )
                    #print 'Creating composite', type_id, node.nodeName, self.Types[type_id].T()

                elif node.nodeName == 'Typedef':  # maybe not preatty and certainly verbose - but it fix everything...
                    #return self.getType( node.getAttribute('type') )
                    self.Types[type_id] = CppType_Typedef( node.getAttribute('name'), self.Nodes[node.getAttribute('context')].getAttribute('demangled'), self.getType( node.getAttribute('type') ) )

                elif node.nodeName in ['Class', 'Struct', 'Union', 'Typedef', 'Enumeration']:
                    # old... #self.Types[type_id] = CppType_Complex( '::' + node.getAttribute('demangled'), '::' + self.Nodes[node.getAttribute('context')].getAttribute('demangled') )

                    self.Types[type_id] = CppType_Complex( node.getAttribute('name'), self.Nodes[node.getAttribute('context')].getAttribute('demangled'), node.nodeName )

                    #self.Types[type_id] = CppType_Complex( node.getAttribute('demangled'), self.Nodes[node.getAttribute('context')].getAttribute('demangled'), node.nodeName )
                    #self.Types[type_id].name = self.Types[type_id].name[len(self.Types[type_id].context)+2:]   # removing namespace prefix
                    #print 'Creating composite', type_id, node.nodeName, self.Types[type_id].name, self.Types[type_id].context

                else: self.Types[type_id] = CppType_Fundamental( 'UnknowType_%s_%s' % (node.nodeName, type_id), node.nodeName) #  + self.getType( node.getAttribute('type') )


                '''
                elif node.nodeName == 'ArrayType':         self.Types[type_id] = '___XQWERTY___Array_%s'  % self.getType( node.getAttribute('type') )
                elif node.nodeName == 'FunctionType':
                    args = []
                    for ch in node.childNodes:
                        if ch.nodeName == 'Argument': args.append( self.getType(ch.getAttribute('type')) )

                    self.Types[type_id] = '___XQWERTY___FunctionType_%s(*)(%s)' % (self.getType( node.getAttribute('returns') ), ', '.join(args) )
                '''
                file_ = node.hasAttribute('file') and self.Files[node.getAttribute('file')]
                line  = node.getAttribute('line')
                self.Types[type_id].file_ = file_
                self.Types[type_id].line  = line

            #print self.Types[type_id].T(), self.Types[type_id].getFile()

        #print type_id, self.Types[type_id].T(), self.Types[type_id].getContext()
        return self.Types[type_id]


    def parse(self, relevantFilesList=None):
        ''' Parse GCC_XML, we do this in two passes, first: namespaces, files and creating Dict of all elements.
            second - parsing actual elements (functions, classes)
        '''
        # Pass №1 Indexing Namespace's and File's
        print 'GccXML parsing, Pass №1...'
        for c in self.dom.documentElement.childNodes:  # dom.documentElement suppose to be dom.getElementsByTagName(GCC_XML')
            if c.nodeName == '#text': continue
            c_id = c.getAttribute('id')
            if c.nodeName in ['FundamentalType', 'PointerType', 'ReferenceType', 'CvQualifiedType', 'Class', 'Struct', 'Union', 'Typedef', 'FunctionType', 'Enumeration', 'ArrayType', 'OffsetType', 'MethodType']: self.TypeNodes[c_id] = c

            elif c.nodeName == 'File':
                self.Files[c_id] = c.getAttribute('name')
                if self.Files[c_id].startswith('./'): self.Files[c_id] = self.Files[c_id][2:]

            #else: self.Nodes[c_id] = c
            self.Nodes[c_id] = c  # we might want to have full reference for simplicity sake...

        print 'GccXML parsing, Pass №2...'
        for c in self.dom.documentElement.childNodes:  # dom.documentElement suppose to be dom.getElementsByTagName(GCC_XML')
            if c.nodeName == '#text': continue
            c_id = c.getAttribute('id')
            if c.nodeName in ['Namespace', 'Class', 'Struct']:
                #if c.hasAttribute('context') : context = self.Namespaces[c.getAttribute('context')].getChildrenContext()
                if c.getAttribute('context') in self.Namespaces : context = self.Namespaces[c.getAttribute('context')].getChildrenContext()
                else: context = ''
                self.Namespaces[c_id] = CppNamespace(c.getAttribute('name'), context)


        relevantFilesList = relevantFilesList or self.Files.values()

        # Pass №3 Actually doing the job...
        print 'GccXML parsing, Pass №3...'
        for node in self.dom.documentElement.childNodes:
            if node.nodeName == '#text': continue
            #if (not node.hasAttribute('file'))  or  (self.Files[node.getAttribute('file')] not in relevantFilesList) : continue

            object_ = None

            if node.nodeName == 'Enumeration':
                context =  '::' + self.Nodes[node.getAttribute('context')].getAttribute('demangled') + '::'
                object_ =  CppEnum(node, context, self)

            if node.nodeName in ['Function', 'OperatorFunction']:
                context =  self.Namespaces[node.getAttribute('context')].getChildrenContext()
                type_ = self.getType( node.getAttribute('returns') )
                object_ =  CppFunction(node, context, type_, self)

            if node.nodeName in ['Class', 'Struct']:
                #if node.nodeName == 'Struct' and node.getAttribute('artificial') == '1': continue

                if node.getAttribute('context') not in self.Namespaces: continue  # Clearly something not right - we probably hitting default/inner GCC classes here...
                context =  self.Namespaces[node.getAttribute('context')].getChildrenContext()
                object_ = CppClass(node, context, self)
                members = node.getAttribute('members')[:-1].split(' ')
                ch_context = object_.getChildrenContext()
                #print node.getAttribute('name'), members
                for m_id in [ x for x in members if x]:
                    m = self.Nodes[m_id]
                    #if m.getAttribute('access') != "public": continue
                    if m.nodeName == 'Constructor':  # and (not m.hasAttribute('artificial') ):
                        object_.constructors.append( CppFunction(m, ch_context, CppType_Fundamental('', ''), self, constructor=True) )
                        object_.constructors[-1].default = m.getAttribute('access') == "public"  and (not m.hasAttribute('artificial'))  # m.hasAttribute('explicit')  # and  ch.getAttribute('is_default_ctor') == '1'

                    if  m.nodeName == 'Destructor': object_.destructor = CppFunction(m, ch_context, None, self)

                    if m.nodeName in ['Method', 'OperatorMethod']:
                        #if m.getAttribute('access') != "public": continue
                        ch_type_ = self.getType( m.getAttribute('returns') )
                        object_.functions.append( CppFunction(m, ch_context, ch_type_, self, memberFunction=True, const=m.hasAttribute('const') ) )

                    if m.nodeName == 'Enumeration' and (not m.getAttribute('name').startswith('$_') ):
                        object_.enums.append( CppEnum(m, ch_context, self) )
                    '''
                    if m.nodeName == 'Field':  # we actually looking for enum here... this it to handle old style: enum {a, b} something:
                        en_node = self.Nodes[ m.getAttribute('type') ]
                        if en_node.nodeName == 'Enumeration':
                            object_.enums.append( CppEnum(en_node, ch_context, self, name=m.getAttribute('name') ) )
                    '''
                    if m.nodeName == 'Field':  # Ok, this could be inner type (typedef, enum etc) or datamemeber
                        object_.requiredTypes.append( self.getType( m.getAttribute('type') ) )

                        # TODO: add some code here to distinguish between datamembers and typedefs
                        object_.dataMembers.append( CppVariable(m, ch_context, self) )




            if object_:
                if not self.Contexts.get(context): self.Contexts[context] = [ object_ ]
                else: self.Contexts[ context ].append( object_ )
                #print '~~~~~', object_.context + object_.name
                self.Objects[ object_.context + object_.name ] = object_
                self.Objects[ node.getAttribute('demangled') ] = object_  # functions could be overloaded...


        print 'GccXML parsing... Done!'


def generateIncludes(incl_list):
    incl_list = list( set( filter(bool, incl_list) ) )
    #print incl_list
    for i, f in enumerate(incl_list):
        if f.endswith('.fwd.hh'):
            if os.path.isfile(f[:-7]+'.hh'): incl_list[i] = f[:-7]+'.hh'

    if 'core/types.hh' in incl_list:  # ＃　→ばか！！！
        incl_list.append('numeric/xyzVector.hh')
        incl_list = list( set( incl_list ) )

    r = ''
    for f in incl_list:
        if f.split('/')[0] in ['utility', 'numeric', 'basic', 'ObjexxFCL', 'core', 'protocols']:
            r += '#include <%s>\n' % f

    return r


def sortObjects(l):
    ''' sort list of C++ objects for wrapping
    '''
    def swap(i, j):
        #print l[i], l[j]
        a = l[i];  l[i]=l[j];  l[j]=a

    f = True
    while f:
        f = False
        for i in range( len(l) ):
            for j in range( i, len(l) ):
                if isinstance(l[i], CppClass) and isinstance(l[j], CppClass) and 'utility::SingletonBase' in l[i].context+l[i].name and 'utility::SingletonBase' in l[j].context+l[j].name:
                    if l[i].context+l[i].name > l[j].context+l[j].name: swap(i, j); f = True;
                    continue

                if isinstance(l[j], CppClass) and 'utility::SingletonBase' in l[j].context+l[j].name: swap(i, j); f = True; continue
                if isinstance(l[i], CppClass) and 'utility::SingletonBase' in l[i].context+l[i].name: continue
                    #if isinstance(l[j], CppClass): # some special cases (default parameter bechavor)

                if isinstance(l[i], CppClass) and isinstance(l[j], CppClass):
                    if l[i].isNeeded(l[j]) and not l[j].isNeeded(l[i]):  swap(i, j); f = True; break
                    #if 'utility::SingletonBase' in l[j].context+l[j].name: break

                if isinstance(l[i], CppClass) and isinstance(l[j], CppEnum): swap(i, j); f = True; break

                if isinstance(l[i], CppFunction) and isinstance(l[j], CppEnum): swap(i, j); f = True; break

                if isinstance(l[i], CppClass) and isinstance(l[j], CppClass): # some special cases (default parameter bechavor)
                    #print l[i].getChildrenContext(), l[j].getChildrenContext()
                    if l[i].getChildrenContext()=='::core::chemical::Atom::' and l[j].getChildrenContext()=='::core::chemical::AtomICoor::': swap(i, j); f = True; break

                    #if 'utility::SingletonBase' in l[i].context+l[i].name  and  'utility::SingletonBase' not in l[j].context+l[j].name: swap(i, j); f = True; break
                    #if 'utility::SingletonBase' in l[i].context+l[i].name  and  'utility::SingletonBase' in l[j].context+l[j].name:
                    #   if l[i].context+l[i].name > l[j].context+l[j].name: swap(i, j); f = True; break



def wrapModule(name, name_spaces, context, relevant_files_list, max_funcion_size, by_hand_beginning='', by_hand_ending='', monolith=False):
    ''' Template for creating one module, and wrapping each elelemnts... (each elements must have .wrap)
    '''
    #print 'by_hand_beginning=%s, by_hand_ending=%s' % (by_hand_beginning, by_hand_ending)

    #print '------ name:', name
    #print '------ name_spaces:', name_spaces
    #print '------ context:', context

    objects = []
    for n in name_spaces:
        if n in context: objects.extend( context[n] )

    if '::utility::' not in name_spaces:
        for n in context:
            if n.startswith('::u'):
                for o in context[n]:
                    if isinstance(o, CppClass):
                        if '::utility::SingletonBase<' in o.context+o.name:
                            print 'Adding SingletonBase instance:', o.context+o.name
                            objects.extend( [o] )

    if name_spaces == ['::core::conformation::']:  # work-aroung for GCCXML template-misplace-bug
        for o in context['::core::graph::']:
            if o.name == 'UpperEdgeGraph<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>':
                objects.append( o );
                #print '######################## ', o

    sortObjects( objects )

    code = []

    module_addon = '#include <utility/py/PyHelper.hh>\n'

    r  = module_addon+'\n'

    s = ''
    includes = []
    prefix_code = ''
    for i in objects:
        force = i.name == 'UpperEdgeGraph<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>'
        #if force: print '--- yay: ', i.name

        if i.file_ in relevant_files_list  or force:
            s += i.wrap('  ')
            prefix_code += i.wrap_prefix_code()

            includes += i.getFileList()


        if len(s)+len(prefix_code)+len(includes) >= max_funcion_size:
            r += 'void %s_partial_%s(void);\n' % (name, len(code))
            code.append( module_addon+'%s\n\n%s\nvoid %s_partial_%s(void)\n{\n%s\n}\n\n' % (generateIncludes(includes), prefix_code, name, len(code), s)  )
            s, prefix_code, includes = '', '', []


    r += by_hand_beginning + by_hand_ending
    #r += '%s\n\n%s\nBOOST_PYTHON_MODULE( %s ) {\n' % (generateIncludes(includes), prefix_code, name)

    module_prefix = '{includes}\n\n{prefix}\nvoid _import_namespace{name}() {{\n' if monolith else '{includes}\n\n{prefix}\nBOOST_PYTHON_MODULE( {name} ) {{\n'

    r += module_prefix.format(includes=generateIncludes(includes), prefix=prefix_code, name=name)

    if by_hand_beginning:
        print '\033[33m\033[1mAdding by_hand_beginning code for %s::%s\033[0m' % (name_spaces, name)
        r += '  __%s_by_hand_beginning__();\n\n' % name_spaces[0].split('::')[-2]

    if by_hand_ending:
        print '\033[33m\033[1mAdding by_hand_ending code for %s::%s\033[0m' % (name_spaces, name)
        s = s + '\n' + '__%s_by_hand_ending__();\n' % name_spaces[0].split('::')[-2]

    #r += '%s\n\n%s\n#ifndef __PYROSETTA_ONE_LIB__\n  BOOST_PYTHON_MODULE( %s ) {\n#else\n  void __wrap%s() {\n#endif\n' % (generateIncludes(includes), prefix_code, name, name_spaces[0].replace('::', '__'))
    for i in range( len(code) ): r += '\n  %s_partial_%s();\n' % (name, i)
    r += '\n' + s + '}\n'

    return code+[r]


def parseAndWrapModule(module_name, namespaces_to_wrap, xml_source, relevant_files_list, ParserType=GccXML, max_funcion_size=1024*1024*1024,
                       by_hand_beginning='', by_hand_ending='', monolith=False):
    print 'Wrapping %s... with file list as: %s...' % (namespaces_to_wrap, relevant_files_list)
    relevant_files_list.append('utility/SingletonBase.hh')  # special case: self self-inheritance template
    for f in relevant_files_list[:] :  relevant_files_list.append('./'+f)

    print 'Parse XML using mini dom...'
    dom = xml.dom.minidom.parse(xml_source)
    #dom =  xml.dom.pulldom.parseString(xml_source)

    #import xml.dom.ext.reader
    #reader = xml.dom.ext.reader.Sax2.reader()
    #dom = reader.fromStream(xml_source)

    print 'Parse XML using mini dom... Done!'
    cxml = ParserType(dom)
    cxml.parse(relevantFilesList=relevant_files_list)
    #return wrapModule(module_name, cxml.Contexts[namespace_to_wrap], cxml.Contexts)
    res = wrapModule(module_name, namespaces_to_wrap, cxml.Contexts, relevant_files_list, max_funcion_size, by_hand_beginning=by_hand_beginning, by_hand_ending=by_hand_ending, monolith=monolith)

    del cxml  # trying to save memory...
    del dom
    gc.collect()

    return res


# --------------------------------------- Monolith related funtions ---------------------------------------
_import_order_ = {
    '' : ['utility', 'numeric', 'basic', 'core', 'protocols'],

    'basic': ['datacache', 'resource_manager'],

    'core': ['graph', 'conformation', 'id', 'io', 'scoring'],
    'core/scoring': ['trie', 'methods', 'func', 'hbonds'],

    'protocols': ['moves', 'jd2', 'jumping', 'environment', 'features', 'evaluation', 'canonical_sampling', 'farna', 'filters',
                  'simple_moves', 'ligand_docking', 'rigid', 'toolbox', 'forge', 'loop_modeling', 'loops', 'wum', 'rosetta_scripts',
                  'simple_filters'],
    'protocols/forge' : ['remodel'],
    'protocols/match' : ['upstream'],
    'protocols/stepwise' : ['screener', 'modeler'],
    #'protocols/loop_modeling' : ['LoopMover'],
}

_import_after_sub_namespaces_ = ['protocols/forge/build'] #'protocols/kinematic_closure']  #, 'protocols/loop_modeling'] 'protocols/rotamer_recovery' 'protocols/stepwise/sampling/rna'

def get_direct_childrens(path, modules):
    ''' Find direct childrens of given module and sort it if possible by propper import order '''
    childrens =  [m for m in modules if path == os.path.split(m.path)[0] ]
    childrens.sort(key=lambda x: x.name )  # we always sort to avoid error because of FS files order

    first = []
    for m in childrens[:]:
        if path in _import_order_  and  m.name in _import_order_[path]:
            first.append(m); childrens.remove(m)


    first.sort(key=lambda x: _import_order_[path].index(x.name) )
    if first: print 'Adjusting import order: {0} --> {1}'.format(path, [m.name for m in first]) #, [m.name for m in childrens])
    return first + childrens


def generate_monolith_module(module, modules, indent='', use_python_code='', python_module_name=''):
    #submodule_cpp_name = '_pyrosetta_{}_submodule_'.format(name)

    name = python_module_name if python_module_name else module.name
    r  = ''
    r += '{{ // {0}\n'.format( python_module_name if python_module_name else module.all_at_once_base )
    r += '  boost::python::scope current;\n'
    r += '  std::string submodule_name( boost::python::extract<const char*>(current.attr("__name__")));\n'
    r += '  submodule_name.append(".{name}");\n'.format(name=name)
    #r += '  std::cout << submodule_name << std::endl;\n'

    r += '  boost::python::object submodule( boost::python::borrowed( PyImport_AddModule(submodule_name.c_str()) ) );\n'

    #r += '  boost::python::object submodule( boost::python::handle<>( boost::python::borrowed( PyImport_AddModule(submodule_name.c_str()) ) ) );\n'

    r += '  current.attr("{name}") = submodule;\n'.format(name=name)
    r += '  submodule.attr("__package__") = "{name}";\n'.format(name=name)
    r += '  boost::python::extract<boost::python::dict>(boost::python::getattr(boost::python::import("sys"),"modules"))()[submodule_name]=submodule;\n'
    #r += '  boost::python::extract<boost::python::dict>(boost::python::getattr(boost::python::import("sys"),"modules"))()["{}"]=submodule;\n'.format(module.path.replace('/','.')) # also append path without root as imported. This should enable realative imports

    #r += '  submodule.attr("__file__") = "<synthetic>";\n'
    #r += '  boost::python::extract<boost::python::dict>(boost::python::getattr(boost::python::import("sys"),"modules"))()[submodule_name]="{name}";'.format(name=name)
    r += '  \n'

    r += '  {\n'  # this inner scope is absolutely nessesary here because we want to make swtiching to submodule context temporary. DO NOT REMOVE!!!!!
    r += '    boost::python::scope submodule_scope(submodule);\n'

    if use_python_code:
        r += '    // Embedding: {0}\n'.format(python_module_name)
        #r += '    std::cout << "{}" << std::endl;\n'.format(use_python_code)
        r += '    boost::python::exec("{0}");\n'.format(use_python_code)
    else:
        self_import = '    _import_namespace{base_name}();\n'.format(base_name=module.all_at_once_base)
        if module.path not in _import_after_sub_namespaces_: r+= self_import
        for m in get_direct_childrens(module.path, modules): r+= generate_monolith_module(m, modules, indent=indent+'  ')
        if module.path in _import_after_sub_namespaces_: r+= self_import; print '↓{0}'.format(module.all_at_once_base)

    r += '  }\n'
    r += '}\n'
    return '\n'.join( [ indent+line if line else line for line in r.split('\n')] )

def generate_monolith_main(root_module, modules, rosetta_library_name, embed_python):
    ''' Geneate main file for monolith build. This function have nothing to do with rest CppParse but its here to keep all generate* funtion together
    '''
    r  = '// Monolith main\n'
    #r += '\n#include <boost/python.hpp>\n'
    r += '\n#include <boost/python.hpp>\n#include <iostream>\n'
    if embed_python:
        r += '#include <boost/python/stl_iterator.hpp>\n'
        r += 'boost::python::object boost_python_dir(boost::python::object object)\n'
        r += '{\n'
        r += '  boost::python::handle<> handle(PyObject_Dir(object.ptr()));\n'
        r += '  return boost::python::object(handle);\n'
        r += '}\n'


    # generate all prototypes
    for m in modules: r += 'void _import_namespace{base_name}();\n'.format(base_name=m.all_at_once_base)

    #r += 'int test_function() { return 42; }\n'
    r += '\n'
    r += 'BOOST_PYTHON_MODULE( {0} )\n{{\n'.format(rosetta_library_name)

    # Disabling Python duplicate warnings
    #r += '  std::string disable_warning("{}");\n'.format("""import warnings\\nwarnings.filterwarnings(\\"ignore\\", \\"to-Python converter for .+ already registered; second conversion method ignored.\\", RuntimeWarning, \\"^rosetta\\\\\\\\.\\")""")
    r += '  std::string disable_warning("{0}");\n'.format("""import warnings\\nwarnings.filterwarnings(\\"ignore\\", \\"to-Python converter for .+ already registered; second conversion method ignored.\\", RuntimeWarning, \\"\\")""")
    #r += '  std::cout << disable_warning << std::endl;\n'
    r += '  boost::python::exec( disable_warning.c_str() );\n'

    #r += 'boost::python::def("test_function", test_function);\n'
    for m in get_direct_childrens('', modules): r+= generate_monolith_module(m, modules, indent='  ')

    # r += '  boost::python::object module = boost::python::import("__main__");\n'
    # r += '  boost::python::object name_space = module.attr("__dict__");\n'
    # r += '  boost::python::exec("print 123\\nimport rosetta3.utility\\nimport utility", name_space, name_space);\n'

    def escape_to_python(s): return s.replace('\\', '\\\\').replace('\n', '\\n').replace('"', '\\"')  # we can't use re.escape beacuse it C++ don't like some of it esacpes...

    # it maybe hard to believe but apparatnly M$ compilers have *hard* limit on length of string literals (WHO ON EARTH WROTE THAT?????)... so we have to construct resulted string on the fly...
    def generate_exec_code(python_code, namespace):
        # in GCC this would be: boost::python::exec("{0}", namespace);\n'.format( escape_to_python(python_code) )
        r = '{\n' \
             '  std::string python_code;\n'

        python_code = python_code.split(' ')
        number_of_lines_to_add = 2048;
        cpp_lines = [ ' '.join(python_code[i:i+number_of_lines_to_add]) for i in range(0, len(python_code), number_of_lines_to_add)]
        for i, line in enumerate(cpp_lines): r += '  python_code += "{0}{1}";\n'.format(' ' if i else '', escape_to_python(line))
        r += '  boost::python::exec(python_code.c_str(), {namespace});\n'.format(namespace=namespace)
        r += '}\n'

        return r


    if embed_python:
        # embed some python files inside our lib
        for fn in sorted(glob.glob(root_module.binding_source_path + '/src/*.py'), key = lambda n: 'z'+n if n.find('__init__')>0 else n ):  # making init last module
            python_module_name = os.path.split(fn)[1][:-3] # removing path and '.py' suffix

            #if python_module_name == '__init__': python_module_name = 'init_module' # we put default python module in to special namespace rosetta.init

            if python_module_name != '__init__':  # skipping because we will embed it without namespace below
                #if python_module_name in ['version']:
                    print 'Embedding:', python_module_name
                    #r+= generate_monolith_module(root_module, modules, indent='  ', use_python_code=escape_to_python(file(fn).read()), python_module_name=python_module_name)
                    python_string = "import sys,imp\n" \
                                    "{python_module_name} = imp.new_module('{python_module_name}')\n" \
                                    "exec \"{code}\" in {python_module_name}.__dict__\n" \
                                    "sys.modules['{python_module_name}']={python_module_name}\n" \
                                    .format(python_module_name=python_module_name, code=escape_to_python( file(fn).read() ) )
                                    # "print \"{python_module_name}\"\n" \
                                    # "setattr(sys.modules[__name__], '{python_module_name}', {python_module_name} )\n" \
                                    #"__name__['{python_module_name}'] = {python_module_name}\n" \

                    # r += '  boost::python::exec("{}");\n'.format( escape_to_python(python_string) )
                    r += '{{ // Embedding: {0}\n'.format(python_module_name)
                    r += '  boost::python::scope current;\n'

                    r += '  std::string submodule_name( boost::python::extract<const char*>(current.attr("__name__")));\n'
                    r += '  submodule_name.append(".{name}");\n'.format(name=python_module_name)
                    #r += '  std::cout << submodule_name << std::endl;\n'

                    r += '  boost::python::object main = boost::python::import("__main__");\n'
                    r += '  boost::python::object main_namespace = main.attr("__dict__");\n'


                    #r += '  boost::python::exec("{0}", main_namespace);\n'.format( escape_to_python(python_string) )
                    r += generate_exec_code(python_string, 'main_namespace')

                    #r += '  boost::python::object submodule = boost::python::extract<boost::python::object>( main_namespace.attr("{python_module_name}") );\n'.format(python_module_name=python_module_name)
                    #r += '  current.attr("{name}") = submodule;\n'.format(name=python_module_name)
                    r += '  current.attr("{name}") = main.attr("{name}");\n'.format(name=python_module_name)
                    r += '  boost::python::extract<boost::python::dict>(boost::python::getattr(boost::python::import("sys"),"modules"))()[submodule_name]=main.attr("{name}");\n'.format(name=python_module_name)
                    r += '}\n'

        init_file = file(root_module.binding_source_path + '/src/__init__.py').read()
        #init_file = init_file.replace('\\', '\\\\').replace('\n', '\\n').replace('"', '\\"')

        r += '{ // Embedding: __init__\n' # in short: we create dummy 'init' submodule, execute __init__ in it and the copy all object in to root namespace
        r += '  boost::python::scope current;\n'

        # # Approach when we copy some of the build-ins in to scope namespace...
        # r += '  boost::python::object main_namespace = current.attr("__dict__");\n'
        # # Now copy some Python build-in defaults in to current scope
        # r += '  boost::python::object submodule = boost::python::import("__main__");\n'
        # r += '  boost::python::object submodule_namespace = submodule.attr("__dict__");\n'
        # #r += '  boost::python::exec("import __builtin__\\n", submodule_namespace, submodule_namespace);\n'
        # #r += "  boost::python::object __builtin__ = boost::python::extract<boost::python::object>( boost::python::eval(\"__builtin__\", submodule_namespace) );\n"
        # #r += '  main_namespace["__builtin__"] = __builtin__;\n'
        # r += '  typedef boost::python::stl_input_iterator<boost::python::str> iterator_type;\n'
        # r += '  for (iterator_type name( submodule_namespace ), end; // for name in dir(object): \n'
        # #r += '  for (iterator_type name(boost_python_dir(__builtin__)), end; // for name in dir(object): \n'
        # r += '    name != end; ++name) {\n'
        # r += '    if( !name->startswith("__") || name->startswith("__import__") || true ) {\n'
        # r += '      std::cout << "Adding:" << boost::python::extract<const char*>(*name) << std::endl;\n'
        # #r += '      current.attr( *name ) = __builtin__.attr( *name );\n'
        # r += '      current.attr( *name ) = submodule_namespace[*name];\n'
        # r += '    } else {\n'
        # r += '      std::cout << "Skipping:" << boost::python::extract<const char*>(*name) << std::endl;\n'
        # r += '    }\n'
        # r += '  }\n'
        # r += '  boost::python::exec("{}", main_namespace, main_namespace);\n'.format(init_file)
        # #r += '  boost::python::exec("{}", main_namespace, main_namespace);\n'.format("import sys\\naaa = 'ABC'\\nprint aaa\\n", init_file)
        # # r += '  boost::python::exec("import __builtin__\\n", submodule_namespace, submodule_namespace);\n'
        # # r += '  boost::python::object __builtin__ = boost::python::extract<boost::python::object>( boost::python::eval("__builtin__", submodule_namespace) );\n'
        # # r += '  main_namespace["__builtin__"] = __builtin__;\n'
        # # r += '  boost::python::object __import__ = boost::python::extract<boost::python::object>( boost::python::eval("__import__", submodule_namespace) );\n'
        # # r += '  main_namespace["__import__"] = __import__;\n'



        # r += '  std::string submodule_name( boost::python::extract<const char*>(current.attr("__name__")));\n'
        # r += '  submodule_name.append(".pyrosetta");\n'

        # r += '  boost::python::exec("{}", init_module);\n'.format("aaa = 'ABC'\\nprint aaa\\n")

        #r += '  Py_Initialize();\n'
        #r += '  boost::python::object main_module = boost::python::import("__main__");\n'
        #r += '  boost::python::object main_namespace = main_module.attr("__dict__");\n'

        #r += '  boost::python::exec("{}", main_namespace, main_namespace);\n'.format("aaa = 'ABC'\\nprint aaa\\n", init_file)

        # r += '  boost::python::exec("{}", main_namespace, main_namespace);\n'.format(init_file)



        r += '  boost::python::object main = boost::python::import("__main__");\n'
        r += '  boost::python::object init_module = main.attr("__dict__");\n'

        #r += '  boost::python::exec("{0}", init_module);\n'.format( escape_to_python(init_file) )
        r += generate_exec_code(init_file, 'init_module')

        #r += '  boost::python::exec("from rosetta.init_module import *\\n", init_module);\n'
        r += '  typedef boost::python::stl_input_iterator<boost::python::str> iterator_type;\n'
        r += '  for (iterator_type name(init_module), end; // for name in dir(object): \n'
        #r += '  for (iterator_type name(boost_python_dir(main_namespace)), end; // for name in dir(object): \n'
        r += '    name != end; ++name) {\n'
        #r += '      std::cout << boost::python::extract<const char*>(*name) << std::endl;\n'
        r += '    if( !name->startswith("__") ) {\n'
        #r += '      std::cout << boost::python::extract<const char*>(*name) << std::endl;\n'
        r += '      current.attr( *name ) = init_module[ *name ];\n'
        r += '    }\n'
        r += '  }\n'
        #r += '  boost::python::extract<boost::python::dict>(boost::python::getattr(boost::python::import("sys"),"modules"))()["rosetta"]=init_module["rosetta"];\n'
        #r += '  current.attr("pyrosetta") = init_module;\n'
        #r += '  boost::python::import(pyrosetta);\n'
        r += '}\n'

    r += '}\n'

    return r


def parseTranslationUnit(TU):
    for i in TU.childNodes:
        print 'nodeName', i.nodeName, i.namespaceURI


def CLANG_XML(clang_xml):
    for n in clang_xml.childNodes:
        print n, dir(n)
        print 'nodeName', n.nodeName
        print '\n\n'

        if n.nodeName == 'TranslationUnit': parseTranslationUnit(n)
    '''
    tr_unit = clang_xml.getElementsByTagName('TranslationUnit')

    for i in tr_unit:
        print i
'''


def main(args):
    print'Starting...'
    #dom = xml.dom.minidom.parse('TestingXML.xml')
    dom = xml.dom.minidom.parse('Example2.xml')

    cl = dom.getElementsByTagName('CLANG_XML')
    print cl, dom.documentElement
    for i in cl:
        CLANG_XML(i)

    cxml = ClangXML(dom)

    cxml.parse(relevantFilesList=['Example2.cpp'])
    print cxml.Contexts

    print wrapModule('test', cxml.Contexts['::NM::'], cxml.Contexts)


def GCC_XML_main(args):
    print'Starting...'
    dom = xml.dom.minidom.parse('Example1.gccxml.xml')

    gxml = GccXML(dom)

    gxml.parse(relevantFilesList=['Example1.cpp'])

    for k in gxml.TypeNodes: gxml.getType(k)
    print 'Types:', gxml.Types
    print 'Files:', gxml.Files

    print wrapModule('test', '::NM::', gxml.Contexts)





class ClangXML:
    def __init__(self, dom):
        self.dom = dom  # Do we even need this???

        self.Contexts = {}

        for i in dom.documentElement.childNodes:  # dom.documentElement suppose to be dom.getElementsByTagName('CLANG_XML')
            if i.nodeName == 'TranslationUnit':  self.TranslationUnit = i
            if i.nodeName == 'ReferenceSection':  self.ReferenceSection = ReferenceSection(i)



    def parse(self, node=None, relevantFilesList=None):
        relevantFilesList = relevantFilesList or self.ReferenceSection.Files.values()
        node = node or self.TranslationUnit
        #print 'relevantFilesList:', relevantFilesList

        if node.nodeName == '#text': return
        if node.hasAttribute('file')  and  self.ReferenceSection.Files[node.getAttribute('file')] not in relevantFilesList: return

        #print '~~~', node

        context =  node.hasAttribute('context') and self.ReferenceSection.getContext( node.getAttribute('context') )
        type_ = node.hasAttribute('type') and self.ReferenceSection.getType( node.getAttribute('type') )

        object_ = None

        if node.nodeName == 'Namespace':
            object_ = CppNamespace(node.getAttribute('name'), context)

        elif node.nodeName == 'Function': object_ =  CppFunction(node, context, type_, self.ReferenceSection)

        elif node.nodeName == 'CXXRecord': object_ = CppClass(node, context, self.ReferenceSection)


        if object_:
            if not self.Contexts.get(context): self.Contexts[context] = [ object_ ]
            else: self.Contexts[ context ].append( object_ )

        #if node.nodeName == 'ReferenceSection':

        if node.nodeName in ['Namespace', 'TranslationUnit'] :  # we only need to automatically parse sub-nodes for namespaces...
            for i in node.childNodes:
                self.parse(i, relevantFilesList=relevantFilesList)


'''
    class method
    def parse_Clang_XML(self, node, context, refSection):
        # Now lets parse all child nodes...
        for ch in node.childNodes:
            if ch.nodeName == '#text': continue  # why on earth this node type is everywhere???
            #print 'Processing: ', ch
            ch_context =  ch.hasAttribute('context') and refSection.getContext( ch.getAttribute('context') )
            ch_type_ = ch.hasAttribute('type') and refSection.getType( ch.getAttribute('type') )

            if ch.nodeName == 'CXXConstructor':
                self.constructors.append( CppFunction(ch, ch_context, ch_type_, refSection) )
                self.constructors[-1].default = ch.hasAttribute('is_default_ctor') and  ch.getAttribute('is_default_ctor') == '1'

            elif ch.nodeName == 'CXXMethod':
                self.functions.append( CppFunction(ch, ch_context, ch_type_, refSection) )

'''

"""
    def __cmp__(self, other):
        #print 'ii? ', other.name, isinstance(other, CppClass)
        if isinstance(other, CppClass):
            print '%s.isBase(%s) --> %s' % (self.name, other.context + other.name, self.isBase(other.context + other.name) )
            print '%s.isBase(%s) --> %s' % (other.name, self.context + self.name, other.isBase(self.context + self.name) )
            if self.isBase(other.context + other.name):
                return 1
            if other.isBase(self.context + self.name): return -1
            '''
                           context+self.name in ['::core::chemical::PatchOperation', '::core::chemical::AddChi']:
                print '__cmp__', self.context+self.name, other.name

            for b in other.bases:
                if self.context + self.name == b.type_.T(): return -1
            for b in self.bases:
                if other.context + other.name == b.type_.T(): return 1
                '''
            return 0
        else:
            return -1
    """


#if __name__ == "__main__": main(sys.argv)
if __name__ == "__main__": GCC_XML_main(sys.argv)

'''
ToDo:
 - investigate why we can't use ::core::scoring::ScoreFunctionFactory::create_score_function ?
 - bp::implicitly_convertible< std::string const &, protocols::simple_moves::MinMover >();
 - docs for enums

'''
