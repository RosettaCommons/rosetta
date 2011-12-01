# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   BuildBuindings.py
## @brief  Build Python buidings for mini
## @author Sergey Lyskov, based on original source from Py++ 'doxygen.py' writen by G.D.
##                        and was distributed under  Boost Software License, Version 1.0.
##                        (See http://www.boost.org/LICENSE_1_0.txt)


class doxygen_doc_extractor:
    """
    Extracts Doxygen styled documentation from source or generates it from description.
    """
    def __init__(self):
        #for caching source
        self.file_name = None
        self.source = None
    #__init__

    def __call__(self, declaration):
        if type(None) == type(declaration.location): return ''  # Seriously u guys...

        if self.file_name != declaration.location.file_name:
            self.file_name = declaration.location.file_name
            self.source = open(declaration.location.file_name).readlines()

        _Code, _Empty, _LineComment, _MultiLineComment, _Terminate = range(5)

        state = _Code
        doc_lines = []
        for lcount in xrange(declaration.location.line-2, -1, -1):
            line = self.source[lcount]
            add_current_line = False

            if state == _Code or state == _Empty:
                #print line.rstrip()[-2:]
                if line.rstrip()[-2:] == "*/":
                    state = _MultiLineComment
                    add_current_line = True

                if line.lstrip()[:2] == "//":
                    state = _LineComment
                    add_current_line = True

            if state == _MultiLineComment:
                if line.lstrip()[:2] == "/*":
                    state = _Terminate
                add_current_line = True

            if state == _LineComment:
                if line.lstrip()[:2] == "//":
                    add_current_line = True
                else:
                    state = _Terminate
                    add_current_line = False

            if state == _Code:
                if self.is_empty(line):
                    state = _Empty
                    add_current_line = False
            elif state == _Empty:
                if self.is_code(line):
                    state = _Terminate
                    add_current_line = False


                #if self.is_code(line):
                    #state == _Terminate
                    #add_current_line = False

            if add_current_line:
                doc_lines.insert(0, self.clear_str(line) )

            if state == _Terminate: break

            '''
            final_str = self.clear_str(line)
            if not find_block_end and self.is_code(line) and lcount<declaration.location.line-2:
                break

            if (not find_block_end) and doc_lines: break

            if final_str:
                if lcount == declaration.location.line-2  and  (not find_block_end):
                    pass
                else:
                    doc_lines.insert(0, final_str)
                    '''
        #except:
        #    pass
        #finally:
        if doc_lines:
            final_doc_lines = [ line.replace("\n","\\n") for line in doc_lines[:-1] ]
            final_doc_lines.append(doc_lines[-1].replace("\n",""))
            return '\"' + ''.join(final_doc_lines) + '\"'
        else:
            return '"%s:%s"' % (declaration.location.file_name.split('/')[-1], declaration.location.line)
            #return '\"\"'
    #__call__()


    def clear_str(self, tmp_str):
        """
        Replace */! by space and \brief, @fn, \param, etc
        """
        clean = lambda tstr, sym, change2 = '': tstr.replace(sym, change2)

        tmp_str = reduce(clean, [tmp_str, '/','*','!',"\\brief","@brief","\\fn","@fn","\\ref","@ref", "\"", "\'", "\\c"])

        #commands list taken form : http://www.stack.nl/~dimitri/doxygen/commands.html
        replacement_list = [
#			"a",
            "addindex",
            "addtogroup",
            "anchor",
            "arg",
            "attention",
            "author",
#			"b",
#			"brief",
            "bug",
#			"c",
            "callgraph",
            "callergraph",
            "category",
            "class",
            ("code","[Code]"),
            "cond",
            "copybrief",
            "copydetails",
            "copydoc",
            "date",
            "def",
            "defgroup",
            "deprecated",
            "details",
            "dir",
            "dontinclude",
            ("dot","[Dot]"),
            "dotfile",
            "e",
            "else",
            "elseif",
            "em",
            ("endcode","[/Code]"),
            "endcond",
            ("enddot","[/Dot]"),
            "endhtmlonly",
            "endif",
            "endlatexonly",
            "endlink",
            "endmanonly",
            "endmsc",
            "endverbatim",
            "endxmlonly",
            "enum",
            "example",
            "exception",
            "extends",
            "f$",
            "f[",
            "f]",
            "f{",
            "f}",
            "file",
#			"fn",
            "headerfile",
            "hideinitializer",
            "htmlinclude",
            "htmlonly",
            "if",
            "ifnot",
            "image",
            "implements",
            "include",
            "includelineno",
            "ingroup",
            "internal",
            "invariant",
            "interface",
            "latexonly",
            "li",
            "line",
            "link",
            "mainpage",
            "manonly",
            "memberof",
            "msc",
#			"n",
            "name",
            "namespace",
            "nosubgrouping",
            "note",
            "overload",
#			"p",
            "package",
            "page",
            "par",
            "paragraph",
            "param",
            "post",
            "pre",
#			"private (PHP only)",
#			"privatesection (PHP only)",
            "property",
#			"protected (PHP only)",
#			"protectedsection (PHP only)",
            "protocol",
#			"public (PHP only)",
#			"publicsection (PHP only)",
#			"ref",
            "relates",
            "relatesalso",
            "remarks",
            "return",
            "retval",
            "sa",
            "section",
            "see",
            "showinitializer",
            "since",
            "skip",
            "skipline",
            "struct",
            "subpage",
            "subsection",
            "subsubsection",
            "test",
            "throw",
            ("todo","TODO"),
            "tparam",
            "typedef",
            "union",
            "until",
            "var",
            "verbatim",
            "verbinclude",
            "version",
            "warning",
            "weakgroup",
            "xmlonly",
            "xrefitem",
#			"$",
#			"@",
#			"\",
#			"&",
#			"~",
#			"<",
#			">",
#			"#",
#			"%",
            ]

        for command in replacement_list:
            try:
                old,new = command
            except ValueError:
                old = command
                new = command.capitalize()+":"
            tmp_str = clean(tmp_str, "@"+old, new)
            tmp_str = clean(tmp_str, "\\"+old, new)

        return tmp_str.lstrip()
    #clean_str()

    def is_code(self, tmp_str):
        """
        Detects if tmp_str is code or not
        """
        try:
            st = tmp_str.lstrip()
            beg = st[:2]
            return st!= '' and beg != "//" and beg != "/*"
        except:
            pass
        return False


    def is_empty(self, tmp_str):
        """
        Detects if tmp_str is empty or not
        """
        return tmp_str.lstrip() == ''


    #is_code()

#class doxygen_doc_extractor
