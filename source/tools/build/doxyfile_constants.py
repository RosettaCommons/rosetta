# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
#
# Doxygen constants.  In separate file to avoid clutter in doxygen.py
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# Each value in the defaults is a tuple made up of
# - The default value
# - The value's type
# - The documentation string

defaults = dict(

    # Project related configuration options
    PROJECT_NAME                = "",
    PROJECT_NUMBER              = "",
    OUTPUT_DIRECTORY            = "",
    OUTPUT_LANGUAGE             = "English",
    USE_WINDOWS_ENCODING        = "NO",
    BRIEF_MEMBER_DESC           = "YES",
    REPEAT_BRIEF                = "YES",
    ABBREVIATE_BRIEF            = "",
    ALWAYS_DETAILED_SEC         = "NO",
    INLINE_INHERITED_MEMB       = "NO",
    FULL_PATH_NAMES             = "NO",
    STRIP_FROM_PATH             = "",
    SHORT_NAMES                 = "NO",
    JAVADOC_AUTOBRIEF           = "NO",
    MULTILINE_CPP_IS_BRIEF      = "NO",
    DETAILS_AT_TOP              = "NO",
    INHERIT_DOCS                = "YES",
    DISTRIBUTE_GROUP_DOC        = "NO",
    TAB_SIZE                    = 8,
    ALIASES                     = "",
    OPTIMIZE_OUTPUT_FOR_C       = "NO",
    OPTIMIZE_OUTPUT_JAVA        = "NO",
    SUBGROUPING                 = "YES",

    # Build related configuration options
    EXTRACT_ALL                 = "NO",
    EXTRACT_PRIVATE             = "NO",
    EXTRACT_STATIC              = "NO",
    EXTRACT_LOCAL_CLASSES       = "YES",
    HIDE_UNDOC_MEMBERS          = "NO",
    HIDE_UNDOC_CLASSES          = "NO",
    HIDE_FRIEND_COMPOUNDS       = "NO",
    HIDE_IN_BODY_DOCS           = "NO",
    INTERNAL_DOCS               = "NO",
    CASE_SENSE_NAMES            = "YES",
    HIDE_SCOPE_NAMES            = "NO",
    SHOW_INCLUDE_FILES          = "YES",
    INLINE_INFO                 = "YES",
    SORT_MEMBER_DOCS            = "YES",
    SORT_BRIEF_DOCS             = "NO",
    SORT_BY_SCOPE_NAME          = "NO",
    GENERATE_TODOLIST           = "YES",
    GENERATE_TESTLIST           = "YES",
    GENERATE_BUGLIST            = "YES",
    GENERATE_DEPRECATEDLIST     = "YES",
    ENABLED_SECTIONS            = "",
    MAX_INITIALIZER_LINES       = 30,
    SHOW_USED_FILES             = "YES",

    # Configuration options related to warning and progress messages
    QUIET                       = "NO",
    WARNINGS                    = "YES",
    WARN_IF_UNDOCUMENTED        = "YES",
    WARN_IF_DOC_ERROR           = "YES",
    WARN_FORMAT                 = '"$file:$line: $text"',
    WARN_LOGFILE                = "",

    # Configuration options related to the input files
    INPUT                       = "",
    FILE_PATTERNS               = "",
    RECURSIVE                   = "NO",
    EXCLUDE                     = "",
    EXCLUDE_SYMLINKS            = "NO",
    EXCLUDE_PATTERNS            = "",
    EXAMPLE_PATH                = "",
    EXAMPLE_PATTERNS            = "",
    EXAMPLE_RECURSIVE           = "NO",
    IMAGE_PATH                  = "",
    INPUT_FILTER                = "",
    FILTER_SOURCE_FILES         = "NO",

    # Configuration options related to source browsing
    SOURCE_BROWSER              = "NO",
    INLINE_SOURCES              = "NO",
    STRIP_CODE_COMMENTS         = "YES",
    REFERENCED_BY_RELATION      = "YES",
    REFERENCES_RELATION         = "YES",
    VERBATIM_HEADERS            = "YES",

    # Configuration options related to the alphabetical class index
    ALPHABETICAL_INDEX          = "NO",
    COLS_IN_ALPHA_INDEX         = 5,
    IGNORE_PREFIX               = "",

    # Configuration options related to the HTML output
    GENERATE_HTML               = "YES",
    HTML_OUTPUT                 = "html",
    HTML_FILE_EXTENSION         = ".html",
    HTML_HEADER                 = "",
    HTML_FOOTER                 = "",
    HTML_STYLESHEET             = "",
    HTML_ALIGN_MEMBERS          = "YES",
    GENERATE_HTMLHELP           = "NO",
    CHM_FILE                    = "",
    HHC_LOCATION                = "",
    GENERATE_CHI                = "NO",
    BINARY_TOC                  = "NO",
    TOC_EXPAND                  = "NO",
    DISABLE_INDEX               = "NO",
    ENUM_VALUES_PER_LINE        = 4,
    GENERATE_TREEVIEW           = "NO",
    TREEVIEW_WIDTH              = 250,

    # Configuration options related to the LaTeX output
    GENERATE_LATEX              = "NO",
    LATEX_OUTPUT                = "latex",
    LATEX_CMD_NAME              = "latex",
    MAKEINDEX_CMD_NAME          = "makeindex",
    COMPACT_LATEX               = "NO",
    PAPER_TYPE                  = "a4wide",
    EXTRA_PACKAGES              = "",
    LATEX_HEADER                = "",
    PDF_HYPERLINKS              = "NO",
    USE_PDFLATEX                = "NO",
    LATEX_BATCHMODE             = "NO",
    LATEX_HIDE_INDICES          = "NO",

    # Configuration options related to the RTF output
    GENERATE_RTF                = "NO",
    RTF_OUTPUT                  = "rtf",
    COMPACT_RTF                 = "NO",
    RTF_HYPERLINKS              = "NO",
    RTF_STYLESHEET_FILE         = "",
    RTF_EXTENSIONS_FILE         = "",

    # Configuration options related to the man page output
    GENERATE_MAN                = "NO",
    MAN_OUTPUT                  = "man",
    MAN_EXTENSION               = ".3",
    MAN_LINKS                   = "NO",
    GENERATE_XML                = "NO",

    # Configuration options related to the XML output
    XML_OUTPUT                  = "xml",
    XML_SCHEMA                  = "",
    XML_DTD                     = "",
    XML_PROGRAMLISTING          = "YES",

    # Configuration options for the AutoGen Definitions output
    GENERATE_AUTOGEN_DEF        = "NO",

    # Configuration options related to the Perl module output
    GENERATE_PERLMOD            = "NO",
    PERLMOD_LATEX               = "NO",
    PERLMOD_PRETTY              = "YES",
    PERLMOD_MAKEVAR_PREFIX      = "",

    # Configuration options related to the preprocessor
    ENABLE_PREPROCESSING        = "YES",
    MACRO_EXPANSION             = "NO",
    EXPAND_ONLY_PREDEF          = "NO",
    SEARCH_INCLUDES             = "YES",
    INCLUDE_PATH                = "",
    INCLUDE_FILE_PATTERNS       = "",
    PREDEFINED                  = "",
    EXPAND_AS_DEFINED           = "",
    SKIP_FUNCTION_MACROS        = "YES",

    # Configuration options related to external references
    TAGFILES                    = "",
    GENERATE_TAGFILE            = "",
    ALLEXTERNALS                = "NO",
    EXTERNAL_GROUPS             = "YES",
    PERL_PATH                   = "/usr/bin/perl",

    # Configuration options related to the 'dot' graphing tool.
    # Everything in this section below HAVE_DOT depends on HAVE_DOT and
    # will default to "NO" if HAVE_DOT is "NO".  The options are set to
    # what they should be if dot is present.
    # By default we build WITHOUT dot.  If the SConscript.doc detects that
    # dot exists on the system it will change these options to use it.
    CLASS_DIAGRAMS              = "YES",
    HIDE_UNDOC_RELATIONS        = "YES",
    HAVE_DOT                    = "NO",
    CLASS_GRAPH                 = "YES",
    COLLABORATION_GRAPH         = "YES",
    UML_LOOK                    = "NO",
    TEMPLATE_RELATIONS          = "NO",
    INCLUDE_GRAPH               = "YES",
    INCLUDED_BY_GRAPH           = "YES",
    CALL_GRAPH                  = "NO",
    GRAPHICAL_HIERARCHY         = "YES",
    DOT_IMAGE_FORMAT            = "png",
    DOT_PATH                    = "",
    DOTFILE_DIRS                = "",
    MAX_DOT_GRAPH_WIDTH         = 1024,
    MAX_DOT_GRAPH_HEIGHT        = 1024,
    MAX_DOT_GRAPH_DEPTH         = 0,
    GENERATE_LEGEND             = "YES",
    DOT_CLEANUP                 = "YES",

    # Configuration options related to the search engine
    SEARCHENGINE                = "NO",
)


comments = dict(

    # ---------------------------------------------------------------------------
    # Project related configuration options
    # ---------------------------------------------------------------------------

    PROJECT_NAME = \
"""The PROJECT_NAME tag is a single word (or a sequence of words surrounded by
quotes) that should identify the project.""",

    PROJECT_NUMBER = \
"""The PROJECT_NUMBER tag can be used to enter a project or revision number.
This could be handy for archiving the generated documentation or if some version
control system is used.""",

    OUTPUT_DIRECTORY = \
"""The OUTPUT_DIRECTORY tag is used to specify the (relative or absolute) base
path where the generated documentation will be put.  If a relative path is
entered, it will be relative to the location where Doxygen was started.  If left
blank the current directory will be used.""",

    OUTPUT_LANGUAGE = \
"""The OUTPUT_LANGUAGE tag is used to specify the language in which all
documentation generated by doxygen is written. Doxygen will use this information
to generate all constant output in the proper language.

The default language is English, other supported languages are: Brazilian,
Catalan, Chinese, Chinese-Traditional, Croatian, Czech, Danish, Dutch, Finnish,
French, German, Greek, Hungarian, Italian, Japanese, Japanese-en (Japanese with
English messages), Korean, Korean-en, Norwegian, Polish, Portuguese, Romanian,
Russian, Serbian, Slovak, Slovene, Spanish, Swedish, and Ukrainian.""",

    USE_WINDOWS_ENCODING = \
"""This tag can be used to specify the encoding used in the generated output.
The encoding is not always determined by the language that is chosen, but also
whether or not the output is meant for Windows or non-Windows users.  In case
there is a difference, setting the USE_WINDOWS_ENCODING tag to YES forces the
Windows encoding (this is the default for the Windows binary), whereas setting
the tag to NO uses a Unix-style encoding (the default for all platforms other
than Windows).""",

    BRIEF_MEMBER_DESC = \
"""If the BRIEF_MEMBER_DESC tag is set to YES (the default) Doxygen will include
brief member descriptions after the members that are listed in the file and
class documentation (similar to JavaDoc).  Set to NO to disable this.""",

    REPEAT_BRIEF = \
"""If the REPEAT_BRIEF tag is set to YES (the default) Doxygen will prepend the
brief description of a member or function before the detailed description. Note:
if both HIDE_UNDOC_MEMBERS and BRIEF_MEMBER_DESC are set to NO, the brief
descriptions will be completely suppressed.""",

    ABBREVIATE_BRIEF = \
'''This tag implements a quasi-intelligent brief description abbreviator that is
used to form the text in various listings. Each string in this list, if found as
the leading text of the brief description, will be stripped from the text and
the result after processing the whole list, is used as the annotated
text. Otherwise, the brief description is used as-is. If left blank, the
following values are used ("$name" is automatically replaced with the name of
the entity): "The $name class" "The $name widget" "The $name file" "is"
"provides" "specifies" "contains" "represents" "a" "an" "the"''',

    ALWAYS_DETAILED_SEC = \
"""If the ALWAYS_DETAILED_SEC and REPEAT_BRIEF tags are both set to YES then
Doxygen will generate a detailed section even if there is only a brief
description.""",

    INLINE_INHERITED_MEMB = \
"""If the INLINE_INHERITED_MEMB tag is set to YES, doxygen will show all
inherited members of a class in the documentation of that class as if those
members were ordinary class members. Constructors, destructors and assignment
operators of the base classes will not be shown.""",

    FULL_PATH_NAMES = \
"""If the FULL_PATH_NAMES tag is set to YES then Doxygen will prepend the full
path before files name in the file list and in the header files. If set to NO
the shortest path that makes the file name unique will be used.""",

    STRIP_FROM_PATH = \
"""If the FULL_PATH_NAMES tag is set to YES then the STRIP_FROM_PATH tag can be
used to strip a user-defined part of the path. Stripping is only done if one of
the specified strings matches the left-hand part of the path. It is allowed to
use relative paths in the argument list. If left blank the directory from which
doxygen is run is used as the path to strip.""",

    SHORT_NAMES = \
"""If the SHORT_NAMES tag is set to YES, doxygen will generate much shorter (but
less readable) file names. This can be useful is your file systems doesn't
support long names like on DOS, Mac, or CD-ROM.""",

    JAVADOC_AUTOBRIEF = \
"""If the JAVADOC_AUTOBRIEF tag is set to YES then Doxygen will interpret the
first line (until the first dot) of a JavaDoc-style comment as the brief
description. If set to NO, the JavaDoc comments will behave just like the
Qt-style comments (thus requiring an explicit @brief command for a brief
description.""",

    MULTILINE_CPP_IS_BRIEF = \
"""The MULTILINE_CPP_IS_BRIEF tag can be set to YES to make Doxygen treat a
multi-line C++ special comment block (i.e. a block of //! or /// comments) as a
brief description. This used to be the default behaviour. The new default is to
treat a multi-line C++ comment block as a detailed description. Set this tag to
YES if you prefer the old behaviour instead.""",

    DETAILS_AT_TOP = \
"""If the DETAILS_AT_TOP tag is set to YES then Doxygen will output the detailed
description near the top, like JavaDoc.  If set to NO, the detailed description
appears after the member documentation.""",

    INHERIT_DOCS = \
"""If the INHERIT_DOCS tag is set to YES (the default) then an undocumented
member inherits the documentation from any documented member that it
re-implements.""",

    DISTRIBUTE_GROUP_DOC = \
"""If member grouping is used in the documentation and the DISTRIBUTE_GROUP_DOC
tag is set to YES, then doxygen will reuse the documentation of the first member
in the group (if any) for the other members of the group. By default all members
of a group must be documented explicitly.""",

    TAB_SIZE = \
 """The TAB_SIZE tag can be used to set the number of spaces in a tab. Doxygen
 uses this value to replace tabs by spaces in code fragments.""",

    ALIASES = \
"""This tag can be used to specify a number of aliases that acts as commands in
the documentation. An alias has the form "name=value". For example adding
"sideeffect=\par Side Effects:\n" will allow you to put the command \sideeffect
(or @sideeffect) in the documentation, which will result in a user-defined
paragraph with heading "Side Effects:". You can put \n's in the value part of an
alias to insert newlines.""",

    OPTIMIZE_OUTPUT_FOR_C = \
"""Set the OPTIMIZE_OUTPUT_FOR_C tag to YES if your project consists of C
sources only. Doxygen will then generate output that is more tailored for C. For
instance, some of the names that are used will be different. The list of all
members will be omitted, etc.""",

    OPTIMIZE_OUTPUT_JAVA = \
"""Set the OPTIMIZE_OUTPUT_JAVA tag to YES if your project consists of Java
sources only. Doxygen will then generate output that is more tailored for
Java. For instance, namespaces will be presented as packages, qualified scopes
will look different, etc.""",

    SUBGROUPING = \
"""Set the SUBGROUPING tag to YES (the default) to allow class member groups of
the same type (for instance a group of public functions) to be put as a subgroup
of that type (e.g. under the Public Functions section). Set it to NO to prevent
subgrouping. Alternatively, this can be done per class using the \nosubgrouping
command.""",

    # ---------------------------------------------------------------------------
    # Build related configuration options
    # ---------------------------------------------------------------------------

    EXTRACT_ALL = \
"""If the EXTRACT_ALL tag is set to YES doxygen will assume all entities in
documentation are documented, even if no documentation was available. Private
class members and static file members will be hidden unless the EXTRACT_PRIVATE
and EXTRACT_STATIC tags are set to YES""",

    EXTRACT_PRIVATE = \
"""If the EXTRACT_PRIVATE tag is set to YES all private members of a class will
be included in the documentation.""",

    EXTRACT_STATIC = \
"""If the EXTRACT_STATIC tag is set to YES all static members of a file will be
included in the documentation.""",

    EXTRACT_LOCAL_CLASSES = \
        """If the EXTRACT_LOCAL_CLASSES tag is set to YES classes (and structs)
defined locally in source files will be included in the documentation.
If set to NO only classes defined in header files are included.""",

    HIDE_UNDOC_MEMBERS = \
"""If the HIDE_UNDOC_MEMBERS tag is set to YES, Doxygen will hide all
undocumented members of documented classes, files or namespaces. If set to NO
(the default) these members will be included in the various overviews, but no
documentation section is generated. This option has no effect if EXTRACT_ALL is
enabled.""",

    HIDE_UNDOC_CLASSES = \
"""If the HIDE_UNDOC_CLASSES tag is set to YES, Doxygen will hide all
undocumented classes that are normally visible in the class hierarchy. If set to
NO (the default) these classes will be included in the various overviews. This
option has no effect if EXTRACT_ALL is enabled.""",

    HIDE_FRIEND_COMPOUNDS = \
"""If the HIDE_FRIEND_COMPOUNDS tag is set to YES, Doxygen will hide all friend
(class|struct|union) declarations. If set to NO (the default) these declarations
will be included in the documentation.""",

    HIDE_IN_BODY_DOCS = \
"""If the HIDE_IN_BODY_DOCS tag is set to YES, Doxygen will hide any
documentation blocks found inside the body of a function. If set to NO (the
default) these blocks will be appended to the
function's detailed documentation block.""",

    INTERNAL_DOCS = \
"""The INTERNAL_DOCS tag determines if documentation that is typed after a
\internal command is included. If the tag is set to NO (the default) then the
documentation will be excluded.
Set it to YES to include the internal documentation.""",

    CASE_SENSE_NAMES = \
"""If the CASE_SENSE_NAMES tag is set to NO then Doxygen will only generate file
names in lower-case letters. If set to YES upper-case letters are also
allowed. This is useful if you have classes or files whose names only differ
in case and if your file system supports case sensitive file names. Windows
users are advised to set this option to NO.""",

    HIDE_SCOPE_NAMES = \
"""If the HIDE_SCOPE_NAMES tag is set to NO (the default) then Doxygen will show
members with their full class and namespace scopes in the documentation. If set
to YES the scope will be hidden.""",

    SHOW_INCLUDE_FILES = \
"""If the SHOW_INCLUDE_FILES tag is set to YES (the default) then Doxygen will
put a list of the files that are included by a file in the documentation of that
file.""",

    INLINE_INFO = \
"""If the INLINE_INFO tag is set to YES (the default) then a tag [inline] is
inserted in the documentation for inline members.""",

    SORT_MEMBER_DOCS = \
"""If the SORT_MEMBER_DOCS tag is set to YES (the default) then doxygen will
sort the (detailed) documentation of file and class members alphabetically by
member name. If set to NO the members will appear in
declaration order.""",

    SORT_BRIEF_DOCS = \
"""If the SORT_BRIEF_DOCS tag is set to YES then doxygen will sort the brief
documentation of file, namespace and class members alphabetically by member
name. If set to NO (the default) the members will appear in
declaration order.""",

    SORT_BY_SCOPE_NAME = \
"""If the SORT_BY_SCOPE_NAME tag is set to YES, the class list will be sorted by
fully-qualified names, including namespaces. If set to NO (the default), the
class list will be sorted only by class name,
not including the namespace part.
Note: This option is not very useful if HIDE_SCOPE_NAMES is set to YES.
Note: This option applies only to the class list, not to the
alphabetical list.""",

    GENERATE_TODOLIST = \
"""The GENERATE_TODOLIST tag can be used to enable (YES) or disable (NO) the
todo list. This list is created by putting \todo commands in the
documentation.""",

    GENERATE_TESTLIST = \
"""The GENERATE_TESTLIST tag can be used to enable (YES) or disable (NO) the
test list. This list is created by putting \test commands in the
documentation.""",

    GENERATE_BUGLIST = \
"""The GENERATE_BUGLIST tag can be used to enable (YES) or disable (NO) the bug
list. This list is created by putting \bug commands in the documentation.""",

    GENERATE_DEPRECATEDLIST = \
"""The GENERATE_DEPRECATEDLIST tag can be used to enable (YES) or disable (NO)
the deprecated list. This list is created by putting \deprecated commands in the
documentation.""",

    ENABLED_SECTIONS = \
"""The ENABLED_SECTIONS tag can be used to enable conditional documentation
sections, marked by \if sectionname ... \endif.""",

    MAX_INITIALIZER_LINES = \
"""The MAX_INITIALIZER_LINES tag determines the maximum number of lines the
initial value of a variable or define consists of for it to appear in the
documentation. If the initializer consists of more lines than specified
here it will be hidden. Use a value of 0 to hide initializers completely.
The appearance of the initializer of individual variables and defines in the
documentation can be controlled using \showinitializer or \hideinitializer
command in the documentation regardless of this setting.""",

    SHOW_USED_FILES = \
"""Set the SHOW_USED_FILES tag to NO to disable the list of files generated at
the bottom of the documentation of classes and structs. If set to YES the list
will mention the files that were used to generate the documentation.""",

    # ---------------------------------------------------------------------------
    # Configuration options related to warning and progress messages
    # ---------------------------------------------------------------------------

    QUIET = \
"""The QUIET tag can be used to turn on/off the messages that are generated by
doxygen. Possible values are YES and NO. If left blank NO is used.""",

    WARNINGS = \
"""The WARNINGS tag can be used to turn on/off the warning messages that are
generated by doxygen. Possible values are YES and NO. If left blank NO is
used.""",

    WARN_IF_UNDOCUMENTED = \
"""If WARN_IF_UNDOCUMENTED is set to YES, then doxygen will generate warnings
for undocumented members. If EXTRACT_ALL is set to YES then this flag will
automatically be disabled.""",

    WARN_IF_DOC_ERROR = \
"""If WARN_IF_DOC_ERROR is set to YES, doxygen will generate warnings for
potential errors in the documentation, such as not documenting some parameters
in a documented function, or documenting parameters that
don't exist or using markup commands wrongly.""",

    WARN_FORMAT = \
"""The WARN_FORMAT tag determines the format of the warning messages that
doxygen can produce. The string should contain the $file, $line, and $text tags,
which will be replaced by the file and line number from which the
warning originated and the warning text.""",

    WARN_LOGFILE = \
"""The WARN_LOGFILE tag can be used to specify a file to which warning and error
messages should be written. If left blank the output is written to stderr.""",

        # ---------------------------------------------------------------------------
        # Configuration options related to the input files
        # ---------------------------------------------------------------------------

    INPUT = \
"""The INPUT tag can be used to specify the files and/or directories that
contain documented source files. You may enter file names like "myfile.cpp" or
directories like "/usr/src/myproject". Separate the files or directories with
spaces.""",

    FILE_PATTERNS = \
"""If the value of the INPUT tag contains directories, you can use the
FILE_PATTERNS tag to specify one or more wildcard pattern (like *.cpp and *.h)
to filter out the source-files in the directories. If left blank the following
patterns are tested: *.c *.cc *.cxx *.cpp *.c++ *.java *.ii *.ixx *.ipp *.i++
*.inl *.h *.hh *.hxx *.hpp *.h++ *.idl *.odl *.cs *.php *.php3 *.inc""",

    RECURSIVE = \
"""The RECURSIVE tag can be used to turn specify whether or not subdirectories
should be searched for input files as well. Possible values are YES and NO. If
left blank NO is used.""",

    EXCLUDE = \
"""The EXCLUDE tag can be used to specify files and/or directories that should
excluded from the INPUT source files. This way you can easily exclude a
subdirectory from a directory tree whose root is specified with the INPUT
tag.""",

    EXCLUDE_SYMLINKS = \
"""The EXCLUDE_SYMLINKS tag can be used select whether or not files or
directories that are symbolic links (a Unix filesystem feature) are excluded
from the input.""",

    EXCLUDE_PATTERNS = \
"""If the value of the INPUT tag contains directories, you can use the
EXCLUDE_PATTERNS tag to specify one or more wildcard patterns to exclude certain
files from those directories.""",

    EXAMPLE_PATH = \
"""The EXAMPLE_PATH tag can be used to specify one or more files or directories
that contain example code fragments that are included (see the \include
command).""",

    EXAMPLE_PATTERNS = \
"""If the value of the EXAMPLE_PATH tag contains directories, you can use the
EXAMPLE_PATTERNS tag to specify one or more wildcard pattern (like *.cpp and
*.h) to filter out the source-files in the directories. If left blank all files
are included.""",

    EXAMPLE_RECURSIVE = \
"""If the EXAMPLE_RECURSIVE tag is set to YES then subdirectories will be
searched for input files to be used with the \include or \dontinclude commands
irrespective of the value of the RECURSIVE tag. Possible values are YES and
NO. If left blank NO is used.""",

    IMAGE_PATH = \
"""The IMAGE_PATH tag can be used to specify one or more files or directories
that contain image that are included in the documentation (see the \image
command).""",

    INPUT_FILTER = \
"""The INPUT_FILTER tag can be used to specify a program that doxygen should
invoke to filter for each input file. Doxygen will invoke the filter program by
executing (via popen()) the command <filter> <input-file>, where <filter> is the
value of the INPUT_FILTER tag, and <input-file> is the name of an input
file. Doxygen will then use the output that the filter program writes to
standard output.""",

    FILTER_SOURCE_FILES = \
"""If the FILTER_SOURCE_FILES tag is set to YES, the input filter (if set using
INPUT_FILTER) will be used to filter the input files when producing source files
to browse (i.e. when SOURCE_BROWSER is set to YES).""",


    # ---------------------------------------------------------------------------
    # Configuration options related to source browsing
    # ---------------------------------------------------------------------------

    SOURCE_BROWSER = \
"""If the SOURCE_BROWSER tag is set to YES then a list of source files will be
generated. Documented entities will be cross-referenced with these
sources. Note: To get rid of all source code in the generated output, make sure
also VERBATIM_HEADERS is set to NO.""",

    INLINE_SOURCES = \
"""Setting the INLINE_SOURCES tag to YES will include the body of functions and
classes directly in the documentation.""",

    STRIP_CODE_COMMENTS = \
"""Setting the STRIP_CODE_COMMENTS tag to YES (the default) will instruct
doxygen to hide any special comment blocks from generated source code
fragments. Normal C and C++ comments will always remain visible.""",

    REFERENCED_BY_RELATION = \
"""If the REFERENCED_BY_RELATION tag is set to YES (the default) then for each
documented function all documented functions referencing it will be listed.""",

    REFERENCES_RELATION = \
"""If the REFERENCES_RELATION tag is set to YES (the default) then for each
documented function all documented entities called/used by that function will be
listed.""",

    VERBATIM_HEADERS = \
"""If the VERBATIM_HEADERS tag is set to YES (the default) then Doxygen will
generate a verbatim copy of the header file for each class for which an include
is specified. Set to NO to disable this.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the alphabetical class index
    # ---------------------------------------------------------------------------

    ALPHABETICAL_INDEX = \
"""If the ALPHABETICAL_INDEX tag is set to YES, an alphabetical index of all
compounds will be generated. Enable this if the project contains a lot of
classes, structs, unions or interfaces.""",

    COLS_IN_ALPHA_INDEX = \
"""If the alphabetical index is enabled (see ALPHABETICAL_INDEX) then the
COLS_IN_ALPHA_INDEX tag can be used to specify the number of columns in which
this list will be split (can be a number in the range [1..20])""",

    IGNORE_PREFIX = \
"""In case all classes in a project start with a common prefix, all classes will
be put under the same header in the alphabetical index. The IGNORE_PREFIX tag
can be used to specify one or more prefixes that should be ignored while
generating the index headers.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the HTML output
    # ---------------------------------------------------------------------------

    GENERATE_HTML = \
"""If the GENERATE_HTML tag is set to YES (the default) Doxygen will generate
HTML output.""",

    HTML_OUTPUT = \
"""The HTML_OUTPUT tag is used to specify where the HTML docs will be put. If a
relative path is entered the value of OUTPUT_DIRECTORY will be put in front of
it. If left blank `html' will be used as the default path.""",

    HTML_FILE_EXTENSION = \
"""The HTML_FILE_EXTENSION tag can be used to specify the file extension for
each generated HTML page (for example: .htm,.php,.asp). If it is left blank
doxygen will generate files with .html extension.""",

    HTML_HEADER = \
"""The HTML_HEADER tag can be used to specify a personal HTML header for each
generated HTML page. If it is left blank doxygen will generate a standard
header.""",

    HTML_FOOTER = \
"""The HTML_FOOTER tag can be used to specify a personal HTML footer for each
generated HTML page. If it is left blank doxygen will generate a standard
footer.""",

    HTML_STYLESHEET = \
"""The HTML_STYLESHEET tag can be used to specify a user-defined cascading style
sheet that is used by each HTML page. It can be used to fine-tune the look of
the HTML output. If the tag is left blank doxygen will generate a default style
sheet. Note that doxygen will try to copy the style sheet file to the HTML
output directory, so don't put your own stylesheet in the HTML output directory
as well, or it will be erased!""",

    HTML_ALIGN_MEMBERS = \
"""If the HTML_ALIGN_MEMBERS tag is set to YES, the members of classes, files or
namespaces will be aligned in HTML using tables. If set to NO a bullet list will
be used.""",

    GENERATE_HTMLHELP = \
"""If the GENERATE_HTMLHELP tag is set to YES, additional index files will be
generated that can be used as input for tools like the Microsoft HTML help
workshop to generate a compressed HTML help file (.chm) of the generated HTML
documentation.""",

    CHM_FILE = \
"""If the GENERATE_HTMLHELP tag is set to YES, the CHM_FILE tag can be used to
specify the file name of the resulting .chm file. You can add a path in front of
the file if the result should not be written to the html output directory.""",

    HHC_LOCATION = \
"""If the GENERATE_HTMLHELP tag is set to YES, the HHC_LOCATION tag can be used
to specify the location (absolute path including file name) of the HTML help
compiler (hhc.exe). If non-empty doxygen will try to run the HTML help compiler
on the generated index.hhp.""",

    GENERATE_CHI = \
"""If the GENERATE_HTMLHELP tag is set to YES, the GENERATE_CHI flag controls if
a separate .chi index file is generated (YES) or that it should be included in
the master .chm file (NO).""",

    BINARY_TOC = \
"""If the GENERATE_HTMLHELP tag is set to YES, the BINARY_TOC flag controls
whether a binary table of contents is generated (YES) or a normal table of
contents (NO) in the .chm file.""",

    TOC_EXPAND = \
"""The TOC_EXPAND flag can be set to YES to add extra items for group members to
the contents of the HTML help documentation and to the tree view.""",

    DISABLE_INDEX = \
"""The DISABLE_INDEX tag can be used to turn on/off the condensed index at top
of each HTML page. The value NO (the default) enables the index and the value
YES disables it.""",

    ENUM_VALUES_PER_LINE = \
"""This tag can be used to set the number of enum values (range [1..20]) that
doxygen will group on one line in the generated HTML documentation.""",

    GENERATE_TREEVIEW = \
"""If the GENERATE_TREEVIEW tag is set to YES, a side panel will be generated
containing a tree-like index structure (just like the one that is generated for
HTML Help). For this to work a browser that supports JavaScript, DHTML, CSS and
frames is required (for instance Mozilla 1.0+, Netscape 6.0+, Internet explorer
5.0+, or Konqueror). Windows users are probably better off using the HTML help
feature.""",

    TREEVIEW_WIDTH = \
"""If the treeview is enabled (see GENERATE_TREEVIEW) then this tag can be used
to set the initial width (in pixels) of the frame in which the tree is shown.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the LaTeX output
    # ---------------------------------------------------------------------------

    GENERATE_LATEX = \
"""If the GENERATE_LATEX tag is set to YES (the default) Doxygen will generate
Latex output.""",

    LATEX_OUTPUT = \
"""The LATEX_OUTPUT tag is used to specify where the LaTeX docs will be put. If
a relative path is entered the value of OUTPUT_DIRECTORY will be put in front of
it. If left blank `latex' will be used as the default path.""",

    LATEX_CMD_NAME = \
"""The LATEX_CMD_NAME tag can be used to specify the LaTeX command name to be
invoked. If left blank `latex' will be used as the default command name.""",

    MAKEINDEX_CMD_NAME = \
"""The MAKEINDEX_CMD_NAME tag can be used to specify the command name to
generate index for LaTeX. If left blank `makeindex' will be used as the default
command name.""",

    COMPACT_LATEX = \
"""If the COMPACT_LATEX tag is set to YES Doxygen generates more compact LaTeX
documents. This may be useful for small projects and may help to save some trees
in general.""",

    PAPER_TYPE = \
"""The PAPER_TYPE tag can be used to set the paper type that is used by the
printer. Possible values are: a4, a4wide, letter, legal and executive. If left
blank a4wide will be used.""",

    EXTRA_PACKAGES = \
"""The EXTRA_PACKAGES tag can be to specify one or more names of LaTeX packages
that should be included in the LaTeX output.""",

    LATEX_HEADER = \
"""The LATEX_HEADER tag can be used to specify a personal LaTeX header for the
generated latex document. The header should contain everything until the first
chapter. If it is left blank doxygen will generate a standard header. Notice:
only use this tag if you know what you are doing!""",

    PDF_HYPERLINKS = \
"""If the PDF_HYPERLINKS tag is set to YES, the LaTeX that is generated is
prepared for conversion to pdf (using ps2pdf). The pdf file will contain links
(just like the HTML output) instead of page references This makes the output
suitable for online browsing using a pdf viewer.""",

    USE_PDFLATEX = \
"""If the USE_PDFLATEX tag is set to YES, pdflatex will be used instead of plain
latex in the generated Makefile. Set this option to YES to get a higher quality
PDF documentation.""",

    LATEX_BATCHMODE = \
"""If the LATEX_BATCHMODE tag is set to YES, doxygen will add the
\\batchmode. command to the generated LaTeX files. This will instruct LaTeX to
keep running if errors occur, instead of asking the user for help. This option
is also used when generating formulas in HTML.""",

    LATEX_HIDE_INDICES = \
"""If LATEX_HIDE_INDICES is set to YES then doxygen will not include the index
chapters (such as File Index, Compound Index, etc.) in the output.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the RTF output
    # ---------------------------------------------------------------------------

    GENERATE_RTF = \
"""If the GENERATE_RTF tag is set to YES Doxygen will generate RTF output The
RTF output is optimized for Word 97 and may not look very pretty with other RTF
readers or editors.""",

    RTF_OUTPUT = \
"""The RTF_OUTPUT tag is used to specify where the RTF docs will be put. If a
relative path is entered the value of OUTPUT_DIRECTORY will be put in front of
it. If left blank `rtf' will be used as the default path.""",

    COMPACT_RTF = \
"""If the COMPACT_RTF tag is set to YES Doxygen generates more compact RTF
documents. This may be useful for small projects and may help to save some trees
in general.""",

    RTF_HYPERLINKS = \
"""If the RTF_HYPERLINKS tag is set to YES, the RTF that is generated will
contain hyperlink fields. The RTF file will contain links (just like the HTML
output) instead of page references. This makes the output suitable for online
browsing using WORD or other programs which support those fields. Note: wordpad
(write) and others do not support links.""",

    RTF_STYLESHEET_FILE = \
"""Load stylesheet definitions from file. Syntax is similar to doxygen's config
file, i.e. a series of assignments. You only have to provide replacements,
missing definitions are set to their default value.""",

    RTF_EXTENSIONS_FILE = \
"""Set optional variables used in the generation of an rtf document. Syntax is
similar to doxygen's config file.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the man page output
    # ---------------------------------------------------------------------------

    GENERATE_MAN = \
"""If the GENERATE_MAN tag is set to YES (the default) Doxygen will generate man
pages""",

    MAN_OUTPUT = \
"""The MAN_OUTPUT tag is used to specify where the man pages will be put. If a
relative path is entered the value of OUTPUT_DIRECTORY will be put in front of
it. If left blank `man' will be used as the default path.""",

    MAN_EXTENSION = \
"""The MAN_EXTENSION tag determines the extension that is added to the generated
man pages (default is the subroutine's section .3)""",

    MAN_LINKS = \
"""If the MAN_LINKS tag is set to YES and Doxygen generates man output, then it
will generate one additional man file for each entity documented in the real man
page(s). These additional files only source the real man page, but without them
the man command would be unable to find the correct page. The default is NO.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the XML output
    # ---------------------------------------------------------------------------

    GENERATE_XML = \
"""If the GENERATE_XML tag is set to YES Doxygen will generate an XML file that
captures the structure of the code including all documentation.""",

    XML_OUTPUT = \
"""The XML_OUTPUT tag is used to specify where the XML pages will be put. If a
relative path is entered the value of OUTPUT_DIRECTORY will be put in front of
it. If left blank `xml' will be used as the default path.""",

    XML_SCHEMA = \
"""The XML_SCHEMA tag can be used to specify an XML schema, which can be used by
a validating XML parser to check the syntax of the XML files.""",

    XML_DTD = \
"""The XML_DTD tag can be used to specify an XML DTD, which can be used by a
validating XML parser to check the syntax of the XML files.""",

    XML_PROGRAMLISTING = \
"""If the XML_PROGRAMLISTING tag is set to YES Doxygen will dump the program
listings (including syntax highlighting and cross-referencing information) to
the XML output. Note that enabling this will significantly increase the size of
the XML output.""",


    # ---------------------------------------------------------------------------
    # Configuration options for the AutoGen Definitions output
    # ---------------------------------------------------------------------------

    GENERATE_AUTOGEN_DEF = \
"""If the GENERATE_AUTOGEN_DEF tag is set to YES Doxygen will generate an
AutoGen Definitions (see autogen.sf.net) file that captures the structure of the
code including all documentation. Note that this feature is still experimental
and incomplete at the moment.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the Perl module output
    # ---------------------------------------------------------------------------

    GENERATE_PERLMOD = \
"""If the GENERATE_PERLMOD tag is set to YES Doxygen will generate a Perl module
file that captures the structure of the code including all documentation. Note
that this feature is still experimental and incomplete at the moment.""",

    PERLMOD_LATEX = \
"""If the PERLMOD_LATEX tag is set to YES Doxygen will generate the necessary
Makefile rules, Perl scripts and LaTeX code to be able to generate PDF and DVI
output from the Perl module output.""",

    PERLMOD_PRETTY = \
"""If the PERLMOD_PRETTY tag is set to YES the Perl module output will be nicely
formatted so it can be parsed by a human reader.  This is useful if you want to
understand what is going on.  On the other hand, if this tag is set to NO the
size of the Perl module output will be much smaller and Perl will parse it just
the same.""",

    PERLMOD_MAKEVAR_PREFIX = \
"""The names of the make variables in the generated doxyrules.make file are
prefixed with the string contained in PERLMOD_MAKEVAR_PREFIX. This is useful so
different doxyrules.make files included by the same Makefile don't overwrite
each other's variables.""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the preprocessor
    # ---------------------------------------------------------------------------
    ENABLE_PREPROCESSING = \
"""If the ENABLE_PREPROCESSING tag is set to YES (the default) Doxygen will
evaluate all C-preprocessor directives found in the sources and include
files.""",

    MACRO_EXPANSION = \
"""If the MACRO_EXPANSION tag is set to YES Doxygen will expand all macro names
in the source code. If set to NO (the default) only conditional compilation will
be performed. Macro expansion can be done in a controlled way by setting
EXPAND_ONLY_PREDEF to YES.""",

    EXPAND_ONLY_PREDEF = \
"""If the EXPAND_ONLY_PREDEF and MACRO_EXPANSION tags are both set to YES then
the macro expansion is limited to the macros specified with the PREDEFINED and
EXPAND_AS_PREDEFINED tags.""",

    SEARCH_INCLUDES = \
"""If the SEARCH_INCLUDES tag is set to YES (the default) the includes files in
the INCLUDE_PATH (see below) will be search if a #include is found.""",

    INCLUDE_PATH = \
"""The INCLUDE_PATH tag can be used to specify one or more directories that
contain include files that are not input files but should be processed by the
preprocessor.""",

    INCLUDE_FILE_PATTERNS = \
"""You can use the INCLUDE_FILE_PATTERNS tag to specify one or more wildcard
patterns (like *.h and *.hpp) to filter out the header-files in the
directories. If left blank, the patterns specified with FILE_PATTERNS will be
used.""",

    PREDEFINED = \
"""The PREDEFINED tag can be used to specify one or more macro names that are
defined before the preprocessor is started (similar to the -D option of
gcc). The argument of the tag is a list of macros of the form: name or
name=definition (no spaces). If the definition and the = are omitted =1 is
assumed.""",

    EXPAND_AS_DEFINED = \
"""If the MACRO_EXPANSION and EXPAND_ONLY_PREDEF tags are set to YES then this
tag can be used to specify a list of macro names that should be expanded. The
macro definition that is found in the sources will be used. Use the PREDEFINED
tag if you want to use a different macro definition.""",

    SKIP_FUNCTION_MACROS = \
"""If the SKIP_FUNCTION_MACROS tag is set to YES (the default) then doxygen's
preprocessor will remove all function-like macros that are alone on a line, have
an all uppercase name, and do not end with a semicolon. Such function macros are
typically used for boiler-plate code, and will confuse the parser if not
removed.""",


    # ---------------------------------------------------------------------------
    # Configuration::additions related to external references
    # ---------------------------------------------------------------------------

    TAGFILES = \
"""The TAGFILES option can be used to specify one or more tagfiles. Optionally
an initial location of the external documentation can be added for each
tagfile. The format of a tag file without this location is as follows: TAGFILES
= file1 file2 ... Adding location for the tag files is done as follows: TAGFILES
= file1=loc1 "file2 = loc2" ... where "loc1" and "loc2" can be relative or
absolute paths or URLs. If a location is present for each tag, the installdox
tool does not have to be run to correct the links.  Note that each tag file must
have a unique name (where the name does NOT include the path) If a tag file is
not located in the directory in which doxygen is run, you must also specify the
path to the tagfile here.""",

    GENERATE_TAGFILE = \
"""When a file name is specified after GENERATE_TAGFILE, doxygen will create a
tag file that is based on the input files it reads.""",

    ALLEXTERNALS = \
"""If the ALLEXTERNALS tag is set to YES all external classes will be listed in
the class index. If set to NO only the inherited external classes will be
listed.""",

    EXTERNAL_GROUPS = \
"""If the EXTERNAL_GROUPS tag is set to YES all external groups will be listed
in the modules index. If set to NO, only the current project's groups will be
listed.""",

    PERL_PATH = \
"""The PERL_PATH should be the absolute path and name of the perl script
interpreter (i.e. the result of `which perl').""",


    # ---------------------------------------------------------------------------
    # Configuration options related to the dot tool
    # ---------------------------------------------------------------------------


    CLASS_DIAGRAMS = \
    """If the
generate a inheritance diagram (in HTML, RTF and LaTeX) for classes with base or
super classes. Setting the tag to NO turns the diagrams off. Note that this
option is superseded by the HAVE_DOT option below. This is only a fallback. It
is recommended to install and use dot, since it yields more powerful graphs.""",

    HIDE_UNDOC_RELATIONS = \
"""If set to YES, the inheritance and collaboration graphs will hide inheritance
and usage relations if the target is undocumented or is not a class.""",

    HAVE_DOT = \
"""If you set the HAVE_DOT tag to YES then doxygen will assume the dot tool is
available from the path. This tool is part of Graphviz, a graph visualization
toolkit from AT&T and Lucent Bell Labs. The other options in this section have
no effect if this option is set to NO (the default)""",

    CLASS_GRAPH = \
"""If the CLASS_GRAPH and HAVE_DOT tags are set to YES then doxygen will
generate a graph for each documented class showing the direct and indirect
inheritance relations. Setting this tag to YES will force the the CLASS_DIAGRAMS
tag to NO.""",

    COLLABORATION_GRAPH = \
"""If the COLLABORATION_GRAPH and HAVE_DOT tags are set to YES then doxygen will
generate a graph for each documented class showing the direct and indirect
implementation dependencies (inheritance, containment, and class references
variables) of the class with other documented classes.""",

    UML_LOOK = \
"""If the UML_LOOK tag is set to YES doxygen will generate inheritance and
collaboration diagrams in a style similar to the OMG's Unified Modeling
Language.""",

    TEMPLATE_RELATIONS = \
"""If set to YES, the inheritance and collaboration graphs will show the
relations between templates and their instances.""",

    INCLUDE_GRAPH = \
"""If the ENABLE_PREPROCESSING, SEARCH_INCLUDES, INCLUDE_GRAPH, and HAVE_DOT
tags are set to YES then doxygen will generate a graph for each documented file
showing the direct and indirect include dependencies of the file with other
documented files.""",

    INCLUDED_BY_GRAPH = \
"""If the ENABLE_PREPROCESSING, SEARCH_INCLUDES, INCLUDED_BY_GRAPH, and HAVE_DOT
tags are set to YES then doxygen will generate a graph for each documented
header file showing the documented files that directly or indirectly include
this file.""",

    CALL_GRAPH = \
"""If the CALL_GRAPH and HAVE_DOT tags are set to YES then doxygen will generate
a call dependency graph for every global function or class method. Note that
enabling this option will significantly increase the time of a run. So in most
cases it will be better to enable call graphs for selected functions only using
the \callgraph command.""",

    GRAPHICAL_HIERARCHY = \
"""If the GRAPHICAL_HIERARCHY and HAVE_DOT tags are set to YES then doxygen will
graphical hierarchy of all classes instead of a textual one.""",

    DOT_IMAGE_FORMAT = \
"""The DOT_IMAGE_FORMAT tag can be used to set the image format of the images
generated by dot. Possible values are png, jpg, or gif If left blank png will be
used.""",

    DOT_PATH = \
"""The tag DOT_PATH can be used to specify the path where the dot tool can be
found. If left blank, it is assumed the dot tool can be found on the path.""",

    DOTFILE_DIRS = \
"""The DOTFILE_DIRS tag can be used to specify one or more directories that
contain dot files that are included in the documentation (see the \dotfile
command).""",

    MAX_DOT_GRAPH_WIDTH = \
"""The MAX_DOT_GRAPH_WIDTH tag can be used  to set the maximum allowed
width (in pixels) of  the graphs generated by  dot. If a graph becomes
larger than this  value,  doxygen will try  to truncate  the graph, so
that it   fits within the   specified   constraint. Beware  that  most
browsers cannot cope with very large images.""",


    MAX_DOT_GRAPH_HEIGHT = \
    """The MAX_DOT_GRAPH_HEIGHT tag can be used to set the maximum allows height
(in pixels) of the graphs generated by dot. If a graph becomes larger than
this value, doxygen will try to truncate the graph, so that it fits within
the specified constraint. Beware that most browsers cannot cope with very
large images.""",

    MAX_DOT_GRAPH_DEPTH = \

"""The MAX_DOT_GRAPH_DEPTH  tag  can be  used  to set the  maximum depth  of the
graphs generated by dot. A depth value of 3 means that only nodes reachable from
the root by following a path via  at most 3 edges will  be shown. Nodes that lay
further from the root node will be  omitted. Note that  setting this option to 1
or 2 may greatly reduce  the computation time needed  for large code bases. Also
note that a graph  may be further truncated if  the graph's image dimensions are
not  sufficient    to     fit  the    graph    (see    MAX_DOT_GRAPH_WIDTH   and
MAX_DOT_GRAPH_HEIGHT).   If  0 is  used for the   depth value (the default), the
graph is not depth-constrained.""",

    GENERATE_LEGEND = \
"""If the GENERATE_LEGEND tag is set to YES (the default) Doxygen will generate
a legend page explaining the meaning of the various boxes and arrows in the dot
generated graphs.""",

    DOT_CLEANUP = \
"""If the DOT_CLEANUP tag is set to YES (the default) Doxygen will remove the
intermediate dot files that are used to generate the various graphs.""",


    # Configuration options related to the search engine

    SEARCHENGINE = \
"""The SEARCHENGINE tag specifies whether or not a search engine should be used.
If set to NO the values of all tags below this one will be ignored.""",

)


header = """Doxyfile

This file describes the settings to be used by the documentation system
doxygen (www.doxygen.org) for a project

All text after a hash (#) is considered a comment and will be ignored
The format is:
      TAG = value [value, ...]
For lists items can also be appended using:
      TAG += value [value, ...]
Values that contain spaces should be placed between quotes (\" \")
"""
