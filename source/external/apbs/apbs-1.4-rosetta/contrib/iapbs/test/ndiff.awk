### -*-awk-*-
### ====================================================================
###  @Awk-file{
###     author          = "Nelson H. F. Beebe",
###     version         = "1.00",
###     date            = "28 January 2000",
###     time            = "08:38:46 MST",
###     filename        = "ndiff.awk",
###     copyright       = "Copyright (c) 2000 Nelson H. F. Beebe. This
###                        code is licensed under the GNU General Public
###                        License, version 2 or later.",
###     address         = "Center for Scientific Computing
###                        University of Utah
###                        Department of Mathematics, 322 INSCC
###                        155 S 1400 E RM 233
###                        Salt Lake City, UT 84112-0090
###                        USA",
###     telephone       = "+1 801 581 5254",
###     FAX             = "+1 801 585 1640, +1 801 581 4148",
###     URL             = "http://www.math.utah.edu/~beebe",
###     checksum        = "10489 476 1904 14057",
###     email           = "beebe@math.utah.edu, beebe@acm.org,
###                        beebe@ieee.org (Internet)",
###     codetable       = "ISO/ASCII",
###     keywords        = "numerical file differencing",
###     supported       = "yes",
###     docstring       = "This program compares two putatively similar
###                        files, ignoring small numeric differences.
###                        Complete documentation can be found in the
###                        accompanying UNIX manual page file,
###                        ndiff.man.
###
###                        Usage:
###                       	awk -f ndiff.awk \
###                       		[-v ABSERR=x] \
###                       		[-v FIELDS=n1a-n1b,n2,n3a-n3b,...] \
###                       		[-v FS=regexp] \
###                       		[-v MINWIDTH=n] \
###                       		[-v QUIET=n] \
###                       		[-v RELERR=x] \
###                       		[-v SILENT=n] \
###                       		infile1 infile2
###
###                        The checksum field above contains a CRC-16
###                        checksum as the first value, followed by the
###                        equivalent of the standard UNIX wc (word
###                        count) utility output of lines, words, and
###                        characters.  This is produced by Robert
###                        Solovay's checksum utility.",
###  }
### ====================================================================

BEGIN \
{
    initialize()

    compare_files(ARGV[1], ARGV[2])

    exit (Ndiff != 0)
}


function abs(a)
{
    ## Return the absolute value of the argument.

    return ((a < 0) ? -a : a)
}


function awkfloat(s)
{
    ## Convert a numeric string to an awk floating-point number, and
    ## return the result as a floating-point number.
    ##
    ## Fortran use has any of E, e, D, d, Q, or q, or even nothing at
    ## all, for the exponent letter, but awk and C only allow E and e.
    ##
    ## Ada usefully permits nonsignificant underscores for
    ## readability: 3.14159265358979323846 and
    ## 3.14159_26535_89793_23846 are equivalent.
    ##
    ## We can safely assume that there are no leading or trailing
    ## whitespace characters, because all strings passed to this
    ## function are the result of splitting lines into
    ## whitespace-delimited fields.

    gsub("_","",s)		# remove Ada-style separators
    gsub("[DdQq]","e",s)	# convert Fortran exponent letters to awk-style
    if (match(s,"[0-9.][-+][0-9]+$")) # then letter-less exponent
	s = substr(s,1,RSTART) "e" substr(s,RSTART+1) # insert exponent letter e
    return (0 + s)		# coerce to a number
}


function compare_all(f1line,f2line,f1parts,f2parts,n, k)
{
    ## Compare all fields in f1line and f2line, assuming that they have
    ## already been split into n parts in f1parts[] and f2parts[].
    ##
    ## If any fields differ, print a diff-style report, and increment
    ## global variable Ndiff,

    for (k = 1; k <= n; ++k)
    {
	if (diff_field(f1parts[k], f2parts[k], k) != 0)
	{
	    report_difference(f1line,f2line,k)
	    return
	}
    }
}


function compare_files(file1,file2, f1line,f2line,f1parts,f2parts,n1,n2)
{
    ## Compare all lines in two files, printing a diff-style report of
    ## differences.  If any numeric differences have been found, print a
    ## one-line report of which matching line had the largest numeric
    ## difference.  Finally, print a diagnostic if the files differ in
    ## length.

    NRLINE = 0
    while (((getline f1line < file1) > 0) && \
	   ((getline f2line < file2) > 0))
    {
	NRLINE++
	n1 = split(f1line,f1parts)
	n2 = split(f2line,f2parts)
	if (n1 == n2)
	{
	    if (N_Fields == 0)
		compare_all(f1line,f2line,f1parts,f2parts,n1)
	    else
		compare_some(f1line,f2line,f1parts,f2parts,n1)
	}
	else
	    report_difference(f1line,f2line,max(n1,n2))
    }
    if (QUIET == 0)
    {
	if (Max_Abserr > 0)
	    printf("### Maximum absolute error in matching lines = %.2e at line %d field %d\n", \
		   Max_Abserr, Max_Abserr_NR, Max_Abserr_NF)
	if (Max_Relerr > 0)
	    printf("### Maximum relative error in matching lines = %.2e at line %d field %d\n", \
		   Max_Relerr, Max_Relerr_NR, Max_Relerr_NF)
    }
    if ((getline f1line < file1) > 0)
	warning("file " file2 " is short")
    if ((getline f2line < file2) > 0)
	warning("file " file1 " is short")
}


function compare_some(f1line,f2line,f1parts,f2parts,n, k,m)
{
    ## Compare selected fields in f1line and f2line, assuming that they
    ## have already been split into n parts in f1parts[] and f2parts[].
    ## The globals (N_Fields, Fields[]) define which fields are to be
    ## compared.
    ##
    ## If any fields differ, print a diff-style report, and increment
    ## global variable Ndiff.

    for (k = 1; (k <= N_Fields) && (k <= n); ++k)
    {
	m = Fields[k]
	if ((m <= n) && (diff_field(f1parts[m], f2parts[m], m) != 0))
	{
	    report_difference(f1line,f2line,m)
	    return
	}
    }
}


function diff_field(field1,field2,nfield)
{
    ## If both fields are identical as strings, return 0.
    ##
    ## Otherwise, if both fields are numeric, return 0 if they are close
    ## enough (as determined by the globals ABSERR and RELERR), or are
    ## both ignorable (as determined by MINWIDTH), and otherwise return
    ## 1.
    ##
    ## Otherwise, return 1.
    ##
    ## The computed absolute and relative errors are saved in global
    ## variables (This_Abserr and This_Relerr) for later use in
    ## diagnostic reports.  These values are always zero for
    ## nonnumeric fields.

    This_Abserr = 0
    This_Relerr = 0

    if (field1 == field2) # handle the commonest, and easiest, case first
	return (0)
    else if ((field1 ~ NUMBER_PATTERN) && (field2 ~ NUMBER_PATTERN))
    {
	## Handle MINWIDTH test while the fields are still strings
	if (ignore(field1) && ignore(field2))
	    return (0)

	## Now coerce both fields to floating-point numbers,
	## converting Fortran-style exponents, if necessary.
	field1 = awkfloat(field1)
	field2 = awkfloat(field2)

	This_Abserr = abs(field1 - field2)
	This_Relerr = maxrelerr(field1,field2)
	if ( ((ABSERR != "") && (This_Abserr > ABSERR)) || \
	     ((RELERR != "") && (This_Relerr > RELERR)) )
	{
	    if (This_Abserr > Max_Abserr)
	    {
		Max_Abserr_NF = nfield
		Max_Abserr_NR = NRLINE
		Max_Abserr = This_Abserr
	    }
	    if (This_Relerr > Max_Relerr)
	    {
		Max_Relerr_NF = nfield
		Max_Relerr_NR = NRLINE
		Max_Relerr = This_Relerr
	    }
	    return (1)
	}
	else
	    return (0)
    }
    else
	return (1)
}


function error(message)
{
    ## Issue an error message and terminate with a failing status code.

    warning("ERROR: " message)
    exit(1)
}


function ignore(field)
{
    ## Return 1 if field is ignorable, because it is shorter than
    ## MINWIDTH and appears to be a real number.  Otherwise, return 0.

    return ((MINWIDTH > 0) && \
	    (length(field) < MINWIDTH) && \
	    (field ~ "[.DdEeQq]"))
}


function initialize( eps)
{
    ## Process command-line options, and initialize global variables.

    Stderr = "/dev/stderr"

    Macheps = machine_epsilon()

    if (ABSERR != "")
	ABSERR = abs(awkfloat(ABSERR)) # coerce to positive number

    if (RELERR != "")
    {
	RELERR = abs(awkfloat(RELERR)) # coerce to positive number
	if (RELERR < Macheps)
	    warning("RELERR = " RELERR " is below machine epsilon " Macheps)
	else if (RELERR >= 1)	# RELERR=nnn means nnn*(machine epsilon)
	    RELERR *= Macheps
    }

    if ((ABSERR == "") && (RELERR == "")) # supply default (see man pages)
	RELERR = max(1.0e-15, 8.0 * Macheps)

	## printf( "RELERR is %15.10f\n", RELERR )
    ## Coerce remaining options to numbers
    MINWIDTH += 0
    QUIET += 0
    SILENT += 0

    Max_Relerr = 0
    Max_Relerr_NR = 0
    Max_Relerr_NF = 0

    Max_Abserr = 0
    Max_Abserr_NR = 0
    Max_Abserr_NF = 0

    This_Abserr = 0
    This_Relerr = 0

    if (FIELDS != "")
	initialize_fields()
    else
	N_Fields = 0

    ## The precise value of this regular expression to match both an
    ## integer and a floating-point number is critical, and documented
    ## in the accompanying manual page: it must match not only the
    ## awk- and C-style -nnn, -n.nnn, and -n.nnne+nn, but also the
    ## Fortran styles -nnn, -n.nnn, -n.D+nn, -.nnD+nn, -nD+nn,
    ## -n.nnnQ+nn, -n.nnnd+nn, and -n.nnn+nnn.  The Fortran forms will
    ## be converted by awkfloat() to awk-form.  Ada permits an
    ## nonsignificant underscore between digits, so we support that
    ## too.

    NUMBER_PATTERN = "^[-+]?([0-9](_?[0-9])*([.]?([0-9](_?[0-9])*)*)?|[.][0-9](_?[0-9])*)([DdEeQq]?[-+]?[0-9](_?[0-9])*)?$"

    Ndiff = 0
    if (ARGC != 3)
	error("Incorrect argument count\n\tUsage: awk -f ndiff.awk [-v ABSERR=x] [-v FIELDS=n1a-n1b,n2,n3a-n3b,...] [-v FS='regexp'] [-v MINWIDTH=n] [-v RELERR=x] infile1 infile2")
}


function initialize_fields( j,k,m,n,numbers,parts)
{
    ## Convert a FIELDS=n1a-n1b,n2,n3a-n3b,... specification to a list
    ## of N_Fields numbers in Fields[].

    N_Fields = 0
    n = split(FIELDS,parts,",")
    for (k = 1; k <= n; ++k)
    {
	m = split(parts[k],numbers,"-+")
	if (m == 1)
	{
	    if (parts[k] !~ "^[0-9]+$")
		error("non-numeric FIELDS value [" parts[k] "]")
	    else if (parts[k] == 0)
		error("zero FIELDS value [" parts[k] "]: fields are numbered from 1")
	    else
		Fields[++N_Fields] = parts[k]
	}
	else if (m == 2)
	{
	    if ((numbers[1] !~ "^[0-9]+$") || \
		(numbers[2] !~ "^[0-9]+$"))
		error("non-numeric FIELDS range [" parts[k] "]")
	    else if ((numbers[1] == 0) || (numbers[2] == 0))
		error("zero value in FIELDS range [" parts[k] "]: fields are numbered from 1")
	    else if (numbers[1] > numbers[2])
		error("bad FIELDS range [" parts[k] "]")
	    else if ((numbers[2] - numbers[1] + 1) > 100)
		error("FIELDS range [" parts[k] "] exceeds 100")
	    else
	    {
		for (j = numbers[1]; j <= numbers[2]; ++j)
		    Fields[++N_Fields] = j
	    }
	}
	else
	    error("bad FIELDS range [" parts[k] "]")
    }
    ## printf("DEBUG: Fields = [")
    ## for (k = 1; k <= N_Fields; ++k)
    ##     printf("%d,", Fields[k])
    ## print "]"
    ## exit(0)
}


function machine_epsilon( x)
{
    ## Tests on these architectures with awk, gawk, mawk, and nawk all
    ## produced identical results:
    ##
    ##		Apple Macintosh PPC G3	Rhapsody 5.5
    ##		DEC Alpha		OSF/1 4.0F
    ##		HP 9000/735		HP-UX 10.01
    ##		IBM PowerPC		AIX 4.2
    ##		Intel Pentium III	GNU/Linux 2.2.12-20smp (Redhat 6.1)
    ##		NeXT Turbostation	Mach 3.3
    ##		SGI Indigo/2		IRIX 5.3
    ##		SGI Origin 200		IRIX 6.5
    ##		Sun SPARC		GNU/Linux 2.2.12-42smp (Redhat 6.1)
    ##		Sun SPARC		Solaris 2.6
    ##		Sun SPARC		Solaris 2.7
    ##
    ##		/usr/local/bin/awk:  2.22045e-16
    ##		/usr/local/bin/gawk: 2.22045e-16
    ##		/usr/local/bin/mawk: 2.22045e-16
    ##		/usr/local/bin/nawk: 2.22045e-16
    ##
    ## Thus, there does not appear to be concern for surprises from
    ## long registers, such as on the Intel x86 architecture.

    x = 1.0
    while ((1.0 + x/2.0) != 1.0)
	x /= 2.0
    return (x)
}


function max(a,b)
{
    ## Return the (numerically or lexicographically) larger of the two
    ## arguments.

    return ((a > b) ? a : b)
}


function maxrelerr(x,y)
{
    ## Return the maximum relative error of two values.

    #x = abs(x + 0)		# coerce to nonnegative numbers
    #y = abs(y + 0)    		# coerce to nonnegative numbers

    ## See the documentation of the -relerr option in ndiff.man for the
    ## explanation of this complex definition:

    if (x == y)
	return (0)
    else if ((x != 0) && (y != 0))
	return (abs(x-y)/min(abs(x),abs(y)))
    else if ((x == 0) && (y != 0))
	return (1)
    else if ((y == 0) && (x != 0))
	return (1)
    else
	return (0)
}


function min(a,b)
{
    ## Return the (numerically or lexicographically) smaller of the two
    ## arguments.

    return ((a < b) ? a : b)
}


function report_difference(f1line,f2line,nfield, emult)
{
    ## Print a diff-style difference of two lines, but also show in
    ## the separator line the field number at which they differ, and
    ## the global absolute and relative errors, if they are nonzero.

    if (SILENT == 0)
    {
	printf("%dc%d\n", NRLINE, NRLINE)
	printf("< %s\n", f1line)
	## if ((This_Abserr != 0) || (This_Relerr != 0))
	## {
	##     emult = This_Relerr / Macheps
	##     if (emult >= 10000)
	##	printf("--- field %d\tabsolute error %.2e\trelative error %.2e\n",
	##	       nfield, This_Abserr, This_Relerr)
	##    else
	##	printf("--- field %d\tabsolute error %.2e\trelative error %.2e [%d*(machine epsilon)]\n",
	##	       nfield, This_Abserr, This_Relerr, int(emult + 0.5))
	##}
	##else
	##    printf("--- field %d\n", nfield)
	printf("> %s\n", f2line)
    }
    Ndiff++
}


function warning(message)
{
    ## Print a warning message on stderr, using emacs
    ## compile-command-style message format.

    if (FNR > 0)
	print FILENAME ":" FNR ":%%" message >Stderr
    else	# special case for diagnostics during initialization
	print message >Stderr
}
