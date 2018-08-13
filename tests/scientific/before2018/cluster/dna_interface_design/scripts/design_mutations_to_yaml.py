#!/usr/bin/python

# ashworth
import sys,string,re,os

inputfiles = []
if len(sys.argv) > 1:
	inputfiles = sys.argv[1:]
else:
	fail = { '_isTestPassed' : False }
	print str(fail)
	sys.exit()

# example of file to be parsed
'''
Statistics for residues that were allowed to mutate:
RMS cutoff for conserved mutable residues: 1.00
-------------------------------------------------------------
AA  : --------Composition-------- ----------Recovery---------
AA  :    nat      f    des      f   cons      f cnsrot      f
-------------------------------------------------------------
ALA :    114  0.054     24  0.011     12  0.105     12  0.105
CYS :     10  0.005      0  0.000      0  0.000      0  0.000
ASP :     97  0.046    149  0.070     45  0.464     29  0.299
GLU :     64  0.030    202  0.095     15  0.234      9  0.141
PHE :     35  0.016     32  0.015     10  0.286      6  0.171
GLY :    161  0.076    130  0.061    107  0.665    107  0.665
HIS :     82  0.039    143  0.067     26  0.317     23  0.280
ILE :     64  0.030     15  0.007      9  0.141      9  0.141
LYS :    231  0.109     78  0.037     26  0.113     14  0.061
LEU :     40  0.019     17  0.008      2  0.050      2  0.050
MET :     25  0.012      3  0.001      0  0.000      0  0.000
ASN :    176  0.083    151  0.071     66  0.375     56  0.318
PRO :     55  0.026     25  0.012     17  0.309     17  0.309
GLN :    105  0.049    146  0.069     31  0.295     21  0.200
ARG :    286  0.135    543  0.256    163  0.570     60  0.210
SER :    255  0.120    208  0.098     82  0.322     70  0.275
THR :    213  0.100    114  0.054     50  0.235     48  0.225
VAL :     57  0.027     12  0.006      4  0.070      4  0.070
TRP :     12  0.006     45  0.021      0  0.000      0  0.000
TYR :     40  0.019     85  0.040     12  0.300      6  0.150
-------------------------------------------------------------
PRT :   2122  1.000   2122  1.000    677  0.319    493  0.232

(also...)

-------------------------------------------------------------
AA  : --------Composition-------- ----------Recovery---------
AA  :    nat      f    des      f   cons      f cnsrot      f
-------------------------------------------------------------
  A :    575  0.251    595  0.260    291  0.506    291  0.506
  C :    571  0.250    556  0.243    274  0.480    274  0.480
  G :    555  0.243    542  0.237    269  0.485    268  0.483
  T :    587  0.257    595  0.260    294  0.501    294  0.501
-------------------------------------------------------------
DNA :   2288  1.000   2288  1.000   1128  0.493   1127  0.493

'''
#                                (cat) :     nat        f     des        (f)    cons        (f)  cnsrot        (f)
re_stats = re.compile('([A-Z0-9-_. ]+) : +[0-9]+ +[0-9.]+ +[0-9]+ +([0-9.]+) +[0-9]+ +([0-9.]+) +[0-9]+ +([0-9.]+)')

stats = {} # for easy lookup for TestPassed criteria

for file in inputfiles:
	if not os.path.exists(file): continue
	for match in re_stats.finditer( open(file).read() ):
		(cat,desfrac,consfrac,cnsrotfrac) = match.groups()
		#print cat,desfrac,consfrac,cnsrotfrac
		for key,val in [ ('%s_DesignedComposition' %cat, desfrac),
		                 ('%s_DesignedConservation' %cat, consfrac),
		                 ('%s_DesignedRotamerConservation' %cat, cnsrotfrac) ]:
			if stats.has_key(key): print 'Warning!: overwriting %s' %key
			stats[key] = val

print '{',
# use old-style python sorting in case of old version
keys = stats.keys()
keys.sort()
for key in keys:
	print '\'%s\' : %s,' %(key,stats[key]),

passed = True
conditions = {
	'PRT_DesignedConservation' : lambda x: x >= 0.3,
	'PRT_DesignedRotamerConservation' : lambda x: x >= 0.2,
	'DNA_DesignedConservation' : lambda x: x >= 0.4,
	'DNA_DesignedRotamerConservation' : lambda x: x >= 0.3,
}

for key,func in conditions.items():
	if not stats.has_key(key):
		passed = False
		break
	if not func(stats[key]):
		passed = False
		break

print '\'_isTestPassed\': %s }' %passed
