Guide To Relax Scripts

When the user supplies a tag (such as "default"),
FastRelax/FastDesign will look for a file with the name:

main/database/sampling/relax_scripts/[tag].(sfxn).[options in alphabetical order].txt

Where elements in () are optional.
For example, <FastRelax relaxscript="legacy"/> will first look for "legacy.ref2015.txt" assuming ref2015 is the scorefunction.
If that does not exist, it will drop the score function and look for a general case: "legacy.txt"


At the time that this was written, the only option available is dualspace.
<FastRelax relaxscript="foo" dualspace="true"/> will first look for "foo.ref2015.dualspace.txt"
and fall back to "foo.dualspace.txt" if need be.


Just one more example:
<Scorefunciton name="sfxn" weights="beta_nov16"/>
<FastRelax relaxscript="foo" scorefxn="sfxn"/>
will look for "foo.beta_nov16.txt" and fall back to "foo.txt"


default.* should always be symlinks that point to the current default script
