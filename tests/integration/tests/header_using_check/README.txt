This integration test checks that you don't have any 'using' statements in headers.

'using' statements in headers are bad because they will apply not just to the header,
but also to any other file that includes the header - either via direct include or 
from an indirect include (a cc file which includes a header which includes another 
header, which includes a header using 'using') This can lead to issues where the 
compiler has difficulty resolving variable names due to conflicts from different 
definitions. (It defeats the purpose of having namespaces in the first place.)

By Rosetta coding convention, 'using namespace' is banned in headers. It's a nuclear
swiss-army chainsaw. Don't do it. This test will always error out as long as 
there is a 'using namespace' in a header.

Specific 'using' statements should only be used if they're a *deliberate* 
attempt to transclude an object from one namespace to another. Do NOT use '
using' statements to avoid typing out the full namespace in headers. 
This test will error out if there is a (non-namespace) using statement in headers, 
but that can be suppressed by placing the offending line (with the grep-style 
filename prefix) in known_usings.txt 
