# integration.addsan and integration.ubsan tests

-----
### integration.addsan
Run GCC's AddressSanitizer on the integration test runs. Failure indicates there are memory issues with the failing tests.


-----
### integration.ubsan
Run GCC's UndefinedBehaviorSanitizer on the integration test runs. Failure indicates the failing tests invoke C++ undefined behavior and hence may have unstable/uncertain results.
