
cd /Users/vmulligan/rosetta_devcopy/Rosetta/main/tests/integration/new/sasa_metric_options

# Do the tests actually exist?
[ -x /Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease ] || exit 1

/Users/vmulligan/rosetta_devcopy/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease  @flags -database /Users/vmulligan/rosetta_devcopy/Rosetta/main/database \
    -parser:protocol test_simple_metric_filter.xml -testing:INTEGRATION_TEST 2>&1 | egrep -vf ../../ignore_list > log
test "${PIPESTATUS[0]}" != '0' && exit 1 || true