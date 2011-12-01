#################################
#
# Generate features database from a sample source
#
# See https://wiki.rosettacommons.org/index.php/FeaturesScientificBenchmark
#
# OUTPUT:
#
#   features_${sample_source_id}_${date_code}.db3
#
#   containing features extracted from the input structures
#
#
################################
################################
# You probably need to set these:

mini=../../../../../..
platform=linuxgccrelease
database=$ROSETTA3_DB

sample_source_id=local_example
# e.g. "top4400", "native_interfaces", "rosetta_r1234_correct_abinitio"

sample_source_description="Local Example"
# e.g. "Rosetta revision 1234 with correct flag an ab initio protocol"

# IMPORTANT REMEBER TO EDIT:
#   - features.xml.TEMPLATE  <- specify which features to extract
#   - flags                  <- specify which how to run Rosetta

################################
################################
# You probably don't need to edit this

date_code=`date '+%y%m%d'`
run_id=${sample_source_id}_${date_code}
database_fname=features_${run_id}.db3
log_fname=ReportToDB.${run_id}.log

sed -e "s@DATABASE_FNAME@${database_fname}@g" -e "s@SAMPLE_SOURCE_DESCRIPTION@${sample_source_description}@g" features.xml.TEMPLATE > features.xml

echo "Generating Featuers Database:"
echo ""
echo ""
echo "writing features.xml from features.xml.TEMPLATE:"
echo "    -RUN_ID ->" ${run_id}
echo "    -SAMPLE_SOURCE_DESCRIPTION ->" ${sample_source_description}
echo ""
echo "writing features databse:"
echo "    -${database_fname}"
echo ""
echo "writing log:"
echo "    -${log_fname}"
echo ""
echo ""

rm -f features_${run_id}.db3

cmd="time $mini/bin/rosetta_scripts.$platform -database ${database} @flags"
echo $cmd "&> ${log_fname}"
$cmd &> ${log_fname}

echo ""
if [ !$? ]; then
  echo "ERROR: Failed to create ${database_fname}. Check ${log_fname} for errors."
else
  echo "Created ${database_fname}."
  echo "Now either move or ${database_fname} somewhere within test/scientific/cluster/features/sample_source/"
  echo "and run analysis scripts."
  echo "Try running \'R CMD BATCH make_plots.R\' in test/scientific/cluster/features"
fi 
###############################