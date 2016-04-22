bash test_clean.sh
mkdir public
cp samples/base-2-length-2-compressed-ore.json public/template.json
./keygen --secparam 20
record00=`./encrypt '["0","00","0"]' | grep -v 'Starting\|Finished\|Generated\|Progress'`
record11=`./encrypt '["1","11","1"]' | grep -v 'Starting\|Finished\|Generated\|Progress'`
./eval '{"L":"'$record00'","R":"'$record11'"}'
