bash test_clean.sh
mkdir public
cp samples/base-2-length-2-compressed-ore.json public/template.json
./keygen -C --secparam 20
record00=`./encrypt -C '["0","00","0"]' | grep -v 'Starting\|Finished\|Generated\|Progress'`
record11=`./encrypt -C '["1","11","1"]' | grep -v 'Starting\|Finished\|Generated\|Progress'`
./eval -C '{"L":"'$record00'","R":"'$record11'"}'
