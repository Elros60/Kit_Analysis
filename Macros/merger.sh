#!/bin/bash
# Load O2 env before running this script
#TARGET_DIR=$1
LOCAL_DIR=$1
#RUN=$3
#PASS=$4

# NEW_DIR=$PASS
# mkdir -p $LOCAL_DIR/$NEW_DIR
cd $LOCAL_DIR/
RUN_LIST=`ls -l $LOCAL_DIR | grep "^d" | awk -F" " '{print $9}'`
echo $RUN_LIST
for RUN in $RUN_LIST
do
	cd $LOCAL_DIR/$RUN
	echo "Start processing $LOCAL_DIR/$RUN"
	hadd -f AnalysisResults.root $LOCAL_DIR/$RUN/*/AnalysisResults.root
done
cd $LOCAL_DIR/