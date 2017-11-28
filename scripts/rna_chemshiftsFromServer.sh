#!/bin/bash

#for f in output/*; do
#	dest=${f/output/chemshifts}
#	dest=${dest/%pdb/ccs}
#	echo "Getting chemshifts for $f"
#	curl -X POST -F pdb=@$f http://50.63.157.7/RAMSEYWebService/upload/ > $dest
#
#done

if [ $# -eq 0 ]; then
	echo "Usage: rna_chemshiftsFromServer.sh <list of pdb-files>"
	exit 0
fi


for f in $@; do
	dest=${f/%pdb/ccs}
	echo -n "Polling RAMSEYWebService for chemshifts: $f -> $dest .. "
	curl --fail --silent -X POST -F pdb=@$f http://50.63.157.7/RAMSEYWebService/upload/ | sed 's/<.*>//g' > $dest
	echo "done"
done
