#!/bin/bash
# Stream only as much as we need: head closes the pipe, curl stops transferring.
url=$1; out=$2; n=$3
curl -s --http1.1 "$url" 2>/dev/null | gzip -dc 2>/dev/null | head -n $((n*4)) | gzip -1 > "$out"
echo "$out $(( $(gzip -dc "$out" | wc -l) / 4 )) reads" >> grab.log
