#!/bin/bash
if [[ -z "$1" ]]
then
        NSNDS=3
    else
        NSNDS=$1
fi
files=$(find ${HOME}/Music/mums -name "*mp4" \
    | gsort -R \
    | head -n${NSNDS} )
echo "Files:"
echo "$files"
ldur=$(sort -rn <( echo "$files" \
    | xargs -I^ sh -c \
    'ffprobe -v error -show_entries stream=duration_ts -of \
     default=noprint_wrappers=1:nokey=1 -i "^";')\
    | head -n1)
echo "Max duration:"
echo "$ldur"
in_s=$(echo "$files" | gawk 'BEGIN {s=""} {s=s"-i "$0" "} END { print s}')
#echo "$in_s"
outfile=/tmp/$RANDOM.wav
fgrph=$(echo "$files" | gawk 'BEGIN {s="";n=0}\
    {s=s"["n":a]apad=whole_len='${ldur}'[pad"n"],";\
        n=n+1; }\
        END { for (m=0;m<n;m++) { s=s"[pad"m"]" }\
        s=s"amix=inputs="n"[aout]";\
        print s }')
echo "$fgrph"
ffmpeg $in_s -filter_complex "${fgrph}"\
    -map "[aout]" -ac 1 -c:a pcm_s16le -f wav $outfile
echo "Wrote to: $outfile"
