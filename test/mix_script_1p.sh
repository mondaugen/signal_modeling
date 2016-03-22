#!/bin/bash
if hash gdate 2>/dev/null; then
    DATE=gdate
else
    DATE=date
fi
if hash gawk 2>/dev/null; then
    AWK=gawk
else
    AWK=awk
fi
if hash gsort 2>/dev/null; then
    SORT=gsort
else
    SORT=sort
fi
SOUNDSPATH=${HOME}/Music/mums
if [[ -z "$1" ]]
then
        NSNDS=3
    else
        NSNDS=$1
fi
if [[ -z "$2" ]]
then
        OUTDIR='/tmp'
    else
        OUTDIR=$2
fi
if [[ -z "$3" ]]
then
    SNDTYPE="mp4"
else
    SNDTYPE=$3
fi
outfile=${OUTDIR}/`${DATE} -Ins`_mix
files=$(cat ${SOUNDSPATH}/single_pitch.INDEX \
    | ${SORT} -R \
    | head -n${NSNDS} | ${AWK} '{print "'${SOUNDSPATH}'/"$0}' )
echo "Files:" >> "${outfile}.INFO"
echo "$files" >> "${outfile}.INFO"
echo "Pitches:" >> "${outfile}.INFO"
echo "${files}" | grep -o -e '[a-gA-G]\#\?[0-9][^/]' \
    | grep -o -e '[a-gA-G]\#\?[0-9]' | python ./ptable.py \
    >> "${outfile}.INFO"
ldur=$(sort -rn <( echo "$files" \
    | xargs -I^ sh -c \
    'ffprobe -v error -show_entries stream=duration_ts -of \
     default=noprint_wrappers=1:nokey=1 -i "^";')\
    | head -n1)
echo "Max duration:"
echo "$ldur"
in_s=$(echo "$files" | gawk 'BEGIN {s=""} {s=s"-i "$0" "} END { print s}')
#echo "$in_s"
fgrph=$(echo "$files" | gawk 'BEGIN {s="";n=0}\
    {s=s"["n":a]apad=whole_len='${ldur}'[pad"n"],";\
        n=n+1; }\
        END { for (m=0;m<n;m++) { s=s"[pad"m"]" }\
        s=s"amix=inputs="n"[aout]";\
        print s }')
echo "$fgrph"
ffmpeg $in_s -filter_complex "${fgrph}"\
    -map "[aout]" -ac 1 -c:a pcm_s16le -f wav ${outfile}.wav
echo "Wrote to: ${outfile}.wav"
