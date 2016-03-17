#!/bin/bash
# Creates file with mixture of N sound files where N is first argument to
# script.
# Writes which files mixed and estimated pitches to outfile, as well as the
# parameters used.
# Results written in directory specified as second argument to script.
# Script must be called from the directory it resides in.
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
outfile=${OUTDIR}/`gdate -Ins`.mf0s
files=$(find ${SOUNDSPATH} -name "*${SNDTYPE}" \
    | gsort -R \
    | head -n${NSNDS} )
echo "Files:" >> "${outfile}"
echo "$files" >> "${outfile}"
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
N=4096
H=256
W=hanning
L=1024
f0min=36
f0max=72
f0d=1
K=40
sgm=1
thresh=-40
Fs=44100
Fmax=20000
P=${NSNDS}
echo "N=$N" >> "${outfile}"
echo "H=$H" >> "${outfile}"
echo "W=$W" >> "${outfile}"
echo "L=$L" >> "${outfile}"
echo "f0min=$f0min" >> "${outfile}"
echo "f0max=$f0max" >> "${outfile}"
echo "f0d=$f0d" >> "${outfile}"
echo "K=$K" >> "${outfile}"
echo "sgm=$sgm" >> "${outfile}"
echo "thresh=$thresh" >> "${outfile}"
echo "Fs=$Fs" >> "${outfile}"
echo "Fmax=$Fmax" >> "${outfile}"
echo "P=$P" >> "${outfile}"
echo "Analysis:" >> "${outfile}"
ffmpeg $in_s -filter_complex "${fgrph}"\
    -map "[aout]" -ac 1 -f f64le pipe:1\
    | octave -q ../f0sp_strm.m \
        $N \
        $H \
        $W \
        $L \
        $f0min \
        $f0max \
        $f0d \
        $K \
        $sgm \
        $thresh \
        $Fs \
        $Fmax \
        $P \
        text >> "${outfile}"
