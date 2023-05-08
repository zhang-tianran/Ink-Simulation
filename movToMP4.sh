#! /bin/bash

# Convert a given .mov file to an .mp4 file,
# optionally speeding it up or muting it

if [[ $# != 2  &&  $# != 3 ]]
then
  echo "Usage: ./movToMP4.sh <mov_filename> <speedup> <mute?>"
  exit 1
fi

FULL_FILEPATH=$1
FILEPATH="${FULL_FILEPATH%.*}"

# Verify that file exists
if [[ ! -e "${FULL_FILEPATH}" ]]
then
  echo "File does not exist"
  exit 1
fi

if [[ $# == 3 && $3 == "mute" ]]
then
  echo "Producing ${FILEPATH}.mp4, sped up by $2x and muted"
  ffmpeg -i "${FILEPATH}.mov" -vcodec h264 -filter:v "setpts=PTS/$2" -an "${FILEPATH}.mp4" < /dev/null > /dev/null 2>&1
else
  echo "Producing ${FILEPATH}.mp4, sped up by $2x"
  ffmpeg -i "${FILEPATH}.mov" -vcodec h264 -filter:v "setpts=PTS/$2" "${FILEPATH}.mp4" < /dev/null > /dev/null 2>&1
fi