#!/bin/sh
# This script builds and installs your Maven project to the given Imagej2/Fiji directory
# It will delete all duplicate .jars which may break your ImageJ installation!
# It's recommended to use this script only on a clean copy of ImageJ which you don't
# use in your daily work
#
# author: Richard Domander (Royal Veterinary College)

if [ "$1" = "" ]
  then
    echo "No ImageJ directory specified"
else
    mvn -Dimagej.app.directory="$1" -Ddelete.other.versions=true clean install
fi
