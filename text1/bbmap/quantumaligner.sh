#!/bin/bash

usage(){
echo "
Written by Brian Bushnell
Last modified May 4, 2025

Description:  Aligns a query sequence to a reference using QuantumAligner.
The sequences can be any characters, but N is a special case.
Outputs the identity, rstart, and rstop positions.
Optionally prints a state space exploration map.
This map can be fed to visualizealignment.sh to make an image.

Usage:
quantumaligner.sh <query> <ref>
quantumaligner.sh <query> <ref> <map>
quantumaligner.sh <query> <ref> <map> <iterations>

Parameters:
query           A literal nucleotide sequence or fasta file.
ref             A literal nucleotide sequence or fasta file.
map             Optional output text file for matrix score space.
                Set to null for benchmarking with no visualization.
iterations      Optional integer for benchmarking multiple iterations.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

#This block allows symlinked shellscripts to correctly set classpath.
pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

calcXmx () {
    # Source the new scripts
    source "$DIR""/memdetect.sh"
    source "$DIR""/javasetup.sh"
    
    parseJavaArgs "--mem=200m" "--mode=fixed" "$@"
    
    # Set environment paths
    setEnvironment
}
calcXmx "$@"

align() {
	local CMD="java $EA $EOOM $SIMD $XMX $XMS -cp $CP aligner.QuantumAligner $@"
	#echo $CMD >&2
	eval $CMD
}

align "$@"
