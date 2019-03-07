#!/bin/bash
declare -i nevents=3000000
ptj="300.0"
particle_dir="C8S_UFO" # "Xquark_UFO", "C8S_UFO", "Coloron_UFO", "ZPB_UFO"
particle_tag="s8" # "qx", "s8", "cp", "zp"
process_type="qq" # "qg" or "gg" or "qq" or "qq" (for qx -> q g)

mg5_path="/home/jfilipek/Programs/MG5_aMC_v2_5_5/bin/mg5_aMC" # path to mg5 executable

# !!! Both have to end with '/' !!!
indir="/home/jfilipek/SampleGeneration/UFO_gen/Model_files/"  # Where UFO model files are
outdir="/home/jfilipek/SampleGeneration/UFO_gen/mg5_out/" # Where output will be generated

tmp_script="tmp_mg5_"$particle_tag"_"$process_type"_script.txt"

# END OF DEFINITIONS. START OF PROGRAM
declare -i runs=$nevents/10000
log_file="tmp_mg5_"$particle_tag"_"$process_type".log"
log_file=$(pwd)"/"$log_file
echo "Number of runs with nevents 10000: $runs"
echo "" > $tmp_script  # Clear file if exists.

# Process tags to actual processes in the madgraph
process="p p > "$particle_tag" j, "$particle_tag" >"
if [[ $process_type == "gg" ]]; then
	process=$process" g g"
elif [[ $process_type == "qq" ]]; then
	process=$process" q q~"
else
	process=$process" q g"
fi

# Start building script.
echo "import model $indir$particle_dir" >> $tmp_script
echo "define q = u c d s t b" >> $tmp_script # Has no effect if "gg"
echo "define q~ = u~ c~ d~ s~ t~ b~" >> $tmp_script

if [[ $particle_tag == "qx" ]]; then
	echo "define qx = ux cx dx sx tx bx" >> $tmp_script	
fi


echo "generate $process" >> $tmp_script
echo $process >> generate_history.log
echo "output "$outdir$particle_tag"_"$process_type >> $tmp_script
echo "launch -i" >> $tmp_script
echo "multi_run $runs" >> $tmp_script
echo "0" >> $tmp_script
echo "set ptj $ptj" >> $tmp_script
echo "0" >> $tmp_script
echo "Starting MG5 for $process"
echo "To follow the process open: tail -f $log_file"
$($mg5_path $tmp_script > $log_file)
# rm tmp_mg5_script.txt
