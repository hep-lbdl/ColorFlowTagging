#!/bin/bash
declare -i nevents=600000
ptj="300.0"
mass="125.0"
width="0.00407"
particle_dir="Xquark_UFO"
particle_tag="qx"
process_type="qg" # "gg" or "qq" or "qg" or "q~g"

# ENVIRONMENT VARIABLES
mg5loc="/phys/users/balbok/Programs/MG5_aMC_v2_5_5/bin/mg5_aMC"
ufoloc="/phys/users/balbok/Projects/UFO/Model_files/$particle_dir"

script_name="tmp_mg5_script_"$particle_tag"_"$process_type".txt"

log_file="tmp_mg5_"$particle_tag"_"$process_type".log"
log_file=$(pwd)"/"$log_file

# END OF DEFINITIONS. START OF PROGRAM
declare -i runs=$nevents/10000
echo "Number of runs with nevents 10000: $runs"
echo "" > $script_name  # Clear file if exists.

# Process tags to actual processes in the madgraph
process="p p > "$particle_tag" j, "$particle_tag" >"
if [[ $process_type = "gg" ]]; then
	process=$process" g g"
elif [[ $process_type = "qq" ]]; then
	process=$process" q q~"
elif [[ $process_type = "qg" ]]; then
	process=$process" q g"
else
	process=$process" q~ g"
fi

# Start building script.
if [[ $particle_tag = "h" ]]; then
	echo "import model heft" >> $script_name
else
	echo "import model $ufoloc" >> $script_name
fi

echo "define q = u c d s t b" >> $script_name # Has no effect if "gg"
echo "define q~ = u~ c~ d~ s~ t~ b~" >> $script_name

if [[ $particle_tag == "qx" ]]; then
	echo "define qx = ux cx dx sx tx bx" >> $script_name
elif [[ $particle_tag == "qx~" ]]; then
	echo "define qx~ = ux~ cx~ dx~ sx~ tx~ bx~" >> $script_name
fi

echo "generate $process" >> $script_name
echo $process >> generate_history.log
echo "output ~/Projects/UFO/mg5_out/"$particle_tag"_"$process_type >> $script_name
echo "launch -i" >> $script_name
echo "multi_run $runs" >> $script_name
echo "0" >> $script_name
echo "set lhe_version 2" >> $script_name
echo "set ptj $ptj" >> $script_name

if [[ $particle_tag == "qx" ]]; then
	echo "set mass mux $mass" >> $script_name
	echo "set mass mdx $mass" >> $script_name
	echo "set mass msx $mass" >> $script_name
	echo "set mass mcx $mass" >> $script_name
	echo "set mass mtx $mass" >> $script_name
	echo "set mass mbx $mass" >> $script_name

	echo "set decay wux $width" >> $script_name
	echo "set decay wdx $width" >> $script_name
	echo "set decay wsx $width" >> $script_name
	echo "set decay wcx $width" >> $script_name
	echo "set decay wtx $width" >> $script_name
	echo "set decay wbx $width" >> $script_name
else
	echo "set mass m"$particle_tag" $mass" >> $script_name
	echo "set decay w"$particle_tag" $width" >> $script_name
fi

echo "0" >> $script_name
echo "Starting MG5 for $process"
echo "To follow the process open: tail -f $log_file"
$($mg5loc $script_name > $log_file)
