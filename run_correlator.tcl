puts "@@@ Starting"

# open the project, don't forget to reset
puts "@@@ Opening project"
#### NOTE the top from set_top must match the name of a function in the code included
#### in this way, Vivado will execute it as the entry point (imagine it as the 'main' of the code)
open_project -reset projCorrelator
# set_top is_in_boundaries_th
# set_top correlator_one
# set_top top_arr_correlator_one
# set_top correlator_stream
# set_top BRAM_to_corr
set_top BRAM_to_corr_nostream
# set_top top_arr_correlator_mult
add_files src/dataformats.cpp
add_files src/matching_LUTs.cpp
add_files src/matcher.cpp
add_files src/algo_parts.cpp
# add_files src/playground.cpp
add_files src/correlator.cpp
# add_files -tb random_tests.cpp
# add_files -tb simple_algo_ref.cpp
add_files -tb correlator_test_stream.cpp

puts "@@@ Opening solution"
# reset the solution
open_solution -reset "solution1"
# part options:
#	xcku9p-ffve900-2-i-EVAL
#	xc7vx690tffg1927-2
#	xcku5p-sfvb784-3-e
#	xcku115-flvf1924-2-i
#	xcvu9p-flga2104-2l-e
set_part {xc7vx690tffg1927-2}
create_clock -period 5 -name default

# do stuff
puts "@@@ C SIM"
csim_design

puts "@@@ C SYNTH"
csynth_design

puts "@@@ C/RTL COSYM"
cosim_design -trace_level all
# export_design -format ip_catalog  -vendor "cern-cms"

## cannot make the part below work, although trying with two different commands
## use the GUI (view waveform after cosim_design)
# puts "@@@ EXPORTING WAVEFORM"
# write_hw_ila_data my_hw_ila_data_file.zip [upload_hw_ila_data hw_ila_1]
# log_wave -r /

# exit Vivado HLS
exit
