#! /bin/bash 


echo "Select folder name"
read fname
mkdir Output/$fname

session=autoSession
window=${session}:0
pane_su2=${window}.0
pane_Calculix=${window}.1
tmux send-keys -t "$pane_su2" C-z '( time ./../Solvers/SU2_fin/bin/SU2_CFD euler_config_coupled.cfg 2>&1 ) | tee su2.log' Enter
tmux select-pane -t "$pane_su2"
tmux select-window -t "$window"
tmux send-keys -t "$pane_Calculix" C-z '( time ccx_preCICE -i flap -precice-participant Calculix 2>&1 ) | tee calculix.log' Enter 
tmux select-pane -t "$pane_Calculix"
tmux select-window -t "$window"
tmux attach-session -t "$session"

tmux detach -s "$session" 

cp flow*.vtk Output/$fname 
cp euler_config_coupled.cfg Output/$fname 
cp flap.inp Output/$fname 
cp calculix.log Output/$fname 
cp su2.log Output/$fname  
cp precice-config.xml Output/$fname 
cp point1.watchpoint.txt Output/$fname 
cp plotDisplacement.sh Output/$fname 
# clean everything 
 ./Allclean

 echo " Copying of results was successfull! Let's hope that simulation went well as well" 
