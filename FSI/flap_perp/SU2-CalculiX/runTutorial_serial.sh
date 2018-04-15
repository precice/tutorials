#! /bin/bash 


ln -sf precice-config-serial.xml precice-config.xml

tmux has-session -t PreciceSU2-Calculix-newtonParallel 2>/dev/null
if [ "$?" -eq 1 ] ; then 
    tmux new-session -d -s 'PreciceSU2-Calculix-newtonParallel'
    tmux split-window -v
fi


session=PreciceSU2-Calculix-newtonParallel
window=${session}:0
pane_su2=${window}.0
pane_Calculix=${window}.1
tmux send-keys -t "$pane_su2" C-z '( time SU2_CFD euler_config_coupled.cfg 2>&1 ) | tee Su2.log' Enter
tmux select-pane -t "$pane_su2"
tmux select-window -t "$window"
tmux send-keys -t "$pane_Calculix" C-z '( time ccx_preCICE -i flap -precice-participant Calculix 2>&1 ) | tee Calculix.log' Enter 
tmux select-pane -t "$pane_Calculix"
tmux select-window -t "$window"
tmux attach-session -t "$session"

tmux detach -s "$session" 

rm -f Output/*
cp flow*.vtk Output
cp euler_config_coupled.cfg Output
cp flap.inp Output
cp Calculix.log Output
cp Su2.log Output
cp precice-config.xml Output
cp point1.watchpoint.txt Output
cp plotDisplacement.sh Output
# clean everything 
rm -f *log
rm -f *vtk
rm -f *txt
rm -f *csv
rm -f *out
rm -f restart_flow*
rm -f forces*
rm -f flap.[^i]*

 echo " Copying of results was successfull! Let's hope that simulation went well as well. Results are in the output
 folder" 
