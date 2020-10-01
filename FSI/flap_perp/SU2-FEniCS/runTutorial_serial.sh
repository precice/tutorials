#! /bin/bash 


tmux has-session -t PreciceSU2-FEniCS-newtonParallel 2>/dev/null
if [ "$?" -eq 1 ] ; then 
    tmux new-session -d -s 'PreciceSU2-FEniCS-newtonParallel'
    tmux split-window -v
fi


session=PreciceSU2-FEniCS-newtonParallel
window=${session}:0
pane_su2=${window}.0
pane_FEniCS=${window}.1
tmux send-keys -t "$pane_su2" C-z '( time $SU2_HOME/SU2_CFD/bin/SU2_CFD Fluid/euler_config_coupled.cfg ) 2>&1  | tee Su2.log' Enter
tmux select-pane -t "$pane_su2"
tmux select-window -t "$window"
tmux send-keys -t "$pane_FEniCS" C-z '( python3 Solid/perp-flap.py ) 2>&1  | tee FEniCS.log' Enter 
tmux select-pane -t "$pane_FEniCS"
tmux select-window -t "$window"
tmux attach-session -t "$session"

tmux detach -s "$session" 

rm -rf Output
mkdir Output
cp flow*.vtk Output
cp euler_config_coupled.cfg Output
cp flap.inp Output
cp FEniCS.log Output
cp Su2.log Output
cp precice-config.xml Output
cp point1.watchpoint.txt Output
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
