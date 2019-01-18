#! /bin/bash

tmux has-session -t PreciceCalculix 2>/dev/null
if [ "$?" -eq 1 ] ; then 
    tmux new-session -d -s 'PreciceCalculix'
    tmux split-window -v
fi

session=PreciceCalculix
window=${session}:0
pane_Calculix1=${window}.0
pane_Calculix2=${window}.1
tmux send-keys -t "$pane_Calculix1" C-z '( time ccx_preCICE -i beam1 -precice-participant Calculix1 ) 2>&1  | tee Calculix1.log' Enter
tmux select-pane -t "$pane_Calculix1"
tmux select-window -t "$window"
tmux send-keys -t "$pane_Calculix2" C-z '( time ccx_preCICE -i beam2 -precice-participant Calculix2 ) 2>&1  | tee Calculix2.log' Enter 
tmux select-pane -t "$pane_Calculix2"
tmux select-window -t "$window"
tmux attach-session -t "$session"

tmux detach -s "$session" 
