#!/bin/bash
# Split-screen dashboard for the LOCAL CNSGA optimizer:
#   top pane  -> local_monitor.sh (eval progress / ETA / best-so-far / stage spread), refreshed every 5s
#   bottom pane -> htop filtered to this user's processes (per-core load, the 14 eval subprocesses, RAM)
#
#   bash dashboard.sh
#
# Detach with Ctrl-b then d (the optimizer keeps running). Re-attach: tmux attach -t linacmon.
# Kill the view (not the run): tmux kill-session -t linacmon.

BASE=/Users/rylandgoldman/Documents/Coding/Cornell2
SESSION=linacmon

tmux has-session -t "$SESSION" 2>/dev/null && exec tmux attach -t "$SESSION"

tmux new-session -d -s "$SESSION"
# top pane: the optimizer dashboard
tmux send-keys -t "$SESSION" "watch -n 5 -t bash $BASE/local_monitor.sh" C-m
# bottom pane: htop, this user's processes, sorted by CPU
tmux split-window -v -t "$SESSION"
tmux send-keys -t "$SESSION" "htop -u $USER -s PERCENT_CPU" C-m
# fix the monitor (top) at 20 rows so its full output always shows; htop takes the rest.
# resize interactively any time: Ctrl-b then hold Ctrl-Up/Ctrl-Down.
tmux resize-pane -t "$SESSION".0 -y 20
tmux select-pane -t "$SESSION".0
exec tmux attach -t "$SESSION"
