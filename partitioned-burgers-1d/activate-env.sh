#!/bin/bash
# source this script to activate the virtual environment
# source activate-env.sh

VENV_DIR=".venv"

if [ -d "$VENV_DIR" ]; then
	: # env already exists
else
	echo "Creating virtual environment..."
	python3 -m venv $VENV_DIR
	# Activate the virtual environment and install packages
	source $VENV_DIR/bin/activate
	pip install -r requirements.txt
fi

source $VENV_DIR/bin/activate