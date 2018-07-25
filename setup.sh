#!/bin/bash
conda env create -f environment.yaml
python download_databases.py --file databases.yaml
