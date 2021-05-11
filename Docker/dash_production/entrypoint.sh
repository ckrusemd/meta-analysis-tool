#!bin/bash

cd /app/
pip install -r requirements.txt
cd /app/src/dash/
gunicorn -w 1 -b :8050 src.dash.app:server