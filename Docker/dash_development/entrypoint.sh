#!bin/bash

cd /app
pip install -r requirements.txt
cd /app/src/dash/
#uvicorn src.dash.app:server --reload  --host 0.0.0.0 --port 8050
python app.py