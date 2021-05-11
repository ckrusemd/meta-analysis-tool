#!bin/bash

cd /app/
/usr/local/bin/python -m pip install --upgrade pip
pip install -r requirements.txt
uvicorn src.api.main:app --host 0.0.0.0 --port 80