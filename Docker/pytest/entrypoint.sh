#!bin/bash

cd /app/
pip install -r requirements.txt

rm /var/www/html/index.nginx-debian.html
pytest -q tests --html=/var/www/html/index.nginx-debian.html --self-contained-html
nginx -g "daemon off;"
