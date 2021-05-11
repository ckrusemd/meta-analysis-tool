Python Data Science Kickstarter 
==============================

<p><small>Based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

Builds on top of the excellent Cookiecutter DS template with services for FastAPI, Dash, MongoDB and pytest-html (creating a nice HTML of collected tests and displaying them through nginx).

I base other projects on this template and invite you to do the same - please raise issues and suggestions.

Docker services:
* FastAPI (:80)
* Plotly Dash Development, enabling hot-reload, (:8051)
* Plotly Dash Production, Gunicorn (:8050)
* MongoDB  (:27017)
* Pytest-html w/ nginx (:8080)

How to run:

```
docker-compose up -d
```
