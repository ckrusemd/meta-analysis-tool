# Meta-Analysis Tool
==============================

Why are systematic reviews still tediously done by hand from A-Z, saved in spreadsheets and intransparently across several software suites?

This project aims to accelerate the process using readily available tools:
* Biopython: We can query the Entrez API to collect abstract information and text from Pubmed.
* MongoDB: We can store the results as JSON to include, label and reproduce abstracts.
* Dash: We can make it all interactive with the Python Dash framework. And follow PRISMA guidelines.
* PythonMeta: We can arrive at the beautiful forest and funnel plots and customize the look.
* Docker: Keep it reproducible, save the images.

Please chip in with issues, suggestions and pull requests if you want to contribute.

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
