from fastapi import Depends, FastAPI

#from .dependencies import get_query_token, get_token_header
from .internal import admin
from .routers import test, mongodb, entrez, mongodb_entrez, mongodb_projects, mongodb_review

app = FastAPI(title="Meta-Analysis Tool API",
    description="API for the Meta-Analysis Tool Application",
    version="0.1",
    #dependencies=[Depends(get_query_token)]
    )


app.include_router(mongodb_review.router)
app.include_router(mongodb_projects.router)
app.include_router(mongodb_entrez.router)
app.include_router(entrez.router)
app.include_router(mongodb.router)
app.include_router(test.router)