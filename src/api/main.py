from fastapi import Depends, FastAPI

#from .dependencies import get_query_token, get_token_header
from .internal import admin
from .routers import test, mongodb

app = FastAPI(title="Python Data Science Kickstarter",
    description="FastAPI Template",
    version="0.1",
    #dependencies=[Depends(get_query_token)]
    )


app.include_router(test.router)
app.include_router(mongodb.router)