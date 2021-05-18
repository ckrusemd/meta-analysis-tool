from fastapi import Body, APIRouter
from pymongo import MongoClient
from pydantic import BaseModel
from loguru import logger
from ..dependencies import Helper
import json

router = APIRouter()

Helper_ = Helper()

################################################
### MODELS
################################################

class EntrezQuery(BaseModel):
    query: str
    email: str

class ProjectModel(BaseModel):
    name: str
    description: str
    query: str
    user_id: int

################################################
### QUERY
################################################

@router.post("/mongodb_projects/clean_collections/", tags=["MongoDB Projects"], summary="Erase all content of Projects-MongoDB collections")
async def clean_mongodb_projects_collections():
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    user_list = db['users'] 
    user_list.delete_many({})
    projects = db['projects']
    projects.delete_many({})
    uid_reviews = db['uid_reviews']
    uid_reviews.delete_many({})
    counts = {  'users' : db['users'].find().count(),
                'projects' : db['projects'].find().count(),
                'uid_reviews' : db['uid_reviews'].find().count()
    }
    return counts

@router.post("/mongodb_projects/generate_collections/", tags=["MongoDB Projects"], summary="Generate empty Projects dictionaries in MongoDB ")
async def generate_empty_mongodb_projects_collections():
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    user_list = db['users'] 
    user_list = []
    projects = db['projects']
    projects = []
    uid_reviews = db['uid_reviews']
    uid_reviews = []
    counts = {  'users' : db['users'].find().count(),
                'projects' : db['projects'].find().count(),
                'uid_reviews' : db['uid_reviews'].find().count()
    }
    return counts

@router.get("/mongodb_projects/projects", tags=["MongoDB Projects"], summary="Get all projects in MongoDB")
async def get_all_mongodb_projects():
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    results = db['projects'].find({},{"_id":0})
    logger.info( results )
    return list( results )

@router.get("/mongodb_projects/users", tags=["MongoDB Projects"], summary="Get all users in MongoDB")
async def get_all_mongodb_users():
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    results = db['users'].find({},{"_id":0})
    logger.info( results )
    return list( results )

@router.get("/mongodb_projects/uid_reviews/{project_id}", tags=["MongoDB Projects"], summary="Get all UID reviews for a project in MongoDB")
async def get_all_mongodb_uid_reviews(project_id: int):
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    result = db['projects'].find_one({ 'uid' : project_id },{"_id":0})
    return result

@router.post("/mongodb_projects/create_project/", tags=["MongoDB Projects"], summary="Create new project")
async def create_project(query_body: ProjectModel = Body(...,example={ 
                                                                        'name': 'Kruse Eiken Vestergaard', 
                                                                        'description' : 'Articles Done By Me',
                                                                        'query' : 'kruse eiken vestergaard',
                                                                        'user_id' : 1  } ) ):

    client = Helper_.getMongoDbClient()
    db = client['Projects']
    query_json = dict( query_body )

    ## Check if exists:
    exists_already = db['projects'].find_one(query_json,{"_id":0})
    if exists_already == None:
        result = db['projects'].insert_one( query_json )
        result = list( db['projects'].find({},{"_id" : 0}) )
        return result
    else:
        return exists_already