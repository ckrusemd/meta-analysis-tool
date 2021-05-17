from fastapi import Body, APIRouter
from pymongo import MongoClient
from pydantic import BaseModel
from loguru import logger
from ..dependencies import Helper

router = APIRouter()

Helper_ = Helper()

################################################
### MODELS
################################################

class EntrezQuery(BaseModel):
    query: str
    email: str

################################################
### QUERY
################################################/

@router.post("/mongodb_review/review/{project_id}/{uid}/{user_id}/{include}", tags=["MongoDB Review"], summary="Review a specific UID")
async def mongodb_review_uid(project_id: int, user_id: int, uid: int, include: bool):
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    uid_reviews = db['uid_reviews']
    uid_reviews_project = uid_reviews[project_id]

    logger.info( list( uid_reviews_project.find() ) )

    ## Check if exists
    if uid_reviews_project.find().count() == 0:
        logger.info( "insert new" )
        ## Insert new
        myquery = { 'uid' : uid, 'user_id': user_id, 'include' : include }
        uid_reviews_project.insert( myquery )
    else:
        logger.info( "update" )
        myquery = { 'uid' : uid, 'user_id': user_id }
        newvalues = { "$set": { 'include' : include  } }
        uid_reviews_project.update(myquery,newvalues)

    logger.info(" Når vi herned? ")
    result = list( uid_reviews_project.find({},{"_id": 0 })  )
    
    return result

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
    return list( db['projects'].find({}) )

@router.get("/mongodb_projects/users", tags=["MongoDB Projects"], summary="Get all projects in MongoDB")
async def get_all_mongodb_projects():
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    return list( db['users'].find({}) )

@router.get("/mongodb_projects/uid_reviews/{project_id}", tags=["MongoDB Projects"], summary="Get all UID reviews for a project in MongoDB")
async def get_all_mongodb_projects(project_id: int):
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    result = db['projects'].find_one({ 'uid' : project_id })
    return result
