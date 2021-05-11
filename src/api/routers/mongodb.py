from fastapi import APIRouter
from pymongo import MongoClient
from loguru import logger

router = APIRouter()


@router.get("/mongodb/testConnection/", tags=["MongoDB"], summary="Test Connection to MongoDB Container (mongodb:27017)")
async def test_mongodb_connection():
    client = MongoClient('mongodb', 27017)
    logger.info(client)
    logger.info(client.list_database_names())
    didConnect = client.admin.command('ismaster')
    logger.info( didConnect )
    return [ { "connection": client.list_database_names() } ]


@router.get("/mongodb/showAllDatabasesAndCollections/", tags=["MongoDB"], summary="Show all Databases, Collection and number of Documents",)
async def show_all_databases_and_collections():
    client = MongoClient('mongodb', 27017)
    databases = client.list_database_names()
    result = {}
    for database in databases:
        result[database] = {}
        collections = client[database].list_collection_names()
        for collection in collections:
            result[database][collection] = client[database][collection].find().count()
    return result