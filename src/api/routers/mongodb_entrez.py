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
################################################


@router.post("/mongodb_entrez/save_query_uids/{project_id}", tags=["MongoDB Entrez"], summary="Perform an Entrez query and save query UIDs to MongoDB")
async def save_query_results_to_mongodb(project_id:int, query_body: EntrezQuery = Body(...,example={ "query": "kruse eiken vestergaard", 'email' : "XXX@YYY.com" } ) ):
    client = Helper_.getMongoDbClient()
    db = client['Entrez']

    # Query
    results_query = Helper_.entrez_search_pubmed(query = query_body.query, email = query_body.email )
    results_query_parsed = Helper_.entrez_search_pubmed_parse_to_uid_list(results_query)
    uid_list = results_query['IdList']
    uid_list = { 'project_id': project_id, 'uid_list' : uid_list }

    # UID reviews
    client = Helper_.getMongoDbClient()
    db = client['Projects']
    project_uids = db['project_uids']

    logger.info( len( list( db['project_uids'].find({"project_id": project_id},{"_id":0}) ) ) )

    # Insert or update
    if len( list( db['project_uids'].find({"project_id": project_id},{"_id":0}) ) )==0:
        logger.info(" Insert ")
        project_uids.insert_one( uid_list )
    else:
        logger.info(" Update ")
        project_uids.delete_many( {"project_id": project_id} )
        project_uids.insert_one( uid_list )


    return list( db['project_uids'].find({},{"_id":0}) )




@router.post("/mongodb_entrez/save_query_results/", tags=["MongoDB Entrez"], summary="Perform an Entrez query and save results to MongoDB")
async def save_query_results_to_mongodb(query_body: EntrezQuery = Body(...,example={ "query": "kruse eiken vestergaard", 'email' : "XXX@YYY.com" } ) ):
    client = Helper_.getMongoDbClient()
    db = client['Entrez']

    # Query
    results_query = Helper_.entrez_search_pubmed(query = query_body.query, email = query_body.email )
    results_query_parsed = Helper_.entrez_search_pubmed_parse_to_uid_list(results_query)
    uid_list = results_query['IdList']
    logger.info(" Query ")

    # UIDs
    results_existing = Helper_.mongodb_find_all_uid()
    results_combined = results_query_parsed + results_existing
    results_combined = list( {v['uid']:v for v in results_combined}.values() )
    db['uid_list'].remove()
    db['uid_list'].insert_many( results_combined )
    logger.info(" UIDs ")


    # Summaries
    results_existing = Helper_.mongodb_find_all_summaries()
    results_summaries = Helper_.entrez_fetch_list_summary(uid_list = uid_list)
    results_combined = results_summaries + results_existing
    results_combined = list( {v['uid']:v for v in results_combined}.values() )
    db['summaries'].remove()
    db['summaries'].insert_many( results_combined )
    logger.info(" Summaries ")


    # Abstracts
    results_existing = Helper_.mongodb_find_all_abstracts()
    results_abstracts = Helper_.entrez_construct_abstract_dict(uids = uid_list)
    results_combined = results_abstracts + results_existing
    results_combined = list( {v['uid']:v for v in results_combined}.values() )
    db['abstracts'].remove()
    db['abstracts'].insert_many( results_combined )
    logger.info(" Abstracts ")


    # Elink
    results_existing = Helper_.mongodb_find_all_elinks()
    results_elinks = Helper_.entrez_fetch_full_text_linkout(uid_list)
    results_combined = results_elinks + results_existing
    results_combined = list( {v['uid']:v for v in results_combined}.values() )
    db['elinks'].remove()
    db['elinks'].insert_many( results_combined )
    logger.info(" Elinks ")


    # Counts
    counts = {  'uid_list' : db['uid_list'].find().count(),
            'summaries' : db['summaries'].find().count(),
            'abstracts' : db['abstracts'].find().count(),
            'elinks' : db['elinks'].find().count()
    }

    return counts

