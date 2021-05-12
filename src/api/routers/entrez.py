from fastapi import Body, APIRouter
from pydantic import BaseModel
from typing import Optional
import json

router = APIRouter()

################################################
### MODELS
################################################

class EntrezQuery(BaseModel):
    query: str
    email: str

class EntrezSingleUID(BaseModel):
    uid: str
    email: str

class EntrezListUID(BaseModel):
    uid_list: list
    email: str

################################################
### QUERY
################################################


@router.post("/entrez/query", tags=["Entrez"], summary="Query the Entrez API")
async def entrez_query(query_body: EntrezQuery = Body(...,example={ "query": "kruse eiken vestergaard", 'email' : "XXX@YYY.com" } ) ):

    def entrez_search_pubmed(query,records_per_query=10,email="XXX@YYY.com"):
        from Bio import Entrez
        Entrez.email = email
        # Search
        handle = Entrez.esearch(db="pubmed",term=query, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        return record

    results = entrez_search_pubmed(query = query_body.query, email = query_body.email )
    return results

################################################
### SUMMARY
################################################

@router.post("/entrez/summary/single", tags=["Entrez"], summary="Summary: Single UID")
async def entrez_summary_single_uid(query_body: EntrezSingleUID = Body(...,example={ "uid": "28197643", 'email' : "XXX@YYY.com" } ) ):

    def entrez_fetch_single_summary(uid,email):
        from Bio import Entrez
        Entrez.email = email
        handle = Entrez.esummary(db="pubmed", id=uid)
        record = Entrez.read(handle)
        return record

    results = entrez_fetch_single_summary(uid = query_body.uid, email = query_body.email )
    return results

@router.post("/entrez/summary/list", tags=["Entrez"], summary="Summary: List UID")
async def entrez_summary_list_of_uid(query_body: EntrezListUID = Body(...,example={ "uid_list": ["28197643","29679305","27848006"], 'email' : "XXX@YYY.com" } ) ):

    def entrez_fetch_list_summary(uid_list,email):
        from Bio import Entrez
        Entrez.email = email
        results = [ Entrez.read( Entrez.esummary(db="pubmed", id=uid) ) for uid in uid_list ]
        return [{'uid':a, 'summary':b} for a,b in zip(uid_list, results)]
        return results

    results = entrez_fetch_list_summary(uid_list = query_body.uid_list, email = query_body.email )
    return results

################################################
### ABSTRACT
################################################

### Helper functions
def flatten_abstract(abstract_xml):
    abstract = ''
    for abstractText in abstract_xml.find_all('abstracttext'):
        if abstractText.get('label') != None:
            abstract = abstract + " " + abstractText.get('label') + ": "
        abstract = abstract + abstractText.text
    return abstract

def entrez_fetch_abstracts(uid,email):
    from Bio import Entrez
    from bs4 import BeautifulSoup as bs
    Entrez.email = email
    handle = Entrez.efetch(db="pubmed", id=uid, rettype='Medline', retmode='xml')
    result = handle.readlines()
    result = b"".join(result)
    bs_content = bs(result, "lxml")
    abstracts = bs_content.find_all('abstract')
    handle.close()
    # Abstract
    results = [ flatten_abstract(abstract) for abstract in abstracts]
    return results

### Endpoints

@router.post("/entrez/abstract/single", tags=["Entrez"], summary="Abstract: List UID")
async def entrez_construct_abstract_dict(query_body: EntrezSingleUID = Body(...,example={ "uid": "28197643", 'email' : "XXX@YYY.com" } ) ):

    results = entrez_fetch_abstracts(uid = query_body.uid, email = query_body.email )
    return results

@router.post("/entrez/abstract/list", tags=["Entrez"], summary="Abstract: List of UIDs")
async def entrez_construct_abstract_dict(query_body: EntrezListUID = Body(...,example={ "uid_list": ["28197643","29679305","27848006"], 'email' : "XXX@YYY.com" } ) ):

    def entrez_construct_abstract_dict(uids,email):
        results = entrez_fetch_abstracts( uids, email )
        return [{'uid':a, 'abstract':b} for a,b in zip(uids, results)]

    results = entrez_construct_abstract_dict(uids = query_body.uid_list, email = query_body.email )
    return results

################################################
### ELINK
################################################

@router.post("/entrez/elink/single", tags=["Entrez"], summary="Elink: Single UID")
async def entrez_query(query_body: EntrezSingleUID = Body(...,example={ "uid": "28197643", 'email' : "XXX@YYY.com" } ) ):

    def entrez_fetch_full_text_linkout(uid):
        import requests
        result = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=" + uid + "&cmd=prlinks&retmode=json").json()
        return result

    results = entrez_fetch_full_text_linkout(uid = query_body.uid )
    return results

@router.post("/entrez/elink/list", tags=["Entrez"], summary="Elink: List UID")
async def entrez_query(query_body: EntrezListUID = Body(...,example={ "uid_list": ["28197643","29679305","27848006"], 'email' : "XXX@YYY.com" } ) ):

    def entrez_fetch_full_text_linkout(uid_list):
        import requests
        query = ",".join(uid_list)
        result = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=" + query + "&cmd=prlinks&retmode=json").json()
        return result

    results = entrez_fetch_full_text_linkout(uid_list = query_body.uid_list )
    return results