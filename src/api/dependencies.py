from fastapi import Header, HTTPException
from pymongo import MongoClient
from loguru import logger

#async def get_token_header(x_token: str = Header(...)):
#    if x_token != "fake-super-secret-token":
#        raise HTTPException(status_code=400, detail="X-Token header invalid")


#async def get_query_token(token: str):
#    if token != "jessica":
#        raise HTTPException(status_code=400, detail="No Jessica token provided")

class Helper:
    def __init__(self):
        self.data = []

    def getMongoDbClient(self):
        return MongoClient('mongodb', 27017)

    def entrez_search_pubmed(self, query,records_per_query=10,email="XXX@YYY.com"):
        from Bio import Entrez
        Entrez.email = email
        # Search
        handle = Entrez.esearch(db="pubmed",term=query, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        return record

    def entrez_search_pubmed_parse_to_uid_list(self, query_result):
        result_parsed = [dict({'uid': value}) for value in query_result['IdList'] ]
        return result_parsed
    
    def entrez_fetch_list_summary(self, uid_list,email="XXX@YYY.com"):
        from Bio import Entrez
        Entrez.email = email
        results = [ Entrez.read( Entrez.esummary(db="pubmed", id=uid) ) for uid in uid_list ]
        return [{'uid':a, 'summary':b} for a,b in zip(uid_list, results)]
        return results

    def flatten_abstract(self, abstract_xml):
        abstract = ''
        for abstractText in abstract_xml.find_all('abstracttext'):
            if abstractText.get('label') != None:
                abstract = abstract + " " + abstractText.get('label') + ": "
            abstract = abstract + abstractText.text
        return abstract

    def entrez_fetch_abstracts(self, uid,email):
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
        results = [ self.flatten_abstract(abstract) for abstract in abstracts]
        return results
        
    def entrez_construct_abstract_dict(self, uids,email="XXX@YYY.com"):
        results = self.entrez_fetch_abstracts( uids, email )
        return [{'uid':a, 'abstract':b} for a,b in zip(uids, results)]

    def entrez_fetch_full_text_linkout(self, uid_list):
        import requests
        query = ",".join(uid_list)
        result = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&id=" + query + "&cmd=prlinks&retmode=json").json()
        results_parsed = result['linksets'][0]['idurllist']
        return [{ 'uid': result['id'], 'objurls' : result['objurls']} for result in results_parsed]


    ### Query databases
    def mongodb_find_all_uid(self):
        client = self.getMongoDbClient()
        return list( client['Entrez']['uid_list'].find({},{"_id" : False }) )

    def mongodb_find_all_summaries(self):
        client = self.getMongoDbClient()
        return list( client['Entrez']['summaries'].find({},{"_id" : False }) )

    def mongodb_find_all_abstracts(self):
        client = self.getMongoDbClient()
        return list( client['Entrez']['abstracts'].find({},{"_id" : False }) )

    def mongodb_find_all_elinks(self):
        client = self.getMongoDbClient()
        return list( client['Entrez']['elinks'].find({},{"_id" : False }) )
        