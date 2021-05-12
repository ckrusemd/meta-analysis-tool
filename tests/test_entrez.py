
import requests
from loguru import logger

def test_entrez_query():
    json_ = {   "query": "kruse eiken vestergaard",  "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/query", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_summary_single():
    json_ = {  "uid": "28197643",  "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/summary/single", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_summary_list():
    json_ = { "uid_list": ["28197643","29679305","27848006"], "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/summary/list", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_abstract_single():
    json_ = {  "uid": "28197643",  "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/abstract/single", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_abstract_list():
    json_ = { "uid_list": ["28197643","29679305","27848006"], "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/abstract/list", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_elink_single():
    json_ = {  "uid": "28197643",  "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/elink/single", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )

def test_entrez_elink_list():
    json_ = { "uid_list": ["28197643","29679305","27848006"], "email": "XXX@YYY.com" }
    response = requests.post("http://api:8080/entrez/elink/list", json = json_)
    assert response.status_code == 200
    assert response.json()
    logger.info( response.json() )