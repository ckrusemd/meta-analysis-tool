from fastapi import Body, APIRouter
from pydantic import BaseModel
from typing import Optional

router = APIRouter()

@router.get("/test/HelloWorld", tags=["Test"], summary="Receive Hello World")
async def hello_world():
    return "Hello World"

@router.get("/test/path/{msg}", tags=["Test"], summary="Test Path Parameter")
async def test_path_parameter(msg):
    return msg

@router.get("/test/path/type/{msg}", tags=["Test"], summary="Test Path Parameter Type")
async def test_path_parameter(msg: str):
    return { msg : isinstance(msg, str) }
    

class HelloWorld(BaseModel):
    msg: str

@router.post("/test/Models", tags=["Test"], response_model=HelloWorld, response_model_include={"msg"}, summary="Test Pydantic Model Body and Response Model")
async def test_model(helloWorld: HelloWorld = Body(...,example={ "msg": "Hello World" } ) ):
    return helloWorld

@router.post("/test/ModelsAndPathparameter/{someParameter}", tags=["Test"], summary="Test Combination of Body and Path Parameter")
async def test_model_and_parameter(helloWorld: HelloWorld = Body(...,example={ "msg": "Hello World" } ), someParameter: str = "Greeting"):
    return { someParameter : helloWorld }