import os
import time
import json
from reasoner_pydantic import Query

from celery import Celery

from .gennifer_api import (
        get_openai_justification,
        get_translator_results,
        parse_translator_results,
        )

celery = Celery(__name__)
celery.conf.broker_url = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379")
celery.conf.result_backend = os.environ.get("CELERY_RESULT_BACKEND", "redis://localhost:6379")
celery.conf.task_routes = {"create_annotation_task": {"queue": 'annotation'}}

@celery.task(name="create_annotation_task")
def create_annotation_task(data, directed):
    # Collect OpenAI justifications
    for e in data:
        source_name = e["source"]["name"]
        target_name = e["target"]["name"]
        e["justification"] = get_openai_justification(source_name, target_name, directed)
    result_obj = get_translator_results(data, directed, None)
    #with open('result.json', 'r') as res_file:
    #    result = Query.parse_obj(json.load(res_file)["fields"]["data"])
    return parse_translator_results(data, result_obj["result"], directed)
