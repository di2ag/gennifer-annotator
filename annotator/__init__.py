import os
import json
import requests_cache

from flask import Flask
from flask_restful import Resource, Api, reqparse
from celery.result import AsyncResult

from .tasks import create_annotation_task

def create_app(test_config=None):
# Initialise environment variables

    # Read secret key
    secret_key_filepath = os.environ.get("SECRET_KEY_FILE", None)
    if secret_key_filepath:
        with open(secret_key_filepath, 'r') as sk_file:
            SECRET_KEY = sk_file.readline().strip()
    else:
        SECRET_KEY = 'dev'

    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
            SECRET_KEY=SECRET_KEY,
            )
    if test_config is None:
        # Load this instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # Load the testing config if passed
        app.config.from_mapping(test_config)

    # Ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # Install a requests cache
    requests_cache.install_cache('annotator_cache')
    
    # Build flask RESTful API
    api = Api(app)

    parser = reqparse.RequestParser()
    parser.add_argument('data', type=list, location='json')
    parser.add_argument('directed', type=bool)

    class Annotate(Resource):
        
        def post(self):
            args = parser.parse_args()
            task = create_annotation_task.delay(args["data"], args["directed"])
            return {"task_id": task.id}, 200

        def get(self, task_id):
            task_result = AsyncResult(task_id)
            result = {
                    "task_id": task_id,
                    "task_status": task_result.status,
                    "task_result": task_result.result,
                    }
            return result, 200

            
    api.add_resource(Annotate, '/run', '/status/<task_id>')

    return app
