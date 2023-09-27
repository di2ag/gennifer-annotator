import os
import json
import time
import requests
import openai
from reasoner_pydantic import Query
from collections import defaultdict

openai.api_key_path = os.environ.get("OPENAI_API_KEY_FILE")

#if openai_api_key_filepath:
#    with open(openai_api_key_filepath, 'r') as openai_file:
#        openai.api_key = openai_file.readline().strip()
ARS_SUBMIT_ENDPOINT = 'https://ars-prod.transltr.io/ars/api/submit'
ARS_MESSAGE_ENDPOINT = 'https://ars-prod.transltr.io/ars/api/messages'
NODE_NORMALIZER_ENDPOINT = 'https://nodenormalization-sri.renci.org/get_normalized_nodes'


def get_openai_justification(source_name, target_name, directed):
    '''
    Function will call the OPENAI end point and ask for a result for this relation.
    '''
    if directed:
        system_content = "You are an expert bioinformatician. You will be provided with a source gene and a target gene, and your task \
                            is to briefly explain a possible genetic regulator relationship that the source gene may have on the target gene. If you do not \
                            think there is a relationship, you should respond with: 'No relationship found.'"
    else:
        system_content = "You are an expert bioinformatician. You will be provided with a two gene names, and your task \
                            is to briefly explain a possible genetic regulatory relationship between these genes. If you do not \
                            think there is a relationship, you should respond with: 'No relationship found.'"

    response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo",
            messages=[
                {
                    "role": "system",
                    "content": system_content
                },
                {
                    "role": "user",
                    "content": f'Source gene name: {source_name}\nTarget gene name: {target_name}.'
                }
            ],
            temperature=0.8,
            max_tokens=256,
            )
    return response['choices'][0]['message']['content']

def add_node(node_id, node_counter, node_map, nodes):
    if node_id not in node_map:
        node_key = f'n{node_counter}'
        node_counter += 1
        node_map[node_id] = node_key
        nodes[node_key] = {
                "ids": [node_id],
                'categories': ['biolink:Gene'],
                'is_set': False
                }
    return node_counter, node_map, nodes

def add_edge(e, node_map, edge_id, edges):
    if e["directed"]:
        predicate = 'biolink:affects'
    else:
        predicate = 'biolink:related_to'
    edge_key = f'e{edge_id}'
    qsubject = node_map[e["source"]["id"]]
    qobject = node_map[e["target"]["id"]]
    edges[edge_key] = {
            'subject': qsubject,
            'object': qobject,
            'predicates': [predicate]
            }
    return edges
'''
def construct_query_graph(data):
    node_counter = 0
    node_map = {}
    # Build query nodes
    nodes = {}
    for d in data:
        node_counter, node_map, nodes = add_node(
                d["source"]["id"],
                node_counter,
                node_map, 
                nodes
                )
        node_counter, node_map, nodes = add_node(
                d["target"]["id"],
                node_counter,
                node_map, 
                nodes
                )
    # Build query edges
    edges = {}
    for edge_id, d in enumerate(data):
        edges = add_edge(d, node_map, edge_id, edges)
    # Construct query
    q = {
            'message': {
                'query_graph': {
                    'nodes': nodes,
                    'edges': edges,
                    }
                }
            }
    return q
'''

def construct_query_graph(data, directed):
    subject_ids = set()
    object_ids = set()
    for d in data:
        subject_ids.add(d["source"]["id"])
        object_ids.add(d["target"]["id"])
    if directed:
        predicate = 'biolink:affects'
    else:
        predicate = 'biolink:related_to'
    # Construct query
    q = {
            'message': {
                'query_graph': {
                    'nodes': {
                        'n00': {
                            "ids": list(subject_ids),
                            'categories': ['biolink:Gene'],
                            'is_set': False
                            },
                        'n01': {
                            "ids": list(object_ids),
                            'categories': ['biolink:Gene'],
                            'is_set': False
                            },
                        },
                    'edges': {
                        "e00": {
                            'subject': 'n00',
                            'object': 'n01',
                            'predicates': [predicate]
                            }
                        }
                    }
                }
            }
    return q

def get_translator_results(data, directed, timeout=None):
    q = construct_query_graph(data, directed)
    start_time = time.time()
    resp = requests.post(ARS_SUBMIT_ENDPOINT, json=q)
    pk = resp.json()["pk"]
    print(f'ARS PK: {pk}')
    result_obj = {
            "query_pk": pk,
            "merged_pk": None,
            "result": None,
            "status": None,
            "message": None,
            }
    # Run until ARS is Done or we timeout
    while True:
        res = requests.get(os.path.join(ARS_MESSAGE_ENDPOINT, pk))
        status = res.json()["fields"]["status"]
        result_obj["status"] = status
        if timeout:
            if time.time() - start_time > timeout:
                result_obj["message"] = 'Timed Out.'
                break
        if status == 'Running':
            time.sleep(10)
            continue
        if status != 'Running' and status != 'Done':
            result_obj['message'] = 'Recieved a strange ARS status, so breaking.'
            break
        result_obj["merged_pk"] = res.json()["fields"]["merged_version"]
        break
    # Check for result and collect.
    if result_obj["merged_pk"]:
       merged_res = requests.get(os.path.join(ARS_MESSAGE_ENDPOINT, result_obj["merged_pk"]))
       result_obj["result"] = Query.parse_obj(merged_res.json()["fields"]["data"])
    return result_obj 

def create_normalization_map(nodes):
    resp = requests.post(NODE_NORMALIZER_ENDPOINT, json={"curies": nodes})
    normed_dict = defaultdict(set)
    for curie, info in resp.json().items():
        normed_dict[curie].add(curie)
        for equiv_id_info in info["equivalent_identifiers"]:
            normed_dict[curie].add(equiv_id_info["identifier"])
    # Create map from any potential normed node to the preferred curie that was passed
    normed_map = {}
    for passed_curie, normed_curie_set in normed_dict.items():
        normed_map[passed_curie] = passed_curie
        for normed_curie in normed_curie_set:
            normed_map[normed_curie] = passed_curie
    return normed_map

def parse_translator_results(data, result, directed):
    all_curies = set.union(set([d["source"]["id"] for d in data]), set([d["target"]["id"] for d in data]))
    normed_map = create_normalization_map(list(all_curies))
    # Create a source, target to subject, object map
    if directed:
        data_map = {(d["source"]["id"], d["target"]["id"]): idx for idx, d in enumerate(data)}
    else:
        data_map = {tuple(sorted([d["source"]["id"], d["target"]["id"]])): idx for idx, d in enumerate(data)}
    for r in result.message.results:
        for analysis in r.analyses:
            for q_edge_key, edge_bindings in analysis.edge_bindings.items():
                for edge_binding in edge_bindings:
                    kedge = result.message.knowledge_graph.edges[edge_binding.id]
                    try:
                        mapped_subject = normed_map[kedge.subject]
                        mapped_object = normed_map[kedge.object]
                    except KeyError:
                        continue
                    if directed:
                        subject_object_tup = (mapped_subject, mapped_object)
                    else:
                        subject_object_tup = tuple(sorted([mapped_subject, mapped_object]))
                    if subject_object_tup not in data_map:
                        continue
                    qualified_predicate = None
                    object_modifier = None
                    object_aspect = None
                    primary_source = None
                    publications = None
                    # Collect qualifiers that matter to gennifer
                    if kedge.qualifiers:
                        for qualifier in kedge.qualifiers:
                            if qualifier.qualifier_type_id == 'biolink:qualified_predicate':
                                qualified_predicate = qualifier.qualifier_value
                            elif qualifier.qualifier_type_id == 'biolink:object_aspect_qualifier':
                                object_aspect = qualifier.qualifier_value
                            elif qualifier.qualifier_type_id == 'biolink:object_modifier_qualifier':
                                object_modifier = qualifier.qualifier_value
                    # Collect sources that matter to gennifer
                    if kedge.sources:
                        for source in kedge.sources:
                            if source.resource_role == 'primary_knowledge_source':
                                primary_source = source.resource_id
                    # Collect attributes that matter to gennifer
                    if kedge.attributes:
                        for attribute in kedge.attributes:
                            if attribute.attribute_type_id == 'biolink:publications':
                                if type(attribute.value) == str:
                                    publications = [attribute.value]
                                else:
                                    publications = list(attribute.value)
                    # Build result object
                    if 'results' not in data[data_map[subject_object_tup]]:
                        data[data_map[subject_object_tup]]["results"] = []
                    data[data_map[subject_object_tup]]["results"].append(
                            {
                                "predicate": kedge.predicate,
                                "qualified_predicate": qualified_predicate,
                                "object_modifier": object_modifier,
                                "object_aspect": object_aspect,
                                "resource_id": analysis.resource_id,
                                "primary_source": primary_source,
                                "publications": publications,
                                }
                            )
    # Add empty results key to edges that no results were found.
    for idx, item in enumerate(data):
        if 'results' not in item:
            data[idx]["results"] = []
    return data
