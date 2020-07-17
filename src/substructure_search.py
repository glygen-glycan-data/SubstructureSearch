
import os
import sys
import time, datetime
import json
import flask
import atexit
import hashlib
import ConfigParser
import multiprocessing, Queue

import pygly.alignment
from pygly.GlycanFormatter import WURCS20Format, GlycoCTFormat


# Default Configuration
flask_API_port = 10980
flask_API_host = "localhost" # "0.0.0.0"

worker_num = 1
max_motif_size = 10
structure_file_path = ""
result_file_path = ""


# Error classes
class SubstructureSearchError(RuntimeError):

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class ParameterError(SubstructureSearchError):

    pass



# Handle parameters and configuration file
if len(sys.argv) > 1:

    kvpara = {}
    for k,v in zip(sys.argv[1::2], sys.argv[2::2]):
        if not k.startswith("-"):
            raise ParameterError("Unknown parameter %s" % k)
        k = k.lstrip("-")
        kvpara[k] = v

    if "c" not in kvpara:
        raise ParameterError("No config file provided")
    else:
        currentdir = os.path.dirname(sys.argv[0])
        currentdirabs = os.path.abspath(currentdir)

        configpath = os.path.join(currentdirabs, kvpara["c"])

        config = ConfigParser.SafeConfigParser()
        config.readfp(open(configpath))

        worker_num = config.get("substructure_search", "cpu_core")
        max_motif_size = config.get("substructure_search", "max_motif_size")

        structure_file_path = config.get("substructure_search", "glycan_set")
        structure_file_path = os.path.join(currentdirabs, structure_file_path)

        flask_API_host = config.get("service", "host")
        flask_API_port = config.get("service", "port")
        # result_file_path = config.get("service", "result_file_log")

        max_motif_size = int(max_motif_size)
        worker_num = int(worker_num)
        flask_API_port = int(flask_API_port)

else:
    raise ParameterError("No config file provided")






# Define functions for flask process
def flask_API_init(shared_resources, flask_API_host, flask_API_port):

    task_queue, result_queue = shared_resources

    results = {}


    app = flask.Flask(__name__)

    # app.config["DEBUG"] = True

    @app.route('/', methods=['GET', 'POST'])
    def home():
        return open(os.path.join(currentdirabs, "index.html")).read()


    @app.route('/date', methods=['GET', 'POST'])
    def getdate():
        return flask.jsonify(datetime.datetime.now())


    @app.route('/queue', methods=['GET', 'POST'])
    def getqueuelength():
        # update_results()
        n = len(filter(lambda x: not x["finished"], results.values()))
        return flask.jsonify(n)


    @app.route('/submit', methods=['GET', 'POST'])
    def submit():

        motif_match_position = "anywhere"
        additional_subst = False
        loose_root_match = False
        query_sequence = ""

        if flask.request.method == "GET":
            para = flask.request.args
        elif flask.request.method == "POST":
            para = flask.request.form
        else:
            return flask.jsonify("METHOD %s is not suppoted" % flask.request.method)

        if "seq" not in para:
            flask.jsonify("Please provide a valid sequence")
        query_sequence = str(para["seq"])

        if "motif_match_position" in para:
            motif_match_position = para["motif_match_position"].lower()
            if motif_match_position not in ["anywhere", "reo", "notre", "fullstructure"]:
                raise ParameterError("motif match position is not recognized")

        if "additional_subst" in para:
            if para["additional_subst"] == 'true':
                additional_subst = True

        if "loose_root_match" in para:
            if para["loose_root_match"] == 'true':
                loose_root_match = True

        tmp = query_sequence + "_" + str(motif_match_position) + "_" + str(additional_subst) + "_" + str(loose_root_match)
        list_id = hashlib.sha256(tmp).hexdigest()

        task = {
            "id": list_id,
            "seq": query_sequence,
            "motif_match_position": motif_match_position,
            "additional_subst": additional_subst,
            "loose_root_match": loose_root_match
        }
        status = {
            "id": list_id,
            "submission_detail": task,
            "finished": False,
            "result": {}
        }

        if list_id in results:
            pass
        else:
            task_queue.put(task)
            results[list_id] = status

        return flask.jsonify({
            "list_id": list_id
        })


    def update_results():

        for i in range(5):
            try:
                res = result_queue.get_nowait()
                results[res["id"]]['finished'] = True
                results[res["id"]]["result"] = res
            except Queue.Empty:
                break
            except KeyError:
                print "Job ID %s is not present" % res["id"]


    @app.route('/retrieve', methods=['GET', 'POST'])
    def retrieve():
        update_results()
        rjson = {"error": "Please provide valid list_id"}
        if flask.request.method == "GET":
            para = flask.request.args
        elif flask.request.method == "POST":
            para = flask.request.form
        else:
            return flask.jsonify({"error": "METHOD %s is not suppoted" % flask.request.method})

        if "list_id" in para:
            list_id = para["list_id"]
            rjson = results.get(list_id, {"error": list_id + " not found"})

        return flask.jsonify(rjson)


    app.run(host=flask_API_host, port=flask_API_port, threaded=False)



# Define functions for computing process
def substructure_search_init(shared_resources, structure_list_file_path):
    task_queue, result_queue = shared_resources

    gp = GlycoCTFormat()
    wp = WURCS20Format()

    motif_match_connected_nodes_cache = pygly.alignment.ConnectedNodesCache()
    mm1 = pygly.alignment.GlyTouCanMotif(connected_nodes_cache=motif_match_connected_nodes_cache)
    mm2 = pygly.alignment.MotifAllowOptionalSub(connected_nodes_cache=motif_match_connected_nodes_cache)
    mm3 = pygly.alignment.MotifLooseRoot(connected_nodes_cache=motif_match_connected_nodes_cache)
    mm4 = pygly.alignment.MotifAllowOptionalSubAndLooseRoot(connected_nodes_cache=motif_match_connected_nodes_cache)


    glycans = {}
    for line in open(structure_list_file_path):
        acc, s = line.strip().split()
        glycans[acc] = wp.toGlycan(s)

    while True:
        task_detail = task_queue.get(block=True)

        seq = task_detail["seq"]
        loose_root_match = task_detail["loose_root_match"]
        jobid = task_detail["id"]
        additional_subst = task_detail["additional_subst"]

        motif_match_position = task_detail["motif_match_position"]

        motif_matcher = mm1
        if additional_subst:
            motif_matcher = mm2
        if loose_root_match:
            motif_matcher = mm3
        if loose_root_match and additional_subst:
            motif_matcher = mm4


        fullstructure = False
        rootOnly = False
        anywhereExceptRoot = False
        if motif_match_position == "anywhere":
            pass
        elif motif_match_position == "reo":
            rootOnly = True
        elif motif_match_position == "notre":
            anywhereExceptRoot = True
        elif motif_match_position == "fullstructure":
            rootOnly = True
            fullstructure = True
        else:
            pass



        matches = []
        error = []
        calculation_start_time = time.time()

        if loose_root_match and not rootOnly:
            error.append("Loose Root Match is only available for reducing-end alignment")

        try:
            if "RES" in seq:
                motif = gp.toGlycan(seq)
            elif "WURCS" in seq:
                motif = wp.toGlycan(seq)
            else:
                raise RuntimeError
        except:
            error.append("Unable to parse")


        if len(error) == 0:
            motif_node_num = len(list(motif.all_nodes()))
            if motif_node_num > max_motif_size:
                error.append("Motif is too big")

        # TODO time out mechanism to avoid running for too long
        for acc, glycan in glycans.items():

            if len(error) != 0:
                break

            if fullstructure:
                if motif_node_num != len(list(glycan.all_nodes())):
                    continue

            if motif_matcher.leq(motif, glycan, rootOnly=rootOnly, anywhereExceptRoot=anywhereExceptRoot):
                matches.append(acc)


        calculation_end_time = time.time()
        calculation_time_cost = calculation_end_time - calculation_start_time


        res = {
            "id": jobid,
            "start time": calculation_start_time,
            "end time": calculation_end_time,
            "time spent": calculation_time_cost,
            "matches": matches,
            "error": error
        }
        result_queue.put(res)




def cleanup():
    for p in worker_processor_pool:
        p.terminate()
    front_end_API_process.terminate()




if __name__ == "__main__":

    task_queue   = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    shared_resources = [task_queue, result_queue, ]

    front_end_API_process = multiprocessing.Process(target=flask_API_init, args=(shared_resources, flask_API_host, flask_API_port))
    front_end_API_process.start()

    worker_processor_pool = []

    for i in range(worker_num):
        worker_processor = multiprocessing.Process(target=substructure_search_init, args=(shared_resources, structure_file_path))
        worker_processor.start()
        worker_processor_pool.append(worker_processor)


    atexit.register(cleanup)
    while True:
        goodbye = not front_end_API_process.is_alive()
        for p in worker_processor_pool:
            if not p.is_alive():
                goodbye = True

        if goodbye:
            cleanup()
            break
        time.sleep(1)




