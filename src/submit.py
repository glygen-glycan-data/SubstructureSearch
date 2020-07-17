import sys
import time
import requests



kvpara = {}
if len(sys.argv) > 1:

    for k,v in zip(sys.argv[1::2], sys.argv[2::2]):
        k = k.lstrip("-")
        kvpara[k] = v

else:
    sys.stdout.write("Please provide sufficient parameters!")
    sys.exit()

if "s" not in kvpara and "seq" not in kvpara:
    sys.stdout.write("Please provide glycan sequence or sequence file!")
    sys.exit()

if "s" in kvpara:
    seq_fp = kvpara["s"]
    seq = open(seq_fp).read().strip()

if "seq" in kvpara:
    seq = kvpara["seq"]

pos = "anywhere"
if "pos" in kvpara:
    pos = kvpara["pos"]

as_flag = 'false'
if "as" in kvpara:
    flag = kvpara["as"]

lr_flag = 'false'
if "lr" in kvpara:
    flag = kvpara["lr"]

main_url = "http://localhost"
if "host" in kvpara:
    main_url = kvpara["host"]

port = 10980
if "port" in kvpara:
    port = kvpara["port"]


if port not in [22, "22"]:
    main_url += ":%s/"%port





params = {
    "seq": seq,
    "motif_match_position": pos,
    "additional_subst": as_flag,
    "loose_root_match": lr_flag
}


try:
    response1 = requests.post(main_url + "submit", params)
    list_id = response1.json()[u'list_id']
except Exception, e:
    sys.stdout.write("Error: has issue connecting to flask API.")
    sys.stdout.write(str(e))
    sys.exit()



params = {"list_id": list_id}
while True:
    time.sleep(1)

    response2 = requests.post(main_url+ "retrieve", params)
    result = response2.json()
    if result[u'finished']:
        if len(result[u"result"][u"error"]) > 0:
            for e in result[u"result"][u"error"]:
                sys.stdout.write("Error: has issue during computing.")
                sys.stdout.write(e)
            sys.exit()
        break

sys.stdout.write("\n".join(sorted(result[u"result"][u'matches'])) + "\n")




