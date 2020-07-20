# User Manual

## Basics
1. Run package/init.sh to start the Flask API and preparing the worker processes.
2. Run package/submit.sh to submit motif.
## Submit.sh
1. Supported Parameters
```
-s glycan sequence(*) file
-seq glycan sequence(*)


-pos motif match position, allowed value ["anywhere", "reo", "notre", "fullstructure"]
-as allow additional substituent (Sulfate or Phosphate), allowed value ["true", "false"]
-lr loose root matching, allowed value ["true", "false"]

-host http://localhost as default
-port 10980 as default

*: GlycoCT and WURCS are supported, and -seq has higher priority than -s
```
Try it out:
```
./submit.sh -s example_wurcs -pos reo      > motif_G00028MO_reducing_end_only.out
./submit.sh -s example_wurcs -pos anywhere > motif_G00028MO_anywhere.out
```


## Custom Options
Located at package/service/service.ini
```
substructure_search.cpu_core: number of CPU is(are) being used for motif alignment.
substructure_search.max_motif_size: maximum motif glycan size by monosaccharide number
substructure_search.glycan_set: glycan structure tsv file name

service.host: default host localhost 
service.port: default port 10980

Note1: there are total (2 + substructure_search.cpu_core) processes running. Primary one, flaks, and deamon processes.
Note2: the larger the substructure_search.max_motif_size and substructure_search.cpu_core is, the more memory it consumes.
```

## GUI
See:
http://"service.host":"service.port".

