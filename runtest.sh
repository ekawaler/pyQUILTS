#!/bin/bash

python quilts.py -—output_dir /samples/em_test —-genome /samples/em_test/genome —-proteome /samples/em_test/proteome —-germline /samples/em_test/germline --somatic /samples/em_test/somatic --junction /samples/em_test/junctions —-junction_file_type mapsplice —-threshB 1 --threshD 1 --threshN 1 --fusion /samples/em_test/fusions