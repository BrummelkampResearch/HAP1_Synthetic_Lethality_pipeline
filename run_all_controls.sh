#!/bin/bash
./analyze_sli.sh --createcontrol --seq=/media/raw_data/HAP1_CONTROL_1.fastq.gz --replicate=1 --name=ControlData-HAP1
./analyze_sli.sh --createcontrol --seq=/media/raw_data/HAP1_CONTROL_2.fastq.gz --replicate=2 --name=ControlData-HAP1
./analyze_sli.sh --createcontrol --seq=/media/raw_data/HAP1_CONTROL_3.fastq.gz --replicate=3 --name=ControlData-HAP1
./analyze_sli.sh --createcontrol --seq=/media/raw_data/HAP1_CONTROL_4.fastq.gz --replicate=4 --name=ControlData-HAP1
