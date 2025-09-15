#!/bin/bash

echo "starting", $HQ_ENTRY >> "already_run.log"
$HQ_ENTRY &
