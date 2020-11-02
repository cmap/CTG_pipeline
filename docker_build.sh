#!/usr/bin/env bash

#change the version number for each new build
docker build -t cmap/ctg:latest -t cmap/ctg:v0.1.0 --rm=true .
