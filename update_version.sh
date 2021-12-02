#!/usr/bin/env bash

#meta.yaml
sed -Ei.bak 's|set version = "[0-9v.]+"|set version = "'"$1"'"|' conda_build/meta.yaml

#setup.py
sed -Ei.bak "s|version='[0-9v.]+'|version='$1'|" setup.py

#Dockerfile
sed -Ei.bak 's|flippyr=[0-9v.]+|flippyr='"$1"'|' Dockerfile
