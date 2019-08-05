#!/bin/bash

if [ -d ./output/$1 ] ; then
  echo "./output/$1 already exist."
else
  mkdir -v ./output/$1
  mkdir -v ./output/$1/log
  mkdir -v ./output/$1/err
  mkdir -v ./output/$1/out
fi
