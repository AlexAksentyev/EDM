#!/bin/bash

find . -type f -name \* \
     ! -name \*.pdf \
     ! -name \*.tex \
     ! -name \*.sty \
     ! -name \*.sh \
     ! -name \*png\
     -print -delete
