#!/bin/bash

find . -type f -name \* \
     ! -name \*.pdf \
     ! -name \*.tex \
     ! -name \*.sty \
     ! -name \*.bib \
     ! -name \*.sh \
     ! -name \*png\
     ! -name makefile
     -print -delete
