#!/bin/bash

find . -name \* \
! -name \*.pdf \
! -name \*.tex \
! -name \*.sty \
! -name \*.sh \
-print -delete
