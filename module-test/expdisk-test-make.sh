#!/bin/sh

gcc ../src/cps/expdisk.c ../src/component.c ../src/utils.c ../src/model.c ../src/spectrum.c ../src/recipe.c expdisk-test.c -o expdisk-test -lm
