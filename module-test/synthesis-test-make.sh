#!/bin/sh

gcc synthesis-test.c ../src/cps/expdisk.c ../src/component.c ../src/model.c ../src/recipe.c ../src/spectrum.c ../src/synthesis.c ../src/utils.c -o synthesis-test -lm
