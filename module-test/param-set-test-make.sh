#!/bin/sh

gcc ../src/cps/sersic.c ../src/cps/expdisk.c ../src/component.c ../src/utils.c ../src/model.c ../src/spectrum.c ../src/recipe.c param-set-test.c -o param-set-test -lm -lgsl -lgslcblas
