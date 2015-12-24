#!/bin/sh

gcc ../src/cps/sersic.c ../src/component.c ../src/utils.c ../src/model.c ../src/spectrum.c ../src/recipe.c sersic-test.c -o sersic-test -lm
