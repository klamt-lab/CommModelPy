#!/bin/bash
# This Unix script generates commodepy's documentation using pdoc3
# In order to run this script correctly, you need pdoc3 first
cd commodepy
pdoc3 commodepy --output-dir ../docs --force --html
