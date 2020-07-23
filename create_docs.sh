#!/bin/bash
# This Unix script generates pyredcom's documentation using pdoc3
# In order to run this script correctly, you need pdoc3 first
cd pyredcom
pdoc3 pyredcom --output-dir ../docs --force --html
