#!/bin/bash
# This Unix script generates commodelpy's documentation using pdoc3
# In order to run this script correctly, you need pdoc3 first
pdoc3 commodelpy --output-dir ./docs --force --html
