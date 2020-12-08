REM This Windows batch script generates commmodelpy's documentation using pdoc3
REM In order to run this script correctly, you need pdoc3 first
pdoc3 commmodelpy --output-dir ../docs --force --html
