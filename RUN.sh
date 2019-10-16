#!/bin/bash

# Ensure script path (from: https://stackoverflow.com/a/246128/9571488)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" > /dev/null 2>&1 && pwd )"

# Run application portal
python $DIR/application/__main__.py
