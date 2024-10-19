#!/bin/bash

ENV_DIR=$(which python) 
ENV_DIR=${ENV_DIR%/bin/python}
ETC_DIR=$ENV_DIR"/etc/conda"
AMETHYST_DIR=${PWD%misc}
echo $AMETHYST_DIR

mkdir -p $ETC_DIR
mkdir -p $ETC_DIR"/activate.d"
mkdir -p $ETC_DIR"/deactivate.d"

# Edit activation script
echo  "export OLD_PYTHONPATH=$PYTHONPATH" > $ETC_DIR"/activate.d/env_vars_2.sh"
echo  "export PYTHONPATH=$PYTHONPATH:"/home/mirto/amethyst >> $ETC_DIR"/activate.d/env_vars_2.sh"

# Edit deactivation script
echo  "export PYTHONPATH=$OLD_PYTHONPATH" > $ETC_DIR"/deactivate.d/env_vars_2.sh"
echo  "export OLD_PYTHONPATH=''" >> $ETC_DIR"/deactivate.d/env_vars_2.sh"

echo "Done"

