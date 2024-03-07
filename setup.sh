#!/bin/bash

# Get the current dir.
if [ -n "$BASH_VERSION" ]; then
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
elif [ -n "$ZSH_VERSION" ]; then
    DIR=${0:a:h}  # https://unix.stackexchange.com/a/115431
else
	echo "Error: Unknown shell; cannot determine path to BridgeRNA2024 local repository"
fi

export BRNA2024_REPO_DIR="${DIR}"
DOCKER_IMAGE="brna2024"

# docker aliases
alias brna2024_docker_build="bash ${BRNA2024_REPO_DIR}/Docker_build.sh"

brna2024_docker_run_func() {

    PIPELINE=${1}
    THREADS=${2}
    WORKDIR=$(realpath ${3})

    docker run \
	    -v ${WORKDIR}:/WORKDIR \
	    ${DOCKER_IMAGE} sh -c "brna2024 ${PIPELINE} --threads ${THREADS} /WORKDIR"
}

alias brna2024_docker_run_interactive="docker run -it"
alias brna2024_docker_run="brna2024_docker_run_func "
