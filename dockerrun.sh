docker run -it rahulbansal2187/cs598ape_hw2 /bin/bash

# docker run -it --security-opt seccomp=unconfined -v `pwd`:/host wsmoses/598ape /bin/bash
# docker run -it -v "$(pwd)":/workspace -w /workspace rahulbansal2187/cs598ape_hw2 /bin/bash
# docker run -it \
#     -v "$(pwd)":/host \
#     --security-opt seccomp=unconfined \
#     rahulbansal2187/cs598ape_hw2 \
#     /bin/bash