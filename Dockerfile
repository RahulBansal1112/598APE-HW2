FROM rahulbansal2187/cs598ape_hw2:latest

WORKDIR /host

# Copy everything in current directory into /host
COPY . /host

CMD ["/bin/bash"]