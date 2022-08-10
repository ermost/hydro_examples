# This is an example for a simple Docker container.
#
# Build this Dockerfile like this (or give it a different name with "-t"):
#
#     $ docker build -f cpp.Dockerfile -t yourname/cpp:latest .
#
# Then you can run it interactively like this:
#
#     $ docker run -it yourname/cpp:latest /bin/bash
#
# In the container you can run the compiled executable:
#
#     $ ./hello_world
#
# You can push the container to your account on Dockerhub like this:
#
#     $ docker push yourname/cpp:latest

FROM ubuntu:latest

WORKDIR /root

RUN apt-get -y update && apt-get install -y clang
