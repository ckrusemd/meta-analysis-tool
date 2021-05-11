FROM ubuntu

RUN apt-get update
RUN apt-get -y install nginx

RUN apt-get update
RUN apt-get install --no-install-recommends -y python3.8 python3-pip python3.8-dev

RUN update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1

EXPOSE 80