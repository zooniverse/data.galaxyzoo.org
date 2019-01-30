FROM node:9

WORKDIR /usr/src/

COPY package.json /usr/src/

RUN npm install

COPY . /usr/src/
