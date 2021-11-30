FROM node:16-alpine

WORKDIR /usr/src/

COPY package.json /usr/src/
COPY package-lock.json /usr/src/

RUN npm ci

COPY . /usr/src/
