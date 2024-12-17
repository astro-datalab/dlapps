FROM gcc:14.2 AS build

# Sync the necessary files for build
WORKDIR /home/fitz/dlapps
COPY . .

# install additional packages
RUN apt-get update && apt-get install -y libpq-dev

# compile the code with make and gcc
RUN make all

# for this stage pull httpd
FROM httpd:2.4

# we'll use the normal working directory
WORKDIR /usr/local/apache2/

# mount point for other CGI scripts
VOLUME /usr/local/apache2/cgi-bin

# Install dependencies that the image doesn't come with but are necessary
# for the application
RUN apt update && \
    apt install -y libpng-dev && \
    apt install -y libpq-dev

# sync the configuration file
COPY ./server/* /usr/local/apache2/conf

# sync the build to this stage
COPY --from=build /home/fitz/dlapps/bin/* /usr/local/apache2/cgi-bin

EXPOSE 80
