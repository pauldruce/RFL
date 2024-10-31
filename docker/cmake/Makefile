# Define variables
TIMESTAMP=$(shell date +%Y%m%d)
IMAGE_NAME=RFL/docker/cmake
TIMESTAMP_FILE=.docker-timestamp

# Main build target
build: $(TIMESTAMP_FILE)

# Timestamp file target with dependencies
$(TIMESTAMP_FILE): Dockerfile ../base/.docker-timestamp
	docker build . -t $(IMAGE_NAME):$(TIMESTAMP)
	docker tag $(IMAGE_NAME):$(TIMESTAMP) $(IMAGE_NAME):latest
	echo $(TIMESTAMP) > $(TIMESTAMP_FILE)
