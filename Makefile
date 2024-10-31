TIMESTAMP=$(shell date +%Y%m%d)
IMAGE_NAME=RFL/rfl
TIMESTAMP_FILE=.docker-timestamp

# Main build target
build: $(TIMESTAMP_FILE)

# Timestamp file target with dependencies
$(TIMESTAMP_FILE): Dockerfile docker/cmake/.docker-timestamp $(wildcard RFL_source/**/*)
	docker build . -t $(IMAGE_NAME):$(TIMESTAMP)
	docker tag $(IMAGE_NAME):$(TIMESTAMP) $(IMAGE_NAME):latest
	echo $(TIMESTAMP) > $(TIMESTAMP_FILE)