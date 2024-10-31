FINAL_IMAGE_NAME=RFL/rfl
BASE_IMAGE_NAME=RFL/rfl-base
CMAKE_IMAGE_NAME=RFL/rfl-cmake
TIMESTAMP=$(shell date +%Y%m%d-%H%M%S)

base-docker:
	docker build --target base -t $(BASE_IMAGE_NAME):$(TIMESTAMP) .

cmake-docker:
	docker build --target cmake -t $(CMAKE_IMAGE_NAME):$(TIMESTAMP) .

rfl-docker:
	docker build --target full -t $(FINAL_IMAGE_NAME):$(TIMESTAMP) .
	docker tag $(FINAL_IMAGE_NAME):$(TIMESTAMP) $(FINAL_IMAGE_NAME):latest

all: base-docker cmake-docker rfl-docker

.PHONY: base-docker cmake-docker rfl-docker all