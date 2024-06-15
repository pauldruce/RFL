
TIMESTAMP=$(shell date +%Y%m%d)

build:
	docker build . -t RFL/docker/base:${TIMESTAMP}
	docker tag RFL/docker/base:${TIMESTAMP} RFL/docker/base:latest