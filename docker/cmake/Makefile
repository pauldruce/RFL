TIMESTAMP=$(shell date +%Y%m%d)


build:
	docker build . -t RFL/docker/cmake:${TIMESTAMP}
	docker tag RFL/docker/cmake:${TIMESTAMP} RFL/docker/cmake:latest