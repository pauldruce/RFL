TIMESTAMP=$(shell date +%Y%m%d)

docker-build:
	docker build . -t RFL/rfl:${TIMESTAMP}
	docker tag RFL/rfl:${TIMESTAMP} RFL/rfl:latest