#!/usr/bin/env sh
set -e -u

cd "$(dirname "$0")"

echo "[systemtests] Building tutorial images..."
docker compose --file docker-compose.tutorial.yaml build

echo "[systemtests] Running tutorial containers..."
docker compose --file docker-compose.tutorial.yaml up

if [ -f docker-compose.field_compare.yaml ]; then
  echo "[systemtests] Running fieldcompare..."
  docker compose --file docker-compose.field_compare.yaml up --exit-code-from field-compare
fi
